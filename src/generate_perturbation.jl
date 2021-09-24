
# Utility to call the cache's organized by symbol
function fill_array_by_symbol_dispatch(f, c, symbols, args...)
    for (i, sym) in enumerate(symbols)
        f(c[i], Val(sym), args...)
    end
end

# The generate_perturbation function calculates the perturbation itself
# It can do used without any derivatives overhead (except, perhaps, extra memory in the cache)
function generate_perturbation(m::PerturbationModel, p_d, p_f, order::Val{1} = Val(1);
                                  cache=SolverCache(m, Val(1), p_d),
                                  settings=PerturbationSolverSettings())
    @assert cache.p_d_symbols == collect(Symbol.(keys(p_d)))

    p = isnothing(p_f) ? p_d : order_vector_by_symbols(merge(p_d, p_f), m.mod.p_symbols)

    # solver type provided to all callbacks
    solver = PerturbationSolver(m, cache, settings)

    ret = calculate_steady_state!(m, cache, settings, p)
    maybe_call_function(settings.calculate_steady_state_callback, ret, m, cache, settings,
                        p)  # before returning
    (ret == :Success) || return FirstOrderPerturbationSolution(ret, m, cache)

    ret = evaluate_first_order_functions!(m, cache, settings, p)
    maybe_call_function(settings.evaluate_functions_callback, ret, m, cache, settings, p,
                        solver)
    (ret == :Success) || return FirstOrderPerturbationSolution(ret, m, cache)

    ret = solve_first_order!(m, cache, settings)
    maybe_call_function(settings.solve_first_order_callback, ret, m, cache, settings)
    (ret == :Success) || return FirstOrderPerturbationSolution(ret, m, cache)
    return FirstOrderPerturbationSolution(:Success, m, cache)
end

# The generate_perturbation function calculates the perturbation itself
# It can do used without any derivatives overhead (except, perhaps, extra memory in the cache)
function generate_perturbation(m::PerturbationModel, p_d, p_f, order::Val{2};
                                   cache=SolverCache(m, Val(2), p_d),
                                   settings=PerturbationSolverSettings())
    @assert cache.p_d_symbols == collect(Symbol.(keys(p_d)))

    p = isnothing(p_f) ? p_d : order_vector_by_symbols(merge(p_d, p_f), m.mod.p_symbols)

    # Calculate the first-order perturbation
    sol_first = generate_perturbation(m, p_d, p_f, Val(1);cache, settings)

    (sol_first.retcode == :Success) || return SecondOrderPerturbationSolution(ret, m, cache)

    # solver type provided to all callbacks
    solver = PerturbationSolver(m, cache, settings)
    ret = evaluate_second_order_functions!(m, cache, settings, p)
    (ret == :Success) || return SecondOrderPerturbationSolution(ret, m, cache)
    maybe_call_function(settings.evaluate_functions_callback, ret, m, cache, settings, p,
                        solver)
    ret = solve_second_order!(m, cache, settings)
    (ret == :Success) || return SecondOrderPerturbationSolution(ret, m, cache)
    return SecondOrderPerturbationSolution(:Success, m, cache)
end

function calculate_steady_state!(m::PerturbationModel, c, settings, p)
    @unpack n_y, n_x = m
    n = n_y + n_x

    (settings.print_level > 2) && println("Calculating steady state")
    try
        if !isnothing(m.mod.ȳ!) && !isnothing(m.mod.x̄!) # use closed form if possible
            m.mod.ȳ!(c.y, p)
            m.mod.x̄!(c.x, p)
            isnothing(m.mod.ȳ_p!) ||
                fill_array_by_symbol_dispatch(m.mod.ȳ_p!, c.y_p, c.p_d_symbols, p)
            isnothing(m.mod.x̄_p!) ||
                fill_array_by_symbol_dispatch(m.mod.x̄_p!, c.x_p, c.p_d_symbols, p)
        elseif !isnothing(m.mod.steady_state!) # use user-provided calculation otherwise
            m.mod.steady_state!(c.y, c.x, p)
        else # fallback is to solve system of equations from user-provided initial condition
            y_0 = zeros(n_y)
            x_0 = zeros(n_x)
            m.mod.ȳ_iv!(y_0, p)
            m.mod.x̄_iv!(x_0, p)
            w_0 = [y_0; x_0]

            if isnothing(m.mod.H̄_w!) # no jacobian
                nlsol = nlsolve((H, w) -> m.mod.H̄!(H, w, p), w_0;
                                DifferentiableStateSpaceModels.nlsolve_options(settings)...)
            else
                J_0 = zeros(n, n)
                F_0 = zeros(n)
                df = OnceDifferentiable((H, w) -> m.mod.H̄!(H, w, p),
                                        (J, w) -> m.mod.H̄_w!(J, w, p), w_0, F_0, J_0)  # TODO: the buffer to use for the w_0 is unclear?
                nlsol = nlsolve(df, w_0;
                                DifferentiableStateSpaceModels.nlsolve_options(settings)...)
            end
            if !converged(nlsol)
                if settings.print_level > 0
                    println("No steady state found\n")
                end
                return :SteadyStateFailure
            end
            settings.print_level > 1 &&
                println("Steady state found in $(nlsol.iterations) iterations\n")
            c.y .= nlsol.zero[1:n_y]
            c.x .= nlsol.zero[(n_y + 1):end]
        end
    catch e
        if !is_linear_algebra_exception(e)
            (settings.print_level > 2) && println("Rethrowing exception")
            rethrow(e)
        else
            settings.print_level == 0 || display(e)
            return :Failure # generic failure
        end
    end
    return :Success
end

function evaluate_first_order_functions!(m, c, settings, p)
    (settings.print_level > 2) && println("Evaluating first-order functions into cache")
    try
        @unpack y, x = c  # Precondition: valid (y, x) steady states

        m.mod.H_yp!(c.H_yp, y, x, p)
        m.mod.H_y!(c.H_y, y, x, p)
        m.mod.H_xp!(c.H_xp, y, x, p)
        m.mod.H_x!(c.H_x, y, x, p)
        m.mod.Γ!(c.Γ, p)
        maybe_call_function(m.mod.Ω!, c.Ω, p) # supports  m.mod.Ω! = nothing
        (length(c.p_d_symbols) > 0) && m.mod.Ψ!(c.Ψ, y, x, p)
    catch e
        if !is_linear_algebra_exception(e)
            (settings.print_level > 2) && println("Rethrowing exception")
            rethrow(e)
        else
            settings.print_level == 0 || display(e)
            return :Failure # generic failure
        end
    end
    return :Success  # no failing code-paths yet.
end

function evaluate_second_order_functions!(m, c, settings, p)
    (settings.print_level > 2) && println("Evaluating second-order functions into cache")
    try
        @unpack y, x = c  # Precondition: valid (y, x) steady states
        (length(c.p_d_symbols) == 0) && m.mod.Ψ!(c.Ψ, y, x, p)  # would have been called otherwise in first_order_functions
        m.mod.Ψ!(c.Ψ, y, x, p)
        m.mod.Ψ_yp!(c.Ψ_yp, y, x, p)
        m.mod.Ψ_y!(c.Ψ_y, y, x, p)
        m.mod.Ψ_xp!(c.Ψ_xp, y, x, p)
        m.mod.Ψ_x!(c.Ψ_x, y, x, p)
    catch e
        if !is_linear_algebra_exception(e)
            (settings.print_level > 2) && println("Rethrowing exception")
            rethrow(e)
        else
            settings.print_level == 0 || display(e)
            return :Failure # generic failure
        end
    end
    return :Success  # no failing code-paths yet.
end

function solve_first_order!(m, c, settings)
    @unpack ϵ_BK, print_level = settings
    @unpack n_x, n_y, n_p, n_ϵ = m
    n = n_x + n_y

    (settings.print_level > 2) && println("Solving first order perturbation")
    try
        A = [c.H_xp c.H_yp]
        B = [c.H_x c.H_y]
        s = schur(complex(A), complex(B)) # Generalized Schur decomposition
        # The generalized eigenvalues λ_i are S_ii / T_ii
        # Following Blanchard-Kahn condition, we reorder the Schur so that
        # S_22 ./ T_22 < 1, ie, the eigenvalues < 1 come last
        # inds = [s.α[i] / s.β[i] >= 1 for i in 1:n]
        inds = abs.(s.α) .>= (1 - ϵ_BK) .* abs.(s.β)
        if sum(inds) != n_x
            # More debugging code???
            if print_level > 0
                @show n_x
                @show sum(inds)
                @show inds
                @show abs.(s.α)
                @show abs.(s.β)
                println("Blanchard-Kahn condition not satisfied\n")
            end
            return :BlanchardKahnFailure
        end

        ordschur!(s, inds)
        # In Julia A = QSZ' and B = QTZ'

        @unpack S, T = s # Extract the Schur components
        Z = s.Z'

        b = 1:n_x
        l = (n_x + 1):n
        g_x = -Z[l, l] \ Z[l, b]
        blob = Z[b, b] .+ Z[b, l] * g_x
        h_x = -blob \ (S[b, b] \ (T[b, b] * blob))
        c.g_x .= real(g_x)
        c.h_x .= real(h_x)

        # fill in Σ, Ω.
        c.Σ .= Symmetric(c.Γ * c.Γ')

        # Q transforms
        c.C_1 .= c.Q * vcat(c.g_x, diagm(ones(n_x)))

        # Stationary Distribution
        c.V = cholesky(Symmetric(lyapd(c.h_x, c.η * c.Σ * c.η')))
        # eta * Gamma
        c.B .= c.η * c.Γ
    catch e
        if !is_linear_algebra_exception(e)
            (settings.print_level > 2) && println("Rethrowing exception")
            rethrow(e)
        else
            settings.print_level == 0 || display(e)
            return :Failure # generic failure
        end
    end
    return :Success
end

function solve_second_order!(m, c, settings)
    @unpack n_x, n_y, n_p, n_ϵ, n_z, η = m
    n = n_x + n_y

    # "Sylvester prep for _xx"
    A = [c.H_y c.H_xp + c.H_yp * c.g_x]
    B = I(n_x * n_x)
    C = [c.H_yp zeros(n, n_x)]
    D = kron(c.h_x, c.h_x)
    E = zeros(n, n_x * n_x)
    R = vcat(c.g_x * c.h_x, c.g_x, c.h_x, I(n_x))
    for i in 1:n
        E[i, :] = -(R' * c.Ψ[i] * R)[:] # (24), flip the sign for gsylv
    end

    # "Sylvester"
    X = gsylv(A, B, C, D, E) # (22)
    c.g_xx .= reshape(X[1:n_y, :], n_y, n_x, n_x)
    c.h_xx .= reshape(X[(n_y + 1):end, :], n_x, n_x, n_x)

    # Linear equations for _σσ
    A_σ = [c.H_yp + c.H_y c.H_xp + c.H_yp * c.g_x]
    C_σ = zeros(n)
    η_sq = η * c.Σ * η'
    H_yp_g = c.H_yp * X[1:n_y, :]
    R_σ = vcat(c.g_x, zeros(n_y, n_x), I(n_x), zeros(n_x, n_x))
    for i in 1:n # (29), flip the sign for (34)
        C_σ[i] -= dot(R_σ' * c.Ψ[i] * R_σ, η_sq)
        C_σ[i] -= dot(H_yp_g[i, :], η_sq)
    end
    X_σ = A_σ \ C_σ # solve (34)
    c.g_σσ .= X_σ[1:n_y]
    c.h_σσ .= X_σ[(n_y + 1):end]

    c.C_0 .= 0.5 * c.Q * vcat(c.g_σσ, zeros(n_x))
    c.C_2 .= zero(eltype(c.C_2))  # reset as we need to use `+=`
    for i in 1:n_z
        for j in 1:n_y
            c.C_2[i, :, :] += 0.5 * c.Q[i, j] * c.g_xx[j, :, :]
        end
    end
    return :Success
end
