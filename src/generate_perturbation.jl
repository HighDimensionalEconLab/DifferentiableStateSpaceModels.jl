
# Utility to call the cache's organized by symbol
function fill_array_by_symbol_dispatch(f, c, symbols, args...)
    for (i, sym) in enumerate(symbols)
        f(c[i], Val(sym), args...)
    end
end

# The generate_perturbation function calculates the perturbation itself
# It can do used without any derivatives overhead (except, perhaps, extra memory in the cache) 
function generate_perturbation(m::PerturbationModel, p_d, p_f, order::Val{1} = Val(1);
                               cache = SolverCache(m, Val(1), p_d),
                               settings = PerturbationSolverSettings())
    @assert cache.p_d_symbols == collect(Symbol.(keys(p_d)))

    p = isnothing(p_f) ? p_d : order_vector_by_symbols(merge(p_d, p_f), m.mod.p_symbols)

    # solver type provided to all callbacks
    ret = calculate_steady_state!(m, cache, settings, p)
    maybe_call_function(settings.calculate_steady_state_callback, ret, m, cache, settings,
                        p)  # before returning
    (ret == :Success) || return FirstOrderPerturbationSolution(ret, m, cache)

    ret = evaluate_first_order_functions!(m, cache, settings, p)
    maybe_call_function(settings.evaluate_functions_callback, ret, m, cache, settings, p)
    (ret == :Success) || return FirstOrderPerturbationSolution(ret, m, cache)

    ret = solve_first_order!(m, cache, settings)
    maybe_call_function(settings.solve_first_order_callback, ret, m, cache, settings)
    (ret == :Success) || return FirstOrderPerturbationSolution(ret, m, cache)
    return FirstOrderPerturbationSolution(:Success, m, cache)
end

# The generate_perturbation function calculates the perturbation itself
# It can do used without any derivatives overhead (except, perhaps, extra memory in the cache)
function generate_perturbation(m::PerturbationModel, p_d, p_f, order::Val{2};
                               cache = SolverCache(m, Val(2), p_d),
                               settings = PerturbationSolverSettings())
    @assert cache.p_d_symbols == collect(Symbol.(keys(p_d)))
    @assert cache.order == Val(2)

    p = isnothing(p_f) ? p_d : order_vector_by_symbols(merge(p_d, p_f), m.mod.p_symbols)

    # Calculate the first-order perturbation
    sol_first = generate_perturbation(m, p_d, p_f, Val(1); cache, settings)

    (sol_first.retcode == :Success) || return SecondOrderPerturbationSolution(sol_first.retcode, m, cache)

    # solver type provided to all callbacks
    ret = evaluate_second_order_functions!(m, cache, settings, p)
    (ret == :Success) || return SecondOrderPerturbationSolution(ret, m, cache)
    maybe_call_function(settings.evaluate_functions_callback, ret, m, cache, settings, p)
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
        buff = c.first_order_solver_buffer
        buff.A .= [c.H_xp c.H_yp]
        buff.B .= [c.H_x c.H_y]
        s = schur!(buff.A, buff.B) # Generalized Schur decomposition, inplace using buffers
        # The generalized eigenvalues λ_i are S_ii / T_ii
        # Following Blanchard-Kahn condition, we reorder the Schur so that
        # S_22 ./ T_22 < 1, ie, the eigenvalues < 1 come last
        # inds = [s.α[i] / s.β[i] >= 1 for i in 1:n]
        inds = abs.(s.α) .>= (1 - ϵ_BK) .* abs.(s.β)
        if sum(inds) != n_x
            if print_level > 0
                # @show (n_x, sum(inds), inds, abs.(s.α), abs.(s.β)) # move to print_level > 1?
                println("Blanchard-Kahn condition not satisfied\n")
            end
            return :BlanchardKahnFailure
        end

        ordschur!(s, inds)
        # In Julia A = QSZ' and B = QTZ'

        b = 1:n_x
        l = (n_x + 1):n

        # Extract the Schur components to real matrices, for inplace factorizations/etc.
        Z = s.Z'

        c.g_x .= real(-Z[l, l] \ Z[l, b])
        blob = Z[b, b] .+ Z[b, l] * c.g_x
        c.h_x .= real(-blob \ (s.S[b, b] \ (s.T[b, b] * blob)))
        # buff.Z .= real(s.Z')
        # buff.Z_ll .= buff.Z[l,l] #preallocate for buffer for inplace LU.  No known structure?

        # Both of these are upper-triangular, helpful for fast linsolve
        # buff.S_bb .= UpperTriangular(real(s.S[b,b]))
        # buff.T_bb .= UpperTriangular(real(s.T[b,b]))

        # TODO: Check if RecursiveFactorization.jl is faster or slower than using MKL/etc. Add toggle
        #Z_ll = lu!(buff.Z_ll)
        # Z_ll = RecursiveFactorization.lu!(buff.Z_ll)
        # c.g_x .= ldiv!(Z_ll, buff.Z[l, b])
        # c.g_x .*= -1

        # The following is an as-inplace-as-possible version of
        # blob = buff.Z[b, b] .+ buff.Z[b, l] * c.g_x
        # c.h_x .= -blob \ (buff.S_bb \ (buff.T_bb * blob))
        # temp = buff.Z[b, b] .+ buff.Z[b, l] * c.g_x        
        # mul!(c.h_x, buff.T_bb, temp)  # doing everything in place
        # ldiv!(buff.S_bb, c.h_x)  # no factorization required since triangular
        # temp_lu = lu!(temp)
        # ldiv!(temp_lu, c.h_x)
        # c.h_x .*= -1

        # fill in Σ, Ω.
        c.Σ .= Symmetric(c.Γ * c.Γ')

        # Q transforms
        mul!(c.C_1, c.Q, vcat(c.g_x, I))

        # Stationary Distribution
        c.η_Σ_sq .= Symmetric(c.η * c.Σ * c.η')  # used in 2nd order as well

        V = cholesky(lyapd(c.h_x, c.η_Σ_sq); check = false) #inplace wouldn't help since allocating for lyapd.  Can use lyapds! perhaps, but would need work and have low payoffs

        # no inplace assignment or copy for cholesky, so reach inside internals for now.  Assumes same uplo flag
        c.V.factors .= V.factors

        # eta * Gamma
        mul!(c.B, c.η, c.Γ)

        # @exfiltrate  # flip on to see intermediate calculations.  TURN OFF BEFORE PROFILING
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

    buff = c.second_order_solver_buffer

    try
        # "Sylvester prep for _xx"
        buff.A .= [c.H_y c.H_xp + c.H_yp * c.g_x]
        buff.C[:, 1:n_y] .= c.H_yp
        ws = IPlusAtKronBWs(n, n, n_x, 2)
        # TODO: Tullio/etc. quadratic form trickier for any of this?
        buff.R .= vcat(c.g_x * c.h_x, c.g_x, c.h_x, c.I_x)
        for i in 1:n
            buff.E[i, :] = vec(buff.R' * c.Ψ[i] * buff.R) # (24), flip the sign for gsylv
        end
        buff.E .*= -1

        # Sylvester
        generalized_sylvester_solver!(buff.A, buff.C, c.h_x, buff.E, 2, ws)
        c.g_xx .= reshape(buff.E[1:n_y, :], n_y, n_x, n_x)
        c.h_xx .= reshape(buff.E[(n_y + 1):end, :], n_x, n_x, n_x)

        # Linear equations for _σσ
        buff.A_σ .= [c.H_yp + c.H_y c.H_xp + c.H_yp * c.g_x]
        C_σ = zeros(n)
        H_yp_g = c.H_yp * buff.E[1:n_y, :]
        buff.R_σ .= vcat(c.g_x, zeros(n_y, n_x), c.I_x, zeros(n_x, n_x))
        for i in 1:n # (29), flip the sign for (34)
            C_σ[i] -= dot(buff.R_σ' * c.Ψ[i] * buff.R_σ, c.η_Σ_sq)
            C_σ[i] -= dot(H_yp_g[i, :], c.η_Σ_sq)
        end
        #A_σ_lu = lu!(buff.A_σ) # modifes the buff.A_σ
        A_σ_lu = RecursiveFactorization.lu!(buff.A_σ)
        ldiv!(A_σ_lu, C_σ) # solve (34) inplace. X_σ = C_σ after modification
        c.g_σσ .= C_σ[1:n_y]
        c.h_σσ .= C_σ[(n_y + 1):end]

        c.C_0 .= 0.5 * c.Q * vcat(c.g_σσ, zeros(n_x))

        # TODO: This looks like a tullio thing
        fill!(c.C_2, zero(eltype(c.C_2)))  # reset as we need to use `+=`
        for i in 1:n_z
            for j in 1:n_y
                c.C_2[i, :, :] += 0.5 * c.Q[i, j] * c.g_xx[j, :, :]
            end
        end
    catch e
        settings.print_level == 0 || display(e)
        return :Failure # generic failure
    end
    #    @exfiltrate  # flip on to see intermediate calculations.  TURN OFF BEFORE PROFILING
    return :Success
end
