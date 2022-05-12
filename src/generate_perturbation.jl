
# Utility to call the cache's organized by symbol
function fill_array_by_symbol_dispatch(f, c, symbols, args...)
    for (i, sym) in enumerate(symbols)
        f(c[i], Val(sym), args...)
    end
end

# create it if not supplied
function create_or_zero_cache(m, cache::Nothing, ::Val{Order}, p_d,
                              zero_cache) where {Order}
    return SolverCache(m,
                       Val(Order),
                       p_d)
end

# otherwise conditionally zero it and return the argument
function create_or_zero_cache(m, cache, ::Val{Order}, p_d, zero_cache) where {Order}
    if zero_cache
        fill_zeros!(cache) # recursively works through cache and sub-types
    end
    return cache
end

# The generate_perturbation function calculates the perturbation itself
# It can do used without any derivatives overhead (except, perhaps, extra memory in the cache) 
function generate_perturbation(m::PerturbationModel, p_d, p_f, order::Val{1} = Val(1);
                               cache = nothing, zero_cache = false,
                               settings = PerturbationSolverSettings())
    c = create_or_zero_cache(m, cache, order, p_d, zero_cache)
    c.p_d_symbols .= collect(Symbol.(keys(p_d)))

    p = isnothing(p_f) ? p_d : order_vector_by_symbols(merge(p_d, p_f), m.mod.m.p_symbols)

    # solver type provided to all callbacks
    ret = calculate_steady_state!(m, c, settings, p)
    maybe_call_function(settings.calculate_steady_state_callback, ret, m, c, settings,
                        p)  # before returning
    (ret == :Success) || return FirstOrderPerturbationSolution(ret, m, c, settings)

    ret = evaluate_first_order_functions!(m, c, settings, p)
    maybe_call_function(settings.evaluate_functions_callback, ret, m, c, settings, p)
    (ret == :Success) || return FirstOrderPerturbationSolution(ret, m, c, settings)

    ret = solve_first_order!(m, c, settings)
    maybe_call_function(settings.solve_first_order_callback, ret, m, c, settings)
    (ret == :Success) || return FirstOrderPerturbationSolution(ret, m, c, settings)
    return FirstOrderPerturbationSolution(:Success, m, c, settings)
end

# The generate_perturbation function calculates the perturbation itself
# It can do used without any derivatives overhead (except, perhaps, extra memory in the cache)
function generate_perturbation(m::PerturbationModel, p_d, p_f, order::Val{2};
                               cache = nothing, zero_cache = false,
                               settings = PerturbationSolverSettings())
    c = create_or_zero_cache(m, cache, order, p_d, zero_cache)
    c.p_d_symbols .= collect(Symbol.(keys(p_d)))
    @assert c.order == Val(2)

    p = isnothing(p_f) ? p_d : order_vector_by_symbols(merge(p_d, p_f), m.mod.m.p_symbols)

    # Calculate the first-order perturbation
    sol_first = generate_perturbation(m, p_d, p_f, Val(1); cache = c, settings,
                                      zero_cache = false)

    (sol_first.retcode == :Success) ||
        return SecondOrderPerturbationSolution(sol_first.retcode, m, c, settings)

    # solver type provided to all callbacks
    ret = evaluate_second_order_functions!(m, c, settings, p)
    (ret == :Success) || return SecondOrderPerturbationSolution(ret, m, c, settings)
    maybe_call_function(settings.evaluate_functions_callback, ret, m, c, settings, p)
    ret = solve_second_order!(m, c, settings)
    (ret == :Success) || return SecondOrderPerturbationSolution(ret, m, c, settings)
    return SecondOrderPerturbationSolution(:Success, m, c, settings)
end

function calculate_steady_state!(m::PerturbationModel, c, settings, p)
    @unpack n_y, n_x = m
    n = n_y + n_x

    (settings.print_level > 2) && println("Calculating steady state")
    try
        if !isnothing(m.mod.m.ȳ!) && !isnothing(m.mod.m.x̄!) # use closed form if possible
            m.mod.m.ȳ!(c.y, p)
            m.mod.m.x̄!(c.x, p)
            isnothing(m.mod.m.ȳ_p!) ||
                fill_array_by_symbol_dispatch(m.mod.m.ȳ_p!, c.y_p, c.p_d_symbols, p)
            isnothing(m.mod.m.x̄_p!) ||
                fill_array_by_symbol_dispatch(m.mod.m.x̄_p!, c.x_p, c.p_d_symbols, p)
        elseif !isnothing(m.mod.m.steady_state!) # use user-provided calculation otherwise
            m.mod.m.steady_state!(c.y, c.x, p)
        else # fallback is to solve system of equations from user-provided initial condition
            y_0 = zeros(n_y)
            x_0 = zeros(n_x)
            m.mod.m.ȳ_iv!(y_0, p)
            m.mod.m.x̄_iv!(x_0, p)
            w_0 = [y_0; x_0]

            if isnothing(m.mod.m.H̄_w!) # no jacobian
                nlsol = nlsolve((H, w) -> m.mod.m.H̄!(H, w, p), w_0;
                                DifferentiableStateSpaceModels.nlsolve_options(settings)...)
            else
                J_0 = zeros(n, n)
                F_0 = zeros(n)
                df = OnceDifferentiable((H, w) -> m.mod.m.H̄!(H, w, p),
                                        (J, w) -> m.mod.m.H̄_w!(J, w, p), w_0, F_0, J_0)  # TODO: the buffer to use for the w_0 is unclear?
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
        if settings.rethrow_exceptions
            rethrow(e)
        elseif e isa LAPACKException || e isa PosDefException
            (settings.print_level > 0) && display(e)
            return :LAPACK_Error
        elseif e isa PosDefException
            (settings.print_level > 0) && display(e)
            return :POSDEF_EXCEPTION
        elseif e isa DomainError
            settings.print_level == 0 || display(e)
            return :Evaluation_Error # function evaluation error
        else
            settings.print_level == 0 || display(e)
            return :FAILURE # generic failure
        end
    end
    return :Success
end

function evaluate_first_order_functions!(m, c, settings, p)
    (settings.print_level > 2) && println("Evaluating first-order functions into cache")
    try
        @unpack y, x = c  # Precondition: valid (y, x) steady states

        m.mod.m.H_yp!(c.H_yp, y, x, p)
        m.mod.m.H_y!(c.H_y, y, x, p)
        m.mod.m.H_xp!(c.H_xp, y, x, p)
        m.mod.m.H_x!(c.H_x, y, x, p)
        m.mod.m.Γ!(c.Γ, p)
        maybe_call_function(m.mod.m.Ω!, c.Ω, p) # supports  m.mod.m.Ω! = nothing
        (length(c.p_d_symbols) > 0) && m.mod.m.Ψ!(c.Ψ, y, x, p)
    catch e
        if settings.rethrow_exceptions
            rethrow(e)
        elseif e isa DomainError
            settings.print_level == 0 || display(e)
            return :Evaluation_Error # function evaluation error
        else
            settings.print_level == 0 || display(e)
            return :FAILURE # generic failure
        end
    end
    return :Success  # no failing code-paths yet.
end

function evaluate_second_order_functions!(m, c, settings, p)
    (settings.print_level > 2) && println("Evaluating second-order functions into cache")
    try
        @unpack y, x = c  # Precondition: valid (y, x) steady states
        (length(c.p_d_symbols) == 0) && m.mod.m.Ψ!(c.Ψ, y, x, p)  # would have been called otherwise in first_order_functions
        m.mod.m.Ψ!(c.Ψ, y, x, p)
        m.mod.m.Ψ_yp!(c.Ψ_yp, y, x, p)
        m.mod.m.Ψ_y!(c.Ψ_y, y, x, p)
        m.mod.m.Ψ_xp!(c.Ψ_xp, y, x, p)
        m.mod.m.Ψ_x!(c.Ψ_x, y, x, p)
    catch e
        if settings.rethrow_exceptions
            rethrow(e)
        elseif e isa DomainError
            settings.print_level == 0 || display(e)
            return :Evaluation_Error # function evaluation error
        else
            settings.print_level == 0 || display(e)
            return :FAILURE # generic failure
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
        (settings.print_level > 3) && println("Calculating schur")
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
            if settings.rethrow_exceptions
                error("Failure of the Blanchard Khan Condition")
            else
                return :Blanchard_Kahn_Failure
            end
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

        # if calculating the ergodic distribution, solve the Lyapunov and use a pivoted cholesky
        if settings.calculate_ergodic_distribution
            try
                # inplace wouldn't help since allocating for lyapd.  Can use lyapds! perhaps, but would need work and have low payoffs
                c.V.mat .= lyapd(c.h_x, c.η_Σ_sq)

                # potentially perturb the covariance to try to get positive definite
                if settings.perturb_covariance > 0.0
                    c.V.mat .+= settings.perturb_covariance * I(size(c.V.mat, 1))  # perturb to ensure it is positive definite
                end
                # Do inplace cholesky and catch error.  Note that Cholesky required for MvNormal construction regardless
                copy!(c.V.chol.factors, c.V.mat) # copy over to the factors for the cholesky and do in place

                # TODO: Later investigate using a pivoted cholesky instead (i.e. Val(true)) because otherwise matrices that are semi-definite will fail
                cholesky!(c.V.chol.factors, Val(false);
                          check = settings.check_posdef_cholesky) # inplace uses V_t with cholesky.  Now V[t]'s chol is upper-UpperTriangular

                # check scale of diagonal to see if it was explosive
                if settings.tol_cholesky > 0 && (norm(c.V.mat, Inf) > settings.tol_cholesky)
                    throw(ErrorException("Failing on norm of covariance matrix"))
                end

            catch e
                if settings.rethrow_exceptions
                    rethrow(e)
                elseif e isa PosDefException
                    (settings.print_level > 0) && display(e)
                    return :POSDEF_EXCEPTION
                else
                    settings.print_level == 0 || display(e)
                    return :GENERAL_CHOLESKY_FAIL
                end
            end
        end
        # eta * Gamma
        mul!(c.B, c.η, c.Γ)

        # @exfiltrate  # flip on to see intermediate calculations.  TURN OFF BEFORE PROFILING
    catch e
        if settings.rethrow_exceptions
            rethrow(e)
        elseif e isa LAPACKException || e isa PosDefException
            (settings.print_level > 0) && display(e)
            return :LAPACK_Error
        elseif e isa PosDefException
            (settings.print_level > 0) && display(e)
            return :POSDEF_EXCEPTION
        else
            settings.print_level == 0 || display(e)
            return :FAILURE # generic failure
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
        buff.B .= c.I_x_2
        buff.C[:, 1:n_y] .= c.H_yp
        kron!(buff.D, c.h_x, c.h_x)
        # TODO: Tullio/etc. quadratic form trickier for any of this?
        buff.R .= vcat(c.g_x * c.h_x, c.g_x, c.h_x, c.I_x)
        for i in 1:n
            buff.E[i, :] = vec(buff.R' * c.Ψ[i] * buff.R) # (24), flip the sign for gsylv
        end
        buff.E .*= -1

        # Sylvester
        # NOTE: Michel's package overwrites buff.E, while the function from MatrixEquations creates a new variable
        if settings.sylvester_solver == :GeneralizedSylvesterSolver
            ws = GeneralizedSylvesterWs(n, n, n_x, 2)
            generalized_sylvester_solver!(buff.A, buff.C, c.h_x, buff.E, 2, ws)
        else # use the old MatrixEquations
            X = gsylv(buff.A, buff.B, buff.C, buff.D, buff.E) # (22)
            buff.E .= X
        end
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
        if settings.rethrow_exceptions
            rethrow(e)
        else
            settings.print_level == 0 || display(e)
            return :Failure # generic failure
        end
    end
    return :Success
end
