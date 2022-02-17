# Fill in the gradients . No defaults since almost always called internally to custom rule
function generate_perturbation_derivatives!(m, p_d, p_f, cache::AbstractSolverCache{1};
                                            settings = PerturbationSolverSettings())
    @assert cache.p_d_symbols == collect(Symbol.(keys(p_d)))

    p = isnothing(p_f) ? p_d : order_vector_by_symbols(merge(p_d, p_f), m.mod.p_symbols)

    # Fill in derivatives
    ret = evaluate_first_order_functions_p!(m, cache, settings, p)
    (ret == :Success) || return ret
    ret = solve_first_order_p!(m, cache, settings)
    maybe_call_function(settings.solve_first_order_p_callback, ret, m, cache, settings)
    return ret
end

function generate_perturbation_derivatives!(m, p_d, p_f, cache::AbstractSolverCache{2};
                                            settings = PerturbationSolverSettings())
    @assert cache.p_d_symbols == collect(Symbol.(keys(p_d)))

    p = isnothing(p_f) ? p_d : order_vector_by_symbols(merge(p_d, p_f), m.mod.p_symbols)

    # Fill in derivatives, first by calling the first-order
    # Fill in derivatives
    ret = evaluate_first_order_functions_p!(m, cache, settings, p)
    (ret == :Success) || return ret
    ret = solve_first_order_p!(m, cache, settings)
    maybe_call_function(settings.solve_first_order_p_callback, ret, m, cache, settings)

    # 2nd Order calculations
    ret = evaluate_second_order_functions_p!(m, cache, settings, p)
    (ret == :Success) || return ret
    ret = solve_second_order_p!(m, cache, settings)
    maybe_call_function(settings.solve_second_order_p_callback, ret, m, cache, settings)
    return ret
end

function evaluate_first_order_functions_p!(m, c, settings, p)
    (settings.print_level > 2) &&
        println("Evaluating first-order function derivatives into cache")
    try
        @unpack y, x = c  # Precondition: valid (y, x) steady states
        isnothing(c.H_p) ||
            fill_array_by_symbol_dispatch(m.mod.H_p!, c.H_p, c.p_d_symbols, y, x, p) #not required if steady_state_p!
        fill_array_by_symbol_dispatch(m.mod.H_yp_p!, c.H_yp_p, c.p_d_symbols, y, x, p)
        fill_array_by_symbol_dispatch(m.mod.H_y_p!, c.H_y_p, c.p_d_symbols, y, x, p)
        fill_array_by_symbol_dispatch(m.mod.H_xp_p!, c.H_xp_p, c.p_d_symbols, y, x, p)
        fill_array_by_symbol_dispatch(m.mod.H_x_p!, c.H_x_p, c.p_d_symbols, y, x, p)
        fill_array_by_symbol_dispatch(m.mod.Γ_p!, c.Γ_p, c.p_d_symbols, p)
        isnothing(c.Ω_p) ||
            fill_array_by_symbol_dispatch(m.mod.Ω_p!, c.Ω_p, c.p_d_symbols, p)
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

function evaluate_second_order_functions_p!(m, c, settings, p)
    (settings.print_level > 2) &&
        println("Evaluating second-order function derivatives into cache")
    try
        @unpack y, x = c  # Precondition: valid (y, x) steady states
        fill_array_by_symbol_dispatch(m.mod.Ψ_p!, c.Ψ_p, c.p_d_symbols, y, x, p)
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

function solve_first_order_p!(m, c, settings)
    @unpack n_x, n_y, n_ϵ = m
    n = n_x + n_y
    n_p = length(c.p_d_symbols)
    (settings.print_level > 2) && println("Solving first order derivatives of perturbation")

    buff = c.first_order_solver_p_buffer
    try
        if isnothing(m.mod.ȳ_p!) && isnothing(m.mod.x̄_p!)
            # Zeroth-order derivatives if not provided
            # Calculating c.y_p, c.x_p
            A_zero = [c.H_y + c.H_yp c.H_x + c.H_xp]
            # H_p is a vector of vectors
            x_zeroth = A_zero \ -hcat(c.H_p...) # (47)
            c.y_p .= x_zeroth[1:n_y, :]
            c.x_p .= x_zeroth[(n_y + 1):n, :]
        end

        # The derivatives
        for i in 1:n_p
            c.Σ_p[i] .= Symmetric(c.Γ_p[i] * c.Γ' + c.Γ * c.Γ_p[i]')
        end

        # Write equation (52) as E + AX + CXD = 0, a generalized Sylvester equation
        # first-order derivatives
        buff.R .= vcat(c.g_x * c.h_x, c.g_x, c.h_x, c.I_x)
        buff.A .= [c.H_y c.H_xp + c.H_yp * c.g_x]

        # i.e. C = [c.H_yp zeros(n, n_x)]
        buff.C[:, 1:n_y] .= c.H_yp
        ws = GeneralizedSylvesterWs(n, n, n_x, 1)

        # Initialize
        for i in 1:n_p
            # p-specific Sylvester preparation
            buff.bar[1:n_y] = c.y_p[i]
            buff.bar[(n_y + 1):(2 * n_y)] = c.y_p[i]
            buff.bar[(2 * n_y + 1):(2 * n_y + n_x)] = c.x_p[i]
            buff.bar[(2 * n_y + n_x + 1):end] = c.x_p[i]

            buff.dH[:, 1:n_y] = c.H_yp_p[i]
            buff.dH[:, (n_y + 1):(2 * n_y)] = c.H_y_p[i]
            buff.dH[:, (2 * n_y + 1):(2 * n_y + n_x)] = c.H_xp_p[i]
            buff.dH[:, (2 * n_y + n_x + 1):end] = c.H_x_p[i]
            
            for j in 1:n
                buff.dH[j, :] += c.Ψ[j] * buff.bar
            end
            mul!(buff.E, buff.dH, buff.R, -1.0, 0.0)

            # solves AXB + CXD = E
            # Sylvester
            generalized_sylvester_solver!(buff.A, buff.C, c.h_x, buff.E, 1, ws)
            c.g_x_p[i] .= buff.E[1:n_y, :]
            c.h_x_p[i] .= buff.E[(n_y + 1):n, :]

            # Q weighted derivatives
            c.C_1_p[i] .= c.Q * vcat(c.g_x_p[i], zeros(n_x, n_x))
            c.A_1_p[i] .= c.h_x_p[i]

            # V derivatives
            tmp = c.h_x_p[i] * Array(c.V) * c.h_x'
            c.V_p[i] .= lyapd(c.h_x, c.η * c.Σ_p[i] * c.η' + tmp + tmp')
            # B derivatives
            c.B_p[i] .= c.η * c.Γ_p[i]
        end

        #@exfiltrate  # flip on to see intermediate calculations.  TURN OFF BEFORE PROFILING        

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

function solve_second_order_p!(m, c, settings)
    @unpack n_x, n_y, n_ϵ, n_z = m
    n = n_x + n_y
    n_p = length(c.p_d_symbols)
    (settings.print_level > 2) &&
        println("Solving second order derivatives of perturbation")
    buff = c.second_order_solver_p_buffer

    try
        # General Prep
        buff.A .= [c.H_y c.H_xp + c.H_yp * c.g_x]
        buff.C[:, 1:n_y] .= c.H_yp
        ws = GeneralizedSylvesterWs(n, n, n_x, 2)
        buff.E .= zeros(n, n_x * n_x)
        buff.R .= vcat(c.g_x * c.h_x, c.g_x, c.h_x, c.I_x)
        buff.A_σ .= [c.H_yp + c.H_y c.H_xp + c.H_yp * c.g_x]
        buff.R_σ .= vcat(c.g_x, zeros(n_y, n_x), c.I_x, zeros(n_x, n_x))

        buff.gh_stack[1:n_y, :] .= reshape(c.g_xx, n_y, n_x * n_x)
        buff.gh_stack[n_y+1:end, :] .= reshape(c.h_xx, n_x, n_x * n_x)
        buff.g_xx_flat .= reshape(c.g_xx, n_y, n_x * n_x)
        for i in 1:n
            for j in 1:n_x
                buff.Ψ_x_sum[j][i] .= c.Ψ_xp[j][i] .+ c.Ψ_x[j][i]
            end
            for j in 1:n_y
                buff.Ψ_y_sum[j][i] .= c.Ψ_yp[j][i] .+ c.Ψ_y[j][i]
            end
        end

        tmp1 = similar(buff.R')
        tmp2 = similar(c.h_x)
        for i in 1:n_p
            # Prep for _xx_p      
            # Compute the total derivatives
            buff.bar[1:n_y] = c.y_p[i]
            buff.bar[(n_y + 1):(2 * n_y)] = c.y_p[i]
            buff.bar[(2 * n_y + 1):(2 * n_y + n_x)] = c.x_p[i]
            buff.bar[(2 * n_y + n_x + 1):end] = c.x_p[i]
            buff.dH[:, 1:n_y] = c.H_yp_p[i]
            buff.dH[:, (n_y + 1):(2 * n_y)] = c.H_y_p[i]
            buff.dH[:, (2 * n_y + 1):(2 * n_y + n_x)] = c.H_xp_p[i]
            buff.dH[:, (2 * n_y + n_x + 1):end] = c.H_x_p[i]
            for j in 1:n
                buff.dH[j, :] += c.Ψ[j] * buff.bar
            end
            for j in 1:n
                buff.dΨ[j] .= c.Ψ_p[i][j]
            end
            for j in 1:n_y
                if (c.y_p[i][j] != 0)
                    for k in 1:n
                        buff.dΨ[k] .+= buff.Ψ_y_sum[j][k] .* c.y_p[i][j]
                    end
                end
            end
            for j in 1:n_x
                if (c.x_p[i][j] != 0)
                    for k in 1:n
                        buff.dΨ[k] .+= buff.Ψ_x_sum[j][k] .* c.x_p[i][j]
                    end
                end
            end

            # Constants: (60)
            buff.R_p[1:n_y, :] = c.g_x_p[i] * c.h_x + c.g_x * c.h_x_p[i]
            buff.R_p[(n_y + 1):(2 * n_y), :] = c.g_x_p[i]
            buff.R_p[(2 * n_y + 1):(2 * n_y + n_x), :] = c.h_x_p[i]
            # Flip the sign of E for Sylvester input
            fill!(buff.E, 0.0)
            for j in 1:n
                mul!(tmp1, buff.R_p', c.Ψ[j])
                mul!(tmp2, tmp1, buff.R)
                buff.E[j, :] .-= vec(tmp2)
                buff.E[j, :] .-= vec(tmp2')
                mul!(tmp1, buff.R', buff.dΨ[j])
                mul!(tmp2, tmp1, buff.R)
                buff.E[j, :] .-= vec(tmp2)
                # buff.E[j, :] .= -vec(buff.R_p' * c.Ψ[j] * buff.R + buff.R' * c.Ψ[j] * buff.R_p + buff.R' * buff.dΨ[j] * buff.R)
            end
            # Constants: (56)
            kron!(buff.kron_h_x, c.h_x, c.h_x)
            buff.E .-= buff.dH[:, 1:n_y] * buff.g_xx_flat * buff.kron_h_x # Plug (57) in (56)
            tmp = c.H_yp * buff.g_xx_flat
            kron!(buff.kron_h_x, c.h_x_p[i], c.h_x)
            mul!(buff.E, tmp, buff.kron_h_x, -1.0, 1.0) # Plug (58) in (56), step 1
            kron!(buff.kron_h_x, c.h_x, c.h_x_p[i]) 
            mul!(buff.E, tmp, buff.kron_h_x, -1.0, 1.0) # Plug (58) in (56), step 2
            buff.E .-= hcat(buff.dH[:, (n_y + 1):(2 * n_y)],
            buff.dH[:, 1:n_y] * c.g_x +
            c.H_yp * c.g_x_p[i] +
            buff.dH[:, (2 * n_y + 1):(2 * n_y + n_x)]) * buff.gh_stack # Plug (59) in (56)

            # Solve the Sylvester equations (56)
            generalized_sylvester_solver!(buff.A, buff.C, c.h_x, buff.E, 2, ws)
            copyto!(c.g_xx_p[i], buff.E[1:n_y, :]) # Reshaping into n_y * n_x * n_x
            copyto!(c.h_xx_p[i], buff.E[(n_y + 1):end, :]) # # Reshaping into n_x * n_x * n_x

            # Prep for _σσ_p
            # Solve _σσ_p
            R_σ_p = vcat(c.g_x_p[i], zeros(n_y + n_x * 2, n_x))
            η_sq_p = c.η * c.Σ_p[i] * c.η'
            C_σ = -(buff.dH[:, 1:n_y] + buff.dH[:, (n_y + 1):(2 * n_y)]) * c.g_σσ
            C_σ -= (buff.dH[:, 1:n_y] * c.g_x + c.H_yp * c.g_x_p[i] + buff.dH[:, (2 * n_y + 1):(2 * n_y + n_x)]) * c.h_σσ
            # C_σ = -hcat(buff.dH[:, 1:n_y] + buff.dH[:, (n_y + 1):(2 * n_y)],
            #             buff.dH[:, 1:n_y] * c.g_x +
            #             c.H_yp * c.g_x_p[i] +
            #             buff.dH[:, (2 * n_y + 1):(2 * n_y + n_x)]) * vcat(c.g_σσ, c.h_σσ) # Plug (65) in (64), flip the sign to solve (64)
            C_σ -= (buff.dH[:, 1:n_y] * buff.g_xx_flat + c.H_yp * buff.E[1:n_y, :]) * vec(c.η_Σ_sq)# (67), 2nd line
            C_σ -= (c.H_yp * buff.g_xx_flat) * vec(η_sq_p) # (67), 3rd line, second part
            for j in 1:n
                mul!(tmp1, R_σ_p', c.Ψ[j])
                mul!(tmp2, tmp1, buff.R_σ)
                C_σ[j] -= dot(tmp2, c.η_Σ_sq)
                C_σ[j] -= dot(tmp2', c.η_Σ_sq)
                mul!(tmp1, buff.R_σ', buff.dΨ[j])
                mul!(tmp2, tmp1, buff.R_σ)
                C_σ[j] -= dot(tmp2, c.η_Σ_sq)
                # C_σ[j] -= dot((R_σ_p' * c.Ψ[j] * buff.R_σ + buff.R_σ' * c.Ψ[j] * R_σ_p + buff.R_σ' * buff.dΨ[j] * buff.R_σ), c.η_Σ_sq) # (67), 1st line
                mul!(tmp1, buff.R_σ', c.Ψ[j])
                mul!(tmp2, tmp1, buff.R_σ)
                C_σ[j] -= dot(tmp2, η_sq_p)
                # C_σ[j] -= dot((buff.R_σ' * c.Ψ[j] * buff.R_σ), η_sq_p) # (67), 3rd line, first part
            end

            # Solve _σσ_p
            X_σ = buff.A_σ \ C_σ # solve (64)
            c.g_σσ_p[:, i] .= X_σ[1:n_y]
            c.h_σσ_p[:, i] .= X_σ[(n_y + 1):end]
            fill!(c.C_2_p[i], 0.0) # reset as we need to use `+=`
            for j in 1:n_z
                for k in 1:n_y
                    if (c.Q[j, k] != 0)
                        c.C_2_p[i][j, :, :] .+= 0.5 * c.Q[j, k] * c.g_xx_p[i][k, :, :]
                    end
                end
            end
            c.A_2_p[i] .= 0.5 * c.h_xx_p[i]
        end

        c.C_0_p .= 0.5 * c.Q * vcat(c.g_σσ_p, zeros(n_x, n_p))
        c.A_0_p .= 0.5 * c.h_σσ_p

    catch e
        settings.print_level == 0 || display(e)
        return :Failure # generic failure
    end
    return :Success
end

function ChainRulesCore.rrule(::typeof(generate_perturbation), m::PerturbationModel,
                              p_d::NamedTuple{DFieldsType,DTupleType}, p_f, order::Val{1};
                              cache = SolverCache(m, Val(1), p_d),
                              settings = PerturbationSolverSettings()) where {DFieldsType,
                                                                              DTupleType}
    (settings.print_level > 2) && println("Calculating generate_perturbation primal ")
    sol = generate_perturbation(m, p_d, p_f, Val(1); cache, settings)
    if (sol.retcode == :Success)
        grad_ret = generate_perturbation_derivatives!(m, p_d, p_f, cache)
        if (grad_ret != :Success)
            sol = FirstOrderPerturbationSolution(:GradientFailure, m, cache)
        end
    end

    function generate_perturbation_pb(Δsol)
        (settings.print_level > 2) && println("Calculating generate_perturbation pullback")
        Δp = (p_d === nothing) ? nothing : zeros(length(p_d))
        if (sol.retcode == :Success) & (p_d !== nothing)
            c = cache # temp to avoid renaming everything
            n_p_d = length(p_d)
            if (~iszero(Δsol.A))
                for i in 1:n_p_d
                    Δp[i] += dot(c.h_x_p[i], Δsol.A)
                end
            end
            if (~iszero(Δsol.g_x))
                for i in 1:n_p_d
                    Δp[i] += dot(c.g_x_p[i], Δsol.g_x)
                end
            end
            if (~iszero(Δsol.C))
                for i in 1:n_p_d
                    Δp[i] += dot(c.C_1_p[i], Δsol.C)
                end
            end
            if (~isnothing(Δsol.x_ergodic))
                if ((Δsol.x_ergodic != NoTangent()) & (Δsol.x_ergodic != ZeroTangent()))
                    for i in 1:n_p_d
                        Δp[i] += dot(c.V_p[i], Δsol.x_ergodic.Σ)
                    end
                end
            end
            if (~iszero(Δsol.Γ))
                for i in 1:n_p_d
                    Δp[i] += dot(c.Γ_p[i], Δsol.Γ)
                end
            end
            if (~iszero(Δsol.B))
                for i in 1:n_p_d
                    Δp[i] += dot(c.B_p[i], Δsol.B)
                end
            end
            if (~isnothing(Δsol.D)) # D is a Distribution object which complicates stuff here
                if ((Δsol.D != NoTangent()) & (Δsol.D != ZeroTangent()))
                    # Only supports diagonal matrices for now.
                    ΔΩ = diag(Δsol.D.Σ) .* c.Ω * 2
                    for i in 1:n_p_d
                        Δp[i] += dot(c.Ω_p[i], ΔΩ)
                    end
                end
            end
            if (~iszero(Δsol.x))
                for i in 1:n_p_d
                    Δp[i] += dot(c.x_p[i], Δsol.x)
                end
            end
            if (~iszero(Δsol.y))
                for i in 1:n_p_d
                    Δp[i] += dot(c.y_p[i], Δsol.y)
                end
            end
        end
        Δp_nt = NamedTuple{DFieldsType,DTupleType}(tuple(Δp...))  # turn tuple into named tuple in the same order
        return NoTangent(), NoTangent(),
               Tangent{NamedTuple{DFieldsType,DTupleType},
                       NamedTuple{DFieldsType,DTupleType}}(Δp_nt), NoTangent(), NoTangent()
    end
    # keep the named tuple the same
    return sol, generate_perturbation_pb
end

function ChainRulesCore.rrule(::typeof(generate_perturbation), m::PerturbationModel,
                              p_d::NamedTuple{DFieldsType,DTupleType}, p_f, order::Val{2};
                              cache = SolverCache(m, Val(2), p_d),
                              settings = PerturbationSolverSettings()) where {DFieldsType,
                                                                              DTupleType}
    (settings.print_level > 2) && println("Calculating generate_perturbation primal ")
    sol = generate_perturbation(m, p_d, p_f, Val(2); cache, settings)
    if (sol.retcode == :Success)
        grad_ret = generate_perturbation_derivatives!(m, p_d, p_f, cache)
        if (grad_ret != :Success)
            sol = SecondOrderPerturbationSolution(:GradientFailure, m, cache)
        end
    end
    function generate_perturbation_pb(Δsol)
        (settings.print_level > 2) && println("Calculating generate_perturbation pullback")
        Δp = (p_d === nothing) ? nothing : zeros(length(p_d))
        if (sol.retcode == :Success) & (p_d !== nothing)
            c = cache # temp to avoid renaming everything
            n_p_d = length(p_d)
            if (~iszero(Δsol.A_1))
                for i in 1:n_p_d
                    Δp[i] += dot(c.A_1_p[i], Δsol.A_1)
                end
            end
            if (~iszero(Δsol.g_x))
                for i in 1:n_p_d
                    Δp[i] += dot(c.g_x_p[i], Δsol.g_x)
                end
            end
            if (~iszero(Δsol.C_1))
                for i in 1:n_p_d
                    Δp[i] += dot(c.C_1_p[i], Δsol.C_1)
                end
            end
            if (~iszero(Δsol.Γ))
                for i in 1:n_p_d
                    Δp[i] += dot(c.Γ_p[i], Δsol.Γ)
                end
            end
            if (~iszero(Δsol.B))
                for i in 1:n_p_d
                    Δp[i] += dot(c.B_p[i], Δsol.B)
                end
            end
            if (~iszero(Δsol.A_2))
                for i in 1:n_p_d
                    Δp[i] += dot(c.A_2_p[i], Δsol.A_2)
                end
            end
            if (~iszero(Δsol.g_xx))
                for i in 1:n_p_d
                    Δp[i] += dot(c.g_xx_p[i], Δsol.g_xx)
                end
            end
            if (~iszero(Δsol.C_2))
                for i in 1:n_p_d
                    Δp[i] += dot(c.C_2_p[i], Δsol.C_2)
                end
            end
            if (~isnothing(Δsol.D)) # D is a Distribution object which complicates stuff here
                if ((Δsol.D != NoTangent()) & (Δsol.D != ZeroTangent()))
                    # Only supports diagonal matrices for now.
                    ΔΩ = diag(Δsol.D.Σ) .* c.Ω * 2
                    for i in 1:n_p_d
                        Δp[i] += dot(c.Ω_p[i], ΔΩ)
                    end
                end
            end
            if (~iszero(Δsol.x))
                for i in 1:n_p_d
                    Δp[i] += dot(c.x_p[i], Δsol.x)
                end
            end
            if (~iszero(Δsol.y))
                for i in 1:n_p_d
                    Δp[i] += dot(c.y_p[i], Δsol.y)
                end
            end
            # Currently the results for the second-order correction terms are not vector of vectors, but just arrays
            if (~iszero(Δsol.C_0))
                Δp += c.C_0_p' * Δsol.C_0
            end
            if (~iszero(Δsol.g_σσ))
                Δp += c.g_σσ_p' * Δsol.g_σσ
            end
            if (~iszero(Δsol.A_0))
                Δp += c.A_0_p' * Δsol.A_0
            end
        end

        Δp_nt = NamedTuple{DFieldsType,DTupleType}(tuple(Δp...)) # turn tuple into named tuple in the same order
        return NoTangent(), NoTangent(),
               Tangent{NamedTuple{DFieldsType,DTupleType},
                       NamedTuple{DFieldsType,DTupleType}}(Δp_nt), NoTangent(), NoTangent()
    end

    return sol, generate_perturbation_pb
end
