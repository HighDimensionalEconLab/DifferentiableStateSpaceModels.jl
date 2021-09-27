# Fill in the gradients . No defaults since almost always called internally to custom rule
function generate_perturbation_derivatives!(m, p_d, p_f, cache::AbstractSolverCache{1};
                                               settings=PerturbationSolverSettings())
    @assert cache.p_d_symbols == collect(Symbol.(keys(p_d)))

    p = isnothing(p_f) ? p_d : order_vector_by_symbols(merge(p_d, p_f), m.mod.p_symbols)

    # Fill in derivatives
    ret = evaluate_first_order_functions_p!(m, cache, settings, p)
    (ret == :Success) || return ret
    ret = solve_first_order_p!(m, cache, settings)
    maybe_call_function(settings.solve_first_order_p_callback, ret, m, cache, settings)
    return (ret == :Success) || return ret
end

function generate_perturbation_derivatives!(m, p_d, p_f, cache::AbstractSolverCache{2};
                                                settings=PerturbationSolverSettings())
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
    return (ret == :Success) || return ret
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
        R = vcat(c.g_x * c.h_x, c.g_x, c.h_x, I(n_x))
        A = [c.H_y c.H_xp + c.H_yp * c.g_x]
        B = Array(I(n_x) * 1.0)
        C = [c.H_yp zeros(n, n_x)]
        D = c.h_x
        AS, CS, Q1, Z1 = schur(A, C)
        BS, DS, Q2, Z2 = schur(B, D) # careful going inplace if passing in c.h_x.  B is a buffer?
        # Initialize
        dH = zeros(2n, n)
        bar = zeros(2n, 1)
        Hstack = zeros(n, 2n)
        for i in 1:n_p
            # p-specific Sylvester preparation
            bar[1:n_y] = c.y_p[i]
            bar[(n_y + 1):(2 * n_y)] = c.y_p[i]
            bar[(2 * n_y + 1):(2 * n_y + n_x)] = c.x_p[i]
            bar[(2 * n_y + n_x + 1):end] = c.x_p[i]

            Hstack[:, 1:n_y] = c.H_yp_p[i]
            Hstack[:, (n_y + 1):(2 * n_y)] = c.H_y_p[i]
            Hstack[:, (2 * n_y + 1):(2 * n_y + n_x)] = c.H_xp_p[i]
            Hstack[:, (2 * n_y + n_x + 1):end] = c.H_x_p[i]
            for j in 1:n
                dH[:, j] = c.Ψ[j] * bar + Hstack[j, :]
            end
            E = -dH'R

            # solves AXB + CXD = E
            # sylvester
            # X = gsylv(A, B, C, D, E)
            Y = adjoint(Q1) * (E * Z2)
            gsylvs!(AS, BS, CS, DS, Y)
            X = Z1 * (Y * adjoint(Q2))
            c.g_x_p[i] .= X[1:n_y, :]
            c.h_x_p[i] .= X[(n_y + 1):n, :]

            # Q weighted derivatives
            c.C_1_p[i] .= c.Q * vcat(c.g_x_p[i], zeros(n_x, n_x))
            c.A_1_p[i] .= c.h_x_p[i]

            # V derivatives
            tmp = c.h_x_p[i] * Array(c.V) * c.h_x'
            c.V_p[i] .= lyapd(c.h_x, c.η * c.Σ_p[i] * c.η' + tmp + tmp')
            # B derivatives
            c.B_p[i] .= c.η * c.Γ_p[i]
        end

        @exfiltrate  # flip on to see intermediate calculations.  TURN OFF BEFORE PROFILING        

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


    # General Prep
    A = [c.H_y c.H_xp + c.H_yp * c.g_x]
    B = Array(I(n_x * n_x) * 1.0)
    C = [c.H_yp zeros(n, n_x)]
    D = kron(c.h_x, c.h_x)
    AS, CS, Q1, Z1 = schur(A, C)
    BS, DS, Q2, Z2 = schur(B, D)
    E = zeros(n, n_x * n_x)
    R = vcat(c.g_x * c.h_x, c.g_x, c.h_x, I(n_x))
    A_σ = [c.H_yp + c.H_y c.H_xp + c.H_yp * c.g_x]
    R_σ = vcat(c.g_x, zeros(n_y, n_x), I(n_x), zeros(n_x, n_x))
    gh_stack = vcat(reshape(c.g_xx, n_y, n_x * n_x), reshape(c.h_xx, n_x, n_x * n_x))
    g_xx_flat = reshape(c.g_xx, n_y, n_x * n_x)

    dH = zeros(n, 2n)
    dΨ = [zeros(2n, 2n) for _ in 1:n]

    for i in 1:n_p
        # Prep for _xx_p
        # Compute the total derivatives
        bar = vcat(c.y_p[i], c.y_p[i], c.x_p[i], c.x_p[i])
        Hstack = hcat(c.H_yp_p[i], c.H_y_p[i], c.H_xp_p[i], c.H_x_p[i])
        for j in 1:n
            dH[j, :] = c.Ψ[j] * bar + Hstack[j, :]
            dΨ[j] .= c.Ψ_p[i][j]
            for k in 1:n_y
                if (c.y_p[i][k] != 0)
                    dΨ[j] += (c.Ψ_yp[k][j] + c.Ψ_y[k][j]) * c.y_p[i][k]
                end
            end
            for k in 1:n_x
                if (c.x_p[i][k] != 0)
                    dΨ[j] += (c.Ψ_xp[k][j] + c.Ψ_x[k][j]) * c.x_p[i][k]
                end
            end
        end

        # Constants: (60)
        R_p = vcat(c.g_x_p[i] * c.h_x + c.g_x * c.h_x_p[i], c.g_x_p[i], c.h_x_p[i],
                   zeros(n_x, n_x))
        # Flip the sign of E for Sylvester input
        for j in 1:n
            E[j, :] .= -((R_p' * c.Ψ[j] * R)[:] +
                         (R' * c.Ψ[j] * R_p)[:] +
                         (R' * dΨ[j] * R)[:])
        end
        # Constants: (56)
        E -= hcat(dH[:, 1:n_y], zeros(n, n_x)) * gh_stack * kron(c.h_x, c.h_x) # Plug (57) in (56)
        E -= hcat(c.H_yp, zeros(n, n_x)) *
             gh_stack *
             (kron(c.h_x_p[i], c.h_x) + kron(c.h_x, c.h_x_p[i])) # Plug (58) in (56)
        E -= hcat(dH[:, (n_y + 1):(2 * n_y)],
                  dH[:, 1:n_y] * c.g_x +
                  c.H_yp * c.g_x_p[i] +
                  dH[:, (2 * n_y + 1):(2 * n_y + n_x)]) * gh_stack # Plug (59) in (56)

        # Solve the Sylvester equations (56)
        # X = gsylv(A, B, C, D, E)
        Y = adjoint(Q1) * (E * Z2)
        gsylvs!(AS, BS, CS, DS, Y)
        X = Z1 * (Y * adjoint(Q2))
        c.g_xx_p[i] .= reshape(X[1:n_y, :], n_y, n_x, n_x)
        c.h_xx_p[i] .= reshape(X[(n_y + 1):end, :], n_x, n_x, n_x)

        # Prep for _σσ_p
        # Solve _σσ_p
        R_σ_p = vcat(c.g_x_p[i], zeros(n_y + n_x * 2, n_x))
        η_sq_p = c.η * c.Σ_p[i] * c.η'
        C_σ = -hcat(dH[:, 1:n_y] + dH[:, (n_y + 1):(2 * n_y)],
                    dH[:, 1:n_y] * c.g_x +
                    c.H_yp * c.g_x_p[i] +
                    dH[:, (2 * n_y + 1):(2 * n_y + n_x)]) * vcat(c.g_σσ, c.h_σσ) # Plug (65) in (64), flip the sign to solve (64)
        C_σ -= (dH[:, 1:n_y] * g_xx_flat + c.H_yp * X[1:n_y, :]) * c.η_Σ_sq[:] # (67), 2nd line
        C_σ -= (c.H_yp * g_xx_flat) * η_sq_p[:] # (67), 3rd line, second part
        for j in 1:n
            C_σ[j] -= dot((R_σ_p' * c.Ψ[j] * R_σ +
                           R_σ' * c.Ψ[j] * R_σ_p +
                           R_σ' * dΨ[j] * R_σ), c.η_Σ_sq) # (67), 1st line
            C_σ[j] -= dot((R_σ' * c.Ψ[j] * R_σ), η_sq_p) # (67), 3rd line, first part
        end

        # Solve _σσ_p
        X_σ = A_σ \ C_σ # solve (64)
        c.g_σσ_p[:, i] .= X_σ[1:n_y]
        c.h_σσ_p[:, i] .= X_σ[(n_y + 1):end]
        fill!(c.C_2_p[i], 0.0) # reset as we need to use `+=`
        for j in 1:n_z
            for k in 1:n_y
                c.C_2_p[i][j, :, :] += 0.5 * c.Q[j, k] * c.g_xx_p[i][k, :, :]
            end
        end
        c.A_2_p[i] .= 0.5 * c.h_xx_p[i]
    end

    c.C_0_p .= 0.5 * c.Q * vcat(c.g_σσ_p, zeros(n_x, n_p))
    c.A_0_p .= 0.5 * c.h_σσ_p
    return :Success
end

# TODO: Hook up adjoints
# CHANGES:
# = Args separate p_d and p_f where p = p_d + p_f solve_first_order
# - Use n_p_d rather than n_p to access inside derivatives.
# - call the solve and then the derivatives separately
# - y_p, x_p, Omega_p all now vectors of vectors.
# - You use the size of p

# function ChainRulesCore.rrule(
#     ::typeof(generate_perturbation),
#     m::AbstractFirstOrderPerturbationModel,
#     p;
#     p_f = nothing,
#     cache = allocate_cache(m),
#     settings = PerturbationSolverSettings(),
# )
#     (settings.print_level > 2) && println("Calculating generate_perturbation primal ")
#     sol = generate_perturbation(m, p; p_f, cache, settings)
#     c = cache # temp to avoid renaming everything

#     function generate_perturbation_pb(Δsol)
#         (settings.print_level > 2) && println("Calculating generate_perturbation pullback")
#         Δp = (p === nothing) ? nothing : zeros(length(p))
#         if (sol.retcode == :Success) & (p !== nothing)
#             if (~iszero(Δsol.A))
#                 for i = 1:sol.n_p
#                     Δp[i] += dot(c.h_x_p[i], Δsol.A)
#                 end
#             end
#             if (~iszero(Δsol.g_x))
#                 for i = 1:sol.n_p
#                     Δp[i] += dot(c.g_x_p[i], Δsol.g_x)
#                 end
#             end
#             if (~iszero(Δsol.C))
#                 for i = 1:sol.n_p
#                     Δp[i] += dot(c.C_1_p[i], Δsol.C)
#                 end
#             end
#             if (~iszero(Δsol.x_ergodic))
#                 for i = 1:sol.n_p
#                     tmp = c.V.L \ c.V_p[i] / c.V.U
#                     tmp[diagind(tmp)] /= 2.0
#                     t1 = c.V.L * LowerTriangular(tmp)
#                     # Cholesky by default stores the U part, but that information is lost when passing back
#                     # from MvNormal. Therefore, we have to extract the components out manually
#                     Δp[i] += dot(t1', UpperTriangular(Δsol.x_ergodic.C.factors)) # dot(c.V_p[i], Δsol.x_ergodic) is the original logic
#                 end
#             end
#             if (~iszero(Δsol.Γ))
#                 for i = 1:sol.n_p
#                     Δp[i] += dot(c.Γ_p[i], Δsol.Γ)
#                 end
#             end
#             if (~iszero(Δsol.B))
#                 for i = 1:sol.n_p
#                     Δp[i] += dot(c.B_p[i], Δsol.B)
#                 end
#             end
#             if (~iszero(Δsol.D))
#                 # Only supports diagonal matrices for now.
#                 Δp += c.Ω_p' * (Δsol.D.σ)
#             end
#             if (~iszero(Δsol.x))
#                 Δp += c.x_p' * Δsol.x
#             end
#             if (~iszero(Δsol.y))
#                 Δp += c.y_p' * Δsol.y
#             end
#         end
#         return nothing, nothing, Δp
#     end
#     # keep the named tuple the same
#     return sol, generate_perturbation_pb
# end

# function ChainRulesCore.rrule(
#     ::typeof(generate_perturbation),
#     m::AbstractSecondOrderPerturbationModel,
#     p;
#     p_f = nothing,
#     cache = SecondOrderSolverCache(m),
#     settings = PerturbationSolverSettings(),
# )
#     (settings.print_level > 2) && println("Calculating generate_perturbation primal")
#     c = cache # temp to avoid renaming everything

#     sol = generate_perturbation(m, p; p_f, cache, settings)

#     function generate_perturbation_pb(Δsol)
#         (settings.print_level > 2) && println("Calculating generate_perturbation pullback")

#         Δp = (p === nothing) ? nothing : zeros(length(p))
#         if (sol.retcode == :Success) & (p !== nothing)
#             if (~iszero(Δsol.A_1))
#                 for i = 1:sol.n_p
#                     Δp[i] += dot(c.A_1_p[i], Δsol.A_1)
#                 end
#             end
#             if (~iszero(Δsol.g_x))
#                 for i = 1:sol.n_p
#                     Δp[i] += dot(c.g_x_p[i], Δsol.g_x)
#                 end
#             end
#             if (~iszero(Δsol.C_1))
#                 for i = 1:sol.n_p
#                     Δp[i] += dot(c.C_1_p[i], Δsol.C_1)
#                 end
#             end
#             if (~iszero(Δsol.Γ))
#                 for i = 1:sol.n_p
#                     Δp[i] += dot(c.Γ_p[i], Δsol.Γ)
#                 end
#             end
#             if (~iszero(Δsol.B))
#                 for i = 1:sol.n_p
#                     Δp[i] += dot(c.B_p[i], Δsol.B)
#                 end
#             end
#             if (~iszero(Δsol.A_2))
#                 for i = 1:sol.n_p
#                     Δp[i] += dot(c.A_2_p[i], Δsol.A_2)
#                 end
#             end
#             if (~iszero(Δsol.g_xx))
#                 for i = 1:sol.n_p
#                     Δp[i] += dot(c.g_xx_p[i], Δsol.g_xx)
#                 end
#             end
#             if (~iszero(Δsol.C_2))
#                 for i = 1:sol.n_p
#                     Δp[i] += dot(c.C_2_p[i], Δsol.C_2)
#                 end
#             end
#             if (~iszero(Δsol.D))
#                 Δp += c.Ω_p' * (Δsol.D.σ)
#             end
#             if (~iszero(Δsol.x))
#                 Δp += c.x_p' * Δsol.x
#             end
#             if (~iszero(Δsol.y))
#                 Δp += c.y_p' * Δsol.y
#             end
#             if (~iszero(Δsol.C_0))
#                 Δp += c.C_0_p' * Δsol.C_0
#             end
#             if (~iszero(Δsol.g_σσ))
#                 Δp += c.g_σσ_p' * Δsol.g_σσ
#             end
#             if (~iszero(Δsol.A_0))
#                 Δp += c.A_0_p' * Δsol.A_0
#             end
#         end
#         return nothing, nothing, Δp
#     end

#     # println("Generated perturbation gradient function")
#     return sol, generate_perturbation_pb
# end