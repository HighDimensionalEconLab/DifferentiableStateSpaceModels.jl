
"""
Solve a model, optionally reusing a cache (which is assumed threadsafe)

$(SIGNATURES)

# Details

Returns the solution
"""
function generate_perturbation(
    m::AbstractFirstOrderPerturbationModel,
    p_d, p_f = nothing; cache = allocate_cache(m),
    settings = PerturbationSolverSettings(),
)

    # solver type provided to all callbacks
    solver = PerturbationSolver(m, cache, settings)

    @timeit_debug "calculate_steady_state" begin
        ret = calculate_steady_state!(m, cache, settings, p, solver)
    end
    maybe_call_function(
        settings.calculate_steady_state_callback,
        ret,
        m,
        cache,
        settings,
        p,
        p_f,
        solver,
    )  # before returning
    (ret == :Success) || return FirstOrderPerturbationSolution(ret, m, cache)

    @timeit_debug "evaluate_functions" begin
        ret = evaluate_functions!(m, cache, settings, p, solver)
    end
    maybe_call_function(
        settings.evaluate_functions_callback,
        ret,
        m,
        cache,
        settings,
        p,
        p_f,
        solver,
    )
    (ret == :Success) || return FirstOrderPerturbationSolution(ret, m, cache)

    @timeit_debug "solve_first_order" begin
        ret = solve_first_order!(m, cache, settings)
    end
    maybe_call_function(settings.solve_first_order_callback, ret, m, cache, settings)
    (ret == :Success) || return FirstOrderPerturbationSolution(ret, m, cache)

    @timeit_debug "solve_first_order_p" begin
        ret = solve_first_order_p!(m, cache, settings)
    end
    maybe_call_function(settings.solve_first_order_p_callback, ret, m, cache, settings)
    (ret == :Success) || return FirstOrderPerturbationSolution(ret, m, cache)
    
    return FirstOrderPerturbationSolution(:Success, m, cache)
end

function ChainRulesCore.rrule(
    ::typeof(generate_perturbation),
    m::AbstractFirstOrderPerturbationModel,
    p;
    p_f = nothing,
    cache = allocate_cache(m),
    settings = PerturbationSolverSettings(),
)
    (settings.print_level > 2) && println("Calculating generate_perturbation primal ")
    sol = generate_perturbation(m, p; p_f, cache, settings)
    c = cache # temp to avoid renaming everything

    function generate_perturbation_pb(Δsol)
        (settings.print_level > 2) && println("Calculating generate_perturbation pullback")
        Δp = (p === nothing) ? nothing : zeros(length(p))
        if (sol.retcode == :Success) & (p !== nothing)
            if (~iszero(Δsol.A))
                for i = 1:sol.n_p
                    Δp[i] += dot(c.h_x_p[i], Δsol.A)
                end
            end
            if (~iszero(Δsol.g_x))
                for i = 1:sol.n_p
                    Δp[i] += dot(c.g_x_p[i], Δsol.g_x)
                end
            end
            if (~iszero(Δsol.C))
                for i = 1:sol.n_p
                    Δp[i] += dot(c.C_1_p[i], Δsol.C)
                end
            end
            if (~iszero(Δsol.x_ergodic))
                for i = 1:sol.n_p
                    tmp = c.V.L \ c.V_p[i] / c.V.U
                    tmp[diagind(tmp)] /= 2.0
                    t1 = c.V.L * LowerTriangular(tmp)
                    # Cholesky by default stores the U part, but that information is lost when passing back
                    # from MvNormal. Therefore, we have to extract the components out manually
                    Δp[i] += dot(t1', UpperTriangular(Δsol.x_ergodic.C.factors)) # dot(c.V_p[i], Δsol.x_ergodic) is the original logic
                end
            end
            if (~iszero(Δsol.Γ))
                for i = 1:sol.n_p
                    Δp[i] += dot(c.Γ_p[i], Δsol.Γ)
                end
            end
            if (~iszero(Δsol.B))
                for i = 1:sol.n_p
                    Δp[i] += dot(c.B_p[i], Δsol.B)
                end
            end
            if (~iszero(Δsol.D))
                # Only supports diagonal matrices for now.
                Δp += c.Ω_p' * (Δsol.D.σ)
            end
            if (~iszero(Δsol.x))
                Δp += c.x_p' * Δsol.x
            end
            if (~iszero(Δsol.y))
                Δp += c.y_p' * Δsol.y
            end
        end
        return nothing, nothing, Δp
    end
    # keep the named tuple the same
    return sol, generate_perturbation_pb
end

function generate_perturbation(
    m::AbstractSecondOrderPerturbationModel,
    p;
    p_f = nothing,
    cache = SecondOrderSolverCache(m),
    settings = PerturbationSolverSettings(),
)

    @unpack use_solution_cache = settings
    # solver type provided to all callbacks
    solver = PerturbationSolver(m, cache, settings)

    ret = calculate_steady_state!(m, cache, settings, p, solver)
    maybe_call_function(
        settings.calculate_steady_state_callback,
        ret,
        m,
        cache,
        settings,
        p,
        p_f,
        solver,
    )  # before returning
    (ret == :Success) || return SecondOrderPerturbationSolution(ret, m, cache)

    # need to evaluate even if perturbations are not required.  Later could separate more cleanly
    ret = evaluate_functions!(m, cache, settings, p, solver)
    maybe_call_function(
        settings.evaluate_functions_callback,
        ret,
        m,
        cache,
        settings,
        p,
        p_f,
        solver,
    )
    (ret == :Success) || return SecondOrderPerturbationSolution(ret, m, cache)

    ret = solve_first_order!(m, cache, settings)
    maybe_call_function(settings.solve_first_order_callback, ret, m, cache, settings)
    (ret == :Success) || return SecondOrderPerturbationSolution(ret, m, cache)

    ret = solve_first_order_p!(m, cache, settings)
    maybe_call_function(settings.solve_first_order_p_callback, ret, m, cache, settings)
    (ret == :Success) || return SecondOrderPerturbationSolution(ret, m, cache)

    ret = solve_second_order!(m, cache, settings)
    maybe_call_function(settings.solve_second_order_callback, ret, m, cache, settings)
    (ret == :Success) || return SecondOrderPerturbationSolution(ret, m, cache)

    ret = solve_second_order_p!(m, cache, settings)
    maybe_call_function(settings.solve_second_order_p_callback, ret, m, cache, settings)
    (ret == :Success) || return SecondOrderPerturbationSolution(ret, m, cache)
    
    return SecondOrderPerturbationSolution(:Success, m, cache)
end

function ChainRulesCore.rrule(
    ::typeof(generate_perturbation),
    m::AbstractSecondOrderPerturbationModel,
    p;
    p_f = nothing,
    cache = SecondOrderSolverCache(m),
    settings = PerturbationSolverSettings(),
)
    (settings.print_level > 2) && println("Calculating generate_perturbation primal")
    c = cache # temp to avoid renaming everything

    sol = generate_perturbation(m, p; p_f, cache, settings)

    function generate_perturbation_pb(Δsol)
        (settings.print_level > 2) && println("Calculating generate_perturbation pullback")

        Δp = (p === nothing) ? nothing : zeros(length(p))
        if (sol.retcode == :Success) & (p !== nothing)
            if (~iszero(Δsol.A_1))
                for i = 1:sol.n_p
                    Δp[i] += dot(c.A_1_p[i], Δsol.A_1)
                end
            end
            if (~iszero(Δsol.g_x))
                for i = 1:sol.n_p
                    Δp[i] += dot(c.g_x_p[i], Δsol.g_x)
                end
            end
            if (~iszero(Δsol.C_1))
                for i = 1:sol.n_p
                    Δp[i] += dot(c.C_1_p[i], Δsol.C_1)
                end
            end
            if (~iszero(Δsol.Γ))
                for i = 1:sol.n_p
                    Δp[i] += dot(c.Γ_p[i], Δsol.Γ)
                end
            end
            if (~iszero(Δsol.B))
                for i = 1:sol.n_p
                    Δp[i] += dot(c.B_p[i], Δsol.B)
                end
            end
            if (~iszero(Δsol.A_2))
                for i = 1:sol.n_p
                    Δp[i] += dot(c.A_2_p[i], Δsol.A_2)
                end
            end
            if (~iszero(Δsol.g_xx))
                for i = 1:sol.n_p
                    Δp[i] += dot(c.g_xx_p[i], Δsol.g_xx)
                end
            end
            if (~iszero(Δsol.C_2))
                for i = 1:sol.n_p
                    Δp[i] += dot(c.C_2_p[i], Δsol.C_2)
                end
            end
            if (~iszero(Δsol.D))
                Δp += c.Ω_p' * (Δsol.D.σ)
            end
            if (~iszero(Δsol.x))
                Δp += c.x_p' * Δsol.x
            end
            if (~iszero(Δsol.y))
                Δp += c.y_p' * Δsol.y
            end
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
        return nothing, nothing, Δp
    end

    # println("Generated perturbation gradient function")
    return sol, generate_perturbation_pb
end

# function calculate_steady_state!(m::AbstractPerturbationModel, c, settings, p, solver)
#     @unpack n_y, n_x, n = m

#     (settings.print_level > 2) && println("Calculating steady state")
#     try
#         if !isnothing(m.mod.ȳ!) && !isnothing(m.mod.x̄!) # use closed form if possible
#             m.mod.ȳ!(c.y, p, solver)
#             m.mod.x̄!(c.x, p, solver)
#             isnothing(m.mod.ȳ_p!) || m.mod.ȳ_p!(c.y_p, p, solver)  # p may be empty
#             isnothing(m.mod.x̄_p!) || m.mod.x̄_p!(c.x_p, p, solver)
#         elseif !isnothing(m.mod.steady_state!) # use user-provided calculation otherwise
#             m.mod.steady_state!(c.y, c.x, p, solver)
#         else # fallback is to solve system of equations from user-provided initial condition
#             y_0 = zeros(n_y)
#             x_0 = zeros(n_x)
#             m.mod.ȳ_iv!(y_0, p, solver)
#             m.mod.x̄_iv!(x_0, p, solver)
#             w_0 = [y_0; x_0]

#             if isnothing(m.mod.H̄_w!) # no jacobian
#                 nlsol = nlsolve(
#                     (H, w) -> m.mod.H̄!(H, w, p, solver),
#                     w_0;
#                     DifferentiableStateSpaceModels.nlsolve_options(solver.settings)...,
#                 )
#             else
#                 J_0 = zeros(n, n)
#                 F_0 = zeros(n)
#                 df = OnceDifferentiable(
#                     (H, w) -> m.mod.H̄!(H, w, p, solver),
#                     (J, w) -> m.mod.H̄_w!(J, w, p, solver),
#                     w_0,
#                     F_0,
#                     J_0,
#                 )  # TODO: the buffer to use for the w_0 is unclear to me same as iv?
#                 nlsol = nlsolve(
#                     df,
#                     w_0;
#                     DifferentiableStateSpaceModels.nlsolve_options(solver.settings)...,
#                 )
#             end
#             if !converged(nlsol)
#                 if settings.print_level > 0
#                     println("No steady state found\n")
#                 end
#                 return :SteadyStateFailure
#             end
#             settings.print_level > 1 &&
#                 println("Steady state found in $(nlsol.iterations) iterations\n")
#             c.y .= nlsol.zero[1:n_y]
#             c.x .= nlsol.zero[(n_y+1):end]
#         end
#     catch e
#         if !is_linear_algebra_exception(e)
#             (settings.print_level > 2) && println("Rethrowing exception")
#             rethrow(e)
#         else
#             settings.print_level == 0 || display(e)
#             return :Failure # generic failure
#         end
#     end
#     return :Success
# end

## Core algorithms

# # Requires valid preallocated solution and filled cache of steady state/etc
# # called by both the first and 2nd order models
# function solve_first_order!(m::AbstractPerturbationModel, c, settings)
#     @unpack ϵ_BK, print_level = settings
#     @unpack n_x, n_y, n_p, n_ϵ, n = m

#     (settings.print_level > 2) && println("Solving first order perturbation")
#     try
#         @timeit_debug "schur" begin
#             A = [c.H_xp c.H_yp]
#             B = [c.H_x c.H_y]
#             s = schur(complex(A), complex(B)) # Generalized Schur decomposition
#         end
#         # The generalized eigenvalues λ_i are S_ii / T_ii
#         # Following Blanchard-Kahn condition, we reorder the Schur so that
#         # S_22 ./ T_22 < 1, ie, the eigenvalues < 1 come last
#         # inds = [s.α[i] / s.β[i] >= 1 for i in 1:n]
#         @timeit_debug "ordschur" begin
#             inds = abs.(s.α) .>= (1 - ϵ_BK) .* abs.(s.β)
#             if sum(inds) != n_x
#                 # More debugging code???
#                 if print_level > 0
#                     @show n_x
#                     @show sum(inds)
#                     @show inds
#                     @show abs.(s.α)
#                     @show abs.(s.β)
#                     println("Blanchard-Kahn condition not satisfied\n")
#                 end
#                 return :BlanchardKahnFailure
#             end

#             ordschur!(s, inds)
#         end
#         # In Julia A = QSZ' and B = QTZ'

#         @timeit_debug "Extracting g_x and h_x" begin

#             @unpack S, T = s # Extract the Schur components
#             Z = s.Z'

#             b = 1:n_x
#             l = (n_x+1):n
#             g_x = -Z[l, l] \ Z[l, b]
#             blob = Z[b, b] .+ Z[b, l] * g_x
#             h_x = -blob \ (S[b, b] \ (T[b, b] * blob))
#             c.g_x .= real(g_x)
#             c.h_x .= real(h_x)
#         end

#         # fill in Σ, Ω.
#         c.Σ .= Symmetric(c.Γ * c.Γ')

#         # And derivatives
#         if (n_p > 0)
#             for i = 1:n_p
#                 c.Σ_p[i] .= Symmetric(c.Γ_p[i] * c.Γ' + c.Γ * c.Γ_p[i]')
#             end
#         end

#         # Q transforms
#         c.C_1 .= c.Q * vcat(c.g_x, diagm(ones(n_x)))

#         # Stationary Distribution
#         c.V = cholesky(Symmetric(lyapd(c.h_x, c.η * c.Σ * c.η')))
#         # eta * Gamma
#         c.B .= c.η * c.Γ
#     catch e
#         if !is_linear_algebra_exception(e)
#             (settings.print_level > 2) && println("Rethrowing exception")
#             rethrow(e)
#         else
#             settings.print_level == 0 || display(e)
#             return :Failure # generic failure
#         end
#     end
#     return :Success
# end

# # Calculate derivatives, requires `solve_first_order!` completion.
# function solve_first_order_p!(m::AbstractPerturbationModel, c, settings)
#     @unpack n_x, n_y, n_p, n_ϵ, n = m
#     (settings.print_level > 2) && println("Solving first order derivatives of perturbation")

#     if n_p == 0
#         return :Success
#     end

#     try
#         if isnothing(m.mod.ȳ_p!) && isnothing(m.mod.x̄_p!)
#             # Zeroth-order derivatives if not provided
#             @timeit_debug "Calculating c.y_p, c.x_p" begin
#                 A_zero = [c.H_y + c.H_yp c.H_x + c.H_xp]
#                 # TODO:  H_p now a vector of vectors
#                 x_zeroth = A_zero \ -c.H_p # (47)
#                 c.y_p .= x_zeroth[1:n_y, :]
#                 c.x_p .= x_zeroth[(n_y+1):n, :]
#             end
#         end

#         # Write equation (52) as E + AX + CXD = 0, a generalized Sylvester equation
#         # first-order derivatives
#         @timeit_debug "General Sylvester preparation" begin
#             R = vcat(c.g_x * c.h_x, c.g_x, c.h_x, I(n_x))
#             A = [c.H_y c.H_xp + c.H_yp * c.g_x]
#             B = Array(I(n_x) * 1.0)
#             C = [c.H_yp zeros(n, n_x)]
#             D = c.h_x
#             AS, CS, Q1, Z1 = schur(A, C)
#             BS, DS, Q2, Z2 = schur(B, D)
#             # Initialize
#             dH = zeros(2n, n)
#             bar = zeros(2n, 1)
#             Hstack = zeros(n, 2n)
#         end
#         for i = 1:n_p
#             @timeit_debug "p-specific Sylvester preparation" begin
#                 bar[1:n_y] = c.y_p[:, i]
#                 bar[(n_y+1):(2*n_y)] = c.y_p[:, i]
#                 bar[(2*n_y+1):(2*n_y+n_x)] = c.x_p[:, i]
#                 bar[(2*n_y+n_x+1):end] = c.x_p[:, i]

#                 Hstack[:, 1:n_y] = c.H_yp_p[i]
#                 Hstack[:, (n_y+1):(2*n_y)] = c.H_y_p[i]
#                 Hstack[:, (2*n_y+1):(2*n_y+n_x)] = c.H_xp_p[i]
#                 Hstack[:, (2*n_y+n_x+1):end] = c.H_x_p[i]
#                 for j = 1:n
#                     dH[:, j] = c.Ψ[j] * bar + Hstack[j, :]
#                 end
#                 E = -dH'R
#             end
#             # solves AXB + CXD = E
#             @timeit_debug "sylvester" begin
#                 # X = gsylv(A, B, C, D, E)
#                 Y = adjoint(Q1) * (E * Z2)
#                 gsylvs!(AS, BS, CS, DS, Y)
#                 X = Z1 * (Y * adjoint(Q2))
#                 c.g_x_p[i] .= X[1:n_y, :]
#                 c.h_x_p[i] .= X[(n_y+1):n, :]
#             end

#             # Q weighted derivatives
#             c.C_1_p[i] .= c.Q * vcat(c.g_x_p[i], zeros(n_x, n_x))
#             c.A_1_p[i] .= c.h_x_p[i]

#             # V derivatives
#             tmp = c.h_x_p[i] * Array(c.V) * c.h_x'
#             c.V_p[i] .= lyapd(c.h_x, c.η * c.Σ_p[i] * c.η' + tmp + tmp')
#             # B derivatives
#             c.B_p[i] .= c.η * c.Γ_p[i]
#         end
#     catch e
#         if !is_linear_algebra_exception(e)
#             (settings.print_level > 2) && println("Rethrowing exception")
#             rethrow(e)
#         else
#             settings.print_level == 0 || display(e)
#             return :Failure # generic failure
#         end
#     end
#     return :Success
# end

# #additional calculations for the 2nd order
# function solve_second_order!(
#     m::AbstractSecondOrderPerturbationModel,
#     c,
#     settings,
# )
#     @unpack n_x, n_y, n_p, n_ϵ, n_z, n, η = m

#     @timeit_debug "Sylvester prep for _xx" begin
#         A = [c.H_y c.H_xp + c.H_yp * c.g_x]
#         B = I(n_x * n_x)
#         C = [c.H_yp zeros(n, n_x)]
#         D = kron(c.h_x, c.h_x)
#         E = zeros(n, n_x * n_x)
#         R = vcat(c.g_x * c.h_x, c.g_x, c.h_x, I(n_x))
#         for i = 1:n
#             E[i, :] = -(R'*c.Ψ[i]*R)[:] # (24), flip the sign for gsylv
#         end
#     end
#     @timeit_debug "Sylvester" begin
#         X = gsylv(A, B, C, D, E) # (22)
#         c.g_xx .= reshape(X[1:n_y, :], n_y, n_x, n_x)
#         c.h_xx .= reshape(X[(n_y+1):end, :], n_x, n_x, n_x)
#     end

#     @timeit_debug "Linear equations for _σσ" begin
#         A_σ = [c.H_yp + c.H_y c.H_xp + c.H_yp * c.g_x]
#         C_σ = zeros(n)
#         η_sq = η * c.Σ * η'
#         H_yp_g = c.H_yp * X[1:n_y, :]
#         R_σ = vcat(c.g_x, zeros(n_y, n_x), I(n_x), zeros(n_x, n_x))
#         for i = 1:n # (29), flip the sign for (34)
#             C_σ[i] -= dot(R_σ' * c.Ψ[i] * R_σ, η_sq)
#             C_σ[i] -= dot(H_yp_g[i, :], η_sq)
#         end
#         X_σ = A_σ \ C_σ # solve (34)
#         c.g_σσ .= X_σ[1:n_y]
#         c.h_σσ .= X_σ[(n_y+1):end]
#     end

#     c.C_0 .= 0.5 * c.Q * vcat(c.g_σσ, zeros(n_x))
#     c.C_2 .= zero(eltype(c.C_2))  # reset as we need to use `+=`
#     for i = 1:n_z
#         for j = 1:n_y
#             c.C_2[i, :, :] += 0.5 * c.Q[i, j] * c.g_xx[j, :, :]
#         end
#     end
#     return :Success
# end

function solve_second_order_p!(
    m::AbstractSecondOrderPerturbationModel,
    c,
    settings,
)

    @unpack n_x, n_y, n_p, n_ϵ, n_z, n, η = m

    @timeit_debug "General Prep" begin
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
        η_sq = η * c.Σ * η'
        gh_stack = vcat(reshape(c.g_xx, n_y, n_x * n_x), reshape(c.h_xx, n_x, n_x * n_x))
        g_xx_flat = reshape(c.g_xx, n_y, n_x * n_x)

        dH = zeros(n, 2n)
        dΨ = [zeros(2n, 2n) for _ = 1:n]
    end

    for i = 1:n_p
        @timeit_debug "Prep for _xx_p" begin
            # Compute the total derivatives
            bar = vcat(c.y_p[:, i], c.y_p[:, i], c.x_p[:, i], c.x_p[:, i])
            Hstack = hcat(c.H_yp_p[i], c.H_y_p[i], c.H_xp_p[i], c.H_x_p[i])
            for j = 1:n
                dH[j, :] = c.Ψ[j] * bar + Hstack[j, :]
                dΨ[j] .= c.Ψ_p[i][j]
                for k = 1:n_y
                    if (c.y_p[k, i] != 0)
                        dΨ[j] += (c.Ψ_yp[k][j] + c.Ψ_y[k][j]) * c.y_p[k, i]
                    end
                end
                for k = 1:n_x
                    if (c.x_p[k, i] != 0)
                        dΨ[j] += (c.Ψ_xp[k][j] + c.Ψ_x[k][j]) * c.x_p[k, i]
                    end
                end
            end

            # Constants: (60)
            R_p = vcat(
                c.g_x_p[i] * c.h_x + c.g_x * c.h_x_p[i],
                c.g_x_p[i],
                c.h_x_p[i],
                zeros(n_x, n_x),
            )
            # Flip the sign of E for Sylvester input
            for j = 1:n
                E[j, :] .= -((R_p'*c.Ψ[j]*R)[:] + (R'*c.Ψ[j]*R_p)[:] + (R'*dΨ[j]*R)[:])
            end
            # Constants: (56)
            E -= hcat(dH[:, 1:n_y], zeros(n, n_x)) * gh_stack * kron(c.h_x, c.h_x) # Plug (57) in (56)
            E -=
                hcat(c.H_yp, zeros(n, n_x)) *
                gh_stack *
                (kron(c.h_x_p[i], c.h_x) + kron(c.h_x, c.h_x_p[i])) # Plug (58) in (56)
            E -=
                hcat(
                    dH[:, (n_y+1):(2*n_y)],
                    dH[:, 1:n_y] * c.g_x +
                    c.H_yp * c.g_x_p[i] +
                    dH[:, (2*n_y+1):(2*n_y+n_x)],
                ) * gh_stack # Plug (59) in (56)
        end

        @timeit_debug "Sylvester for _xx_p" begin
            # Solve the Sylvester equations (56)
            # X = gsylv(A, B, C, D, E)
            Y = adjoint(Q1) * (E * Z2)
            gsylvs!(AS, BS, CS, DS, Y)
            X = Z1 * (Y * adjoint(Q2))
            c.g_xx_p[i] .= reshape(X[1:n_y, :], n_y, n_x, n_x)
            c.h_xx_p[i] .= reshape(X[(n_y+1):end, :], n_x, n_x, n_x)
        end

        @timeit_debug "Prep for _σσ_p" begin
            # Solve _σσ_p
            R_σ_p = vcat(c.g_x_p[i], zeros(n_y + n_x * 2, n_x))
            η_sq_p = η * c.Σ_p[i] * η'
            C_σ =
                -hcat(
                    dH[:, 1:n_y] + dH[:, (n_y+1):(2*n_y)],
                    dH[:, 1:n_y] * c.g_x +
                    c.H_yp * c.g_x_p[i] +
                    dH[:, (2*n_y+1):(2*n_y+n_x)],
                ) * vcat(c.g_σσ, c.h_σσ) # Plug (65) in (64), flip the sign to solve (64)
            C_σ -= (dH[:, 1:n_y] * g_xx_flat + c.H_yp * X[1:n_y, :]) * η_sq[:] # (67), 2nd line
            C_σ -= (c.H_yp * g_xx_flat) * η_sq_p[:] # (67), 3rd line, second part
            for j = 1:n
                C_σ[j] -= dot(
                    (R_σ_p' * c.Ψ[j] * R_σ + R_σ' * c.Ψ[j] * R_σ_p + R_σ' * dΨ[j] * R_σ),
                    η_sq,
                ) # (67), 1st line
                C_σ[j] -= dot((R_σ' * c.Ψ[j] * R_σ), η_sq_p) # (67), 3rd line, first part
            end
        end

        @timeit_debug "Solve _σσ_p" begin
            X_σ = A_σ \ C_σ # solve (64)
            c.g_σσ_p[:, i] .= X_σ[1:n_y]
            c.h_σσ_p[:, i] .= X_σ[(n_y+1):end]

        end
        fill!(c.C_2_p[i], 0.0) # reset as we need to use `+=`
        for j = 1:n_z
            for k = 1:n_y
                c.C_2_p[i][j, :, :] += 0.5 * c.Q[j, k] * c.g_xx_p[i][k, :, :]
            end
        end
        c.A_2_p[i] .= 0.5 * c.h_xx_p[i]
    end

    c.C_0_p .= 0.5 * c.Q * vcat(c.g_σσ_p, zeros(n_x, n_p))
    c.A_0_p .= 0.5 * c.h_σσ_p
    return :Success
end

# Later, can specialize based on the cache later, but not necessary for now
function evaluate_functions!(m, c::AbstractFirstOrderSolverCache, settings, p, solver)
    (settings.print_level > 2) && println("Evaluating first-order functions into cache")
    try
        @unpack y, x = c  # Precondition: valid (y, x) steady states

        m.mod.H_yp!(c.H_yp, y, x, p, solver)
        m.mod.H_y!(c.H_y, y, x, p, solver)
        m.mod.H_xp!(c.H_xp, y, x, p, solver)
        m.mod.H_x!(c.H_x, y, x, p, solver)
        m.mod.Γ!(c.Γ, p, solver)
        maybe_call_function(m.mod.Ω!, c.Ω, p, solver) # supports  m.mod.Ω! = nothing

        if m.mod.n_p > 0
            if !isnothing(c.H_p)  # not required if steady_state_p! there
                m.mod.H_p!(c.H_p, y, x, p, solver)
            end
            m.mod.H_yp_p!(c.H_yp_p, y, x, p, solver)
            m.mod.H_y_p!(c.H_y_p, y, x, p, solver)
            m.mod.H_xp_p!(c.H_xp_p, y, x, p, solver)
            m.mod.H_x_p!(c.H_x_p, y, x, p, solver)
            m.mod.Γ_p!(c.Γ_p, p, solver)
            m.mod.Ψ!(c.Ψ, y, x, p, solver)
            maybe_call_function(m.mod.Ω_p!, c.Ω_p, p, solver)
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
    return :Success  # no failing code-paths yet.
end

function evaluate_functions!(m, c::AbstractSecondOrderSolverCache, settings, p, solver)
    (settings.print_level > 2) && println("Evaluating second-order functions into cache")
    try
        @unpack y, x = c  # Precondition: valid (y, x) steady states

        # m.mod.H_yp!(c.H_yp, y, x, p, solver)
        # m.mod.H_y!(c.H_y, y, x, p, solver)
        # m.mod.H_xp!(c.H_xp, y, x, p, solver)
        # m.mod.H_x!(c.H_x, y, x, p, solver)
        # m.mod.Γ!(c.Γ, p, solver)
        # m.mod.Ψ!(c.Ψ, y, x, p, solver)
        # m.mod.Ψ_yp!(c.Ψ_yp, y, x, p, solver)
        # m.mod.Ψ_y!(c.Ψ_y, y, x, p, solver)
        # m.mod.Ψ_xp!(c.Ψ_xp, y, x, p, solver)
        # m.mod.Ψ_x!(c.Ψ_x, y, x, p, solver)
        maybe_call_function(m.mod.Ω!, c.Ω, p, solver) # supports  m.mod.Ω! = nothing
        if m.mod.n_p > 0
            if !isnothing(c.H_p)  # not required if steady_state_p! there
                m.mod.H_p!(c.H_p, y, x, p, solver)
            end
            m.mod.H_yp_p!(c.H_yp_p, y, x, p, solver)
            m.mod.H_y_p!(c.H_y_p, y, x, p, solver)
            m.mod.H_xp_p!(c.H_xp_p, y, x, p, solver)
            m.mod.H_x_p!(c.H_x_p, y, x, p, solver)
            m.mod.Γ_p!(c.Γ_p, p, solver)
            m.mod.Ψ_p!(c.Ψ_p, y, x, p, solver)
            maybe_call_function(m.mod.Ω_p!, c.Ω_p, p, solver)
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
    return :Success  # no failing code-paths yet.
end
