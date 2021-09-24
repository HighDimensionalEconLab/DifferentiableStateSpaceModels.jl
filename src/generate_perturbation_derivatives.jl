# Fill in the gradients
function first_order_perturbation_derivatives!(m, p_d, p_f, cache; settings =PerturbationSolverSettings())
    @assert cache.p_d_symbols == collect(Symbol.(keys(p_d)))

    p = isnothing(p_f) ? p_d : order_vector_by_symbols(merge(p_d, p_f), m.mod.p_symbols)

    # solver type provided to all callbacks
    solver = PerturbationSolver(m, cache, settings)

    # Fill in derivatives
    ret = evaluate_first_order_functions_p!(m, cache, settings, p)
    (ret == :Success) || return ret
    ret = solve_first_order_p!(m, cache, settings)
    maybe_call_function(settings.solve_first_order_p_callback, ret, m, cache, settings)
    (ret == :Success) || return ret
end

function second_order_perturbation_derivatives!(m, p_d, p_f, cache; settings =PerturbationSolverSettings())
    @assert cache.p_d_symbols == collect(Symbol.(keys(p_d)))

    p = isnothing(p_f) ? p_d : order_vector_by_symbols(merge(p_d, p_f), m.mod.p_symbols)

    # solver type provided to all callbacks
    solver = PerturbationSolver(m, cache, settings)

    # Fill in derivatives
    ret = evaluate_first_order_functions_p!(m, cache, settings, p)
    (ret == :Success) || return ret
    ret = evaluate_second_order_functions_p!(m, cache, settings, p)
    (ret == :Success) || return ret
    ret = solve_first_order_p!(m, cache, settings)
    maybe_call_function(settings.solve_first_order_p_callback, ret, m, cache, settings)
    (ret == :Success) || return ret    
end


function evaluate_first_order_functions_p!(m, c, settings, p)
    (settings.print_level > 2) && println("Evaluating first-order function derivatives into cache")
    try        
        @unpack y, x = c  # Precondition: valid (y, x) steady states
        isnothing(c.H_p) || fill_array_by_symbol_dispatch(m.mod.H_p!, c.H_p, c.p_d_symbols,  y, x, p) #not required if steady_state_p!
        fill_array_by_symbol_dispatch(m.mod.H_yp_p!, c.H_yp_p, c.p_d_symbols,  y, x, p)
        fill_array_by_symbol_dispatch(m.mod.H_y_p!, c.H_y_p, c.p_d_symbols,  y, x, p)
        fill_array_by_symbol_dispatch(m.mod.H_xp_p!, c.H_xp_p, c.p_d_symbols,  y, x, p)
        fill_array_by_symbol_dispatch(m.mod.H_x_p!, c.H_x_p, c.p_d_symbols,  y, x, p)
        fill_array_by_symbol_dispatch(m.mod.Γ_p!, c.Γ_p, c.p_d_symbols,p)
        isnothing(c.Ω_p) || fill_array_by_symbol_dispatch(m.mod.Ω_p!, c.Ω_p, c.p_d_symbols,p)
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
    (settings.print_level > 2) && println("Evaluating second-order function derivatives into cache")
    try        
        @unpack y, x = c  # Precondition: valid (y, x) steady states
        fill_array_by_symbol_dispatch(m.mod.Ψ_p!, c.Ψ_p, c.p_d_symbols,  y, x, p)
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


function solve_first_order_p!(m,c,settings)
    @unpack n_x, n_y, n_ϵ = m
    n = n_x + n_y
    n_p = length(c.p_d_symbols)
    (settings.print_level > 2) && println("Solving first order derivatives of perturbation")


    try
        if isnothing(m.mod.ȳ_p!) && isnothing(m.mod.x̄_p!)
            # Zeroth-order derivatives if not provided
            # Calculating c.y_p, c.x_p
            A_zero = [c.H_y + c.H_yp c.H_x + c.H_xp]
            # H_p is a vector of vectors
            x_zeroth = A_zero \ -hcat(c.H_p...) # (47)
            c.y_p .= x_zeroth[1:n_y, :]
            c.x_p .= x_zeroth[(n_y+1):n, :]
        end

        # Write equation (52) as E + AX + CXD = 0, a generalized Sylvester equation
        # first-order derivatives
        R = vcat(c.g_x * c.h_x, c.g_x, c.h_x, I(n_x))
        A = [c.H_y c.H_xp + c.H_yp * c.g_x]
        B = Array(I(n_x) * 1.0)
        C = [c.H_yp zeros(n, n_x)]
        D = c.h_x
        AS, CS, Q1, Z1 = schur(A, C)
        BS, DS, Q2, Z2 = schur(B, D)
        # Initialize
        dH = zeros(2n, n)
        bar = zeros(2n, 1)
        Hstack = zeros(n, 2n)
        for i = 1:n_p
            # p-specific Sylvester preparation
            bar[1:n_y] = c.y_p[i]
            bar[(n_y+1):(2*n_y)] = c.y_p[i]
            bar[(2*n_y+1):(2*n_y+n_x)] = c.x_p[i]
            bar[(2*n_y+n_x+1):end] = c.x_p[i]

            Hstack[:, 1:n_y] = c.H_yp_p[i]
            Hstack[:, (n_y+1):(2*n_y)] = c.H_y_p[i]
            Hstack[:, (2*n_y+1):(2*n_y+n_x)] = c.H_xp_p[i]
            Hstack[:, (2*n_y+n_x+1):end] = c.H_x_p[i]
            for j = 1:n
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
            c.h_x_p[i] .= X[(n_y+1):n, :]

            # Q weighted derivatives
            c.C_1_p[i] .= c.Q * vcat(c.g_x_p[i], zeros(n_x, n_x))
            c.A_1_p[i] .= c.h_x_p[i]

            # V derivatives
            tmp = c.h_x_p[i] * Array(c.V) * c.h_x'
            c.V_p[i] .= lyapd(c.h_x, c.η * c.Σ_p[i] * c.η' + tmp + tmp')
            # B derivatives
            c.B_p[i] .= c.η * c.Γ_p[i]
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


function solve_second_order_p!(m, c, settings)
        # The derivatives
        for i = 1:n_p
            c.Σ_p[i] .= Symmetric(c.Γ_p[i] * c.Γ' + c.Γ * c.Γ_p[i]')
        end
end