# Fill in the gradients
function first_order_perturbation_derivatives!(m, p_d, p_f, cache; settings =PerturbationSolverSettings())
    @assert cache.p_d_symbols == collect(Symbol.(keys(p_d)))

    p = isnothing(p_f) ? p_d : order_vector_by_symbols(merge(p_d, p_f), m.mod.p_symbols)

    # solver type provided to all callbacks
    solver = PerturbationSolver(m, cache, settings)

    # Fill in derivatives
    ret = evaluate_first_order_functions_p!(m, cache, settings, p)
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