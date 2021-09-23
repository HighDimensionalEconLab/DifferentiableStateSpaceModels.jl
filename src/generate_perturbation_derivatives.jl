# Fill in the gradients
function evaluate_first_order_functions_p!(m, c, settings, p)
    (settings.print_level > 2) && println("Evaluating first-order function derivatives into cache")
    try
        # TODO: Symbol dispatching
        @unpack y, x = c  # Precondition: valid (y, x) steady states
        if !isnothing(c.H_p)  # not required if steady_state_p! there
            m.mod.H_p!(c.H_p, y, x, p)
        end
        m.mod.H_yp_p!(c.H_yp_p, y, x, p)
        m.mod.H_y_p!(c.H_y_p, y, x, p)
        m.mod.H_xp_p!(c.H_xp_p, y, x, p)
        m.mod.H_x_p!(c.H_x_p, y, x, p)
        m.mod.Γ_p!(c.Γ_p, p)
        maybe_call_function(m.mod.Ω_p!, c.Ω_p, p)
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