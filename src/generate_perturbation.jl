
# Utility to call the cache's organized by symbol
function fill_array_by_symbol_dispatch(f, c, symbols, args...)
    for (i, sym) in enumerate(symbols)
        f(c[i], Val(sym), args...)
    end
end

# The generate_perturbation function calculates the perturbation itself
# It can do used without any derivatives overhead (except, perhaps, extra memory in the cache)
function first_order_perturbation(m::PerturbationModel, p_d, p_f=nothing;
                                  cache=SolverCache(m, Val(1), collect(Symbol.(keys(p_d)))),
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

    # @timeit_debug "solve_first_order" begin
    #     ret = solve_first_order!(m, cache, settings)
    # end
    # maybe_call_function(settings.solve_first_order_callback, ret, m, cache, settings)
    # (ret == :Success) || return FirstOrderPerturbationSolution(ret, m, cache)

    # @timeit_debug "solve_first_order_p" begin
    #     ret = solve_first_order_p!(m, cache, settings)
    # end
    # maybe_call_function(settings.solve_first_order_p_callback, ret, m, cache, settings)
    # (ret == :Success) || return FirstOrderPerturbationSolution(ret, m, cache)

    return FirstOrderPerturbationSolution(:Success, m, cache)
end

# The generate_perturbation function calculates the perturbation itself
# It can do used without any derivatives overhead (except, perhaps, extra memory in the cache)
function second_order_perturbation(m::PerturbationModel, p_d, p_f=nothing;
                                   cache=SolverCache(m, Val(2), collect(Symbol.(keys(p_d)))),
                                   settings=PerturbationSolverSettings())
    @assert cache.p_d_symbols == collect(Symbol.(keys(p_d)))

    p = isnothing(p_f) ? p_d : order_vector_by_symbols(merge(p_d, p_f), m.mod.p_symbols)

    # solver type provided to all callbacks
    solver = PerturbationSolver(m, cache, settings)
    ret = calculate_steady_state!(m, cache, settings, p)
    maybe_call_function(settings.calculate_steady_state_callback, ret, m, cache, settings,
                        p, solver)
    (ret == :Success) || return SecondOrderPerturbationSolution(ret, m, cache)

    ret = evaluate_first_order_functions!(m, cache, settings, p)
    (ret == :Success) || return SecondOrderPerturbationSolution(ret, m, cache)
    ret = evaluate_second_order_functions!(m, cache, settings, p)
    (ret == :Success) || return SecondOrderPerturbationSolution(ret, m, cache)
    maybe_call_function(settings.evaluate_functions_callback, ret, m, cache, settings, p,
                        solver)

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
            isnothing(m.mod.ȳ_p!) || fill_array_by_symbol_dispatch(m.mod.ȳ_p!, c.y_p, c.p_d_symbols, p)
            isnothing(m.mod.x̄_p!) || fill_array_by_symbol_dispatch(m.mod.x̄_p!, c.x_p, c.p_d_symbols, p)
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