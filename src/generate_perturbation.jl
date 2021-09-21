# The generate_perturbation function calculates the perturbation itself
# It can do used without any derivatives overhead (except, perhaps, extra memory in the cache)
function first_order_perturbation(
    m::PerturbationModel,
    p_d, p_f = nothing; cache = SolverCache(m, Val(1), length(p_d)),
    settings = PerturbationSolverSettings()
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