using DifferentiableStateSpaceModels: all_equal_struct, get_threadsafe_cache
@testset "Checks cache reuse first-order" begin
    settings = PerturbationSolverSettings(; print_level = 2)
    settings_no_cache = PerturbationSolverSettings(; print_level = 2, use_solution_cache = false)
    m = @include_example_module(Examples.rbc_observables_benchmark)
    p = [0.4, 0.96]
    p_f = [0.7, 0.1, 0.01, 0.0001]
    cache = allocate_cache(m)
    sol_1 = generate_perturbation(m, p; p_f, settings, cache)
    c_1 = get_threadsafe_cache(cache,m, p, p_f).c |>deepcopy
    sol_2 = generate_perturbation(m, p; p_f, settings, cache)
    c_2 = get_threadsafe_cache(cache,m, p, p_f).c |>deepcopy
    sol_3 = generate_perturbation(m, p; p_f, settings = settings_no_cache, cache)
    c_3 = get_threadsafe_cache(cache,m, p, p_f).c |>deepcopy
    @test all_equal_struct(sol_1, sol_2)
    @test all_equal_struct(sol_2, sol_3)
    @test all_equal_struct(c_1, c_2)
    @test all_equal_struct(c_2, c_3)
    cache_0 = allocate_cache(m)  # try with a different cache
    generate_perturbation(m, p; p_f, settings, cache = cache_0)
    c_0 = get_threadsafe_cache(cache,m, p, p_f).c |>deepcopy # before calling anything
    @test all_equal_struct(c_0, c_1)

    # change the `p` and use the cache to check
    p = [0.41, 0.96]
    sol_1 = generate_perturbation(m, p; p_f, settings, cache)
    c_1 = get_threadsafe_cache(cache,m, p, p_f).c |>deepcopy
    cache_0 = allocate_cache(m)  # try with a different cache
    sol_0 = generate_perturbation(m, p; p_f, settings, cache = cache_0)
    c_0 = get_threadsafe_cache(cache,m, p, p_f).c |>deepcopy # before calling anything
    @test all_equal_struct(c_0, c_1)
    @test all_equal_struct(sol_0, sol_1)
end


@testset "Checks cache reuse second-order" begin
    settings = PerturbationSolverSettings(; print_level = 2)
    settings_no_cache = PerturbationSolverSettings(; print_level = 2, use_solution_cache = false)
    m = @include_example_module(Examples.rbc_observables_benchmark, 2)
    p = [0.4, 0.96]
    p_f = [0.7, 0.1, 0.01, 0.0001]
    cache = allocate_cache(m)
    sol_1 = generate_perturbation(m, p; p_f, settings, cache)
    c_1 = get_threadsafe_cache(cache,m, p, p_f).c |>deepcopy
    sol_2 = generate_perturbation(m, p; p_f, settings, cache)
    c_2 = get_threadsafe_cache(cache,m, p, p_f).c |>deepcopy
    sol_3 = generate_perturbation(m, p; p_f, settings = settings_no_cache, cache)
    c_3 = get_threadsafe_cache(cache,m, p, p_f).c |>deepcopy
    @test all_equal_struct(sol_1, sol_2)
    @test all_equal_struct(sol_2, sol_3)
    @test all_equal_struct(c_1, c_2)
    @test all_equal_struct(c_2, c_3)
    cache_0 = allocate_cache(m)  # try with a different cache
    generate_perturbation(m, p; p_f, settings, cache = cache_0)
    c_0 = get_threadsafe_cache(cache,m, p, p_f).c |>deepcopy # before calling anything
    @test all_equal_struct(c_0, c_1)

    # change the `p` and use the cache to check
    p = [0.41, 0.96]
    sol_1 = generate_perturbation(m, p; p_f, settings, cache)
    c_1 = get_threadsafe_cache(cache,m, p, p_f).c |>deepcopy
    cache_0 = allocate_cache(m)  # try with a different cache
    sol_0 = generate_perturbation(m, p; p_f, settings, cache = cache_0)
    c_0 = get_threadsafe_cache(cache,m, p, p_f).c |>deepcopy # before calling anything
    @test all_equal_struct(c_0, c_1)
    @test all_equal_struct(sol_0, sol_1)
end
