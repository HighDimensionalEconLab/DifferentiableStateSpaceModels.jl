using DifferentiableStateSpaceModels, BenchmarkTools, Zygote, Random

const RBC_LIKELIHOODS = BenchmarkGroup()
const RBC_LIKELIHOODS["kalman"] = BenchmarkGroup()
const RBC_LIKELIHOODS["joint"] = BenchmarkGroup()
const RBC_LIKELIHOODS["second_order_joint"] = BenchmarkGroup()
const RBC_LIKELIHOODS["perturbation"] = BenchmarkGroup()
const RBC_LIKELIHOODS["solve"] = BenchmarkGroup()

const m_rbc = @include_example_module(Examples.rbc_observables_benchmark)
const m_second_rbc = @include_example_module(Examples.rbc_observables_benchmark, 2)

function generate_benchmark_data(
    m::FirstOrderPerturbationModel,
    p,
    p_f,
    T,
    x0,
    ϵ = [randn(m.n_ϵ) for _ = 1:T],
)
    return (
        solve(generate_perturbation(m, p; p_f), x0, (0, T), LTI(); noise = ϵ).z[2:end],
        ϵ,
    )
end

function generate_benchmark_data(
    m::SecondOrderPerturbationModel,
    p,
    p_f,
    T,
    x0,
    ϵ = [randn(m.n_ϵ) for _ = 1:T],
)
    return (
        solve(generate_perturbation(m, p; p_f), x0, (0, T), QTI(); noise = ϵ).z[2:end],
        ϵ,
    )
end
const p_rbc = [0.4, 0.96]
const p_f_rbc = [0.7, 0.1, 0.01, 0.0001]
const T = 100
const x0_rbc = zeros(m_rbc.n_x)
const z_rbc, ϵ_rbc = generate_benchmark_data(m_rbc, p_rbc, p_f_rbc, T, x0_rbc)
const x0_second_rbc = zeros(m_rbc.n_x)
const z_second_rbc, ϵ_second_rbc =
    generate_benchmark_data(m_second_rbc, p_rbc, p_f_rbc, T, x0_second_rbc)
const cache_rbc = allocate_cache(m_rbc)
const cache_second_rbc = allocate_cache(m_second_rbc)
const settings_skip_hash = PerturbationSolverSettings(use_solution_cache = false)  #reuses cache storage but doesn't check cache


# conditional likelihood
function benchmark_likelihood_kalman(p, p_f, z, m, cache)
    sol = generate_perturbation(m, p; p_f, cache)
    return solve(sol, sol.x_ergodic, (0, length(z)); observables = z).logpdf
end

# Just the perturbation DifferentiableStateSpaceModels
# WE do not expect the reuse of the memory in the cache to save much time for small models like this
RBC_LIKELIHOODS["perturbation"]["first_order_perturbation"] =
    @benchmarkable generate_perturbation($m_rbc, $p_rbc; p_f = $p_f_rbc, cache) setup =
        (cache = allocate_cache(m_rbc))
RBC_LIKELIHOODS["perturbation"]["first_order_perturbation_reused_cache"] =
    @benchmarkable generate_perturbation($m_rbc, $p_rbc; p_f = $p_f_rbc, cache = $cache_rbc)  # shouldn't do many calculations  since we didn't change the p or p_f
RBC_LIKELIHOODS["perturbation"]["first_order_perturbation_reused_cache_no_hash"] =
    @benchmarkable generate_perturbation(
        $m_rbc,
        $p_rbc;
        p_f = $p_f_rbc,
        cache = $cache_rbc,
        settings = $settings_skip_hash,
    )  #reuse memory but no skipping perturbation
RBC_LIKELIHOODS["perturbation"]["second_order_perturbation"] =
    @benchmarkable generate_perturbation($m_second_rbc, $p_rbc; p_f = $p_f_rbc, cache) setup =
        (cache = allocate_cache(m_second_rbc))
RBC_LIKELIHOODS["perturbation"]["second_order_perturbation_reused_cache"] =
    @benchmarkable generate_perturbation(
        $m_second_rbc,
        $p_rbc;
        p_f = $p_f_rbc,
        cache = $cache_second_rbc,
    )  # shouldn't do any significant calculations  since we didn't change the p or p_f
RBC_LIKELIHOODS["perturbation"]["second_order_perturbation_reused_cache_no_hash"] =
    @benchmarkable generate_perturbation(
        $m_second_rbc,
        $p_rbc;
        p_f = $p_f_rbc,
        cache = $cache_second_rbc,
        settings = $settings_skip_hash,
    )  #reuse memory but no skipping perturbation


# Testing the solve isolated sequential 
const sol_rbc = generate_perturbation(m_rbc, p_rbc; p_f = p_f_rbc)
const sol_second_rbc = generate_perturbation(m_second_rbc, p_rbc; p_f = p_f_rbc)
const tspan = (0, T)
RBC_LIKELIHOODS["solve"]["kalman"] =
    @benchmarkable solve($sol_rbc, $sol_rbc.x_ergodic, $tspan; observables = $z_rbc).logpdf
# RBC_LIKELIHOODS["solve"]["kalman_gradient"] = @benchmarkable gradient((A, B) -> solve(A, B, $sol_rbc.x_ergodic, $tspan; h = $sol_rbc.C, D = $sol_rbc.D, observables = $z_rbc).logpdf, $sol_rbc.A, $sol_rbc.B)


# for testing the solve on its own
test_evolution(x, p, t) = p.A * x
test_volatility(x, p, t) = p.B
test_observation(x, p, t) = p.C * x
const p_test = (; sol_rbc.A, sol_rbc.B, sol_rbc.C)
RBC_LIKELIHOODS["solve"]["first_order_joint"] = @benchmarkable solve(
    test_evolution,
    test_volatility,
    $x0_rbc,
    $tspan,
    $p_test;
    h = test_observation,
    D = $sol_rbc.D,
    observables = $z_rbc,
    noise = $ϵ_rbc,
).logpdf
RBC_LIKELIHOODS["solve"]["first_order_joint_gradient"] = @benchmarkable gradient(
    p_test ->
        solve(
            test_evolution,
            test_volatility,
            $x0_rbc,
            $tspan,
            p_test;
            h = test_observation,
            D = $sol_rbc.D,
            observables = $z_rbc,
            noise = $ϵ_rbc,
        ).logpdf,
    $p_test,
)
RBC_LIKELIHOODS["solve"]["first_order_joint_gradient_ϵ"] = @benchmarkable gradient(
    (p_test, ϵ) ->
        solve(
            test_evolution,
            test_volatility,
            $x0_rbc,
            $tspan,
            p_test;
            h = test_observation,
            D = $sol_rbc.D,
            observables = $z_rbc,
            noise = ϵ,
        ).logpdf,
    $p_test,
    $ϵ_rbc,
)
RBC_LIKELIHOODS["solve"]["second_order_joint"] = @benchmarkable solve(
    $sol_second_rbc,
    $x0_second_rbc,
    $tspan;
    observables = $z_second_rbc,
    noise = $ϵ_second_rbc,
).logpdf
RBC_LIKELIHOODS["solve"]["second_order_joint_gradient_ϵ"] = @benchmarkable gradient(
    ϵ ->
        solve(
            $sol_second_rbc,
            $x0_second_rbc,
            $tspan;
            observables = $z_second_rbc,
            noise = ϵ,
        ).logpdf,
    $ϵ_second_rbc,
)

# likelihood and gradient with cache reallocated each time
RBC_LIKELIHOODS["kalman"]["likelihood"] =
    @benchmarkable benchmark_likelihood_kalman($p_rbc, $p_f_rbc, $z_rbc, $m_rbc, cache) setup =
        (cache = allocate_cache(m_rbc))
RBC_LIKELIHOODS["kalman"]["likelihood_gradient"] = @benchmarkable gradient(
    p -> benchmark_likelihood_kalman(p, $p_f_rbc, $z_rbc, $m_rbc, cache),
    $p_rbc,
) setup = (cache = allocate_cache(m_rbc))

# filled cache
generate_perturbation(m_rbc, p_rbc; p_f = p_f_rbc, cache = cache_rbc)  # solve to fill it
RBC_LIKELIHOODS["kalman"]["likelihood_reused_cache"] =
    @benchmarkable benchmark_likelihood_kalman($p_rbc, $p_f_rbc, $z_rbc, $m_rbc, $cache_rbc)

# gradient of kalman with a cache already filled
const cache_grad_rbc = allocate_cache(m_rbc)
generate_perturbation(m_rbc, p_rbc; p_f = p_f_rbc, cache = cache_grad_rbc)
RBC_LIKELIHOODS["kalman"]["likelihood_gradient_reused_cache"] = @benchmarkable gradient(
    p -> benchmark_likelihood_kalman(p, $p_f_rbc, $z_rbc, $m_rbc, $cache_grad_rbc),
    $p_rbc,
)

# joint-likelihood
function benchmark_joint_likelihood(p, ϵ, x0, p_f, z, m, cache)
    sol = generate_perturbation(m, p; p_f, cache)
    # (sol.retcode == :Success) || return 0.0
    return solve(sol, x0, (0, length(z)); noise = ϵ, observables = z).logpdf
end

# likelihood and gradient with cache reallocated each time
RBC_LIKELIHOODS["joint"]["likelihood"] = @benchmarkable benchmark_joint_likelihood(
    $p_rbc,
    $ϵ_rbc,
    $x0_rbc,
    $p_f_rbc,
    $z_rbc,
    $m_rbc,
    cache,
) setup = (cache = allocate_cache(m_rbc))
RBC_LIKELIHOODS["joint"]["likelihood_gradient_p"] = @benchmarkable gradient(
    p ->
        benchmark_joint_likelihood(p, $ϵ_rbc, $x0_rbc, $p_f_rbc, $z_rbc, $m_rbc, cache),
    $p_rbc,
) setup = (cache = allocate_cache(m_rbc))
RBC_LIKELIHOODS["joint"]["likelihood_gradient_p_ϵ"] = @benchmarkable gradient(
    (p, ϵ) ->
        benchmark_joint_likelihood(p, ϵ, $x0_rbc, $p_f_rbc, $z_rbc, $m_rbc, cache),
    $p_rbc,
    $ϵ_rbc,
) setup = (cache = allocate_cache(m_rbc))
RBC_LIKELIHOODS["joint"]["likelihood_gradient_ϵ"] = @benchmarkable gradient(
    ϵ ->
        benchmark_joint_likelihood($p_rbc, ϵ, $x0_rbc, $p_f_rbc, $z_rbc, $m_rbc, cache),
    $ϵ_rbc,
) setup = (cache = allocate_cache(m_rbc))


# gradient with a cache already filled
const cache_rbc_joint_grad = allocate_cache(m_rbc)
generate_perturbation(m_rbc, p_rbc; p_f = p_f_rbc, cache = cache_rbc_joint_grad)
RBC_LIKELIHOODS["joint"]["likelihood_reused"] = @benchmarkable benchmark_joint_likelihood(
    $p_rbc,
    $ϵ_rbc,
    $x0_rbc,
    $p_f_rbc,
    $z_rbc,
    $m_rbc,
    $cache_grad_rbc,
)
RBC_LIKELIHOODS["joint"]["likelihood_gradient_p_reused"] = @benchmarkable gradient(
    p -> benchmark_joint_likelihood(
        p,
        $ϵ_rbc,
        $x0_rbc,
        $p_f_rbc,
        $z_rbc,
        $m_rbc,
        $cache_grad_rbc,
    ),
    $p_rbc,
)
RBC_LIKELIHOODS["joint"]["likelihood_gradient_p_ϵ_reused"] = @benchmarkable gradient(
    (p, ϵ) -> benchmark_joint_likelihood(
        p,
        ϵ,
        $x0_rbc,
        $p_f_rbc,
        $z_rbc,
        $m_rbc,
        $cache_grad_rbc,
    ),
    $p_rbc,
    $ϵ_rbc,
)
RBC_LIKELIHOODS["joint"]["likelihood_gradient_ϵ_reused"] = @benchmarkable gradient(
    ϵ -> benchmark_joint_likelihood(
        $p_rbc,
        ϵ,
        $x0_rbc,
        $p_f_rbc,
        $z_rbc,
        $m_rbc,
        $cache_grad_rbc,
    ),
    $ϵ_rbc,
)

# Second order
# likelihood and gradient with cache reallocated each time
RBC_LIKELIHOODS["second_order_joint"]["likelihood"] =
    @benchmarkable benchmark_joint_likelihood(
        $p_rbc,
        $ϵ_second_rbc,
        $x0_second_rbc,
        $p_f_rbc,
        $z_second_rbc,
        $m_second_rbc,
        cache,
    ) setup = (cache = allocate_cache(m_second_rbc))
RBC_LIKELIHOODS["second_order_joint"]["likelihood_gradient_p"] = @benchmarkable gradient(
    p -> benchmark_joint_likelihood(
        p,
        $ϵ_second_rbc,
        $x0_second_rbc,
        $p_f_rbc,
        $z_second_rbc,
        $m_second_rbc,
        cache,
    ),
    $p_rbc,
) setup = (cache = allocate_cache(m_second_rbc))
RBC_LIKELIHOODS["second_order_joint"]["likelihood_gradient_p_ϵ"] = @benchmarkable gradient(
    (p, ϵ) -> benchmark_joint_likelihood(
        p,
        ϵ,
        $x0_second_rbc,
        $p_f_rbc,
        $z_second_rbc,
        $m_second_rbc,
        cache,
    ),
    $p_rbc,
    $ϵ_second_rbc,
) setup = (cache = allocate_cache(m_second_rbc))
RBC_LIKELIHOODS["second_order_joint"]["likelihood_gradient_ϵ"] = @benchmarkable gradient(
    ϵ -> benchmark_joint_likelihood(
        $p_rbc,
        ϵ,
        $x0_second_rbc,
        $p_f_rbc,
        $z_second_rbc,
        $m_second_rbc,
        cache,
    ),
    $ϵ_second_rbc,
) setup = (cache = allocate_cache(m_second_rbc))


# gradient with a cache already filled
const cache_second_grad_rbc = allocate_cache(m_second_rbc)
generate_perturbation(m_second_rbc, p_rbc; p_f = p_f_rbc, cache = cache_second_grad_rbc)
RBC_LIKELIHOODS["second_order_joint"]["likelihood_reused"] =
    @benchmarkable benchmark_joint_likelihood(
        $p_rbc,
        $ϵ_second_rbc,
        $x0_second_rbc,
        $p_f_rbc,
        $z_second_rbc,
        $m_second_rbc,
        $cache_second_grad_rbc,
    )
RBC_LIKELIHOODS["second_order_joint"]["likelihood_gradient_p_reused"] =
    @benchmarkable gradient(
        p -> benchmark_joint_likelihood(
            p,
            $ϵ_second_rbc,
            $x0_second_rbc,
            $p_f_rbc,
            $z_second_rbc,
            $m_second_rbc,
            $cache_second_grad_rbc,
        ),
        $p_rbc,
    )
RBC_LIKELIHOODS["second_order_joint"]["likelihood_gradient_p_ϵ_reused"] =
    @benchmarkable gradient(
        (p, ϵ) -> benchmark_joint_likelihood(
            p,
            ϵ,
            $x0_second_rbc,
            $p_f_rbc,
            $z_second_rbc,
            $m_second_rbc,
            $cache_second_grad_rbc,
        ),
        $p_rbc,
        $ϵ_second_rbc,
    )
RBC_LIKELIHOODS["second_order_joint"]["likelihood_gradient_ϵ_reused"] =
    @benchmarkable gradient(
        ϵ -> benchmark_joint_likelihood(
            $p_rbc,
            ϵ,
            $x0_second_rbc,
            $p_f_rbc,
            $z_second_rbc,
            $m_second_rbc,
            $cache_second_grad_rbc,
        ),
        $ϵ_second_rbc,
    )

# return the RBC_LIKELIHOODS benchmarks to the parent
RBC_LIKELIHOODS
