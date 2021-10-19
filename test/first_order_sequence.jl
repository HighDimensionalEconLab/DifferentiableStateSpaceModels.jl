using DifferentiableStateSpaceModels, LinearAlgebra, Test, Zygote
using DifferentiableStateSpaceModels.Examples
using FiniteDiff: finite_difference_gradient

@testset "Sequence Simulation, 1st order" begin
    m = @include_example_module(Examples.rbc_observables)
    p_f = (ρ=0.2, δ=0.02, σ=0.01, Ω_1=0.01)
    p_d = (α=0.5, β=0.95)

    c = SolverCache(m, Val(1), p_d)
    sol = generate_perturbation(m, p_d, p_f; cache = c)

    T = 9
    eps_value = [[0.22], [0.01], [0.14], [0.03], [0.15], [0.21], [0.22], [0.05], [0.18]]
    x0 = zeros(m.n_x)
    simul = solve(
        DifferentiableStateSpaceModels.dssm_evolution,
        DifferentiableStateSpaceModels.dssm_volatility,
        x0,
        (0, T),
        sol;
        h = DifferentiableStateSpaceModels.dssm_observation,
        noise = eps_value,
    )
    @test simul.z[2:end] ≈ [[-0.0014843113235688628, 0.0], [-0.001672969226342977, -0.013660616212663025], [-0.0025907902434120613, -0.016424018091334355], [-0.002808352075911724, -0.02507880927301923], [-0.003749820417827369, -0.027731843802628737], [-0.00514124775142971, -0.036595970922849726], [-0.006607765317710925, -0.05006822400565087], [-0.006885970365182856, -0.06457803532896965], [-0.00789049635407183, -0.06822941930528102]]
    @test simul.z ≈ solve(sol, x0, (0, T), DifferentiableStateSpaceModels.LTI(); noise = eps_value).z

    # inference
    @inferred solve(
        DifferentiableStateSpaceModels.dssm_evolution,
        DifferentiableStateSpaceModels.dssm_volatility,
        x0,
        (0, T),
        sol;
        h = DifferentiableStateSpaceModels.dssm_observation,
        noise = eps_value,
    )
    # @inferred generate_perturbation(m, p_d, p_f)
end

@testset "Gradients, generate_perturbation + simulation, 1st order" begin
    m = @include_example_module(Examples.rbc_observables)
    p_f = (ρ=0.2, δ=0.02, σ=0.01, Ω_1=0.01)
    p_d = (α=0.5, β=0.95)
    p_d_input = [0.5, 0.95]

    T = 9
    ϵ_mat = [0.22, 0.01, 0.14, 0.03, 0.15, 0.21, 0.22, 0.05, 0.18]
    x0 = zeros(m.n_x)

    function sum_test_joint_first(p_d_input, ϵ_mat, x0, T, p_f, m; kwargs...)
        p_d = (α=p_d_input[1], β=p_d_input[2])
        sol = generate_perturbation(m, p_d, p_f; kwargs...)
        ϵ = map(i -> ϵ_mat[i:i], 1:T)
        simul = solve(
            DifferentiableStateSpaceModels.dssm_evolution,
            DifferentiableStateSpaceModels.dssm_volatility,
            x0,
            (0, T),
            sol;
            h = DifferentiableStateSpaceModels.dssm_observation,
            noise = ϵ,
        )
        return sum(sum(simul.z))
    end
    settings = PerturbationSolverSettings()
    # The function I wrote here is not type stable
    # @inferred sum_test_joint_first(p_d_input, ϵ_mat, x0, T, p_f, m; settings)
    res_zygote = gradient(
        (p_d_input, ϵ_mat) -> sum_test_joint_first(p_d_input, ϵ_mat, x0, T, p_f, m; settings),
        p_d_input,
        ϵ_mat
    )
    # p
    res_finite = finite_difference_gradient(
       p_d_input -> sum_test_joint_first(p_d_input, ϵ_mat, x0, T, p_f, m; settings),
       p_d_input,
    )
    @test res_zygote[1] ≈ res_finite
    # ϵ
    res_finite = finite_difference_gradient(
        ϵ_mat -> sum_test_joint_first(p_d_input, ϵ_mat, x0, T, p_f, m; settings),
        ϵ_mat,
    )
    @test res_zygote[2] ≈ res_finite

    # inference
    # @inferred sum_test_joint_first(p_d_input, ϵ_mat, x0, T, p_f, m; settings)
end

function kalman_test(p, p_f, m, cache, z, tspan)
    sol = generate_perturbation(m, p; p_f, cache)
    return solve(sol, sol.x_ergodic, tspan; observables = z).logpdf
end

@testset "Kalman filter and its gradient" begin
    m = @include_example_module(Examples.rbc_observables_benchmark)

    p_f = [0.2, 0.02, 0.01, sqrt(0.01)]
    p = [0.5, 0.95]
    z = [
        [-0.6949847708598687, -0.8456988740809867],
        [-0.7804117657996692, 0.07781473603479207],
        [-1.1363021614363802, -2.41253450179418],
        [-0.2140813001516194, -0.10914617826240575],
        [-1.0365874981404577, 0.9869373465251516],
        [-0.7321498641416826, 0.012293325072265942],
        [-0.054809260599132194, -1.8233591236618099],
        [0.5407452466493482, -0.9773559802938866],
        [1.3968232347532277, -2.139194998843768],
        [1.3427856576886874, -0.3100476471887863],
    ]

    tspan = (0, length(z))
    cache = allocate_cache(m)  # reuse the cache
    sol = generate_perturbation(m, p; p_f)
    res = kalman_test(p, p_f, m, cache, z, tspan)
    @test res ≈ -805.7558351251781
    res = gradient(p -> kalman_test(p, p_f, m, cache, z, tspan), p)
    res_finite = finite_difference_gradient(p -> kalman_test(p, p_f, m, cache, z, tspan), p)
    @test res_finite ≈ res[1]

    # inference
    @inferred solve(sol, sol.x_ergodic, tspan; observables = z)
    @inferred kalman_test(p, p_f, m, cache, z, tspan)
end

function likelihood_test_joint_first(p, p_f, ϵ, x0, m, tspan, z)
    sol = generate_perturbation(m, p; p_f)
    return solve(
        dssm_evolution,
        dssm_volatility,
        x0,
        tspan,
        sol;
        observables = z,
        h = dssm_observation,
        sol.D,
        noise = ϵ,
    ).logpdf
end
likelihood_test_joint_first_sol(p, p_f, ϵ, x0, m, tspan, z) =
    solve(generate_perturbation(m, p; p_f), x0, tspan; observables = z, noise = ϵ).logpdf

@testset "Gradients, generate_perturbation + likelihood, 1st order" begin
    m = @include_example_module(Examples.rbc_observables_benchmark)
    p_f = [0.2, 0.02, 0.01, sqrt(0.01)]
    p = [0.5, 0.95]
    ϵ = [[0.22], [0.01], [0.14], [0.03], [0.15], [0.21], [0.22], [0.05], [0.18]]
    z = [
        [-0.6949847708598687, -0.8456988740809867],
        [-0.7804117657996692, 0.07781473603479207],
        [-1.1363021614363802, -2.41253450179418],
        [-0.2140813001516194, -0.10914617826240575],
        [-1.0365874981404577, 0.9869373465251516],
        [-0.7321498641416826, 0.012293325072265942],
        [-0.054809260599132194, -1.8233591236618099],
        [0.5407452466493482, -0.9773559802938866],
        [1.3968232347532277, -2.139194998843768],
    ]
    x0 = zeros(m.n_x)
    tspan = (0, length(z))

    res = gradient((p, ϵ) -> likelihood_test_joint_first(p, p_f, ϵ, x0, m, tspan, z), p, ϵ)
    @test res[1] ≈ [303.7133186356109, 553.6149537473261]
    @test res[2] ≈ [
        [40.62454806083384],
        [39.38899479341156],
        [25.095297618483304],
        [26.06697625612332],
        [33.10959536324157],
        [31.484308705831474],
        [19.172319198105615],
        [11.464791870737214],
        [-0.9477420442978448],
    ]

    #lifting the solution object directly.

    res2 = gradient(
        (p, ϵ) -> likelihood_test_joint_first_sol(p, p_f, ϵ, x0, m, tspan, z),
        p,
        ϵ,
    )
    @test res[1] ≈ res2[1]
    @test res[2] ≈ res2[2]

    # inference
    sol = generate_perturbation(m, p; p_f)
    @inferred solve(
        dssm_evolution,
        dssm_volatility,
        x0,
        tspan,
        sol;
        observables = z,
        h = dssm_observation,
        sol.D,
        noise = ϵ,
    )
    @inferred likelihood_test_joint_first(p, p_f, ϵ, x0, m, tspan, z)
end

function minimal_likelihood_test_kalman_first(A, B, C, D, u0, noise, observables)
    return solve(
        A,
        B,
        C,
        D,
        u0,
        (0, length(observables)),
        LTILikelihood();
        noise = nothing,
        observables,
    ).logpdf
end

@testset "FVGQ20 Kalman likelhood derivative in 1st order" begin
    path = joinpath(pkgdir(DifferentiableStateSpaceModels), "test", "data")
    file_prefix = "FVGQ20"
    A = Matrix(DataFrame(CSV.File(joinpath(path, "$(file_prefix)_A.csv"), header = false)))
    B = Matrix(DataFrame(CSV.File(joinpath(path, "$(file_prefix)_B.csv"), header = false)))
    C = Matrix(DataFrame(CSV.File(joinpath(path, "$(file_prefix)_C.csv"), header = false)))
    D_raw =
        Matrix(DataFrame(CSV.File(joinpath(path, "$(file_prefix)_D.csv"), header = false)))
    D = Turing.TuringDiagMvNormal(zero(vec(D_raw)), vec(D_raw))
    observables_raw = Matrix(
        DataFrame(
            CSV.File(joinpath(path, "$(file_prefix)_observables.csv"), header = false),
        ),
    )
    noise_raw = Matrix(
        DataFrame(CSV.File(joinpath(path, "$(file_prefix)_noise.csv"), header = false)),
    )
    observables = [observables_raw[i, :] for i = 1:size(observables_raw, 1)]
    noise = [noise_raw[i, :] for i = 1:size(noise_raw, 1)]
    u0_raw = Matrix(
        DataFrame(CSV.File(joinpath(path, "$(file_prefix)_ergodic.csv"), header = false)),
    )
    u0 = DistributionsAD.TuringDenseMvNormal(
        zeros(size(u0_raw, 1)),
        cholesky(Symmetric(u0_raw)),
    )

    minimal_likelihood_test_kalman_first(A, B, C, D, u0, noise, observables)

    res = gradient(minimal_likelihood_test_kalman_first, A, B, C, D, u0, noise, observables)

    # Some tests
    @test finite_difference_gradient(
        A -> minimal_likelihood_test_kalman_first(A, B, C, D, u0, noise, observables),
        A,
    ) ≈ res[1] rtol = 1E-3
    @test finite_difference_gradient(
        B -> minimal_likelihood_test_kalman_first(A, B, C, D, u0, noise, observables),
        B,
    ) ≈ res[2] rtol = 1E-3
    @test finite_difference_gradient(
        C -> minimal_likelihood_test_kalman_first(A, B, C, D, u0, noise, observables),
        C,
    ) ≈ res[3] rtol = 1E-3

    # missing test for the D = MvNormal.  No finite_difference_gradient support yet.

    # @test finite_difference_gradient(u0 -> minimal_likelihood_test_kalman_first(A, B, C, D, u0, noise, observables), u0) ≈ res[5] rtol=1E-7

    observables_grad = finite_difference_gradient(
        observables_mat -> minimal_likelihood_test_kalman_first(
            A,
            B,
            C,
            D,
            u0,
            noise,
            [observables_mat[i, :] for i = 1:size(observables_mat, 1)],
        ),
        observables_raw,
    )
    @test [observables_grad[i, :] for i = 1:size(observables_raw, 1)] ≈ res[7] rtol = 1E-5

    # inference
    @inferred solve(
        A,
        B,
        C,
        D,
        u0,
        (0, length(observables)),
        LTILikelihood();
        noise = nothing,
        observables,
    )
end

function minimal_likelihood_test_joint_first(A, B, C, D, u0, noise, observables)
    return solve(
        A,
        B,
        C,
        D,
        u0,
        (0, length(observables)),
        LTILikelihood();
        noise,
        observables,
    ).logpdf
end

@testset "FVGQ20 joint likelhood derivative in 1st order" begin
    path = joinpath(pkgdir(DifferentiableStateSpaceModels), "test", "data")
    file_prefix = "FVGQ20"
    A = Matrix(DataFrame(CSV.File(joinpath(path, "$(file_prefix)_A.csv"), header = false)))
    B = Matrix(DataFrame(CSV.File(joinpath(path, "$(file_prefix)_B.csv"), header = false)))
    C = Matrix(DataFrame(CSV.File(joinpath(path, "$(file_prefix)_C.csv"), header = false)))
    D_raw =
        Matrix(DataFrame(CSV.File(joinpath(path, "$(file_prefix)_D.csv"), header = false)))
    D = MvNormal(vec(D_raw))
    observables_raw = Matrix(
        DataFrame(
            CSV.File(joinpath(path, "$(file_prefix)_observables.csv"), header = false),
        ),
    )
    noise_raw = Matrix(
        DataFrame(CSV.File(joinpath(path, "$(file_prefix)_noise.csv"), header = false)),
    )
    observables = [observables_raw[i, :] for i = 1:size(observables_raw, 1)]
    noise = [noise_raw[i, :] for i = 1:size(noise_raw, 1)]
    u0 = zeros(size(A, 1))

    res = gradient(minimal_likelihood_test_joint_first, A, B, C, D, u0, noise, observables)

    # Some tests
    @test finite_difference_gradient(
        A -> minimal_likelihood_test_joint_first(A, B, C, D, u0, noise, observables),
        A,
    ) ≈ res[1] rtol = 1E-5
    @test finite_difference_gradient(
        B -> minimal_likelihood_test_joint_first(A, B, C, D, u0, noise, observables),
        B,
    ) ≈ res[2] rtol = 1E-5
    @test finite_difference_gradient(
        C -> minimal_likelihood_test_joint_first(A, B, C, D, u0, noise, observables),
        C,
    ) ≈ res[3] rtol = 1E-5

    # missing test for the D = MvNormal.  No finite_difference_gradient support yet.

    @test finite_difference_gradient(
        u0 -> minimal_likelihood_test_joint_first(A, B, C, D, u0, noise, observables),
        u0,
    ) ≈ res[5] rtol = 1E-7

    noise_grad = finite_difference_gradient(
        noise_mat -> minimal_likelihood_test_joint_first(
            A,
            B,
            C,
            D,
            u0,
            [noise_mat[i, :] for i = 1:size(noise_mat, 1)],
            observables,
        ),
        noise_raw,
    )
    @test [noise_grad[i, :] for i = 1:size(noise_raw, 1)] ≈ res[6] rtol = 1E-7

    observables_grad = finite_difference_gradient(
        observables_mat -> minimal_likelihood_test_joint_first(
            A,
            B,
            C,
            D,
            u0,
            noise,
            [observables_mat[i, :] for i = 1:size(observables_mat, 1)],
        ),
        observables_raw,
    )
    @test [observables_grad[i, :] for i = 1:size(observables_raw, 1)] ≈ res[7] rtol = 1E-7

    # inference
    @inferred solve(
        A,
        B,
        C,
        D,
        u0,
        (0, length(observables)),
        LTILikelihood();
        noise,
        observables,
    )
end
