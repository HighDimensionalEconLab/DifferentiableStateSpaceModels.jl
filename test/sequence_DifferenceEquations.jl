# Test the downstream `solve` in DifferenceEquations.jl

using DifferentiableStateSpaceModels, DifferenceEquations, LinearAlgebra, Test, Zygote
using Distributions, DistributionsAD
using DifferentiableStateSpaceModels.Examples

@testset "Sequence Simulation, 1st order" begin
    m = @include_example_module(Examples.rbc_observables)
    p_f = (ρ = 0.2, δ = 0.02, σ = 0.01, Ω_1 = 0.1)
    p_d = (α = 0.5, β = 0.95)

    c = SolverCache(m, Val(1), p_d)
    sol = generate_perturbation(m, p_d, p_f; cache = c)

    T = 9
    eps_value = [[0.22], [0.01], [0.14], [0.03], [0.15], [0.21], [0.22], [0.05], [0.18]]
    x0 = zeros(m.n_x)

    # General simulation
    problem = StateSpaceProblem(
        DifferentiableStateSpaceModels.dssm_evolution,
        DifferentiableStateSpaceModels.dssm_volatility,
        DifferentiableStateSpaceModels.dssm_observation,
        x0,
        (0, T),
        sol,
        noise = eps_value
    )    
    simul = DifferenceEquations.solve(problem, NoiseConditionalFilter())
    @test simul.z[2:end] ≈
          [[-0.0014843113235688628, 0.0],
           [-0.001672969226342977, -0.013660616212663025],
           [-0.0025907902434120613, -0.016424018091334355],
           [-0.002808352075911724, -0.02507880927301923],
           [-0.003749820417827369, -0.027731843802628737],
           [-0.00514124775142971, -0.036595970922849726],
           [-0.006607765317710925, -0.05006822400565087],
           [-0.006885970365182856, -0.06457803532896965],
           [-0.00789049635407183, -0.06822941930528102]]
        
    # Linear specific problem
    linear_problem = LinearStateSpaceProblem(
        sol.A,
        sol.B,
        sol.C,
        x0,
        (0, T),
        noise = eps_value
    )
    simul_linear = DifferenceEquations.solve(linear_problem, NoiseConditionalFilter())
    @test simul_linear.z[2:end] ≈
            [[-0.0014843113235688628, 0.0],
             [-0.001672969226342977, -0.013660616212663025],
             [-0.0025907902434120613, -0.016424018091334355],
             [-0.002808352075911724, -0.02507880927301923],
             [-0.003749820417827369, -0.027731843802628737],
             [-0.00514124775142971, -0.036595970922849726],
             [-0.006607765317710925, -0.05006822400565087],
             [-0.006885970365182856, -0.06457803532896965],
             [-0.00789049635407183, -0.06822941930528102]]

    # inference
    @inferred DifferenceEquations.solve(problem, NoiseConditionalFilter())
    @inferred DifferenceEquations.solve(linear_problem, NoiseConditionalFilter())
    @inferred generate_perturbation(m, p_d, p_f; cache = c)
end

function likelihood_test_joint_first_general(p_d_input, p_f, ϵ, x0, m, tspan, z)
    p_d = (α = p_d_input[1], β = p_d_input[2])
    sol = generate_perturbation(m, p_d, p_f, Val(1))
    problem = StateSpaceProblem(
        DifferentiableStateSpaceModels.dssm_evolution,
        DifferentiableStateSpaceModels.dssm_volatility,
        DifferentiableStateSpaceModels.dssm_observation,
        x0,
        tspan,
        sol,
        noise = ϵ,
        obs_noise = sol.D,
        observables = z
    )
    return DifferenceEquations.solve(problem, NoiseConditionalFilter()).loglikelihood
end

function likelihood_test_joint_first_linear(p_d_input, p_f, ϵ, x0, m, tspan, z)
    p_d = (α = p_d_input[1], β = p_d_input[2])
    sol = generate_perturbation(m, p_d, p_f, Val(1))
    problem = LinearStateSpaceProblem(
        sol.A,
        sol.B,
        sol.C,
        x0,
        tspan,
        noise = ϵ,
        obs_noise = sol.D,
        observables = z
    )
    return DifferenceEquations.solve(problem, NoiseConditionalFilter()).loglikelihood
end

@testset "Gradients, generate_perturbation + likelihood, 1st order" begin
    m = @include_example_module(Examples.rbc_observables)
    p_f = (ρ = 0.2, δ = 0.02, σ = 0.01, Ω_1 = 0.1)
    p_d = (α = 0.5, β = 0.95)
    p_d_input = [0.5, 0.95]
    ϵ = [[0.22], [0.01], [0.14], [0.03], [0.15], [0.21], [0.22], [0.05], [0.18]]
    z = [[-0.6949847708598687, -0.8456988740809867],
         [-0.7804117657996692, 0.07781473603479207],
         [-1.1363021614363802, -2.41253450179418],
         [-0.2140813001516194, -0.10914617826240575],
         [-1.0365874981404577, 0.9869373465251516],
         [-0.7321498641416826, 0.012293325072265942],
         [-0.054809260599132194, -1.8233591236618099],
         [0.5407452466493482, -0.9773559802938866],
         [1.3968232347532277, -2.139194998843768]]
    x0 = zeros(m.n_x)
    tspan = (0, length(z))

    res = gradient((p_d_input, ϵ) -> likelihood_test_joint_first_general(p_d_input, p_f, ϵ, x0, m, tspan, z), p_d_input, ϵ)
    @test res[1] ≈ [303.7133186356109, 553.6149537473261]
    @test res[2] ≈ [[40.62454806083384], [39.38899479341156], [25.095297618483304],
                    [26.06697625612332], [33.10959536324157], [31.484308705831474],
                    [19.172319198105615], [11.464791870737214], [-0.9477420442978448]]

    res = gradient((p_d_input, ϵ) -> likelihood_test_joint_first_linear(p_d_input, p_f, ϵ, x0, m, tspan, z), p_d_input, ϵ)
    @test res[1] ≈ [303.7133186356109, 553.6149537473261]
    @test res[2] ≈ [[40.62454806083384], [39.38899479341156], [25.095297618483304],
                    [26.06697625612332], [33.10959536324157], [31.484308705831474],
                    [19.172319198105615], [11.464791870737214], [-0.9477420442978448]]

end

@testset "Sequence Simulation, 2nd order" begin
    m = @include_example_module(Examples.rbc_observables)
    p_f = (ρ = 0.2, δ = 0.02, σ = 0.01, Ω_1 = 0.1)
    p_d = (α = 0.5, β = 0.95)

    c = SolverCache(m, Val(1), p_d)
    sol = generate_perturbation(m, p_d, p_f; cache = c)

    T = 9
    eps_value = [[0.22], [0.01], [0.14], [0.03], [0.15], [0.21], [0.22], [0.05], [0.18]]
    obs_noise = [zeros(2) for _ in 1:T] # there is no observation noises
    x0 = zeros(m.n_x)

    problem = StateSpaceProblem(
        DifferentiableStateSpaceModels.dssm_evolution,
        DifferentiableStateSpaceModels.dssm_volatility,
        DifferentiableStateSpaceModels.dssm_observation,
        x0,
        (0, T),
        sol,
        noise = DefinedNoise(eps_value),
        obs_noise = DefinedNoise(obs_noise)
    )
    
    m = @include_example_module(Examples.rbc_observables)
    p_f = (ρ = 0.2, δ = 0.02, σ = 0.01, Ω_1 = 0.1)
    p_d = (α = 0.5, β = 0.95)

    c = SolverCache(m, Val(2), p_d)
    sol = generate_perturbation(m, p_d, p_f, Val(2); cache = c)

    T = 9
    eps_value = [[0.22], [0.01], [0.14], [0.03], [0.15], [0.21], [0.22], [0.05], [0.18]]
    obs_noise = [zeros(2) for _ in 1:T] # there is no observation noises
    x0 = zeros(m.n_x)
    problem = StateSpaceProblem(
        DifferentiableStateSpaceModels.dssm_evolution,
        DifferentiableStateSpaceModels.dssm_volatility,
        DifferentiableStateSpaceModels.dssm_observation,
        [x0; x0],
        (0, T),
        sol,
        noise = DefinedNoise(eps_value),
        obs_noise = DefinedNoise(obs_noise)
    )
    simul = DifferenceEquations.solve(problem, NoiseConditionalFilter())
    @test simul.z[2:end] ≈ [[-0.0014120420256672264, -7.824904812715083e-5],
                            [-0.001607843339241866, -0.013798593509356017],
                            [-0.0025317633821568975, -0.016632915260522855],
                            [-0.002755836083102567, -0.02534820495624617],
                            [-0.0037026560860619877, -0.028065833613511802],
                            [-0.005097998299327853, -0.036982697654476426],
                            [-0.006567668012302603, -0.05049239849879983],
                            [-0.006851208817373075, -0.06503101829364497],
                            [-0.007859486666066144, -0.06873403795558215]]
    # Compare with old sequence.jl stuff
    @test simul.z ≈
          DifferentiableStateSpaceModels.solve(sol, x0, (0, T), DifferentiableStateSpaceModels.QTI(); noise = eps_value).z
    @inferred DifferenceEquations.solve(problem, NoiseConditionalFilter())
end

function likelihood_test_joint_second(p_d_input, p_f, ϵ, x0, m, tspan, z)
    p_d = (α = p_d_input[1], β = p_d_input[2])
    sol = generate_perturbation(m, p_d, p_f, Val(2))
    problem = StateSpaceProblem(
        DifferentiableStateSpaceModels.dssm_evolution,
        DifferentiableStateSpaceModels.dssm_volatility,
        DifferentiableStateSpaceModels.dssm_observation,
        x0,
        tspan,
        sol,
        noise = DefinedNoise(ϵ),
        obs_noise = sol.D,
        observables = z
    )
    return DifferenceEquations.solve(problem, NoiseConditionalFilter(); vectype = Zygote.Buffer).loglikelihood
end

@testset "Gradients, generate_perturbation + likelihood, 2nd order" begin
    m = @include_example_module(Examples.rbc_observables)
    p_f = (ρ = 0.2, δ = 0.02, σ = 0.01, Ω_1 = 0.1)
    p_d = (α = 0.5, β = 0.95)
    p_d_input = [0.5, 0.95]
    ϵ = [[0.22], [0.01], [0.14], [0.03], [0.15], [0.21], [0.22], [0.05], [0.18]]
    z = [[-0.6949847708598687, -0.8456988740809867],
         [-0.7804117657996692, 0.07781473603479207],
         [-1.1363021614363802, -2.41253450179418],
         [-0.2140813001516194, -0.10914617826240575],
         [-1.0365874981404577, 0.9869373465251516],
         [-0.7321498641416826, 0.012293325072265942],
         [-0.054809260599132194, -1.8233591236618099],
         [0.5407452466493482, -0.9773559802938866],
         [1.3968232347532277, -2.139194998843768]]
    x0 = zeros(2 * m.n_x)
    tspan = (0, length(z))

    res = gradient((p_d_input, ϵ) -> likelihood_test_joint_second(p_d_input, p_f, ϵ, x0, m,
                                                                  tspan, z), p_d_input, ϵ)
    @test res[1] ≈ [305.5874661276336, 559.166700806099]
    @test res[2] ≈ [[40.5141940179588], [39.32706019833505], [25.02785099195666],
                    [26.010688843169483], [33.01985483763039], [31.381238099783715],
                    [19.106378855992403], [11.441562042277948], [-0.9454627257067805]]


    # inferred
    # @inferred likelihood_test_joint_second(p_d_input, p_f, ϵ, x0, m, tspan, z)
    # @inferred likelihood_test_joint_second_sol(p_d_input, p_f, ϵ, x0, m, tspan, z)
end
