### This file is just a convenient file to horserace the speed of gradient evaluation through different backends and models.
### It is orthogonal to the main source or test files.

using DifferentiableStateSpaceModels, Symbolics, LinearAlgebra, Test
using CSV, DataFrames, Zygote
using DifferenceEquations
using DifferentiableStateSpaceModels.Examples

## FVGQ model and preparation
# FVGQ model
isdefined(Main, :FVGQ20) || include(joinpath(pkgdir(DifferentiableStateSpaceModels),
"test/generated_models/FVGQ20.jl"))
m = PerturbationModel(Main.FVGQ20)
p_d = (β = 0.998, h = 0.97, ϑ = 1.17, κ = 9.51, α = 0.21, θp = 0.82, χ = 0.63,
γR = 0.77, γy = 0.19, γΠ = 1.29, Πbar = 1.01, ρd = 0.12, ρφ = 0.93, ρg = 0.95,
g_bar = 0.3, σ_A = exp(-3.97), σ_d = exp(-1.51), σ_φ = exp(-2.36),
σ_μ = exp(-5.43), σ_m = exp(-5.85), σ_g = exp(-3.0), Λμ = 3.4e-3, ΛA = 2.8e-3)
p_f = (δ = 0.025, ε = 10, ϕ = 0, γ2 = 0.001, Ω_ii = sqrt(1e-5))
# FVGQ data
path = joinpath(pkgdir(DifferentiableStateSpaceModels), "test", "data")
file_prefix = "FVGQ20"
A = Matrix(DataFrame(CSV.File(joinpath(path, "$(file_prefix)_A.csv"); header = false)))
B = Matrix(DataFrame(CSV.File(joinpath(path, "$(file_prefix)_B.csv"); header = false)))
C = Matrix(DataFrame(CSV.File(joinpath(path, "$(file_prefix)_C.csv"); header = false)))
observables_raw = Matrix(DataFrame(CSV.File(joinpath(path, "$(file_prefix)_observables.csv"); header = false)))
noise_raw = Matrix(DataFrame(CSV.File(joinpath(path, "$(file_prefix)_noise.csv"); header = false)))
observables = [observables_raw[i, :] for i in 1:size(observables_raw, 1)]
noise = [noise_raw[i, :] for i in 1:size(noise_raw, 1)]
u0 = zeros(size(A, 1))
tspan = (0, length(observables))

# FVGQ first-order model solution
c = SolverCache(m, Val(1), p_d)
sol = generate_perturbation(m, p_d, p_f; cache = c)
generate_perturbation_derivatives!(m, p_d, p_f, c)

function likelihood_test_joint_first(p_d, p_f, noise, u0, m, tspan, observables)
	sol = generate_perturbation(m, p_d, p_f)
	return DifferentiableStateSpaceModels.solve(DifferentiableStateSpaceModels.dssm_evolution,
				DifferentiableStateSpaceModels.dssm_volatility, u0, tspan, sol;
				observables, h = DifferentiableStateSpaceModels.dssm_observation,
				sol.D, noise).logpdf
end

function likelihood_test_joint_first_customAD(p_d, p_f, noise, u0, m, tspan, observables)
	sol = generate_perturbation(m, p_d, p_f)
	return DifferentiableStateSpaceModels.solve(DifferentiableStateSpaceModels.dssm_evolution,
				DifferentiableStateSpaceModels.dssm_volatility, u0, tspan, sol,
                DifferentiableStateSpaceModels.GeneralLikelihood();
				observables, h = DifferentiableStateSpaceModels.dssm_observation,
				sol.D, noise).logpdf
end

function likelihood_test_joint_first_sol(p_d, p_f, noise, u0, m, tspan, observables)
	sol = generate_perturbation(m, p_d, p_f)
	return DifferentiableStateSpaceModels.solve(sol.A, sol.B, sol.C, sol.D, u0, tspan, DifferentiableStateSpaceModels.LTILikelihood(); observables, noise).logpdf
end

function likelihood_test_joint_first_DE(p_d, p_f, noise, u0, m, tspan, observables)
    sol = generate_perturbation(m, p_d, p_f, Val(1))
    problem = LinearStateSpaceProblem(sol.A, sol.B, sol.C, u0, tspan,
                                      noise = noise, obs_noise = sol.D, observables = observables)
    return DifferenceEquations.solve(problem, NoiseConditionalFilter()).loglikelihood
end

f = (p_d, noise) -> likelihood_test_joint_first(p_d, p_f, noise, u0, m, tspan, observables)
res = gradient(f, p_d, noise)

g = (p_d, noise) -> likelihood_test_joint_first_sol(p_d, p_f, noise, u0, m, tspan, observables)
res = gradient(g, p_d, noise)

h = (p_d, noise) -> likelihood_test_joint_first_DE(p_d, p_f, noise, u0, m, tspan, observables)
res = gradient(h, p_d, noise)

q = (p_d, noise) -> likelihood_test_joint_first_customAD(p_d, p_f, noise, u0, m, tspan, observables)
res = gradient(q, p_d, noise)

function likelihood_test_kalman(p_d, p_f, m, tspan, observables)
    sol = generate_perturbation(m, p_d, p_f, Val(1))
    return DifferentiableStateSpaceModels.solve(sol, sol.x_ergodic, tspan; observables).logpdf
end

function likelihood_test_kalman_DE(p_d, p_f, m, tspan, observables)
    sol = generate_perturbation(m, p_d, p_f, Val(1))
    problem = LinearStateSpaceProblem(sol.A, sol.B, sol.C, sol.x_ergodic, tspan,
                                      noise = nothing, obs_noise = sol.D, observables = observables)
    # Solve with Kalman filter
    return DifferenceEquations.solve(problem, KalmanFilter()).loglikelihood
end

f = p_d -> likelihood_test_kalman(p_d, p_f, m, tspan, observables)
res = gradient(f, p_d)

g = p_d -> likelihood_test_kalman_DE(p_d, p_f, m, tspan, observables)
res = gradient(g, p_d)

# FVGQ second-order model solution
c = SolverCache(m, Val(2), p_d)
sol = generate_perturbation(m, p_d, p_f, Val(2); cache = c)
generate_perturbation_derivatives!(m, p_d, p_f, c)

function likelihood_test_joint_second(p_d, p_f, noise, u0, m, tspan, observables)
    sol = generate_perturbation(m, p_d, p_f, Val(2))
    return DifferentiableStateSpaceModels.solve(DifferentiableStateSpaceModels.dssm_evolution,
                 DifferentiableStateSpaceModels.dssm_volatility, [u0; u0], tspan, sol;
                 observables, h = DifferentiableStateSpaceModels.dssm_observation,
                 sol.D, noise).logpdf
end

function likelihood_test_joint_second_customAD(p_d, p_f, noise, u0, m, tspan, observables)
	sol = generate_perturbation(m, p_d, p_f, Val(2))
	return DifferentiableStateSpaceModels.solve(DifferentiableStateSpaceModels.dssm_evolution,
				DifferentiableStateSpaceModels.dssm_volatility, [u0; u0], tspan, sol,
                DifferentiableStateSpaceModels.GeneralLikelihood();
				observables, h = DifferentiableStateSpaceModels.dssm_observation,
				sol.D, noise).logpdf
end

function likelihood_test_joint_second_sol(p_d, p_f, noise, u0, m, tspan, observables)
	sol = generate_perturbation(m, p_d, p_f, Val(2))
	return DifferentiableStateSpaceModels.solve(sol, u0, tspan; observables, noise).logpdf
end

function likelihood_test_joint_second_DE(p_d, p_f, noise, u0, m, tspan, observables)
    sol = generate_perturbation(m, p_d, p_f, Val(2))
    problem = StateSpaceProblem(
        DifferentiableStateSpaceModels.dssm_evolution,
        DifferentiableStateSpaceModels.dssm_volatility,
        DifferentiableStateSpaceModels.dssm_observation,
        [u0; u0],
        tspan,
        sol,
        noise = DefinedNoise(noise),
        obs_noise = sol.D,
        observables = observables
    )
    return DifferenceEquations.solve(problem, NoiseConditionalFilter(); vectype = Zygote.Buffer).loglikelihood
end

f = (p_d, noise) -> likelihood_test_joint_second(p_d, p_f, noise, u0, m, tspan, observables)
res = gradient(f, p_d, noise)

g = (p_d, noise) -> likelihood_test_joint_second_sol(p_d, p_f, noise, u0, m, tspan, observables)
res = gradient(g, p_d, noise)

h = (p_d, noise) -> likelihood_test_joint_second_DE(p_d, p_f, noise, u0, m, tspan, observables)
res = gradient(h, p_d, noise)

q = (p_d, noise) -> likelihood_test_joint_second_customAD(p_d, p_f, noise, u0, m, tspan, observables)
res = gradient(q, p_d, noise)
