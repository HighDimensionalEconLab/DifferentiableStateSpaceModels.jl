using DifferentiableStateSpaceModels, DifferenceEquations, LinearAlgebra, Test, Zygote
using Distributions, DelimitedFiles
using DifferentiableStateSpaceModels.Examples
using FiniteDiff: finite_difference_gradient
using ChainRulesTestUtils

# Raw usage without helper for DSSM solution objects
function kalman_test_raw(p_d_input, p_f, m, z)
    p_d = (α = p_d_input[1], β = p_d_input[2])
    sol = generate_perturbation(m, p_d, p_f)
    linear_problem = LinearStateSpaceProblem(sol.A, sol.B, sol.x_ergodic, (0, size(z, 2));
                                             sol.C, observables_noise = sol.D,
                                             u0_prior = sol.x_ergodic, noise = nothing,
                                             observables = z)
    return solve(linear_problem, KalmanFilter()).logpdf
end

function kalman_test(p_d_input, p_f, m, z)
    p_d = (α = p_d_input[1], β = p_d_input[2])
    sol = generate_perturbation(m, p_d, p_f)
    linear_problem = LinearStateSpaceProblem(sol, sol.x_ergodic, (0, size(z, 2));
                                             observables = z)
    return solve(linear_problem, KalmanFilter()).logpdf
end

# Can't have inside of the testset until #117 fixed
#@testset "Kalman filter and its gradient" begin
const m = @include_example_module(Examples.rbc_observables) # const required due to #117 bug

p_f = (ρ = 0.2, δ = 0.02, σ = 0.01, Ω_1 = 0.1)
p_d = (α = 0.5, β = 0.95)
p_d_input = [0.5, 0.95]
# Data
z = [-0.6949847708598687 -0.8456988740809867;
     -0.7804117657996692 0.07781473603479207;
     -1.1363021614363802 -2.41253450179418;
     -0.2140813001516194 -0.10914617826240575;
     -1.0365874981404577 0.9869373465251516;
     -0.7321498641416826 0.012293325072265942;
     -0.054809260599132194 -1.8233591236618099;
     0.5407452466493482 -0.9773559802938866;
     1.3968232347532277 -2.139194998843768;
     1.3427856576886874 -0.3100476471887863]' |> collect

@test kalman_test_raw(p_d, p_f, m, z) ≈ -805.7558351251781
@test kalman_test(p_d, p_f, m, z) ≈ -805.7558351251781

# Not really needed if test_rrule working
# res = gradient(p_d_input -> kalman_test(p_d_input, p_f, m, z), p_d_input)
# @test res[1] ≈ [1581.9868695425796, 2731.7261992816493]
test_rrule(Zygote.ZygoteRuleConfig(), p_d_input -> kalman_test(p_d_input, p_f, m, z),
           p_d_input; rrule_f = rrule_via_ad, check_inferred = false)
# end

function likelihood_test_joint_first_raw(p_d_input, p_f, ϵ, x0, m, z)
    p_d = (α = p_d_input[1], β = p_d_input[2])
    sol = generate_perturbation(m, p_d, p_f, Val(1))
    problem = LinearStateSpaceProblem(sol.A, sol.B, x0, (0, size(z, 2)); sol.C, noise = ϵ,
                                      observables_noise = sol.D, observables = z)
    return solve(problem, DirectIteration()).logpdf
end

function likelihood_test_joint_first(p_d_input, p_f, ϵ, x0, m, z)
    p_d = (α = p_d_input[1], β = p_d_input[2])
    sol = generate_perturbation(m, p_d, p_f, Val(1))
    problem = LinearStateSpaceProblem(sol, x0, (0, size(z, 2)); noise = ϵ, observables = z)
    return solve(problem, DirectIteration()).logpdf
end

# Can't put in testset until #117 fixed
# @testset "Gradients, generate_perturbation + likelihood, 1st order" begin
# m = @include_example_module(Examples.rbc_observables) # put back in after #117 is fixed
p_f = (ρ = 0.2, δ = 0.02, σ = 0.01, Ω_1 = 0.1)
p_d = (α = 0.5, β = 0.95)
p_d_input = [0.5, 0.95]
ϵ = reshape([0.22 0.01 0.14 0.03 0.15 0.21 0.22 0.05 0.18], 1, 9)
z = [-0.6949847708598687 -0.8456988740809867;
     -0.7804117657996692 0.07781473603479207;
     -1.1363021614363802 -2.41253450179418;
     -0.2140813001516194 -0.10914617826240575;
     -1.0365874981404577 0.9869373465251516;
     -0.7321498641416826 0.012293325072265942;
     -0.054809260599132194 -1.8233591236618099;
     0.5407452466493482 -0.9773559802938866;
     1.3968232347532277 -2.139194998843768]' |> collect
x0 = zeros(m.n_x)

res = gradient((p_d_input, ϵ) -> likelihood_test_joint_first_raw(p_d_input, p_f, ϵ, x0, m,
                                                                 z), p_d_input, ϵ)
@test res[1] ≈ [303.7133186356109, 553.6149537473261]
@test res[2] ≈
      [40.62454806083384 39.38899479341156 25.095297618483304 26.06697625612332 33.10959536324157 31.484308705831474 19.172319198105615 11.464791870737214 -0.9477420442978448]
res = gradient((p_d_input, ϵ) -> likelihood_test_joint_first(p_d_input, p_f, ϵ, x0, m,
                                                             z), p_d_input, ϵ)
@test res[1] ≈ [303.7133186356109, 553.6149537473261]
@test res[2] ≈
      [40.62454806083384 39.38899479341156 25.095297618483304 26.06697625612332 33.10959536324157 31.484308705831474 19.172319198105615 11.464791870737214 -0.9477420442978448]
test_rrule(Zygote.ZygoteRuleConfig(),
           (p_d_input, ϵ) -> likelihood_test_joint_first(p_d_input, p_f, ϵ, x0, m,
                                                         z), p_d_input, ϵ;
           rrule_f = rrule_via_ad, check_inferred = false)
# end