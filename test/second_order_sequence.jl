using DifferentiableStateSpaceModels, DifferenceEquations, LinearAlgebra, Test, Zygote
using DifferentiableStateSpaceModels.Examples
using DelimitedFiles
using FiniteDiff: finite_difference_gradient
using ChainRulesTestUtils

function likelihood_test_joint_second_raw(p_d_input, p_f, ϵ, x0, m, z; kwargs...)
    p_d = (α = p_d_input[1], β = p_d_input[2])
    sol = generate_perturbation(m, p_d, p_f, Val(2))
    problem = QuadraticStateSpaceProblem(sol.A_0, sol.A_1, sol.A_2, sol.B, x0,
                                         (0, size(z, 2)); observables_noise = sol.D,
                                         sol.C_0, sol.C_1,
                                         sol.C_2,
                                         noise = ϵ, observables = z, kwargs...)
    return solve(problem, DirectIteration()).logpdf
end

function likelihood_test_joint_second(p_d_input, p_f, ϵ, x0, m, z; kwargs...)
    p_d = (α = p_d_input[1], β = p_d_input[2])
    sol = generate_perturbation(m, p_d, p_f, Val(2))
    problem = QuadraticStateSpaceProblem(sol, x0,
                                         (0, size(z, 2)); observables = z, noise = ϵ,
                                         kwargs...)
    return solve(problem, DirectIteration()).logpdf
end

# @testset "Sequence likelihood and gradient 2nd order" begin
m = @include_example_module(Examples.rbc_observables)
p_f = (ρ = 0.2, δ = 0.02, σ = 0.01, Ω_1 = 0.1)
p_d = (α = 0.5, β = 0.95)
z = [-0.6949847708598687 -0.8456988740809867;
     -0.7804117657996692 0.07781473603479207;
     -1.1363021614363802 -2.41253450179418;
     -0.2140813001516194 -0.10914617826240575;
     -1.0365874981404577 0.9869373465251516;
     -0.7321498641416826 0.012293325072265942;
     -0.054809260599132194 -1.8233591236618099;
     0.5407452466493482 -0.9773559802938866;
     1.3968232347532277 -2.139194998843768]' |> collect
ϵ = reshape([0.22 0.01 0.14 0.03 0.15 0.21 0.22 0.05 0.18], 1, 9)
x0 = zeros(m.n_x)

c = SolverCache(m, Val(2), p_d)
sol = generate_perturbation(m, p_d, p_f, Val(2); cache = c)

@test likelihood_test_joint_second_raw(p_d, p_f, ϵ, x0, m, z) ≈
      -1077.8016532874813  # Regression test, should verify!  Couldn't find existing number to check against
@test likelihood_test_joint_second(p_d, p_f, ϵ, x0, m, z) ≈
      -1077.8016532874813  # Regression test, should verify!  Couldn't find existing number to 
@inferred likelihood_test_joint_second_raw(p_d, p_f, ϵ, x0, m, z)
@inferred likelihood_test_joint_second(p_d, p_f, ϵ, x0, m, z)

@test likelihood_test_joint_second(p_d, p_f, ϵ, x0, m, z) ≈ -1077.8016532874813 # regression test, couldn't find independent value to compare against

res = gradient((args...) -> likelihood_test_joint_second(args..., m, z), p_d, p_f, ϵ,
               x0)
@test [values(res[1])...] ≈ [305.5874661276336, 559.166700806099]
@test isnothing(res[2])
@test res[3] ≈
      reshape([40.5141940179588 39.32706019833505 25.02785099195666 26.010688843169483 33.01985483763039 31.381238099783715 19.106378855992403 11.441562042277948 -0.9454627257067805],
              1, 9)
@test res[4] ≈ [-562.7501911667255; -4462.244956879686] # Regression test, not sure of the correct number
# CRTU is preferred but FiniteDifferences is having trouble at this point
# end
