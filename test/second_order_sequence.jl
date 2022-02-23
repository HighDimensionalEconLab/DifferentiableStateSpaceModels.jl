using DifferentiableStateSpaceModels, LinearAlgebra, Test, Zygote
using DifferentiableStateSpaceModels.Examples
using DelimitedFiles
using FiniteDiff: finite_difference_gradient
using ChainRulesTestUtils

function likelihood_test_joint_second(p_d_input, p_f, ϵ, x0, m, z; kwargs...)
    p_d = (α = p_d_input[1], β = p_d_input[2])
    sol = generate_perturbation(m, p_d, p_f, Val(2))
    problem = QuadraticStateSpaceProblem(sol.A_0, sol.A_1, sol.A_2, sol.B, sol.C_0, sol.C_1,
                                         sol.C_2, x0, (0, size(z, 2)); obs_noise = sol.D,
                                         noise = ϵ, observables = z, kwargs...)
    return solve(problem, NoiseConditionalFilter(); save_everystep = false).loglikelihood
end

@testset "Sequence likelihood and gradient 2nd order" begin
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

    problem = QuadraticStateSpaceProblem(sol.A_0, sol.A_1, sol.A_2, sol.B, sol.C_0, sol.C_1,
                                         sol.C_2, x0, (0, size(z, 2)); obs_noise = sol.D,
                                         noise = ϵ, observables = z)
    @test solve(problem, NoiseConditionalFilter(); save_everystep = false).loglikelihood ≈
          -1077.8016532874813  # Regression test, should verify!  Couldn't find existing number to check against
    @inferred solve(problem, NoiseConditionalFilter(); save_everystep = false)

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
end
# @testset "Gradients, generate_perturbation + simulation, 2nd order" begin
#     m = @include_example_module(Examples.rbc_observables)
#     p_f = (ρ=0.2, δ=0.02, σ=0.01, Ω_1=0.1)
#     p_d = (α=0.5, β=0.95)
#     p_d_input = [0.5, 0.95]

#     T = 9
#     ϵ_mat = [0.22, 0.01, 0.14, 0.03, 0.15, 0.21, 0.22, 0.05, 0.18]
#     x0 = zeros(m.n_x)

#     function sum_test_joint_second(p_d_input, ϵ_mat, x0, T, p_f, m; kwargs...)
#         p_d = (α=p_d_input[1], β=p_d_input[2])
#         sol = generate_perturbation(m, p_d, p_f, Val(2); kwargs...)
#         ϵ = map(i -> ϵ_mat[i:i], 1:T)
#         simul = solve(
#             DifferentiableStateSpaceModels.dssm_evolution,
#             DifferentiableStateSpaceModels.dssm_volatility,
#             [x0; x0],
#             (0, T),
#             sol;
#             h = DifferentiableStateSpaceModels.dssm_observation,
#             noise = ϵ,
#         )
#         return sum(sum(simul.z))
#     end

#     settings = PerturbationSolverSettings()
#     # @inferred sum_test_joint_second(p_d_input, ϵ_mat, x0, T, p_f, m; settings)
#     res_zygote = gradient(
#         (p_d_input, ϵ_mat) -> sum_test_joint_second(p_d_input, ϵ_mat, x0, T, p_f, m; settings),
#         p_d_input,
#         ϵ_mat,
#     )
#     res_finite_p = finite_difference_gradient(
#         p_d_input -> sum_test_joint_second(p_d_input, ϵ_mat, x0, T, p_f, m; settings),
#         p_d_input,
#     )
#     res_finite_ϵ = finite_difference_gradient(
#         ϵ_mat -> sum_test_joint_second(p_d_input, ϵ_mat, x0, T, p_f, m; settings),
#         ϵ_mat,
#     )
#     @test res_zygote[2] ≈ res_finite_ϵ
#     @test isapprox(res_zygote[1][1], res_finite_p[1]; rtol = 1e-5)
#     @test isapprox(res_zygote[1][2], res_finite_p[2]; rtol = 1e-5)
# end

# @testset "FVGQ20 joint likelhood derivative in 2nd order" begin
#     path = joinpath(pkgdir(DifferentiableStateSpaceModels), "test", "data")
#     file_prefix = "FVGQ20"
#     A_0_raw = readdlm(joinpath(pkgdir(DifferentiableStateSpaceModels),
#                                "test/data/$(file_prefix)_A_0.csv"))
#     A_0 = vec(A_0_raw)
#     A_1 = readdlm(joinpath(pkgdir(DifferentiableStateSpaceModels),
#                            "test/data/$(file_prefix)_A_1.csv"))
#     A_2_raw = readdlm(joinpath(pkgdir(DifferentiableStateSpaceModels),
#                                "test/data/$(file_prefix)_A_2.csv"))
#     A_2 = reshape(A_2_raw, length(A_0), length(A_0), length(A_0))
#     B = readdlm(joinpath(pkgdir(DifferentiableStateSpaceModels),
#                          "test/data/$(file_prefix)_B.csv"))
#     C_0_raw = readdlm(joinpath(pkgdir(DifferentiableStateSpaceModels),
#                                "test/data/$(file_prefix)_C_0.csv"))
#     C_0 = vec(C_0_raw)
#     C_1 = readdlm(joinpath(pkgdir(DifferentiableStateSpaceModels),
#                            "test/data/$(file_prefix)_C_1.csv"))
#     C_2_raw = readdlm(joinpath(pkgdir(DifferentiableStateSpaceModels),
#                                "test/data/$(file_prefix)_C_2.csv"))
#     C_2 = reshape(C_2_raw, length(C_0), length(A_0), length(A_0))
#     D_raw = readdlm(joinpath(pkgdir(DifferentiableStateSpaceModels),
#                              "test/data/$(file_prefix)_D.csv"))
#     D = MvNormal(Diagonal(map(abs2, vec(D_raw))))
#     observables_raw = Matrix(DataFrame(CSV.File(joinpath(path,
#                                                          "$(file_prefix)_observables.csv");
#                                                 header = false)))
#     noise_raw = readdlm(joinpath(pkgdir(DifferentiableStateSpaceModels),
#                                  "test/data/$(file_prefix)_noise.csv"))
#     observables = [observables_raw[i, :] for i in 1:size(observables_raw, 1)]
#     noise = [noise_raw[i, :] for i in 1:size(noise_raw, 1)]
#     u0 = zeros(length(A_0))

#     res = gradient(minimal_likelihood_test_joint_second, A_0, A_1, A_2, B, C_0, C_1, C_2, D,
#                    u0, noise, observables)

#     # Some tests
#     @test finite_difference_gradient(A_0 -> minimal_likelihood_test_joint_second(A_0, A_1,
#                                                                                  A_2, B,
#                                                                                  C_0, C_1,
#                                                                                  C_2, D, u0,
#                                                                                  noise,
#                                                                                  observables),
#                                      A_0) ≈ res[1] rtol = 1E-5
#     @test finite_difference_gradient(A_1 -> minimal_likelihood_test_joint_second(A_0, A_1,
#                                                                                  A_2, B,
#                                                                                  C_0, C_1,
#                                                                                  C_2, D, u0,
#                                                                                  noise,
#                                                                                  observables),
#                                      A_1) ≈ res[2] rtol = 1E-5
#     # I didn't add the tests for A_2 and C_2 because of their dimensions
#     @test finite_difference_gradient(B -> minimal_likelihood_test_joint_second(A_0, A_1,
#                                                                                A_2, B, C_0,
#                                                                                C_1, C_2, D,
#                                                                                u0, noise,
#                                                                                observables),
#                                      B) ≈ res[4] rtol = 1E-5
#     @test finite_difference_gradient(C_0 -> minimal_likelihood_test_joint_second(A_0, A_1,
#                                                                                  A_2, B,
#                                                                                  C_0, C_1,
#                                                                                  C_2, D, u0,
#                                                                                  noise,
#                                                                                  observables),
#                                      C_0) ≈ res[5] rtol = 1E-5
#     @test finite_difference_gradient(C_1 -> minimal_likelihood_test_joint_second(A_0, A_1,
#                                                                                  A_2, B,
#                                                                                  C_0, C_1,
#                                                                                  C_2, D, u0,
#                                                                                  noise,
#                                                                                  observables),
#                                      C_1) ≈ res[6] rtol = 1E-5
#     # missing test for the D = MvNormal.  No finite_difference_gradient support yet.

#     @test finite_difference_gradient(u0 -> minimal_likelihood_test_joint_second(A_0, A_1,
#                                                                                 A_2, B, C_0,
#                                                                                 C_1, C_2, D,
#                                                                                 u0, noise,
#                                                                                 observables),
#                                      u0) ≈ res[9] rtol = 1E-7

#     noise_grad = finite_difference_gradient(noise_mat -> minimal_likelihood_test_joint_second(A_0,
#                                                                                               A_1,
#                                                                                               A_2,
#                                                                                               B,
#                                                                                               C_0,
#                                                                                               C_1,
#                                                                                               C_2,
#                                                                                               D,
#                                                                                               u0,
#                                                                                               [noise_mat[i,
#                                                                                                          :]
#                                                                                                for i in
#                                                                                                    1:size(noise_mat,
#                                                                                                           1)],
#                                                                                               observables),
#                                             noise_raw)
#     @test [noise_grad[i, :] for i in 1:size(noise_raw, 1)] ≈ res[10] rtol = 1E-7

#     observables_grad = finite_difference_gradient(observables_mat -> minimal_likelihood_test_joint_second(A_0,
#                                                                                                           A_1,
#                                                                                                           A_2,
#                                                                                                           B,
#                                                                                                           C_0,
#                                                                                                           C_1,
#                                                                                                           C_2,
#                                                                                                           D,
#                                                                                                           u0,
#                                                                                                           noise,
#                                                                                                           [observables_mat[i,
#                                                                                                                            :]
#                                                                                                            for i in
#                                                                                                                1:size(observables_mat,
#                                                                                                                       1)]),
#                                                   observables_raw)
#     @test [observables_grad[i, :] for i in 1:size(observables_raw, 1)] ≈ res[11] rtol = 1E-7

#     # inference
#     @inferred solve(A_0, A_1, A_2, B, C_0, C_1, C_2, D, u0, (0, length(observables)),
#                     DifferentiableStateSpaceModels.QTILikelihood(); noise, observables)
# end
