using DifferentiableStateSpaceModels,
    ModelingToolkit,
    SparseArrays,
    LinearAlgebra,
    Parameters,
    Test,
    TimerOutputs,
    BenchmarkTools
using DifferentiableStateSpaceModels: all_equal_struct

# @testset "Dense RBC Constructor" begin
#     #p_d = [α, β]
#     # m = @include_example_module(Examples.rbc_observables_benchmark)

#     # # bookkeeping tests
#     # @test m.n == 4
#     # @test m.n_y == 2
#     # @test m.n_x == 2
#     # @test m.n_p == 2
#     # @test m.n_ϵ == 1
#     # @test m.n_z == 2
#     # @test m.η == reshape([0; -1], 2, m.n_ϵ)
#     # @test m.Q == [1.0 0.0 0.0 0.0; 0.0 0.0 1.0 0.0]

#     # # function tests (steady state)
#     # p_f = [0.2, 0.02, 0.01, 0.01]
#     # p = [0.5, 0.95]
#     # c = allocate_cache(m)
#     # @inferred allocate_cache(m)
#     # sol = generate_perturbation(m, p; p_f, cache = c)
#     # @inferred generate_perturbation(m, p; p_f, cache = c)
#     # @test sol.retcode == :Success

#     # y = zeros(m.n_y)
#     # x = zeros(m.n_x)

#     # # m.mod.ȳ!(y, p)
#     # # m.mod.x̄!(x, p)
#     # # @test y ≈ [5.936252888048733, 6.884057971014498]
#     # # @test x ≈ [47.39025414828825, 0.0]

#     # m.mod.H_yp!(c.H_yp, y, x, p)
#     # @test c.H_yp ≈ [0.028377570562199098 0.0; 0.0 0.0; 0.0 0.0; 0.0 0.0]

#     # m.mod.H_y!(c.H_y, y, x, p)
#     # @test c.H_y ≈ [-0.0283775705621991 0.0; 1.0 -1.0; 0.0 1.0; 0.0 0.0]

#     # m.mod.H_xp!(c.H_xp, y, x, p)
#     # @test c.H_xp ≈ [
#     #     0.00012263591151906127 -0.011623494029190608
#     #     1.0 0.0
#     #     0.0 0.0
#     #     0.0 1.0
#     # ]

#     # m.mod.H_x!(c.H_x, y, x, p)
#     # @test c.H_x ≈ [
#     #     0.0 0.0
#     #     -0.98 0.0
#     #     -0.07263157894736837 -6.884057971014498
#     #     0.0 -0.2
#     # ]

#     # m.mod.H_yp_p!(c.H_yp_p, y, x, p)
#     # @test c.H_yp_p ≈ [
#     #     [0.011471086498795562 0.0; 0.0 0.0; 0.0 0.0; 0.0 0.0],
#     #     [0.029871126907577997 0.0; 0.0 0.0; 0.0 0.0; 0.0 0.0],
#     # ]

#     # m.mod.H_y_p!(c.H_y_p, y, x, p)
#     # @test c.H_y_p ≈
#     #       [[0.0 0.0; 0.0 0.0; 0.0 0.0; 0.0 0.0], [0.0 0.0; 0.0 0.0; 0.0 0.0; 0.0 0.0]]

#     # m.mod.H_xp_p!(c.H_xp_p, y, x, p)
#     # @test c.H_xp_p ≈ [
#     #     [
#     #         0.000473180436623283 -0.06809527035753198
#     #         0.0 0.0
#     #         0.0 0.0
#     #         0.0 0.0
#     #     ],
#     #     [
#     #         0.00012909043317795924 -0.01223525687283222
#     #         0.0 0.0
#     #         0.0 0.0
#     #         0.0 0.0
#     #     ],
#     # ]

#     # m.mod.H_x_p!(c.H_x_p, y, x, p)
#     # @test c.H_x_p ≈ [
#     #     [
#     #         0.0 0.0
#     #         0.0 0.0
#     #         -0.4255060477077458 -26.561563542978472
#     #         0.0 0.0
#     #     ],
#     #     [0.0 0.0; 0.0 0.0; 0.0 0.0; 0.0 0.0],
#     # ]

#     # m.mod.Ψ!(c.Ψ, y, x, p)
#     # @test c.Ψ[1] ≈ [
#     #     -0.009560768753410337 0.0 0.0 0.0 -2.0658808482697935e-5 0.0019580523687917364 0.0 0.0
#     #     0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
#     #     0.0 0.0 0.009560768753410338 0.0 0.0 0.0 0.0 0.0
#     #     0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
#     #     -2.0658808482697935e-5 0.0 0.0 0.0 -3.881681383327978e-6 0.00012263591151906127 0.0 0.0
#     #     0.0019580523687917364 0.0 0.0 0.0 0.00012263591151906127 -0.011623494029190608 0.0 0.0
#     #     0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
#     #     0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
#     # ]
#     # @test c.Ψ[2] ≈ zeros(8, 8)
#     # @test c.Ψ[3] ≈ [
#     #     0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
#     #     0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
#     #     0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
#     #     0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
#     #     0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
#     #     0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
#     #     0.0 0.0 0.0 0.0 0.0 0.0 0.0007663134567721225 -0.07263157894736837
#     #     0.0 0.0 0.0 0.0 0.0 0.0 -0.07263157894736837 -6.884057971014498
#     # ]
#     # @test c.Ψ[4] ≈ zeros(8, 8)

#     # m.mod.Γ!(c.Γ, p)
#     # @test c.Γ ≈ [0.01]

#     # m.mod.Γ_p!(c.Γ_p, p)
#     # @test c.Γ_p ≈ [[0.0], [0.0]]

#     # m.mod.H_p!(c.H_p, y, x, p)
#     # @test c.H_p ≈ [
#     #     -0.06809527035753199 -0.1773225633743801
#     #     0.0 0.0
#     #     -26.561563542978472 0.0
#     #     0.0 0.0
#     # ]

#     # solution tests
#     @test c.y ≈ [5.936252888048733, 6.884057971014498]
#     @test c.x ≈ [47.39025414828824, 0.0]
#     @test c.y_p ≈ [
#         55.78596896689701 76.10141579073955
#         66.89124302798608 105.01995379122064
#     ]
#     @test c.x_p ≈ [555.2637030544529 1445.9269000240533; 0.0 0.0]
#     @test c.g_x ≈ [
#         0.0957964300241661 0.6746869652586178
#         0.07263157894736878 6.884057971014507
#     ]
#     @test c.h_x ≈ [0.9568351489232028 6.209371005755889; -1.5076865909646354e-18 0.2]
#     @test c.g_x_p ≈ [
#         [
#             -0.12465264193058262 5.596211904442805
#             -1.2823781479976832e-15 66.89124302798608
#         ],
#         [
#             -1.6946742377792863 -0.8343618226192915
#             -1.1080332409972313 105.01995379122064
#         ],
#     ]
#     @test c.h_x_p ≈ [
#         [0.12465264193058134 61.29503112354326; 0.0 0.0],
#         [0.586640996782055 105.85431561383992; 0.0 0.0],
#     ]
#     @test c.Σ ≈ [1e-4]
#     @test c.Σ_p ≈ [[0.0], [0.0]]

#     @test c.Ω ≈ [0.01, 0.01]
#     @test c.Ω_p ≈ [0.0 0.0; 0.0 0.0]

#     @test c.η == reshape([0; -1], 2, m.n_ϵ)
#     @test c.B ≈ [0.0; -0.01]
#     @test c.B_p ≈ [[0.0; 0.0], [0.0; 0.0]]
#     @test c.Q ≈ [1.0 0.0 0.0 0.0; 0.0 0.0 1.0 0.0]
#     @test c.C_1 ≈ [0.0957964300241661 0.6746869652585828; 1.0 0.0]
#     @test c.C_1_p ≈ [
#         [-0.12465264193057919 5.596211904442171; 0.0 0.0],
#         [-1.6946742377792825 -0.8343618226202246; 0.0 0.0],
#     ]
#     @test Array(c.V) ≈ [
#         0.07005411173180227 0.00015997603451513485
#         0.00015997603451513485 0.00010416666666666667
#     ]
#     @test c.V_p ≈ [
#         [1.584528257999749 0.0015841155991973127; 0.0015841155991973127 0.0],
#         [3.336643488330957 0.002750404724942799; 0.002750404724942799 0.0],
#     ]

#     # Should be identical to c
#     @test sol.n == 4
#     @test sol.n_y == 2
#     @test sol.n_x == 2
#     @test sol.n_p == 2
#     @test sol.n_ϵ == 1
#     @test c.y ≈ sol.y
#     @test c.x ≈ sol.x
#     @test sol.n_z == 2
#     @test c.g_x ≈ sol.g_x
#     @test c.h_x ≈ sol.A
#     @test c.B ≈ sol.B
#     @test c.Ω ≈ sol.D.σ
#     @test c.Q ≈ sol.Q
#     @test c.η ≈ sol.η
#     @test Array(c.V) ≈ Array(sol.x_ergodic.C)  # The `cov` not working on TuringMvNormal quite yet.
#     @test sol.retcode == :Success
# end

# @testset "Dense RBC Module with Module" begin
#     H, mod_vals, model_name = Examples.rbc()
#     model_cache_path =
#         save_first_order_module(H; model_name, overwrite_model_cache = true, mod_vals...)
#     model_module = include(model_cache_path)
#     m = FirstOrderPerturbationModel(model_module)
#     # bookkeeping tests
#     @test m.n == 4
#     @test m.n_y == 2
#     @test m.n_x == 2
#     @test m.n_p == 2
#     @test m.n_ϵ == 1
#     @test m.n_z == 4
#     @test m.η == reshape([0; -1], m.n_x, m.n_ϵ)

#     # function tests (steady state)
#     p_f = [0.2, 0.02, 0.01]
#     p = [0.5, 0.95]
#     y = zeros(m.n_y)
#     x = zeros(m.n_x)
#     @inferred allocate_cache(m)
#     c = allocate_cache(m)

#     sol = generate_perturbation(m, p; p_f, cache = c)

#     m.mod.ȳ!(y, p)
#     m.mod.x̄!(x, p)
#     @test y ≈ [5.936252888048733, 6.884057971014498]
#     @test x ≈ [47.39025414828825, 0.0]

#     m.mod.H_yp!(c.H_yp, y, x, p)
#     @test c.H_yp ≈ [0.028377570562199098 0.0; 0.0 0.0; 0.0 0.0; 0.0 0.0]

#     m.mod.H_y!(c.H_y, y, x, p)
#     @test c.H_y ≈ [-0.0283775705621991 0.0; 1.0 -1.0; 0.0 1.0; 0.0 0.0]

#     m.mod.H_xp!(c.H_xp, y, x, p)
#     @test c.H_xp ≈ [
#         0.00012263591151906127 -0.011623494029190608
#         1.0 0.0
#         0.0 0.0
#         0.0 1.0
#     ]

#     m.mod.H_x!(c.H_x, y, x, p)
#     @test c.H_x ≈ [
#         0.0 0.0
#         -0.98 0.0
#         -0.07263157894736837 -6.884057971014498
#         0.0 -0.2
#     ]

#     m.mod.H_yp_p!(c.H_yp_p, y, x, p)
#     @test c.H_yp_p ≈ [
#         [0.011471086498795562 0.0; 0.0 0.0; 0.0 0.0; 0.0 0.0],
#         [0.029871126907577997 0.0; 0.0 0.0; 0.0 0.0; 0.0 0.0],
#     ]

#     m.mod.H_y_p!(c.H_y_p, y, x, p)
#     @test c.H_y_p ≈
#           [[0.0 0.0; 0.0 0.0; 0.0 0.0; 0.0 0.0], [0.0 0.0; 0.0 0.0; 0.0 0.0; 0.0 0.0]]

#     m.mod.H_xp_p!(c.H_xp_p, y, x, p)
#     @test c.H_xp_p ≈ [
#         [
#             0.000473180436623283 -0.06809527035753198
#             0.0 0.0
#             0.0 0.0
#             0.0 0.0
#         ],
#         [
#             0.00012909043317795924 -0.01223525687283222
#             0.0 0.0
#             0.0 0.0
#             0.0 0.0
#         ],
#     ]

#     m.mod.H_x_p!(c.H_x_p, y, x, p)
#     @test c.H_x_p ≈ [
#         [
#             0.0 0.0
#             0.0 0.0
#             -0.4255060477077458 -26.561563542978472
#             0.0 0.0
#         ],
#         [0.0 0.0; 0.0 0.0; 0.0 0.0; 0.0 0.0],
#     ]

#     m.mod.Ψ!(c.Ψ, y, x, p)
#     @test c.Ψ[1] ≈ [
#         -0.009560768753410337 0.0 0.0 0.0 -2.0658808482697935e-5 0.0019580523687917364 0.0 0.0
#         0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
#         0.0 0.0 0.009560768753410338 0.0 0.0 0.0 0.0 0.0
#         0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
#         -2.0658808482697935e-5 0.0 0.0 0.0 -3.881681383327978e-6 0.00012263591151906127 0.0 0.0
#         0.0019580523687917364 0.0 0.0 0.0 0.00012263591151906127 -0.011623494029190608 0.0 0.0
#         0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
#         0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
#     ]
#     @test c.Ψ[2] ≈ zeros(8, 8)
#     @test c.Ψ[3] ≈ [
#         0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
#         0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
#         0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
#         0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
#         0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
#         0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
#         0.0 0.0 0.0 0.0 0.0 0.0 0.0007663134567721225 -0.07263157894736837
#         0.0 0.0 0.0 0.0 0.0 0.0 -0.07263157894736837 -6.884057971014498
#     ]
#     @test c.Ψ[4] ≈ zeros(8, 8)

#     m.mod.Γ!(c.Γ, p)
#     @test c.Γ ≈ [0.01]

#     m.mod.Γ_p!(c.Γ_p, p)
#     @test c.Γ_p ≈ [[0.0], [0.0]]

#     m.mod.H_p!(c.H_p, y, x, p)
#     @test c.H_p ≈ [
#         -0.06809527035753199 -0.1773225633743801
#         0.0 0.0
#         -26.561563542978472 0.0
#         0.0 0.0
#     ]

#     # solution tests
#     @test sol.retcode == :Success

#     @test sol.y ≈ [5.936252888048733, 6.884057971014498]
#     @test sol.x ≈ [47.39025414828824, 0.0]
#     @test c.y_p ≈ [
#         55.78596896689701 76.10141579073955
#         66.89124302798608 105.01995379122064
#     ]
#     @test c.x_p ≈ [555.2637030544529 1445.9269000240533; 0.0 0.0]
#     @test c.g_x ≈ [
#         0.0957964300241661 0.6746869652586178
#         0.07263157894736878 6.884057971014507
#     ]
#     @test c.h_x ≈ [0.9568351489232028 6.209371005755889; -1.5076865909646354e-18 0.2]
#     @test c.g_x_p ≈ [
#         [
#             -0.12465264193058262 5.596211904442805
#             -1.2823781479976832e-15 66.89124302798608
#         ],
#         [
#             -1.6946742377792863 -0.8343618226192915
#             -1.1080332409972313 105.01995379122064
#         ],
#     ]
#     @test c.h_x_p ≈ [
#         [0.12465264193058134 61.29503112354326; 0.0 0.0],
#         [0.586640996782055 105.85431561383992; 0.0 0.0],
#     ]
#     @test c.Σ ≈ [1e-4]
#     @test c.Σ_p ≈ [[0.0], [0.0]]
# end

# @testset "Dense RBC Constructor with Observables and Module" begin
#     H, mod_vals, model_name = Examples.rbc_observables()
#     model_cache_path =
#         save_first_order_module(H; model_name, overwrite_model_cache = true, mod_vals...)
#     model_module = include(model_cache_path)
#     m = FirstOrderPerturbationModel(model_module)

#     # bookkeeping tests
#     @test m.n == 4
#     @test m.n_y == 2
#     @test m.n_x == 2
#     @test m.n_p == 4
#     @test m.n_ϵ == 1
#     @test m.n_z == 2
#     @test m.η == reshape([0; -1], 2, m.n_ϵ)
#     @test m.Q == [1.0 0.0 0.0 0.0; 0.0 0.0 1.0 0.0]

#     # function tests (steady state)
#     p = [0.5, 0.95, 0.2, 0.011]
#     p_f = [0.02, 0.01, 0.012]

#     # solution tests
#     c = allocate_cache(m)
#     sol = generate_perturbation(m, p; p_f, cache = c)

#     @test c.y ≈ [5.936252888048733, 6.884057971014498]
#     @test c.x ≈ [47.39025414828824, 0.0]
#     @test c.y_p ≈ [
#         55.78596896689701 76.10141579073955 0.0 0.0
#         66.89124302798608 105.01995379122064 0.0 0.0
#     ]
#     @test c.x_p ≈ [555.2637030544529 1445.9269000240533 0.0 0.0; 0.0 0.0 0.0 0.0]
#     @test c.g_x ≈ [
#         0.0957964300241661 0.6746869652586178
#         0.07263157894736878 6.884057971014507
#     ]
#     @test c.h_x ≈ [0.9568351489232028 6.209371005755889; -1.5076865909646354e-18 0.2]
#     @test c.g_x_p ≈ [
#         [
#             -0.12465264193058262 5.596211904442805
#             -1.2823781479976832e-15 66.89124302798608
#         ],
#         [
#             -1.6946742377792863 -0.8343618226192915
#             -1.1080332409972313 105.01995379122064
#         ],
#         [1.3921362894665956e-18 0.29450084693461975; 0.0 0.0],
#         zeros(2, 2),
#     ]
#     @test c.h_x_p ≈ [
#         [0.12465264193058134 61.29503112354326; 0.0 0.0],
#         [0.586640996782055 105.85431561383992; 0.0 0.0],
#         [-1.3921362894665956e-18 -0.29450084693461975; 4.7271045362244914e-18 1.0],
#         zeros(2, 2),
#     ]
#     @test c.Σ ≈ [1e-4]
#     @test c.Σ_p ≈ [[0.0], [0.0], [0.0], [0.0]]

#     @test c.Ω ≈ [0.011, 0.012]
#     @test c.Ω_p ≈ [0.0 0.0 0.0 1.0; 0.0 0.0 0.0 0.0]

#     @test c.η == reshape([0; -1], 2, m.n_ϵ)
#     @test c.Q ≈ [1.0 0.0 0.0 0.0; 0.0 0.0 1.0 0.0]

#     # Should be identical to c
#     @test sol.n == 4
#     @test sol.n_y == 2
#     @test sol.n_x == 2
#     @test sol.n_p == 4
#     @test sol.n_ϵ == 1
#     @test sol.n_z == 2
#     @test c.y ≈ sol.y
#     @test c.x ≈ sol.x
#     @test c.g_x ≈ sol.g_x
#     @test c.h_x ≈ sol.A
#     @test c.B ≈ sol.B
#     @test c.Ω ≈ sol.D.σ
#     @test c.Q ≈ sol.Q
#     @test c.η ≈ sol.η
#     @test sol.retcode == :Success
# end

# @testset "Dense RBC with Timer" begin
#     m = @include_example_module(Examples.rbc)
#     p_f = [0.2, 0.02, 0.01]
#     p = [0.5, 0.95]
#     generate_perturbation(m, p; p_f) # warmup

#     # display timer information
#     reset_timer!()
#     generate_perturbation(m, p; p_f)
#     print_timer()
# end

# @testset "Dense RBC, Empty Theta" begin
#     m = @include_example_module(Examples.rbc_empty_p)

#     # bookkeeping tests
#     @test m.n == 4
#     @test m.n_y == 2
#     @test m.n_x == 2
#     @test m.n_p == 0
#     @test m.n_ϵ == 1
#     @test m.η == reshape([0; -1], 2, m.n_ϵ)

#     # function tests (steady state)
#     p_f = [0.5, 0.95, 0.2, 0.02, 0.01]
#     c = allocate_cache(m)
#     sol = generate_perturbation(m, nothing; p_f, cache = c)
    
#     y = zeros(m.n_y)
#     x = zeros(m.n_x)
#     m.mod.ȳ!(y, nothing, p_f, nothing)
#     m.mod.x̄!(x, nothing, p_f, nothing)
#     @test y ≈ [5.936252888048733, 6.884057971014498]
#     @test x ≈ [47.39025414828825, 0.0]

#     m.mod.H_yp!(c.H_yp, y, x, nothing, p_f, nothing)
#     @test c.H_yp ≈ [0.028377570562199098 0.0; 0.0 0.0; 0.0 0.0; 0.0 0.0]

#     m.mod.H_y!(c.H_y, y, x, nothing, p_f, nothing)
#     @test c.H_y ≈ [-0.0283775705621991 0.0; 1.0 -1.0; 0.0 1.0; 0.0 0.0]

#     m.mod.H_xp!(c.H_xp, y, x, nothing, p_f, nothing)
#     @test c.H_xp ≈ [
#         0.00012263591151906127 -0.011623494029190608
#         1.0 0.0
#         0.0 0.0
#         0.0 1.0
#     ]

#     m.mod.H_x!(c.H_x, y, x, nothing, p_f, nothing)

#     @test c.H_x ≈ [
#         0.0 0.0
#         -0.98 0.0
#         -0.07263157894736837 -6.884057971014498
#         0.0 -0.2
#     ]
#     @test isnothing(m.mod.Ψ!)  #Not generated if theta is empty

#     m.mod.Γ!(c.Γ, nothing, p_f, nothing)
#     @test c.Γ ≈ [0.01]
#     @test isnothing(m.mod.H_p!)

#     # solution tests
#     @test c.y ≈ [5.936252888048733, 6.884057971014498]
#     @test c.x ≈ [47.39025414828824, 0.0]
#     @test c.y_p == Array{Float64}(undef, m.n_y, 0)
#     @test c.x_p == Array{Float64}(undef, m.n_x, 0)
#     @test c.g_x ≈ [
#         0.0957964300241661 0.6746869652586178
#         0.07263157894736878 6.884057971014507
#     ]
#     @test c.h_x ≈ [0.9568351489232028 6.209371005755889; -1.5076865909646354e-18 0.2]
#     @test c.g_x_p == Array{Array{Float64,2},1}(undef, 0)
#     @test c.h_x_p == Array{Array{Float64,2},1}(undef, 0)
#     @test c.Σ ≈ [1e-4]
# end

@testset "Dense RBC, Empty p_f" begin
    m = @include_example_module(Examples.rbc_empty_p_f)    
    
    # bookkeeping tests
    @test m.n == 4
    @test m.n_y == 2
    @test m.n_x == 2
    @test m.n_p == 5
    @test m.n_ϵ == 1
    @test m.η == reshape([0; -1], 2, m.n_ϵ)

    # function tests (steady state)
    p = [0.5, 0.95, 0.2, 0.02, 0.01]
    c = allocate_cache(m)
    sol = generate_perturbation(m, p; cache = c)

    y = zeros(m.n_y)
    x = zeros(m.n_x)
    m.mod.ȳ!(y, p)
    m.mod.x̄!(x, p)
    @test y ≈ [5.936252888048733, 6.884057971014498]
    @test x ≈ [47.39025414828825, 0.0]

    m.mod.H_yp!(c.H_yp, y, x, p)
    @test c.H_yp ≈ [0.028377570562199098 0.0; 0.0 0.0; 0.0 0.0; 0.0 0.0]

    m.mod.H_y!(c.H_y, y, x, p)
    @test c.H_y ≈ [-0.0283775705621991 0.0; 1.0 -1.0; 0.0 1.0; 0.0 0.0]

    m.mod.H_xp!(c.H_xp, y, x, p)
    @test c.H_xp ≈ [
        0.00012263591151906127 -0.011623494029190608
        1.0 0.0
        0.0 0.0
        0.0 1.0
    ]

    m.mod.H_x!(c.H_x, y, x, p)
    @test c.H_x ≈ [
        0.0 0.0
        -0.98 0.0
        -0.07263157894736837 -6.884057971014498
        0.0 -0.2
    ]

    m.mod.Ψ!(c.Ψ, y, x, p)
    @test c.Ψ[1] ≈ [
        -0.009560768753410337 0.0 0.0 0.0 -2.0658808482697935e-5 0.0019580523687917364 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.009560768753410338 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        -2.0658808482697935e-5 0.0 0.0 0.0 -3.881681383327978e-6 0.00012263591151906127 0.0 0.0
        0.0019580523687917364 0.0 0.0 0.0 0.00012263591151906127 -0.011623494029190608 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
    ]
    @test c.Ψ[2] ≈ zeros(8, 8)
    @test c.Ψ[3] ≈ [
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0007663134567721225 -0.07263157894736837
        0.0 0.0 0.0 0.0 0.0 0.0 -0.07263157894736837 -6.884057971014498
    ]
    @test c.Ψ[4] ≈ zeros(8, 8)

    m.mod.Γ!(c.Γ, p)
    @test c.Γ ≈ [0.01]

    m.mod.Γ_p!(c.Γ_p, p)
    @test c.Γ_p ≈ [[0.0], [0.0], [0.0], [0.0], [1.0]]

    m.mod.H_p!(c.H_p, y, x, p)
    @test c.H_p ≈ [
        -0.06809527035753199 -0.1773225633743801 0.0 0.16003361344537806 0.0
        0.0 0.0 0.0 47.39025414828824 0.0
        -26.561563542978472 0.0 0.0 0.0 0.0
        0.0 0.0 -0.0 0.0 0.0
    ]

    # solution tests
    @test c.y ≈ [5.936252888048733, 6.884057971014498]
    @test c.x ≈ [47.39025414828824, 0.0]
    @test c.y_p ≈ [
        55.78596896689701 76.10141579073955 0.0 -116.07178189943077 0.0
        66.89124302798608 105.01995379122064 0.0 -94.78050829657676 0.0
    ]
    @test c.x_p ≈ [
        555.2637030544529 1445.9269000240533 0.0 -1304.9490272717098 0.0
        0.0 0.0 0.0 0.0 0.0
    ]
    @test c.g_x ≈ [
        0.0957964300241661 0.6746869652586178
        0.07263157894736878 6.884057971014507
    ]
    @test c.h_x ≈ [0.9568351489232028 6.209371005755889; -1.5076865909646354e-18 0.2]
    @test c.g_x_p ≈ [
        [
            -0.12465264193058262 5.596211904442805
            -1.2823781479976832e-15 66.89124302798608
        ],
        [
            -1.6946742377792863 -0.8343618226192915
            -1.1080332409972313 105.01995379122064
        ],
        [4.640454298222595e-19 0.2945008469346586; 0.0 0.0],
        [
            0.6277268890968761 -5.0369653468355455
            1.0000000000000024 -94.78050829657676
        ],
        [0.0 0.0; 0.0 0.0],
    ]
    @test c.h_x_p ≈ [
        [0.12465264193058134 61.29503112354326; 0.0 0.0],
        [0.586640996782055 105.85431561383992; 0.0 0.0],
        [
            -4.640454298222595e-19 -0.2945008469346586
            1.5757015120748295e-18 1.0
        ],
        [-0.6277268890968737 -89.74354294974118; 0.0 0.0],
        [0.0 0.0; 0.0 0.0],
    ]
    @test c.Σ ≈ [1e-4]
    @test c.Σ_p ≈ [[0.0], [0.0], [0.0], [0.0], [0.02]]

    # Call again reusing the cache but not the solution
    settings = PerturbationSolverSettings(; use_solution_cache = false)
    c = allocate_cache(m)
    sol = generate_perturbation(m, p; cache = c, settings)

    @test c.y ≈ [5.936252888048733, 6.884057971014498]
    @test c.x ≈ [47.39025414828824, 0.0]
    @test c.y_p ≈ [
        55.78596896689701 76.10141579073955 0.0 -116.07178189943077 0.0
        66.89124302798608 105.01995379122064 0.0 -94.78050829657676 0.0
    ]
    @test c.x_p ≈ [
        555.2637030544529 1445.9269000240533 0.0 -1304.9490272717098 0.0
        0.0 0.0 0.0 0.0 0.0
    ]
    @test c.g_x ≈ [
        0.0957964300241661 0.6746869652586178
        0.07263157894736878 6.884057971014507
    ]
    @test c.h_x ≈ [0.9568351489232028 6.209371005755889; -1.5076865909646354e-18 0.2]
    @test c.g_x_p ≈ [
        [
            -0.12465264193058262 5.596211904442805
            -1.2823781479976832e-15 66.89124302798608
        ],
        [
            -1.6946742377792863 -0.8343618226192915
            -1.1080332409972313 105.01995379122064
        ],
        [4.640454298222595e-19 0.2945008469346586; 0.0 0.0],
        [
            0.6277268890968761 -5.0369653468355455
            1.0000000000000024 -94.78050829657676
        ],
        [0.0 0.0; 0.0 0.0],
    ]
    @test c.h_x_p ≈ [
        [0.12465264193058134 61.29503112354326; 0.0 0.0],
        [0.586640996782055 105.85431561383992; 0.0 0.0],
        [
            -4.640454298222595e-19 -0.2945008469346586
            1.5757015120748295e-18 1.0
        ],
        [-0.6277268890968737 -89.74354294974118; 0.0 0.0],
        [0.0 0.0; 0.0 0.0],
    ]
    @test c.Σ ≈ [1e-4]
    @test c.Σ_p ≈ [[0.0], [0.0], [0.0], [0.0], [0.02]]

    # Call again reusing the cache
    sol = generate_perturbation(m, p; cache = c, settings)

    @test c.y ≈ [5.936252888048733, 6.884057971014498]
    # should  add more...
end

@testset "Dense RBC, Empty p_f, Multiple Shocks" begin
    m = @include_example_module(Examples.rbc_empty_p_f_multiple_shocks)        

    # bookkeeping tests
    @test m.n == 4
    @test m.n_y == 2
    @test m.n_x == 2
    @test m.n_p == 5
    @test m.n_ϵ == 2
    @test m.η == [0 -1; 0 -1]

    # function tests (steady state)
    p = [0.5, 0.95, 0.2, 0.02, 0.01]
    c = allocate_cache(m)
    sol = generate_perturbation(m, p; cache = c)

    y = zeros(m.n_y)
    x = zeros(m.n_x)
    m.mod.ȳ!(y, p)
    m.mod.x̄!(x, p)
    @test y ≈ [5.936252888048733, 6.884057971014498]
    @test x ≈ [47.39025414828825, 0.0]

    m.mod.H_yp!(c.H_yp, y, x, p)
    @test c.H_yp ≈ [0.028377570562199098 0.0; 0.0 0.0; 0.0 0.0; 0.0 0.0]

    m.mod.H_y!(c.H_y, y, x, p)
    @test c.H_y ≈ [-0.0283775705621991 0.0; 1.0 -1.0; 0.0 1.0; 0.0 0.0]

    m.mod.H_xp!(c.H_xp, y, x, p)
    @test c.H_xp ≈ [
        0.00012263591151906127 -0.011623494029190608
        1.0 0.0
        0.0 0.0
        0.0 1.0
    ]

    m.mod.H_x!(c.H_x, y, x, p)
    @test c.H_x ≈ [
        0.0 0.0
        -0.98 0.0
        -0.07263157894736837 -6.884057971014498
        0.0 -0.2
    ]

    m.mod.Ψ!(c.Ψ, y, x, p)
    @test c.Ψ[1] ≈ [
        -0.009560768753410337 0.0 0.0 0.0 -2.0658808482697935e-5 0.0019580523687917364 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.009560768753410338 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        -2.0658808482697935e-5 0.0 0.0 0.0 -3.881681383327978e-6 0.00012263591151906127 0.0 0.0
        0.0019580523687917364 0.0 0.0 0.0 0.00012263591151906127 -0.011623494029190608 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
    ]
    @test c.Ψ[2] ≈ zeros(8, 8)
    @test c.Ψ[3] ≈ [
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0007663134567721225 -0.07263157894736837
        0.0 0.0 0.0 0.0 0.0 0.0 -0.07263157894736837 -6.884057971014498
    ]
    @test c.Ψ[4] ≈ zeros(8, 8)

    m.mod.Γ!(c.Γ, p)
    @test c.Γ ≈ [0.01 0.0; 0.01 0.01]

    m.mod.Γ_p!(c.Γ_p, p)
    @test c.Γ_p ≈ [
        [0.0 0.0; 0.0 0.0],
        [0.0 0.0; 0.0 0.0],
        [0.0 0.0; 0.0 0.0],
        [0.0 0.0; 0.0 0.0],
        [1.0 0.0; 1.0 1.0],
    ]

    # solution tests
    @test c.y ≈ [5.936252888048733, 6.884057971014498]
    @test c.x ≈ [47.39025414828824, 0.0]
    @test c.y_p ≈ [
        55.78596896689701 76.10141579073955 0.0 -116.07178189943077 0.0
        66.89124302798608 105.01995379122064 0.0 -94.78050829657676 0.0
    ]
    @test c.x_p ≈ [
        555.2637030544529 1445.9269000240533 0.0 -1304.9490272717098 0.0
        0.0 0.0 0.0 0.0 0.0
    ]
    @test c.g_x ≈ [
        0.0957964300241661 0.6746869652586178
        0.07263157894736878 6.884057971014507
    ]
    @test c.h_x ≈ [0.9568351489232028 6.209371005755889; -1.5076865909646354e-18 0.2]
    @test c.g_x_p ≈ [
        [
            -0.12465264193058262 5.596211904442805
            -1.2823781479976832e-15 66.89124302798608
        ],
        [
            -1.6946742377792863 -0.8343618226192915
            -1.1080332409972313 105.01995379122064
        ],
        [4.640454298222595e-19 0.2945008469346586; 0.0 0.0],
        [
            0.6277268890968761 -5.0369653468355455
            1.0000000000000024 -94.78050829657676
        ],
        [0.0 0.0; 0.0 0.0],
    ]
    @test c.h_x_p ≈ [
        [0.12465264193058134 61.29503112354326; 0.0 0.0],
        [0.586640996782055 105.85431561383992; 0.0 0.0],
        [
            -4.640454298222595e-19 -0.2945008469346586
            1.5757015120748295e-18 1.0
        ],
        [-0.6277268890968737 -89.74354294974118; 0.0 0.0],
        [0.0 0.0; 0.0 0.0],
    ]
    @test c.Σ ≈ [0.0001 0.0001; 0.0001 0.0002]
    @test c.Σ_p ≈ [
        [0.0 0.0; 0.0 0.0],
        [0.0 0.0; 0.0 0.0],
        [0.0 0.0; 0.0 0.0],
        [0.0 0.0; 0.0 0.0],
        [0.02 0.02; 0.02 0.04],
    ]
end

# @testset "Initial Conditions + Reordering" begin
    # m = @include_example_module(Examples.rbc_solve_steady_state)       
    
    # # bookkeeping tests
    # @test m.n == 4
    # @test m.n_y == 2
    # @test m.n_x == 2
    # @test m.n_p == 2
    # @test m.n_ϵ == 1
    # @test m.η == reshape([0; -1], 2, m.n_ϵ)

    # p_f = [0.2, 0.02, 0.01]
    # p = [0.5, 0.95]
    # c = allocate_cache(m)
    # sol = generate_perturbation(m, p; p_f, cache = c)

    # y = zeros(m.n_y)
    # x = zeros(m.n_x)
    # m.mod.ȳ_iv!(y, p)
    # m.mod.x̄_iv!(x, p)
    # @test y ≈ [5.936252888048733, 6.884057971014498]
    # @test x ≈ [47.39025414828825, 0.0]

#     # solution tests
#     @test c.y ≈ [5.936252888048733, 6.884057971014498]
#     @test c.x ≈ [47.39025414828824, 0.0]
#     @test c.y_p ≈ [
#         55.78596896689701 76.10141579073955
#         66.89124302798608 105.01995379122064
#     ]
#     @test c.x_p ≈ [555.2637030544529 1445.9269000240533; 0.0 0.0]
#     @test c.g_x ≈ [
#         0.0957964300241661 0.6746869652586178
#         0.07263157894736878 6.884057971014507
#     ]
#     @test c.h_x ≈ [0.9568351489232028 6.209371005755889; -1.5076865909646354e-18 0.2]
#     @test c.g_x_p ≈ [
#         [
#             -0.12465264193058262 5.596211904442805
#             -1.2823781479976832e-15 66.89124302798608
#         ],
#         [
#             -1.6946742377792863 -0.8343618226192915
#             -1.1080332409972313 105.01995379122064
#         ],
#     ]
#     @test c.h_x_p ≈ [
#         [0.12465264193058134 61.29503112354326; 0.0 0.0],
#         [0.586640996782055 105.85431561383992; 0.0 0.0],
#     ]
#     @test c.Σ ≈ [1e-4]
#     @test c.Σ_p ≈ [[0.0], [0.0]]
#     @test c.y ≈ sol.y
#     @test c.x ≈ sol.x
#     @test c.g_x ≈ sol.g_x
#     @test c.h_x ≈ sol.A
#     @test c.B ≈ sol.B
#     @test c.Ω === sol.D  # nothing
#     @test c.Q ≈ sol.Q
#     @test c.η ≈ sol.η
# end

# @testset "Initial Conditions Not at Steady State" begin
#     m = @include_example_module(Examples.rbc_solve_steady_state_different_iv)        
    
#     p_f = [0.2, 0.02, 0.01]
#     p = [0.5, 0.95]
#     sol = generate_perturbation(m, p; p_f)
#     @test sol.y ≈ [5.936252888048733, 6.884057971014498]
#     @test sol.x ≈ [47.39025414828825, 0.0]
# end

# @testset "Reordered SS tests" begin
#     m = @include_example_module(Examples.rbc_reorder_ss)      

#     # bookkeeping tests
#     @test m.n == 4
#     @test m.n_y == 2
#     @test m.n_x == 2
#     @test m.n_p == 2
#     @test m.n_ϵ == 1
#     @test m.η == reshape([0; -1], 2, m.n_ϵ)

#     # function tests (steady state)
#     p_f = [0.2, 0.02, 0.01]
#     p = [0.5, 0.95]
#     c = allocate_cache(m)
#     sol = generate_perturbation(m, p; p_f, cache = c)
#     y = zeros(m.n_y)
#     x = zeros(m.n_x)
#     m.mod.ȳ!(y, p)
#     m.mod.x̄!(x, p)
#     @test y ≈ [5.936252888048733, 6.884057971014498]
#     @test x ≈ [47.39025414828825, 0.0]

#     m.mod.H_yp!(c.H_yp, y, x, p)
#     @test c.H_yp ≈ [0.028377570562199098 0.0; 0.0 0.0; 0.0 0.0; 0.0 0.0]

#     m.mod.H_y!(c.H_y, y, x, p)
#     @test c.H_y ≈ [-0.0283775705621991 0.0; 1.0 -1.0; 0.0 1.0; 0.0 0.0]

#     m.mod.H_xp!(c.H_xp, y, x, p)
#     @test c.H_xp ≈ [
#         0.00012263591151906127 -0.011623494029190608
#         1.0 0.0
#         0.0 0.0
#         0.0 1.0
#     ]

#     m.mod.H_x!(c.H_x, y, x, p)
#     @test c.H_x ≈ [
#         0.0 0.0
#         -0.98 0.0
#         -0.07263157894736837 -6.884057971014498
#         0.0 -0.2
#     ]

#     m.mod.H_yp_p!(c.H_yp_p, y, x, p)
#     @test c.H_yp_p ≈ [
#         [0.011471086498795562 0.0; 0.0 0.0; 0.0 0.0; 0.0 0.0],
#         [0.029871126907577997 0.0; 0.0 0.0; 0.0 0.0; 0.0 0.0],
#     ]

#     m.mod.H_y_p!(c.H_y_p, y, x, p)
#     @test c.H_y_p ≈
#           [[0.0 0.0; 0.0 0.0; 0.0 0.0; 0.0 0.0], [0.0 0.0; 0.0 0.0; 0.0 0.0; 0.0 0.0]]

#     m.mod.H_xp_p!(c.H_xp_p, y, x, p)
#     @test c.H_xp_p ≈ [
#         [
#             0.000473180436623283 -0.06809527035753198
#             0.0 0.0
#             0.0 0.0
#             0.0 0.0
#         ],
#         [
#             0.00012909043317795924 -0.01223525687283222
#             0.0 0.0
#             0.0 0.0
#             0.0 0.0
#         ],
#     ]

#     m.mod.H_x_p!(c.H_x_p, y, x, p)
#     @test c.H_x_p ≈ [
#         [
#             0.0 0.0
#             0.0 0.0
#             -0.4255060477077458 -26.561563542978472
#             0.0 0.0
#         ],
#         [0.0 0.0; 0.0 0.0; 0.0 0.0; 0.0 0.0],
#     ]

#     m.mod.Ψ!(c.Ψ, y, x, p)
#     @test c.Ψ[1] ≈ [
#         -0.009560768753410337 0.0 0.0 0.0 -2.0658808482697935e-5 0.0019580523687917364 0.0 0.0
#         0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
#         0.0 0.0 0.009560768753410338 0.0 0.0 0.0 0.0 0.0
#         0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
#         -2.0658808482697935e-5 0.0 0.0 0.0 -3.881681383327978e-6 0.00012263591151906127 0.0 0.0
#         0.0019580523687917364 0.0 0.0 0.0 0.00012263591151906127 -0.011623494029190608 0.0 0.0
#         0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
#         0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
#     ]
#     @test c.Ψ[2] ≈ zeros(8, 8)
#     @test c.Ψ[3] ≈ [
#         0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
#         0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
#         0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
#         0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
#         0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
#         0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
#         0.0 0.0 0.0 0.0 0.0 0.0 0.0007663134567721225 -0.07263157894736837
#         0.0 0.0 0.0 0.0 0.0 0.0 -0.07263157894736837 -6.884057971014498
#     ]
#     @test c.Ψ[4] ≈ zeros(8, 8)

#     m.mod.Γ!(c.Γ, p)
#     @test c.Γ ≈ [0.01]

#     m.mod.Γ_p!(c.Γ_p, p)
#     @test c.Γ_p ≈ [[0.0], [0.0]]

#     # solution tests
#     @test c.y ≈ [5.936252888048733, 6.884057971014498]
#     @test c.x ≈ [47.39025414828824, 0.0]
#     @test c.y_p ≈ [
#         55.78596896689701 76.10141579073955
#         66.89124302798608 105.01995379122064
#     ]
#     @test c.x_p ≈ [555.2637030544529 1445.9269000240533; 0.0 0.0]
#     @test c.g_x ≈ [
#         0.0957964300241661 0.6746869652586178
#         0.07263157894736878 6.884057971014507
#     ]
#     @test c.h_x ≈ [0.9568351489232028 6.209371005755889; -1.5076865909646354e-18 0.2]
#     @test c.g_x_p ≈ [
#         [
#             -0.12465264193058262 5.596211904442805
#             -1.2823781479976832e-15 66.89124302798608
#         ],
#         [
#             -1.6946742377792863 -0.8343618226192915
#             -1.1080332409972313 105.01995379122064
#         ],
#     ]
#     @test c.h_x_p ≈ [
#         [0.12465264193058134 61.29503112354326; 0.0 0.0],
#         [0.586640996782055 105.85431561383992; 0.0 0.0],
#     ]
#     @test c.Σ ≈ [1e-4]
#     @test c.Σ_p ≈ [[0.0], [0.0]]
# end

@testset "Schur Decomposition Failure" begin
    m = @include_example_module(Examples.rbc)
    p_f = [0.2, 0.02, 0.01]
    p = [100.0, 0.0]  # garbage parameters
    settings = PerturbationSolverSettings(; print_level = 0)
    sol = generate_perturbation(m, p; p_f, settings)
    @test sol.retcode == :Failure
end

@testset "BK Condition Failure" begin
    m = @include_example_module(Examples.rbc_empty_p_f)

    p = [0.5, 0.95, 1.01, 0.02, 0.01]  # rho > 1.
    settings = PerturbationSolverSettings(; print_level = 0)
    sol = generate_perturbation(m, p; settings)
    @test sol.retcode == :BlanchardKahnFailure
end

# @testset "Steady State Failure" begin
#     m = @include_example_module(Examples.rbc_solve_steady_state_different_iv)       
#     p_f = [0.2, 0.02, 0.01]
#     p = [0.5, 0.95]
#     settings = PerturbationSolverSettings(; print_level = 0, nlsolve_iterations = 2)  # simulate failure by insufficient iterationrs
#     sol = generate_perturbation(m, p; p_f, settings)
#     @test sol.retcode == :SteadyStateFailure
# end

@testset "Callbacks" begin
    calculate_steady_state_callback_triggered = false
    function calculate_steady_state_callback(ret, m, cache, settings, p, p_f, solver)
        calculate_steady_state_callback_triggered = true

        return nothing
    end
    evaluate_functions_callback_triggered = false
    function evaluate_functions_callback(ret, m, cache, settings, p, p_f, solver)
        evaluate_functions_callback_triggered = true
        return nothing
    end
    solve_first_order_callback_triggered = false
    function solve_first_order_callback(ret, m, cache, settings)
        solve_first_order_callback_triggered = true
        return nothing
    end
    solve_first_order_p_callback_triggered = false
    function solve_first_order_p_callback(ret, m, cache, settings)
        solve_first_order_p_callback_triggered = true
        return nothing
    end

    settings = PerturbationSolverSettings(;
        print_level = 0,
        calculate_steady_state_callback,
        evaluate_functions_callback,
        solve_first_order_callback,
        solve_first_order_p_callback,
    )

    m = @include_example_module(Examples.rbc_observables)
    p = [0.5, 0.95, 0.2, 0.011]
    p_f = [0.02, 0.01, 0.012]
    c = allocate_cache(m)
    sol = generate_perturbation(m, p; p_f, cache = c, settings)
    @test calculate_steady_state_callback_triggered == true
    @test evaluate_functions_callback_triggered == true
    @test solve_first_order_callback_triggered == true
    @test solve_first_order_p_callback_triggered == true
end

# TODO: Not completely confident on test coverage here
# Removing hash checks for now.  Can reenable later

# @testset "Checks on hashing and recalculations" begin
#     calculate_steady_state_callback_triggered = false
#     function calculate_steady_state_callback(ret, m, cache, settings, p, p_f, solver)
#         calculate_steady_state_callback_triggered = true
#         return nothing
#     end
#     evaluate_functions_callback_triggered = false
#     function evaluate_functions_callback(ret, m, cache, settings, p, p_f, solver)
#         evaluate_functions_callback_triggered = true
#         return nothing
#     end
#     solve_first_order_callback_triggered = false
#     function solve_first_order_callback(ret, m, cache, settings)
#         solve_first_order_callback_triggered = true
#         return nothing
#     end
#     solve_first_order_p_callback_triggered = false
#     function solve_first_order_p_callback(ret, m, cache, settings)
#         solve_first_order_p_callback_triggered = true
#         return nothing
#     end

#     settings = PerturbationSolverSettings(;
#         print_level = 0,
#         calculate_steady_state_callback,
#         evaluate_functions_callback,
#         solve_first_order_callback,
#         solve_first_order_p_callback,
#     )

#     m = @include_example_module(Examples.rbc_observables)

#     p = [0.5, 0.95, 0.2, 0.011]
#     p_f = [0.02, 0.01, 0.012]
#     base_cache = allocate_cache(m)
#     sol = generate_perturbation(m, p; p_f, cache = base_cache, settings)
#     @test calculate_steady_state_callback_triggered == true
#     @test evaluate_functions_callback_triggered == true
#     @test solve_first_order_callback_triggered == true
#     @test solve_first_order_p_callback_triggered == true

#     # Resolve only changing the Omega_1.  Shouldn't require recalculation of anything important.
#     cache = deepcopy(base_cache)
#     base_c = deepcopy(base_cache)
#     p_change_Ω_1 = [0.5, 0.95, 0.2, 0.013]

#     @test base_c.p_ss_hash == hash(
#         DifferentiableStateSpaceModels.get_hash_subset(m.select_p_ss_hash, p_change_Ω_1),
#     )
#     @test hash(p) != hash(p_change_Ω_1)
#     @test base_c.p_hash != hash(p_change_Ω_1)

#     # Recalculating
#     calculate_steady_state_callback_triggered = false
#     evaluate_functions_callback_triggered = false
#     solve_first_order_callback_triggered = false
#     solve_first_order_p_callback_triggered = false
#     sol = generate_perturbation(m, p_change_Ω_1; p_f, settings, cache)
#     @test calculate_steady_state_callback_triggered == false
#     @test evaluate_functions_callback_triggered == true
#     @test solve_first_order_callback_triggered == false
#     @test solve_first_order_p_callback_triggered == false

#     # Resolve only changing the σ.  No need to recalculate the steady state - flagged to recalculate the perturbation, even if not strictly needed
#     cache = deepcopy(base_cache)
#     p_f_change_σ = [0.02, 0.015, 0.012]
#     @test base_c.p_f_ss_hash == hash(
#         DifferentiableStateSpaceModels.get_hash_subset(m.select_p_f_ss_hash, p_f_change_σ),
#     )
#     @test base_c.p_f_perturbation_hash != hash(
#         DifferentiableStateSpaceModels.get_hash_subset(
#             m.select_p_f_perturbation_hash,
#             p_f_change_σ,
#         ),
#     )
#     calculate_steady_state_callback_triggered = false
#     evaluate_functions_callback_triggered = false
#     solve_first_order_callback_triggered = false
#     solve_first_order_p_callback_triggered = false
#     sol = generate_perturbation(m, p; p_f = p_f_change_σ, settings, cache)
#     @test calculate_steady_state_callback_triggered == false
#     @test evaluate_functions_callback_triggered == true
#     @test solve_first_order_callback_triggered == true
#     @test solve_first_order_p_callback_triggered == true
# end
