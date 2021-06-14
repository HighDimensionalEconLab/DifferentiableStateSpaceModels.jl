using DifferentiableStateSpaceModels, ModelingToolkit, SparseArrays, LinearAlgebra,
      Parameters, Test, TimerOutputs, BenchmarkTools

@testset "Dense RBC 2nd Order" begin
    m = @include_example_module(Examples.rbc_observables_benchmark, 2)

    # bookkeeping tests
    @test m.n == 4
    @test m.n_y == 2
    @test m.n_x == 2
    @test m.n_p == 2
    @test m.n_ϵ == 1
    @test m.n_z == 2
    @test m.η == reshape([0; -1], 2, m.n_ϵ)
    @test m.Q == [1.0 0.0 0.0 0.0; 0.0 0.0 1.0 0.0]

    # Solution and Cache test
    p_f = [0.2, 0.02, 0.01, 0.01]
    p = [0.5, 0.95]
    @inferred allocate_cache(m)    
    @inferred generate_perturbation(m, p; p_f)

    cache = allocate_cache(m)
    @inferred generate_perturbation(m, p; p_f, cache)
    sol = generate_perturbation(m, p; p_f, cache)
    @unpack c = get_threadsafe_cache(cache, m, p, p_f)

    # Cache not in solution
    @test c.H_yp ≈ [0.028377570562199098 0.0; 0.0 0.0; 0.0 0.0; 0.0 0.0]
    @test c.H_y ≈ [-0.0283775705621991 0.0; 1.0 -1.0; 0.0 1.0; 0.0 0.0]
    @test c.H_xp ≈ [
        0.00012263591151906127 -0.011623494029190608
        1.0 0.0
        0.0 0.0
        0.0 1.0
    ]
    @test c.H_x ≈ [
        0.0 0.0
        -0.98 0.0
        -0.07263157894736837 -6.884057971014498
        0.0 -0.2
    ]
    @test c.H_yp_p ≈ [[0.011471086498795562 0.0; 0.0 0.0; 0.0 0.0; 0.0 0.0],
           [0.029871126907577997 0.0; 0.0 0.0; 0.0 0.0; 0.0 0.0]]
    @test c.H_y_p ≈
          [[0.0 0.0; 0.0 0.0; 0.0 0.0; 0.0 0.0], [0.0 0.0; 0.0 0.0; 0.0 0.0; 0.0 0.0]]
    @test c.H_xp_p ≈ [[
               0.000473180436623283 -0.06809527035753198
               0.0 0.0
               0.0 0.0
               0.0 0.0
           ], [
               0.00012909043317795924 -0.01223525687283222
               0.0 0.0
               0.0 0.0
               0.0 0.0
           ]]
    @test c.H_x_p ≈ [[
               0.0 0.0
               0.0 0.0
               -0.4255060477077458 -26.561563542978472
               0.0 0.0
           ], [0.0 0.0; 0.0 0.0; 0.0 0.0; 0.0 0.0]]
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
    @test c.Γ_p ≈ [[0.0], [0.0]]
    @test c.H_p ≈ [
        -0.06809527035753199 -0.1773225633743801
        0.0 0.0
        -26.561563542978472 0.0
        0.0 0.0
    ]
    @test c.Ψ_yp[1] ≈ [[
               0.0048317190660755365 0.0 0.0 0.0 6.960218465183541e-6 -0.0006596930439853133 0.0 0.0
               0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
               0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
               0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
               6.960218465183541e-6 0.0 0.0 0.0 6.53894208439611e-7 -2.0658808482697952e-5 0.0 0.0
               -0.0006596930439853133 0.0 0.0 0.0 -2.0658808482697952e-5 0.0019580523687917372 0.0 0.0
               0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
               0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
           ], zeros(8, 8), zeros(8, 8), zeros(8, 8)]

    @test c.Ψ_yp[2] ≈ [zeros(8, 8), zeros(8, 8), zeros(8, 8), zeros(8, 8)]
    @test c.Ψ_y[1] ≈ [[
               0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
               0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
               0.0 0.0 -0.0048317190660755365 0.0 0.0 0.0 0.0 0.0
               0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
               0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
               0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
               0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
               0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
           ], zeros(8, 8), zeros(8, 8), zeros(8, 8)]

    @test c.Ψ_y[2] ≈ [zeros(8, 8), zeros(8, 8), zeros(8, 8), zeros(8, 8)]
    @test c.Ψ_xp[1] ≈ [[
               6.9602184651835404e-6 0.0 0.0 0.0 6.53894208439611e-7 -2.0658808482697952e-5 0.0 0.0
               0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
               0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
               0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
               6.53894208439611e-7 0.0 0.0 0.0 2.0477213369556226e-7 -3.881681383327981e-6 0.0 0.0
               -2.0658808482697952e-5 0.0 0.0 0.0 -3.881681383327981e-6 0.00012263591151906138 0.0 0.0
               0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
               0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
           ], zeros(8, 8), zeros(8, 8), zeros(8, 8)]
    @test c.Ψ_xp[2] ≈ [[
               -0.0006596930439853132 0.0 0.0 0.0 -2.0658808482697952e-5 0.0019580523687917372 0.0 0.0
               0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
               0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
               0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
               -2.0658808482697952e-5 0.0 0.0 0.0 -3.881681383327981e-6 0.00012263591151906138 0.0 0.0
               0.0019580523687917372 0.0 0.0 0.0 0.00012263591151906138 -0.011623494029190612 0.0 0.0
               0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
               0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
           ], zeros(8, 8), zeros(8, 8), zeros(8, 8)]
    @test c.Ψ_x[1] ≈ [zeros(8, 8), zeros(8, 8),
           [
               0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
               0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
               0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
               0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
               0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
               0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
               0.0 0.0 0.0 0.0 0.0 0.0 -2.4255412970806024e-5 0.000766313456772123
               0.0 0.0 0.0 0.0 0.0 0.0 0.000766313456772123 -0.07263157894736838
           ], zeros(8, 8)]
    @test c.Ψ_x[2] ≈ [zeros(8, 8), zeros(8, 8),
           [
               0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
               0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
               0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
               0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
               0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
               0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
               0.0 0.0 0.0 0.0 0.0 0.0 0.000766313456772123 -0.07263157894736838
               0.0 0.0 0.0 0.0 0.0 0.0 -0.07263157894736838 -6.884057971014497
           ], zeros(8, 8)]
    @test c.Ψ_p[1] ≈ [[
               -0.0038647566790457792 0.0 0.0 0.0 -7.971028956261657e-5 0.011471086498795566 0.0 0.0
               0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
               0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
               0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
               -7.971028956261657e-5 0.0 0.0 0.0 -1.238935629209052e-5 0.00047318043662328334 0.0 0.0
               0.011471086498795566 0.0 0.0 0.0 0.00047318043662328334 -0.068095270357532 0.0 0.0
               0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
               0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
           ], zeros(8, 8),
           [
               0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
               0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
               0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
               0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
               0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
               0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
               0.0 0.0 0.0 0.0 0.0 0.0 0.002956756561550659 -0.42550604770774586
               0.0 0.0 0.0 0.0 0.0 0.0 -0.42550604770774586 -26.56156354297847
           ], zeros(8, 8)]

    @test c.Ψ_p[2] ≈ [[
               -0.010063967108852991 0.0 0.0 0.0 -2.1746114192313637e-5 0.0020611077566228815 0.0 0.0
               0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
               0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
               0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
               -2.1746114192313637e-5 0.0 0.0 0.0 -4.085980403503138e-6 0.00012909043317795935 0.0 0.0
               0.0020611077566228815 0.0 0.0 0.0 0.00012909043317795935 -0.012235256872832223 0.0 0.0
               0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
               0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
           ], zeros(8, 8), zeros(8, 8), zeros(8, 8)]
    @test c.y_p ≈ [
        55.78596896689701 76.10141579073955
        66.89124302798608 105.01995379122064
    ]
    @test c.x_p ≈ [555.2637030544529 1445.9269000240533; 0.0 0.0]
    @test c.g_x_p ≈ [[
               -0.12465264193058262 5.596211904442805
               -1.2823781479976832e-15 66.89124302798608
           ],
           [
               -1.6946742377792863 -0.8343618226192915
               -1.1080332409972313 105.01995379122064
           ]]

    # Solution tests
    @test c.h_x_p ≈ [[0.12465264193058134 61.29503112354326; 0.0 0.0],
           [0.586640996782055 105.85431561383992; 0.0 0.0]]
    @test c.Γ ≈ reshape([0.01], 1, 1)
    @test c.Γ_p ≈ [[0.0], [0.0]]
    @test c.Σ ≈ [1e-4]
    @test c.Σ_p ≈ [[0.0], [0.0]]

    @test c.g_xx[:, :, 1] ≈ [
        -0.000371083339499955 0.005130472630563616
        -0.0007663134567721218 0.07263157894736846
    ]
    @test c.g_xx[:, :, 2] ≈
          [0.005130472630563628 0.6265410073784347; 0.07263157894736846 6.8840579710144985]
    @test c.h_xx[:, :, 1] ≈ [-0.0003952301172721667 0.06750110631680481; 0.0 0.0]
    @test c.h_xx[:, :, 2] ≈ [0.06750110631680481 6.257516963636062; 0.0 0.0]
    @test c.g_σσ ≈ [0.0001564980962543059, 0.0]
    @test c.h_σσ ≈ [-0.0001564980962543059, 0.0]

    @test c.g_xx_p[1] ≈ cat([
                  0.005919879353027914 0.0028179768532923975
                  0.010511393863734047 1.0255059778850425e-15
              ],
              [
                  0.0028179768532926456 5.321141188794126
                  1.0255059778850425e-15 66.89124302798595
              ], dims = 3)
    @test c.g_xx_p[2] ≈ cat([
                  0.017520602529482888 -0.1729419353982087
                  0.03507155408568066 -1.1080332409972282
              ],
              [
                  -0.17294193539820824 -0.804951434933394
                  -1.1080332409972282 105.01995379122047
              ], dims = 3)
    @test c.h_xx_p[1] ≈ cat([0.0045915145107061316 -0.0028179768532913723; 0.0 0.0],
              [-0.0028179768532916203 61.570101839191814; 0.0 0.0], dims = 3)
    @test c.h_xx_p[2] ≈ cat([0.017550951556197774 -0.9350913055990194; 0.0 0.0],
              [-0.9350913055990198 105.82490522615385; 0.0 0.0], dims = 3)
    @test c.g_σσ_p ≈ [0.001363945590429837 0.0035331253439556264; 0.0 0.0]
    @test c.h_σσ_p ≈ [-0.001363945590429837 -0.0035331253439556264; 0.0 0.0]
    @test c.Ω_p ≈ [0.0 0.0; 0.0 0.0]
    @test c.B_p ≈ [[0.0; 0.0], [0.0; 0.0]]
    @test c.C_1_p ≈ [[-0.12465264193057919 5.596211904442171; 0.0 0.0], [-1.6946742377792825 -0.8343618226202246; 0.0 0.0]]
    @test c.C_2_p[1] ≈ cat([0.002959939676513957 0.0014089884266461996; 0.0 0.0], [0.001408988426646328 2.660570594397063; 0.0 0.0], dims = 3)
    @test c.C_2_p[2] ≈ cat([0.008760301264741449 -0.08647096769910435; 0.0 0.0], [-0.08647096769910412 -0.402475717466697; 0.0 0.0], dims = 3)
    @test c.C_0_p ≈ [0.0006819727952149194 0.0017665626719778132; 0.0 0.0]

    # Solution
    @test sol.y ≈ [5.936252888048733, 6.884057971014498]
    @test sol.x ≈ [47.39025414828824, 0.0]

    @test sol.A_0 ≈ 0.5 * c.h_σσ
    @test sol.A_1 ≈ [0.9568351489232028 6.209371005755889; -1.5076865909646354e-18 0.2]
    @test sol.A_2 ≈ 0.5 * c.h_xx
    @test sol.B ≈ [0.0; -0.01]
    @test sol.C_0 ≈ [7.824904812715295e-5, 0.0]
    @test sol.C_1 ≈ [0.0957964300241661 0.6746869652585828; 1.0 0.0]
    @test sol.C_2 ≈ cat([-0.0001855416697499775 0.002565236315281808; 0.0 0.0], [0.002565236315281814 0.31327050368921733; 0.0 0.0], dims = 3)
    @test sol.D.σ ≈ [0.01, 0.01]
    @test sol.Γ ≈ [0.01]
    @test sol.g_x ≈ [
        0.0957964300241661 0.6746869652586178
        0.07263157894736878 6.884057971014507
    ]
    @test sol.n == 4
    @test sol.n_y == 2
    @test sol.n_x == 2
    @test sol.n_p == 2
    @test sol.n_ϵ == 1
    @test sol.n_z == 2
    @test sol.η == reshape([0; -1], 2, m.n_ϵ)
    @test sol.Q ≈ [1.0 0.0 0.0 0.0; 0.0 0.0 1.0 0.0]
    @test c.g_x ≈ sol.g_x
    @test c.h_x ≈ sol.A_1
    @test c.Q ≈ sol.Q
    @test c.η ≈ sol.η
    @test c.g_σσ ≈ sol.g_σσ
    @test c.g_xx ≈ sol.g_xx

    @test c.B ≈ sol.B
    @test sol.retcode == :Success
end

@testset "Dense RBC 2nd Order, sigma derivatives" begin
    function rbc_with_sigma()
        H, nt = Examples.rbc()
        p = [nt.p; nt.p_f[3]]
        p_f = [nt.p_f[1]; nt.p_f[2]]
        return H, merge(nt, (; p, p_f)), "rbc_with_sigma"
    end
    H, mod_vals = rbc_with_sigma()
    m = SecondOrderPerturbationModel(H; mod_vals...)
    p_f = [0.2, 0.02]
    p = [0.5, 0.95, 0.01]
    cache = allocate_cache(m)
    sol = generate_perturbation(m, p; p_f, cache)
    @unpack c = get_threadsafe_cache(cache, m, p, p_f)

    @test c.g_σσ_p ≈
          [0.001363945590429837 0.0035331253439556264 0.03129961925086118; 0.0 0.0 0.0]
    @test c.h_σσ_p ≈
          [-0.001363945590429837 -0.0035331253439556264 -0.03129961925086118; 0.0 0.0 0.0]
end

@testset "Dense RBC 2nd Module Generate" begin
    H, mod_vals, model_name = Examples.rbc()
    model_cache_path = save_second_order_module(H; overwrite_model_cache = true,
                                                model_name = model_name * "_second_test",
                                                mod_vals...)
    model_module = include(model_cache_path)
    m = SecondOrderPerturbationModel(model_module)
    # function tests (steady state)
    p_f = [0.2, 0.02, 0.01]
    p = [0.5, 0.95]

    # sol tests
    cache = allocate_cache(m)
    sol = generate_perturbation(m, p; p_f, cache)
    @unpack c = get_threadsafe_cache(cache, m, p, p_f)

    @test c.y ≈ [5.936252888048733, 6.884057971014498]
    @test c.x ≈ [47.39025414828824, 0.0]
    @test c.y_p ≈ [
        55.78596896689701 76.10141579073955
        66.89124302798608 105.01995379122064
    ]
    @test c.x_p ≈ [555.2637030544529 1445.9269000240533; 0.0 0.0]
    @test c.g_x ≈ [
        0.0957964300241661 0.6746869652586178
        0.07263157894736878 6.884057971014507
    ]
    @test c.h_x ≈ [0.9568351489232028 6.209371005755889; -1.5076865909646354e-18 0.2]
    @test c.g_x_p ≈ [[
               -0.12465264193058262 5.596211904442805
               -1.2823781479976832e-15 66.89124302798608
           ],
           [
               -1.6946742377792863 -0.8343618226192915
               -1.1080332409972313 105.01995379122064
           ]]
    @test c.h_x_p ≈ [[0.12465264193058134 61.29503112354326; 0.0 0.0],
           [0.586640996782055 105.85431561383992; 0.0 0.0]]
    @test c.Σ ≈ [1e-4]
    @test c.Σ_p ≈ [[0.0], [0.0]]

    @test c.g_xx[:, :, 1] ≈ [
        -0.000371083339499955 0.005130472630563616
        -0.0007663134567721218 0.07263157894736846
    ]
    @test c.g_xx[:, :, 2] ≈
          [0.005130472630563628 0.6265410073784347; 0.07263157894736846 6.8840579710144985]
    @test c.h_xx[:, :, 1] ≈ [-0.0003952301172721667 0.06750110631680481; 0.0 0.0]
    @test c.h_xx[:, :, 2] ≈ [0.06750110631680481 6.257516963636062; 0.0 0.0]
    @test c.g_σσ ≈ [0.0001564980962543059, 0.0]
    @test c.h_σσ ≈ [-0.0001564980962543059, 0.0]

    @test c.g_xx_p[1] ≈ cat([
                  0.005919879353027914 0.0028179768532923975
                  0.010511393863734047 1.0255059778850425e-15
              ],
              [
                  0.0028179768532926456 5.321141188794126
                  1.0255059778850425e-15 66.89124302798595
              ], dims = 3)
    @test c.g_xx_p[2] ≈ cat([
                  0.017520602529482888 -0.1729419353982087
                  0.03507155408568066 -1.1080332409972282
              ],
              [
                  -0.17294193539820824 -0.804951434933394
                  -1.1080332409972282 105.01995379122047
              ], dims = 3)
    @test c.h_xx_p[1] ≈ cat([0.0045915145107061316 -0.0028179768532913723; 0.0 0.0],
              [-0.0028179768532916203 61.570101839191814; 0.0 0.0], dims = 3)
    @test c.h_xx_p[2] ≈ cat([0.017550951556197774 -0.9350913055990194; 0.0 0.0],
              [-0.9350913055990198 105.82490522615385; 0.0 0.0], dims = 3)
    @test c.g_σσ_p ≈ [0.001363945590429837 0.0035331253439556264; 0.0 0.0]
    @test c.h_σσ_p ≈ [-0.001363945590429837 -0.0035331253439556264; 0.0 0.0]

    @test c.Ω === nothing
    @test c.Ω_p === nothing
    @test sol.n == 4
    @test sol.n_y == 2
    @test sol.n_x == 2
    @test sol.n_p == 2
    @test sol.n_ϵ == 1
    @test sol.n_z == 4
    @test c.η == reshape([0; -1], 2, m.n_ϵ)
    @test c.Q == I
    @test c.y ≈ sol.y
    @test c.x ≈ sol.x
    @test c.g_x ≈ sol.g_x
    @test c.h_x ≈ sol.A_1
    @test c.B ≈ sol.B
    @test c.Ω === sol.D  # should be nothing
    @test c.Q ≈ sol.Q
    @test c.η ≈ sol.η
    @test c.g_σσ ≈ sol.g_σσ
    @test c.h_σσ ≈ sol.A_0 * 2
    @test c.g_xx ≈ sol.g_xx
    @test c.h_xx ≈ sol.A_2 * 2
    @test sol.retcode == :Success
end

@testset "Checks on hashing, recalculations, callbacks" begin
    calculate_steady_state_callback_triggered = false
    function calculate_steady_state_callback(ret, m, c, settings, p, p_f, solver)
        calculate_steady_state_callback_triggered = true
        return nothing
    end
    evaluate_functions_callback_triggered = false
    function evaluate_functions_callback(ret, m, c, settings, p, p_f, solver)
        evaluate_functions_callback_triggered = true
        return nothing
    end
    solve_first_order_callback_triggered = false
    function solve_first_order_callback(ret, m, c, settings)
        solve_first_order_callback_triggered = true
        return nothing
    end
    solve_first_order_p_callback_triggered = false
    function solve_first_order_p_callback(ret, m, c, settings)
        solve_first_order_p_callback_triggered = true
        return nothing
    end
    solve_second_order_callback_triggered = false
    function solve_second_order_callback(ret, m, c, settings)
        solve_second_order_callback_triggered = true
        return nothing
    end
    solve_second_order_p_callback_triggered = false
    function solve_second_order_p_callback(ret, m, c, settings)
        solve_second_order_p_callback_triggered = true
        return nothing
    end

    # note that callbacks
    settings = PerturbationSolverSettings(; calculate_steady_state_callback,
                                                     evaluate_functions_callback,
                                                     solve_first_order_callback,
                                                     solve_first_order_p_callback,
                                                     solve_second_order_callback,
                                                     solve_second_order_p_callback)

    m = @include_example_module(Examples.rbc_observables, 2)

    p = [0.5, 0.95, 0.2, 0.011]
    p_f = [0.02, 0.01, 0.012]
    base_cache = allocate_cache(m)
    base_c = get_threadsafe_cache(base_cache, m, p, p_f).c

    sol = generate_perturbation(m, p; p_f, settings, cache = base_cache)
    @test calculate_steady_state_callback_triggered == true
    @test evaluate_functions_callback_triggered == true
    @test solve_first_order_callback_triggered == true
    @test solve_first_order_p_callback_triggered == true
    @test solve_second_order_callback_triggered == true
    @test solve_second_order_p_callback_triggered == true

    # Resolve only changing the Omega_1.  Shouldn't require recalculation of anything important.

    cache = deepcopy(base_cache)
    @unpack c = get_threadsafe_cache(cache, m, p, p_f)
    p_change_Ω_1 = [0.5, 0.95, 0.2, 0.013]
    @test base_c.p_ss_hash ==
          hash(DifferentiableStateSpaceModels.get_hash_subset(m.select_p_ss_hash,
                                                              p_change_Ω_1))
    @test hash(p) != hash(p_change_Ω_1)
    @test base_c.p_hash != hash(p_change_Ω_1)
    calculate_steady_state_callback_triggered = false
    evaluate_functions_callback_triggered = false
    solve_first_order_callback_triggered = false
    solve_first_order_p_callback_triggered = false
    solve_second_order_callback_triggered = false
    solve_second_order_p_callback_triggered = false

    sol = generate_perturbation(m, p_change_Ω_1; p_f, settings, cache)
    @test calculate_steady_state_callback_triggered == false
    @test evaluate_functions_callback_triggered == true
    @test solve_first_order_callback_triggered == false
    @test solve_first_order_p_callback_triggered == false
    @test solve_second_order_callback_triggered == false
    @test solve_second_order_p_callback_triggered == false

    # Resolve only changing the σ.  No need to recalculate the steady state - flagged to recalculate the perturbation, even if not strictly needed
    cache = deepcopy(base_cache)
    p_f_change_σ = [0.02, 0.015, 0.012]
    @test base_c.p_f_ss_hash ==
          hash(DifferentiableStateSpaceModels.get_hash_subset(m.select_p_f_ss_hash,
                                                              p_f_change_σ))
    @test base_c.p_f_perturbation_hash !=
          hash(DifferentiableStateSpaceModels.get_hash_subset(m.select_p_f_perturbation_hash,
                                                              p_f_change_σ))
    calculate_steady_state_callback_triggered = false
    evaluate_functions_callback_triggered = false
    solve_first_order_callback_triggered = false
    solve_first_order_p_callback_triggered = false
    solve_second_order_callback_triggered = false
    solve_second_order_p_callback_triggered = false

    sol = generate_perturbation(m, p; p_f = p_f_change_σ, settings, cache)
    @test calculate_steady_state_callback_triggered == false
    @test evaluate_functions_callback_triggered == true
    @test solve_first_order_callback_triggered == true
    @test solve_first_order_p_callback_triggered == true
    @test solve_second_order_callback_triggered == true
    @test solve_second_order_p_callback_triggered == true
end

@testset "Caching" begin
m = @include_example_module(Examples.rbc_observables, 2)

settings = PerturbationSolverSettings(;print_level = 3)
p = [0.5, 0.95, 0.2, 0.011]
p_f = [0.02, 0.01, 0.012]
cache = allocate_cache(m)

sol = generate_perturbation(m, p; p_f, settings, cache)
for a in (:p_hash, :p_f_hash)
    @show getfield(cache.caches[1], a)
end

p = [0.51, 0.95, 0.2, 0.011]
p_f = [0.02, 0.01, 0.012]
sol = generate_perturbation(m, p; p_f, settings, cache)
for a in (:p_hash, :p_f_hash)
    @show getfield(cache.caches[1], a)
end
end
