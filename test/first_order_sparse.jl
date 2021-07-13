@testset "Sparse RBC " begin
    H, mod_vals = Examples.rbc_observables_benchmark()
    m = FirstOrderPerturbationModel(H; functions_type = SparseFunctions(), mod_vals...)

    # bookkeeping tests
    @test m.n == 4
    @test m.n_y == 2
    @test m.n_x == 2
    @test m.n_p == 2
    @test m.n_ϵ == 1
    @test m.n_z == 2
    @test m.η == reshape([0; -1], 2, m.n_ϵ)
    @test m.Q == [1.0 0.0 0.0 0.0; 0.0 0.0 1.0 0.0]

    # function tests (steady state)
    p_f = [0.2, 0.02, 0.01, 0.01]
    p = [0.5, 0.95]
    cache = allocate_cache(m)
    sol = generate_perturbation(m, p; p_f, cache)
    @unpack c = get_threadsafe_cache(cache, m, p, p_f)

    y = zeros(m.n_y)
    x = zeros(m.n_x)

    m.ȳ!(y, p, p_f, nothing)
    m.x̄!(x, p, p_f, nothing)
    @test y ≈ [5.936252888048733, 6.884057971014498]
    @test x ≈ [47.39025414828825, 0.0]

    m.H_yp!(c.H_yp, y, x, p, p_f, nothing)
    @test c.H_yp ≈ [0.028377570562199098 0.0; 0.0 0.0; 0.0 0.0; 0.0 0.0]

    m.H_y!(c.H_y, y, x, p, p_f, nothing)
    @test c.H_y ≈ [-0.0283775705621991 0.0; 1.0 -1.0; 0.0 1.0; 0.0 0.0]

    m.H_xp!(c.H_xp, y, x, p, p_f, nothing)
    @test c.H_xp ≈ [
        0.00012263591151906127 -0.011623494029190608
        1.0 0.0
        0.0 0.0
        0.0 1.0
    ]

    m.H_x!(c.H_x, y, x, p, p_f, nothing)
    @test c.H_x ≈ [
        0.0 0.0
        -0.98 0.0
        -0.07263157894736837 -6.884057971014498
        0.0 -0.2
    ]

    m.H_yp_p!(c.H_yp_p, y, x, p, p_f, nothing)
    @test c.H_yp_p ≈ [
        [0.011471086498795562 0.0; 0.0 0.0; 0.0 0.0; 0.0 0.0],
        [0.029871126907577997 0.0; 0.0 0.0; 0.0 0.0; 0.0 0.0],
    ]

    m.H_y_p!(c.H_y_p, y, x, p, p_f, nothing)
    @test c.H_y_p ≈
          [[0.0 0.0; 0.0 0.0; 0.0 0.0; 0.0 0.0], [0.0 0.0; 0.0 0.0; 0.0 0.0; 0.0 0.0]]

    m.H_xp_p!(c.H_xp_p, y, x, p, p_f, nothing)
    @test c.H_xp_p ≈ [
        [
            0.000473180436623283 -0.06809527035753198
            0.0 0.0
            0.0 0.0
            0.0 0.0
        ],
        [
            0.00012909043317795924 -0.01223525687283222
            0.0 0.0
            0.0 0.0
            0.0 0.0
        ],
    ]

    m.H_x_p!(c.H_x_p, y, x, p, p_f, nothing)
    @test c.H_x_p ≈ [
        [
            0.0 0.0
            0.0 0.0
            -0.4255060477077458 -26.561563542978472
            0.0 0.0
        ],
        [0.0 0.0; 0.0 0.0; 0.0 0.0; 0.0 0.0],
    ]

    m.Ψ!(c.Ψ, y, x, p, p_f, nothing)
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

    m.Γ!(c.Γ, p, p_f, nothing)
    @test c.Γ ≈ [0.01]

    m.Γ_p!(c.Γ_p, p, p_f, nothing)
    @test c.Γ_p ≈ [[0.0], [0.0]]

    m.H_p!(c.H_p, y, x, p, p_f, nothing)
    @test c.H_p ≈ [
        -0.06809527035753199 -0.1773225633743801
        0.0 0.0
        -26.561563542978472 0.0
        0.0 0.0
    ]

    # solution tests
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
    @test c.g_x_p ≈ [
        [
            -0.12465264193058262 5.596211904442805
            -1.2823781479976832e-15 66.89124302798608
        ],
        [
            -1.6946742377792863 -0.8343618226192915
            -1.1080332409972313 105.01995379122064
        ],
    ]
    @test c.h_x_p ≈ [
        [0.12465264193058134 61.29503112354326; 0.0 0.0],
        [0.586640996782055 105.85431561383992; 0.0 0.0],
    ]
    @test c.Σ ≈ [1e-4]
    @test c.Σ_p ≈ [[0.0], [0.0]]

    @test c.Ω ≈ [0.01, 0.01]
    @test c.Ω_p ≈ [0.0 0.0; 0.0 0.0]

    @test c.η == reshape([0; -1], 2, m.n_ϵ)
    @test c.B ≈ [0.0; -0.01]
    @test c.B_p ≈ [[0.0; 0.0], [0.0; 0.0]]
    @test c.Q ≈ [1.0 0.0 0.0 0.0; 0.0 0.0 1.0 0.0]
    @test c.C_1 ≈ [0.0957964300241661 0.6746869652585828; 1.0 0.0]
    @test c.C_1_p ≈ [
        [-0.12465264193057919 5.596211904442171; 0.0 0.0],
        [-1.6946742377792825 -0.8343618226202246; 0.0 0.0],
    ]
    @test Array(c.V) ≈ [
        0.07005411173180227 0.00015997603451513485
        0.00015997603451513485 0.00010416666666666667
    ]
    @test c.V_p ≈ [
        [1.584528257999749 0.0015841155991973127; 0.0015841155991973127 0.0],
        [3.336643488330957 0.002750404724942799; 0.002750404724942799 0.0],
    ]

    # Should be identical to c
    @test sol.n == 4
    @test sol.n_y == 2
    @test sol.n_x == 2
    @test sol.n_p == 2
    @test sol.n_ϵ == 1
    @test sol.n_z == 2
    @test c.y ≈ sol.y
    @test c.x ≈ sol.x
    @test c.g_x ≈ sol.g_x
    @test c.h_x ≈ sol.A
    @test c.B ≈ sol.B
    @test c.Ω ≈ sol.D.σ  # should be nothing
    @test c.Q ≈ sol.Q
    @test c.η ≈ sol.η
    @test Array(c.V) ≈ Array(sol.x_ergodic.C)
    @test sol.retcode == :Success
end

@testset "Sparse RBC Module with Module" begin
    H, mod_vals, model_name = Examples.rbc()
    model_cache_path = save_first_order_module(
        H;
        functions_type = SparseFunctions(),
        model_name = "sparse_" * model_name,
        overwrite_model_cache = true,
        mod_vals...,
    )
    model_module = include(model_cache_path)
    m = FirstOrderPerturbationModel(model_module)
    cache = allocate_cache(m)
    p_f = [0.2, 0.02, 0.01]
    p = [0.5, 0.95]
    @unpack c = get_threadsafe_cache(cache, m, p, p_f)

    # bookkeeping tests
    @test m.n == 4
    @test m.n_y == 2
    @test m.n_x == 2
    @test m.n_p == 2
    @test m.n_ϵ == 1
    @test m.η == reshape([0; -1], m.n_x, m.n_ϵ)
    @test m.n_z == 4

    # function tests (steady state)

    y = zeros(m.n_y)
    x = zeros(m.n_x)

    m.ȳ!(y, p, p_f, nothing)
    m.x̄!(x, p, p_f, nothing)
    @test y ≈ [5.936252888048733, 6.884057971014498]
    @test x ≈ [47.39025414828825, 0.0]

    m.H_yp!(c.H_yp, y, x, p, p_f, nothing)
    @test c.H_yp ≈ [0.028377570562199098 0.0; 0.0 0.0; 0.0 0.0; 0.0 0.0]

    m.H_y!(c.H_y, y, x, p, p_f, nothing)
    @test c.H_y ≈ [-0.0283775705621991 0.0; 1.0 -1.0; 0.0 1.0; 0.0 0.0]

    m.H_xp!(c.H_xp, y, x, p, p_f, nothing)
    @test c.H_xp ≈ [
        0.00012263591151906127 -0.011623494029190608
        1.0 0.0
        0.0 0.0
        0.0 1.0
    ]

    m.H_x!(c.H_x, y, x, p, p_f, nothing)
    @test c.H_x ≈ [
        0.0 0.0
        -0.98 0.0
        -0.07263157894736837 -6.884057971014498
        0.0 -0.2
    ]

    m.H_yp_p!(c.H_yp_p, y, x, p, p_f, nothing)
    @test c.H_yp_p ≈ [
        [0.011471086498795562 0.0; 0.0 0.0; 0.0 0.0; 0.0 0.0],
        [0.029871126907577997 0.0; 0.0 0.0; 0.0 0.0; 0.0 0.0],
    ]

    m.H_y_p!(c.H_y_p, y, x, p, p_f, nothing)
    @test c.H_y_p ≈
          [[0.0 0.0; 0.0 0.0; 0.0 0.0; 0.0 0.0], [0.0 0.0; 0.0 0.0; 0.0 0.0; 0.0 0.0]]

    m.H_xp_p!(c.H_xp_p, y, x, p, p_f, nothing)
    @test c.H_xp_p ≈ [
        [
            0.000473180436623283 -0.06809527035753198
            0.0 0.0
            0.0 0.0
            0.0 0.0
        ],
        [
            0.00012909043317795924 -0.01223525687283222
            0.0 0.0
            0.0 0.0
            0.0 0.0
        ],
    ]

    m.H_x_p!(c.H_x_p, y, x, p, p_f, nothing)
    @test c.H_x_p ≈ [
        [
            0.0 0.0
            0.0 0.0
            -0.4255060477077458 -26.561563542978472
            0.0 0.0
        ],
        [0.0 0.0; 0.0 0.0; 0.0 0.0; 0.0 0.0],
    ]

    m.Ψ!(c.Ψ, y, x, p, p_f, nothing)
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

    m.Γ!(c.Γ, p, p_f, nothing)
    @test c.Γ ≈ [0.01]

    m.Γ_p!(c.Γ_p, p, p_f, nothing)
    @test c.Γ_p ≈ [[0.0], [0.0]]

    m.H_p!(c.H_p, y, x, p, p_f, nothing)
    @test c.H_p ≈ [
        -0.06809527035753199 -0.1773225633743801
        0.0 0.0
        -26.561563542978472 0.0
        0.0 0.0
    ]
    #
    # solution tests
    @inferred allocate_cache(m)
    cache = allocate_cache(m)
    @inferred generate_perturbation(m, p; p_f)
    @inferred generate_perturbation(m, p; p_f, cache)
    @unpack c = get_threadsafe_cache(cache, m, p, p_f)
    sol = generate_perturbation(m, p; p_f, cache)

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
    @test c.g_x_p ≈ [
        [
            -0.12465264193058262 5.596211904442805
            -1.2823781479976832e-15 66.89124302798608
        ],
        [
            -1.6946742377792863 -0.8343618226192915
            -1.1080332409972313 105.01995379122064
        ],
    ]
    @test c.h_x_p ≈ [
        [0.12465264193058134 61.29503112354326; 0.0 0.0],
        [0.586640996782055 105.85431561383992; 0.0 0.0],
    ]
    @test c.Σ ≈ [1e-4]
    @test c.Σ_p ≈ [[0.0], [0.0]]
    @test isnothing(c.Ω)
    @test isnothing(c.Ω_p)

    @test c.η == reshape([0; -1], 2, m.n_ϵ)
    @test c.Q == I

    # Should be identical to c
    @test sol.n == 4
    @test sol.n_y == 2
    @test sol.n_x == 2
    @test sol.n_p == 2
    @test sol.n_ϵ == 1
    @test sol.n_z == 4
    @test c.y ≈ sol.y
    @test c.x ≈ sol.x
    @test c.g_x ≈ sol.g_x
    @test c.h_x ≈ sol.A
    @test c.B ≈ sol.B
    @test c.Ω === sol.D  # should be nothing
    @test c.Q ≈ sol.Q
    @test c.η ≈ sol.η
    @test sol.retcode == :Success
end

@testset "Sparse RBC Module with include helper" begin
    m = @include_example_module(Examples.rbc_observables, 1, SparseFunctions())
    p = [0.5, 0.95, 0.2, 0.011]
    p_f = [0.02, 0.01, 0.012]
    sol = generate_perturbation(m, p; p_f)

    @test sol.y ≈ [5.936252888048733, 6.884057971014498]
    @test sol.x ≈ [47.39025414828824, 0.0]
end
