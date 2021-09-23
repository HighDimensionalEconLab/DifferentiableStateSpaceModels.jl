using DifferentiableStateSpaceModels, Symbolics, Test
using DifferentiableStateSpaceModels.Examples
using DifferentiableStateSpaceModels: arrange_vector_from_symbols, fill_array_by_symbol_dispatch

@testset "Second Order Construction" begin
    m = @include_example_module(Examples.rbc_observables)

    # bookkeeping tests
    @test m.n_y == 2
    @test m.n_x == 2
    @test m.n_p == 6
    @test m.n_ϵ == 1
    @test m.n_z == 2
    @test m.η == reshape([0; -1], 2, m.n_ϵ)
    @test m.Q == [1.0 0.0 0.0 0.0; 0.0 0.0 1.0 0.0]

    # function tests (steady state)
    # Basic Steady State
    p_f = (ρ = 0.2, δ = 0.02, σ = 0.01, Ω_1 = 0.01)
    p_d = (α = 0.5, β = 0.95)
    sol = second_order_perturbation(m, p_d, p_f)
    @inferred second_order_perturbation(m, p_d, p_f)
    @test sol.y ≈ [5.936252888048733, 6.884057971014498]
    @test sol.x ≈ [47.39025414828825, 0.0]
    @test sol.retcode == :Success

    # Call all variables differentiated
    sol = second_order_perturbation(m, merge(p_d, p_f))
    @test sol.y ≈ [5.936252888048733, 6.884057971014498]
    @test sol.x ≈ [47.39025414828825, 0.0]

    # With a prebuilt cache
    @test_throws AssertionError second_order_perturbation(m, p_d, p_f;cache =  SolverCache(m, Val(2), 1)) # Wrong cache size
    c = SolverCache(m, Val(2), length(p_d))
    sol = second_order_perturbation(m, p_d, p_f;cache = c)
    @test sol.y ≈ [5.936252888048733, 6.884057971014498]
    @test sol.x ≈ [47.39025414828825, 0.0]
    @inferred second_order_perturbation(m, p_d, p_f;cache = c)
end

@testset "Second Order Function Evaluation" begin
    m = @include_example_module(Examples.rbc_observables)
    p_f = (ρ = 0.2, δ = 0.02, σ = 0.01, Ω_1 = 0.01)
    p_d = (α = 0.5, β = 0.95)
    p_d_symbols =collect(Symbol.(keys(p_d)))  #The order of derivatives in p_d
    c = SolverCache(m, Val(2), length(p_d))

    # Create parameter vector in the same ordering the internal algorithms would
    p = arrange_vector_from_symbols(merge(p_d, p_f), m.mod.p_symbols)

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

    # The derivative ones dispatch by the derivative symbol
    fill_array_by_symbol_dispatch(m.mod.H_x_p!, c.H_x_p, p_d_symbols, y, x, p)
    @test c.H_x_p ≈ [
        [
            0.0 0.0
            0.0 0.0
            -0.4255060477077458 -26.561563542978472
            0.0 0.0
        ],
        [0.0 0.0; 0.0 0.0; 0.0 0.0; 0.0 0.0],
    ]
    fill_array_by_symbol_dispatch(m.mod.H_yp_p!, c.H_yp_p, p_d_symbols, y, x, p)
    @test c.H_yp_p ≈ [
        [0.011471086498795562 0.0; 0.0 0.0; 0.0 0.0; 0.0 0.0],
        [0.029871126907577997 0.0; 0.0 0.0; 0.0 0.0; 0.0 0.0],
    ]

    fill_array_by_symbol_dispatch(m.mod.H_y_p!, c.H_y_p, p_d_symbols, y, x, p)
    @test c.H_y_p ≈
        [[0.0 0.0; 0.0 0.0; 0.0 0.0; 0.0 0.0], [0.0 0.0; 0.0 0.0; 0.0 0.0; 0.0 0.0]]


    fill_array_by_symbol_dispatch(m.mod.H_xp_p!, c.H_xp_p, p_d_symbols, y, x, p)
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

    fill_array_by_symbol_dispatch(m.mod.Γ_p!, c.Γ_p, p_d_symbols, p)

    @test c.Γ_p ≈ [[0.0], [0.0]]

    fill_array_by_symbol_dispatch(m.mod.H_p!, c.H_p, p_d_symbols, y, x, p)

    @test c.H_p ≈ [
        [-0.06809527035753199, 0.0, -26.561563542978472, 0.0],
        [ -0.1773225633743801, 0.0, 0.0, 0.0]]


    # Second order checks!
    m.mod.Ψ_yp!(c.Ψ_yp, y, x, p)
    @test c.Ψ_yp[1] ≈ [
        [
            0.0048317190660755365 0.0 0.0 0.0 6.960218465183541e-6 -0.0006596930439853133 0.0 0.0
            0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
            0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
            0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
            6.960218465183541e-6 0.0 0.0 0.0 6.53894208439611e-7 -2.0658808482697952e-5 0.0 0.0
            -0.0006596930439853133 0.0 0.0 0.0 -2.0658808482697952e-5 0.0019580523687917372 0.0 0.0
            0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
            0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        ],
        zeros(8, 8),
        zeros(8, 8),
        zeros(8, 8),
    ]
    @test c.Ψ_yp[2] ≈ [zeros(8, 8), zeros(8, 8), zeros(8, 8), zeros(8, 8)]

    m.mod.Ψ_y!(c.Ψ_y, y, x, p)
    @test c.Ψ_y[1] ≈ [
        [
            0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
            0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
            0.0 0.0 -0.0048317190660755365 0.0 0.0 0.0 0.0 0.0
            0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
            0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
            0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
            0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
            0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        ],
        zeros(8, 8),
        zeros(8, 8),
        zeros(8, 8),
    ]

    @test c.Ψ_y[2] ≈ [zeros(8, 8), zeros(8, 8), zeros(8, 8), zeros(8, 8)]

    m.mod.Ψ_xp!(c.Ψ_xp, y, x, p)
    @test c.Ψ_xp[1] ≈ [
        [
            6.9602184651835404e-6 0.0 0.0 0.0 6.53894208439611e-7 -2.0658808482697952e-5 0.0 0.0
            0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
            0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
            0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
            6.53894208439611e-7 0.0 0.0 0.0 2.0477213369556226e-7 -3.881681383327981e-6 0.0 0.0
            -2.0658808482697952e-5 0.0 0.0 0.0 -3.881681383327981e-6 0.00012263591151906138 0.0 0.0
            0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
            0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        ],
        zeros(8, 8),
        zeros(8, 8),
        zeros(8, 8),
    ]
    @test c.Ψ_xp[2] ≈ [
        [
            -0.0006596930439853132 0.0 0.0 0.0 -2.0658808482697952e-5 0.0019580523687917372 0.0 0.0
            0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
            0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
            0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
            -2.0658808482697952e-5 0.0 0.0 0.0 -3.881681383327981e-6 0.00012263591151906138 0.0 0.0
            0.0019580523687917372 0.0 0.0 0.0 0.00012263591151906138 -0.011623494029190612 0.0 0.0
            0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
            0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        ],
        zeros(8, 8),
        zeros(8, 8),
        zeros(8, 8),
    ]
    m.mod.Ψ_x!(c.Ψ_x, y, x, p)
    @test c.Ψ_x[1] ≈ [
        zeros(8, 8),
        zeros(8, 8),
        [
            0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
            0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
            0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
            0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
            0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
            0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
            0.0 0.0 0.0 0.0 0.0 0.0 -2.4255412970806024e-5 0.000766313456772123
            0.0 0.0 0.0 0.0 0.0 0.0 0.000766313456772123 -0.07263157894736838
        ],
        zeros(8, 8),
    ]
    @test c.Ψ_x[2] ≈ [
        zeros(8, 8),
        zeros(8, 8),
        [
            0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
            0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
            0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
            0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
            0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
            0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
            0.0 0.0 0.0 0.0 0.0 0.0 0.000766313456772123 -0.07263157894736838
            0.0 0.0 0.0 0.0 0.0 0.0 -0.07263157894736838 -6.884057971014497
        ],
        zeros(8, 8),
    ]
    
    fill_array_by_symbol_dispatch(m.mod.Ψ_p!, c.Ψ_p, p_d_symbols, y, x, p)    
    @test c.Ψ_p[1] ≈ [
        [
            -0.0038647566790457792 0.0 0.0 0.0 -7.971028956261657e-5 0.011471086498795566 0.0 0.0
            0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
            0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
            0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
            -7.971028956261657e-5 0.0 0.0 0.0 -1.238935629209052e-5 0.00047318043662328334 0.0 0.0
            0.011471086498795566 0.0 0.0 0.0 0.00047318043662328334 -0.068095270357532 0.0 0.0
            0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
            0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        ],
        zeros(8, 8),
        [
            0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
            0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
            0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
            0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
            0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
            0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
            0.0 0.0 0.0 0.0 0.0 0.0 0.002956756561550659 -0.42550604770774586
            0.0 0.0 0.0 0.0 0.0 0.0 -0.42550604770774586 -26.56156354297847
        ],
        zeros(8, 8),
    ]

    @test c.Ψ_p[2] ≈ [
        [
            -0.010063967108852991 0.0 0.0 0.0 -2.1746114192313637e-5 0.0020611077566228815 0.0 0.0
            0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
            0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
            0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
            -2.1746114192313637e-5 0.0 0.0 0.0 -4.085980403503138e-6 0.00012909043317795935 0.0 0.0
            0.0020611077566228815 0.0 0.0 0.0 0.00012909043317795935 -0.012235256872832223 0.0 0.0
            0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
            0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        ],
        zeros(8, 8),
        zeros(8, 8),
        zeros(8, 8),
    ]
end
