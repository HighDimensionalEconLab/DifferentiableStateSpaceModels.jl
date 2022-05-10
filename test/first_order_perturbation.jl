using DifferentiableStateSpaceModels, Symbolics, LinearAlgebra, Zygote, Test, BenchmarkTools
using DifferentiableStateSpaceModels.Examples
using DifferentiableStateSpaceModels: order_vector_by_symbols,
                                      fill_array_by_symbol_dispatch, all_fields_equal

# # # Use while testing internals
m = @include_example_module(Examples.rbc_observables)
# Basic Steady State
p_f = (ρ = 0.2, δ = 0.02, σ = 0.01, Ω_1 = 0.01)
p_d = (α = 0.5, β = 0.95)
c = SolverCache(m, Val(1), p_d)
sol = generate_perturbation(m, p_d, p_f, Val(1); cache = c) # manually passing in order
# generate_perturbation_derivatives!(m, p_d, p_f, c)
# ex = DifferentiableStateSpaceModels.exfiltrated    

@testset "Construction" begin
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
    sol = generate_perturbation(m, p_d, p_f, Val(1)) # manually passing in first derivative, but not really required
    sol = generate_perturbation(m, p_d, p_f) # Default is first-order
    @inferred generate_perturbation(m, p_d, p_f)
    @inferred generate_perturbation(m, p_d, p_f, Val(1))
    @test sol.y ≈ [5.936252888048733, 6.884057971014498]
    @test sol.x ≈ [47.39025414828825, 0.0]
    @test sol.retcode == :Success

    # Call all variables differentiated
    sol = generate_perturbation(m, merge(p_d, p_f), nothing)
    @test sol.y ≈ [5.936252888048733, 6.884057971014498]
    @test sol.x ≈ [47.39025414828825, 0.0]

    # With a prebuilt cache
    c = SolverCache(m, Val(1), p_d)
    sol = generate_perturbation(m, p_d, p_f; cache = c)
    @test sol.y ≈ [5.936252888048733, 6.884057971014498]
    @test sol.x ≈ [47.39025414828825, 0.0]
    @inferred generate_perturbation(m, p_d, p_f; cache = c)
end

@testset "Function Evaluation" begin
    m = @include_example_module(Examples.rbc_observables)
    p_f = (ρ = 0.2, δ = 0.02, σ = 0.01, Ω_1 = 0.01)
    p_d = (α = 0.5, β = 0.95)
    p_d_symbols = collect(Symbol.(keys(p_d)))
    c = SolverCache(m, Val(1), p_d)

    # Create parameter vector in the same ordering the internal algorithms would
    p = order_vector_by_symbols(merge(p_d, p_f), m.mod.m.p_symbols)

    y = zeros(m.n_y)
    x = zeros(m.n_x)

    m.mod.m.ȳ!(y, p)
    m.mod.m.x̄!(x, p)
    @test y ≈ [5.936252888048733, 6.884057971014498]
    @test x ≈ [47.39025414828825, 0.0]

    m.mod.m.H_yp!(c.H_yp, y, x, p)
    @test c.H_yp ≈ [0.028377570562199098 0.0; 0.0 0.0; 0.0 0.0; 0.0 0.0]

    m.mod.m.H_y!(c.H_y, y, x, p)
    @test c.H_y ≈ [-0.0283775705621991 0.0; 1.0 -1.0; 0.0 1.0; 0.0 0.0]

    m.mod.m.H_xp!(c.H_xp, y, x, p)
    @test c.H_xp ≈ [0.00012263591151906127 -0.011623494029190608
                    1.0 0.0
                    0.0 0.0
                    0.0 1.0]

    m.mod.m.H_x!(c.H_x, y, x, p)
    @test c.H_x ≈ [0.0 0.0
                   -0.98 0.0
                   -0.07263157894736837 -6.884057971014498
                   0.0 -0.2]

    m.mod.m.Ψ!(c.Ψ, y, x, p)
    @test c.Ψ[1] ≈
          [-0.009560768753410337 0.0 0.0 0.0 -2.0658808482697935e-5 0.0019580523687917364 0.0 0.0
           0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
           0.0 0.0 0.009560768753410338 0.0 0.0 0.0 0.0 0.0
           0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
           -2.0658808482697935e-5 0.0 0.0 0.0 -3.881681383327978e-6 0.00012263591151906127 0.0 0.0
           0.0019580523687917364 0.0 0.0 0.0 0.00012263591151906127 -0.011623494029190608 0.0 0.0
           0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
           0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]
    @test c.Ψ[2] ≈ zeros(8, 8)
    @test c.Ψ[3] ≈ [0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
                    0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
                    0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
                    0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
                    0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
                    0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
                    0.0 0.0 0.0 0.0 0.0 0.0 0.0007663134567721225 -0.07263157894736837
                    0.0 0.0 0.0 0.0 0.0 0.0 -0.07263157894736837 -6.884057971014498]
    @test c.Ψ[4] ≈ zeros(8, 8)

    m.mod.m.Γ!(c.Γ, p)
    @test c.Γ ≈ [0.01]

    m.mod.m.Ω!(c.Ω, p)
    @test c.Ω ≈ [0.01, 0.01]

    # The derivative ones dispatch by the derivative symbol
    fill_array_by_symbol_dispatch(m.mod.m.H_x_p!, c.H_x_p, p_d_symbols, y, x, p)
    @test c.H_x_p ≈ [[0.0 0.0
                      0.0 0.0
                      -0.4255060477077458 -26.561563542978472
                      0.0 0.0], [0.0 0.0
                                 0.0 0.0;
                                 0.0 0.0;
                                 0.0 0.0]]
    fill_array_by_symbol_dispatch(m.mod.m.H_yp_p!, c.H_yp_p, p_d_symbols, y, x, p)
    @test c.H_yp_p ≈ [[0.011471086498795562 0.0; 0.0 0.0; 0.0 0.0; 0.0 0.0],
                      [0.029871126907577997 0.0; 0.0 0.0; 0.0 0.0; 0.0 0.0]]

    fill_array_by_symbol_dispatch(m.mod.m.H_y_p!, c.H_y_p, p_d_symbols, y, x, p)
    @test c.H_y_p ≈ [[0.0 0.0; 0.0 0.0; 0.0 0.0; 0.0 0.0], [0.0 0.0
                                                            0.0 0.0;
                                                            0.0 0.0; 0.0 0.0]]

    fill_array_by_symbol_dispatch(m.mod.m.H_xp_p!, c.H_xp_p, p_d_symbols, y, x, p)
    @test c.H_xp_p ≈ [[0.000473180436623283 -0.06809527035753198
                       0.0 0.0
                       0.0 0.0
                       0.0 0.0], [0.00012909043317795924 -0.01223525687283222
                                  0.0 0.0
                                  0.0 0.0
                                  0.0 0.0]]

    fill_array_by_symbol_dispatch(m.mod.m.Γ_p!, c.Γ_p, p_d_symbols, p)

    @test c.Γ_p ≈ [[0.0], [0.0]]

    fill_array_by_symbol_dispatch(m.mod.m.H_p!, c.H_p, p_d_symbols, y, x, p)

    @test c.H_p ≈ [[-0.06809527035753199, 0.0, -26.561563542978472, 0.0],
                   [-0.1773225633743801, 0.0, 0.0, 0.0]]

    fill_array_by_symbol_dispatch(m.mod.m.Ω_p!, c.Ω_p, p_d_symbols, p)
    @test c.Ω_p ≈ [[0.0, 0.0], [0.0, 0.0]]
end

@testset "Evaluation into cache and check solution" begin
    m = @include_example_module(Examples.rbc_observables)
    p_f = (ρ = 0.2, δ = 0.02, σ = 0.01, Ω_1 = 0.01)
    p_d = (α = 0.5, β = 0.95)
    c = SolverCache(m, Val(1), p_d)
    sol = generate_perturbation(m, p_d, p_f; cache = c)
    # Create parameter vector in the same ordering the internal algorithms would

    @test c.y ≈ [5.936252888048733, 6.884057971014498]
    @test c.x ≈ [47.39025414828825, 0.0]
    @test c.H_yp ≈ [0.028377570562199098 0.0; 0.0 0.0; 0.0 0.0; 0.0 0.0]
    @test c.H_y ≈ [-0.0283775705621991 0.0; 1.0 -1.0; 0.0 1.0; 0.0 0.0]
    @test c.H_xp ≈ [0.00012263591151906127 -0.011623494029190608
                    1.0 0.0
                    0.0 0.0
                    0.0 1.0]
    @test c.H_x ≈ [0.0 0.0
                   -0.98 0.0
                   -0.07263157894736837 -6.884057971014498
                   0.0 -0.2]
    @test c.Ψ[1] ≈
          [-0.009560768753410337 0.0 0.0 0.0 -2.0658808482697935e-5 0.0019580523687917364 0.0 0.0
           0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
           0.0 0.0 0.009560768753410338 0.0 0.0 0.0 0.0 0.0
           0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
           -2.0658808482697935e-5 0.0 0.0 0.0 -3.881681383327978e-6 0.00012263591151906127 0.0 0.0
           0.0019580523687917364 0.0 0.0 0.0 0.00012263591151906127 -0.011623494029190608 0.0 0.0
           0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
           0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]
    @test c.Ψ[2] ≈ zeros(8, 8)
    @test c.Ψ[3] ≈ [0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
                    0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
                    0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
                    0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
                    0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
                    0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
                    0.0 0.0 0.0 0.0 0.0 0.0 0.0007663134567721225 -0.07263157894736837
                    0.0 0.0 0.0 0.0 0.0 0.0 -0.07263157894736837 -6.884057971014498]
    @test c.Ψ[4] ≈ zeros(8, 8)
    @test c.Γ ≈ [0.01]
    @test c.Ω ≈ [0.01, 0.01]

    # solution tests
    @test c.y ≈ [5.936252888048733, 6.884057971014498]
    @test c.x ≈ [47.39025414828824, 0.0]
    @test c.g_x ≈ [0.0957964300241661 0.6746869652586178
                   0.07263157894736878 6.884057971014507]
    @test c.h_x ≈ [0.9568351489232028 6.209371005755889; -1.5076865909646354e-18 0.2]
    @test c.Σ ≈ [1e-4]
    @test c.Ω ≈ [0.01, 0.01]
    @test c.η == reshape([0; -1], 2, m.n_ϵ)
    @test c.B ≈ [0.0; -0.01]
    @test c.Q ≈ [1.0 0.0 0.0 0.0; 0.0 0.0 1.0 0.0]
    @test c.C_1 ≈ [0.0957964300241661 0.6746869652585828; 1.0 0.0]
    @test Array(c.V) ≈ [0.07005411173180227 0.00015997603451513485
                        0.00015997603451513485 0.00010416666666666667]

    # Check the solution type matches these all
    fields_to_compare = (:y, :x, :g_x, :B, :Q, :η, :Γ)
    @test all_fields_equal(c, sol, fields_to_compare)
    @test c.h_x ≈ sol.A
    @test c.C_1 ≈ sol.C
    @test c.V ≈ sol.x_ergodic.Σ # Covariance matrix in MvNormal
    @test c.Ω ≈ sqrt.(diag(sol.D.Σ))
end

@testset "Evaluate Derivatives into cache" begin
    m = @include_example_module(Examples.rbc_observables)
    p_f = (ρ = 0.2, δ = 0.02, σ = 0.01, Ω_1 = 0.01)
    p_d = (α = 0.5, β = 0.95)
    p_d_symbols = collect(Symbol.(keys(p_d)))
    c = SolverCache(m, Val(1), p_d)
    sol = generate_perturbation(m, p_d, p_f; cache = c)
    @test :Success == generate_perturbation_derivatives!(m, p_d, p_f, c) # Solves and fills the cache
    @inferred generate_perturbation_derivatives!(m, p_d, p_f, c)

    @test c.H_x_p ≈ [[0.0 0.0
                      0.0 0.0
                      -0.4255060477077458 -26.561563542978472
                      0.0 0.0], [0.0 0.0
                                 0.0 0.0;
                                 0.0 0.0; 0.0 0.0]]

    @test c.H_yp_p ≈ [[0.011471086498795562 0.0; 0.0 0.0; 0.0 0.0; 0.0 0.0],
                      [0.029871126907577997 0.0; 0.0 0.0; 0.0 0.0; 0.0 0.0]]
    @test c.H_y_p ≈ [[0.0 0.0; 0.0 0.0; 0.0 0.0; 0.0 0.0], [0.0 0.0
                                                            0.0 0.0;
                                                            0.0 0.0; 0.0 0.0]]
    @test c.H_xp_p ≈ [[0.000473180436623283 -0.06809527035753198
                       0.0 0.0
                       0.0 0.0
                       0.0 0.0], [0.00012909043317795924 -0.01223525687283222
                                  0.0 0.0
                                  0.0 0.0
                                  0.0 0.0]]
    @test c.Γ_p ≈ [[0.0], [0.0]]
    @test c.H_p ≈ [[-0.06809527035753199, 0.0, -26.561563542978472, 0.0],
                   [-0.1773225633743801, 0.0, 0.0, 0.0]]
    @test c.Ω_p ≈ [[0.0, 0.0], [0.0, 0.0]]

    # solution
    @test hcat(c.y_p...) ≈ [55.78596896689701 76.10141579073955
                            66.89124302798608 105.01995379122064]
    @test hcat(c.x_p...) ≈ [555.2637030544529 1445.9269000240533; 0.0 0.0]
    @test c.g_x_p ≈ [[-0.12465264193058262 5.596211904442805
                      -1.2823781479976832e-15 66.89124302798608],
                     [-1.6946742377792863 -0.8343618226192915
                      -1.1080332409972313 105.01995379122064]]
    @test c.h_x_p ≈ [[0.12465264193058134 61.29503112354326; 0.0 0.0],
                     [0.586640996782055 105.85431561383992; 0.0 0.0]]
    @test c.Σ_p ≈ [[0.0], [0.0]]
    @test hcat(c.Ω_p...) ≈ [0.0 0.0; 0.0 0.0]
    @test c.B_p ≈ [[0.0; 0.0], [0.0
                                0.0]]
    @test c.C_1_p ≈ [[-0.12465264193057919 5.596211904442171; 0.0 0.0],
                     [-1.6946742377792825 -0.8343618226202246; 0.0 0.0]]
    @test c.V_p ≈ [[1.584528257999749 0.0015841155991973127; 0.0015841155991973127 0.0],
                   [3.336643488330957 0.002750404724942799; 0.002750404724942799 0.0]]
end

@testset "Steady State with Initial Conditions" begin
    m = @include_example_module(Examples.rbc_solve_steady_state)
    p_f = (ρ = 0.2, δ = 0.02, σ = 0.01, Ω_1 = 0.01)
    p_d = (α = 0.5, β = 0.95)

    sol = generate_perturbation(m, p_d, p_f)
    @inferred generate_perturbation(m, p_d, p_f)
    @test sol.y ≈ [5.936252888048733, 6.884057971014498]
    @test sol.x ≈ [47.39025414828825, 0.0]
    @test sol.retcode == :Success
end

@testset "Steady State with Initial Conditions" begin
    m = @include_example_module(Examples.rbc_solve_steady_state)
    p_f = (ρ = 0.2, δ = 0.02, σ = 0.01, Ω_1 = 0.01)
    p_d = (α = 0.5, β = 0.95)

    # Starting at steady state in optimizer, so should be 0 or 1 iteration
    settings = PerturbationSolverSettings(; print_level = 0)
    sol = generate_perturbation(m, p_d, p_f; settings)
    @inferred generate_perturbation(m, p_d, p_f)
    @test sol.y ≈ [5.936252888048733, 6.884057971014498]
    @test sol.x ≈ [47.39025414828825, 0.0]
    @test sol.retcode == :Success

    # At an initial condition outside of the steady state
    m = @include_example_module(Examples.rbc_solve_steady_state_different_iv)
    p_f = (ρ = 0.2, δ = 0.02, σ = 0.01, Ω_1 = 0.01)
    p_d = (α = 0.5, β = 0.95)

    sol = generate_perturbation(m, p_d, p_f; settings)
    @test sol.y ≈ [5.936252888048733, 6.884057971014498]
    @test sol.x ≈ [47.39025414828825, 0.0]
    @test sol.retcode == :Success

    # Change the order of the parameters/etc.
    p_f = (σ = 0.01, δ = 0.02, Ω_1 = 0.01)
    p_d = (α = 0.5, ρ = 0.2, β = 0.95)
    sol = generate_perturbation(m, p_d, p_f; settings)  # no p_f and rearranged
    @test sol.y ≈ [5.936252888048733, 6.884057971014498]
    @test sol.x ≈ [47.39025414828825, 0.0]
    @test sol.retcode == :Success

    m = @include_example_module(Examples.rbc_solve_steady_state_different_iv)
    p_f = (ρ = 0.2, δ = 0.02, σ = 0.01, Ω_1 = 0.01)
    p_d = (α = 0.5, β = 0.95)

    settings = PerturbationSolverSettings(; print_level = 0, nlsolve_iterations = 2)  # simulate failure by insufficient iterationrs
    sol = generate_perturbation(m, p_d, p_f; settings)
    @test sol.retcode == :SteadyStateFailure
end

@testset "Construction and solution with no Omega" begin
    m = @include_example_module(Examples.rbc)
    p_f = (ρ = 0.2, δ = 0.02, σ = 0.01)
    p_d = (α = 0.5, β = 0.95)
    sol = generate_perturbation(m, p_d, p_f)
    @inferred generate_perturbation(m, p_d, p_f)
    @test sol.y ≈ [5.936252888048733, 6.884057971014498]
    @test sol.x ≈ [47.39025414828825, 0.0]
    @test sol.retcode == :Success

    c = SolverCache(m, Val(1), p_d)
    sol = generate_perturbation(m, p_d, p_f; cache = c)
    @test :Success == generate_perturbation_derivatives!(m, p_d, p_f, c)  # Solves and fills the cache
    @inferred generate_perturbation(m, p_d, p_f; cache = c)
    @inferred generate_perturbation_derivatives!(m, p_d, p_f, c)

    @inferred generate_perturbation(m, p_d, p_f; cache = c)
    @test sol.retcode == :Success

    @test sol.y ≈ [5.936252888048733, 6.884057971014498]
    @test sol.x ≈ [47.39025414828824, 0.0]
    @test hcat(c.y_p...) ≈ [55.78596896689701 76.10141579073955
                            66.89124302798608 105.01995379122064]
    @test hcat(c.x_p...) ≈ [555.2637030544529 1445.9269000240533; 0.0 0.0]
    @test c.g_x ≈ [0.0957964300241661 0.6746869652586178
                   0.07263157894736878 6.884057971014507]
    @test c.h_x ≈ [0.9568351489232028 6.209371005755889; -1.5076865909646354e-18 0.2]
    @test c.g_x_p ≈ [[-0.12465264193058262 5.596211904442805
                      -1.2823781479976832e-15 66.89124302798608],
                     [-1.6946742377792863 -0.8343618226192915
                      -1.1080332409972313 105.01995379122064]]
    @test c.h_x_p ≈ [[0.12465264193058134 61.29503112354326; 0.0 0.0],
                     [0.586640996782055 105.85431561383992; 0.0 0.0]]
    @test c.Σ ≈ [1e-4]
    @test c.Σ_p ≈ [[0.0], [0.0]]
end

@testset "Evaluate rbc_observables_separate_variance derivatives into cache" begin
    m = @include_example_module(Examples.rbc_observables_separate_variance)

    p_f = (δ = 0.02, σ = 0.01, Ω_2 = 0.012)
    p_d = (α = 0.5, β = 0.95, ρ = 0.2, Ω_1 = 0.011)
    c = SolverCache(m, Val(1), p_d)
    sol = generate_perturbation(m, p_d, p_f; cache = c)
    @test :Success == generate_perturbation_derivatives!(m, p_d, p_f, c)

    @test c.y ≈ [5.936252888048733, 6.884057971014498]
    @test c.x ≈ [47.39025414828824, 0.0]
    @test hcat(c.y_p...) ≈ [55.78596896689701 76.10141579073955 0.0 0.0
                            66.89124302798608 105.01995379122064 0.0 0.0]
    @test hcat(c.x_p...) ≈ [555.2637030544529 1445.9269000240533 0.0 0.0; 0.0 0.0 0.0 0.0]
    @test c.g_x ≈ [0.0957964300241661 0.6746869652586178
                   0.07263157894736878 6.884057971014507]
    @test c.h_x ≈ [0.9568351489232028 6.209371005755889; -1.5076865909646354e-18 0.2]
    @test c.g_x_p ≈ [[-0.12465264193058262 5.596211904442805
                      -1.2823781479976832e-15 66.89124302798608],
                     [-1.6946742377792863 -0.8343618226192915
                      -1.1080332409972313 105.01995379122064],
                     [1.3921362894665956e-18 0.29450084693461975; 0.0 0.0], zeros(2, 2)]
    @test c.h_x_p ≈ [[0.12465264193058134 61.29503112354326; 0.0 0.0],
                     [0.586640996782055 105.85431561383992; 0.0 0.0],
                     [-1.3921362894665956e-18 -0.29450084693461975; 4.7271045362244914e-18 1.0],
                     zeros(2, 2)]
    @test c.Σ ≈ [1e-4]
    @test c.Σ_p ≈ [[0.0], [0.0], [0.0], [0.0]]

    @test c.Ω ≈ [0.011, 0.012]
    @test hcat(c.Ω_p...) ≈ [0.0 0.0 0.0 1.0; 0.0 0.0 0.0 0.0]

    @test c.η == reshape([0; -1], 2, m.n_ϵ)
    @test c.Q ≈ [1.0 0.0 0.0 0.0; 0.0 0.0 1.0 0.0]

    # Should be identical to c
    @test sol.n_y == 2
    @test sol.n_x == 2
    @test sol.n_ϵ == 1
    @test sol.n_z == 2
    @test c.y ≈ sol.y
    @test c.x ≈ sol.x
    @test c.g_x ≈ sol.g_x
    @test c.h_x ≈ sol.A
    @test c.B ≈ sol.B
    @test c.Ω ≈ sqrt.(diag(sol.D.Σ))
    @test c.Q ≈ sol.Q
    @test c.η ≈ sol.η
    @test sol.retcode == :Success
end

@testset "Construction and solution no p_f" begin
    m = @include_example_module(Examples.rbc)
    p_f = nothing
    p_d = (α = 0.5, β = 0.95, ρ = 0.2, δ = 0.02, σ = 0.01)
    c = SolverCache(m, Val(1), p_d)
    sol = generate_perturbation(m, p_d, p_f; cache = c)
    generate_perturbation_derivatives!(m, p_d, p_f, c)
    @inferred generate_perturbation(m, p_d, p_f; cache = c)
    @inferred generate_perturbation_derivatives!(m, p_d, p_f, c)

    @test c.y ≈ [5.936252888048733, 6.884057971014498]
    @test c.x ≈ [47.39025414828825, 0.0]
    @test c.H_yp ≈ [0.028377570562199098 0.0; 0.0 0.0; 0.0 0.0; 0.0 0.0]
    @test c.H_y ≈ [-0.0283775705621991 0.0; 1.0 -1.0; 0.0 1.0; 0.0 0.0]
    @test c.H_xp ≈ [0.00012263591151906127 -0.011623494029190608
                    1.0 0.0
                    0.0 0.0
                    0.0 1.0]
    @test c.H_x ≈ [0.0 0.0
                   -0.98 0.0
                   -0.07263157894736837 -6.884057971014498
                   0.0 -0.2]
    @test c.Ψ[1] ≈
          [-0.009560768753410337 0.0 0.0 0.0 -2.0658808482697935e-5 0.0019580523687917364 0.0 0.0
           0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
           0.0 0.0 0.009560768753410338 0.0 0.0 0.0 0.0 0.0
           0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
           -2.0658808482697935e-5 0.0 0.0 0.0 -3.881681383327978e-6 0.00012263591151906127 0.0 0.0
           0.0019580523687917364 0.0 0.0 0.0 0.00012263591151906127 -0.011623494029190608 0.0 0.0
           0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
           0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]
    @test c.Ψ[2] ≈ zeros(8, 8)
    @test c.Ψ[3] ≈ [0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
                    0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
                    0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
                    0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
                    0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
                    0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
                    0.0 0.0 0.0 0.0 0.0 0.0 0.0007663134567721225 -0.07263157894736837
                    0.0 0.0 0.0 0.0 0.0 0.0 -0.07263157894736837 -6.884057971014498]
    @test c.Ψ[4] ≈ zeros(8, 8)
    @test c.Γ ≈ [0.01]
    @test c.Γ_p ≈ [[0.0], [0.0], [0.0], [0.0], [1.0]]
    @test hcat(c.H_p...) ≈
          [-0.06809527035753199 -0.1773225633743801 0.0 0.16003361344537806 0.0
           0.0 0.0 0.0 47.39025414828824 0.0
           -26.561563542978472 0.0 0.0 0.0 0.0
           0.0 0.0 -0.0 0.0 0.0]

    # solution tests
    @test c.y ≈ [5.936252888048733, 6.884057971014498]
    @test c.x ≈ [47.39025414828824, 0.0]
    @test hcat(c.y_p...) ≈ [55.78596896689701 76.10141579073955 0.0 -116.07178189943077 0.0
                            66.89124302798608 105.01995379122064 0.0 -94.78050829657676 0.0]
    @test hcat(c.x_p...) ≈ [555.2637030544529 1445.9269000240533 0.0 -1304.9490272717098 0.0
                            0.0 0.0 0.0 0.0 0.0]
    @test c.g_x ≈ [0.0957964300241661 0.6746869652586178
                   0.07263157894736878 6.884057971014507]
    @test c.h_x ≈ [0.9568351489232028 6.209371005755889; -1.5076865909646354e-18 0.2]
    @test c.g_x_p ≈ [[-0.12465264193058262 5.596211904442805
                      -1.2823781479976832e-15 66.89124302798608],
                     [-1.6946742377792863 -0.8343618226192915
                      -1.1080332409972313 105.01995379122064],
                     [4.640454298222595e-19 0.2945008469346586; 0.0 0.0],
                     [0.6277268890968761 -5.0369653468355455
                      1.0000000000000024 -94.78050829657676], [0.0 0.0
                                                               0.0 0.0]]
    @test c.h_x_p ≈ [[0.12465264193058134 61.29503112354326; 0.0 0.0],
                     [0.586640996782055 105.85431561383992; 0.0 0.0],
                     [-4.640454298222595e-19 -0.2945008469346586
                      1.5757015120748295e-18 1.0], [-0.6277268890968737 -89.74354294974118
                                                    0.0 0.0], [0.0 0.0
                                                               0.0 0.0]]
    @test c.Σ ≈ [1e-4]
    @test c.Σ_p ≈ [[0.0], [0.0], [0.0], [0.0], [0.02]]
end

@testset "Construction and solution, multiple shocks" begin
    m = @include_example_module(Examples.rbc_multiple_shocks)
    p_f = nothing
    p_d = (α = 0.5, β = 0.95, ρ = 0.2, δ = 0.02, σ = 0.01)
    c = SolverCache(m, Val(1), p_d)
    sol = generate_perturbation(m, p_d, p_f; cache = c)
    generate_perturbation_derivatives!(m, p_d, p_f, c)
    @inferred generate_perturbation(m, p_d, p_f; cache = c)
    @inferred generate_perturbation_derivatives!(m, p_d, p_f, c)

    @test sol.y ≈ [5.936252888048733, 6.884057971014498]
    @test sol.x ≈ [47.39025414828825, 0.0]
    @test c.H_yp ≈ [0.028377570562199098 0.0; 0.0 0.0; 0.0 0.0; 0.0 0.0]
    @test c.H_y ≈ [-0.0283775705621991 0.0; 1.0 -1.0; 0.0 1.0; 0.0 0.0]
    @test c.H_xp ≈ [0.00012263591151906127 -0.011623494029190608
                    1.0 0.0
                    0.0 0.0
                    0.0 1.0]
    @test c.H_x ≈ [0.0 0.0
                   -0.98 0.0
                   -0.07263157894736837 -6.884057971014498
                   0.0 -0.2]
    @test c.Ψ[1] ≈
          [-0.009560768753410337 0.0 0.0 0.0 -2.0658808482697935e-5 0.0019580523687917364 0.0 0.0
           0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
           0.0 0.0 0.009560768753410338 0.0 0.0 0.0 0.0 0.0
           0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
           -2.0658808482697935e-5 0.0 0.0 0.0 -3.881681383327978e-6 0.00012263591151906127 0.0 0.0
           0.0019580523687917364 0.0 0.0 0.0 0.00012263591151906127 -0.011623494029190608 0.0 0.0
           0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
           0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]
    @test c.Ψ[2] ≈ zeros(8, 8)
    @test c.Ψ[3] ≈ [0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
                    0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
                    0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
                    0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
                    0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
                    0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
                    0.0 0.0 0.0 0.0 0.0 0.0 0.0007663134567721225 -0.07263157894736837
                    0.0 0.0 0.0 0.0 0.0 0.0 -0.07263157894736837 -6.884057971014498]
    @test c.Ψ[4] ≈ zeros(8, 8)
    @test c.Γ ≈ [0.01 0.0; 0.01 0.01]
    @test c.Γ_p ≈ [[0.0 0.0; 0.0 0.0], [0.0 0.0
                                        0.0 0.0], [0.0 0.0
                                                   0.0 0.0], [0.0 0.0
                                                              0.0 0.0], [1.0 0.0
                                                                         1.0 1.0]]

    # solution tests
    @test c.y ≈ [5.936252888048733, 6.884057971014498]
    @test c.x ≈ [47.39025414828824, 0.0]
    @test hcat(c.y_p...) ≈ [55.78596896689701 76.10141579073955 0.0 -116.07178189943077 0.0
                            66.89124302798608 105.01995379122064 0.0 -94.78050829657676 0.0]
    @test hcat(c.x_p...) ≈ [555.2637030544529 1445.9269000240533 0.0 -1304.9490272717098 0.0
                            0.0 0.0 0.0 0.0 0.0]
    @test c.g_x ≈ [0.0957964300241661 0.6746869652586178
                   0.07263157894736878 6.884057971014507]
    @test c.h_x ≈ [0.9568351489232028 6.209371005755889; -1.5076865909646354e-18 0.2]
    @test c.g_x_p ≈ [[-0.12465264193058262 5.596211904442805
                      -1.2823781479976832e-15 66.89124302798608],
                     [-1.6946742377792863 -0.8343618226192915
                      -1.1080332409972313 105.01995379122064],
                     [4.640454298222595e-19 0.2945008469346586; 0.0 0.0],
                     [0.6277268890968761 -5.0369653468355455
                      1.0000000000000024 -94.78050829657676], [0.0 0.0
                                                               0.0 0.0]]
    @test c.h_x_p ≈ [[0.12465264193058134 61.29503112354326; 0.0 0.0],
                     [0.586640996782055 105.85431561383992; 0.0 0.0],
                     [-4.640454298222595e-19 -0.2945008469346586
                      1.5757015120748295e-18 1.0], [-0.6277268890968737 -89.74354294974118
                                                    0.0 0.0], [0.0 0.0
                                                               0.0 0.0]]
    @test c.Σ ≈ [0.0001 0.0001; 0.0001 0.0002]
    @test c.Σ_p[1] ≈ zeros(2, 2)
    @test c.Σ_p[2] ≈ zeros(2, 2)
    @test c.Σ_p[3] ≈ zeros(2, 2)
    @test c.Σ_p[4] ≈ zeros(2, 2)
    @test c.Σ_p[5] ≈ [0.02 0.02; 0.02 0.04]
end

@testset "Schur Decomposition Failure" begin
    m = @include_example_module(Examples.rbc)
    p_f = (ρ = 0.2, δ = 0.02, σ = 0.01)
    p_d = (α = 100.0, β = 0.0) # garbage parameters
    settings = PerturbationSolverSettings(; print_level = 0)
    sol = generate_perturbation(m, p_d, p_f; settings)
    @test sol.retcode == :LAPACK_Error
end

# Version with rethrowing exceptions
@testset "Schur Decomposition Failure Exceptions" begin
    m = @include_example_module(Examples.rbc)
    p_f = (ρ = 0.2, δ = 0.02, σ = 0.01)
    p_d = (α = 100.0, β = 0.0) # garbage parameters
    settings = PerturbationSolverSettings(; print_level = 0, rethrow_exceptions = true)
    @test_throws LAPACKException generate_perturbation(m, p_d, p_f; settings)
end

@testset "BK Condition Failure" begin
    m = @include_example_module(Examples.rbc)
    p_f = nothing
    p_d = (α = 0.5, β = 0.95, ρ = 1.01, δ = 0.02, σ = 0.01) # rho > 1
    settings = PerturbationSolverSettings(; print_level = 0)
    sol = generate_perturbation(m, p_d, p_f; settings)
    @test sol.retcode == :Blanchard_Kahn_Failure
end

@testset "BK Condition Failure Exceptions" begin
    m = @include_example_module(Examples.rbc)
    p_f = nothing
    p_d = (α = 0.5, β = 0.95, ρ = 1.01, δ = 0.02, σ = 0.01) # rho > 1
    settings = PerturbationSolverSettings(; print_level = 0, rethrow_exceptions = true)
    @test_throws ErrorException generate_perturbation(m, p_d, p_f; settings)
end

@testset "Function evaluation Failure" begin
    m = @include_example_module(Examples.rbc)
    p_f = nothing
    p_d = (α = -0.5, β = 0.95, ρ = 0.2, δ = 0.02, σ = 0.01)
    settings = PerturbationSolverSettings(; print_level = 0)
    sol = generate_perturbation(m, p_d, p_f; settings)
    @test sol.retcode == :Evaluation_Error
end

@testset "Function evaluation Failure Exceptions" begin
    m = @include_example_module(Examples.rbc)
    p_f = nothing
    p_d = (α = -0.5, β = 0.95, ρ = 0.2, δ = 0.02, σ = 0.01)
    settings = PerturbationSolverSettings(; print_level = 0, rethrow_exceptions = true)
    @test_throws DomainError generate_perturbation(m, p_d, p_f; settings)
end

@testset "Ergodic distribution failure" begin
    m = @include_example_module(Examples.rbc)
    p_f = nothing
    p_d = (α = 0.5, β = 0.95, ρ = 0.2, δ = 0.02, σ = 10000)
    settings = PerturbationSolverSettings(; print_level = 0)
    sol = generate_perturbation(m, p_d, p_f; settings)
    @test sol.retcode != :Success
end

@testset "Callbacks" begin
    calculate_steady_state_callback_triggered = false
    function calculate_steady_state_callback(ret, m, cache, settings, p)
        calculate_steady_state_callback_triggered = true

        return nothing
    end
    evaluate_functions_callback_triggered = false
    function evaluate_functions_callback(ret, m, cache, settings, p)
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

    settings = PerturbationSolverSettings(; print_level = 0,
                                          calculate_steady_state_callback,
                                          evaluate_functions_callback,
                                          solve_first_order_callback,
                                          solve_first_order_p_callback)

    m = @include_example_module(Examples.rbc_observables)
    p_f = (ρ = 0.2, δ = 0.02, σ = 0.01, Ω_1 = 0.01)
    p_d = (α = 0.5, β = 0.95)

    c = SolverCache(m, Val(1), p_d)
    sol = generate_perturbation(m, p_d, p_f; settings, cache = c)
    generate_perturbation_derivatives!(m, p_d, p_f, c, ; settings)

    @test calculate_steady_state_callback_triggered == true
    @test evaluate_functions_callback_triggered == true
    @test solve_first_order_callback_triggered == true
    @test solve_first_order_p_callback_triggered == true
end

@testset "First Order Pullback inference" begin
    m = @include_example_module(Examples.rbc_observables)
    p_f = (ρ = 0.2, δ = 0.02, σ = 0.01, Ω_1 = 0.01)
    p_d = (α = 0.5, β = 0.95)

    c = SolverCache(m, Val(1), p_d)
    sol = generate_perturbation(m, p_d, p_f; cache = c)

    _, pb = Zygote.pullback(generate_perturbation, m, p_d, p_f, Val(1))
    @inferred pb(sol)
end