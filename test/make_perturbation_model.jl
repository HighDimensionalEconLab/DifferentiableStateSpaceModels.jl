using DifferentiableStateSpaceModels, Symbolics, Test

# Setup for test
@testset "Basic construction of model" begin
    ∞ = Inf
    @variables α, β, ρ, δ, σ, Ω_1, Ω_2
    @variables t::Integer, k(..), z(..), c(..), q(..)

    x = [k, z]
    y = [c, q]
    p = [α, β, ρ, δ, σ]

    H = [1 / c(t) - (β / c(t + 1)) * (α * exp(z(t + 1)) * k(t + 1)^(α - 1) + (1 - δ)),
         c(t) + k(t + 1) - (1 - δ) * k(t) - q(t), q(t) - exp(z(t)) * k(t)^α,
         z(t + 1) - ρ * z(t)]

    steady_states = [k(∞) ~ (((1 / β) - 1 + δ) / α)^(1 / (α - 1)), z(∞) ~ 0,
                     c(∞) ~ (((1 / β) - 1 + δ) / α)^(α / (α - 1)) -
                            δ * (((1 / β) - 1 + δ) / α)^(1 / (α - 1)),
                     q(∞) ~ (((1 / β) - 1 + δ) / α)^(α / (α - 1))]

    steady_states_iv = [k(∞) ~ (((1 / β) - 1 + δ) / α)^(1 / (α - 1)), z(∞) ~ 0,
                        c(∞) ~ (((1 / β) - 1 + δ) / α)^(α / (α - 1)) -
                               δ * (((1 / β) - 1 + δ) / α)^(1 / (α - 1)),
                        q(∞) ~ (((1 / β) - 1 + δ) / α)^(α / (α - 1))]

    n_ϵ = 1
    n_z = 2
    n_x = length(x)
    n_y = length(y)
    n_p = length(p)
    Γ = reshape([σ], n_ϵ, n_ϵ)
    η = reshape([0; -1], n_x, n_ϵ) # η is n_x * n_ϵ matrix

    Q = zeros(n_z, n_x + n_y)
    Q[1, 1] = 1.0
    Q[2, 3] = 1.0

    Ω = [Ω_1, Ω_2]

    model_name = "rbc_temp"
    verbose = true
    save_ip = true
    save_oop = false
    max_order = 2
    skipzeros = false
    fillzeros = false

    module_cache_path = make_perturbation_model(H; model_name, t, y, x, p, steady_states,
                                                steady_states_iv, Γ, Ω, η, Q,
                                                overwrite_model_cache = true, verbose,
                                                max_order, save_ip, save_oop, skipzeros,
                                                fillzeros)
    make_perturbation_model(H; model_name, t, y, x, p, steady_states, steady_states_iv, Γ,
                            Ω, η, Q, overwrite_model_cache = true)

    make_perturbation_model(H; model_name, t, y, x, p, steady_states, steady_states_iv, Γ,
                            Ω, η, Q, overwrite_model_cache = false, verbose, max_order,
                            save_ip, save_oop, skipzeros, fillzeros)
    #module_cache_path = join_path(default_model_cache_location(), "rbc_temp.jl")

    # Load the module
    include(module_cache_path)

    # Test the construction
    m = PerturbationModel(Main.rbc_temp)
    @test m.n_y == n_y
    @test m.max_order == max_order
    @test m.mod.n_z == n_z
    # Note that this is inherently dynamic and cannot be inferred, so @inferred PerturbationModel(Main.rbc_observables) would fail

    c = SolverCache(m, Val(2), [:a, :b, :c]) # the exact symbol names won't matter for inference
    @inferred SolverCache(m, Val(2), [:a, :b, :c])
    @inferred SolverCache(m, Val(2), [:a])  # less differentiated parameters shouldn't matter
    @inferred SolverCache(m, Val(1), [:a, :b, :c])  # lower order shouldn't break inference either

    # Test convenience macro.  Use dict or named tuple of args, without the model_name
    args = (; t, y, x, p, steady_states, steady_states_iv, Γ, Ω, η, Q,
            overwrite_model_cache = true, verbose, max_order, save_ip, save_oop, skipzeros,
            fillzeros)

    # If the model_nam.jl doesn't exist it creates it, inclues it, and uses that to construct
    # If it exists but the modeule name isn't loaded, then it includes it and uses the constructor
    # If the module name is loaded, it just calls the constructor
    m = @make_and_include_perturbation_model("rbc_conv2", H, args)
    m = PerturbationModel(Main.rbc_conv2) # After creation, can use module constructor or simply call this function again 
    m = @make_and_include_perturbation_model("rbc_conv2", H, args) # should not regenerate
    m = @make_and_include_perturbation_model("rbc_temp", H, args) # already loaded from previous include
end

@testset "Loading Examples" begin
    # Tests for loading examples
    H, mod_vals = DifferentiableStateSpaceModels.Examples.rbc()
    make_perturbation_model(H; model_name = "rbc_test_2", overwrite_model_cache = true,
                            verbose = true, mod_vals...)

    m = @include_example_module(DifferentiableStateSpaceModels.Examples.rbc)
    @test m.n_y == 2

    # Load second one
    m = @include_example_module(DifferentiableStateSpaceModels.Examples.rbc)
    @test m.n_y == 2

    m = @include_example_module(DifferentiableStateSpaceModels.Examples.rbc_observables)
    @test m.n_y == 2

    m = @include_example_module(DifferentiableStateSpaceModels.Examples.rbc_solve_steady_state)
    @test m.n_y == 2

    m = @include_example_module(DifferentiableStateSpaceModels.Examples.rbc_multiple_shocks)
    @test m.n_y == 2

    m = @include_example_module(DifferentiableStateSpaceModels.Examples.rbc_observables_separate_variance)
    @test m.n_y == 2
end
