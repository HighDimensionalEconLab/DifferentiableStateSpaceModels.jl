using DifferentiableStateSpaceModels, Symbolics, LinearAlgebra, Test, Zygote
using ChainRulesTestUtils
using DifferentiableStateSpaceModels.Examples

@testset "SGUext First Order" begin
    m = @include_example_module(Examples.sguext)
    p_d = (γ = 2.0, ω = 1.455, ρ = 0.42, σe = 0.0129, δ = 0.1, ψ = 0.000742, α = 0.32,
           ϕ = 0.028, β = 1.0 / (1.0 + 0.04), r_w = 0.04, d_bar = 0.7442, ρ_u = 0.2,
           σu = 0.003, ρ_v = 0.4, σv = 0.1)
    p_f = nothing

    c = SolverCache(m, Val(1), p_d)
    settings = PerturbationSolverSettings(; print_level = 0)
    sol = generate_perturbation(m, p_d, p_f; cache = c, settings)
    generate_perturbation_derivatives!(m, p_d, p_f, c)

    # TODO check that these are actually the correct results, currently only a regression test
    # Also, does not check the derivative details in the cache
    @test sol.retcode == :Success
end

# No D
function test_first_order(p_d, p_f, m)
    sol = generate_perturbation(m, p_d, p_f)#, Val(1); cache = c) # manually passing in order
    return sum(sol.y) + sum(sol.x) + sum(sol.A) + sum(sol.B) + sum(sol.C) +
           sum(sol.x_ergodic.Σ.mat)
end
# Gradients.  Can't put in a testset until #117 fixed
#@testset "SGU 1st order Gradients" begin
const m_sgu = @include_example_module(Examples.sguext)
p_d = (γ = 2.0, ω = 1.455, ρ = 0.42, σe = 0.0129, δ = 0.1, ψ = 0.000742, α = 0.32,
       ϕ = 0.028, β = 1.0 / (1.0 + 0.04), r_w = 0.04, d_bar = 0.7442, ρ_u = 0.2,
       σu = 0.003, ρ_v = 0.4, σv = 0.1)
p_f = nothing

test_first_order(p_d, p_f, m_sgu)
gradient((args...) -> test_first_order(args..., p_f, m_sgu), p_d)

test_rrule(Zygote.ZygoteRuleConfig(),
           (args...) -> test_first_order(args..., p_f, m_sgu), p_d;
           rrule_f = rrule_via_ad,
           check_inferred = false, rtol = 1e-7)  # note the rtol is not the default, but it is good enough

function test_second_order_no_D(p_d, p_f, m)
    sol = generate_perturbation(m, p_d, p_f, Val(2))
    return sum(sol.y) + sum(sol.x) + sum(sol.A_0) + +sum(sol.A_1) + sum(sol.A_2) +
           sum(sol.B) + sum(sol.C_0) + sum(sol.C_1) + sum(sol.C_2) +
           sum(sol.g_xx) + sum(sol.g_σσ) + sum(sol.g_x)
end
test_second_order_no_D(p_d, p_f, m_sgu)
gradient((args...) -> test_second_order_no_D(args..., p_f, m_sgu), p_d)

test_rrule(Zygote.ZygoteRuleConfig(),
           (args...) -> test_second_order_no_D(args..., p_f, m_sgu), p_d;
           rrule_f = rrule_via_ad,
           check_inferred = false, rtol = 1e-8)  # note the rtol is not the default, but it is good enough
# end

@testset "SGUext Second Order" begin
    m = @include_example_module(Examples.sguext)
    p_d = (γ = 2.0, ω = 1.455, ρ = 0.42, σe = 0.0129, δ = 0.1, ψ = 0.000742, α = 0.32,
           ϕ = 0.028, β = 1.0 / (1.0 + 0.04), r_w = 0.04, d_bar = 0.7442, ρ_u = 0.2,
           σu = 0.003, ρ_v = 0.4, σv = 0.1)
    p_f = nothing

    c = SolverCache(m, Val(2), p_d)
    sol = generate_perturbation(m, p_d, p_f, Val(2); cache = c)
    generate_perturbation_derivatives!(m, p_d, p_f, c)  # Solves and fills the cache

    # TODO check that these are actually the correct results, currently only a regression test
    # Also, does not check the derivative details in the cache
    @test sol.retcode == :Success
end