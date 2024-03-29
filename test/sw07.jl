using DifferentiableStateSpaceModels, Symbolics, LinearAlgebra, Test, Zygote, Statistics
using ChainRulesTestUtils
using FiniteDiff
using FiniteDiff: finite_difference_derivative, finite_difference_gradient,
                  finite_difference_jacobian, finite_difference_hessian

isdefined(Main, :SW07) || include(joinpath(pkgdir(DifferentiableStateSpaceModels),
                                           "test/generated_models/SW07.jl"))
const m_sw = PerturbationModel(Main.SW07)

# @testset "Dense SW First Order" begin
#     m = @include_example_module(Examples.SW07)

#     p = [10, 0.51, 10, 0, 0.7, 0.742, 0, 0, 0.24, 0.2696, 6.0144, 0.025, 1.5, 0.6361, 1.5,
#          0.3243, 0.8087, 0.47, 0.6, 1.9423, 1.5, 1.488, 0.2347, 0.0593, 0.8762, 0.9977,
#          0.5799, 0.9957, 0.7165, 0, 0, 0, 0.3982, 0.18]
#     generate_perturbation(m, p)
#     reset_timer!()
#     sol = generate_perturbation(m, p)
#     print_timer()
#     @test sol.retcode == :Success
# end
@testset "Dense SW First Order" begin
    p_d = (ε_w = 10.0, ρ_ga = 0.51, ε_p = 10.0, l_bar = 0.0, Π_bar = 0.7, B = 0.742,
           μ_w = 0.0,
           μ_p = 0.0, α = 0.24, ψ = 0.2696, φ = 6.0144, δ = 0.025, σ_c = 1.5, λ = 0.6361,
           ϕ_p = 1.5, ι_w = 0.3243, ξ_w = 0.8087, ι_p = 0.47, ξ_p = 0.6, σ_l = 1.9423,
           ϕ_w = 1.5, r_π = 1.488, r_Δy = 0.2347, r_y = 0.0593, ρ = 0.8762, ρ_a = 0.9977,
           ρ_b = 0.5799, ρ_g = 0.9957, ρ_i = 0.7165, ρ_r = 0.0, ρ_p = 0.0, ρ_w = 0.0,
           γ_bar = 0.3982, gy_ss = 0.18, se_a = 0.4618, se_b = 1.8513, se_g = 0.6090,
           se_i = 0.6017, se_m = 0.2397, se_π = 0.1455, se_w = 0.2089)
    p_f = (Ω_ii = sqrt(1e-5),)

    c = SolverCache(m_sw, Val(1), p_d)
    sol = generate_perturbation(m_sw, p_d, p_f; cache = c)
    # The following numbers have been compared with the Dynare results
    @test sol.C[:, 1:6] ≈
          [-0.28813568284395347 -0.7118643171560466 -0.0004453884466803 -0.17884070071357855 0.012233819522172222 0.48611470010516994;
           -0.28168902504971266 0.2816890250497125 -0.002016285719923244 0.04626339754117198 0.017034324795422214 -0.32512931888395524;
           -0.3665018568648009 0.36650185686480125 0.004360055105974288 0.003264911328657694 -0.00045589630688400444 -0.041922405285856536;
           -0.040739750405921434 0.04073975040592146 0.000771707665569567 0.005420514459266393 -0.0006905315403806023 0.01174656579355143;
           -0.057581376815551866 0.057581376815551984 0.0016828778070008797 -0.023261941823524417 -0.0022078632769032902 0.011646660084958956;
           0.15435193444138412 -0.1543519344413842 0.020207601741540456 -0.04757202750310869 -0.07053731527422884 0.11980533590670855;
           -0.1977437516488323 0.1977437516488322 -0.00046210279747159976 -0.1976714626296963 0.008849526885981281 0.3442945770530765]
end

# Check gradients
p_d = (ε_w = 10.0, ρ_ga = 0.51, ε_p = 10.0, l_bar = 0.0, Π_bar = 0.7, B = 0.742, μ_w = 0.0,
       μ_p = 0.0, α = 0.24, ψ = 0.2696, φ = 6.0144, δ = 0.025, σ_c = 1.5, λ = 0.6361,
       ϕ_p = 1.5, ι_w = 0.3243, ξ_w = 0.8087, ι_p = 0.47, ξ_p = 0.6, σ_l = 1.9423,
       ϕ_w = 1.5, r_π = 1.488, r_Δy = 0.2347, r_y = 0.0593, ρ = 0.8762, ρ_a = 0.9977,
       ρ_b = 0.5799, ρ_g = 0.9957, ρ_i = 0.7165, ρ_r = 0.0, ρ_p = 0.0, ρ_w = 0.0,
       γ_bar = 0.3982, gy_ss = 0.18, se_a = 0.4618, se_b = 1.8513, se_g = 0.6090,
       se_i = 0.6017, se_m = 0.2397, se_π = 0.1455, se_w = 0.2089)
p_f = (Ω_ii = sqrt(1e-5),)
settings = PerturbationSolverSettings(; tol_cholesky = 1e9)  # or zero
function test_first_order(p_d, p_f, m, settings)
    sol = generate_perturbation(m, p_d, p_f; settings)
    return (sum(sol.y) + sum(sol.x) + sum(sol.A) + sum(sol.C) + sum(sol.D) +
            sum(sol.g_x) + sum(sol.B))
end

test_first_order(p_d, p_f, m_sw, settings)
function test_first_order_closure(p_d)
    return test_first_order(p_d, p_f, m_sw, settings)
end
gradient(test_first_order_closure, p_d)

# Amount of failurs likely due to random projections being too large?  But not enough to mess up results
test_rrule(Zygote.ZygoteRuleConfig(), test_first_order_closure, p_d; rrule_f = rrule_via_ad,
           check_inferred = false, rtol = 1e-7)

#The ergodic distribution is causing us trouble...
function test_first_order_other(p_d, p_f, m, settings)
    sol = generate_perturbation(m, p_d, p_f; settings)
    return mean(sol.x_ergodic_var)  # These are causing trouble
end

# Numerical error on the ergodic distribution higher.
test_rrule(Zygote.ZygoteRuleConfig(),
           p_d -> test_first_order_other(p_d, p_f, m_sw, settings), p_d;
           rrule_f = rrule_via_ad,
           check_inferred = false, rtol = 1e-6)

grad_val = gradient(p_d -> test_first_order_other(p_d, p_f, m_sw, settings), p_d)

test_eps = sqrt(eps())
@test grad_val[1].α ≈
      (test_first_order_other(merge(p_d, (; α = p_d.α + test_eps)), p_f, m_sw, settings) -
       test_first_order_other(merge(p_d, (; α = p_d.α - test_eps)), p_f, m_sw, settings)) /
      (2 * test_eps) rtol = 1e-5  # PASSES BUT NOT ACCEPTABLE PRECISION.  MASKING ERROR?

@test grad_val[1].λ ≈
      (test_first_order_other(merge(p_d, (; λ = p_d.λ + test_eps)), p_f, m_sw, settings) -
       test_first_order_other(merge(p_d, (; λ = p_d.λ - test_eps)), p_f, m_sw, settings)) /
      (2 * test_eps) rtol = 1e-4 # PASSES BUT NOT ACCEPTABLE PRECISION.  MASKING ERROR?