using DifferentiableStateSpaceModels, Symbolics, LinearAlgebra, Zygote, Test, Distributions
using DifferentiableStateSpaceModels.Examples
using DifferentiableStateSpaceModels: order_vector_by_symbols,
                                      fill_array_by_symbol_dispatch, all_fields_equal
using ChainRulesTestUtils
function test_first_order(p_d, p_f, m, c)
    sol = generate_perturbation(m, p_d, p_f, Val(1); cache = c)
    return sum(sol.y) + sum(sol.x) + sum(sol.A) + sum(sol.B) + sum(sol.C) +
           sum(cov(sol.D)) + sum(sol.x_ergodic.Σ.mat) + sum(sol.g_x)
end

# Some sort of inference issues.  Trouble putting in function and need the `const` for now
# see #117
# @testset "grad_tests" begin
const m_grad = @include_example_module(Examples.rbc_observables)  # const fixes current bug.  Can't move inside
const c = SolverCache(m_grad, Val(1), Val(2))
p_f = (ρ = 0.2, δ = 0.02, σ = 0.01, Ω_1 = 0.01)
p_d = (α = 0.5, β = 0.95)
p_d_2 = (α = 0.54, β = 0.92)
test_first_order(p_d, p_f, m_grad, c)
grad_test_1 = gradient((args...) -> test_first_order(args..., p_f, m_grad, c), p_d)
grad_test_2 = gradient((args...) -> test_first_order(args..., p_f, m_grad, c), p_d)
grad_test_3 = gradient((args...) -> test_first_order(args..., p_f, m_grad, c), p_d_2)

fd_eps = sqrt(eps())
@test grad_test_1[1].α ≈
      (test_first_order(merge(p_d, (; α = p_d.α + fd_eps)), p_f, m_grad, c) -
       test_first_order(merge(p_d, (; α = p_d.α - fd_eps)), p_f, m_grad, c)) / (2 * fd_eps)
@test grad_test_2[1].α ≈
      (test_first_order(merge(p_d, (; α = p_d.α + fd_eps)), p_f, m_grad, c) -
       test_first_order(merge(p_d, (; α = p_d.α - fd_eps)), p_f, m_grad, c)) / (2 * fd_eps)
@test grad_test_1[1].β ≈
      (test_first_order(merge(p_d, (; β = p_d.β + fd_eps)), p_f, m_grad, c) -
       test_first_order(merge(p_d, (; β = p_d.β - fd_eps)), p_f, m_grad, c)) / (2 * fd_eps)
@test grad_test_2[1].β ≈
      (test_first_order(merge(p_d, (; β = p_d.β + fd_eps)), p_f, m_grad, c) -
       test_first_order(merge(p_d, (; β = p_d.β - fd_eps)), p_f, m_grad, c)) / (2 * fd_eps)
@test grad_test_3[1].α ≈
      (test_first_order(merge(p_d_2, (; α = p_d_2.α + fd_eps)), p_f, m_grad, c) -
       test_first_order(merge(p_d_2, (; α = p_d_2.α - fd_eps)), p_f, m_grad, c)) /
      (2 * fd_eps)
@test grad_test_3[1].β ≈
      (test_first_order(merge(p_d_2, (; β = p_d_2.β + fd_eps)), p_f, m_grad, c) -
       test_first_order(merge(p_d_2, (; β = p_d_2.β - fd_eps)), p_f, m_grad, c)) /
      (2 * fd_eps)
@test grad_test_1[1].β ≈
      (test_first_order(merge(p_d, (; β = p_d.β + fd_eps)), p_f, m_grad, c) -
       test_first_order(merge(p_d, (; β = p_d.β - fd_eps)), p_f, m_grad, c)) / (2 * fd_eps)

# some weird CRTU behavior.
test_rrule(Zygote.ZygoteRuleConfig(),
           (args...) -> test_first_order(args..., p_f, m_grad, c), p_d;
           rrule_f = rrule_via_ad,
           check_inferred = false)
# end

# Slows things down and didn't catch any bugs
# # Example with more cross-gradients
# #@testset "grad_cross_tests" begin
# const m_grad_cross = @include_example_module(Examples.rbc_cross_gradients)  # const fixes current bug.  Can't move inside
# p_f = (; p_5 = 0.006011187045517736, p_6 = -0.001608593944662794)
# p_d = (p_1 = 0.47263197503003485, p_2 = 0.9456023867955173, p_3 = 0.18364072167470738,
#        p_4 = 0.01155202369386235)
# @test test_first_order(p_d, p_f, m_grad_cross) ≈ 77.13512914470338
# gradient((args...) -> test_first_order(args..., p_f, m_grad_cross), p_d)
# test_rrule(Zygote.ZygoteRuleConfig(),
#            (args...) -> test_first_order(args..., p_f, m_grad_cross), p_d;
#            rrule_f = rrule_via_ad,
#            check_inferred = false)

# # random tangents for p_5 and p_6 causing trouble with gradients
# # manually checking with crude central finitediff
# p_f = (;)
# p_d = (p_1 = 0.47263197503003485, p_2 = 0.9456023867955173, p_3 = 0.18364072167470738,
#        p_4 = 0.01155202369386235, p_5 = 0.006011187045517736, p_6 = -0.001608593944662794)
# all_grads = gradient((args...) -> test_first_order(args..., p_f, m_grad_cross), p_d)

# fd_eps = sqrt(eps())
# @test all_grads[1].p_5 ≈
#       (test_first_order(merge(p_d, (; p_5 = p_d.p_5 + fd_eps)), p_f, m_grad_cross) -
#        test_first_order(merge(p_d, (; p_5 = p_d.p_5 - fd_eps)), p_f, m_grad_cross)) /
#       (2 * fd_eps)
# @test all_grads[1].p_6 ≈
#       (test_first_order(merge(p_d, (; p_6 = p_d.p_6 + fd_eps)), p_f, m_grad_cross) -
#        test_first_order(merge(p_d, (; p_6 = p_d.p_6 - fd_eps)), p_f, m_grad_cross)) /
#       (2 * fd_eps) rtol = 1e-5
# #end
