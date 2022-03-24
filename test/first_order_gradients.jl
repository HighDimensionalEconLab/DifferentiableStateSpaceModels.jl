using DifferentiableStateSpaceModels, Symbolics, LinearAlgebra, Zygote, Test, Distributions
using DifferentiableStateSpaceModels.Examples
using DifferentiableStateSpaceModels: order_vector_by_symbols,
                                      fill_array_by_symbol_dispatch, all_fields_equal
using ChainRulesTestUtils
function test_first_order(p_d, p_f, m)
    sol = generate_perturbation(m, p_d, p_f)#, Val(1); cache = c) # manually passing in order
    return sum(sol.y) + sum(sol.x) + sum(sol.A) + sum(sol.B) + sum(sol.C) +
           sum(cov(sol.D)) + sum(sol.x_ergodic.Σ.mat) + sum(sol.g_x)
end

# Some sort of inference issues.  Trouble putting in function and need the `const` for now
# see #117
#@testset "grad_tests" begin
const m_grad = @include_example_module(Examples.rbc_observables)  # const fixes current bug.  Can't move inside
p_f = (ρ = 0.2, δ = 0.02, σ = 0.01, Ω_1 = 0.01)
p_d = (α = 0.5, β = 0.95)
test_first_order(p_d, p_f, m_grad)
gradient((args...) -> test_first_order(args..., p_f, m_grad), p_d)

test_rrule(Zygote.ZygoteRuleConfig(),
           (args...) -> test_first_order(args..., p_f, m_grad), p_d;
           rrule_f = rrule_via_ad,
           check_inferred = false)
#end