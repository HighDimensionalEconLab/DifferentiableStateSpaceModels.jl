using DifferentiableStateSpaceModels, Symbolics, LinearAlgebra, Zygote, Test, Distributions
using DifferentiableStateSpaceModels.Examples
using DifferentiableStateSpaceModels: order_vector_by_symbols,
                                      fill_array_by_symbol_dispatch, all_fields_equal
using ChainRulesTestUtils
function test_second_order(p_d, p_f, m)
    sol = generate_perturbation(m, p_d, p_f, Val(2))
    return sum(sol.y) + sum(sol.x) + sum(sol.A_0) + +sum(sol.A_1) + sum(sol.A_2) +
           sum(sol.B) + sum(sol.C_0) + sum(sol.C_1) + sum(sol.C_2) + sum(sol.D) +
           sum(sol.g_xx) + sum(sol.g_σσ) + sum(sol.g_x)
end

# Some sort of inference issues.  Trouble putting in function and need the `const` for now
# see #117
#@testset "grad_tests" begin
const m_grad_2 = @include_example_module(Examples.rbc_observables)  # const fixes current bug.  Can't move inside
p_f = (ρ = 0.2, δ = 0.02, σ = 0.01, Ω_1 = 0.01)
p_d = (α = 0.5, β = 0.95)
test_second_order(p_d, p_f, m_grad_2)
gradient((args...) -> test_second_order(args..., p_f, m_grad), p_d)

test_rrule(Zygote.ZygoteRuleConfig(),
           (args...) -> test_second_order(args..., p_f, m_grad_2), p_d;
           rrule_f = rrule_via_ad,
           check_inferred = false)

#end