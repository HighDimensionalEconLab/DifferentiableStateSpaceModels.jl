#Benchmarking of RBC variants
using DifferentiableStateSpaceModels, BenchmarkTools
using DifferentiableStateSpaceModels.Examples
const RBC = BenchmarkGroup()
const RBC["rbc_observables"] = BenchmarkGroup()

# Prepare the model for the RBC observables
const m_rbc_observables = @include_example_module(Examples.rbc_observables)
# Basic Steady State
const p_f_rbc_observables = (ρ=0.2, δ=0.02, σ=0.01, Ω_1=0.01)
const p_d_rbc_observables = (α=0.5, β=0.95)
const c_rbc_observables_1 = SolverCache(m_rbc_observables, Val(1), p_d_rbc_observables)
const c_rbc_observables_2 = SolverCache(m_rbc_observables, Val(2), p_d_rbc_observables)

# The derivatives caches need to be already filled in case they are done out of order
const c_rbc_observables_1_p = SolverCache(m_rbc_observables, Val(1), p_d_rbc_observables)
const c_rbc_observables_2_p = SolverCache(m_rbc_observables, Val(2), p_d_rbc_observables)
generate_perturbation(m_rbc_observables, p_d_rbc_observables, p_f_rbc_observables, Val(1); cache = c_rbc_observables_1_p)
generate_perturbation(m_rbc_observables, p_d_rbc_observables, p_f_rbc_observables, Val(2); cache = c_rbc_observables_2_p)

# rbc_observables
const RBC["rbc_observables"]["first_order"] =
    @benchmarkable generate_perturbation($m_rbc_observables, $p_d_rbc_observables, $p_f_rbc_observables, Val(1); cache = $c_rbc_observables_1)
const RBC["rbc_observables"]["first_order_p"] =
    @benchmarkable generate_perturbation_derivatives!($m_rbc_observables, $p_d_rbc_observables, $p_f_rbc_observables, $c_rbc_observables_1_p)
const RBC["rbc_observables"]["second_order"] =
    @benchmarkable generate_perturbation($m_rbc_observables, $p_d_rbc_observables, $p_f_rbc_observables, Val(2); cache = $c_rbc_observables_2)
const RBC["rbc_observables"]["second_order_p"] =
    @benchmarkable generate_perturbation_derivatives!($m_rbc_observables, $p_d_rbc_observables, $p_f_rbc_observables, $c_rbc_observables_2_p)
const RBC["rbc_observables"]["second_order_no_cache"] =
    @benchmarkable generate_perturbation($m_rbc_observables, $p_d_rbc_observables, $p_f_rbc_observables, Val(2)) # just to see how long the cache takes
const RBC["rbc_observables"]["first_order_solver_cache"] = @benchmarkable SolverCache($m_rbc_observables, Val(1), $p_d_rbc_observables)
const RBC["rbc_observables"]["second_order_solver_cache"] = @benchmarkable SolverCache($m_rbc_observables, Val(2), $p_d_rbc_observables)

# return for the suite
RBC

# Uncomment to run on its own, otherwise nests into full benchmarking
# results = run(RBC; verbose = true)
# results["rbc_observables"] # etc.
