#Benchmarking of SGU variants
using DifferentiableStateSpaceModels, BenchmarkTools
using DifferentiableStateSpaceModels.Examples
const SGU = BenchmarkGroup()
const SGU["sgu"] = BenchmarkGroup()

# Prepare the model for the sgu observables
const m_sgu = @include_example_module(Examples.sgu)
# Basic Steady State
const p_f_sgu = nothing
const p_d_sgu = (γ = 2.0, ω = 1.455, ρ = 0.42, σe = 0.0129, δ = 0.1, ψ = 0.000742, α = 0.32,
                 ϕ = 0.028, β = 1.0 / (1.0 + 0.04), r_w = 0.04, d_bar = 0.7442)
const c_sgu_1 = SolverCache(m_sgu, Val(1), p_d_sgu)
const c_sgu_2 = SolverCache(m_sgu, Val(2), p_d_sgu)

# The derivatives caches need to be already filled in case they are done out of order
const c_sgu_1_p = SolverCache(m_sgu, Val(1), p_d_sgu)
const c_sgu_2_p = SolverCache(m_sgu, Val(2), p_d_sgu)
generate_perturbation(m_sgu, p_d_sgu, p_f_sgu, Val(1); cache = c_sgu_1_p)
generate_perturbation(m_sgu, p_d_sgu, p_f_sgu, Val(2); cache = c_sgu_2_p)

# Can't hurt to  call the derivatives to precompile as well, eventhough they will be written over
generate_perturbation_derivatives!(m_sgu, p_d_sgu, p_f_sgu, c_sgu_1_p)
generate_perturbation_derivatives!(m_sgu, p_d_sgu, p_f_sgu, c_sgu_2_p)

# sgu
const SGU["sgu"]["first_order"] = @benchmarkable generate_perturbation($m_sgu, $p_d_sgu,
                                                                       $p_f_sgu, Val(1);
                                                                       cache = $c_sgu_1)
const SGU["sgu"]["first_order_p"] = @benchmarkable generate_perturbation_derivatives!($m_sgu,
                                                                                      $p_d_sgu,
                                                                                      $p_f_sgu,
                                                                                      $c_sgu_1_p)
const SGU["sgu"]["second_order"] = @benchmarkable generate_perturbation($m_sgu, $p_d_sgu,
                                                                        $p_f_sgu, Val(2);
                                                                        cache = $c_sgu_2)
const SGU["sgu"]["second_order_p"] = @benchmarkable generate_perturbation_derivatives!($m_sgu,
                                                                                       $p_d_sgu,
                                                                                       $p_f_sgu,
                                                                                       $c_sgu_2_p)
const SGU["sgu"]["second_order_no_cache"] = @benchmarkable generate_perturbation($m_sgu,
                                                                                 $p_d_sgu,
                                                                                 $p_f_sgu,
                                                                                 Val(2)) # just to see how long the cache takes
const SGU["sgu"]["first_order_solver_cache"] = @benchmarkable SolverCache($m_sgu, Val(1),
                                                                          $p_d_sgu)
const SGU["sgu"]["second_order_solver_cache"] = @benchmarkable SolverCache($m_sgu, Val(2),
                                                                           $p_d_sgu)

# return for the suite
SGU

# Uncomment to run on its own, otherwise nests into full benchmarking
# results = run(SGU; verbose = true)
# results["sgu"] # etc.
