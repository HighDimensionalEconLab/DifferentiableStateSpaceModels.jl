#Benchmarking of FVGQ20 variants
#loads the files from the unit test folder, where they should be manually updated.
using DifferentiableStateSpaceModels, BenchmarkTools
using DifferentiableStateSpaceModels.Examples
isdefined(Main, :FVGQ20) || include(joinpath(pkgdir(DifferentiableStateSpaceModels),
                                             "test/generated_models/FVGQ20.jl"))
const FVGQ = BenchmarkGroup()
const FVGQ["fvgq"] = BenchmarkGroup()

# Prepare the model for the FVGQ20 observables
const m_fvgq = PerturbationModel(Main.FVGQ20)
# Basic Steady State
p_d_fvgq = (β = 0.998, h = 0.97, ϑ = 1.17, κ = 9.51, α = 0.21, θp = 0.82, χ = 0.63,
            γR = 0.77, γy = 0.19, γΠ = 1.29, Πbar = 1.01, ρd = 0.12, ρφ = 0.93, ρg = 0.95,
            g_bar = 0.3, σ_A = exp(-3.97), σ_d = exp(-1.51), σ_φ = exp(-2.36),
            σ_μ = exp(-5.43), σ_m = exp(-5.85), σ_g = exp(-3.0), Λμ = 3.4e-3, ΛA = 2.8e-3)
const p_f_fvgq = (δ = 0.025, ε = 10, ϕ = 0, γ2 = 0.001, Ω_ii = sqrt(1e-5))
const c_fvgq_1 = SolverCache(m_fvgq, Val(1), p_d_fvgq)
const c_fvgq_2 = SolverCache(m_fvgq, Val(2), p_d_fvgq)

# The derivatives caches need to be already filled in case they are done out of order
const c_fvgq_1_p = SolverCache(m_fvgq, Val(1), p_d_fvgq)
const c_fvgq_2_p = SolverCache(m_fvgq, Val(2), p_d_fvgq)
generate_perturbation(m_fvgq, p_d_fvgq, p_f_fvgq, Val(1); cache = c_fvgq_1_p)
generate_perturbation(m_fvgq, p_d_fvgq, p_f_fvgq, Val(2); cache = c_fvgq_2_p)

# Can't hurt to  call the derivatives to precompile as well, eventhough they will be written over
generate_perturbation_derivatives!(m_fvgq, p_d_fvgq, p_f_fvgq, c_fvgq_1_p)
generate_perturbation_derivatives!(m_fvgq, p_d_fvgq, p_f_fvgq, c_fvgq_2_p)

# fvgq
const FVGQ["fvgq"]["first_order"] = @benchmarkable generate_perturbation($m_fvgq, $p_d_fvgq,
                                                                         $p_f_fvgq, Val(1);
                                                                         cache = $c_fvgq_1)
const FVGQ["fvgq"]["first_order_p"] = @benchmarkable generate_perturbation_derivatives!($m_fvgq,
                                                                                        $p_d_fvgq,
                                                                                        $p_f_fvgq,
                                                                                        $c_fvgq_1_p)
const FVGQ["fvgq"]["second_order"] = @benchmarkable generate_perturbation($m_fvgq,
                                                                          $p_d_fvgq,
                                                                          $p_f_fvgq, Val(2);
                                                                          cache = $c_fvgq_2)
const FVGQ["fvgq"]["second_order_p"] = @benchmarkable generate_perturbation_derivatives!($m_fvgq,
                                                                                         $p_d_fvgq,
                                                                                         $p_f_fvgq,
                                                                                         $c_fvgq_2_p)
const FVGQ["fvgq"]["second_order_no_cache"] = @benchmarkable generate_perturbation($m_fvgq,
                                                                                   $p_d_fvgq,
                                                                                   $p_f_fvgq,
                                                                                   Val(2)) # just to see how long the cache takes

# return for the suite
FVGQ

# Uncomment to run on its own, otherwise nests into full benchmarking
# results = run(fvgq; verbose = true)
# results["fvgq"] # etc.
