using DifferentiableStateSpaceModels, BenchmarkTools, Zygote, Random

const MODELS = BenchmarkGroup()

MODELS["rbc"] = BenchmarkGroup()
MODELS["SGU"] = BenchmarkGroup()
MODELS["SW07"] = BenchmarkGroup()
MODELS["generation"] = BenchmarkGroup()


## RBC variations
const p_f = [0.2, 0.02, 0.01]
const p = [0.5, 0.95]
const p_all = [p; p_f]
const settings_anderson = PerturbationSolverSettings(nlsolve_method = :anderson)

# Baseline, reallocating the cache each time.
rbc_model = @include_example_module(Examples.rbc, 1, DenseFunctions())
MODELS["rbc"]["rbc_model"] = @benchmarkable generate_perturbation($rbc_model, $p; p_f = $p_f)

rbc_second_model = @include_example_module(Examples.rbc, 2, DenseFunctions())
MODELS["rbc"]["rbc_second_model"] = @benchmarkable generate_perturbation($rbc_second_model, $p;
                                                                p_f = $p_f)
#TODO: try without cache reallocation

# Sparse
rbc_sparse_model = @include_example_module(Examples.rbc, 1, SparseFunctions())
MODELS["rbc"]["rbc_sparse_model"] = @benchmarkable generate_perturbation($rbc_sparse_model, $p;
                                                                p_f = $p_f)

# Steady State
rbc_steady_state_model = @include_example_module(Examples.rbc_solve_steady_state, 1, DenseFunctions())
MODELS["rbc"]["rbc_steady_state_model"] = @benchmarkable generate_perturbation($rbc_steady_state_model, $p;
                                                                p_f = $p_f)
MODELS["rbc"]["rbc_steady_state_model_anderson"] = @benchmarkable generate_perturbation($rbc_steady_state_model,
                                                                              $p;
                                                                              p_f = $p_f,
                                                                              settings = $settings_anderson)

# Multiple Shocks
rbc_empty_p_f_multiple_shocks_model = @include_example_module(Examples.rbc_empty_p_f_multiple_shocks, 1, DenseFunctions())
MODELS["rbc"]["rbc_empty_p_f_multiple_shocks_model"] = @benchmarkable generate_perturbation($rbc_empty_p_f_multiple_shocks_model,
                                                                $p_all; p_f = nothing)

# SGU parameters
const p_sgu = [2.0, 1.455, 0.42, 0.0129, 0.1, 0.000742, 0.32, 0.028, 1.0 / (1.0 + 0.04),
               0.04, 0.7442]   #From Cesa-Bianchi (2012)
const p_f_sgu = []
sgu_model = @include_example_module(Examples.sgusmallopen, 1, DenseFunctions())
MODELS["SGU"]["sgu_model"] = @benchmarkable generate_perturbation($sgu_model, $p_sgu;
                                                                p_f = $p_f_sgu)

sgu_sparse_model = @include_example_module(Examples.sgusmallopen, 1, SparseFunctions())
MODELS["SGU"]["sgu_sparse_model"] = @benchmarkable generate_perturbation($sgu_sparse_model, $p_sgu;
                                                                p_f = $p_f_sgu)

sgu_second_model = @include_example_module(Examples.sgusmallopen, 2, DenseFunctions())
MODELS["SGU"]["sgu_second_model"] = @benchmarkable generate_perturbation($sgu_second_model, $p_sgu;
                                                                p_f = $p_f_sgu)

## Smets and Wouters 2007
const p_sw = [10, 0.51, 10, 0, 0.7, 0.742, 0, 0, 0.24, 0.2696, 6.0144, 0.025, 1.5, 0.6361,
              1.5, 0.3243, 0.8087, 0.47, 0.6, 1.9423, 1.5, 1.488, 0.2347, 0.0593, 0.8762,
              0.9977, 0.5799, 0.9957, 0.7165, 0, 0, 0, 0.3982, 0.18]
const p_f_sw = []

# Baseline
sw07_model = @include_example_module(Examples.SW07, 1, DenseFunctions())
MODELS["SW07"]["sw07_model"] = @benchmarkable generate_perturbation($sw07_model, $p_sw;
                                                                 p_f = $p_f_sw)

# Sparse
sw07_sparse_model = @include_example_module(Examples.SW07, 1, SparseFunctions())
MODELS["SW07"]["sw07_sparse_model"] = @benchmarkable generate_perturbation($sw07_sparse_model, $p_sw;
                                                                 p_f = $p_f_sw)

# Second-order
sw07_second_model = @include_example_module(Examples.SW07, 2, DenseFunctions())
MODELS["SW07"]["sw07_second_model"] = @benchmarkable generate_perturbation($sw07_second_model, $p_sw;
                                                                 p_f = $p_f_sw)
## KS parameters
# const p_ks = [0.95, 3.0, 1.0/3.0, 0.2, 0.5, 0.95]  # β, γ, α, δ, sigs, ρz
# const p_ks = [1.0, 0.0, 2.0] #smoother, zlo, zhi, wlo, whi

## Model Creation Benchmarking
H_gen_rbc, mod_vals_gen_rbc = Examples.rbc()
MODELS["generation"]["rbc"] = @benchmarkable FirstOrderPerturbationModel($H_gen_rbc;
                                                                                       $mod_vals_gen_rbc...)



# Return the group for the suite
MODELS
