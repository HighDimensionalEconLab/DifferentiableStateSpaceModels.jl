using DifferentiableStateSpaceModels
using Test, LinearAlgebra

# The BLAS threads is still an issue in Julia 1.7
# This has no effect with MKL
DifferentiableStateSpaceModels.set_blas_threads()

println(
   "Running Testsuite with Threads.nthreads() = $(Threads.nthreads()) BLAS.vendor = $(BLAS.vendor()), and BLAS.num_threads = $(BLAS.get_num_threads()) \n",
)

# Delete the .function_cache
# e.g. ENV["DSSM_TEST_DELETE_CACHE"] = "false" environment variable to turn off, can be global
get(ENV, "DSSM_TEST_DELETE_CACHE", "true") == "true" &&
   rm(default_model_cache_location(), force = true, recursive = true)

include("make_perturbation_model.jl")
include("first_order_perturbation.jl")
include("second_order_perturbation.jl")
include("first_order_sequence.jl")
include("second_order_sequence.jl")
# include("rbc_estimation.jl")
include("sgu.jl")
include("FVGQ20.jl")
include("symbolic_utils.jl")
include("utils.jl")
