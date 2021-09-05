using DifferentiableStateSpaceModels
# using Test, SparseArrays, ModelingToolkit, Parameters, LinearAlgebra
using Test, LinearAlgebra, Parameters, ModelingToolkit, SparseArrays, TimerOutputs

println(
    "Running Testsuite with Threads.nthreads() = $(Threads.nthreads()) BLAS.vendor = $(BLAS.vendor())\n",
)
# See https://github.com/JuliaLang/julia/issues/33409
if (BLAS.vendor() == :openblas64)
    blas_num_threads = min(4, Int64(round(Sys.CPU_THREADS / 2)))  # even lower?
    println("Setting BLAS threads = $blas_num_threads")
    BLAS.set_num_threads(blas_num_threads)
end

# Delete the .function_cache
# e.g. ENV["DSSM_TEST_DELETE_CACHE"] = "false" environment variable to turn off, can be global
get(ENV, "DSSM_TEST_DELETE_CACHE", "true") == "true" &&
    rm(default_model_cache_location(), force = true, recursive = true)


include("first_order_dense.jl")
#include("mtk_generation.jl")
include("second_order_dense.jl")
include("first_order_sequence.jl")
include("second_order_sequence.jl")
include("rbc_estimation.jl")
#include("mtk_utils.jl")
include("utils.jl")
