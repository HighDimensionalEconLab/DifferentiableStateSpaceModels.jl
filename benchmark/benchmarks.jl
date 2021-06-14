using DifferentiableStateSpaceModels, BenchmarkTools
using Test, SparseArrays, ModelingToolkit, Parameters, LinearAlgebra, Random

println("Creating Benchmarking with Threads.nthreads()=$(Threads.nthreads()) BLAS.vendor = $(BLAS.vendor())\n")
# See https://github.com/JuliaLang/julia/issues/33409
if (BLAS.vendor() == :openblas64 && !haskey(ENV, "OPENBLAS_NUM_THREADS"))
    blas_num_threads = min(4, Int64(round(Sys.CPU_THREADS / 2)))  # or lower?
    println("Setting BLAS threads = $blas_num_threads")
    BLAS.set_num_threads(blas_num_threads)
end

# Setting miniumum number of evalations to avoid compilation
BenchmarkTools.DEFAULT_PARAMETERS.evals = 3  # at least 3

# Benchmark groups
const SUITE = BenchmarkGroup()
SUITE["utilities"] = include(pkgdir(DifferentiableStateSpaceModels) * "/benchmark/utilities.jl")
SUITE["rbc_likelihoods"] = include(pkgdir(DifferentiableStateSpaceModels) * "/benchmark/rbc_likelihoods.jl")
SUITE["models"] = include(pkgdir(DifferentiableStateSpaceModels) * "/benchmark/models.jl")


# See README.md
# Can execute with:
# results = run(SUITE; verbose = true)

#Then extract with things like:
# median(results["rbc_likelihoods"])
# results["rbc_likelihoods"]["joint"]["likelihood"]
# Or save stuff, make a change, and then use `judge(m1, m2)`
