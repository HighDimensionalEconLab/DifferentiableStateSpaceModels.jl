using MKL
using DifferentiableStateSpaceModels, BenchmarkTools
using Test, Parameters, LinearAlgebra, Random

using DifferentiableStateSpaceModels
using Test, LinearAlgebra

# The BLAS threads is still an issue in Julia 1.7
# This has no effect with MKL
# See https://github.com/JuliaLang/julia/issues/33409
# Default even lower?
function set_blas_threads(openblas_threads = min(4, Int64(round(Sys.CPU_THREADS / 2))))
    if (BLAS.vendor() == :openblas64)
        println("Setting openblas64 threads = $openblas_threads")
        BLAS.set_num_threads(openblas_threads)
    end
end
set_blas_threads()

println("Running Testsuite with Threads.nthreads() = $(Threads.nthreads()) BLAS.vendor = $(BLAS.vendor()), and BLAS.num_threads = $(BLAS.get_num_threads()) \n")

# Setting miniumum number of evalations to avoid compilation
BenchmarkTools.DEFAULT_PARAMETERS.evals = 3  # at least 3

# Benchmark groups
const SUITE = BenchmarkGroup()
SUITE["utilities"] = include(pkgdir(DifferentiableStateSpaceModels) *
                             "/benchmark/utilities.jl")
SUITE["rbc"] = include(pkgdir(DifferentiableStateSpaceModels) * "/benchmark/rbc.jl")
SUITE["sgu"] = include(pkgdir(DifferentiableStateSpaceModels) * "/benchmark/sgu.jl")
SUITE["fvgq"] = include(pkgdir(DifferentiableStateSpaceModels) * "/benchmark/fvgq.jl")
#SUITE["rbc_likelihoods"] =
# include(pkgdir(DifferentiableStateSpaceModels) * "/benchmark/rbc_likelihoods.jl")

# See README.md
# Can execute with:
# results = run(SUITE; verbose = true)

#Then extract with things like:
# median(results["rbc_likelihoods"])
# results["rbc_likelihoods"]["joint"]["likelihood"]
# Or save stuff, make a change, and then use `judge(m1, m2)`
