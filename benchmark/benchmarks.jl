using DifferentiableStateSpaceModels, BenchmarkTools
using Test, Parameters, LinearAlgebra, Random

using DifferentiableStateSpaceModels
using Test, LinearAlgebra

# The BLAS threads is still an issue in Julia 1.7
# This has no effect with MKL
DifferentiableStateSpaceModels.set_blas_threads()

println(
    "Running Testsuite with Threads.nthreads() = $(Threads.nthreads()) BLAS.vendor = $(BLAS.vendor()), and BLAS.num_threads = $(BLAS.get_num_threads()) \n",
)

# Setting miniumum number of evalations to avoid compilation
BenchmarkTools.DEFAULT_PARAMETERS.evals = 3  # at least 3

# Benchmark groups
const SUITE = BenchmarkGroup()
SUITE["utilities"] =
    include(pkgdir(DifferentiableStateSpaceModels) * "/benchmark/utilities.jl")
SUITE["rbc"] = include(pkgdir(DifferentiableStateSpaceModels) * "/benchmark/rbc.jl")
#SUITE["rbc_likelihoods"] =
    # include(pkgdir(DifferentiableStateSpaceModels) * "/benchmark/rbc_likelihoods.jl")

# See README.md
# Can execute with:
# results = run(SUITE; verbose = true)

#Then extract with things like:
# median(results["rbc_likelihoods"])
# results["rbc_likelihoods"]["joint"]["likelihood"]
# Or save stuff, make a change, and then use `judge(m1, m2)`
