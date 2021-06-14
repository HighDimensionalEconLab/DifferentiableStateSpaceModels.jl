using DifferentiableStateSpaceModels
using Test, SparseArrays, ModelingToolkit, Parameters, LinearAlgebra
println("Running Testsuite with Threads.nthreads()=$(Threads.nthreads()) BLAS.vendor = $(BLAS.vendor())\n")
# See https://github.com/JuliaLang/julia/issues/33409
if (BLAS.vendor() == :openblas64)
    blas_num_threads = min(4, Int64(round(Sys.CPU_THREADS / 2)))  # even lower?
    println("Setting BLAS threads = $blas_num_threads")
    BLAS.set_num_threads(blas_num_threads)
end

include("sgu.jl")
include("sw07.jl")
