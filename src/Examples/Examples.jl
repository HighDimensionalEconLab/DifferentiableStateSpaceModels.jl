# Examples used in docs, tests, and benchmarking
module Examples
using Parameters, ModelingToolkit
using DifferentiableStateSpaceModels
using DifferentiableStateSpaceModels: default_model_cache_location
using Base: isdefined, joinpath, isfile

include("rbc.jl")
include("sgusmallopen.jl")
include("SW07.jl")
#include("krusellsmith.jl")
end

