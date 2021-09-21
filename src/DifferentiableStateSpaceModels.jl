module DifferentiableStateSpaceModels

using Logging
using MatrixEquations
using TensorCast
using DistributionsAD
using ChainRulesCore
using DocStringExtensions
using Distributions
using Parameters
using MacroTools
using LinearAlgebra
using Zygote
using NLsolve
using Symbolics
using SymbolicUtils
using Latexify
using LaTeXStrings
using StructArrays
using TimerOutputs


export make_perturbation_model, default_model_cache_location, PerturbationModel, SolverCache, @include_example_module

export solve # will be replaced by SciML

include("utils.jl")
include("symbolic_utils.jl")
include("make_perturbation_model.jl")
include("types.jl")
include("generate_perturbation.jl")
include("generate_perturbation_derivatives.jl")
#include("sequence.jl")
include("Examples/Examples.jl")

end # module
