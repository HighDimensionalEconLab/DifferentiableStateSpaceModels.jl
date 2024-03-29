module DifferentiableStateSpaceModels

using ChainRulesCore
using GeneralizedSylvesterSolver
using Logging
using Latexify
using LaTeXStrings
using LinearAlgebra
using MacroTools
using MatrixEquations
using NLsolve
using Parameters
using RecursiveFactorization # Or perhaps MKL will dominate this on all hardware.
using StructArrays
using Symbolics
using SymbolicUtils
using DifferenceEquations
# Flip off when not debugging
# using Infiltrator

export make_perturbation_model, default_model_cache_location, PerturbationModel,
       SolverCache, PerturbationSolverSettings, @include_example_module,
       @make_and_include_perturbation_model, generate_perturbation,
       generate_perturbation_derivatives!, model_H_latex, model_steady_states_iv_latex,
       model_steady_states_latex, verify_steady_state, irf

# export solve # will be replaced by SciML

include("utils.jl")
include("symbolic_utils.jl")
include("make_perturbation_model.jl")
include("types.jl")
include("generate_perturbation.jl")
include("generate_perturbation_derivatives.jl")
include("sequential.jl")
include("Examples/Examples.jl")

end # module
