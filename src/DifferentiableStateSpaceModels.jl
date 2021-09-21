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


export generate_perturbation_model, default_model_cache_location, PerturbationModel, SolverCache


# OLD STUFF
export FirstOrderPerturbationModel,
    generate_perturbation,
    dssm_evolution,
    dssm_volatility,
    dssm_observation,
    FirstOrderPerturbationSolution,
    FirstOrderSolverCache,
    @include_example_module,
    connect_markov_variables,
    save_first_order_module,
    SecondOrderPerturbationModel,
    SecondOrderPerturbationSolution,
    SecondOrderSolverCache,
    save_second_order_module,
    PerturbationSolverSettings,
    Examples,
    default_model_cache_location,
    allocate_cache,
    AbstractFirstOrderPerturbationModel,
    AbstractSecondOrderPerturbationModel,
    LTI,
    LTILikelihood,
    QTI,
    QTILikelihood,
    save_model_results

export solve # will be replaced by SciML soon

include("utils.jl")
include("symbolic_utils.jl")
include("generate_perturbation_model.jl")
include("types.jl")
#include("module_constructor.jl")
#include("perturbation.jl")
#include("sequence.jl")
#include("Examples/Examples.jl")

end # module
