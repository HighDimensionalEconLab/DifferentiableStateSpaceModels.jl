module DifferentiableStateSpaceModels

# Removed by Cameron Pfiffer, July 25th 2021.
# May need to be reincluded if it turns out that testing
# was insufficient.
# using RecursiveArrayTools
# using ForwardDiff
# using Turing
# using Dates
# using Tullio
# using TensorOperations

using Logging
using MatrixEquations
using GeneralizedGenerated
using TensorCast
using DistributionsAD
using ChainRulesCore
using DocStringExtensions
using Distributions
using ModelingToolkit
using Parameters
using MacroTools
using LinearAlgebra
using Zygote
using NLsolve
using SparseArrays
using TimerOutputs

using ModelingToolkit:
    build_function, hessian, SerialForm, MultithreadedForm, DistributedForm, Term

export FirstOrderPerturbationModel,
    DenseFunctions,
    SparseFunctions,
    generate_perturbation,
    dssm_evolution,
    dssm_volatility,
    dssm_observation,
    FirstOrderPerturbationSolution,
    FirstOrderSolverCache,
    @make_markov,
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
    get_threadsafe_cache,
    AbstractFirstOrderPerturbationModel,
    ThreadLocalCache,
    AbstractSecondOrderPerturbationModel,
    LTI,
    LTILikelihood,
    QTI,
    QTILikelihood,
    save_model_results

export solve # will be replaced by SciML soon

include("utils.jl")
include("mtk_utils.jl")
include("types.jl")
include("cache_utils.jl")
include("mtk_constructor.jl")
include("module_constructor.jl")
include("perturbation.jl")
include("sequence.jl")
include("Examples/Examples.jl")

end # module
