using DifferentiableStateSpaceModels
using Test, LinearAlgebra

# Delete the .function_cache
rm(default_model_cache_location(); force = true, recursive = true)

@time include("make_perturbation_model.jl")
@time include("first_order_perturbation.jl")
@time include("second_order_perturbation.jl")
@time include("first_order_sequence.jl")
get(ENV, "CI", "false") == "false" && @time include("second_order_sequence.jl")
get(ENV, "CI", "false") == "false" && @time include("sgu.jl")
get(ENV, "CI", "false") == "false" && @time include("FVGQ20.jl")
@time include("symbolic_utils.jl")
@time include("utils.jl")