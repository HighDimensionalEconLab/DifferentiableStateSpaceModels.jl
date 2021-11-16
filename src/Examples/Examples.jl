# Examples used in docs, tests, and benchmarking
module Examples
using Parameters, Symbolics
using DifferentiableStateSpaceModels

include("rbc.jl")
include("sgu.jl")
include("FVGQ20.jl")
end

# Macro to load example can't be in that module.  Needs to be macro due to world-age issues

"""
Return model from a stored module, creating as necessary .  If the module is
- already included, it just returns the `Module`
- saved to disk but not already included, it does `include`
- not saved to disk, it saves it, includes it

The specification of a `model_generator` has a sepcification like `H, mod_vals = model_generator()`
where `H` is an symbolics expression the `mod_vals` can be splatted into the `PerturbationModel` constructor.

The name of the stored `module` is the name of `model_generator` function and the order. 

$(SIGNATURES)

# Examples
```julia-repl
julia> m = @include_example_module(Examples.rbc_observables_benchmark)  # stores or loads `.function_cache/rbc_observables.jl` 
"""
macro include_example_module(model_generator)
    return quote
        function_name = "$(Base.nameof($(model_generator)))"  # strips out module name
        local model_name = "$($model_generator)"
        local model_filename = "$(model_name).jl"
        local model_cache_path = joinpath(DifferentiableStateSpaceModels.default_model_cache_location(),
                                          model_filename)
        if !isfile(model_cache_path)
            local H, mod_vals = $model_generator()
            make_perturbation_model(H; model_name, mod_vals...)
            include(model_cache_path)
        elseif !isdefined(Main, Symbol(model_name))
            include(model_cache_path)
        end
        local mod = getfield(Main, Symbol(model_name))
        DifferentiableStateSpaceModels.PerturbationModel(mod)
    end |> esc
end
