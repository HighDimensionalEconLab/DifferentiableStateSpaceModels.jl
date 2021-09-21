

# Prefer to be inside the DifferentiableStateSpaceModels/Examples, but having trouble
# Must be a macro to avoid world-age issues
"""
Return model from a stored module, creating as necessary .  If the module is
- already included, it just returns the `Module`
- saved to disk but not already included, it does `include`
- not saved to disk, it saves it, includes it

The specification of a `model_generator` has a sepcification like `H, mod_vals = model_generator()`
where `H` is an MTK expression the `mod_vals` can be splatted into the `FirstOrderPerturbationModel` or `SecondOrderPerturbationModel` constructor.

The name of the stored `module` is the name of `model_generator` function and the order. 

$(SIGNATURES)

# Examples
```julia-repl
julia> m = @include_example_module(Examples.rbc_observables_benchmark)  # stores or loads `.function_cache/rbc_observables_1.jl` 

julia> m = @include_example_module(Examples.rbc_observables_benchmark, 2) # second order, loads `.function_cache/rbc_observables_2.jl` 
```
"""
macro include_example_module(
    model_generator,
    order = 1
)
    quote
        function_name = "$(Base.nameof($(model_generator)))"  # strips out module name
        local model_name = "$($model_generator)_$($order)"
        local model_filename = "$(model_name).jl"
        local save_module_function =
            ($order == 1) ? DifferentiableStateSpaceModels.save_first_order_module :
            DifferentiableStateSpaceModels.save_second_order_module
        local model_cache_path = joinpath(
            DifferentiableStateSpaceModels.default_model_cache_location(),
            model_filename,
        )
        if !isfile(model_cache_path)
            local H , mod_vals = $model_generator()
            save_module_function(
                H;
                model_name,
                mod_vals...,
            )
            include(model_cache_path)
        elseif !isdefined(Main, Symbol(model_name))
            include(model_cache_path)
        end
        local mod = getfield(Main, Symbol(model_name))
        ($order == 1) ? DifferentiableStateSpaceModels.FirstOrderPerturbationModel(mod) :
        DifferentiableStateSpaceModels.SecondOrderPerturbationModel(mod)
    end |> esc
end
