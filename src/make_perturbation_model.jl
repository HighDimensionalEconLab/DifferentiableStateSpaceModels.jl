# Helpers
function default_model_cache_location()
    return joinpath(pkgdir(DifferentiableStateSpaceModels), ".function_cache")
end

# TODO Add docstring
function make_perturbation_model(
    H;
    t,
    y,
    x,
    steady_states = nothing,
    steady_states_iv = nothing,    
    Γ,
    Ω = nothing,
    η,
    Q = I,
    p = nothing,
    model_name,
    model_cache_location = default_model_cache_location(),
    overwrite_model_cache = false,
    verbose = true,
    max_order = 2,
    save_ip = true,
    save_oop = false, # only does inplace by default
    skipzeros = false,
    fillzeros = false
)
    @assert max_order ∈ [1,2]
    @assert save_ip || save_oop
    @assert skipzeros == false # currently broken in symbolics otherwise?

    module_cache_path = joinpath(model_cache_location, model_name * ".jl")

    # only load cache if the module isn't already loaded in memory
    if (isdefined(Main, Symbol(model_name)) && !overwrite_model_cache)
        verbose && printstyled("Using existing module $model_name\n", color = :cyan)
        return module_cache_path = module_cache_path
    end

    # if path already exists
    if (ispath(module_cache_path) && !overwrite_model_cache)
        # path exists and not overwriting
        verbose &&
            printstyled("Model already generated at $module_cache_path\n", color = :cyan)
        return module_cache_path    
    end        

    n_y = length(y)
    n_x = length(x)
    n = n_y + n_x
    n_p = length(p)
    @assert n_p > 0  # code written to have at least one parameter
    n_ϵ = size(η, 2)
    n_z = (Q == I) ? n : size(Q, 1)

    # TODO: error check that p, y, x has no overlap
    # Get the markovian variables and create substitutions
    y_subs = StructArray(make_substitutions.(t, y))
    x_subs = StructArray(make_substitutions.(t, x))
    y = y_subs.var
    x = x_subs.var
    y_p = y_subs.var_p
    x_p = x_subs.var_p
    y_ss = y_subs.var_ss
    x_ss = x_subs.var_ss
    subs = vcat(x_subs, y_subs)
    all_to_markov = vcat(subs.markov_t, subs.markov_tp1, subs.markov_inf)
    all_to_var = vcat(subs.tp1_to_var, subs.inf_to_var)
    # Helper to take the [z(∞) ~ expr] and become [z => expr] after substitutions
    equations_to_dict(equations) = Dict(Symbol(substitute(substitute(eq.lhs, all_to_markov), all_to_var)) => Num(substitute(eq.rhs, all_to_markov)) for eq in equations)

    # create functions in correct order
    # TODO: Errors if missing any in expressions
    ȳ = arrange_vector_from_symbols(equations_to_dict(steady_states), y_subs.symbol)
    x̄ = arrange_vector_from_symbols(equations_to_dict(steady_states), x_subs.symbol)
    ȳ_iv = arrange_vector_from_symbols(equations_to_dict(steady_states_iv), y_subs.symbol)
    x̄_iv = arrange_vector_from_symbols(equations_to_dict(steady_states_iv), x_subs.symbol)

    # Get any latex generated stuff we wish for pretty display of the model
    H_latex = latexify(H)
    steady_states_latex = latexify(steady_states)
    steady_states_iv_latex = latexify(steady_states_iv)

    # steady state requiers differentiation after substitution, and wrt [y; x]
    H = substitute.(H, Ref(all_to_markov))
    H̄ = deepcopy(H)
    H̄ = substitute_and_simplify(H̄, all_to_var)

    # Derivatives utilities return nothing if either argument nothing
    H̄_w = nested_differentiate(H̄, [y; x]) # differentiate post-substitution wrt w = [y;x], force a dense one
    H_yp = nested_differentiate(H, y_p)
    H_y = nested_differentiate(H, y)
    H_xp = nested_differentiate(H, x_p)
    H_x = nested_differentiate(H, x)
    Ψ = [Symbolics.hessian(f, [y_p; y; x_p; x]; simplify=true) for f in H]  # need for 1st order derivatives
    Ψ_yp = (max_order < 2) ? nothing : nested_differentiate(Ψ, y_p)
    Ψ_y = (max_order < 2) ? nothing : nested_differentiate(Ψ, y)
    Ψ_xp = (max_order < 2) ? nothing : nested_differentiate(Ψ, x_p)
    Ψ_x = (max_order < 2) ? nothing : nested_differentiate(Ψ, x)


    # The parameter derivatives are maps for dispatching by Symbol
    # utility function substitutes/simplifies because these aren't themselves differentiated
    differentiate_to_dict(f, p) =  Dict([Symbol(p_val) => substitute_and_simplify(nested_differentiate(f, p_val), all_to_var) for p_val in p])
    H_p = differentiate_to_dict(H, p)
    Γ_p = differentiate_to_dict(Γ, p)
    Ω_p = differentiate_to_dict(Ω, p)
    ȳ_p = differentiate_to_dict(ȳ, p)
    x̄_p = differentiate_to_dict(x̄, p)
    H_yp_p = differentiate_to_dict(H_yp, p)
    H_xp_p = differentiate_to_dict(H_xp, p)
    H_y_p = differentiate_to_dict(H_y, p)
    H_x_p = differentiate_to_dict(H_x, p)
    Ψ_p = (max_order < 2) ? nothing : differentiate_to_dict(Ψ, p)

    # apply substitutions and simplify.  nothing stays nothing
    H = substitute_and_simplify(H, all_to_markov)
    H_yp = substitute_and_simplify(H_yp, all_to_var)
    H_xp = substitute_and_simplify(H_xp, all_to_var)
    H_x = substitute_and_simplify(H_x, all_to_var)
    H_y = substitute_and_simplify(H_y, all_to_var)
    Ψ = substitute_and_simplify(Ψ, all_to_var)
    Ψ_yp = substitute_and_simplify(Ψ_yp, all_to_var)
    Ψ_y = substitute_and_simplify(Ψ_y, all_to_var)
    Ψ_xp = substitute_and_simplify(Ψ_xp, all_to_var)
    Ψ_x = substitute_and_simplify(Ψ_x, all_to_var)

    # Generate all functions, and rename using utility
    function build_named_function(f, name, args...; symbol_dispatch = nothing)
        expr = isnothing(symbol_dispatch) ? build_function(f, args...;skipzeros,fillzeros) : build_function(f, nothing, args...;skipzeros,fillzeros)
        ip_function = name_symbolics_function(expr[1], Symbol(name); inplace = false, symbol_dispatch)
        oop_function = name_symbolics_function(expr[2], Symbol(name*"!"); inplace = true, symbol_dispatch)
        return ip_function, oop_function
    end
    build_named_function(::Nothing, name, args...; symbol_dispatch = nothing) = nothing

    Γ_expr = build_named_function(Γ, "Γ", p)
    Ω_expr = build_named_function(Ω, "Ω", p)
    H_expr = build_named_function(H, "H", y_p, y, y_ss, x_p, x, x_ss, p)
    H_yp_expr = build_named_function(H_yp, "H_yp", y, x, p)
    H_y_expr = build_named_function(H_y, "H_y", y, x, p)
    H_xp_expr = build_named_function(H_xp, "H_xp", y, x, p)
    H_x_expr = build_named_function(H_x, "H_x", y, x, p)
    Ψ_expr = build_named_function(Ψ, "Ψ", y, x, p)
    Ψ_yp_expr = (max_order < 2) ? nothing : build_named_function(Ψ_yp, "Ψ_yp", y, x, p)
    Ψ_y_expr = (max_order < 2) ? nothing : build_named_function(Ψ_y, "Ψ_y", y, x, p)
    Ψ_xp_expr = (max_order < 2) ? nothing : build_named_function(Ψ_xp, "Ψ_xp", y, x, p)
    Ψ_x_expr = (max_order < 2) ? nothing : build_named_function(Ψ_x, "Ψ_x", y, x, p)
    H̄_expr = build_named_function(H̄, "H̄", [y; x], p)
    H̄_w_expr = build_named_function(H̄_w, "H̄_w", [y; x], p)
    ȳ_iv_expr = build_named_function(ȳ_iv, "ȳ_iv", p)
    x̄_iv_expr = build_named_function(x̄_iv, "x̄_iv", p)
    ȳ_expr = build_named_function(ȳ, "ȳ", p)
    x̄_expr = build_named_function(x̄, "x̄", p)

    # Derivatives dispatch by symbol, and contain dummy "nothing" argument to replace
    build_function_to_dict(f, name, args...) =  Dict(key => build_named_function(f[key], name, args...;symbol_dispatch = key) for key in keys(f))

    Γ_p_expr = build_function_to_dict(Γ_p, "Γ_p", p)
    Ω_p_expr = build_function_to_dict(Ω_p, "Ω_p", p)
    H_yp_p_expr = build_function_to_dict(H_yp_p, "H_yp_p", y, x, p)
    H_y_p_expr = build_function_to_dict(H_y_p, "H_y_p", y, x, p)
    H_xp_p_expr = build_function_to_dict(H_xp_p, "H_xp_p", y, x, p)
    H_x_p_expr = build_function_to_dict(H_x_p, "H_x_p", y, x, p)
    H_p_expr = build_function_to_dict(H_p, "H_p", y, x, p)
    ȳ_p_expr = build_function_to_dict(ȳ_p, "ȳ_p", p)
    x̄_p_expr = build_function_to_dict(x̄_p, "x̄_p", p)
    Ψ_p_expr = (max_order < 2) ? nothing : build_function_to_dict(Ψ_p, "Ψ_p", y, x, p)

    verbose && printstyled("Done Building Model\n", color = :cyan)

    # Separate filenames for different orders and function types.  For example
    mkpath(model_cache_location)
    mkpath(joinpath(model_cache_location,model_name))
    # module_cache_path has the core module stuff and includes the others
    zero_order_oop_path = joinpath(model_cache_location, model_name, "zero_order_oop.jl")
    zero_order_ip_path = joinpath(model_cache_location, model_name, "zero_order_ip.jl")
    first_order_oop_path = joinpath(model_cache_location, model_name, "first_order_oop.jl")
    first_order_ip_path = joinpath(model_cache_location, model_name, "first_order_ip.jl")
    second_order_oop_path = joinpath(model_cache_location, model_name, "second_order_oop.jl")
    second_order_ip_path = joinpath(model_cache_location, model_name, "second_order_ip.jl")


    # Basic definitions are independent of the order
    open(module_cache_path, "w") do io
        write(io, "module $(model_name)\n")
        write(
            io,
            "using LinearAlgebra, SymbolicUtils, LaTeXStrings\n",
        )  # SymbolicUtils used in the generated functions
        # TODO: Use Val{...} instead?  Useful for any dispatching?
        write(io, "const max_order = $max_order\n")
        write(io, "const n_y = $n_y\n")
        write(io, "const n_x = $n_x\n")
        write(io, "const n_p = $n_p\n")
        write(io, "const n_ϵ = $n_ϵ\n")
        write(io, "const n_z = $n_z\n")
        if n_ϵ == 1
            write(io, "const η = reshape($η, $n_x, $n_ϵ)\n")
        else
            write(io, "const η = $η\n")
        end
        write(io, "const Q = $Q\n")
        write(io, "const has_Ω = $(isnothing(Ω))\n")
        write(io, "# Display definitions\n")
        write(io, "const x_symbols = $(x_subs.symbol)\n")
        write(io, "const y_symbols = $(y_subs.symbol)\n")
        write(io, "const u_symbols = $([y_subs.symbol; x_subs.symbol])\n")
        write(io, "const p_symbols = $(Symbol.(p))\n")
        write(io, "const H_latex = L\"$H_latex\"\n")
        write(io, "const steady_states_latex = L\"$steady_states_latex\"\n")
        write(io, "const steady_states_iv_latex = L\"$steady_states_iv_latex\"\n")

        write(io, "# Function definitions\n")
        save_oop && write(io, "include(\"$model_name/zero_order_oop.jl\")\n")
        save_ip && write(io, "include(\"$model_name/zero_order_ip.jl\")\n")
        save_oop && write(io, "include(\"$model_name/first_order_oop.jl\")\n")
        save_ip && write(io, "include(\"$model_name/first_order_ip.jl\")\n")
        save_oop && max_order > 1 && write(io, "include(\"$model_name/second_order_oop.jl\")\n")
        save_ip && max_order > 1 && write(io, "include(\"$model_name/second_order_ip.jl\")\n")
        return write(io, "end\n") # end module
    end

    # Zero order includes steady state calculations and derivatives
    save_ip && open(zero_order_ip_path, "w") do io
        write(io, string(Γ_expr[2]) * "\n\n")
        if isnothing(Ω)
            write(io, "const Ω! = nothing\n\n")
        else            
            write(io, string(Ω_expr[2]) * "\n\n")
        end
        write(io, string(H̄_expr[2]) * "\n\n")
        write(io, string(H̄_w_expr[2]) * "\n\n")
        write(io, string(ȳ_iv_expr[2]) * "\n\n")
        write(io, string(x̄_iv_expr[2]) * "\n\n")
        write(io, string(ȳ_expr[2]) * "\n\n")
        write(io, string(x̄_expr[2]) * "\n\n")
        write(io, "const steady_state! = nothing\n\n")
        foreach(fun -> write(io, string(fun[2]) * "\n\n"), values(Γ_p_expr))
        if isnothing(Ω)
            write(io, "const Ω_p = nothing\n\n")
        else        
            foreach(fun -> write(io, string(fun[2]) * "\n\n"), values(Ω_p_expr))
        end
        foreach(fun -> write(io, string(fun[2]) * "\n\n"), values(ȳ_p_expr))
        foreach(fun -> write(io, string(fun[2]) * "\n\n"), values(x̄_p_expr))
    end
    save_oop && open(zero_order_oop_path, "w") do io
        write(io, string(Γ_expr[1]) * "\n\n")
        if isnothing(Ω)
            write(io, "const Ω = nothing\n\n")
        else        
            write(io, string(Ω_expr[1]) * "\n\n")
        end
        write(io, string(H̄_expr[1]) * "\n\n")
        write(io, string(H̄_w_expr[1]) * "\n\n")
        write(io, string(ȳ_iv_expr[1]) * "\n\n")
        write(io, string(x̄_iv_expr[1]) * "\n\n")
        write(io, string(ȳ_expr[1]) * "\n\n")
        write(io, string(x̄_expr[1]) * "\n\n")
        write(io, "const steady_state = nothing\n\n")
        foreach(fun -> write(io, string(fun[1]) * "\n\n"), values(Γ_p_expr))
        if isnothing(Ω)
            write(io, "const Ω_p = nothing\n\n")
        else        
            foreach(fun -> write(io, string(fun[1]) * "\n\n"), values(Ω_p_expr))
        end
        foreach(fun -> write(io, string(fun[1]) * "\n\n"), values(ȳ_p_expr))
        foreach(fun -> write(io, string(fun[1]) * "\n\n"), values(x̄_p_expr))
    end

    # First order perturbations + d/dp
    save_ip && open(first_order_ip_path, "w") do io
        write(io, string(H_expr[2]) * "\n\n")
        write(io, string(H_yp_expr[2]) * "\n\n")
        write(io, string(H_y_expr[2]) * "\n\n")
        write(io, string(H_xp_expr[2]) * "\n\n")
        write(io, string(H_x_expr[2]) * "\n\n")
        write(io, string(Ψ_expr[2]) * "\n\n")
        foreach(fun -> write(io, string(fun[2]) * "\n\n"), values(H_yp_p_expr))        
        foreach(fun -> write(io, string(fun[2]) * "\n\n"), values(H_y_p_expr))        
        foreach(fun -> write(io, string(fun[2]) * "\n\n"), values(H_xp_p_expr))        
        foreach(fun -> write(io, string(fun[2]) * "\n\n"), values(H_x_p_expr))        
        foreach(fun -> write(io, string(fun[2]) * "\n\n"), values(H_p_expr))
    end
    save_oop && open(first_order_oop_path, "w") do io
        write(io, string(H_expr[1]) * "\n\n")
        write(io, string(H_yp_expr[1]) * "\n\n")
        write(io, string(H_y_expr[1]) * "\n\n")
        write(io, string(H_xp_expr[1]) * "\n\n")
        write(io, string(H_x_expr[1]) * "\n\n")
        write(io, string(Ψ_expr[1]) * "\n\n")
        foreach(fun -> write(io, string(fun[1]) * "\n\n"), values(H_yp_p_expr))        
        foreach(fun -> write(io, string(fun[1]) * "\n\n"), values(H_y_p_expr))        
        foreach(fun -> write(io, string(fun[1]) * "\n\n"), values(H_xp_p_expr))        
        foreach(fun -> write(io, string(fun[1]) * "\n\n"), values(H_x_p_expr))        
        foreach(fun -> write(io, string(fun[1]) * "\n\n"), values(H_p_expr))
    end

    # Second order perturbations + d/dp
    max_order > 1 && save_ip && open(second_order_ip_path, "w") do io
        write(io, string(Ψ_yp_expr[2]) * "\n\n")
        write(io, string(Ψ_y_expr[2]) * "\n\n")
        write(io, string(Ψ_xp_expr[2]) * "\n\n")
        write(io, string(Ψ_x_expr[2]) * "\n\n")
        foreach(fun -> write(io, string(fun[2]) * "\n\n"), values(Ψ_p_expr))        
    end

    max_order > 1 && save_oop && open(second_order_oop_path, "w") do io
        write(io, string(Ψ_yp_expr[1]) * "\n\n")
        write(io, string(Ψ_y_expr[1]) * "\n\n")
        write(io, string(Ψ_xp_expr[1]) * "\n\n")
        write(io, string(Ψ_x_expr[1]) * "\n\n")
        foreach(fun -> write(io, string(fun[1]) * "\n\n"), values(Ψ_p_expr))        
    end
    verbose && printstyled("Saved $model_name to $module_cache_path\n", color = :cyan)

    # if the module was included, gets the module name, otherwise returns nothing
    return module_cache_path
end
