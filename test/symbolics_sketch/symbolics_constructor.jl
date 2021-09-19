using Symbolics, SymbolicUtils, MacroTools, StructArrays, Test, LinearAlgebra

# TO UTILITIES
function make_substitutions(f_var)
    t = first(@variables t)  # harcoded index for now
    sym_name = f_var.f.name
    sym_name_p = Symbol(string(sym_name) * "_p")
    sym_name_ss = Symbol(string(sym_name) * "_ss")
    names = @variables $sym_name $sym_name_p $sym_name_ss
    return (symbol = sym_name,
            var = names[1],
            var_p = names[2],
            markov_t = f_var(t) => names[1],
            markov_tp1 = f_var(t+1) => names[2],
            markov_inf = f_var(Inf) => names[3],
            tp1_to_var = names[2] => names[1],
            inf_to_var = names[3] => names[1])
end

# Extracts from named tuple, dictionary, etc. tp create a new vector in the order of "symbols"
arrange_vector_from_symbols(x, symbols) = [x[sym] for sym in symbols]

substitute_and_simplify(f::Num, subs; simplify=true) = simplify ? Symbolics.simplify(substitute(f, subs)) : Symbolics.substitute(f, subs)
substitute_and_simplify(f::AbstractArray, subs; simplify = true) = substitute_and_simplify.(f, Ref(subs); simplify)

#Variations of differentiate depending which create matrices, vectors of matrices, etc.
#Recursion isn't quite right because differentiating a vector gives a matrix rather than a vector of vectors.
#Later could try Array{Num, 3} instead if algorithms can be organized appropriately - at which point recursion to tensors makes more sense.
nested_differentiate(f::Vector{Num}, x::Vector{Num}; simplify = true) = [expand_derivatives(Differential(var)(f_val), simplify) for f_val in f, var in x]

nested_differentiate(f::Matrix{Num}, x::Vector{Num}; simplify = true) = [expand_derivatives.(Differential(var).(f), simplify) for var in x]

nested_differentiate(f::Vector{<:Real}, x::Num; simplify = true) = [expand_derivatives(Differential(x)(f_val), simplify) for f_val in f]

nested_differentiate(f::Matrix{Num}, x::Num; simplify = true) = expand_derivatives.(Differential(x).(f), simplify)

nested_differentiate(f::Vector{Matrix{Num}}, x::Num; simplify = true) = [expand_derivatives.(Differential(x).(f_val), simplify) for f_val in f]  #e.g. d psi for a variable

nested_differentiate(f::Vector{Matrix{Num}}, x::Vector{Num}; simplify = true) = [nested_differentiate(f, x_val;simplify) for x_val in x]

nested_differentiate(::Nothing, x) = nothing
nested_differentiate(f, ::Nothing) = nothing

#e.g. d psi for a vector

# Names the expression, and optionally repalces the first argument (after the out) with a dispatch by symbol
function name_symbolics_function(expr, name;inplace = false, symbol_dispatch=nothing, striplines = true)
    # add name for dispatching, and an argument for the derivative
    expr_dict = splitdef(expr)
    expr_dict[:name] = name # Must be a symbol
    
    #if replacing first parameter to dispatching on a symbol
    if !isnothing(symbol_dispatch)
        dispatch_position = inplace ? 2 : 1
        expr_dict[:args][dispatch_position] = :(::Val{$(QuoteNode(symbol_dispatch))})
    end
    named_expr = combinedef(expr_dict)
    return striplines ? MacroTools.striplines(named_expr) : named_expr
end

######### For Test Setup

# Setup for test
const ∞ = Inf
@variables α, β, ρ, δ, σ, Ω_1, Ω_2
@variables t::Integer, k(..), z(..), c(..), q(..)

x = [k, z]
y = [c, q]
p = [α, β, ρ, δ, σ]

H = [
    1 / c(t) - (β / c(t+1)) * (α * exp(z(t+1)) * k(t+1)^(α - 1) + (1 - δ)),
    c(t) + k(t+1) - (1 - δ) * k(t) - q(t),
    q(t) - exp(z(t)) * k(t)^α,
    z(t+1) - ρ * z(t),
]

steady_states = [k(∞) ~ (((1 / β) - 1 + δ) / α)^(1 / (α - 1)),
z(∞) ~ 0,
c(∞) ~ (((1 / β) - 1 + δ) / α)^(α / (α - 1)) -
        δ * (((1 / β) - 1 + δ) / α)^(1 / (α - 1)),
q(∞) ~ (((1 / β) - 1 + δ) / α)^(α / (α - 1)),
]

steady_states_iv = [k(∞) ~ (((1 / β) - 1 + δ) / α)^(1 / (α - 1)),
z(∞) ~ 0,c(∞) ~ (((1 / β) - 1 + δ) / α)^(α / (α - 1)) -
        δ * (((1 / β) - 1 + δ) / α)^(1 / (α - 1)),
q(∞) ~ (((1 / β) - 1 + δ) / α)^(α / (α - 1)),
]

n_ϵ = 1
n_z = 2
n_x = length(x)
n_y = length(y)
n_p = length(p)
Γ = reshape([σ], n_ϵ, n_ϵ)
η = reshape([0; -1], n_x, n_ϵ) # η is n_x * n_ϵ matrix

Q = zeros(n_z, n_x + n_y)
Q[1, 1] = 1.0
Q[2, 3] = 1.0

Ω = [Ω_1, Ω_2]

model_name = "rbc_obervables"
overwrite_model_cache = true
verbose = true
model_cache_location = "./test/symbolics_sketch"

###### INSIDE FUNCTION

    model_cache_path = joinpath(model_cache_location, model_name * ".jl")

    # only load cache if the module isn't already loaded in memory
    if (isdefined(Main, Symbol(model_name)) && !overwrite_model_cache)
        verbose && printstyled("Using existing module $model_name\n", color = :cyan)
        return model_cache_path = model_cache_path
    end

    # if path already exists
    if (ispath(model_cache_path) && !overwrite_model_cache)
        # path exists and not overwriting
        verbose &&
            printstyled("Model already generated at $model_cache_path\n", color = :cyan)
        return model_cache_path    
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
    y_subs = StructArray(make_substitutions.(y))
    x_subs = StructArray(make_substitutions.(x))
    y = y_subs.var
    x = x_subs.var
    y_p = y_subs.var_p
    x_p = x_subs.var_p
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
    Ψ = [Symbolics.hessian(f, [y_p; y; x_p; x]; simplify=true) for f in H]
    Ψ_yp = nested_differentiate(Ψ, y_p)
    Ψ_y = nested_differentiate(Ψ, y)
    Ψ_xp = nested_differentiate(Ψ, x_p)
    Ψ_x = nested_differentiate(Ψ, x)

    
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
    Ψ_p = differentiate_to_dict(Ψ, p)

    # apply substitutions and simplify
    H_yp = substitute_and_simplify(H_yp, all_to_var)
    H_xp = substitute_and_simplify(H_xp, all_to_var)
    H_x = substitute_and_simplify(H_x, all_to_var)
    H_y = substitute_and_simplify(H_y, all_to_var)
    Ψ = substitute_and_simplify(Ψ, all_to_var)
    Ψ_yp = substitute_and_simplify(Ψ_yp, all_to_var)
    Ψ_y = substitute_and_simplify(Ψ_y, all_to_var)
    Ψ_xp = substitute_and_simplify(Ψ_xp, all_to_var)
    Ψ_x = substitute_and_simplify(Ψ_x, all_to_var)

    # Generate all functions and to save to files.
    Γ_expr = build_dssm_function(Γ, p, p_f)
    Γ_p_expr = build_dssm_function(Γ_p, p, p_f)
    Ω_expr = build_dssm_function(Ω, p, p_f)
    Ω_p_expr = build_dssm_function(Ω_p, p, p_f)
    H_expr = build_dssm_function(H, y_p, y, y_ss, x_p, x, x_ss, p, p_f)
    H_yp_expr = build_dssm_function(H_yp_sub, y, x, p, p_f)
    H_y_expr = build_dssm_function(H_y_sub, y, x, p, p_f)
    H_xp_expr = build_dssm_function(H_xp_sub, y, x, p, p_f)
    H_x_expr = build_dssm_function(H_x_sub, y, x, p, p_f)
    H_yp_p_expr = build_dssm_function(H_yp_p_sub, y, x, p, p_f)
    H_y_p_expr = build_dssm_function(H_y_p_sub, y, x, p, p_f)
    H_xp_p_expr = build_dssm_function(H_xp_p_sub, y, x, p, p_f)
    H_x_p_expr = build_dssm_function(H_x_p_sub, y, x, p, p_f)
    H_p_expr = build_dssm_function(H_p_sub, y, x, p, p_f)
    Ψ_expr = build_dssm_function(Ψ_sub, y, x, p, p_f)
    Ψ_p_expr = build_dssm_function(Ψ_p_sub, y, x, p, p_f)
    Ψ_yp_expr = build_dssm_function(Ψ_yp_sub, y, x, p, p_f)
    Ψ_y_expr = build_dssm_function(Ψ_y_sub, y, x, p, p_f)
    Ψ_xp_expr = build_dssm_function(Ψ_xp_sub, y, x, p, p_f)
    Ψ_x_expr = build_dssm_function(Ψ_x_sub, y, x, p, p_f)
    H̄_expr = build_dssm_function(H̄_sub, [y; x], p, p_f)
    H̄_w_expr = build_dssm_function(H̄_w, [y; x], p, p_f)
    ȳ_iv_expr = build_dssm_function(ȳ_iv, p, p_f)
    x̄_iv_expr = build_dssm_function(x̄_iv, p, p_f)
    ȳ_expr = build_dssm_function(ȳ, p, p_f)
    x̄_expr = build_dssm_function(x̄, p, p_f)
    ȳ_p_expr = build_dssm_function(ȳ_p_sub, p, p_f)
    x̄_p_expr = build_dssm_function(x̄_p_sub, p, p_f)

    verbose && printstyled("Done Building Model\n", color = :cyan)

    mkpath(model_cache_location)
    open(model_cache_path, "w") do io
        write(io, "module $(model_name)\n")
        write(
            io,
            "using LinearAlgebra, DifferentiableStateSpaceModels, ModelingToolkit\n",
        )
        write(io, "const n_y = $n_y\n")
        write(io, "const n_x = $n_x\n")
        write(io, "const n = $n\n")
        write(io, "const n_p = $n_p\n")
        write(io, "const n_ϵ = $n_ϵ\n")
        write(io, "const n_z = $n_z\n")
        if n_ϵ == 1
            write(io, "const η = reshape($η, $n_x, $n_ϵ)\n")
        else
            write(io, "const η = $η\n")
        end
        write(io, "const Q = $Q\n")
        write(io, "const Γ! = $(Γ_expr)\n")
        write(io, "const Γ_p! = $(Γ_p_expr)\n")
        write(io, "const Ω! = $(Ω_expr)\n")
        write(io, "const Ω_p! = $(Ω_p_expr)\n")
        write(io, "const H! = $(H_expr)\n")
        write(io, "const H_yp! = $(H_yp_expr)\n")
        write(io, "const H_y! = $(H_y_expr)\n")
        write(io, "const H_xp! = $(H_xp_expr)\n")
        write(io, "const H_x! = $(H_x_expr)\n")
        write(io, "const H_yp_p! = $(H_yp_p_expr)\n")
        write(io, "const H_y_p! = $(H_y_p_expr)\n")
        write(io, "const H_xp_p! = $(H_xp_p_expr)\n")
        write(io, "const H_x_p! = $(H_x_p_expr)\n")
        write(io, "const H_p! = $(H_p_expr)\n")
        write(io, "const Ψ! = $(Ψ_expr)\n")
        write(io, "const Ψ_p! = $(Ψ_p_expr)\n")
        write(io, "const Ψ_yp! = $(Ψ_yp_expr)\n")
        write(io, "const Ψ_y! = $(Ψ_y_expr)\n")
        write(io, "const Ψ_xp! = $(Ψ_xp_expr)\n")
        write(io, "const Ψ_x! = $(Ψ_x_expr)\n")
        write(io, "const H̄! = $(H̄_expr)\n")
        write(io, "const H̄_w! = $(H̄_w_expr)\n")
        write(io, "const ȳ_iv! = $(ȳ_iv_expr)\n")
        write(io, "const x̄_iv! = $(x̄_iv_expr)\n")
        write(io, "const ȳ! = $(ȳ_expr)\n")
        write(io, "const x̄! = $(x̄_expr)\n")
        write(io, "const ȳ_p! = $(ȳ_p_expr)\n")
        write(io, "const x̄_p! = $(x̄_p_expr)\n")
        write(io, "const steady_state! = nothing\n")
        return write(io, "end\n") # end module
    end
    verbose && printstyled("Saved $model_name to $model_cache_path\n", color = :cyan)
    return model_cache_path
