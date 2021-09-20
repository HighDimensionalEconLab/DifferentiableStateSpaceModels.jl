function make_substitutions(t, f_var)
    #t = first(@variables t)  # harcoded index for now
    sym_name = f_var.f.name
    sym_name_p = Symbol(string(sym_name) * "_p")
    sym_name_ss = Symbol(string(sym_name) * "_ss")
    names = @variables $sym_name $sym_name_p $sym_name_ss
    return (symbol = sym_name,
            var = names[1],
            var_p = names[2],
            var_ss = names[3],
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
substitute_and_simplify(f::Nothing, subs; simplify=true) = nothing


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