using Symbolics, SymbolicUtils, MacroTools, StructArrays, Test

function make_substitutions(t, f_var)
    sym_name = f_var.f.name
    sym_name_p = Symbol(string(sym_name) * "_p")
    sym_name_ss = Symbol(string(sym_name) * "_ss")
    names = @variables $sym_name $sym_name_p $sym_name_ss
    return (symbol = sym_name,
            var = names[1],
            markov_t = f_var(t) => names[1],
            markov_tp1 = f_var(t+1) => names[2],
            markov_inf = f_var(Inf) => names[3],
            tp1_to_var = names[2] => names[1],
            inf_to_var = names[3] => names[1])
end


# Extracts from named tuple, dictionary, etc. tp create a new vector in the order of "symbols"
arrange_vector_from_symbols(x, symbols) = [x[sym] for sym in symbols]

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

# Tests
const ∞ = Inf
@variables α, β, ρ, δ, σ
@variables t::Integer, k(..), z(..), c(..), q(..)
x = [k, z]
y = [c, q]
p = [α, β, ρ, δ, σ]
subs_x = StructArray(make_substitutions.(t, x))
subs_y = StructArray(make_substitutions.(t, y))
subs = vcat(subs_x, subs_y)
p_symbols = Symbol.(p)

p_val = (α = 0.1, β = 0.5, ρ = 0.1, δ = 1.9, σ= 1.9)
p_val_2 = (ρ = 0.1, α = 0.1, σ= 1.9, β = 0.5, δ = 1.9)
p_vec = arrange_vector_from_symbols(p_val, p_symbols)
p_vec_2 = arrange_vector_from_symbols(p_val_2, p_symbols)
@test p_vec ≈ p_vec_2


# Test for the created functions
@variables x, y
func = [x, y]
func_2 = [x^2, y]
u = [x,y]
u_val = [0.1, 0.6] # [x,y]
# No symbolic dispatching
ex, ex_ip = build_function(func, u;linenumbers = false)  # the nothing placeholder is for the Val{symbol} dispatch
named_ex = name_symbolics_function(ex, :my_func)
named_ex_ip = name_symbolics_function(ex_ip, :my_func!; inplace = true)
eval(named_ex)
eval(named_ex_ip)

@test my_func(u_val) ≈ [0.1, 0.6]
out = zeros(2)
my_func!(out, u_val)
@test out ≈ [0.1, 0.6]

# Symbolic dispatching
ex2, ex2_ip = build_function(func_2, nothing, u;linenumbers = false)  # the nothing placeholder is for the Val{symbol} dispatch
dispatch_by = :a
named_ex2 = name_symbolics_function(ex2, :func_symb; symbol_dispatch = dispatch_by)
named_ex2_ip = name_symbolics_function(ex2_ip, :func_symb!; inplace = true, symbol_dispatch = dispatch_by)

# Use them
eval(named_ex2)
eval(named_ex2_ip)

@test func_symb(Val(dispatch_by), u_val) ≈ [0.01, 0.6]
out = zeros(2)
func_symb!(out, Val(dispatch_by), u_val)
@test out ≈ [0.01, 0.6]


 

##########
const ∞ = Inf
@variables α, β, ρ, δ, σ
@variables t::Integer, k(..), z(..), c(..), q(..)
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
# @show Symbolics.get_variables(x[1])
x = [k, z]
y = [c, q]
subs_x = StructArray(make_substitutions.(t, x))
subs_y = StructArray(make_substitutions.(t, y))
subs = vcat(subs_x, subs_y)
subs_all_to_markov = vcat(subs.markov_t, subs.markov_tp1, subs.markov_inf)
subs_all_to_var = vcat(subs.tp1_to_var, subs.inf_to_var)
H_markov = substitute.(H, Ref(subs_all_to_markov))

# Sort by the symbols
steady_states_dict = Dict(Symbol(substitute(substitute(eq.lhs, subs_all_to_markov), subs_all_to_var)) => substitute(eq.rhs, subs_all_to_markov) for eq in steady_states)
steady_state_vector = arrange_vector_from_symbols(steady_states_dict, subs.symbol)