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


vector_from_symbols(x, symbols) = [x[sym] for sym in symbols]

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
p_vec = vector_from_symbols(p_val, p_symbols)
p_vec_2 = vector_from_symbols(p_val_2, p_symbols)
@test p_vec ≈ p_vec_2


# 

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

q(::Val{:a}, x) = x
q(::Val{:b}, x) = x^2

function m(sym)
    return q(Val(:a), 0.1)
end
m(:a)

q(Val(:a), 3)

#####################
@variables x, y
func = [x, y]
u = [x,y]

ex1, ex2 = build_function(func, nothing, u;linenumbers = false)  # the nothing placeholder is for the Val{symbol} dispatch

# add name for dispatching, and an argument for the derivative
ex_dict = splitdef(ex1)
ex_dict_2 = splitdef(ex2)
ex_dict[:name] = :my_func
ex_dict_2[:name] = :my_func!
ex_dict[:args][1] = :(::Val{:a})
ex_dict_2[:args][2] = :(::Val{:a})  # second position for these
ex1_mod = combinedef(ex_dict)
ex2_mod = combinedef(ex_dict_2)
