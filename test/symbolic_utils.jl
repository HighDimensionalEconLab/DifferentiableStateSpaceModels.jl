using Test
using Symbolics, StructArrays
using DifferentiableStateSpaceModels: order_vector_by_symbols, make_substitutions

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
p_vec = order_vector_by_symbols(p_val, p_symbols)
p_vec_2 = order_vector_by_symbols(p_val_2, p_symbols)
@test p_vec ≈ p_vec_2

# Check for skipzeros with symbolics.jl
@variables x,y
f = [0, x]
f_expr = build_function(f, [x,y];skipzeros = true, fillzeros=true, expression = Val{false})

out = Vector{Float64}(undef, 2)
u = [5.0, 3.1]
out_oop = f_expr[1](u)
f_expr[2](out, u)
@test_broken out_oop ≈ out
