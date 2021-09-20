using DifferentiableStateSpaceModels, Symbolics

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

model_name = "rbc_observables"
overwrite_model_cache = true
verbose = true
save_ip = true
save_oop = false
max_order = 2
skipzeros = false
fillzeros = false

module_cache_path = generate_perturbation_model(H;model_name, t,y,x,p, steady_states, steady_states_iv,Γ,Ω,η,Q,overwrite_model_cache,verbose,max_order,save_ip,save_oop,skipzeros,fillzeros)

generate_perturbation_model(H;model_name, t,y,x,p, steady_states, steady_states_iv,Γ,Ω,η,Q,overwrite_model_cache = false,verbose,max_order,save_ip,save_oop,skipzeros,fillzeros)