# DifferentiableStateSpaceModels

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://HighDimensionalEconLab.github.io/DifferentiableStateSpaceModels.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://HighDimensionalEconLab.github.io/DifferentiableStateSpaceModels.jl/dev)
[![Build Status](https://github.com/HighDimensionalEconLab/DifferentiableStateSpaceModels.jl/workflows/CI/badge.svg)](https://github.com/HighDimensionalEconLab/DifferentiableStateSpaceModels.jl/actions)
[![Coverage](https://codecov.io/gh/HighDimensionalEconLab/DifferentiableStateSpaceModels.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/HighDimensionalEconLab/DifferentiableStateSpaceModels.jl)

## Development and Benchmarking
See [development.md](development.md) for contributing code and running benchmarks

## Defining Models
Install this package with `] add DifferentiableStateSpaceModels` then the full code to create the RBC model is

```julia
using DifferentiableStateSpaceModels, Symbolics, Test

∞ = Inf
@variables α, β, ρ, δ, σ, Ω_1, Ω_2
@variables t::Integer, k(..), z(..), c(..), q(..)

x = [k, z]
y = [c, q]
p = [α, β, ρ, δ, σ]

H = [1 / c(t) - (β / c(t + 1)) * (α * exp(z(t + 1)) * k(t + 1)^(α - 1) + (1 - δ)),
     c(t) + k(t + 1) - (1 - δ) * k(t) - q(t),
     q(t) - exp(z(t)) * k(t)^α,
     z(t + 1) - ρ * z(t)]

steady_states = [k(∞) ~ (((1 / β) - 1 + δ) / α)^(1 / (α - 1)),
                 z(∞) ~ 0,
                 c(∞) ~ (((1 / β) - 1 + δ) / α)^(α / (α - 1)) -
                        δ * (((1 / β) - 1 + δ) / α)^(1 / (α - 1)),
                 q(∞) ~ (((1 / β) - 1 + δ) / α)^(α / (α - 1))]

Γ = reshape([σ], 1, 1)
η = reshape([0; -1], length(x), 1) # η is n_x * n_ϵ matrix

Q = zeros(2, length(x) + length(y))
Q[1, 1] = 1.0
Q[2, 3] = 1.0

Ω = [Ω_1, Ω_2]

# Arguments for creating the model
args = (; t, y, x, p, steady_states, Γ, Ω, η, Q)

# Generates the files and includes if required.  If model already created, then just loads
m = @make_and_include_perturbation_model("my_model", H, args) # Convenience macro
```

After generation of the model, they can be included as any other julia files in your code (e.g. `include(joinpath(pkgdir(DifferentiableStateSpaceModels), ".function_cache","my_model.jl"))` or moved somewhere more convenient.

After inclusion either throught the `@make_and_include_perturbation_mode` or direction inclusion, you can create a model with `m = PerturbationModel(Main.my_model)`

## Solving Perturbations
Assuming the above model was created and loaded in one way or another

```julia
m = @make_and_include_perturbation_model("my_model", H, args) # or m = PerturbationModel(Main.my_model)
p_d = (α=0.5, β=0.95)  # Differentiated parameters
p_f = (ρ=0.2, δ=0.02, σ=0.01, Ω_1=0.01)
sol = generate_perturbation(m, p_d, p_f) # Default is first-order

# Query the solution
sol.y ≈ [5.936252888048733, 6.884057971014498]
sol.x ≈ [47.39025414828825, 0.0]
sol.retcode == :Success
```

## Derivatives of the Perturbation Solvers

The perturbation solver fills a cache for values which can be used for derivatives.

For example,
```julia
using Zygote
function f(params; m, p_f)
    p_d = (α=params[1], β=params[2])  # Differentiated parameters
    sol = generate_perturbation(m, p_d, p_f) # Default is first-order.
    return sum(sol.A) # An ad-hoc example: reducing the law-of-motion matrix into one number
end

# To call it
m = PerturbationModel(Main.my_model)
p_f = (ρ=0.2, δ=0.02, σ=0.01, Ω_1=0.01)
param_val = [0.5, 0.95] # as a vector, but not required
f(param_val; m, p_f) # Function works on its own, calculating perturbation
# Query the solution
f(param_val; m, p_f) ≈ 7.366206154679124

# But you can also gets its gradient with Zygote/etc.
gradient(params -> f(params; m, p_f), param_val)
# Result check
gradient(params -> f(params; m, p_f), param_val)[1] ≈ [61.41968376547458, 106.44095661062319]
```

## Example Usage for State-Space Model Estimation

With [Turing](https://turing.ml/stable/), a Julia package for Bayesian inference with
probabilistic programming, we can estimate State-Space models using Hamiltonian Monte-Carlo methods.

Following the RBC example above, we show the model estimation with simulated data.

```julia
# Synthetic data generation
p_d = (α=0.5, β=0.95) # Pseudo-true parameters
p_f = (ρ=0.2, δ=0.02, σ=0.01, Ω_1=0.1)
sol = generate_perturbation(m, p_d, p_f, Val(1))
sol_second = generate_perturbation(m, p_d, p_f, Val(2))

T = 20 # Length of synthetic data
ϵ = [randn(model_rbc.n_ϵ) for _ = 1:T] # Shocks
x0 = zeros(model_rbc.n_x) # Initial state
fake_z = solve(sol, x0, (0, T), DifferentiableStateSpaceModels.LTI(); noise = ϵ).z # First-order synthetic data
fake_z_second = solve(sol_second, x0, (0, T), DifferentiableStateSpaceModels.QTI(); noise = ϵ).z # Second-order synthetic data
```

Estimate the RBC model in first-order with Kalman filter derived marginal likelihood:
```julia
using Turing
using Turing: @addlogprob!
Turing.setadbackend(:zygote)

# Specify the Turing model
@model function rbc_kalman(z, m, p_f)
    # List the priors of the parameters
    α ~ Uniform(0.2, 0.8)
    β ~ Uniform(0.5, 0.99)
    p_d = (α = α, β = β) # Put all the parameters into a NamedTuple
    sol = generate_perturbation(m, p_d, p_f, Val(1))
    if !(sol.retcode == :Success)
        @addlogprob! -Inf
        return
    end
    @addlogprob! solve(sol, sol.x_ergodic, (0, length(z)); observables = z).logpdf
end

turing_model = rbc_kalman(fake_z, model_rbc, p_f)
n_samples = 1000 # Number of total samples to draw
n_adapts = 100 # Number of adapts for the No-U-Turn sampler
δ = 0.65 # Target acceptance rate of the No-U-Turn sampler
chain = sample(turing_model, NUTS(n_adapts, δ), n_samples; progress = true)
```

Estimate the RBC model in second-order with joint likelihood approach:
```julia
using Turing
using Turing: @addlogprob!
Turing.setadbackend(:zygote)

@model function rbc_second(z, m, p_f, x0 = zeros(m.n_x))
    α ~ Uniform(0.2, 0.8)
    β ~ Uniform(0.5, 0.99)
    p_d = (α = α, β = β)
    T = length(z)
    ϵ_draw ~ MvNormal(T, 1.0)
    ϵ = map(i -> ϵ_draw[((i-1)*m.n_ϵ+1):(i*m.n_ϵ)], 1:T)
    sol = generate_perturbation(m, p_d, p_f, Val(2))
    if !(sol.retcode == :Success)
        @addlogprob! -Inf
        return
    end
    @addlogprob! solve(sol, x0, (0, T); noise = ϵ, observables = z).logpdf
end

turing_model = rbc_joint(fake_z_second, model_rbc, p_f, c)
n_samples = 1000
n_adapts = 100
δ = 0.65
max_depth = 5
chain = sample(turing_model, NUTS(n_adapts, δ; max_depth), n_samples; progress = true)
```
