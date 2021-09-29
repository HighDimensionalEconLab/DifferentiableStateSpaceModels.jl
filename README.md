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
p_d = (α=0.5, β=0.95)  # differentiated parameters
p_f = (ρ=0.2, δ=0.02, σ=0.01, Ω_1=0.01)
sol = generate_perturbation(m, p_d, p_f) # Default is first-order

# query the solution
sol.y ≈ [5.936252888048733, 6.884057971014498]
sol.x ≈ [47.39025414828825, 0.0]
sol.retcode == :Success
```

## Derivatives of the Perturbation Solvers
TBD finish and test.

The perturbation solver fills a cache for values which can be used for derivatives.

For example,
```julia
using Zygote
function f(params;m, p_f)
    p_d = (α=params[1], β=params[2])  # differentiated parameters
    sol = generate_perturbation(m, p_d, p_f) # Default is first-order
    return # TODO, Stuff with sol
end

# To call it
m = PerturbationModel(Main.my_model)
p_f = (ρ=0.2, δ=0.02, σ=0.01, Ω_1=0.01)
param_val = [0.5, 0.95] # as a vector, but not required
f(param_val;m, p_f) # Function works on its own, calculating perturbation

# But you can also gets its gradient with Zygote/etc.
gradient(params -> f(params; m, p_f), param_val)
```

## Example Usage for HMC Estimation

TBD
