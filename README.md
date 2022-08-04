# DifferentiableStateSpaceModels

[![Build Status](https://github.com/HighDimensionalEconLab/DifferentiableStateSpaceModels.jl/workflows/CI/badge.svg)](https://github.com/HighDimensionalEconLab/DifferentiableStateSpaceModels.jl/actions)
[![Coverage](https://codecov.io/gh/HighDimensionalEconLab/DifferentiableStateSpaceModels.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/HighDimensionalEconLab/DifferentiableStateSpaceModels.jl)

## Development and Benchmarking
See [development.md](development.md) for contributing code and running benchmarks

## Defining Models
Install this package with `] add DifferentiableStateSpaceModels` then the full code to create the RBC model is

```julia
using DifferentiableStateSpaceModels, Symbolics, Test

∞ = Inf
@variables α, β, ρ, δ, σ, Ω_1
@variables t::Integer, k(..), z(..), c(..), q(..)

x = [k, z]
y = [c, q]
p = [α, β, ρ, δ, σ, Ω_1]

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

Ω = [Ω_1, Ω_1]

# Generates the files and includes if required.  If the model is already created, then just loads
# Saves as ".function_cache/my_model.jl"
m = @make_and_include_perturbation_model("my_model", H, (; t, y, x, p, steady_states, Γ, Ω, η, Q)) 
```

After generation of the model, they can be included as any other julia files in your code (e.g. `include(joinpath(pkgdir(DifferentiableStateSpaceModels), ".function_cache","my_model.jl"))` or moved somewhere more convenient.

Inclusion through the `@make_and_include_perturbation_model` creates the model automatically; after direct inclusion through a julia file, you can create a model with `m = PerturbationModel(Main.my_model)`.

## Solving Perturbations
Assuming the above model was created and loaded in one way or another as `m`

```julia
p_d = (α=0.5, β=0.95)  # Differentiated parameters
p_f = (ρ=0.2, δ=0.02, σ=0.01, Ω_1=0.01)
sol = generate_perturbation(m, p_d, p_f) # Default is first-order

# Query the solution
@assert sol.y ≈ [5.936252888048733, 6.884057971014498]
@assert sol.x ≈ [47.39025414828825, 0.0]
sol.retcode == :Success
```

## Derivatives of the Perturbation Solvers

The perturbation solver fills a cache for values used for calculating derivatives.

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
display(f(param_val; m, p_f)) # Function works on its own, calculating perturbation
# Query the solution
@assert f(param_val; m, p_f) ≈ 7.366206154679124

# But you can also get its gradient with Zygote/etc.
display(gradient(params -> f(params; m, p_f), param_val))
# Result check
gradient(params -> f(params; m, p_f), param_val)[1] ≈ [61.41968376547458, 106.44095661062319]
```

## Using Solution Cache

You can also pass your own cache into the solver, which can then be reused later to avoid costly memory re-allocations from repeated runs.
```julia
using Zygote
function f2(p_d; m, p_f, cache)
    sol = generate_perturbation(m, p_d, p_f; cache) # Default is first-order.
    return sum(sol.A) # An ad-hoc example: reducing the law-of-motion matrix into one number
end

# To call it
m = PerturbationModel(Main.my_model)
p_d = (α=0.5, β=0.95)  # Differentiated parameters
p_f = (ρ=0.2, δ=0.02, σ=0.01, Ω_1=0.01)
cache = SolverCache(m, Val(1), p_d)
display(f2(p_d; m, p_f, cache)) # Function works on its own, calculating perturbation
# Query the solution
@assert f2(p_d; m, p_f, cache) ≈ 7.366206154679124

# But you can also get its gradient with Zygote/etc.
gradient(params -> f2(params; m, p_f, cache), p_d)[1]
```

## Simulating Data with DifferenceEquations

Given the model solutions, we can then use [DifferenceEquations.jl](https://github.com/SciML/DifferenceEquations.jl) to simulate an ensemble of trajectories.
```julia
using DifferenceEquations, Distributions, Plots, DiffEqBase
p_f = (ρ = 0.2, δ = 0.02, σ = 0.01, Ω_1 = 0.01) # Fixed parameters
p_d = (α = 0.5, β = 0.95) # Pseudo-true values
sol = generate_perturbation(m, p_d, p_f) # Solution to the first-order RBC
sol_2 = generate_perturbation(m, p_d, p_f, Val(2)) # Solution to the second-order RBC

# Simulate multiple trajectories with T observations
T = 20
trajectories = 40

# draw from ergodic for each simulation
x_iv = MvNormal(sol.x_ergodic_var) # draw initial conditions from the ergodic distribution
problem = LinearStateSpaceProblem(sol, x_iv, (0, T))

# Solve multiple trajectories and plot an ensemble
ensemble_results = solve(EnsembleProblem(problem), DirectIteration(), EnsembleThreads();
                 trajectories)
plot(EnsembleSummary(ensemble_results))
```

## Example Usage for State-Space Model Estimation

With [Turing](https://turing.ml/stable/), a Julia package for Bayesian inference with
probabilistic programming, we can estimate State-Space models using Hamiltonian Monte-Carlo methods.

You can check out this Jupyter notebook for [Full Example Usage for RBC Estimation](notebooks/rbc_example.ipynb).
