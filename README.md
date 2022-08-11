# DifferentiableStateSpaceModels

[![Build Status](https://github.com/HighDimensionalEconLab/DifferentiableStateSpaceModels.jl/workflows/CI/badge.svg)](https://github.com/HighDimensionalEconLab/DifferentiableStateSpaceModels.jl/actions)
[![Coverage](https://codecov.io/gh/HighDimensionalEconLab/DifferentiableStateSpaceModels.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/HighDimensionalEconLab/DifferentiableStateSpaceModels.jl)

**Warning:**  This package is a proof of concept.  While the code will remain working with proper use of a [Julia manifest](https://pkgdocs.julialang.org/v1/), you should not rely on this for research unless you are prepared to modify the source and help maintain it.

For a more complete example with the code below and Bayesian estimation, see the [rbc estimation notebook](notebooks/estimate_rbc.ipynb).

## Development and Benchmarking
See [development.md](development.md) for contributing code and running benchmarks

# Model Class
The model follows [Schmitt-Grohe and Uribe (2004)](http://www.columbia.edu/~mu2166/2nd_order/2nd_order.pdf) timing convention.  The system takes a nonlinear expectational difference equation including all first-order conditions for decisions and the system evolution equations,

$$
\mathbb{E}_{t}\mathcal{H}\left(y',y,x',x;p\right)=0
$$

where $y$ are the control variables, $x$ are the states, and $p$ is a vector of deep parameters of interest.  Expectations are taken over forward-looking variables and an underlying random process $\epsilon'$.


In addition, we consider an observation equation - which might be noisy, for
$$
z = Q \cdot \begin{bmatrix}y &x\end{bmatrix}^{\top} + \nu
$$
where $\nu$ may or may not be normally distributed but $\mathbb{E}(\nu) = 0$ and $\mathbb{V}(\nu) = \Omega(p) \Omega(p)^{\top}$.

Assume that there is a non-stochastic steady state of this problem as $y_{ss}, x_{ss}$.

## Perturbation Solution
Define the deviation from the non-stochastic steady state as $\hat{x} \equiv x - x_{ss}, \hat{y} \equiv y - y_{ss},$ and $\hat{z} \equiv z - z_{ss}$.

The solution finds the first or second order perturbation around that non-stochastic steady state, and yields

$$
x' = h(x; p) + \eta \ \Gamma(p)\ \epsilon'
$$

where $\eta$ describes how shocks affect the law of motion and $\mathbb{E}(\epsilon') = 0$.  Frequently this would be organized such that $\mathbb{V}(\epsilon)= I$, but that is not required.  In addition, it could instead be interpreted as for $x' = h(x; p) + \eta \ \epsilon'$ with $\mathbb{V}(\epsilon') = \Gamma(p) \Gamma(p)^{\top}$.

and with the policy equation,

$$
\hat{y} = g(\hat{x}; p)
$$


$$
y = g(x; p)
$$

and finally, substitution in for the observation equation

$$
z= Q \begin{bmatrix} g(x;p) \\ x\end{bmatrix} + \nu
$$

## First Order Solutions
Perturbation approximates the above model equations, where $h$ and $g$ are not available explicitly, by a Taylor expansion around the steady state. For example, in the case of the 1st order model the solution finds

$$
\hat{x}' = A(p)\ \hat{x} + B(p) \epsilon'
$$

and

$$
\hat{y} = g_x(p) \ \hat{x}
$$

and 

$$
\hat{z} = C(p)\ \hat{x} + \nu
$$
where $C(p)\equiv Q \begin{bmatrix} g_x(p) \\ I\end{bmatrix}$, $B(p) \equiv \eta \Gamma(p)$, and $\mathbb{V}(v) = D(\nu) D(p)^{\top}$.  Normality of $\nu$ or $\epsilon'$ is not required in general.

This is a linear state-space model and if the priors and shocks are Gaussian, a marginal likelihood can be evaluated with classic methods such as a Kalman Filter.  The output of the perturbation can be used manually, or in conjunction with [DifferenceEquations.jl](https://github.com/SciML/DifferenceEquations.jl).

Second-order solutions are defined similarly.  See the estimation notebook for more details.
## Gradients
All of the above use standard solution methods.  The primary contribution of this package is that all of these model elements are **differentiable**. Hence, these gradeients can be composed for use in applications such as optimization, gradient-based estimation methods, and with  [DifferenceEquations.jl](https://github.com/SciML/DifferenceEquations.jl)  which provides differentiable simulations and likelihoods for state-space models.  That is, if we think of a perturbation solver mapping $p$ to solutions (e.g. in first order $\mathbf{P}(p) \to (A, B, C, D)$, then we can find the gradients $\partial_p \mathbf{P}(p), \partial_p A(p)$ etc.  Or, when these gradients are available for use with reverse-mode auto-differentiation, it can take "wobbles" in $A, B, C, D$ and go back to the "wiggles" of the underlying $p$ through $\mathbf{P}$.  See [ChainRules.jl](https://juliadiff.org/ChainRulesCore.jl/stable/maths/propagators.html) for more details on AD.  
# Examples

## Model Primitives
Models are defined using a Dynare-style DSL using [Symbolics.jl](https://github.com/JuliaSymbolics/Symbolics.jl).  The list of primitives are:
1. The list of variables for the controls $y$, state $x$, and deep parameters $p$.
2. The set of equations $H$ as a function of $p, y(t), y(t+1), x(t),$ and $x(t+1)$.  No $t-1$ timing is allowed.
3. The loading of shocks $\eta$ as a fixed matrix of constants
4. The shock covariance Cholesky factor $\Gamma$ as a function of parameters $p$
5. The observation equation $Q$ as a fixed matrix.
6. The Cholesky factor of the observation errors, $\Omega$ as a function of parameters $p$.  At this point only a diagonal matrix is supported.
7. Either the steady state equations for all of $y$ and $x$ in closed form as a function of $p$, or initial conditions for the nonlinear solution to solve for the steady state as functions of $p$
## Defining Models
Install this package with `] add DifferentiableStateSpaceModels`, then the full code to create the RBC model is

```julia
∞ = Inf
@variables α, β, ρ, δ, σ, Ω_1
@variables t::Integer, k(..), z(..), c(..), q(..)

x = [k, z] # states
y = [c, q] # controls
p = [α, β, ρ, δ, σ, Ω_1] # parameters

H = [1 / c(t) - (β / c(t + 1)) * (α * exp(z(t + 1)) * k(t + 1)^(α - 1) + (1 - δ)),
     c(t) + k(t + 1) - (1 - δ) * k(t) - q(t),
     q(t) - exp(z(t)) * k(t)^α,
     z(t + 1) - ρ * z(t)]  # system of model equations

# analytic solutions for the steady state.  Could pass initial values and run solver and use initial values with steady_states_iv
steady_states = [k(∞) ~ (((1 / β) - 1 + δ) / α)^(1 / (α - 1)),
                 z(∞) ~ 0,
                 c(∞) ~ (((1 / β) - 1 + δ) / α)^(α / (α - 1)) -
                        δ * (((1 / β) - 1 + δ) / α)^(1 / (α - 1)),
                 q(∞) ~ (((1 / β) - 1 + δ) / α)^(α / (α - 1))]


Γ = [σ;;] # matrix for the 1 shock.  The [;;] notation just makes it a matrix rather than vector in julia
η = [0; -1;;] # η is n_x * n_ϵ matrix.  The [;; ]notation just makes it a matrix rather than vector in julia

# observation matrix.  order is "y" then "x" variables, so [c,q,k,z] in this example
Q = [1.0 0  0   0; # select c as first "z" observable
     0   0  1.0 0] # select k as second "z" observable

# diagonal cholesky of covariance matrix for observation noise (so these are standard deviations).  Non-diagonal observation noise not currently supported
Ω = [Ω_1, Ω_1]

# Generates the files and includes if required.  If the model is already created, then just loads
overwrite_model_cache  = true
model_rbc = @make_and_include_perturbation_model("rbc_notebook_example", H, (; t, y, x, p, steady_states, Γ, Ω, η, Q, overwrite_model_cache))
```

After generation of the model, they can be included as any other julia files in your code (e.g. `include(joinpath(pkgdir(DifferentiableStateSpaceModels), ".function_cache","my_model.jl"))`) or moved somewhere more convenient.

Inclusion through the `@make_and_include_perturbation_model` creates the model automatically; after direct inclusion through a julia file, you can create a model with `m = PerturbationModel(Main.my_model)`.

## Solving Perturbations

Assuming the above model was created and loaded in one way or another as `m`
```julia
p_f = (ρ = 0.2, δ = 0.02, σ = 0.01, Ω_1 = 0.01) # Fixed parameters
p_d = (α = 0.5, β = 0.95) # Pseudo-true values
sol = generate_perturbation(model_rbc, p_d, p_f) # Solution to the first-order RBC
sol_2 = generate_perturbation(model_rbc, p_d, p_f, Val(2)); # Solution to the second-order RBC
@show sol.retcode, sol_2.retcode, verify_steady_state(m, p_d, p_f) # the final call checks that the analytically provided steady-state solution is correct
```

The perturbation solution (in the canonical form described in the top section) can be queried from the resulting solution.  A few examples for the first order solution are below,
```julia
@show sol.y, sol.x  # steady states y_ss and x_ss  These are the values such that y ≡ ŷ + sol.y and x ≡ x̂ + sol.x
@show sol.g_x # the policy
@show sol.A, sol.B # the evolution equation of the state, so that x̂' = A x̂ + B ϵ
@show sol.C, sol.D; # the evolution equation of the state, so that z = C x̂ + ν  with variance of ν as D D'.
@show sol.x_ergodic_var; # covariance matrix of the ergodic distribution of x̂, which is mean zero since x̂ ≡ x - x_ss
```

## Functions of Perturbation Solutions (and their Derivatives)
The core feature of this library is to enable gradients of the perturbation solutions with respect to parameters (i.e., anything in the `p_d` vector).  To show this, we will construct a function which uses the resulting law of motion and finds the gradient of the results with respect to this value.
```julia
function IRF(p_d, ϵ_0; m, p_f, steps)
    sol = generate_perturbation(m, p_d, p_f) # First-order perturbation by default, pass Val(2) as additional argument to do 2nd order.
    x = sol.B * ϵ_0 # start after applying impulse with the model's shock η and Γ(p)
    for _ in 1:steps
        # A note: you cannot use mutating expressions here with most AD code.  i.e. x .= sol.A * x  wouldn't work
        # For more elaborate simuluations, you would want to use DifferenceEquations.jl in practice
        x = sol.A * x # iterate forward using first-order observation equation        
    end    
    return [0, 1]' * sol.C * x # choose the second observable using the model's C observation equation since first-order
end

m = model_rbc  # ensure notebook executed above
p_f = (ρ=0.2, δ=0.02, σ=0.01, Ω_1=0.01) # not differentiated 
p_d = (α=0.5, β=0.95) # different parameters
steps = 10 # steps ahead to forecast
ϵ_0 = [1.0] # shock size
IRF(p_d, ϵ_0; m, p_f, steps) # Function works on its own, calculating perturbation
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

However, the real benefit is that this function can itself be differentiated, to find gradients with respect to the deep parameters `p_d` and the shock `ϵ_0`
```julia
# Using the Zygote auto-differentiation library already loaded above
p_d = (α=0.5, β=0.95) # different parameters
ϵ_0 = [1.0] # shock size
IRF_grad = gradient((p_d, ϵ_0) -> IRF(p_d, ϵ_0; m, p_f, steps), p_d, ϵ_0) # "closes" over the m, p_f, and steps to create a new function, and differentiates it with respect to other arguments
```
## Simulating Data with DifferenceEquations
## Simulating with DifferenceEquations.jl
The manual iteration of the state-space model from the perturbation solution is possible, but can be verbose and difficult to achieve efficiency for gradients.  One benefit of this package is that it creates state-space models in a form consistent with [DifferenceEquations.jl](https://github.com/SciML/DifferenceEquations.jl) which can be easily simulated, visualized, and estimated.

To do this, we will calculate a perturbation solution then simulate it for various `x_0` drawn from the ergodic solution.

```julia
p_f = (ρ = 0.2, δ = 0.02, σ = 0.01, Ω_1 = 0.01) # Fixed parameters
p_d = (α = 0.5, β = 0.95) # Pseudo-true values
m = model_rbc  # ensure notebook executed above
sol = generate_perturbation(m, p_d, p_f) # Solution to the first-order RBC

# Simulate T observations from a random initial condition
T = 20

# draw from ergodic distribution for the initial condition
x_iv = MvNormal(sol.x_ergodic_var)
problem = LinearStateSpaceProblem(sol, x_iv, (0, T))
plot(solve(problem))
```

The `LinearStateSpaceProblem` type is automatically constructed from the underlying perturbation.  However, we can override any of these options, or pass in our own noise rather than simulate it for a particular experiment

```julia
noise = Matrix([1.0; zeros(T-1)]') # the ϵ shocks are "noise" in DifferenceEquations for SciML compatibility
x_iv = [0.0, 0.0]  # can pass in a single value rather than a distribution 
problem = LinearStateSpaceProblem(sol, x_iv, (0, T); noise)
plot(solve(problem))
```

To demonstrate the composition of gradients between DifferenceEquations and DifferentiableStateSpaceModels lets adapt this function to simulates an impulse with fixed noise and looks at the final observable

```julia
function last_observable(p_d, noise, x_iv; m, p_f, T)
    sol = generate_perturbation(m, p_d, p_f)
    problem = LinearStateSpaceProblem(sol, x_iv, (0, T);noise, observables_noise = nothing)  # removing observation noise
    return solve(problem).z[end][2]  # return 2nd argument of last observable
end
T = 100
noise = Matrix([1.0; zeros(T-1)]') # the ϵ shocks are "noise" in DifferenceEquations for SciML compatibility
x_iv = [0.0, 0.0]  # can pass in a single value rather than a distribution
p_f = (ρ = 0.2, δ = 0.02, σ = 0.01, Ω_1 = 0.01) # Fixed parameters
p_d = (α = 0.5, β = 0.95) # Pseudo-true values
m = model_rbc  # ensure notebook executed above
last_observable(p_d, noise, x_iv; m, p_f, T)
```

And, as before, we can calculate gradients with respect to the underlying `p_d` parameters, but also with respect to the noise which will demonstrate a key benefit of these methods, as they can let us do a joint likelihood of the latent variables in cases where they cannot be easily marginalized out (e.g., non-Gaussian or nonlinear).  Note that the dimensionality of this gradient is over 100.

```julia
gradient((p_d, noise, x_iv) -> last_observable(p_d, noise, x_iv; m, p_f, T), p_d, noise, x_iv)
```

Finally, we can use a simple utility functions to investigate an IRF.
```julia
p_f = (ρ = 0.2, δ = 0.02, σ = 0.01, Ω_1 = 0.01) # Fixed parameters
p_d = (α = 0.5, β = 0.95) # Pseudo-true values
m = model_rbc
sol = generate_perturbation(model_rbc, p_d, p_f)
ϵ0 = [1.0]
T = 10
sim = irf(sol, ϵ0, T)
plot(sim)
```

## Sneak Peak at SciML Compatible Functionality
Finally, there are a variety of features of SciML which are supported.  For example, parallel simulations of ensembles and associated summary statistics.

```julia
# Simulate multiple trajectories with T observations
trajectories = 40
x_iv = MvNormal(sol.x_ergodic_var)
problem = LinearStateSpaceProblem(sol, x_iv, (0, T))

# Solve multiple trajectories and plot an ensemble
ensemble_results = solve(EnsembleProblem(problem), DirectIteration(), EnsembleThreads();
                 trajectories)
summ = EnsembleSummary(ensemble_results)  # see SciML documentation.  Calculates median and other quantles automatically.
summ.med # median values for the "x" simulated ensembles

plot(summ, fillalpha= 0.2) # plots by default show the median and quantiles of both variables.  Modifying transparency as an example
```

## Calculate sequence of observables
We can use the underlying state-space model to easily simulate states and observables

```julia
# Simulate T observations
T = 20

p_f = (ρ = 0.2, δ = 0.02, σ = 0.01, Ω_1 = 0.01) # Fixed parameters
p_d = (α = 0.5, β = 0.95) # Pseudo-true values
sol = generate_perturbation(model_rbc, p_d, p_f) # Solution to the first-order RBC

x_iv = MvNormal(sol.x_ergodic_var) # draw initial conditions from the ergodic distribution
problem = LinearStateSpaceProblem(sol, x_iv, (0, T))
sim = solve(problem, DirectIteration())
ϵ = sim.W # store the underlying noise in the simulation

# Collapse to simulated observables as a matrix  - as required by current DifferenceEquations.jl likelihood
# see https://github.com/SciML/DifferenceEquations.jl/issues/55 for direct support of this datastructure
z_rbc = hcat(sim.z...) 
```
