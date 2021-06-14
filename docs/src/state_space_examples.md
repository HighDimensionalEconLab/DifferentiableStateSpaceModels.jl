# State Space Examples


This tutorial will introduce you to the functionality for solving discrete-time models with an observation equation. Other
introductions can be found by [checking out DiffEqTutorials.jl](https://github.com/JuliaDiffEq/DiffEqTutorials.jl). The [Linear State Space Examples](@ref linear_state_space_examples) provides specialized algorithms when the models are linear and linear-quadratic.


!!! note

    This tutorial assumes you have read the [Ordinary Differential Equations tutorial](@ref ode_example) and the mathematics of the [State Spaces Models](@ref state_space_types)


## Example 1: General State Space Model

The canonical form of the nonlinear model with additive noise is

$$
u_{n+1} = f(u_n,p,t_n) + g(u_n,p,t_n) w_{n+1}
$$
with
$$
z_n = h(u_n, p, t_n) +  v_n
$$

In additional, we will tend to use $v_n \sim D$.

Frequently, the $D = N(0, R)$ and $w_{n+1} \sim N(0,I)$.

Creating a `StateSpaceProblem` and simulating it for a simple, linear equation.

```julia
using DifferentialEquations, LinearAlgebra, Distributions, Random, Plots

f(u, p, t) = [p[1] 0.0; 0.1   0.7] * u
g(u, p, t) = Diagonal([0.1, p[2]])
h(u, p, t) = [0.5 0.5] * u
R = [p[3]]  # one observable

# Simulate data
p = [0.8, 0.05, 0.01]
T = 10
u₀=[0.0, 0.1]
tspan = (0, T)
D = MvNormal(0, R)
prob = StateSpaceProblem(f,g,u₀,tspan,p; h = h, D = D)

# TODO:  TEMPORARY TESTING CODE!  REPLACE WITH `solve` BELOW WHEN FUNCTIONAL
draw_from_prior(u0::AbstractVector) = u0  # degenerate prior
draw_from_prior(u0) = rand(u0) # otherwise, assume it is a distribution
perfect_observation(u,p,t) = u
function simulate_example(f, g, u0_prior, p, tspan; c = perfect_observation, d = nothing, noise = nothing)

    T = tspan[2] - tspan[1]  # hardcoded right now

    # Draw from prior (or use as a degeenrate initial condition)
    u0 = draw_from_prior(u0_prior)
    B_0 = g(u0, p, 0)
    C_0 = c(u0, p, 0)
    n_u, n_w = size(B_0)
    n_z = length(C_0)

    w = isnothing(noise) ? [randn(n_w) for _ in 1:T] : noise

    u = u0
    v = isnothing(d) ? zeros(n_z) : rand(MvNormal(zeros(n_z), d(u, p, 0)))
    z =  c(u, p, 0) .+ v
    us = [u]
    zs = [z]
    vs = [v]
    for t in 1:T
        u = f(u, p, t-1) .+ g(u, p, t-1) * w[t]
        v = isnothing(d) ? zeros(n_z) : rand(MvNormal(zeros(n_z), d(u, p, t)))
        z =  c(u, p, t) .+ v
        push!(us, u)
        push!(zs, z)
        push!(vs, v)
    end
    t = 0:T  # hardcoded
    t_w = 1:T  # none time 0.
    return (u=DiffEqArray(us, t), z=DiffEqArray(zs,t), v=DiffEqArray(v,t), w=DiffEqArray(w,t_w))
end
u, z, v, w = simulate_example(f, g, u₀, p,tspan; c = c, d = d)
#u, z, v, w = simulate_example(f, g, MvNormal(u₀, Diagonal([0.01, 0.02])), p,tspan; c = c, d = d)
#sol = z
## END REMOVE
```

We can `solve` the model to simulate a path.  Since we have not provided the $w_t$ or $v_t$ sequence, it will simulate it using the default Gaussian draws.

```julia
sol = solve(prob)  # default algorithm is just iteration
@show sol[1]  # This is the observation at the first time period.
@show sol[1,:]  # the observation of the first value for all periods

plot(sol)  # or plot(sol.z)
```

The `u` state is hidden. To access the simulated values,

```julia
plot(sol.u)
```

!!! note

    Since the output of the state space model is the observables, `sol[i,j]` refers to `sol.z[i,j]` instead of `sol.u[i,j]`

Assuming that you chose to save the `noise` and the `observational_noise` (i.e. `solve(prob; save_noise = true, save_observational_noise = true`) are the defaults, then you can also access them through

```julia
plot(sol.W)  # noise on the evolution equation
plot(sol.V)  # observational noise, if it exists.
````


While not necessary here, if other algorithms are available they are passed as with other SciML solve arguments (e.g. the default is `solve(prob, AdditiveNoiseIteration())`).

## Example 2: Simulating Given Noise

We can simulate the problem given a particular noise process.  For now, we will assume that there is no observational noise - i.e., using the default `D = nothing`.


```julia
prob = StateSpaceProblem(f,g,u₀,tspan,p; h = h, noise = W)  # removed observational noise, fixing `W`
```

We can simulate this given an external `noise` argument which replaces the default random gaussian used in `solve`.  Here, we provide an impulse at the 2nd time period.
```julia
W = [[0.0, 0.0]]
W = vcat([[0.0, 0.0]],
         [[1.0, 1.0]], # impulse at 2nd period
         [zeros(2) for _ in 3:T])
sol = solve(prob)

plot(sol)  # dynamics with impulse at 2nd period
```

## Example 3: Differentiating the Solution Given Noise

The `StateSpace` model and its `solve` support differentiation.  For example, we can simulate a sequence of observables using our `p` and a random `W`.

```julia
W = [rand(2) for _ in 1:T]
prob = StateSpaceProblem(f,g,u₀,tspan,p; h = h, noise = W)  # no observational noise
z = solve(prob)  # simulated observables
```

Then, we will define a function which stakes in the `W` and `p` and calculates the sum-of-squares deviation from the data.

```julia
function sum_squared_loss(z, u₀, p, W)
   prob = StateSpaceProblem(f,g,u₀,tspan,p; h = h, noise = W)
   sol = solve(prob)
   return norm(sol .- z, 2)  # sum of squares of deviations.
end
sum_squared_loss(z, u₀, p, W)
```
As expected, the deviation is zero since we fed in the exact initial condition, paramere, and noise vectors.

Now let us see the gradient of this function as we change `p` and/or `W` using Zygote's reverse-mode AD.

```julia
using Zygote
grad_ss_loss = gradient((p, W) -> sum_squared_loss(z, u₀, p, W), p, W)
```
In the above, we have bound the initial condition and the data to make it a function of only `p` and `W` and then differentiated the loss function with respect to the parameters and noise.

## Example 4: Likelihood and Gradient Given Noise

Instead of an ad-hoc loss function above, we might want to calculate the likelihood of the data given the solution to the model.  This could be done by using the `d(..)` for each period and comparing the log-likelihood of the observables, for example.

```julia

D = MvNormal(zeros(1), [p[3]])  # one observable
function manual_logpdf_given_noise(z, u₀, p, W)
   prob = StateSpaceProblem(f,g,u₀,tspan,p; h = h, noise = W)  # not simulating with observation error
   sol = solve(prob)
   loglik = 0.0
   D = MvNormal(zeros(2), R)
   for (n, z) in enumerate(sol)
      loglik += logpdf(D, z[n] .- sol[n])  # log pdf of observable
    end
    return loglik
end
grad_manual_logpdf = gradient((p, W) -> manual_logpdf_given_noise(z, u₀, p, W), p, W)
```

Alteratively, we can leave the observation error in the model itself, and use a filter to calculate the log likelihood.
```julia
u₀_prior = u₀ # degenerate prior
prob = StateSpaceProblem(f,g,u₀_prior,tspan,p; h = h, D = D, noise = W, observables = z) # gaussian observation error, u₀ constant prior

# The NoiseConditionalFilter calculates the likelihood given the noise
sol = solve(prob, NoiseConditionalFilter())  # runs the filter given the prior
sol.logpdf  # log likelihood, additional return type from a filter algorithm
```
Which can also be differentiated,
```julia
function logpdf_given_noise(z, u₀, p, W)
   prob = StateSpaceProblem(f,g,u₀,tspan,p; h = h, D=D, noise = W, observables = z)  # Gaussian noise
   return solve(prob, NoiseConditionalFilter()).logpdf
end
gradient((p, W) -> logpdf_given_noise(z, u₀, p, W), p, W)
```
Note that when we are creating the StateSpaceProblem, the `u₀` for filtering is interpreted as a prior (and a degenerate prior if it is constant).

## Example 5: Ensemble

!!! note

    Since the `y` is the observable in the state space, all ensemble features use the `sol = sol.z`.  For now, a limitation is that length of `y` and `u` must be the same.  To access the underlying `u`, add it to the observation equation, or implement an `output_func` callback.

First, we will make a small change to our model to add a second observable.

```julia

f(u, p, t) = [p[1] 0.0; 0.1   0.7] * u
g(u, p, t) = Diagonal([0.1, p[2]])
c(u, p, t) = [0.5 0.5; 0.0 1.0] * u
d(u, p, t) = [p[3]; 0.01]  # two observables

# Simulate data
p = [0.8, 0.05, 0.01]
T = 10
u₀=[0.0, 0.1]
tspan = (0, T)
D = MvNormal(zeros(2), Diagonal([p[3]; 0.01]))
prob = StateSpaceProblem(f,g,u₀,tspan,p; h = h, D = D)
sol = solve(prob)
plot(sol, vars = [1])  # Plot the first variable.
```

Rather than plotting a single simulation, we can calculate [ensembles](https://diffeq.sciml.ai/stable/features/ensemble/) in parallel,
```julia
ensemble_prob = EnsembleProblem(prob)  # could add output_func, reduction, u_init, etc.
sol = solve(ensemble_prob,EnsembleThreads(),trajectories=1000)
```

Given the ensemble solution, summary statistics and plotting are standard

```julia
summ = EnsembleSummary(sol)  # the mean, 5th, and 95th quantiles by default
plot(summ)
```

## Example 6: Nonlinear Filtering
**NOTE: SPECULATIVE, FUTURE VERSIONS!**  
While the previous examples calculate the likelihood conditioning on the `w_t` process, to calculate the marginal likelihood over all `w_t, v_t` shocks we need to choose a a nonlinear filter given a prior on the initial state.

!!! note

    See [linear state space examples](@ref linear_state_space_examples) for cases with exact likelihoods given the Kalman Filter.


Here, we will setup a `StateSpaceProblem` with a prior, and calculate the likeihood of the observables using the `UnscentedKalmanFilter` (which will be exact in this case since we provide a linear gaussian model).

TODO: Do we need to pass in the jacobian for the UKF?

```julia
f(u, p, t) = [p[1] 0.0; 0.1   0.7] * u
g(u, p, t) = Diagonal([0.1, p[2]])
h(u, p, t) = [0.5 0.5] * u


# Simulate data
p = [0.8, 0.05, 0.01]
T = 10
u₀_prior = MvNormal([0.0, 0.1], Diagonal([0.01, 0.01])
D = MvNormal([0.0], [p[3]])  # one observable
tspan = (0, T)

prob = StateSpaceProblem(f,g,u₀_prior,tspan,p; h = h, D = D, observables = z)  # prior for initial condition
# Simulate some observables, where u0 is drawn from the prior
z = solve(prob)

# calculate the likelihood of those observations with the UKF
sol = solve(prob, UnscentedKalmanFilter(); save_everystep = true)
@show sol.logpdf
```
Or, we can use the `sol` to extract the sequence of posteriors.

The default from `sol[i]` is the mean of the i'th posterior, where `sol.posterior[i]` gives the posterior distribution, which has a `cov` since it is of type `MvNormal` here.
```julia
plot(sol.u)  #  posterior mean?  Or smoothed?  Careful?
plot(sol.t, [cov(posterior) for posterior in sol])  # posterior covariance

# TODO: or add a recipe of some sort?  posterior mean and the 5th/95th quantiles around it?
plot(sol)
```
To run the filter without storing the intermediate values (e.g. if you only need the log likelihood), use the standard `solve` options - where `save_everystep = false` is the default.
```julia
sol = solve(prob, UnscentedKalmanFilter(); save_everystep = false, observables = z)
sol.t  # Only the beginning and end.
```

If we wish to run a smoother, we can replace `UnscentedKalmanFilter()` with `UnscentedKalmanSmoother()`.

As before, we are able to differentiate the filter.
```julia
function logpdf_UKF(z, u₀, p)
   prob = StateSpaceProblem(f,g,u₀,tspan,p; h = h, D=D)  # Gaussian noise
   return solve(prob,UnscentedKalmanFilter(); observables = z, save_everystep = false).logpdf
end
gradient(p -> logpdf_UKF(z, u₀, p), p)
```
Note that we have bound the initial prior but could have parameterized and differentiated that as well.