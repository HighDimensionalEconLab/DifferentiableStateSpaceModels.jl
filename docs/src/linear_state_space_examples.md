# Linear State Space Examples

This tutorial provides additional features for linear models

!!! note

    This tutorial assumes you have read the [Ordinary Differential Equations tutorial](@ref ode_example), the mathematics of the [State Spaces Models](@ref state_space_types), and the more general [State Space Examples](@ref state_space_examples) tutorial.
    
    In addition, this package uses the `MatrixFreeOperator` in [DiffEqOperators.jl](https://github.com/SciML/DiffEqOperators.jl) to wrap the general linear equations.


The canonical form of the linear model is

$$
u_{n+1} = A(p, t_n) u_n + B(p,t_n) w_{n+1}
$$
with
$$
z_n = C(p, t_n) u_n +  v_n
$$

where we will tend to use $v_n \sim N(0, R)$ and $w_{n+1} \sim N(0,I)$.


However, we will implement  on the the linear-time invariant (LTI) version,
$$
u_{n+1} = A u_n +  B w_{n+1}
$$
with
$$
z_n = C u_n +  v_n
$$
and $v_n \sim N(0, R)$ where $A, B, C$ and $R$ may be parameterized by $p$.

## Example 1: Linear (and Time-Invariant) State Space Model

Creating a `LinearStateSpaceProblem` and simulating it for a simple, linear equation.

```julia
using DifferentialEquations, LinearAlgebra, Distributions, Random, Plots
A = [0.8 0.0; 0.1 0.7]
B = Diagonal([0.1, 0.5])
C = [0.5 0.5] # one observable
R = [0.01]

# Simulate data
T = 10
u₀=[0.0, 0.1]
tspan = (0, T)

prob = LinearStateSpaceProblem(A,B,u₀,tspan; C = C, R = R)
```

We can `solve` the model to simulate a path.  Since we have not provided the $w_t$ or $v_t$ sequence, it will simulate it using the default Gaussian draws.  The use of the algorithm `LinearGaussian()` is a specialization

```julia
sol = solve(prob, LinearGaussian())  # default algorithm is linear-gaussian iteration
@show sol[1]  # This is the observation at the first time period.
@show sol[1,:]  # the observation of the first value for all periods

plot(sol)  # or plot(sol.z)
```

The `u` state is not-observable in the primary output.  To access the simulated values,

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

## Example 2: Kalman Filter for LTI System

Here, we will setup a `LinearStateSpaceProblem` with a prior, and calculate the likelihood of the observables using the `KalmanFilter` (which will be exact in this case since we provide a linear gaussian model).

```julia
p = [0.8, 0.05, 0.01]
T = 10
u₀_prior = MvNormal([0.0, 0.1], Diagonal([0.01, 0.01])
tspan = (0, T)

prob = LinearStateSpaceProblem(A, B, u₀_prior, tspan; C = C, R = R)  # prior for initial condition
# Simulate some observables, where u0 is drawn from the prior
sol_sim = solve(prob, LinearGaussian())
z = sol_sim
```

Then, attach the observables, and calculate the likelihood,

```julia
prob = LinearStateSpaceProblem(A, B, u₀_prior, tspan; C = C, R = R, observables = z)
sol = solve(prob, KalmanFilter(); save_everystep = true)
@show sol.logpdf
```

Or, we can use the `sol` to extract the sequence of posteriors.

```julia
# Or to extract the posteriors
plot(sol.t, [mean(posterior) for posterior in sol.posteriors])
plot(sol.t, [cov(posterior) for posterior in  sol.posteriors])  # posterior covariance

# TODO: add recipe of some sort?  posterior mean and the 5th/95th quantiles around it?
plot(sol)
```

To run the filter without storing the intermediate values (e.g. if you only need the log likelihood), use the standard `solve` options - where `save_everystep = false` is the default.

```julia
sol = solve(prob, KalmanFilter(); save_everystep = false, save_posteriors = false)
sol.t  # does not save all time periods and saves no posteriors
```

Note: If we wish to run a smoother, we can replace the algorithm with `KalmanSmoother()`, which will update the `u` and `posteriors` on its back-pass.


## Example 3: Differentiating the Kalman Filter

We are then able to differentiate the filter.
```julia
function logpdf_KF(z, u₀, p)
    # parameterize matrices
    A = [p[1] 0.0; 0.1 0.7]
    B = Diagonal([0.1, p[2]])
    C = [0.5 0.5]
    R = [p[3]]  # one observable

    prob = LinearStateSpaceProblem(A, B, u₀,tspan; C=C, R=R)  # Gaussian noise
    return solve(prob, KalmanFilter(), save_everystep = false).logpdf
end
gradient(p -> logpdf_KF(z, u₀, p), p)
```

Note that we have bound the initial prior but could have parameterized and differentiated that as well.
