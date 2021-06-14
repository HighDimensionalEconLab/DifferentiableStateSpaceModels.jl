using BenchmarkTools, Zygote, LinearAlgebra, RecursiveArrayTools, DistributionsAD, Turing,
      Parameters
using Zygote: @adjoint

struct LinearStateSpaceModel{T1,T2,T3,T4,T5,T6}
    A::T1
    B::T2
    C::T3
    D::T4
    x_mean::T5
    x_variance::T6
end
@adjoint function LinearStateSpaceModel(A, B, C, D, x_mean, x_variance)
    return LinearStateSpaceModel(A, B, C, D, x_mean, x_variance), Δ -> (Δ.A, Δ.B, Δ.C, Δ.D, Δ.x_mean, Δ.x_variance)
end

function Kalman(mod::LinearStateSpaceModel, obs)
    T = length(obs)
    T1 = promote_type(eltype(mod.A), eltype(mod.B), eltype(mod.C), eltype(mod.D))
    z = Vector{Vector{T1}}(undef, T)
    V = Vector{Matrix{T1}}(undef, T)
    cur_x = deepcopy(mod.x_mean)
    cur_P = deepcopy(mod.x_variance)
    for i in 1:T
        # Kalman iteration
        cur_x = mod.A * cur_x
        cur_P = mod.A * cur_P * mod.A' + mod.B * mod.B'
        z[i] = mod.C * cur_x
        V[i] = mod.C * cur_P * mod.C + mod.D * mod.D'
        V[i] = (V[i] + V[i]') / 2.0 # make sure V is symmetric -- Hermitian form
        cur_x = cur_x + cur_P' * mod.C' / V[i] * (obs[i] - z[i])
        cur_P = cur_P - cur_P' * mod.C' / V[i] * mod.C * cur_P
    end
    return z, V
end

function _Kalman_zygote(mod::LinearStateSpaceModel, obs)
    T = length(obs)
    z = Zygote.Buffer(Vector{Vector{Float64}}(undef, T))
    V = Zygote.Buffer(Vector{Matrix{Float64}}(undef, T))
    cur_x = Zygote.Buffer(Vector{Vector{Float64}}(undef, T))
    cur_P = Zygote.Buffer(Vector{Matrix{Float64}}(undef, T))
    for i in 1:T
        # Kalman iteration
        if (i > 1)
            cur_x[i] = mod.A * cur_x[i - 1]
            cur_P[i] = mod.A * cur_P[i - 1] * mod.A' + mod.B * mod.B'
        else
            cur_x[1] = mod.A * mod.x_mean
            cur_P[1] = mod.A * mod.x_variance * mod.A' + mod.B * mod.B'
        end
        z[i] = mod.C * cur_x[i]
        V[i] = mod.C * cur_P[i] * mod.C' + mod.D * mod.D'
        V[i] = (V[i] + V[i]') / 2.0 # make sure V is symmetric -- Hermitian form
        cur_x[i] = cur_x[i] + cur_P[i]' * mod.C' / V[i] * (obs[i] - z[i])
        cur_P[i] = cur_P[i] - cur_P[i]' * mod.C' / V[i] * mod.C * cur_P[i]
    end
    copy(cur_x)
    copy(cur_P)
    return copy(z), copy(V)
end

@adjoint Kalman(mod, obs) = Zygote.pullback(_Kalman_zygote, mod, obs)

function Kalman_Likelihood_Only(mod::LinearStateSpaceModel, obs)
    T = length(obs)
    T1 = promote_type(eltype(mod.A), eltype(mod.B), eltype(mod.C), eltype(mod.D))
    loglik = zero(T1)
    cur_x = deepcopy(mod.x_mean)
    cur_P = deepcopy(mod.x_variance)
    for i in 1:T
        # Kalman iteration
        cur_x = mod.A * cur_x
        cur_P = mod.A * cur_P * mod.A' + mod.B * mod.B'
        z = mod.C * cur_x
        V = mod.C * cur_P * mod.C + mod.D * mod.D'
        V = (V + V') / 2.0 # make sure V is symmetric -- Hermitian form
        loglik += logpdf(MvNormal(z, V), obs[i])
        cur_x = cur_x + cur_P' * mod.C' / V * (obs[i] - z)
        cur_P = cur_P - cur_P' * mod.C' / V * mod.C * cur_P
    end
    return loglik
end

function _Kalman_Likelihood_Only_zygote(mod::LinearStateSpaceModel, obs)
    T = length(obs)
    # z = Zygote.Buffer(Vector{Vector{Float64}}(undef, T))
    # V = Zygote.Buffer(Vector{Matrix{Float64}}(undef, T))
    cur_x = Zygote.Buffer(Vector{Vector{Float64}}(undef, T))
    cur_P = Zygote.Buffer(Vector{Matrix{Float64}}(undef, T))
    loglik = 0.0
    for i in 1:T
        # Kalman iteration
        if (i > 1)
            cur_x[i] = mod.A * cur_x[i - 1]
            cur_P[i] = mod.A * cur_P[i - 1] * mod.A' + mod.B * mod.B'
        else
            cur_x[1] = mod.A * mod.x_mean
            cur_P[1] = mod.A * mod.x_variance * mod.A' + mod.B * mod.B'
        end
        z = mod.C * cur_x[i]
        V = mod.C * cur_P[i] * mod.C' + mod.D * mod.D'
        # V[i] = (V[i] + V[i]') / 2.0 # make sure V is symmetric -- Hermitian form
        loglik += logpdf(MvNormal(z, V), obs[i])
        cur_x[i] = cur_x[i] + cur_P[i]' * mod.C' / V * (obs[i] - z)
        cur_P[i] = cur_P[i] - cur_P[i]' * mod.C' / V * mod.C * cur_P[i]
    end
    copy(cur_x)
    copy(cur_P)
    return loglik
end

@adjoint function Kalman_Likelihood_Only(mod, obs)
    return Zygote.pullback(_Kalman_Likelihood_Only_zygote, mod, obs)
end

@adjoint diagm(x::AbstractVector) = diagm(x), dy -> (diag(dy),)
@adjoint function diagm(pr::Pair)
    return diagm(pr), dy -> ((first = nothing, second = diag(dy, first(pr))))
end
Zygote.refresh()

function mymodelbuilder(p, x_mean, x_var)
    N = div(length(p), 4)
    p_A = p[1:N]
    p_B = p[(N + 1):(2N)]
    p_C = p[(2N + 1):(3N)]
    p_D = p[(3N + 1):(4N)]
    return LinearStateSpaceModel(diagm(p_A), diagm(p_B), diagm(p_C), diagm(p_D), x_mean,
                                 x_var)
end

# then call it with something like
function generate_test_data(N, T, x_0_prod = 1.1, gamma = 0.8, shock_var = 0.01,
                            obs_coeff = 1.0, obs_var = 0.01, x_mean_val = 1.0,
                            x_var_val = 0.02)
    p_A = fill(gamma, N)
    p_B = fill(shock_var, N)
    p_C = fill(obs_coeff, N)
    p_D = fill(obs_var, N)
    p = [p_A; p_B; p_C; p_D]
    x_mean = fill(x_mean_val, N)
    x_var = Symmetric(diagm(fill(x_var_val, N)))
    mod = mymodelbuilder(p, x_mean, x_var)
    # starting point
    x = Vector{Vector{Float64}}(undef, T)
    obs = Vector{Vector{Float64}}(undef, T)
    # Solving deterministically, but slightly off prior mean.
    x[1] = x_0_prod * x_mean
    obs[1] = mod.C * x[1]
    for i in 2:T
        x[i] = mod.A * x[i - 1]
        obs[i] = mod.C * x[i]
    end
    return (; p, x, obs, N, T, x_mean, x_var)
end

function Kalman_likelihood(p, x_mean, x_var, obs)
    mod = mymodelbuilder(p, x_mean, x_var)
    z, V = Kalman(mod, obs)
    loglik = zero(eltype(p))
    for i in eachindex(z, V)
        loglik += logpdf(MvNormal(z[i], V[i]), obs[i])
    end
    return loglik
end

function Kalman_likelihood_option1(p, x_mean, x_var, obs)
    mod = mymodelbuilder(p, x_mean, x_var)
    loglik = Kalman_Likelihood_Only(mod, obs)
    return loglik
end

function benchmark_kalman(dat; calculate_forward = true)
    @unpack N, T, p, obs, x_mean, x_var = dat
    println("Primal Calculation")
    @btime Kalman_likelihood($p, $x_mean, $x_var, $obs)
    @btime Kalman_likelihood_option1($p, $x_mean, $x_var, $obs)
    println("Zygote Reverse | only p")
    @btime gradient(p -> Kalman_likelihood(p, $x_mean, $x_var, $obs), $p)
    @btime gradient(p -> Kalman_likelihood_option1(p, $x_mean, $x_var, $obs), $p)
    println("Zygote Reverse | all")
    @btime gradient((p, x_mean, x_var) -> Kalman_likelihood(p, x_mean, x_var, $obs), $p,
                    $x_mean, $x_var)
    @btime gradient((p, x_mean, x_var) -> Kalman_likelihood_option1(p, x_mean, x_var, $obs),
                    $p, $x_mean, $x_var)
    if calculate_forward
        println("Zygote Forward | only p")
        @btime Zygote.forward_jacobian(p -> Kalman_likelihood(p, $x_mean, $x_var, $obs), $p)
        @btime Zygote.forward_jacobian(p -> Kalman_likelihood_option1(p, $x_mean, $x_var,
                                                                      $obs), $p)
    end
    return nothing
end

## Usage
T = 200
n_1 = 5
n_2 = 50
n_3 = 100
dat_1 = generate_test_data(n_1, T)
dat_2 = generate_test_data(n_2, T)
dat_3 = generate_test_data(n_3, T)
benchmark_kalman(dat_1)
benchmark_kalman(dat_2)
benchmark_kalman(dat_3, calculate_forward = false)
