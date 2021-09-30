

# NOTE: this model file is run from `/deps/generate_models.jl` or in the build phase.

# Model as in https://www.sas.upenn.edu/~jesusfv/ARE_Estimation.pdf

function FVGQ20()
    # Parameter values
    # δ = 0.025
    # ε = 10
    # ϕ = 0
    # γ2 = 0.001

    # Estimated parameters, taken from FV(2010), Table 3, p. 38, median estimate parameters
    # β = 0.998
    # h = 0.97
    # ϑ = 1.17
    # κ = 9.51
    # α = 0.21

    # θp = 0.82
    # χ = 0.63
    # γR = 0.77
    # γy = 0.19

    # γΠ = 1.29
    # Πbar = 1.01
    # ρd  = 0.12
    # ρφ = 0.93
    # σ_A = -3.97
    # σ_d = -1.51

    # σ_φ = -2.36
    # σ_μ = -5.43
    # σ_m = -5.85
    # Λμ = 3.4e-3
    # ΛA = 2.8e-3

    # 7 variable relabelings
    # 7 observations
    # 38 equations

    ∞ = Inf
    # 28 parameters. 23 of them to estimate, 5 parameters as constants
    @variables δ, ε, ϕ, γ2, β, h, ϑ, κ, α, θp, χ, γR, γy, γΠ, Πbar, ρd, ρφ, ρg, g_bar, σ_A, σ_d, σ_φ, σ_μ, σ_m, σ_g, Λμ, ΛA, Ω_ii
    # States and Controls, 14 + 24 = 38
    @variables t::Integer, c_m(..), x_m(..), Π_m(..), w_m(..), R_m(..), yd_m(..), vp_m(..), k(..), d(..), φ(..), μ_I(..), μ_A(..), mshock(..), g(..), c(..), x(..), Π(..), w(..), R(..), yd(..), vp(..), Πstar(..), u(..), q(..), λ(..), μ_z(..), r(..), l(..), mc(..), g1(..), g2(..), ob_μz(..), ob_Π(..), ob_R(..), ob_dw(..), ob_dy(..), ob_l(..), ob_μ(..)

    # Precompute some parameters, and steady state values
    ΛYd = (ΛA + α * Λμ) / (1 - α)
    Λx = exp(ΛYd)
    μ_z_val = exp(ΛYd)
    μ_I_val = exp(Λμ)
    γ1 = μ_z_val * μ_I_val / β - (1 - δ)
    Rbar = Πbar * μ_z_val / β
    Πstar_val = ((1 - θp * Πbar^(-(1 - ε) * (1 - χ))) / (1 - θp))^(1 / (1 - ε))
    vp_val = (1 - θp) / (1 - θp * Πbar^((1 - χ) * ε)) * Πstar_val^(-ε)
    mc_val = (1 - β * θp * Πbar^((1 - χ) * ε)) / (1 - β * θp * Πbar^(-(1 - ε) * (1 - χ))) * (ε - 1) / ε * Πstar_val
    w_val = (1 - α) * (mc_val * (α / γ1)^α)^(1 / (1 - α))
    k_val = α / (1 - α) * w_val / γ1 * μ_z_val * μ_I_val
    yd_val = (exp(ΛA) / μ_z_val * k_val^α - ϕ) / vp_val
    x_val = (μ_z_val * μ_I_val - (1 - δ)) / (μ_z_val * μ_I_val) * k_val
    c_val = (1 - g_bar) * yd_val - x_val
    λ_val = (1 - h * β / μ_z_val) / (1 - h / μ_z_val) / c_val
    ψ = w_val * λ_val
    g1_val = λ_val * yd_val / (1 - β * θp * Πbar^(-(1 - ε) * (1 - χ))) * (ε - 1) / ε * Πstar_val
    g2_val = λ_val * yd_val / (1 - β * θp * Πbar^(-(1 - ε) * (1 - χ))) * Πstar_val

    x_sym = [c_m, x_m, Π_m, w_m, R_m, yd_m, vp_m, k, d, φ, μ_I, μ_A, mshock, g]
    y_sym = [c, x, Π, w, R, yd, vp, Πstar, u, q, λ, μ_z, r, l, mc, g1, g2, ob_μz, ob_Π, ob_R, ob_dw, ob_dy, ob_l, ob_μ]
    # Merge p and p_f
    p = [β, h, ϑ, κ, α, θp, χ, γR, γy, γΠ, Πbar, ρd, ρφ, ρg, g_bar, σ_A, σ_d, σ_φ, σ_μ, σ_m, σ_g, Λμ, ΛA, δ, ε, ϕ, γ2, Ω_ii]

    H = [
         # Household Problem
         d(t) * (c(t) - h * c_m(t) * μ_z(t)^(-1))^(-1) - h * β * d(t+1) * (c(t+1) * μ_z(t+1) - h * c(t))^(-1) - λ(t),
         ψ * φ(t) * l(t)^ϑ - w(t) * λ(t),
         λ(t) - β * λ(t+1) * μ_z(t+1)^(-1) * R(t) / Π(t+1),
         γ1 + γ2 * (u(t) - 1) - r(t),
         q(t) - β * λ(t+1) / λ(t) / μ_z(t+1) / μ_I(t+1) * ((1 - δ) * q(t+1) + r(t+1) * u(t+1) - (γ1 * (u(t+1) - 1) + γ2 / 2 * (u(t+1) - 1)^2)),
         q(t) * (1 - (κ / 2 * (x(t) / x_m(t) * μ_z(t) - Λx)^2) - (κ * (x(t) / x_m(t) * μ_z(t) - Λx) * x(t) / x_m(t) * μ_z(t))) + β * q(t+1) * λ(t+1) / λ(t) / μ_z(t+1) * κ * (x(t+1) / x(t) * μ_z(t+1) - Λx) * (x(t+1) / x(t) * μ_z(t+1))^2 - 1,

         # Firm Problem
         λ(t) * mc(t) * yd(t) + β * θp * (Π(t)^χ / Π(t+1))^(-ε) * g1(t+1) - g1(t),
         λ(t) * Πstar(t) * yd(t) + β * θp * (Π(t)^χ / Π(t+1))^(1 - ε) * Πstar(t) / Πstar(t+1) * g2(t+1) - g2(t),
         ε * g1(t) - (ε - 1) * g2(t),
         u(t) * k(t) / l(t) - α / (1 - α) * w(t) / r(t) * μ_z(t) * μ_I(t),
         (1 / (1 - α))^(1 - α) * (1 / α)^α * w(t)^(1 - α) * r(t)^α - mc(t),

         # Price Evolution
         θp * (Π_m(t)^χ / Π(t))^(1 - ε) + (1 - θp) * Πstar(t)^(1 - ε) - 1,

         # Taylor Rule
         R(t) / Rbar - (R_m(t) / Rbar)^γR * ((Π(t) / Πbar)^γΠ * ((yd(t) / yd_m(t) * μ_z(t)) / exp(ΛYd))^γy)^(1 - γR) * mshock(t),

         # Market Clear
         c(t) + x(t) + g(t) + μ_z(t)^(-1) * μ_I(t)^(-1) * (γ1 * (u(t) - 1) + γ2 / 2 * (u(t) - 1)^2) * k(t) - yd(t),
         (μ_A(t) / μ_z(t) * (u(t) * k(t))^α * l(t)^(1 - α) - ϕ) / vp(t) - yd(t),
         θp * (Π_m(t)^χ / Π(t))^(-ε) * vp_m(t) + (1 - θp) * Πstar(t)^(-ε) - vp(t),
         k(t+1) * μ_z(t) * μ_I(t) - (1 - δ) * k(t) - μ_z(t) * μ_I(t) * (1 - κ / 2 * (x(t) / x_m(t) * μ_z(t) - Λx)^2) * x(t),

         # Variable Mapping
         c(t) - c_m(t+1), x(t) - x_m(t+1), Π(t) - Π_m(t+1), w(t) - w_m(t+1), R(t) - R_m(t+1), yd(t) - yd_m(t+1), vp(t) - vp_m(t+1),

         # Shock Evolution
         μ_z(t) - μ_A(t)^(1 / (1 - α)) * μ_I(t)^(α / (1 - α)),
         log(d(t+1)) - ρd * log(d(t)),
         log(φ(t+1)) - ρφ * log(φ(t)),
         log(μ_I(t+1)) - Λμ,
         log(μ_A(t+1)) - ΛA,
         log(mshock(t+1)) - 0,
         log(g(t+1)) - ρg * log(g(t)) - (1 - ρg) * log(g_bar * yd_val),

         # Observation equations Π, R, dw, dy, l, μ_I(-1)
         ob_μz(t) - log(μ_z(t)), ob_Π(t) - log(Π(t)), ob_R(t) - log(R(t)), ob_dw(t) - (log(w(t)) - log(w_m(t))), ob_dy(t) - (log(yd(t)) - log(yd_m(t))), ob_l(t) - log(l(t)), ob_μ(t) - log(μ_I(t))
    ]

    steady_states = [c_m(∞) ~ c_val, x_m(∞) ~ x_val, Π_m(∞) ~ Πbar, w_m(∞) ~ w_val, R_m(∞) ~ Rbar, yd_m(∞) ~ yd_val, vp_m(∞) ~ vp_val, k(∞) ~ k_val, d(∞) ~ 1.0, φ(∞) ~ 1.0, μ_I(∞) ~ exp(Λμ), μ_A(∞) ~ exp(ΛA), mshock(∞) ~ 1.0, g(∞) ~ g_bar * yd_val, c(∞) ~ c_val, x(∞) ~ x_val, Π(∞) ~ Πbar, w(∞) ~ w_val, R(∞) ~ Rbar, yd(∞) ~ yd_val, vp(∞) ~ vp_val, Πstar(∞) ~ Πstar_val, u(∞) ~ 1.0, q(∞) ~ 1.0, λ(∞) ~ λ_val, μ_z(∞) ~ exp(ΛYd), r(∞) ~ γ1, l(∞) ~ 1.0, mc(∞) ~ mc_val, g1(∞) ~ g1_val, g2(∞) ~ g2_val, ob_μz(∞) ~ ΛYd, ob_Π(∞) ~ log(Πbar),ob_R(∞) ~ log(Rbar), ob_dw(∞) ~ 0.0, ob_dy(∞) ~ 0.0, ob_l(∞) ~ 0.0, ob_μ(∞) ~ Λμ]

    n_ϵ = 6
    n_x = length(x_sym)
    n_y = length(y_sym)
    n_z = 6
    Γ = zeros(Num, n_ϵ, n_ϵ) # Make sure it is not a float64 matrix
    Γ[1, 1] = σ_d
    Γ[2, 2] = σ_φ
    Γ[3, 3] = σ_μ
    Γ[4, 4] = σ_A
    Γ[5, 5] = σ_m
    Γ[6, 6] = σ_g
    η_mat = zeros(n_x, n_ϵ)
    # The 9-14 states are: d φ μ_I μ_A mshock g
    η_mat[9, 1] = 1
    η_mat[10, 2] = 1
    η_mat[11, 3] = 1
    η_mat[12, 4] = 1
    η_mat[13, 5] = 1
    η_mat[14, 6] = 1
    # Define Q
    Q = zeros(n_z, n_x + n_y)
    Q[:, (n_y - 5):n_y] = [1.0 0.0 0.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0 0.0 0.0; 0.0 0.0 1.0 0.0 0.0 0.0; 0.0 0.0 0.0 1.0 0.0 0.0; 0.0 0.0 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 0.0 0.0 1.0]
    # Define Ω
    Ω = fill(Ω_ii, n_z)

	return H, (; t, x = x_sym, y = y_sym, p, steady_states, Γ, η = η_mat, Q, Ω), "FVGQ20"

end

p_d = (β=0.998, h=0.97, ϑ=1.17, κ=9.51, α=0.21, θp=0.82, χ=0.63, γR=0.77, γy=0.19, γΠ=1.29, Πbar=1.01, ρd=0.12, ρφ=0.93, ρg=0.95, g_bar=0.3, σ_A=exp(-3.97), σ_d=exp(-1.51), σ_φ=exp(-2.36), σ_μ=exp(-5.43), σ_m=exp(-5.85), σ_g=exp(-3.0), Λμ=3.4e-3, ΛA=2.8e-3)
p_f = (δ=0.025, ε=10, ϕ=0, γ2=0.001, Ω_ii=sqrt(1e-5))
