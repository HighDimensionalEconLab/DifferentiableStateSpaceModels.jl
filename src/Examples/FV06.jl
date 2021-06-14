# Model as in https://www.sas.upenn.edu/~jesusfv/benchmark_DSGE.pdf

function FV06()
    # Parameter values
    # δ = 0.025
    # ε = 10
    # η = 10
    # ϕ = 0
    # γ2 = 0.001

    # Estimated parameters, taken from FV(2010), Table 3, p. 38, median estimate parameters
    # β = 0.998
    # h = 0.97
    # ψ = 8.92
    # ϑ = 1.17
    # κ = 9.51
    # α = 0.21

    # θp = 0.82
    # χ = 0.63
    # θw = 0.68
    # χw = 0.62
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

    # p_FV = [0.025, 10, 10, 0, 0.001]
    # p_FV = [0.998, 0.97, 8.92, 1.17, 9.51, 0.21, 0.82, 0.63, 0.68, 0.62, 0.77, 0.19, 1.29, 1.01, 0.12, 0.93, -3.97, -1.51, -2.36, -5.43, -5.85, 3.4e-3, 2.8e-3]

    # 8 variable relabelings
    # 36 equations

    # 28 parameters. 23 of them to estimate, 5 parameters as constants
    @parameters δ ε η ϕ γ2 β h ψ ϑ κ α θp χ θw χw γR γy γΠ Πbar ρd ρφ σ_A σ_d σ_φ σ_μ σ_m Λμ ΛA
    # States, 14
    @variables c_m x_m Π_m w_m R_m yd_m vp_m vw_m k d φ μ_I μ_A mshock
    # Controls, 22
    @variables c x Π w R yd vp vw wstar Πstarw Πstar u q f λ μ_z r l ld mc g1 g2
    # Markovize everything
    @make_markov c_m x_m Π_m w_m R_m yd_m vp_m vw_m mshock d φ μ_I μ_A k c x Π w R yd vp vw wstar Πstarw Πstar u q f λ μ_z r l ld mc g1 g2

    # Precompute some parameters
    ΛYd = (ΛA + α * Λμ) / (1 - α)
    Λx = exp(ΛYd)
    μ_z_value = exp(ΛYd)
    μ_I_value = exp(Λμ)
    γ1 = μ_z_value * μ_I_value / β - (1 - δ)
    Rbar = 1 + (Πbar * μ_z_value / β - 1)

    x_sym = [c_m, x_m, Π_m, w_m, R_m, yd_m, vp_m, vw_m, k, d, φ, μ_I, μ_A, mshock]
    y_sym = [c, x, Π, w, R, yd, vp, vw, wstar, Πstarw, Πstar, u, q, f, λ, μ_z, r, l, ld, mc,
             g1, g2]
    p = [β, h, ψ, ϑ, κ, α, θp, χ, θw, χw, γR, γy, γΠ, Πbar, ρd, ρφ, σ_A, σ_d, σ_φ, σ_μ, σ_m,
         Λμ, ΛA]
    p_f = [δ, ε, η, ϕ, γ2]

    H = [
         # Household Problem
         d * (c - h * c_m * μ_z^(-1))^(-1) - h * β * d_p * (c_p * μ_z_p - h * c)^(-1) - λ,
         λ - β * λ_p * μ_z_p^(-1) * R / Π_p, γ1 + γ2 * (u - 1) - r,
         q -
         β * λ_p / λ / μ_z_p / μ_I_p *
         ((1 - δ) * q_p + r_p * u_p - (γ1 * (u_p - 1) + γ2 / 2 * (u_p - 1)^2)),
         q *
         (1 - (κ / 2 * (x / x_m * μ_z - Λx)^2) - (κ * (x / x_m * μ_z - Λx) * x / x_m * μ_z)) +
         β * q_p * λ_p / λ / μ_z_p * κ * (x_p / x * μ_z_p - Λx) * (x_p / x * μ_z_p)^2 - 1,
         (η - 1) / η * wstar^(1 - η) * λ * w^η * ld +
         β * θw * (Π^χw / Π_p)^(1 - η) * (wstar_p / wstar * μ_z_p)^(η - 1) * f_p - f,
         ψ * d * φ * Πstarw^(-η * (1 + ϑ)) * ld^(1 + ϑ) +
         β *
         θw *
         (Π^χw / Π_p)^(-η * (1 + ϑ)) *
         (wstar_p / wstar * μ_z_p)^(η * (1 + ϑ)) *
         f_p - f,

         # Firm Problem
         λ * mc * yd + β * θp * (Π^χ / Π_p)^(-ε) * g1_p - g1,
         λ * Πstar * yd + β * θp * (Π^χ / Π_p)^(1 - ε) * Πstar / Πstar_p * g2_p - g2,
         ε * g1 - (ε - 1) * g2, u * k / ld - α / (1 - α) * w / r * μ_z * μ_I,
         (1 / (1 - α))^(1 - α) * (1 / α)^α * w^(1 - α) * r^α - mc,

         # Wage and Price Evolution
         θw * (Π_m^χw / Π)^(1 - η) * (w_m / w / μ_z)^(1 - η) + (1 - θw) * Πstarw^(1 - η) -
         1, θp * (Π_m^χ / Π)^(1 - ε) + (1 - θp) * Πstar^(1 - ε) - 1,

         # Taylor Rule
         R / Rbar -
         (R_m / Rbar)^γR *
         ((Π / Πbar)^γΠ * ((yd / yd_m * μ_z) / exp(ΛYd))^γy)^(1 - γR) *
         exp(mshock),

         # Market Clear
         c + x + μ_z^(-1) * μ_I^(-1) * (γ1 * (u - 1) + γ2 / 2 * (u - 1)^2) * k - yd,
         (μ_A / μ_z * (u * k)^α * ld^(1 - α) - ϕ) / vp - yd, l - vw * ld,
         θp * (Π_m^χ / Π)^(-ε) * vp_m + (1 - θp) * Πstar^(-ε) - vp,
         θw * (w_m / w * μ_z^(-1) * Π_m^χw / Π)^(-η) * vw_m + (1 - θw) * (Πstarw)^(-η) - vw,
         k_p * μ_z * μ_I - (1 - δ) * k -
         μ_z * μ_I * (1 - κ / 2 * (x / x_m * μ_z - Λx)^2) * x,

         # Definition of Π*_w
         Πstarw - wstar / w,

         # Variable Mapping
         c - c_m_p, x - x_m_p, Π - Π_m_p, w - w_m_p, R - R_m_p, yd - yd_m_p, vp - vp_m_p,
         vw - vw_m_p,

         # Shock Evolution
         μ_z - μ_A^(1 / (1 - α)) * μ_I^(α / (1 - α)), log(d_p) - ρd * log(d),
         log(φ_p) - ρφ * log(φ), log(μ_I_p) - Λμ, log(μ_A_p) - ΛA, mshock_p - 0]

    x̄_iv = [c_m_ss ~ 0.5, x_m_ss ~ 0.1, Π_m_ss ~ 1.0, w_m_ss ~ 1.01, R_m_ss ~ 1.01,
             yd_m_ss ~ 0.5, vp_m_ss ~ 1.0, vw_m_ss ~ 1.0, k_ss ~ 2.0, d_ss ~ 1.0,
             φ_ss ~ 1.0, μ_I_ss ~ exp(Λμ), μ_A_ss ~ exp(ΛA), mshock_ss ~ 0]
    ȳ_iv = [c_ss ~ 0.5, x_ss ~ 0.1, Π_ss ~ 1.0, w_ss ~ 1.01, R_ss ~ 1.01, yd_ss ~ 0.5,
             vp_ss ~ 1.0, vw_ss ~ 1.0, wstar_ss ~ 1.1, Πstarw_ss ~ 1.01, Πstar_ss ~ 1.01,
             u_ss ~ 1.0, q_ss ~ 1.0, f_ss ~ 2.5, λ_ss ~ 2.5, μ_z_ss ~ exp(ΛYd), r_ss ~ γ1,
             l_ss ~ 0.3, ld_ss ~ 0.3, mc_ss ~ 1.0, g1_ss ~ 7.0, g2_ss ~ 7.0]

    n_ϵ = 5
    n_x = length(x_sym)
    n_y = length(y_sym)
    n_p = length(p)
    n_p = length(p_f)
    Γ = zeros(Term, n_ϵ, n_ϵ) # Make sure it is not a float64 matrix
    Γ[1, 1] = exp(σ_d)
    Γ[2, 2] = exp(σ_φ)
    Γ[3, 3] = exp(σ_μ)
    Γ[4, 4] = exp(σ_A)
    Γ[5, 5] = exp(σ_m)
    η = zeros(n_x, n_ϵ)
    # The 10-14 states are: d φ μ_I μ_A mshock
    η[10, 1] = 1
    η[11, 2] = 1
    η[12, 3] = 1
    η[13, 4] = 1
    η[14, 5] = 1

    return H,
           (; x = x_sym, y = y_sym, x̄ = nothing, ȳ = nothing, Γ, η, p_f, p, x̄_iv, ȳ_iv),
           "FV06"

end
