#sguext.jl
#Code to implement extended version of Schmitt-Grohe and Uribe (2003, JIE) Model 2
#Shamelessly copied from Dynare .mod file of Cesa-Bianchi (2012)
#https://sites.google.com/site/ambropo/dynarecodes
#and from Dynare .mod file of Pfeifer (2019)
#https://github.com/JohannesPfeifer/DSGE_mod/tree/master/SGU_2003

function sguext()
    ∞ = Inf
    # Parameters
    @variables γ, ω, ρ, σe, δ, ψ, α, ϕ, β, r_w, d_bar, ρ_u, σu, ρ_v, σv
    # x and y
    @variables t::Integer, d(..), c(..), h(..), GDP(..), i(..), k(..), a(..), λ(..), tb(..),
               ca(..), riskpremium(..), r(..), Ld(..), Lr(..), Lk(..), ζ(..), μ(..)

    x = [a, ζ, μ, Ld, Lk, Lr]
    y = [d, c, h, GDP, i, k, λ, tb, ca, riskpremium, r]
    p = [γ, ω, ρ, σe, δ, ψ, α, ϕ, β, r_w, d_bar, ρ_u, σu, ρ_v, σv]

    # Model equations
    # Rewritten from Dynare to Schmitt-Grohe and Uribe timing conventions
    H = [d(t) - (1 + exp(Lr(t))) * Ld(t) + exp(GDP(t)) - exp(c(t)) - exp(i(t)) -
         (ϕ / 2) * (exp(k(t)) - exp(Lk(t)))^2, #Debt evolution
         exp(GDP(t)) - exp(a(t)) * (exp(Lk(t))^α) * (exp(h(t))^(1 - α)), #Production function
         exp(k(t)) - exp(i(t)) - (1 - δ) * exp(Lk(t)), #Capital evolution
         exp(λ(t)) - β * exp(μ(t)) * (1 + exp(r(t))) * exp(λ(t + 1)), #Euler equation
         (exp(c(t)) - ((exp(h(t))^ω) / ω))^(-γ) - exp(λ(t)), #Marginal utility
         ((exp(c(t)) - ((exp(h(t))^ω) / ω))^(-γ)) * (exp(h(t))^(ω - 1)) -
         exp(λ(t)) * (1 - α) * exp(GDP(t)) / exp(h(t)), #Labor FOC
         exp(λ(t)) * (1 + ϕ * (exp(k(t)) - exp(Lk(t)))) -
         β *
         exp(μ(t)) * exp(λ(t + 1)) *
         (α * exp(GDP(t + 1)) / exp(k(t)) + 1 - δ + ϕ * (exp(k(t + 1)) - exp(k(t)))), #Investment FOC
         a(t + 1) - ρ * a(t), #TFP evolution
         exp(r(t)) - r_w - riskpremium(t), #Interest rate
         riskpremium(t) - ψ * (exp(d(t) - d_bar) - 1) - ζ(t), #Risk premium
         tb(t) - 1 +
         ((exp(c(t)) + exp(i(t)) + (ϕ / 2) * (exp(k(t)) - exp(Lk(t)))^2) / exp(GDP(t))), #Trade Balance
         ca(t) - (1 / exp(GDP(t))) * (Ld(t) - d(t)), #Current Account
         Ld(t + 1) - d(t), #Auxiliary lagged variables
         Lk(t + 1) - k(t), Lr(t + 1) - r(t),
         ζ(t+1) - ρ_u * ζ(t), #Risk premium shock
         μ(t+1) - ρ_v * μ(t)] #Demand shock

    #Define representations for steady state

    hstar = ((1 - α) * (α / (r_w + δ))^(α / (1 - α)))^(1 / (ω - 1))
    kstar = hstar / (((r_w + δ) / α)^(1 / (1 - α)))
    istar = δ * kstar
    GDPstar = (kstar^α) * (hstar^(1 - α))
    cstar = GDPstar - istar - r_w * d_bar
    tbstar = GDPstar - cstar - istar
    λstar = (cstar - ((hstar^ω) / ω))^(-γ)

    #Steady state values

    steady_states = [a(∞) ~ 0, Ld(∞) ~ d_bar, Lk(∞) ~ log(kstar), Lr(∞) ~ log((1 - β) / β),
                     d(∞) ~ d_bar, c(∞) ~ log(cstar), h(∞) ~ log(hstar),
                     GDP(∞) ~ log(GDPstar), i(∞) ~ log(istar), k(∞) ~ log(kstar),
                     λ(∞) ~ log(λstar), tb(∞) ~ log(tbstar), ca(∞) ~ 0, riskpremium(∞) ~ 0,
                     r(∞) ~ log((1 - β) / β),
                     ζ(∞) ~ 0,
                     μ(∞) ~ 0]
    n_ϵ = 3
    Γ = zeros(Num, n_ϵ, n_ϵ) # make sure it is not a float64 matrix
    Γ[1, 1] = σe
    Γ[2, 2] = σu
    Γ[3, 3] = σv
    η = zeros(n_x, n_ϵ) # η is n_x * n_ϵ matrix
    η[1, 1] = -1 # e
    η[2, 2] = -1 # u
    η[3, 3] = -1 # v


    # Add some observation matrix for Q and Ω?
    return H, (; t, x, y, p, steady_states, Γ, η), "sguext"
end