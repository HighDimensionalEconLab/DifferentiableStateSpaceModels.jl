#sgusmallopen.jl
#Code to implement Schmitt-Grohe and Uribe (2003, JIE) Model 2
#Shamelessly copied from Dynare .mod file of Cesa-Bianchi (2012)
#https://sites.google.com/site/ambropo/dynarecodes
#and from Dynare .mod file of Pfeifer (2019)
#https://github.com/JohannesPfeifer/DSGE_mod/tree/master/SGU_2003

function sgusmallopen()

    @parameters γ ω ρ σe δ ψ α ϕ β r_w d_bar
    @variables d c h GDP i k a λ tb ca riskpremium r Ld Lr Lk
    @make_markov d c h GDP i k a λ tb ca riskpremium r Ld Lr Lk

    x = [a, Ld, Lk, Lr]
    y = [d, c, h, GDP, i, k, λ, tb, ca, riskpremium, r]
    p = [γ, ω, ρ, σe, δ, ψ, α, ϕ, β, r_w, d_bar] #Differentiated parameters
    p_f = [] #non-differentiated parameters

    #Model equations
    #Rewritten from Dynare to Schmitt-Grohe and Uribe timing conventions
    H = [
        d - (1 + exp(Lr)) * Ld + exp(GDP) - exp(c) - exp(i) -
        (ϕ / 2) * (exp(k) - exp(Lk))^2,            #Debt evolution
        exp(GDP) - exp(a) * (exp(Lk)^α) * (exp(h)^(1 - α)),                                                   #Production function
        exp(k) - exp(i) - (1 - δ) * exp(Lk),                                                                  #Capital evolution
        exp(λ) - β * (1 + exp(r)) * exp(λ_p),                                                                 #Euler equation
        (exp(c) - ((exp(h)^ω) / ω))^(-γ) - exp(λ),                                                      #Marginal utility
        ((exp(c) - ((exp(h)^ω) / ω))^(-γ)) * (exp(h)^(ω - 1)) -
        exp(λ) * (1 - α) * exp(GDP) / exp(h),                  #Labor FOC
        exp(λ) * (1 + ϕ * (exp(k) - exp(Lk))) -
        β * exp(λ_p) * (α * exp(GDP_p) / exp(k) + 1 - δ + ϕ * (exp(k_p) - exp(k))),  #Investment FOC
        a_p - ρ * a,                                                                                      #TFP evolution
        exp(r) - r_w - riskpremium,                                                                     #Interest rate
        riskpremium - ψ * (exp(d - d_bar) - 1),                                                           #Risk premium
        tb - 1 + ((exp(c) + exp(i) + (ϕ / 2) * (exp(k) - exp(Lk))^2) / exp(GDP)),                               #Trade Balance
        ca - (1 / exp(GDP)) * (Ld - d),                                                                       #Current Account
        Ld_p - d,                                                                                       #Auxiliary lagged variables
        Lk_p - k,
        Lr_p - r,
    ]

    #Define representations for steady state

    hstar = ((1 - α) * (α / (r_w + δ))^(α / (1 - α)))^(1 / (ω - 1))
    kstar = hstar / (((r_w + δ) / α)^(1 / (1 - α)))
    istar = δ * kstar
    GDPstar = (kstar^α) * (hstar^(1 - α))
    cstar = GDPstar - istar - r_w * d_bar
    tbstar = GDPstar - cstar - istar
    λstar = (cstar - ((hstar^ω) / ω))^(-γ)

    #Steady state values

    x̄ = [a_ss ~ 0, Ld_ss ~ d_bar, Lk_ss ~ log(kstar), Lr_ss ~ log((1 - β) / β)]
    ȳ = [
        d_ss ~ d_bar,
        c_ss ~ log(cstar),
        h_ss ~ log(hstar),
        GDP_ss ~ log(GDPstar),
        i_ss ~ log(istar),
        k_ss ~ log(kstar),
        λ_ss ~ log(λstar),
        tb_ss ~ log(tbstar),
        ca_ss ~ 0,
        riskpremium_ss ~ 0,
        r_ss ~ log((1 - β) / β),
    ]

    ## TO DO: Finish the rest
    x̄_iv = deepcopy(x̄)
    ȳ_iv = deepcopy(ȳ)

    n_ϵ = 1
    n_x = length(x)
    n_y = length(y)
    n_p = length(p)
    n_p = length(p_f)
    Γ = reshape([σe], n_ϵ, n_ϵ)
    η = reshape([-1, 0, 0, 0], n_x, n_ϵ) # η is n_x * n_ϵ matrix

    return H, (; x, y, x̄, ȳ, Γ, η, p_f, p, x̄_iv, ȳ_iv), "sgusmallopen"
end
