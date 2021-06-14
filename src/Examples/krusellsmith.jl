#krusellsmith.jl
# MTK representation of Krusell Smith (1998) model as in Childers & Dogra (2019)
using SpecialFunctions

function krusellsmith(order = 1, functions_type = DenseFunctions(), simplify = true,
                      parallel = ModelingToolkit.SerialForm())

    nw = 35
    ns = 8

    @parameters β γ α δ sigs ρz smoother slo shi zlo zhi wlo whi
    @variables ℓ[1:nw] m[1:nw] k z Lℓ[1:nw] Lm[1:nw]
    @make_markov ℓ m k z Lℓ Lm

    x = cat(Lm, Lℓ, k, z, dims = 1)
    y = cat(m, ℓ, dims = 1)
    p = [β, γ, α, δ, sigs, ρz] #Differentiated parameters
    p_f = [smoother, zlo, zhi] #non-differentiated parameters

    #Setup
    # Cross-section distribution ranges
    slo = 0.5
    shi = 1.5
    wlo = 0.5
    whi = 15.0

    # make grids
    wgrid = LinRange(wlo, whi, nw)           #Evenly spaced grid
    wweights = ((whi - wlo) / nw) * ones(nw)          #quadrature weights

    sgrid = LinRange(slo, shi, ns)   # sgrid is grid for "skill"
    sweights = ((shi - slo) / ns) * ones(ns)  #quadrature weights

    #Generate densities manually because Distributions.jl doesn't allow MTK objects
    lognpdf(x, mu, sig) = exp.(-(log.(x) .- mu) .^ 2 / 2 * sig^2) ./ (x * sig * sqrt(2 * π))
    @register lognpdf(x, mu, sig)
    logncdf(x, mu, sig) = 0.5 .+ 0.5 * erf((log.(x) .- mu) ./ sqrt(2) * sig)  #calls erf from SpecialFunctions.jl
    @register logncdf(x, mu, sig)
    function trunc_lognpdf(x, UP, LO, mu, sig)
        return (x .- LO .>= -eps()) .* (x .- UP .<= eps()) .* lognpdf.(x, mu, sig) ./
               (logncdf(UP, mu, sig) - logncdf(LO, mu, sig))
    end
    @register trunc_lognpdf(x, UP, LO, mu, sig)

    g = trunc_lognpdf.(sgrid, shi, slo, 0.0, sigs) # density g evaluated on sgrid
    L = sum(sweights .* g .* sgrid)                   # Aggregate Labor endowment
    #L = simplify(L1) #This gets used a lot, so simplifying now is helpful

    MPK(Z, K) = α * Z * (K^(α - 1.0)) * (L^(1.0 - α)) + 1.0 - δ   # Marginal product of capital + 1-δ
    MPL(Z, K) = (1.0 - α) * Z * (K^α) * (L^(-α))  # Marginal product of labor

    function mollifier(x; zhi = zhi, zlo = zlo, smoother = smoother)
        In = 0.443993816237631
        temp = (-1.0 .+ 2.0 * (x .- zlo) ./ (zhi - zlo)) ./ smoother
        res = ((zhi .- zlo) / 2.0) * exp.(min.(-1.0 ./ (1.0 .- temp .^ 2), 2.0))
        res = (res ./ (In * smoother))
        res = res .* (x .> zlo) .* (x .< zhi)
        return res
    end
    # TODO: This changed at some point after our implementation.
    @register mollifier(x; zhi = zhi, zlo = zlo, smoother = smoother)

    ## Derivative of mollifier fucntion: not needed with automatic derivatives
    # function dmollifier(x, zhi, zlo, smoother)
    #     temp = (-1.0 .+ 2.0 * (x .- zlo) ./ (zhi .- zlo)) ./ smoother
    #     dy  = -(2 * temp ./ ((1 .- temp.^2).^2)) .* (2 ./ (smoother*(zhi .- zlo))) .* mollifier(x, zhi=zhi, zlo=zlo, smoother=smoother)
    # end

    cfunc(ell) = min.(ell .^ (-1.0 / γ), wgrid)
    function mollificand(ell, zp, kp)
        return (repeat(wgrid', nw, ns) .-
                MPK(zp, kp) * repeat(wgrid - cfunc(ell), 1, nw * ns)) / MPL(zp, kp) .-
               kron(sgrid', ones(nw, nw))
    end
    function KFmollificand(ell, zp, kp)
        return (repeat(wgrid, 1, nw * ns) -
                MPK(zp, kp) * repeat(wgrid' - cfunc(ell)', nw, ns)) / MPL(zp, kp) .-
               kron(sgrid', ones(nw, nw))
    end

    #Model equations
    H = [
        ℓ .-
        β *
        (MPK(z_p, k_p) / MPL(z_p, k_p)) *
        mollifier(mollificand(ℓ, z_p, k_p)) *
        kron(sweights .* g, wweights .* (cfunc(ℓ_p) .^ (-γ))) #Euler
        m .-
        mollifier(KFmollificand(Lℓ, z, k)) *
        kron(sweights .* (g / MPL(z, k)), wweights .* Lm) #Kolmogorov Forward
        k_p .- sum((wgrid .- cfunc(ℓ)) .* (m .* wweights)) #market clearing
        z_p .- ρz * z #Aggregate shock
        Lm_p .- m #Auxiliary lagged variable
        Lℓ_p .- ℓ #Auxiliary lagged variable
    ]

    n_ϵ = 1
    n_x = length(x)
    n_y = length(y)
    n_p = length(p)
    n_p = length(p_f)
    Γ = reshape([sigs], n_ϵ, n_ϵ)
    η = reshape([zeros(2 * nw + 1); -1], n_x, n_ϵ) # η is n_x * n_ϵ matrix

    return H, (; x, y, Γ, η, p_f, p, functions_type, simplify, parallel, order)
end
