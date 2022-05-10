import Base.deepcopy_internal
# A wrapper for Module.  The only purpose is a specialization of deepcopy since otherwise the "mod" property in the PerturbationModule brakes multithreaded MCMC
struct ModuleWrapper
    m::Module
end
function deepcopy_internal(x::ModuleWrapper, stackdict::IdDict)
    if haskey(stackdict, x)
        return stackdict[x]::ModuleWrapper
    end
    y = ModuleWrapper(x.m)
    stackdict[x] = y
    return y
end

# Model Types
# Model Types.  The template args are required for inference for cache/perturbation solutions
struct PerturbationModel{MaxOrder,HasΩ,T1,T2}
    mod::ModuleWrapper

    # Could extract from type, but here for simplicity
    max_order::Int64
    n_y::Int64
    n_x::Int64
    n_p::Int64
    n_ϵ::Int64
    n_z::Int64

    # Variations
    has_Ω::Bool
    η::T1
    Q::T2
end

# Construct from a module.  Inherently type unstable, so use function barrier from return type
function PerturbationModel(mod)
    return PerturbationModel{mod.max_order,mod.has_Ω,typeof(mod.η),typeof(mod.Q)}(ModuleWrapper(mod),
                                                                                  mod.max_order,
                                                                                  mod.n_y,
                                                                                  mod.n_x,
                                                                                  mod.n_p,
                                                                                  mod.n_ϵ,
                                                                                  mod.n_z,
                                                                                  mod.has_Ω,
                                                                                  mod.η,
                                                                                  mod.Q)
end

# TODO: Add in latex stuff for the mod.H_latex, 
function Base.show(io::IO, ::MIME"text/plain", m::PerturbationModel)
    return print(io,
                 "Perturbation Model: n_y = $(m.n_y), n_x = $(m.n_x), n_p = $(m.n_p), n_ϵ = $(m.n_ϵ), n_z = $(m.n_z)\n y = $(m.mod.m.y_symbols) \n x = $(m.mod.m.x_symbols) \n p = $(m.mod.m.p_symbols)")
end

# Buffers for the solvers to reduce allocations
# General rule for cache vs. buffers
# 1. If something should be used in multiple parts of the algorithm, put it in the cache
# 2. Otherwise, use the buffers, which you can "trash" with inplace operations as required
# 3. The cache should never be modified after it has been filled in a given sequence of events, buffers can be

struct FirstOrderSolverBuffers
    A::Matrix{Complex{Float64}}
    B::Matrix{Complex{Float64}}
    Z::Matrix{Float64}  # real version, transposed relative to the schur
    Z_ll::Matrix{Float64}
    S_bb::UpperTriangular{Float64,Matrix{Float64}}
    T_bb::UpperTriangular{Float64,Matrix{Float64}}
end
function FirstOrderSolverBuffers(n_y, n_x, n_p_d, n_ϵ, n_z)
    return FirstOrderSolverBuffers(zeros(Complex{Float64}, n_x + n_y, n_x + n_y),
                                   zeros(Complex{Float64}, n_x + n_y, n_x + n_y),
                                   zeros(n_x + n_y, n_x + n_y),
                                   zeros(n_y, n_y),
                                   UpperTriangular(zeros(n_x, n_x)),
                                   UpperTriangular(zeros(n_x, n_x)))
end
struct FirstOrderDerivativeSolverBuffers
    R::Matrix{Float64}
    A::Matrix{Float64}
    C::Matrix{Float64}
    D::Matrix{Float64}
    E::Matrix{Float64}
    dH::Matrix{Float64}
    bar::Matrix{Float64}
end
function FirstOrderDerivativeSolverBuffers(n_y, n_x, n_p_d, n_ϵ, n_z)
    return FirstOrderDerivativeSolverBuffers(zeros(2 * (n_x + n_y), n_x),
                                             zeros(n_x + n_y, n_x + n_y),
                                             zeros(n_x + n_y, n_x + n_y),
                                             zeros(n_x, n_x),
                                             zeros(n_x + n_y, n_x),
                                             zeros(n_x + n_y, 2 * (n_x + n_y)),
                                             zeros(2 * (n_x + n_y), 1))
end
struct SecondOrderSolverBuffers
    A::Matrix{Float64}
    B::Matrix{Float64}
    C::Matrix{Float64}
    D::Matrix{Float64}
    E::Matrix{Float64}
    R::Matrix{Float64}
    A_σ::Matrix{Float64}
    R_σ::Matrix{Float64}
end
function SecondOrderSolverBuffers(n_y, n_x, n_p_d, n_ϵ, n_z)
    return SecondOrderSolverBuffers(zeros(n_x + n_y, n_x + n_y),
                                    zeros(n_x^2, n_x^2), zeros(n_x + n_y, n_x + n_y),
                                    zeros(n_x^2, n_x^2),
                                    zeros(n_x + n_y, n_x^2), zeros(2 * (n_x + n_y), n_x),
                                    zeros(n_x + n_y, n_x + n_y),
                                    zeros(2 * (n_x + n_y), n_x))
end
struct SecondOrderDerivativeSolverBuffers
    A::Matrix{Float64}
    B::Matrix{Float64}
    C::Matrix{Float64}
    D::Matrix{Float64}
    E::Matrix{Float64}
    R::Matrix{Float64}
    dH::Matrix{Float64}
    dΨ::Vector{Matrix{Float64}}
    gh_stack::Matrix{Float64}
    g_xx_flat::Matrix{Float64}
    Ψ_x_sum::Vector{Vector{Matrix{Float64}}}
    Ψ_y_sum::Vector{Vector{Matrix{Float64}}}
    bar::Matrix{Float64}
    kron_h_x::Matrix{Float64}
    R_p::Matrix{Float64}
    A_σ::Matrix{Float64}
    R_σ::Matrix{Float64}
end

function SecondOrderDerivativeSolverBuffers(n_y, n_x, n_p_d, n_ϵ, n_z)
    return SecondOrderDerivativeSolverBuffers(zeros(n_x + n_y, n_x + n_y),
                                              zeros(n_x^2, n_x^2),
                                              zeros(n_x + n_y, n_x + n_y),
                                              zeros(n_x^2, n_x^2),
                                              zeros(n_x + n_y, n_x^2),
                                              zeros(2 * (n_x + n_y), n_x),
                                              zeros(n_x + n_y, 2 * (n_x + n_y)),
                                              [zeros(2 * (n_x + n_y), 2 * (n_x + n_y))
                                               for _ in 1:(n_x + n_y)],
                                              zeros(n_x + n_y, n_x^2),
                                              zeros(n_y, n_x^2),
                                              [[zeros(2 * (n_x + n_y), 2 * (n_x + n_y))
                                                for _ in 1:(n_x + n_y)]
                                               for _ in 1:n_x],
                                              [[zeros(2 * (n_x + n_y), 2 * (n_x + n_y))
                                                for _ in 1:(n_x + n_y)]
                                               for _ in 1:n_y],
                                              zeros(2 * (n_x + n_y), 1),
                                              zeros(n_x^2, n_x^2),
                                              zeros(2 * (n_x + n_y), n_x),
                                              zeros(n_x + n_y, n_x + n_y),
                                              zeros(2 * (n_x + n_y), n_x))
end

# The cache if for both 1st and 2nd order
# Constructors set values to nothing as appropriate
abstract type AbstractSolverCache{Order} end
Base.@kwdef struct SolverCache{Order,ΩType,Ω_pType,QType,ηType,g_σσType,g_xxType} <:
                   AbstractSolverCache{Order}
    order::Val{Order}  # allows inference in construction
    p_d_symbols::Vector{Symbol}
    H::Vector{Float64}
    H_yp::Matrix{Float64}
    H_y::Matrix{Float64}
    H_xp::Matrix{Float64}
    H_x::Matrix{Float64}
    H_yp_p::Vector{Matrix{Float64}}
    H_y_p::Vector{Matrix{Float64}}
    H_xp_p::Vector{Matrix{Float64}}
    H_x_p::Vector{Matrix{Float64}}
    H_p::Vector{Vector{Float64}}
    Γ::Matrix{Float64}
    Γ_p::Vector{Matrix{Float64}}
    Σ::Symmetric{Float64,Matrix{Float64}}
    Σ_p::Vector{Symmetric{Float64,Matrix{Float64}}}
    Ω::ΩType
    Ω_p::Ω_pType
    Ψ::Vector{Matrix{Float64}}

    # Used in solution
    x::Vector{Float64}
    y::Vector{Float64}
    y_p::Vector{Vector{Float64}}
    x_p::Vector{Vector{Float64}}
    g_x::Matrix{Float64}
    h_x::Matrix{Float64}
    g_x_p::Vector{Matrix{Float64}}
    h_x_p::Vector{Matrix{Float64}}
    B::Matrix{Float64}
    B_p::Vector{Matrix{Float64}}
    Q::QType
    η::ηType
    A_1_p::Vector{Matrix{Float64}}
    C_1::Matrix{Float64}
    C_1_p::Vector{Matrix{Float64}}
    V::PDMats.PDMat{Float64,Matrix{Float64}}
    V_p::Vector{Matrix{Float64}}
    η_Σ_sq::Symmetric{Float64,Matrix{Float64}}

    # Additional for 2nd order
    Ψ_p::Union{Nothing,Vector{Vector{Matrix{Float64}}}}
    Ψ_yp::Union{Nothing,Vector{Vector{Matrix{Float64}}}}
    Ψ_y::Union{Nothing,Vector{Vector{Matrix{Float64}}}}
    Ψ_xp::Union{Nothing,Vector{Vector{Matrix{Float64}}}}
    Ψ_x::Union{Nothing,Vector{Vector{Matrix{Float64}}}}
    g_xx::g_xxType
    h_xx::Union{Nothing,Array{Float64,3}}
    g_σσ::g_σσType
    h_σσ::Union{Nothing,Vector{Float64}}
    g_xx_p::Union{Nothing,Vector{Array{Float64,3}}}
    h_xx_p::Union{Nothing,Vector{Array{Float64,3}}}
    g_σσ_p::Union{Nothing,Matrix{Float64}}
    h_σσ_p::Union{Nothing,Matrix{Float64}}

    # Additional for solution type 2nd order
    A_0_p::Union{Nothing,Matrix{Float64}}
    A_2_p::Union{Nothing,Vector{Array{Float64,3}}}
    C_0::Union{Nothing,Vector{Float64}}
    C_2::Union{Nothing,Array{Float64,3}}
    C_0_p::Union{Nothing,Matrix{Float64}}
    C_2_p::Union{Nothing,Vector{Array{Float64,3}}}

    # Buffers for additional calculations
    first_order_solver_buffer::FirstOrderSolverBuffers
    first_order_solver_p_buffer::FirstOrderDerivativeSolverBuffers
    second_order_solver_buffer::Union{Nothing,SecondOrderSolverBuffers}
    second_order_solver_p_buffer::Union{Nothing,SecondOrderDerivativeSolverBuffers}
    I_x::Matrix{Float64}  #dense identity matrices
    I_x_2::Matrix{Float64}
    zeros_x_x::Matrix{Float64}
    zeros_y_x::Matrix{Float64}
end

function SolverCache(::Val{Order}, ::Val{HasΩ}, N_p_d, N_y, N_x, N_ϵ, N_z, Q,
                     η) where {Order,HasΩ}
    return SolverCache(; order = Val(Order), p_d_symbols = Vector{Symbol}(undef, N_p_d),
                       H = zeros(N_x + N_y),
                       H_yp = zeros(N_x + N_y, N_y), H_y = zeros(N_x + N_y, N_y),
                       H_xp = zeros(N_x + N_y, N_x), H_x = zeros(N_x + N_y, N_x),
                       Γ = zeros(N_ϵ, N_ϵ), Ω = !HasΩ ? nothing : zeros(N_z),
                       Ψ = [zeros(2(N_x + N_y), 2(N_x + N_y)) for i in 1:(N_x + N_y)],
                       H_p = [zeros(N_x + N_y) for i in 1:N_p_d],
                       H_yp_p = [zeros(N_x + N_y, N_y) for i in 1:N_p_d],
                       H_y_p = [zeros(N_x + N_y, N_y) for i in 1:N_p_d],
                       H_xp_p = [zeros(N_x + N_y, N_x) for i in 1:N_p_d],
                       H_x_p = [zeros(N_x + N_y, N_x) for i in 1:N_p_d],
                       Γ_p = [zeros(N_ϵ, N_ϵ) for i in 1:N_p_d],
                       Ω_p = !HasΩ ? nothing : [zeros(N_z) for i in 1:N_p_d],
                       x = zeros(N_x), y = zeros(N_y), y_p = [zeros(N_y) for i in 1:N_p_d],
                       x_p = [zeros(N_x) for i in 1:N_p_d], g_x = zeros(N_y, N_x),
                       h_x = zeros(N_x, N_x), g_x_p = [zeros(N_y, N_x) for _ in 1:N_p_d],
                       h_x_p = [zeros(N_x, N_x) for _ in 1:N_p_d],
                       Σ = Symmetric(zeros(N_ϵ, N_ϵ)), η_Σ_sq = Symmetric(zeros(N_x, N_x)),
                       Σ_p = [Symmetric(zeros(N_ϵ, N_ϵ)) for _ in 1:N_p_d], Q, η,
                       B = zeros(N_x, N_ϵ), B_p = [zeros(N_x, N_ϵ) for _ in 1:N_p_d],
                       C_1 = zeros(N_z, N_x), C_1_p = [zeros(N_z, N_x) for _ in 1:N_p_d],
                       A_1_p = [zeros(N_x, N_x) for _ in 1:N_p_d],
                       V = PDMat{Float64,Matrix{Float64}}(N_x,
                                                          Matrix{Float64}(undef, N_x, N_x),
                                                          Cholesky{Float64,Matrix{Float64}}(Matrix{Float64}(undef,
                                                                                                            N_x,
                                                                                                            N_x),
                                                                                            'U',
                                                                                            0)),
                       V_p = [zeros(N_x, N_x) for _ in 1:N_p_d],

                       # Stuff for 2nd order
                       Ψ_p = (Order == 1) ? nothing :
                             [[zeros(2 * (N_x + N_y), 2 * (N_x + N_y))
                               for _ in 1:(N_x + N_y)] for _ in 1:N_p_d],
                       Ψ_yp = (Order == 1) ? nothing :
                              [[zeros(2 * (N_x + N_y), 2 * (N_x + N_y))
                                for _ in 1:(N_x + N_y)] for _ in 1:N_y],
                       Ψ_y = (Order == 1) ? nothing :
                             [[zeros(2 * (N_x + N_y), 2 * (N_x + N_y))
                               for _ in 1:(N_x + N_y)] for _ in 1:N_y],
                       Ψ_xp = (Order == 1) ? nothing :
                              [[zeros(2 * (N_x + N_y), 2 * (N_x + N_y))
                                for _ in 1:(N_x + N_y)] for _ in 1:N_x],
                       Ψ_x = (Order == 1) ? nothing :
                             [[zeros(2 * (N_x + N_y), 2 * (N_x + N_y))
                               for _ in 1:(N_x + N_y)] for _ in 1:N_x],
                       g_xx = (Order == 1) ? nothing : zeros(N_y, N_x, N_x),
                       h_xx = (Order == 1) ? nothing : zeros(N_x, N_x, N_x),
                       g_σσ = (Order == 1) ? nothing : zeros(N_y),
                       h_σσ = (Order == 1) ? nothing : zeros(N_x),
                       g_xx_p = (Order == 1) ? nothing :
                                [zeros(N_y, N_x, N_x) for _ in 1:N_p_d],
                       h_xx_p = (Order == 1) ? nothing :
                                [zeros(N_x, N_x, N_x) for _ in 1:N_p_d],
                       g_σσ_p = (Order == 1) ? nothing : zeros(N_y, N_p_d),
                       h_σσ_p = (Order == 1) ? nothing : zeros(N_x, N_p_d),
                       A_0_p = (Order == 1) ? nothing : zeros(N_x, N_p_d),
                       A_2_p = (Order == 1) ? nothing :
                               [zeros(N_x, N_x, N_x) for _ in 1:N_p_d],
                       C_0 = (Order == 1) ? nothing : zeros(N_z),
                       C_0_p = (Order == 1) ? nothing : zeros(N_z, N_p_d),
                       C_2 = (Order == 1) ? nothing : zeros(N_z, N_x, N_x),
                       C_2_p = (Order == 1) ? nothing :
                               [zeros(N_z, N_x, N_x) for _ in 1:N_p_d],

                       # buffers for algorithms
                       first_order_solver_buffer = FirstOrderSolverBuffers(N_y, N_x, N_p_d,
                                                                           N_ϵ, N_z),
                       first_order_solver_p_buffer = FirstOrderDerivativeSolverBuffers(N_y,
                                                                                       N_x,
                                                                                       N_p_d,
                                                                                       N_ϵ,
                                                                                       N_z),
                       second_order_solver_buffer = (Order == 1) ? nothing :
                                                    SecondOrderSolverBuffers(N_y, N_x,
                                                                             N_p_d, N_ϵ,
                                                                             N_z),
                       second_order_solver_p_buffer = (Order == 1) ? nothing :
                                                      SecondOrderDerivativeSolverBuffers(N_y,
                                                                                         N_x,
                                                                                         N_p_d,
                                                                                         N_ϵ,
                                                                                         N_z),
                       I_x = Matrix{Float64}(I(N_x)), I_x_2 = Matrix{Float64}(I(N_x^2)),
                       zeros_x_x = zeros(N_x, N_x), zeros_y_x = zeros(N_y, N_x))
end

function SolverCache(m::PerturbationModel{MaxOrder,HasΩ,T1,T2}, ::Val{Order},
                     p_d) where {Order,MaxOrder,HasΩ,T1,T2}
    return SolverCache(Val(Order),
                       Val(HasΩ),
                       length(p_d),
                       m.n_y,
                       m.n_x,
                       m.n_ϵ,
                       m.n_z,
                       m.Q,
                       m.η)
end
Base.@kwdef struct PerturbationSolverSettings{T1,T2,T3,T4,T5,T6}
    rethrow_exceptions::Bool = false  # rethrows all exceptions to aid in debugging/etc.  Otherwise just uses 
    print_level::Int64 = 1  # 0 is no output at all
    ϵ_BK::Float64 = 1e-6 # For checking Blanchard-Kahn condition
    tol_cholesky::Float64 = 1e9 # for checking norm of covariance matrix, etc.
    check_posdef_cholesky::Bool = true
    perturb_covariance::Float64 = eps() # perturb the covariance matrix to ensure it is postive definite
    nlsolve_method::Symbol = :trust_region
    nlsolve_iterations::Int64 = 1000
    nlsolve_show_trace::Bool = false
    nlsolve_ftol::Float64 = 1e-8
    use_solution_cache::Bool = true
    evaluate_functions_callback::T1 = nothing
    calculate_steady_state_callback::T2 = nothing
    solve_first_order_callback::T3 = nothing
    solve_first_order_p_callback::T4 = nothing
    solve_second_order_callback::T5 = nothing
    solve_second_order_p_callback::T6 = nothing
    sylvester_solver::Symbol = :MatrixEquations
end

function nlsolve_options(s::PerturbationSolverSettings)
    return (method = s.nlsolve_method, iterations = s.nlsolve_iterations,
            show_trace = s.nlsolve_show_trace, ftol = s.nlsolve_ftol)
end

# State Space types
abstract type AbstractPerturbationSolution end
abstract type AbstractFirstOrderPerturbationSolution <: AbstractPerturbationSolution end
abstract type AbstractSecondOrderPerturbationSolution <: AbstractPerturbationSolution end

# For this, all are dense due to schur decomposition
# All are dense due to schur decomposition
struct FirstOrderPerturbationSolution{T1<:AbstractVector,T2<:AbstractVector,
                                      T3<:AbstractMatrix,T4<:AbstractMatrix,
                                      T5<:AbstractMatrix,
                                      T6<:Union{Nothing,Distribution},
                                      T7<:Union{Nothing,AbstractMatrix,
                                                UniformScaling},
                                      T8<:AbstractMatrix,T9<:AbstractMatrix,
                                      T10<:Distribution,T11<:AbstractMatrix} <:
       AbstractFirstOrderPerturbationSolution
    retcode::Symbol
    x_symbols::Vector{Symbol}
    y_symbols::Vector{Symbol}
    p_symbols::Vector{Symbol}
    p_d_symbols::Vector{Symbol}
    u_symbols::Vector{Symbol}
    # TODO: differentiated parameter ordering?
    n_y::Int64
    n_x::Int64
    n_p::Int64
    n_ϵ::Int64
    n_z::Int64
    y::T1
    x::T2
    g_x::T3
    A::T4
    B::T5
    C::T9  # i.e. Q * g_x
    D::T6  # current a matrix or nothing, later could make more general
    Q::T7  # can be nothing
    η::T8
    x_ergodic::T10
    Γ::T11
end

maybe_diagonal(x::AbstractVector) = MvNormal(Diagonal(abs2.(x)))
maybe_diagonal(x) = x # otherwise, just return raw.  e.g. nothing

function FirstOrderPerturbationSolution(retcode, m::PerturbationModel, c::SolverCache)
    return FirstOrderPerturbationSolution(retcode,
                                          m.mod.m.x_symbols,
                                          m.mod.m.y_symbols,
                                          m.mod.m.p_symbols,
                                          c.p_d_symbols,
                                          m.mod.m.u_symbols,
                                          m.n_y,
                                          m.n_x,
                                          m.n_p,
                                          m.n_ϵ,
                                          m.n_z,
                                          c.y,
                                          c.x,
                                          c.g_x,
                                          c.h_x,
                                          c.B,
                                          c.C_1,
                                          maybe_diagonal(c.Ω),
                                          c.Q,
                                          c.η,
                                          MvNormal(zeros(m.n_x), c.V), # already has cholesky taken
                                          c.Γ)
end

struct SecondOrderPerturbationSolution{T1<:AbstractVector,T2<:AbstractVector,
                                       T3<:AbstractMatrix,T4<:AbstractMatrix,
                                       T5<:AbstractMatrix,
                                       T6<:Union{Nothing,Distribution},
                                       T7<:Union{Nothing,AbstractMatrix,
                                                 UniformScaling},
                                       T8<:AbstractMatrix,T9<:AbstractMatrix,
                                       T10,T11<:AbstractArray,T12,
                                       T13<:AbstractVector,T14<:AbstractMatrix,
                                       T15<:AbstractVector,
                                       T16<:AbstractArray} <:
       AbstractSecondOrderPerturbationSolution
    retcode::Symbol
    x_symbols::Vector{Symbol}
    y_symbols::Vector{Symbol}
    p_symbols::Vector{Symbol}
    p_d_symbols::Vector{Symbol}
    u_symbols::Vector{Symbol}
    n_y::Int64
    n_x::Int64
    n_p::Int64
    n_ϵ::Int64
    n_z::Int64
    y::T1
    x::T2
    g_x::T3
    B::T5
    D::T6
    Q::T7  # can be nothing
    η::T8
    Γ::T9

    g_xx::T10
    g_σσ::T12
    A_0::T13
    A_1::T4
    A_2::T11

    C_0::T15
    C_1::T14
    C_2::T16
end

function SecondOrderPerturbationSolution(retcode, m::PerturbationModel, c::SolverCache)
    return SecondOrderPerturbationSolution(retcode,
                                           m.mod.m.x_symbols,
                                           m.mod.m.y_symbols,
                                           m.mod.m.p_symbols,
                                           c.p_d_symbols,
                                           m.mod.m.u_symbols,
                                           m.n_y,
                                           m.n_x,
                                           m.n_p,
                                           m.n_ϵ,
                                           m.n_z,
                                           c.y,
                                           c.x,
                                           c.g_x,
                                           c.B,
                                           maybe_diagonal(c.Ω),
                                           c.Q,
                                           c.η,
                                           c.Γ,
                                           c.g_xx,
                                           c.g_σσ,
                                           0.5 * c.h_σσ,
                                           c.h_x,
                                           0.5 * c.h_xx,
                                           c.C_0,
                                           c.C_1,
                                           c.C_2)
end
