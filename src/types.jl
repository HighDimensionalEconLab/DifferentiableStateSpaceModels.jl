
# Model Types.  The template args are required for inference for cache/perturbation solutions
struct PerturbationModel{MaxOrder,N_y,N_x,N_ϵ,N_z,N_p,HasΩ,T1,T2}
    mod::Module

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
    return PerturbationModel{mod.max_order,mod.n_y,mod.n_x,mod.n_ϵ,mod.n_z,mod.n_p,
                             mod.has_Ω,typeof(mod.η),typeof(mod.Q)}(mod, mod.max_order,
                                                                    mod.n_y, mod.n_x,
                                                                    mod.n_p, mod.n_ϵ,
                                                                    mod.n_z, mod.has_Ω,
                                                                    mod.η, mod.Q)
end

# TODO: Add in latex stuff for the mod.H_latex, 
function Base.show(io::IO, ::MIME"text/plain", m::PerturbationModel) where {T}
    return print(io,
                 "Perturbation Model: n_y = $(m.n_y), n_x = $(m.n_x), n_p = $(m.n_p), n_ϵ = $(m.n_ϵ), n_z = $(m.n_z)\n y = $(m.mod.y_symbols) \n x = $(m.mod.x_symbols) \n p = $(m.mod.p_symbols)")
end

# The cache if for both 1st and 2nd order
# Constructors set values to nothing as appropriate
Base.@kwdef mutable struct SolverCache{Order,MatrixType,MatrixType2,MatrixType3,MatrixType4,
                                       MatrixType5,VectorType,VectorType2,
                                       VectorOfVectorType,VectorOfMatrixType,
                                       VectorOfMatrixType2,VectorOfMatrixType3,
                                       VectorOrNothingType,VectorOfVectorOrNothingType,
                                       MatrixScalingOrNothingType,SymmetricMatrixType,
                                       SymmetricVectorOfMatrixType,
                                       VectorOfVectorOfMatrixType,ThreeTensorType,
                                       CholeskyType,ChangeVarianceType,
                                       VectorOfThreeTensorType}
    order::Val{Order}  # allows inference in construction
    p_d_symbols::Vector{Symbol}
    H::VectorType
    H_yp::MatrixType
    H_y::MatrixType
    H_xp::MatrixType
    H_x::MatrixType
    H_yp_p::VectorOfMatrixType
    H_y_p::VectorOfMatrixType
    H_xp_p::VectorOfMatrixType
    H_x_p::VectorOfMatrixType
    H_p::VectorOfVectorType
    Γ::MatrixType2
    Γ_p::VectorOfMatrixType2
    Σ::SymmetricMatrixType
    Σ_p::SymmetricVectorOfMatrixType
    Ω::VectorOrNothingType
    Ω_p::VectorOfVectorOrNothingType
    Ψ::VectorOfMatrixType

    # Used in solution
    x::VectorType
    y::VectorType
    y_p::VectorOfVectorType
    x_p::VectorOfVectorType
    g_x::MatrixType4
    h_x::MatrixType4
    g_x_p::VectorOfMatrixType3
    h_x_p::VectorOfMatrixType3
    B::MatrixType2
    B_p::VectorOfMatrixType2
    Q::MatrixScalingOrNothingType
    η::MatrixType3
    A_1_p::VectorOfMatrixType3
    C_1::MatrixType4
    C_1_p::VectorOfMatrixType3
    V::CholeskyType
    V_p::ChangeVarianceType

    # Additional for 2nd order
    Ψ_p::VectorOfVectorOfMatrixType
    Ψ_yp::VectorOfVectorOfMatrixType
    Ψ_y::VectorOfVectorOfMatrixType
    Ψ_xp::VectorOfVectorOfMatrixType
    Ψ_x::VectorOfVectorOfMatrixType
    g_xx::ThreeTensorType
    h_xx::ThreeTensorType
    g_σσ::VectorType2
    h_σσ::VectorType2
    g_xx_p::VectorOfThreeTensorType
    h_xx_p::VectorOfThreeTensorType
    g_σσ_p::MatrixType5
    h_σσ_p::MatrixType5

    # Additional for solution type 2nd order
    A_0_p::MatrixType5
    A_2_p::VectorOfThreeTensorType
    C_0::VectorType2
    C_2::ThreeTensorType
    C_0_p::MatrixType5
    C_2_p::VectorOfThreeTensorType
end

# The Val(2), etc. for the order required for inference to function
# Note that the n_p_d is the number of differentiated parameters to allocate for

function SolverCache(m::PerturbationModel{MaxOrder,N_y,N_x,N_ϵ,N_z,N_p,HasΩ,T1,T2},
                     ::Val{Order},
                     p_d_symbols) where {Order,MaxOrder,N_y,N_x,N_ϵ,N_z,N_p,HasΩ,T1,T2}
    n_p_d = length(p_d_symbols)
    return SolverCache(; order=Val(Order), p_d_symbols, H=zeros(N_x + N_y),
                       H_yp=zeros(N_x + N_y, N_y), H_y=zeros(N_x + N_y, N_y),
                       H_xp=zeros(N_x + N_y, N_x), H_x=zeros(N_x + N_y, N_x),
                       Γ=zeros(N_ϵ, N_ϵ), Ω=!HasΩ ? nothing : zeros(N_z),
                       Ψ=[zeros(2(N_x + N_y), 2(N_x + N_y)) for i in 1:(N_x + N_y)],
                       H_p=[zeros(N_x + N_y) for i in 1:n_p_d],
                       H_yp_p=[zeros(N_x + N_y, N_y) for i in 1:n_p_d],
                       H_y_p=[zeros(N_x + N_y, N_y) for i in 1:n_p_d],
                       H_xp_p=[zeros(N_x + N_y, N_x) for i in 1:n_p_d],
                       H_x_p=[zeros(N_x + N_y, N_x) for i in 1:n_p_d],
                       Γ_p=[zeros(N_ϵ, N_ϵ) for i in 1:n_p_d],
                       Ω_p=!HasΩ ? nothing : [zeros(N_z) for i in 1:n_p_d], x=zeros(N_x),
                       y=zeros(N_y), y_p=[zeros(N_y) for i in 1:n_p_d],
                       x_p=[zeros(N_x) for i in 1:n_p_d], g_x=zeros(N_y, N_x),
                       h_x=zeros(N_x, N_x), g_x_p=[zeros(N_y, N_x) for _ in 1:n_p_d],
                       h_x_p=[zeros(N_x, N_x) for _ in 1:n_p_d],
                       Σ=Symmetric(zeros(N_ϵ, N_ϵ)),
                       Σ_p=[Symmetric(zeros(N_ϵ, N_ϵ)) for _ in 1:n_p_d], m.Q, m.η,
                       B=zeros(N_x, N_ϵ), B_p=[zeros(N_x, N_ϵ) for _ in 1:n_p_d],
                       C_1=zeros(N_z, N_x), C_1_p=[zeros(N_z, N_x) for _ in 1:n_p_d],
                       A_1_p=[zeros(N_x, N_x) for _ in 1:n_p_d], V=cholesky(Array(I(N_x))),
                       V_p=[zeros(N_x, N_x) for _ in 1:n_p_d],

                       # Stuff for 2nd order
                       Ψ_p=(Order == 1) ? nothing :
                           [[zeros(2 * (N_x + N_y), 2 * (N_x + N_y)) for _ in 1:(N_x + N_y)]
                            for _ in 1:n_p_d],
                       Ψ_yp=(Order == 1) ? nothing :
                            [[zeros(2 * (N_x + N_y), 2 * (N_x + N_y))
                              for _ in 1:(N_x + N_y)] for _ in 1:N_y],
                       Ψ_y=(Order == 1) ? nothing :
                           [[zeros(2 * (N_x + N_y), 2 * (N_x + N_y)) for _ in 1:(N_x + N_y)]
                            for _ in 1:N_y],
                       Ψ_xp=(Order == 1) ? nothing :
                            [[zeros(2 * (N_x + N_y), 2 * (N_x + N_y))
                              for _ in 1:(N_x + N_y)] for _ in 1:N_x],
                       Ψ_x=(Order == 1) ? nothing :
                           [[zeros(2 * (N_x + N_y), 2 * (N_x + N_y)) for _ in 1:(N_x + N_y)]
                            for _ in 1:N_x],
                       g_xx=(Order == 1) ? nothing : zeros(N_y, N_x, N_x),
                       h_xx=(Order == 1) ? nothing : zeros(N_x, N_x, N_x),
                       g_σσ=(Order == 1) ? nothing : zeros(N_y),
                       h_σσ=(Order == 1) ? nothing : zeros(N_x),
                       g_xx_p=(Order == 1) ? nothing :
                              [zeros(N_y, N_x, N_x) for _ in 1:n_p_d],
                       h_xx_p=(Order == 1) ? nothing :
                              [zeros(N_x, N_x, N_x) for _ in 1:n_p_d],
                       g_σσ_p=(Order == 1) ? nothing : zeros(N_y, n_p_d),
                       h_σσ_p=(Order == 1) ? nothing : zeros(N_x, n_p_d),
                       A_0_p=(Order == 1) ? nothing : zeros(N_x, n_p_d),
                       A_2_p=(Order == 1) ? nothing :
                             [zeros(N_x, N_x, N_x) for _ in 1:n_p_d],
                       C_0=(Order == 1) ? nothing : zeros(N_z),
                       C_0_p=(Order == 1) ? nothing : zeros(N_z, n_p_d),
                       C_2=(Order == 1) ? nothing : zeros(N_z, N_x, N_x),
                       C_2_p=(Order == 1) ? nothing :
                             [zeros(N_z, N_x, N_x) for _ in 1:n_p_d])
end
Base.@kwdef struct PerturbationSolverSettings{T1,T2,T3,T4,T5,T6}
    print_level::Int64 = 1  # 0 is no output at all
    ϵ_BK::Float64 = 1e-6 # For checking Blanchard-Kahn condition
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
end

function nlsolve_options(s::PerturbationSolverSettings)
    return (method=s.nlsolve_method, iterations=s.nlsolve_iterations,
            show_trace=s.nlsolve_show_trace, ftol=s.nlsolve_ftol)
end

# For callbacks
struct PerturbationSolver{T1,T2,T3}
    model::T1
    cache::T2
    settings::T3
end

# State Space types
abstract type AbstractPerturbationSolution end
abstract type AbstractFirstOrderPerturbationSolution <: AbstractPerturbationSolution end

# For this, all are dense due to schur decomposition
# All are dense due to schur decomposition
Base.@kwdef struct FirstOrderPerturbationSolution{T1<:AbstractVector,T2<:AbstractVector,
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
    D::T6  # current a matrix or nothing, later could make more general
    Q::T7  # can be nothing
    η::T8
    C::T9  # i.e. Q * g_x
    x_ergodic::T10
    Γ::T11
end

maybe_diagonal(x::AbstractVector) = TuringDiagMvNormal(zero(x), x)
maybe_diagonal(x) = x # otherwise, just return raw.  e.g. nothing

function FirstOrderPerturbationSolution(retcode, m::PerturbationModel, c::SolverCache)
    return FirstOrderPerturbationSolution(; retcode, m.mod.x_symbols, m.mod.y_symbols,
                                          m.mod.u_symbols, m.mod.p_symbols, c.p_d_symbols,
                                          m.n_x, m.n_y, m.n_p, m.n_ϵ, m.n_z, c.Q, c.η, c.y,
                                          c.x, c.B, D=maybe_diagonal(c.Ω), c.g_x, A=c.h_x,
                                          C=c.C_1,
                                          x_ergodic=TuringDenseMvNormal(zeros(m.n_x), c.V),
                                          c.Γ)
end

Base.@kwdef struct SecondOrderPerturbationSolution{T1<:AbstractVector,T2<:AbstractVector,
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
                   AbstractPerturbationSolution
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

    C_1::T14
    C_0::T15
    C_2::T16
end

function SecondOrderPerturbationSolution(retcode, m::PerturbationModel, c::SolverCache)
    return SecondOrderPerturbationSolution(; retcode, m.mod.x_symbols, m.mod.y_symbols,
                                           m.mod.u_symbols, m.mod.p_symbols, c.p_d_symbols,
                                           m.n_x, m.n_y, m.n_p, m.n_ϵ, m.n_z, c.Q, c.η, c.y,
                                           c.x, c.B, D=maybe_diagonal(c.Ω), c.Γ, c.g_x,
                                           A_1=c.h_x, c.g_xx, A_2=0.5 * c.h_xx, c.g_σσ,
                                           A_0=0.5 * c.h_σσ, c.C_1, c.C_0, c.C_2)
end
