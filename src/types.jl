
# Model Types

struct PerturbationModel{MaxOrder, N_y, N_x, N_ϵ, N_z, N_p, HasΩ, T1, T2}
    mod::Module
    η::T1
    Q::T2

    # Could extract from type, but here for simplicity
    max_order::Int64
    n_y::Int64
    n_x::Int64
    n_p::Int64
    n_ϵ::Int64
    n_z::Int64
end


# Construct from a module.  Inherently type unstable, so use after a function barrier
function PerturbationModel(mod)
    return PerturbationModel{mod.max_order, mod.n_y, mod.n_x, mod.n_ϵ, mod.n_z, mod.n_p, isnothing(mod.Ω!), typeof(mod.η), typeof(mod.Q)}(mod, mod.η, mod.Q, mod.max_order, mod.n_y, mod.n_x, mod.n_p, mod.n_ϵ, mod.n_z)
end

# TODO: Add in latex stuff for the mod.H_latex, 
function Base.show(
    io::IO,
    ::MIME"text/plain",
    m::PerturbationModel,
) where {T}
    return print(
        io,
        "Perturbation Model: n_y = $(m.n_y), n_x = $(m.n_x), n_p = $(m.n_p), n_ϵ = $(m.n_ϵ), n_z = $(m.n_z)\n y = $(m.mod.y_symbols) \n x = $(m.mod.x_symbols) \n p = $(m.mod.p_symbols)",        
    )
end    

# abstract type AbstractFirstOrderPerturbationModel <:
#               AbstractPerturbationModel end
# abstract type AbstractSecondOrderPerturbationModel <:
#               AbstractPerturbationModel end

# function Base.show(
#     io::IO,
#     ::MIME"text/plain",
#     m::AbstractFirstOrderPerturbationModel,
# ) where {T}
#     return print(
#         io,
#         "First Order Model n_y = $(m.n_y), n_x = $(m.n_x), n_p = $(m.n_p), n_ϵ = $(m.n_ϵ), n_z = $(m.n_z)\n",
#     )
# end
# function Base.show(
#     io::IO,
#     ::MIME"text/plain",
#     m::AbstractSecondOrderPerturbationModel,
# ) where {T}
#     return print(
#         io,
#         "Second Order Model n_y = $(m.n_y), n_x = $(m.n_x), n_p = $(m.n_p), n_ϵ = $(m.n_ϵ), n_z = $(m.n_z)\n",
#     )
# end
# #TODO The reason this has all of the structure rather than just a reference to the module is so
# # (1) it can hold runtime genreated functions at some point and (2) otherwise inference is tricky
# # In particular, would need to have anything that can potentially change the return type as something inferable
# # such as anything which could be "nothing", etc.  This could be done, but would take some effort.
# # or, maybe it doesn't matter if the `generate_perturbation` is type stable.
# Base.@kwdef struct FirstOrderPerturbationModel{
#     T1<:AbstractMatrix,
#     T2,
#     T3,
#     T4,
#     T5,
#     T6,
#     T7,
#     T8,
#     T9,
#     T10,
#     T11,
#     T12,
#     T13,
#     T14,
#     T15,
#     T16,
#     T17,
#     T18,
#     T19,
#     T20,
#     T21,
#     T22,
#     T23,
#     T24,
#     T25,
#     T26,
# } <: AbstractFirstOrderPerturbationModel
#     # bookkeeping
#     n::Int64
#     n_y::Int64
#     n_x::Int64
#     n_p::Int64
#     n_ϵ::Int64
#     n_z::Int64
#     η::T1
#     Q::T2

#     # functions
#     Γ!::T3
#     Γ_p!::T4
#     Ω!::T5
#     Ω_p!::T6

#     # 0th order
#     H!::T7
#     # 1st order solutions.
#     H_yp!::T8
#     H_y!::T9
#     H_xp!::T10
#     H_x!::T11
#     H_yp_p!::T12
#     H_y_p!::T13
#     H_xp_p!::T14
#     H_x_p!::T15
#     H_p!::T16
#     Ψ!::T17 # need for the estimation of 1st order

#     # steady state
#     H̄!::T18
#     H̄_w!::T19
#     ȳ_iv!::T20
#     x̄_iv!::T21
#     ȳ!::T22
#     x̄!::T23
#     ȳ_p!::T24
#     x̄_p!::T25
#     steady_state!::T26
# end

# Base.@kwdef struct SecondOrderPerturbationModel{
#     T1<:AbstractMatrix,
#     T2,
#     T3,
#     T4,
#     T5,
#     T6,
#     T7,
#     T8,
#     T9,
#     T10,
#     T11,
#     T12,
#     T13,
#     T14,
#     T15,
#     T16,
#     T17,
#     T18,
#     T19,
#     T20,
#     T21,
#     T22,
#     T23,
#     T24,
#     T25,
#     T26,
#     T27,
#     T28,
#     T29,
#     T30,
#     T31,
# } <: AbstractSecondOrderPerturbationModel
#     # bookkeeping
#     n::Int64
#     n_y::Int64
#     n_x::Int64
#     n_p::Int64
#     n_ϵ::Int64
#     n_z::Int64
#     η::T1
#     Q::T2

#     # functions
#     Γ!::T3
#     Γ_p!::T4
#     Ω!::T5
#     Ω_p!::T6

#     # 0th order
#     H!::T7
#     # 1st order solutions.
#     H_yp!::T8
#     H_y!::T9
#     H_xp!::T10
#     H_x!::T11
#     H_yp_p!::T12
#     H_y_p!::T13
#     H_xp_p!::T14
#     H_x_p!::T15
#     H_p!::T16
#     Ψ!::T17 # need for the estimation of 1st order

#     # steady state
#     H̄!::T18
#     H̄_w!::T19
#     ȳ_iv!::T20
#     x̄_iv!::T21
#     ȳ!::T22
#     x̄!::T23
#     ȳ_p!::T24
#     x̄_p!::T25
#     steady_state!::T26

#     # 2nd order drivatives
#     Ψ_p!::T27
#     Ψ_yp!::T28
#     Ψ_y!::T29
#     Ψ_xp!::T30
#     Ψ_x!::T31
# end

# ## Structures to hold the reusable, mutable cache for the solvers
# abstract type AbstractSolverCache end
# abstract type AbstractFirstOrderSolverCache <: AbstractSolverCache end
# abstract type AbstractSecondOrderSolverCache <: AbstractSolverCache end


# allocate_cache(m::AbstractFirstOrderPerturbationModel) = FirstOrderSolverCache(m)
# allocate_cache(m::AbstractSecondOrderPerturbationModel) = SecondOrderSolverCache(m)

# Base.@kwdef mutable struct FirstOrderSolverCache{
#     MatrixType<:AbstractMatrix,
#     MatrixType2<:AbstractMatrix,
#     MatrixType3<:AbstractMatrix,
#     MatrixType4<:AbstractMatrix,
#     VectorType<:AbstractVector,
#     VectorOfMatrixType<:AbstractVector{<:AbstractMatrix},
#     VectorOfMatrixType2<:AbstractVector{<:AbstractMatrix},
#     VectorOfMatrixType3<:AbstractVector{<:AbstractMatrix},
#     VectorOrNothingType<:Union{Nothing,AbstractVector},
#     MatrixOrNothingType<:Union{Nothing,AbstractMatrix},
#     MatrixScalingOrNothingType<:Union{Nothing,AbstractMatrix,UniformScaling},
#     SymmetricMatrixType<:AbstractMatrix,
#     SymmetricVectorOfMatrixType<:AbstractVector{<:AbstractMatrix},
#     CholeskyType<:Cholesky,
#     ChangeVarianceType<:AbstractVector{<:AbstractMatrix},
# } <: AbstractFirstOrderSolverCache

#     H::VectorType
#     H_yp::MatrixType
#     H_y::MatrixType
#     H_xp::MatrixType
#     H_x::MatrixType
#     H_yp_p::VectorOfMatrixType
#     H_y_p::VectorOfMatrixType
#     H_xp_p::VectorOfMatrixType
#     H_x_p::VectorOfMatrixType
#     H_p::MatrixType
#     Γ::MatrixType2
#     Γ_p::VectorOfMatrixType2
#     Σ::SymmetricMatrixType
#     Σ_p::SymmetricVectorOfMatrixType
#     Ω::VectorOrNothingType
#     Ω_p::MatrixOrNothingType
#     Ψ::VectorOfMatrixType

#     # Used in solution
#     x::VectorType
#     y::VectorType
#     y_p::MatrixType4  # usually dense
#     x_p::MatrixType4
#     g_x::MatrixType4
#     h_x::MatrixType4
#     g_x_p::VectorOfMatrixType3  # usually vector of dense
#     h_x_p::VectorOfMatrixType3
#     B::MatrixType2
#     B_p::VectorOfMatrixType2
#     Q::MatrixScalingOrNothingType
#     η::MatrixType3 # might not be floating points
#     A_1_p::VectorOfMatrixType3
#     C_1::MatrixType4
#     C_1_p::VectorOfMatrixType3
#     V::CholeskyType
#     V_p::ChangeVarianceType
# end

# # initialized to zero rather than undef since the MTK genreated functions don't replace zeros
# # This is the dense matrix default behavior.  Other algorithm types can add additional constructors
# function FirstOrderSolverCache(m::AbstractFirstOrderPerturbationModel)
#     @unpack n_x, n_y, n, n_p, n_ϵ, n_z = m

#     return FirstOrderSolverCache(;
#         H = zeros(n),
#         H_yp = zeros(n, n_y),
#         H_y = zeros(n, n_y),
#         H_xp = zeros(n, n_x),
#         H_x = zeros(n, n_x),
#         Γ = zeros(n_ϵ, n_ϵ),
#         Ω = isnothing(m.Ω!) ? nothing : zeros(n_z),
#         Ψ = [zeros(2n, 2n) for i = 1:n],
#         H_p = zeros(n, n_p),
#         H_yp_p = [zeros(n, n_y) for i = 1:n_p],
#         H_y_p = [zeros(n, n_y) for i = 1:n_p],
#         H_xp_p = [zeros(n, n_x) for i = 1:n_p],
#         H_x_p = [zeros(n, n_x) for i = 1:n_p],
#         Γ_p = [zeros(n_ϵ, n_ϵ) for i = 1:n_p],
#         Ω_p = isnothing(m.Ω_p!) ? nothing : zeros(n_z, n_p),
#         x = zeros(n_x),
#         y = zeros(n_y),
#         y_p = zeros(n_y, n_p),
#         x_p = zeros(n_x, n_p),
#         g_x = zeros(n_y, n_x),
#         h_x = zeros(n_x, n_x),
#         g_x_p = [zeros(n_y, n_x) for _ = 1:n_p],
#         h_x_p = [zeros(n_x, n_x) for _ = 1:n_p],
#         Σ = Symmetric(zeros(n_ϵ, n_ϵ)),
#         Σ_p = [Symmetric(zeros(n_ϵ, n_ϵ)) for _ = 1:n_p],
#         m.Q,
#         m.η,
#         C_1 = zeros(n_z, n_x),
#         A_1_p = [zeros(n_x, n_x) for _ = 1:n_p],
#         C_1_p = [zeros(n_z, n_x) for _ = 1:n_p],
#         V = cholesky(Array(I(n_x))),
#         V_p = [zeros(n_x, n_x) for _ = 1:n_p],
#         B = zeros(n_x, n_ϵ),
#         B_p = [zeros(n_x, n_ϵ) for _ = 1:n_p],
#     )
# end

# Base.@kwdef mutable struct SecondOrderSolverCache{
#     MatrixType<:AbstractMatrix,
#     MatrixType2<:AbstractMatrix,
#     MatrixType3<:AbstractMatrix,
#     MatrixType4<:AbstractMatrix,
#     VectorType<:AbstractVector,
#     VectorOfMatrixType<:AbstractVector{<:AbstractMatrix},
#     VectorOfMatrixType2<:AbstractVector{<:AbstractMatrix},
#     VectorOfMatrixType3<:AbstractVector{<:AbstractMatrix},
#     VectorOrNothingType<:Union{Nothing,AbstractVector},
#     MatrixOrNothingType<:Union{Nothing,AbstractMatrix},
#     MatrixScalingOrNothingType<:Union{Nothing,AbstractMatrix,UniformScaling},
#     SymmetricMatrixType<:AbstractMatrix,
#     SymmetricVectorOfMatrixType<:AbstractVector{<:AbstractMatrix},
#     VectorOfVectorOfMatrixType<:AbstractVector{<:AbstractVector},
#     ThreeTensorType<:Array{<:Number,3},
#     CholeskyType<:Cholesky,
#     ChangeVarianceType<:AbstractVector{<:AbstractMatrix},
#     VectorOfThreeTensorType<:AbstractVector{<:Array{<:Number,3}},
# } <: AbstractSecondOrderSolverCache
#     H::VectorType
#     H_yp::MatrixType
#     H_y::MatrixType
#     H_xp::MatrixType
#     H_x::MatrixType
#     H_yp_p::VectorOfMatrixType
#     H_y_p::VectorOfMatrixType
#     H_xp_p::VectorOfMatrixType
#     H_x_p::VectorOfMatrixType
#     H_p::MatrixType
#     Γ::MatrixType2
#     Γ_p::VectorOfMatrixType2
#     Σ::SymmetricMatrixType
#     Σ_p::SymmetricVectorOfMatrixType
#     Ω::VectorOrNothingType
#     Ω_p::MatrixOrNothingType
#     Ψ::VectorOfMatrixType

#     # Used in solution
#     x::VectorType
#     y::VectorType
#     y_p::MatrixType4  # usually dense
#     x_p::MatrixType4
#     g_x::MatrixType4
#     h_x::MatrixType4
#     g_x_p::VectorOfMatrixType3  # usually vector of dense
#     h_x_p::VectorOfMatrixType3
#     B::MatrixType2
#     B_p::VectorOfMatrixType2
#     Q::MatrixScalingOrNothingType
#     η::MatrixType3 # might not be floating points

#     # Additional for 2nd order
#     Ψ_p::VectorOfVectorOfMatrixType
#     Ψ_yp::VectorOfVectorOfMatrixType
#     Ψ_y::VectorOfVectorOfMatrixType
#     Ψ_xp::VectorOfVectorOfMatrixType
#     Ψ_x::VectorOfVectorOfMatrixType
#     g_xx::ThreeTensorType
#     h_xx::ThreeTensorType
#     g_σσ::VectorType
#     h_σσ::VectorType
#     g_xx_p::VectorOfThreeTensorType
#     h_xx_p::VectorOfThreeTensorType
#     g_σσ_p::MatrixType4
#     h_σσ_p::MatrixType4

#     A_1_p::VectorOfMatrixType3
#     A_0_p::MatrixType4
#     A_2_p::VectorOfThreeTensorType
#     C_1::MatrixType4
#     C_0::VectorType
#     C_2::ThreeTensorType
#     C_1_p::VectorOfMatrixType3
#     C_0_p::MatrixType4
#     C_2_p::VectorOfThreeTensorType

#     V::CholeskyType
#     V_p::ChangeVarianceType
# end

# function SecondOrderSolverCache(m::AbstractSecondOrderPerturbationModel)
#     @unpack n_x, n_y, n, n_p, n_ϵ, n_z = m

#     return SecondOrderSolverCache(;
#         H = zeros(n),
#         H_yp = zeros(n, n_y),
#         H_y = zeros(n, n_y),
#         H_xp = zeros(n, n_x),
#         H_x = zeros(n, n_x),
#         Γ = zeros(n_ϵ, n_ϵ),
#         Ω = isnothing(m.Ω!) ? nothing : zeros(n_z),
#         Ψ = [zeros(2n, 2n) for i = 1:n],
#         H_p = zeros(n, n_p),
#         H_yp_p = [zeros(n, n_y) for i = 1:n_p],
#         H_y_p = [zeros(n, n_y) for i = 1:n_p],
#         H_xp_p = [zeros(n, n_x) for i = 1:n_p],
#         H_x_p = [zeros(n, n_x) for i = 1:n_p],
#         Γ_p = [zeros(n_ϵ, n_ϵ) for i = 1:n_p],
#         Ω_p = isnothing(m.Ω_p!) ? nothing : zeros(n_z, n_p),
#         x = zeros(n_x),
#         y = zeros(n_y),
#         y_p = zeros(n_y, n_p),
#         x_p = zeros(n_x, n_p),
#         g_x = zeros(n_y, n_x),
#         h_x = zeros(n_x, n_x),
#         g_x_p = [zeros(n_y, n_x) for _ = 1:n_p],
#         h_x_p = [zeros(n_x, n_x) for _ = 1:n_p],
#         Σ = Symmetric(zeros(n_ϵ, n_ϵ)),
#         Σ_p = [Symmetric(zeros(n_ϵ, n_ϵ)) for _ = 1:n_p],
#         m.Q,
#         m.η,
#         B = zeros(n_x, n_ϵ),
#         B_p = [zeros(n_x, n_ϵ) for _ = 1:n_p],
#         g_xx = zeros(n_y, n_x, n_x),
#         h_xx = zeros(n_x, n_x, n_x),
#         g_σσ = zeros(n_y),
#         h_σσ = zeros(n_x),
#         Ψ_yp = [[zeros(2n, 2n) for _ = 1:n] for _ = 1:n_y],
#         Ψ_y = [[zeros(2n, 2n) for _ = 1:n] for _ = 1:n_y],
#         Ψ_xp = [[zeros(2n, 2n) for _ = 1:n] for _ = 1:n_x],
#         Ψ_x = [[zeros(2n, 2n) for _ = 1:n] for _ = 1:n_x],
#         Ψ_p = [[zeros(2n, 2n) for _ = 1:n] for _ = 1:n_p],
#         g_xx_p = [zeros(n_y, n_x, n_x) for _ = 1:n_p],
#         h_xx_p = [zeros(n_x, n_x, n_x) for _ = 1:n_p],
#         g_σσ_p = zeros(n_y, n_p),
#         h_σσ_p = zeros(n_x, n_p),
#         C_1 = zeros(n_z, n_x),
#         C_1_p = [zeros(n_z, n_x) for _ = 1:n_p],
#         C_0 = zeros(n_z),
#         C_0_p = zeros(n_z, n_p),
#         C_2 = zeros(n_z, n_x, n_x),
#         C_2_p = [zeros(n_z, n_x, n_x) for _ = 1:n_p],
#         A_0_p = zeros(n_x, n_p),
#         A_1_p = [zeros(n_x, n_x) for _ = 1:n_p],
#         A_2_p = [zeros(n_x, n_x, n_x) for _ = 1:n_p],
#         V = cholesky(Array(I(n_x))),
#         V_p = [zeros(n_x, n_x) for _ = 1:n_p],
#     )
# end


# Base.@kwdef struct PerturbationSolverSettings{T1,T2,T3,T4,T5,T6}
#     print_level::Int64 = 1  # 0 is no output at all
#     ϵ_BK::Float64 = 1e-6 # For checking Blanchard-Kahn condition
#     nlsolve_method::Symbol = :trust_region
#     nlsolve_iterations::Int64 = 1000
#     nlsolve_show_trace::Bool = false
#     nlsolve_ftol::Float64 = 1e-8
#     use_solution_cache::Bool = true
#     evaluate_functions_callback::T1 = nothing
#     calculate_steady_state_callback::T2 = nothing
#     solve_first_order_callback::T3 = nothing
#     solve_first_order_p_callback::T4 = nothing
#     solve_second_order_callback::T5 = nothing
#     solve_second_order_p_callback::T6 = nothing
# end

# function nlsolve_options(s::PerturbationSolverSettings)
#     return (
#         method = s.nlsolve_method,
#         iterations = s.nlsolve_iterations,
#         show_trace = s.nlsolve_show_trace,
#         ftol = s.nlsolve_ftol,
#     )
# end

# struct PerturbationSolver{T1,T2,T3}
#     model::T1
#     cache::T2
#     settings::T3
# end

# # State Space types
# abstract type AbstractPerturbationSolution end
# abstract type AbstractFirstOrderPerturbationSolution <: AbstractPerturbationSolution end

# # For this, all are dense due to schur decomposition
# # All are dense due to schur decomposition
# Base.@kwdef struct FirstOrderPerturbationSolution{
#     T1<:AbstractVector,
#     T2<:AbstractVector,
#     T3<:AbstractMatrix,
#     T4<:AbstractMatrix,
#     T5<:AbstractMatrix,
#     T6<:Union{Nothing,Distribution},
#     T7<:Union{Nothing,AbstractMatrix,UniformScaling},
#     T8<:AbstractMatrix,
#     T9<:AbstractMatrix,
#     T10<:Distribution,
#     T11<:AbstractMatrix,
# } <: AbstractFirstOrderPerturbationSolution

#     retcode::Symbol
#     n::Int64
#     n_y::Int64
#     n_x::Int64
#     n_p::Int64
#     n_ϵ::Int64
#     n_z::Int64
#     y::T1
#     x::T2
#     g_x::T3
#     A::T4
#     B::T5
#     D::T6  # current a matrix or nothing, later could make more general
#     Q::T7  # can be nothing
#     η::T8
#     C::T9  # i.e. Q * g_x
#     x_ergodic::T10
#     Γ::T11
# end


# maybe_diagonal(x::AbstractVector) = TuringDiagMvNormal(zero(x), x)
# maybe_diagonal(x) = x # otherwise, just return raw.  e.g. nothing

# function FirstOrderPerturbationSolution(
#     retcode,
#     m::AbstractFirstOrderPerturbationModel,
#     c::FirstOrderSolverCache,
# )
#     return FirstOrderPerturbationSolution(;
#         retcode,
#         m.n_x,
#         m.n_y,
#         m.n_p,
#         m.n_ϵ,
#         m.n,
#         m.n_z,
#         c.Q,
#         c.η,
#         c.y,
#         c.x,
#         c.B,
#         D = maybe_diagonal(c.Ω),
#         c.g_x,
#         A = c.h_x,
#         C = c.C_1,
#         x_ergodic = TuringDenseMvNormal(zeros(m.n_x), c.V),
#         c.Γ,
#     )
# end

# Base.@kwdef struct SecondOrderPerturbationSolution{
#     T1<:AbstractVector,
#     T2<:AbstractVector,
#     T3<:AbstractMatrix,
#     T4<:AbstractMatrix,
#     T5<:AbstractMatrix,
#     T6<:Union{Nothing,Distribution},
#     T7<:Union{Nothing,AbstractMatrix,UniformScaling},
#     T8<:AbstractMatrix,
#     T9<:AbstractMatrix,
#     T10,
#     T11<:AbstractArray,
#     T12,
#     T13<:AbstractVector,
#     T14<:AbstractMatrix,
#     T15<:AbstractVector,
#     T16<:AbstractArray,
# } <: AbstractPerturbationSolution

#     retcode::Symbol
#     n::Int64
#     n_y::Int64
#     n_x::Int64
#     n_p::Int64
#     n_ϵ::Int64
#     n_z::Int64
#     y::T1
#     x::T2
#     g_x::T3
#     B::T5
#     D::T6
#     Q::T7  # can be nothing
#     η::T8
#     Γ::T9

#     g_xx::T10
#     g_σσ::T12
#     A_0::T13
#     A_1::T4
#     A_2::T11

#     C_1::T14
#     C_0::T15
#     C_2::T16
# end

# function SecondOrderPerturbationSolution(
#     retcode,
#     m::AbstractSecondOrderPerturbationModel,
#     c::SecondOrderSolverCache,
# )
#     return SecondOrderPerturbationSolution(;
#         retcode,
#         m.n_x,
#         m.n_y,
#         m.n_p,
#         m.n_ϵ,
#         m.n,
#         m.n_z,
#         c.Q,
#         c.η,
#         c.y,
#         c.x,
#         c.B,
#         D = maybe_diagonal(c.Ω),
#         c.Γ,
#         c.g_x,
#         A_1 = c.h_x,
#         c.g_xx,
#         A_2 = 0.5 * c.h_xx,
#         c.g_σσ,
#         A_0 = 0.5 * c.h_σσ,
#         c.C_1,
#         c.C_0,
#         c.C_2,
#     )
# end
