using DifferentiableStateSpaceModels, BenchmarkTools, Test
# function test_first_order(p_d, p_f, m)
#     return generate_perturbation(m, p_d, p_f)#, Val(1); cache = c) # manually passing in order 
# end
# Some sort of inference issues.  Trouble putting in function and need the `const` for now
#@testset "grad_tests" begin
const m_other = @include_example_module(DifferentiableStateSpaceModels.Examples.rbc_observables)  # const fixes current bug.  Can't move inside
# const m = @include_example_module(DifferentiableStateSpaceModels.Examples.rbc_observables)  # const fixes current bug.  Can't move inside
m = PerturbationModel{1,1,1,1,1,1,true,Nothing,Nothing}(DifferentiableStateSpaceModels.ModuleWrapper(Main.rbc_observables),
                                                        1, 1, 1, 1, 1, 1, true, nothing,
                                                        nothing)

p_d = (α = 0.5, β = 0.95)
SolverCache(m, Val(1), p_d)
@inferred SolverCache(m, Val(1), p_d)
@btime SolverCache(m, Val(1), Val(2))

# Base.@kwdef struct SimpleSolverCache{Order,MatrixType}
#     order::Val{Order}  # allows inference in construction
#     H::MatrixType
# end

# function SimpleSolverCache(m::PerturbationModel{MaxOrder,N_y,N_x,N_ϵ,N_z,N_p,HasΩ,T1,T2},
#                            ::Val{Order},
#                            p_d) where {Order,MaxOrder,N_y,N_x,N_ϵ,N_z,N_p,HasΩ,T1,T2}
#     return SimpleSolverCache(; order = Val(Order), H = zeros(N_x, N_x))
# end

# SimpleSolverCache(m, Val(1), p_d)
# @btime SimpleSolverCache(m, Val(1), p_d)