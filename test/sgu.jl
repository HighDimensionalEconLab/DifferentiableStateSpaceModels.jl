using DifferentiableStateSpaceModels, ModelingToolkit, SparseArrays, LinearAlgebra,
      Parameters, Test, TimerOutputs, BenchmarkTools

@testset "Dense SGU First Order" begin
    m = @include_example_module(Examples.sgusmallopen)

    p = [2.0, 1.455, 0.42, 0.0129, 0.1, 0.000742, 0.32, 0.028, 1.0 / (1.0 + 0.04), 0.04,
         0.7442]   #From Cesa-Bianchi (2012)
    generate_perturbation(m, p)
    reset_timer!()
    sol = generate_perturbation(m, p)
    print_timer()
    @test sol.retcode == :Success
end

@testset "Dense SGU Second Order" begin
    m = @include_example_module(Examples.sgusmallopen, 2)

    p = [2.0, 1.455, 0.42, 0.0129, 0.1, 0.000742, 0.32, 0.028, 1.0 / (1.0 + 0.04), 0.04,
         0.7442]   #From Cesa-Bianchi (2012)
    generate_perturbation(m, p)
    reset_timer!()
    sol = generate_perturbation(m, p)
    print_timer()
    @test sol.retcode == :Success
end

# @testset "Sparse SGU First Order" begin
#     H, mod_vals, model_name = Examples.sgusmallopen()
#     model = FirstOrderPerturbationModel(H; functions_type = SparseFunctions(),
#                                                    mod_vals...)

#     p = [2.0, 1.455, 0.42, 0.0129, 0.1, 0.000742, 0.32, 0.028, 1.0 / (1.0 + 0.04), 0.04,
#          0.7442]   #From Cesa-Bianchi (2012)
#     generate_perturbation(model, p)
#     reset_timer!()
#     sol = generate_perturbation(model, p)
#     print_timer()
#     @test sol.retcode == :Success
#     # TODO: Need smoke tests
# end
