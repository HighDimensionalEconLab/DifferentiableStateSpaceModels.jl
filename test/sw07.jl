using DifferentiableStateSpaceModels,
    SparseArrays,
    LinearAlgebra,
    Parameters,
    Test,
    TimerOutputs,
    BenchmarkTools

@testset "Dense SW First Order" begin
    m = @include_example_module(Examples.SW07)

    p = [
        10,
        0.51,
        10,
        0,
        0.7,
        0.742,
        0,
        0,
        0.24,
        0.2696,
        6.0144,
        0.025,
        1.5,
        0.6361,
        1.5,
        0.3243,
        0.8087,
        0.47,
        0.6,
        1.9423,
        1.5,
        1.488,
        0.2347,
        0.0593,
        0.8762,
        0.9977,
        0.5799,
        0.9957,
        0.7165,
        0,
        0,
        0,
        0.3982,
        0.18,
    ]
    generate_perturbation(m, p)
    reset_timer!()
    sol = generate_perturbation(m, p)
    print_timer()
    @test sol.retcode == :Success
end
