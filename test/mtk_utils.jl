using DifferentiableStateSpaceModels, ModelingToolkit, SparseArrays, Test,
      GeneralizedGenerated
using DifferentiableStateSpaceModels: match_equation_lhs, check_markov_name_convention,
                                      check_ss_name_convention, make_markov_variable,
                                      make_ss_variable

@testset "Utilities Tests" begin
    @variables a b c α β γ
    eq = a ~ 2 + b + c
    @test match_equation_lhs(a, eq) == true
    @test match_equation_lhs(b, eq) == false

    @test isequal(Variable(:a_p), make_markov_variable(a))
    @test isequal(Variable(:α_p), make_markov_variable(α))
    @test isequal(Variable(:a_ss), make_ss_variable(a))
    @test isequal(Variable(:α_ss), make_ss_variable(α))

    @variables abc_p _p abc_ss _ss αβγ_p αβγ_ss

    @test_throws ArgumentError check_markov_name_convention(abc_p)
    @test_throws ArgumentError check_markov_name_convention(αβγ_p)
    @test_throws ArgumentError check_markov_name_convention(_p)
    @test_throws ArgumentError check_ss_name_convention(abc_ss)
    @test_throws ArgumentError check_ss_name_convention(αβγ_ss)
    @test_throws ArgumentError check_ss_name_convention(_ss)

    sol = connect_markov_variables([a, c])
    @test isequal(Variable(:a_ss), sol[1][3])
    @test isequal(Variable(:c_p), sol[2][2])
end