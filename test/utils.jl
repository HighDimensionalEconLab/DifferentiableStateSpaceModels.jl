using DifferentiableStateSpaceModels, Test, LinearAlgebra
using DifferentiableStateSpaceModels: vech, inv_vech
@testset "vech" begin
    A = [1 2; 3 4]
    @test vech(A) == [1, 3, 4]
    @test vech(LowerTriangular(A)) == [1, 3, 4]
    @test_throws DimensionMismatch vech([1 2 4; 3 4 5])

    @test inv_vech([1, 3, 4], 2) == LowerTriangular(A)
    @test inv_vech([1, 3, 4]) == LowerTriangular(A)
    @test_throws ErrorException inv_vech([1, 3, 4, 5])
    @test_throws ErrorException inv_vech([1, 3, 4], 3)

    B = LowerTriangular(rand(10, 10))
    @test inv_vech(vech(B)) ≈ B
end
