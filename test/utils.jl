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


@testset "thread local cache" begin
    cache = ThreadLocalCache([1,2,3])
    @test cache.caches.count == 1
    Threads.@threads for _ in 1:(Threads.nthreads()*2) # > number of julia threads
        c = cache()
        c[1] = Threads.threadid()
    end
    @test cache.caches.count == Threads.nthreads()
    @show collect(cache.caches)
end