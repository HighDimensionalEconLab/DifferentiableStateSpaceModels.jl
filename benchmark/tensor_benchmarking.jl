using LinearAlgebra, Tullio, Test, ForwardDiff, Random, Zygote, BenchmarkTools
using LoopVectorization, TensorOperations, TensorCast, CUDA
using Zygote: @adjoint
# 2-tensors
# from Michael Abbott
# donate to ChainRules
@adjoint function dot(x,A,y)
    Ay = A * y
    z = adjoint(x) * Ay
    z, dz -> begin
        dx = lmul!(dz, conj!(Ay))
        dA = dz .* x .* adjoint(y)
        dy = lmul!(dz, adjoint(A) * x)
        (dx, dA, dy)
    end
end


quad_form_raw(x, A::AbstractMatrix) =  x' * A * x
quad_form_ten(x, A::AbstractArray{<:Number,2}) = @tullio c := x[j] * A[j,k] * x[k]
quad_form_dot(x, A::AbstractMatrix) = dot(x, A, x)

# 3-tensors
quad_form_raw(x, A::AbstractArray{<:Number,3}) = [x' * A[i,:,:] * x for i in 1:size(A,1)]
quad_form_mul(x, A::AbstractArray{<:Number,3}) = @matmul c[l] := sum(j) (@matmul [l,j] := sum(k) A[l,j,k] * x[k]) * x[j]
quad_form_avx(x, A::AbstractArray{<:Number,3}) = @tullio c[l] := x[j] * A[l,j,k] * x[k] tensor=false;
quad_form_ten(x, A::AbstractArray{<:Number,3}) = @tullio c[l] := x[j] * A[l,j,k] * x[k];
quad_form_base(x, A::AbstractArray{<:Number,3}) = @tullio c[l] := x[j] * A[l,j,k] * x[k] tensor=false avx=false;

#tests 
scalar_objective(f, x, A) = norm(f(x, A))

function print_benchmarks(N, use_cuda = true)
    A_2 = rand(N,N)
    A_3 = rand(N,N,N)
    x = rand(N)
    # 2-tensors
    printstyled("Raw Julia 2-tensor for N = $N\n", color=:red)
    @btime scalar_objective(quad_form_raw, $x, $A_2)
    @btime gradient(scalar_objective, quad_form_raw, $x, $A_2)

    printstyled("Dot Julia 2-tensor for N = $N\n", color=:red)
    @btime scalar_objective(quad_form_dot, $x, $A_2)
    @btime gradient(scalar_objective, quad_form_dot, $x, $A_2)    

    printstyled("Tullio with Tensor Julia 2-tensor for N = $N\n", color=:red)
    @btime scalar_objective(quad_form_ten, $x, $A_2)
    @btime gradient(scalar_objective, quad_form_ten, $x, $A_2)

    # three tensor
    printstyled("Raw Julia 3-tensor for N = $N\n", color=:red)
    @btime scalar_objective(quad_form_raw, $x, $A_3)
    @btime gradient(scalar_objective, quad_form_raw, $x, $A_3)

    printstyled("TensorCast Julia 3-tensor for N = $N\n", color=:red)
    @btime scalar_objective(quad_form_mul, $x, $A_3)
    @btime gradient(scalar_objective, quad_form_mul, $x, $A_3)

    printstyled("AVX Julia 3-tensor for N = $N\n", color=:red)
    @btime scalar_objective(quad_form_avx, $x, $A_3)
    @btime gradient(scalar_objective, quad_form_avx, $x, $A_3)

    printstyled("Tullio with Tensor Julia 3-tensor for N = $N\n", color=:red)
    @btime scalar_objective(quad_form_ten, $x, $A_3)
    @btime gradient(scalar_objective, quad_form_ten, $x, $A_3)

    printstyled("Tullio Base Julia 3-tensor for N = $N\n", color=:red)
    @btime scalar_objective(quad_form_base, $x, $A_3)
    @btime gradient(scalar_objective, quad_form_base, $x, $A_3)


    if use_cuda
        A_3 = CuArray(A_3)
        x = CuArray(x)
        
        printstyled("CUDA TensorCast Julia 3-tensor for N = $N\n", color=:red)
        @btime CUDA.@sync blocking=false scalar_objective(quad_form_mul, $x, $A_3)
        @btime CUDA.@sync blocking=false gradient(scalar_objective, quad_form_mul, $x, $A_3)
    
        # Not sure why it is not working?
        # printstyled("CUDA Tullio  Julia 3-tensor for N = $N\n", color=:red)
        # @btime CUDA.@sync blocking=false scalar_objective(quad_form_ten, $x, $A_3)
        # @btime CUDA.@sync blocking=false gradient(scalar_objective, quad_form_ten, $x, $A_3)        
    end
end

Zygote.refresh()
# try for different sizes
print_benchmarks(100, false)
#print_benchmarks(100)
#print_benchmarks(200)