function vech(A::AbstractMatrix{T}) where {T}
    m = LinearAlgebra.checksquare(A)
    v = Vector{T}(undef, (m * (m + 1)) >> 1)
    k = 0
    for j = 1:m, i = j:m
        @inbounds v[k+=1] = A[i, j]
    end
    return v
end

function inv_vech(
    v::AbstractVector{T},
    n = Int(round(1 / 2 * (sqrt(1 + 8 * length(v) - 1)))),
) where {T}
    n * (n + 1) / 2 == length(v) || error("length(v) != n(n+1)/2")
    A = LowerTriangular(zeros(T, n, n))
    indices = [0; cumsum(n:-1:1)]
    for j = 1:n, i = j:n
        @inbounds A[i, j] = v[indices[j]+i-(j-1)]
    end
    return A
end

# inefficient comparision for unit testing.  Not for high performance code
all_equal_struct(x1::T, x2::T) where {T} =
    all(getfield.(Ref(x1), fieldnames(T)) .== getfield.(Ref(x2), fieldnames(T)))
