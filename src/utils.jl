function vech(A::AbstractMatrix{T}) where {T}
    m = LinearAlgebra.checksquare(A)
    v = Vector{T}(undef, (m * (m + 1)) >> 1)
    k = 0
    for j in 1:m, i in j:m
        @inbounds v[k += 1] = A[i, j]
    end
    return v
end

function inv_vech(v::AbstractVector{T},
                  n = Int(round(1 / 2 * (sqrt(1 + 8 * length(v) - 1))))) where {T}
    n * (n + 1) / 2 == length(v) || error("length(v) != n(n+1)/2")
    A = LowerTriangular(zeros(T, n, n))
    indices = [0; cumsum(n:-1:1)]
    for j in 1:n, i in j:n
        @inbounds A[i, j] = v[indices[j] + i - (j - 1)]
    end
    return A
end

# inefficient comparision for unit testing.  Not for high performance code
# Does not require the same type
function all_fields_equal(x1::T1, x2::T2, fields) where {T1,T2}
    return all(getfield.(Ref(x1), fields) .== getfield.(Ref(x2), fields))
end

# Helpers for generation 

# conditionally call to support functions = nothing
maybe_call_function(f, args...) = f(args...)
maybe_call_function(::Nothing, args...) = nothing

# Rethrow if an exception is non-LAPACK -- we shouldn't fail silently for non-LAPACK errors
# The reason that LAPACK errors are simple rethrown is that they typically mean model is in
# a non-convergence region and new parameters should be drawn
function is_linear_algebra_exception(e)
    if e isa LAPACKException
        return true
    elseif e isa PosDefException
        return true
    else
        return false
    end
end