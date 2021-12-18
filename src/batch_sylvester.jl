function batch_sylvds!(A::AbstractMatrix{T1}, B::AbstractMatrix{T1}, C::AbstractVector{Matrix{T1}}) where T1<:LinearAlgebra.BlasReal
    n_p = length(C)
    W = [similar(A, size(A, 1), 2) for _ in 1:n_p]

    m, n = LinearAlgebra.checksquare(A,B)
    (size(C[1], 1) == m && size(C[1], 2) == n) || throw(DimensionMismatch("C matrices must be of size $m by $n"))
    (m, 2) == size(W[1]) || throw(DimensionMismatch("W matrices must be of size $m by 2"))
    ONE = one(T1)
 
    # determine the structure of the real Schur form of A
    ba, pa = MatrixEquations.sfstruct(A)
    bb, pb = MatrixEquations.sfstruct(B)
    
    G = [Matrix{T1}(undef, 2, 2) for _ in 1:n_p]
    WA = Matrix{T1}(undef, 2, 2)
    Xw = Matrix{T1}(undef, 4, 4)

    j = 1
    for ll = 1:pb
        dl = bb[ll]
        dll = 1:dl
        il1 = 1:j-1
        j1 = j + dl - 1
        l = j:j1
        i = m
        for kk = pa:-1:1
            dk = ba[kk]
            dkk = 1:dk
            i1 = i - dk + 1
            k = i1:i
            for p in 1:n_p
                Ckl = view(C[p], k, l)
                y = view(G[p], 1:dk, 1:dl)
                copyto!(y, Ckl)

                if kk < pa
                    ir = i+1:m
                    W1 = view(WA, dkk, dll)
                    mul!(W1, view(A, k, ir), view(C[p], ir, l))
                    mul!(y, W1, view(B, l, l), -ONE, ONE)
                end

                if ll > 1
                    ic = i1:m
                    mul!(view(W[p], k, dll), view(C[p], k, il1), view(B, il1, l))
                    mul!(y, view(A, k, ic), view(W[p], ic, dll), -ONE, ONE)
                end

                MatrixEquations.sylvd2!(false, false, y, dk, dl, view(A, k, k), view(B, l, l), Xw)  
                copyto!(Ckl, y)
            end
            i -= dk
        end
        j += dl
    end
    return C
end
