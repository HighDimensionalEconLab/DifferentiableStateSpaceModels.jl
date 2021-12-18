function sylvds!(A::AbstractMatrix{T1}, B::AbstractMatrix{T1}, C::AbstractMatrix{T1}, W::AbstractMatrix{T1} = similar(A,size(A,1),2)) where T1<:BlasReal
    m, n = LinearAlgebra.checksquare(A,B)
    (size(C,1) == m && size(C,2) == n ) || throw(DimensionMismatch("C must be an $m x $n matrix"))
    (m, 2) == size(W) || throw(DimensionMismatch("W must be an $m x 2 matrix"))
    ONE = one(T1)
 
    # determine the structure of the real Schur form of A
    ba, pa = MatrixEquations.sfstruct(A)
    bb, pb = MatrixEquations.sfstruct(B)
    
    G = Matrix{T1}(undef,2,2)
    WA = Matrix{T1}(undef,2,2)
 
    Xw = Matrix{T1}(undef,4,4)
    if !adjA && !adjB
       # """
       # The (K,L)th block of X is determined starting from
       # bottom-left corner column by column by
 
       #            A(K,K)*X(K,L)*B(L,L) + X(K,L) = C(K,L) - R(K,L)
 
       # where
       #                        M
       #            R(K,L) = { SUM [A(K,J)*X(J,L)] } * B(L,L) +
       #                      J=K+1
       #                        M             L-1
       #                       SUM { A(K,J) * SUM [X(J,I)*B(I,L)] }.
       #                       J=K            I=1
       # """
       j = 1
       for ll = 1:pb
           dl = bb[ll]
           dll = 1:dl
           il1 = 1:j-1
           j1 = j+dl-1
           l = j:j1
           i = m
           for kk = pa:-1:1
               dk = ba[kk]
               dkk = 1:dk
               i1 = i-dk+1
               k = i1:i
               Ckl = view(C,k,l)
               y = view(G,1:dk,1:dl)
               copyto!(y,Ckl)
               if kk < pa
                  ir = i+1:m
                  W1 = view(WA,dkk,dll)
                  mul!(W1,view(A,k,ir),view(C,ir,l))
                  mul!(y,W1,view(B,l,l),-ONE,ONE)
               end
               if ll > 1
                  ic = i1:m
                  mul!(view(W,k,dll),view(C,k,il1),view(B,il1,l))
                  mul!(y,view(A,k,ic),view(W,ic,dll),-ONE,ONE)
               end
               sylvd2!(adjA,adjB,y,dk,dl,view(A,k,k),view(B,l,l),Xw)  
               copyto!(Ckl,y)
            i -= dk
           end
           j += dl
       end