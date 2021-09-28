using Tullio, BenchmarkTools, LinearAlgebra

A = rand(100,100)
A_transpose = Array(A')
B = rand(100,100)
C = zeros(100,100)
C_vec = vec(copy(C))
swapspace = zeros(100,100)

function quadform_1!(C, A, B)
    C .= A' * B * A
end

function quadform_2!(C, swapspace, A, B)
    mul!(C, mul!(swapspace, A', B), A)
end

@benchmark quadform_1!($C, $A, $B)
@benchmark quadform_2!($C, $swapspace, $A, $A_transpose, $B)
