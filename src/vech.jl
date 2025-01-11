### Convert symmetric matrix to vector
vech(A::Matrix{Float64}) = A[tril(trues(size(A)))]

### Convert vector to symmetric matrix
function vech2vec_matrix(n::Int)
    G = zeros(Int,n*n,Int(n*(n+1)/2))
    for i=1:n
        for j=1:i
            G[(j-1)*n+i,Int((j-1)*(2*n-j)/2+i)] = 1
            G[(i-1)*n+j,Int((j-1)*(2*n-j)/2+i)] = 1
        end
    end
    return G
end
function inv_vech(x::Vector{Float64})
    n::Int = Int((-1+sqrt(1+8*length(x)))/2)
    return reshape(vech2vec_matrix(n)*x,(n,n))
end