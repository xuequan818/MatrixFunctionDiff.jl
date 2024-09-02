using PermutationSymmetricTensors: SymmetricTensor

@inline function Base.getindex(A::SymmetricTensor{T, N, dim}, I::Int64...) where {T, dim, N}
    boundscheck_ex1 = :(@boundscheck ((I[1]>length(A) || I[1]<1) && throw(BoundsError(A, I))))
    if dim == 1 && length(I) == 1 
        index_ex = :(@inbounds A.data[I[1]])
        return :($boundscheck_ex1; $index_ex)
    end
    if dim == 1 && length(I) == 2 
        check_ex = :(I[2] == 1 || throw(DimensionMismatch("This $dim-dimensional symmetric tensor is being indexed with $(length(I)) indices.")))
        index_ex = :(@inbounds A.data[I[1]])
        return :($boundscheck_ex1; $check_ex; $index_ex)
    end
    if dim > 1 && length(I) == 1 
        if big(N)^dim > typemax(Int)
            index_ex = :(@inbounds A[Int128(I[1])])
            return :($boundscheck_ex1; $index_ex)
        else
            index_ex = :(@inbounds A[CartesianIndices(A)[I[1]]])
            return :($boundscheck_ex1; $index_ex)
        end
    end
    if length(I) != dim
        return :( throw(DimensionMismatch("This $dim-dimensional symmetric tensor is being indexed with $(length(I)) indices.")))
    end
    ex = :(I2 = sort(SVector(I), rev=true))
    ex1 = :(@boundscheck (I2[1]>N || I2[end]<1) && throw(BoundsError(A, I))) 
    ex2 = :(ind = 0; lin_ind=A.linear_indices)
    for i in 1:dim
        ex2 = :($ex2; @inbounds ind += lin_ind[$i][I2[$i]])
    end
    ex3 = :(@inbounds A.data[ind])
    return ex = :($ex; $ex1; $ex2; $ex3)
end