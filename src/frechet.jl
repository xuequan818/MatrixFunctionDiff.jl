"""
    mat_fun_frechet(f, eigs, Ψ::AbstractMatrix, hs::Vector{AbstractMatrix})
    mat_fun_frechet(f, H::AbstractMatrix, hs::Vector{AbstractMatrix})

Return the n-th order Fréchet derivative `d^nf(H)hs[1]…hs[n]`, assuming `f` is called as `f(x)`.
"""
@inline function mat_fun_frechet(f::Function, eigs::Vector{Float64},
                         Ψ::AbstractMatrix, 
                         hs::Vector{V}) where {V<:AbstractMatrix}
    N = length(eigs)
    order = length(hs)
    DD_F = DD_tensor(f, eigs, order)
    hs = map(x -> Ψ' * x * Ψ, hs)

    # F_1 =  h_1 ∘ Λ^{0,1}
    if order == 1
        return Ψ * (hs[1] .* DD_F) * Ψ'
    end

    N0 = N^(order - 2)
    DD_F = reshape(DD_F, N, N, N, N0)

    T = promote_type(eltype(Ψ), eltype(V))
    val = zeros(T, N, N)
    h12 = zeros(T, N, N, N)

    # loop for the permutations
    pert = collect(permutations(1:order))
    for p in pert 
        @views hp = hs[p]

        # compute {h}^{1,2}_{:,i,:} := (h_1)_{:,i}(h_2)_{i,:}
        @. h12 = zero(T)
        @views for i = 1:N
            h12[:, i, :] = hp[1][:, i] * transpose(hp[2][i, :])
        end

        # compute {F}^{0,2,…,n}_{:,:,j_3,…,j_n} := ∑_{i=1}^N ({h}^{1,2} ∘ Λ^{0,1,…,n}_{:,:,:,j_3,…,j_n})_{:,i,:}
        hF = zeros(T, N, N, N0)
        @views for i = 1:N0
            hF[:, :, i] = dropdims(sum(h12 .* DD_F[:, :, :, i], dims=2); dims=2)
        end
        hF = reshape(hF, ntuple(x -> N, order))
        hF = permutedims(hF, tuple(2:order..., 1))

        # compute the recursion 
        # {F}^{m,…,n,0}_{:,j_m,…,j_n} = ∑_{i=1}^N (h_m ∘ {F}^{m-1,…,n,0}_{:,:,j_m,…,j_n})_{i,:}
        for k = 3:order
            Nk = N^(order + 1 - k)
            DD_Fk = reshape(hF, N, N, Nk)
            hF = zeros(T, N, Nk)
            @views for i = 1:Nk
                hF[:, i] = dropdims(sum(hp[k] .* DD_Fk[:, :, i], dims=1), dims=1)
            end
        end

        val += hF 
    end

    Ψ * transpose(val) * Ψ'
end

@inline function mat_fun_frechet(f::Function, H::AbstractMatrix,
                         hs::Vector{V}) where {V<:AbstractMatrix}
    # only support hermitian matrices                        
    @assert norm(H - H') < eps(real(eltype(H)))
    for h in hs
        @assert norm(h - h') < eps(real(eltype(h)))
    end

    # diagonalize H
    # compute the full eigen pairs
    eigs, Ψ = eigen(H)

    # compute the frechet derivative by eigen pairs
    mat_fun_frechet(f, eigs, Ψ, hs)
end

# Generate the divided difference tensor 
# DD_F = f[λ_i0, λ_i2, ..., λ_in]
# By the permutation symmetry of the divided difference, 
# we just calculate the irreducible vals
function DD_tensor(f::Function, eigs::Vector{T}, order::Integer) where {T}
    N = length(eigs)
    dim = order + 1
    eigs_take(ind) = map(x -> eigs[x], ind)

    DD_F_sym_index = find_full_indices(N, dim)
    DD_F_sym_val = @. divided_difference(f, eigs_take(DD_F_sym_index))
    DD_F = SymmetricTensor(DD_F_sym_val, Val(N), Val(dim))
    DD_F_ARRAY = zeros(eltype(DD_F_sym_val), ntuple(x->N,dim))
    for (i, vali) in enumerate(DD_F)
        DD_F_ARRAY[i] = vali
    end

    return DD_F_ARRAY
end