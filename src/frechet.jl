@doc row"""
    mat_fun_frechet(f, eigs, Ψ::AbstractMatrix, hs::Vector{AbstractMatrix})
    mat_fun_frechet(f, H::AbstractMatrix, hs::Vector{AbstractMatrix})

Return the n-th order Fréchet derivative `d^nf(H)hs_1…hs_n`, assuming `f` is called as `f(x)`.
Use array operations to efficiently compute the Fréchet derivative.
For simplicity, just consider the no permutation case and define
```math
(F_n)_{kℓ}:=∑_{i_1,⋯,i_{n-1}=1}^N(h_1)_{k,i_1}⋯ (h_n)_{i_{n-1},ℓ}Λ^{0,1,…,n-1,n}_{k,i_1,…,i_{n-1},ℓ},
```
where $Λ^{0,…,n}_{i_0,…,i_n} := f[λ_{i_0},⋯,λ_{i_n}]$.
It is immediately to obtain that $F_1 =  h_1 Λ^{0,1}$, 
and $F_2 = ∑_{i=1}^N (\mathfrak{h}^{1,2} ∘ Λ^{0,1,2})_{:,i,:}$ with $\mathfrak{h}^{1,2}_{:,i,:} := (h_1)_{:,i}(h_2)_{i,:}$.
For $n ≥ 3$, first compute 
```math
\mathfrak{F}^{0,2,…,n}_{:,:,j_3,…,j_n} := ∑_{i=1}^N (\mathfrak{h}^{1,2} ∘ Λ^{0,1,…,n}_{:,:,:,j_3,…,j_n})_{:,i,:}
```
and permute such that the $0$-dimension is at the end $
\mathfrak{F}^{2,…,n,0}$. 
Then for $m ≥ 3$, there is the recursion
```math
\mathfrak{F}^{m,…,n,0}_{:,j_m,…,j_n} = ∑_{i=1}^N (h_m ∘ \mathfrak{F}^{m-1,…,n,0}_{:,:,j_m,…,j_n})_{i,:}
```
and $F_n = (\mathfrak{F}^{n,0})^T$.
"""
function mat_fun_frechet(f::Function, eigs::Vector{Float64},
                         Ψ::AbstractMatrix, 
                         hs::Vector{V}) where {V<:AbstractMatrix}
    N = length(eigs)
    order = length(hs)
    DD_F = DD_tensor(f, eigs, order)
    hs = map(x -> Ψ' * x * Ψ, hs)

    if order == 1
        @. DD_F = DD_F .* hs[1]
        return Ψ * DD_F * Ψ'
    end

    N0 = N^(order - 2)
    DD_F = reshape(DD_F, N, N, N, N0)

    T = promote_type(eltype(Ψ), eltype(V))
    val = zeros(T, N, N)
    h12 = zeros(T, N, N, N)
    pert = collect(permutations(1:order))
    for p in pert
        @views hp = hs[p]

        @. h12 = zero(TT)
        @views for i = 1:N
            h12[:, i, :] = hp[1][:, i] * transpose(hp[2][i, :])
        end

        hF = zeros(TT, N, N, N0)
        @views for i = 1:N0
            hF[:, :, i] = dropdims(sum(h12 .* DD_F[:, :, :, i], dims=2); dims=2)
        end
        hF = reshape(hF, ntuple(x -> N, order))

        if order > 2
            hF = permutedims(hF, tuple(2:order..., 1))
        end

        for k = 3:order
            Nk = N^(order + 1 - k)
            DD_Fk = reshape(hF, N, N, Nk)
            hF = zeros(T, N, Nk)
            @views for i = 1:Nk
                hF[:, i] = dropdims(sum(hp[k] .* DD_Fk[:, :, i], dims=1), dims=1)
            end
        end

        # directly add hF^{n,0} by the symmetry of frechet derivative  
        val += hF 
    end

    Ψ * val * Ψ'
end

function mat_fun_frechet(f::Function, H::AbstractMatrix,
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
function DD_tensor(f::Function, eigs::Vector{T}, order::Int64) where {T}
    N = length(eigs)
    dim = order + 1
    eigs_take(ind) = map(x -> eigs[x], ind)

    DD_F_sym_index = find_full_indices(N, dim)
    DD_F_sym_val = @. divided_difference(f, eigs_take(DD_F_sym_index))
    DD_F = SymmetricTensor(DD_F_sym_val, Val(N), Val(dim))

    return DD_F
end