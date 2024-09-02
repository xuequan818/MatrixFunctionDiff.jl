module FrechetTest

using Test
using MatrixFunctionDiff
using DividedDifferences
using LinearAlgebra
using Combinatorics

# compute the frechet derivative by element loop
function frechet_naive(f::Function, eigs::Vector{Float64},
					   Ψ::AbstractMatrix, 
					   hs::Vector{V}) where {V<:AbstractMatrix}
    N = length(eigs)
    order = length(hs)
    hs = map(x -> Ψ' * x * Ψ, hs)

    pert = collect(permutations(1:order))
    T = promote_type(eltype(Ψ), eltype(V))
    val = zeros(T, N, N)
    for i = 1:N, j = 1:N
        ktr = ones(Int, order - 1)
        for k = 1:N^(order-1)
            kind = vcat(i, ktr, j)
            hval = zero(T)
            for p in pert
                hval += prod(l -> hs[p[l]][kind[l], kind[l+1]], 1:order)
            end
            ΛF = divided_difference(f, eigs[kind])
            val[i, j] += hval * ΛF

            if length(ktr) > 0
                ktr[1] += 1
                for ℓ = 1:order-2
                    if ktr[ℓ] == N + 1
                        ktr[ℓ] = 1
                        ktr[ℓ+1] += 1
                    end
                end
            end
        end
    end

    return Ψ * val * Ψ'
end

const N = 10
const Funs = [DividedDifferences.heaviside, x -> 1.0 / (1.0 + exp(10 * (x - 0.0)))]

@testset "Real symmetric matrix" begin
    X = rand(N, N);
    H = 0.5 * (X + X');
	eigs, Ψ = eigen(H);
	for order in 1:1
		h = [rand(N, N) for i = 1:order]
		hs = [0.5 * (x + x') for x in h]
		for f in Funs
			fd_test = frechet_naive(f, eigs, Ψ, hs);
			fd = mat_fun_frechet(f, eigs, Ψ, hs);
			@test isapprox(fd_test, fd)
		end
	end
end

@testset "Complex hermitian matrix" begin
    X = rand(ComplexF64, N, N)
    H = 0.5 * (X + X')
    eigs, Ψ = eigen(H)
    for order in 1:1
        h = [rand(ComplexF64, N, N) for i = 1:order]
        hs = [0.5 * (x + x') for x in h]
        for f in Funs
            fd_test = frechet_naive(f, eigs, Ψ, hs)
            fd = mat_fun_frechet(f, eigs, Ψ, hs)
            @test isapprox(fd_test, fd)
        end
    end
end

end