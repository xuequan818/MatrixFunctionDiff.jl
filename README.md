[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://xuequan818.github.io/MatrixFunctionDiff.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://xuequan818.github.io/MatrixFunctionDiff.jl/dev/)
[![Build Status](https://github.com/xuequan818/MatrixFunctionDiff.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/xuequan818/MatrixFunctionDiff.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/xuequan818/MatrixFunctionDiff.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/xuequan818/MatrixFunctionDiff.jl)

# MatrixFunctionDiff.jl

A julia package for computing **Fréchet derivatives** of scalar functions with matricies as variables. The higher order Fréchet derivatives are first written in a formula similar to the [Daleskii-Krein theorem](https://www.ams.org/books/trans2/047/), and then computed by the [divided difference](https://github.com/xuequan818/DividedDifferences.jl).

#### Usage

```julia
julia> using MatrixFunctionDiff
julia> using DividedDifferences

julia> f = heaviside; # Heaviside step function

julia> order = 2; # 2nd order Fréchet derivative

julia> H = 0.1 .* reshape(collect(-4:11), 4, 4) + diagm(collect(-2:1))
4×4 Matrix{Float64}:
 -2.4   0.0  0.4  0.8
 -0.3  -0.9  0.5  0.9
 -0.2   0.2  0.6  1.0
 -0.1   0.3  0.7  2.1

julia> h = map(x -> Array(reshape(x, 4, 4)) .+ diagm(0.1*ones(4)), [-0.2:0.05:0.55, -1:0.2:2]);

julia> mat_fun_frechet(f, H, h)
4×4 Matrix{Float64}:
 0.00483228  -0.00766712  -0.0402159   -0.0588928
 0.0098505   -0.0152527   -0.0825338   -0.0813227
 0.0412279   -0.0186076    0.00722674   0.00971879
 0.0189937   -0.0209538    0.00616146   0.00319365
```

MatrixFunctionDiff can also compute the Fréchet derivatives by the eigenpairs of `H`:

```julia
using LinearAlgebra

julia> eigs, Ψ = eigen(H);

julia> mat_fun_frechet(f, eigs, Ψ, h)
4×4 Matrix{Float64}:
 0.00483228  -0.00766712  -0.0402159   -0.0588928
 0.0098505   -0.0152527   -0.0825338   -0.0813227
 0.0412279   -0.0186076    0.00722674   0.00971879
 0.0189937   -0.0209538    0.00616146   0.00319365
```

For more details, please see the [documentation](https://xuequan818.github.io/MatrixFunctionDiff.jl/dev/).