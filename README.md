[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://xuequan818.github.io/MatrixFunctionDiff.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://xuequan818.github.io/MatrixFunctionDiff.jl/dev/)
[![Build Status](https://github.com/xuequan818/MatrixFunctionDiff.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/xuequan818/MatrixFunctionDiff.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/xuequan818/MatrixFunctionDiff.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/xuequan818/MatrixFunctionDiff.jl)

# MatrixFunctionDiff.jl

A julia package for computing **Fréchet derivatives** of scalar functions with matricies as variables. The higher order Fréchet derivatives are first written in a formula similar to the [Daleskii-Krein theorem](https://www.ams.org/books/trans2/047/), and then computed by the [divided difference](https://github.com/xuequan818/DividedDifferences.jl).

#### Usage

```julia
julia> using MatrixFunctionDiff
julia> using DividedDifferences

julia> f = DividedDifferences.heaviside; # Heaviside step function
julia> order = 2; # 2nd order Fréchet derivative

```

For more details, please see the [documentation](https://).