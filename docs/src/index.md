# MatrixFunctionDiff.jl

A julia package for computing **Fréchet derivatives** of scalar functions with matricies as variables. The higher order Fréchet derivatives are first written in a formula similar to the [Daleskii-Krein theorem](https://www.ams.org/books/trans2/047/), and then computed by the [divided difference](https://github.com/xuequan818/DividedDifferences.jl).


---
## Installation
MatrixFunctionDiff.jl is a registered package, so it can be installed by running

```julia
julia> Pkg.add("MatrixFunctionDiff")
```

## Related packages
- [DividedDifferences.jl](https://github.com/xuequan818/DividedDifferences.jl): Divided difference for Julia
