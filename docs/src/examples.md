# Examples

The [`mat_fun_frechet`](@ref) function computes the higher order Fréchet derivatives ${{\rm d}}^{n}f(H)hs[1]\cdots hs[n]$. It requires that the matrix `H` be a real symmetric or complex hermitian, the function `f` be a scalar function, and `hs` be a vector consisting of real symmetric or complex hermitian matrices. The order of the derivative is equal to the length of `hs`.

```julia
julia> using MatrixFunctionDiff

julia> f(x) = 1.0 / (1.0 + exp(100 * (x - 0.1))); # Fermi Dirac function

julia> order = 3; # 3rd order Fréchet derivative

julia> N = 10; # dimension of matrix

julia> X = rand(N, N);

julia> H = 0.5 * (X + X');

julia> h = [rand(N, N) for i = 1:order];

julia> hs = [0.5 * (x + x') for x in h];

julia> mat_fun_frechet(f, H, hs);
```

MatrixFunctionDiff can also compute the Fréchet derivatives by the eigenpairs of `H`:
```julia
using LinearAlgebra

julia> eigs, Ψ = eigen(H);

julia> mat_fun_frechet(f, eigs, Ψ, hs);
```

