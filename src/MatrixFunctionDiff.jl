module MatrixFunctionDiff

using Combinatorics
using DividedDifferences
using LinearAlgebra
using PermutationSymmetricTensors

export mat_fun_frechet
include("symtensor.jl")
include("frechet.jl")

end # module MatrixFunctionDiff
