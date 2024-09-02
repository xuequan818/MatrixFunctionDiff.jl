module MatrixFunctionDiff

using Combinatorics
using DividedDifferences
using LinearAlgebra
using PermutationSymmetricTensors
import PermutationSymmetricTensors: getindex

export mat_fun_frechet
include("frechet.jl")

end # module MatrixFunctionDiff
