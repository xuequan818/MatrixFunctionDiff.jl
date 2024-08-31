using MatrixFunctionDiff
using Test

@testset "MatrixFunctionDiff.jl" begin
    t0 = time()
    @testset "Fréchet Derivative" begin
        println("##### Testing Fréchet derivative functionality...")
        t = @elapsed include("FrechetTest.jl")
        println("##### done (took $t seconds).")
    end
    println("##### Running all MatrixFunctionDiff tests took $(time() - t0) seconds.")
end
