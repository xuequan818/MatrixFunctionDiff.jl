using MatrixFunctionDiff
using Documenter

DocMeta.setdocmeta!(MatrixFunctionDiff, :DocTestSetup, :(using MatrixFunctionDiff); recursive=true)

makedocs(;
    modules=[MatrixFunctionDiff],
    authors="Xue Quan <xuequan818@gmail.com>",
    sitename="MatrixFunctionDiff.jl",
    format=Documenter.HTML(;
        canonical="https://xuequan818.github.io/MatrixFunctionDiff.jl",
        edit_link="master",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Theory" => "theory.md"
    ],
)

deploydocs(;
    repo="github.com/xuequan818/MatrixFunctionDiff.jl",
    devbranch="main",
)
