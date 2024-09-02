using MatrixFunctionDiff
using Documenter

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
        "Background" => "background.md",
        "Examples" => "examples.md",
        "API" => "api.md"
    ],
)

deploydocs(;
    repo="github.com/xuequan818/MatrixFunctionDiff.jl",
    devbranch="master",
)
