using TriangularSolve
using Documenter

DocMeta.setdocmeta!(TriangularSolve, :DocTestSetup, :(using TriangularSolve); recursive=true)

makedocs(;
    modules=[TriangularSolve],
    authors="chriselrod <elrodc@gmail.com> and contributors",
    repo="https://github.com/JuliaSIMD/TriangularSolve.jl/blob/{commit}{path}#{line}",
    sitename="TriangularSolve.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://JuliaSIMD.github.io/TriangularSolve.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/JuliaSIMD/TriangularSolve.jl",
)
