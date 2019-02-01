using Documenter

makedocs(
    sitename = "AtomicLevels",
    pages = [
        "Home" => "index.md",
    ],
    assets = ["assets/latex.js"],
)

deploydocs(repo = "github.com/JuliaAtoms/AtomicLevels.jl.git")
