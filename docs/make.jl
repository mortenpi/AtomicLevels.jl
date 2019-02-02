using Documenter
using AtomicLevels

makedocs(
    modules = [AtomicLevels],
    sitename = "AtomicLevels",
    pages = [
        "Home" => "index.md",
        "Orbitals" => "orbitals.md",
        "Configurations" => "configurations.md",
        "CSFs" => "csfs.md",
        "Other utilities" => "utilities.md",
        "Internals" => "internals.md",
    ],
)

deploydocs(repo = "github.com/JuliaAtoms/AtomicLevels.jl.git")
