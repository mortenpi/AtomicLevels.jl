module AtomicLevels

using UnicodeFun
using Formatting
using Parameters
using BlockBandedMatrices
using WignerSymbols
macro hi_str(s)
    parse(HalfInteger, s)
end

if VERSION < v"1.1-DEV"
    isnothing(::Nothing) = true
    isnothing(::Any) = false
end

include("common.jl")
include("parity.jl")
include("orbitals.jl")
include("configurations.jl")
include("excited_configurations.jl")
include("terms.jl")
include("allchoices.jl")
include("jj_terms.jl")
include("couple_terms.jl")
include("csfs.jl")
include("jj2lsj.jl")
include("levels.jl")

module Utils
include("utils/print_states.jl")
end

end # module
