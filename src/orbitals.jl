using UnicodeFun
using Printf

Orbital{I<:Integer,S<:AbstractString} = Tuple{I,I,I,S} # (n,ell,occ,*/c/i)
Config = Vector{Orbital}

degeneracy(ell::Integer) = 2ell + 1
degeneracy(o::Orbital) = 2*degeneracy(o[2])
filled(o::Orbital) = (o[3] == degeneracy(o))

parity(orbital::Orbital) = (-1)^(orbital[3]*orbital[2])
parity(config::Config) = mapreduce(o -> parity(o), *, config)

import Base: open
open(config::Config) =
    filter(config) do o
        o[4] != "c"
    end
closed(config::Config) =
    filter(config) do o
        o[4] == "c"
    end
import Base: fill
fill(c::Config) = [Orbital((o[1], o[2], degeneracy(o), o[4])) for o in c]

nelc(c::Config) = sum([o[3] for o in c])

import Base.isless
isless(o1::Orbital, o2::Orbital) = o1[1] < o2[1] || (o1[1] == o2[1] && o1[2] < o2[2])

noble_gases = Dict("He" => "1s2",
                   "Ne" => "[He] 2s2 2p6",
                   "Ar" => "[Ne] 3s2 3p6",
                   "Kr" => "[Ar] 3d10 4s2 4p6",
                   "Xe" => "[Kr] 4d10 5s2 5p6",
                   "Rn" => "[Xe] 4f14 5d10 6s2 6p6")
function core_ref_set(symbol::AbstractString, s::AbstractString="c")
    if symbol ∉ keys(noble_gases)
        error("Unknown core, $(symbol)")
    end
    s == "" && (s = "c")
    sort(ref_set_list(join(map(o -> "$o$s", split(noble_gases[symbol])), " ")))
end

function ref_set_list(ref_set::AbstractString)
    m = match(r"\[([a-zA-Z]{2})\]([ci*]{0,1})(.*)", ref_set)
    core = []
    if m != nothing
        core = core_ref_set(m[1], m[2])
        ref_set = strip(m[3])
    end
    orbs = map(split(ref_set)) do orb
        m = match(r"([0-9]+)([a-z])([0-9]*)([ci*]{0,1})", orb)
        n = parse(Int, m[1])
        ell_i = something(findfirst(isequal(m[2][1]), ells), 0) - 1
        ell_i >= n && error("Invalid orbital $(m[1])$(m[2])")
        Orbital((n,
                 ell_i,
                 m[3] != "" ? parse(Int, m[3]) : 1,
                 m[4] != "" ? m[4] : "*"))
    end
    Config([core;orbs])
end

macro c_str(s)
    ref_set_list(s)
end

import Base.show, Base.string

function latex(o::Orbital)
    s = if o[4] == "c"
       "^c"
    elseif o[4] == "i"
        "^i"
    else
        ""
    end
    occ = o[3] > 1 ? "^{$(o[3])}" : ""
    "$(o[1])$(ells[o[2]+1])$(occ)$(s)"
end
show(io::IO, o::Orbital) = print(io, to_latex(latex(o)))

function show(io::IO, ::MIME"text/latex", o::Orbital, wrap = true)
    wrap && print(io, "\$")
    print(io, "\\mathrm{$(latex(o))}")
    wrap && print(io, "\$")
end

function string(o::Orbital)
    m =  o[4] == "*" ? "" : o[4]
    "$(o[1])$(ells[o[2]+1])$(o[3] > 1 ? o[3] : "")$m"
end

function strip_core(c::Config)
    core_str = ""
    c1 = sort(copy(c))
    s = c1[1][4]

    for core in reverse(["He", "Ne", "Ar", "Kr", "Xe", "Rn"])
        core_rs = core_ref_set(core, s)
        stripped = all(map(o -> o ∈ c, core_rs)) ? sort(filter(o -> o ∉ core_rs, c)) : c
        if length(stripped) < length(c1) && length(stripped) > 0
            core_str = core
            c1 = stripped
            break
        end
    end

    core_str, occursin(s, "ci") ? s : "", c1
end

function latex(c::Config)
    core,cs,c = strip_core(c)
    core_str = if core == ""
        ""
    elseif cs == ""
        "[$(core)] "
    else
        @sprintf("[%s]^{%s} ", core, cs)
    end

    @sprintf("%s%s", core_str,
             join(map(latex, c), " "))
end
show(io::IO, c::Config) = print(io, to_latex(latex(c)))

function show(io::IO, m::MIME"text/latex", c::Config, wrap = true)
    wrap && print(io, "\$")
    print(io, "\\mathrm{$(latex(c))}")
    wrap && print(io, "\$")
end

function string(c::Config)
    core,cs,c = strip_core(c)
    core_str = if core == ""
        ""
    elseif cs == ""
        "$(core)_"
    else
        @sprintf("%s_%s_",
                 core == "" ? "" : "$(core)",
                 cs == "" ? "" : cs)
    end
    @sprintf("%s%s", core_str,
             join(map(string, c), "_"))
end

export Orbital, Config, degeneracy, filled, parity, open, closed, fill, nelc,
ref_set_list, @c_str, show, string
