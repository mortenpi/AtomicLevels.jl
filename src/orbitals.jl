typealias Orbital{I<:Integer,S<:AbstractString} Tuple{I,I,I,S} # (n,ell,occ,*/c/i)
typealias Config Vector{Orbital}

degeneracy(ell::Integer) = 2ell + 1
degeneracy(o::Orbital) = 2*degeneracy(o[2])
filled(o::Orbital) = (o[3] == degeneracy(o))

import Base.parity
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

import Base.isless
isless(o1::Orbital, o2::Orbital) = o1[1] < o2[1] || (o1[1] == o2[1] && o1[2] < o2[2])

noble_gases = Dict("He" => "1s2",
                   "Ne" => "[He] 2s2 2p6",
                   "Ar" => "[Ne] 3s2 3p6",
                   "Kr" => "[Ar] 3d10 4s2 4p6",
                   "Xe" => "[Kr] 4d10 5s2 5p6",
                   "Rn" => "[Xe] 4f14 5d10 6s2 6p6")

function ref_set_list(ref_set::AbstractString)
    m = match(r"\[([a-zA-Z]{2})\]([ci*]{0,1})(.*)", ref_set)
    core = []
    if m != nothing
        if m[1] ∉ keys(noble_gases)
            error("Unknown core, $m[1]")
        end
        s = m[2] != "" ? m[2] : "c"
        core_ref_set = join(map(o -> "$o$s", split(noble_gases[m[1]])), " ")
        core = sort(ref_set_list(core_ref_set))
        ref_set = strip(m[3])
    end
    orbs = map(split(ref_set)) do orb
        m = match(r"([0-9]+)([a-z])([0-9]*)([ci*]{0,1})", orb)
        n = parse(Int, m[1])
        ell_i = searchindex(ells, m[2]) - 1
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

import Base.print, Base.show, Base.writemime, Base.string

function print(io::IO, o::Orbital)
    print(io, o[1], ells[o[2]+1])
    o[3]>1 && print(io,superscript(o[3]))
    if o[4] == "c"
        print(io, "ᶜ")
    elseif o[4] == "i"
        print(io, "ⁱ")
    end
end
show(io::IO, o::Orbital) = print(io, o)

function writemime(io::IO, ::MIME"text/latex", o::Orbital, wrap = true)
    m =  o[4] == "*" ? "" : o[4]
    wrap && print(io, "\$")
    print(io, "$(o[1])\\mathrm{$(ells[o[2]+1])}^{$(o[3])$m}")
    wrap && print(io, "\$")
end

function string(o::Orbital)
    m =  o[4] == "*" ? "" : o[4]
    "$(o[1])$(ells[o[2]+1])$(o[3] > 1 ? o[3] : "")$m"
end

function print(io::IO, c::Config)
    for o in c
        print(io, o)
    end
end
show(io::IO, c::Config) = print(io, c)

function writemime(io::IO, m::MIME"text/latex", c::Config, wrap = true)
    wrap && print(io, "\$")
    for o in c
        writemime(io, m, o, false)
    end
    wrap && print(io, "\$")
end

string(c::Config) = join([string(o) for o in c], "_")

export Orbital, Config, degeneracy, filled, parity, open, closed, fill,
ref_set_list, @c_str, print, show, string
