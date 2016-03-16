typealias Orbital{I<:Integer,S<:AbstractString} Tuple{I,I,I,S} # (n,ell,occ,*/c/i)
typealias Config Vector{Orbital}

degeneracy(o::Orbital) = 2(2o[2] + 1)

import Base.parity
parity(orbital::Orbital) = (-1)^(orbital[3]*orbital[2])
parity(config::Config) = mapreduce(o -> parity(o), *, config)

function ref_set_list(ref_set::AbstractString)
    orbs = map(split(ref_set)) do orb
        m = match(r"([0-9]+)([a-z])([0-9]*)([ci*]{0,1})", orb)
        Orbital((parse(Int, m[1]),
                 searchindex(ells, m[2]) - 1,
                 m[3] != "" ? parse(Int, m[3]) : 1,
                 m[4] != "" ? m[4] : "*"))
    end
    Config(orbs)
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
    "$(o[1])$(ells[o[2]+1])$(o[3])$m"
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

export Orbital, Config, degeneracy, parity, @c_str, print, show, string
