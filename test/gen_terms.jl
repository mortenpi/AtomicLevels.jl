using AtomicLevels
using Base.Test

import AtomicLevels: ells, ref_set_list

function test_config{A<:Number,B<:Number}(c::Config, ts::Vector{Tuple{A,B}})
    cts = sort(terms(c))
    p = parity(c)
    ts = sort([Term(t...,p) for t in ts])
    ccts = copy(cts)
    for t in cts
        if t in ts
            ccts = filter!(e -> e != t, ccts)
            ts = filter!(e -> e != t, ts)
        end
    end
    if length(ts) == 0
        println(c, ", success, terms: $cts")
    else
        println(c, ", fail, terms: $cts, missing: ", ts, ", should not be there: ", ccts)
    end
    @test length(ts) == 0
end

function test_config{A<:Number,B<:Number}(cs::Tuple{Config,Config}, ts::Vector{Tuple{A,B}})
    map(cs) do c
        test_config(c, ts)
    end
end

function get_config(c::AbstractString)
    n = searchindex(ells, c[1])
    conf = ref_set_list("$n$c")
    ell = n - 1
    g = 2(2ell + 1)
    N = length(c)>1 ? parse(Int, c[2:end]) : 1
    if g != N && g/2 != N
        conf = (conf,ref_set_list("$n$(c[1])$(g-N)"))
    end
    conf
end

function test_config{S<:AbstractString}(c::S, t::S)
    conf = get_config(c)
    m = match(r"([0-9]+)([A-Z])", t)
    test_config(conf, [(searchindex(ells, lowercase(m[2]))-1, (parse(Int, m[1])-1)//2)])
end

function test_config{S<:AbstractString}(c::S, ts::Vector{S})
    conf = get_config(c)

    p1 = r"([0-9]+)\(((?:[A-Z][0-9]*)+)\)"
    p2 = r"([0-9]+)([A-Z])"
    toL = s -> searchindex(ells, lowercase(s)) - 1
    ts = map(ts) do t
        m = match(p1, t)
        if m != nothing
            S = (parse(Int, m[1]) - 1)//2
            m = matchall(r"([A-Z])([0-9]*)", m[2])
            [Tuple{Int,Rational{Int}}((toL(mm[1]),S)) for mm in m]
        else
            m = match(p2, t)
            S = (parse(Int, m[1]) - 1)//2
            Tuple{Int,Rational{Int}}((toL(m[2]), S))
        end
    end
    test_config(conf, vcat(ts...))
end

# Table taken from Cowan, p. 110
# Numbers following term symbols indicate the amount of times
# different terms with the same (L,S) occur.
configs = [("s", "2S")
           ("s2", "1S")
           ("p", "2P")
           ("p2", ["1(SD)", "3P"])
           ("p3", ["2(PD)", "4S"])
           ("d", "2D")
           ("d2", ["1(SDG)", "3(PF)"])
           ("d3", ["2(PD2FGH)", "4(PF)"])
           ("d4", ["1(S2D2FG2I)", "3(P2DF2GH)", "5D"])
           ("d5", ["2(SPD2F2G2HI)", "4(PDFG)", "6S"])
           ("f", "2F")
           ("f2", ["1(SDGI)", "3(PFH)"])
           ("f3", ["2(PD2F2G2H2IKL)", "4(SDFGI)"])
           ("f4", ["1(S2D4FG4H2I3KL2N)", "3(P3D2F4G3H4I2K2LM)", "5(SDFGI)"])
           ("f5", ["2(P4D5F7G5H7I5K5L3M2NO)", "4(SP2D3F4G4H3I3K2LM)", "6(PFH)"])
           ("f6", ["1(S4PD6F4G8H4I7K3L4M2N2Q)", "3(P6D5F9G7H9I6K6L4M4NO)", "5(SPD3F2G3H2I2KL)", "7F"])
           ("f7", ["2(S2P5D7F10G10H9I9K7L5M4N2OQ)", "4(S2P2D6F5G7H5I5K3L3MN)", "6(PDFGHI)", "8S"])]

for c in configs
    test_config(c...)
end
