using AtomicLevels
using UnicodeFun
using Test

@testset "Terms" begin
    @testset "Construction" begin
        @test T"1S" == Term(0, 0, 1)
        @test T"1Se" == Term(0, 0, 1)
        @test T"1So" == Term(0, 0, -1)
        @test T"2So" == Term(0, 1//2, -1)
        @test T"4P" == Term(1, 3//2, 1)
        @test T"3D" == Term(2, 1, 1)
        @test T"3Do" == Term(2, 1, -1)
        @test T"1[54]" == Term(54, 0, 1)
        @test T"1[3/2]" == Term(3//2, 0, 1)
        @test T"2[3/2]o" == Term(3//2, 1//2, -1)
        
        @test_throws ArgumentError AtomicLevels.term_string("1[4/3]")
        @test_throws ArgumentError AtomicLevels.term_string("1[43/]")
        @test_throws ArgumentError AtomicLevels.term_string("1[/43/]")
        @test_throws ArgumentError AtomicLevels.term_string("1[/43]")
        @test_throws ArgumentError AtomicLevels.term_string("P")
        @test_throws ArgumentError AtomicLevels.term_string("asdf")
    end

    @testset "Properties" begin
        @test multiplicity(T"1S") == 1
        @test multiplicity(T"1So") == 1
        @test multiplicity(T"2S") == 2
        @test multiplicity(T"2So") == 2
        @test multiplicity(T"3S") == 3
        @test multiplicity(T"3So") == 3

        @test weight(T"1S") == 1
        @test weight(T"1P") == 3
        @test weight(T"2S") == 2
        @test weight(T"2P") == 6

        @test T"1S" > T"1So"
        @test T"1So" < T"1S"
        @test T"1S" < T"2S"
        @test T"1S" < T"1P"
        @test T"1S" < T"3S"
        @test T"1P" < T"3S"
    end

    @testset "Pretty printing" begin
        map([T"1S" => "¹S",
             T"2So" => "²Sᵒ",
             T"4[3/2]" => "⁴[3/2]"]) do (t,s)
            @test "$(t)" == s
        end
    end

    @testset "Orbital terms" begin
        function test_single_orbital_terms(orb::Orbital{I,R}, occ::I, ts::Vector{Term{I,R,LT}}) where {I<:Integer,R<:Rational{I},LT}
            cts = sort(terms(orb, occ))
            ts = sort(ts)
            ccts = copy(cts)
            for t in cts
                if t in ts
                    @test count(e -> e == t, ccts) == count(e -> e == t, ts)
                    ccts = filter!(e -> e != t, ccts)
                    ts = filter!(e -> e != t, ts)
                end
            end
            if length(ts) != 0
                println("fail, terms: ", join(string.(cts), ", "))
                println("missing: ", join(string.(ts), ", "))
                println("should not be there: ", join(string.(ccts), ", "))
                println("==========================================================")
            end
            @test length(ts) == 0
        end

        function test_single_orbital_terms(orb::Orbital{I,R}, occs::Tuple{I,I}, ts::Vector{Term{I,R,LT}}) where {I<:Integer,R<:Rational{I},LT}
            map(occs) do occ
                test_single_orbital_terms(orb, occ, ts)
            end
        end

        function get_orbital(o::AbstractString)
            n = something(findfirst(isequal(o[1]), AtomicLevels.spectroscopic), 0)
            ℓ = n - 1
            orb = Orbital(n, ℓ)
            g = 2(2ℓ + 1)
            occ = length(o)>1 ? parse(Int, o[2:end]) : 1
            if g != occ && g/2 != occ
                occ = (occ, g-occ)
            end
            orb, occ
        end

        function test_orbital_terms(o_ts::Pair{<:ST,<:ST}) where {ST<:AbstractString}
            orb,occ = get_orbital(o_ts[1])
            m = match(r"([0-9]+)([A-Z])", o_ts[2])
            L = something(findfirst(isequal(lowercase(m[2])[1]), AtomicLevels.spectroscopic), 0)-1
            S = (parse(Int, m[1])-1)//2
            test_single_orbital_terms(orb, occ, [Term(L, S, parity(orb)^occ[1])])
        end

        function test_orbital_terms(o_ts::Pair{<:ST,<:Vector{ST}}) where {ST<:AbstractString}
            orb,occ = get_orbital(o_ts[1])
            
            p = parity(orb)^occ[1]
            ts = o_ts[2]

            p1 = r"([0-9]+)\(((?:[A-Z][0-9]*)+)\)"
            p2 = r"([0-9]+)([A-Z])"
            toL = s -> something(findfirst(isequal(lowercase(s)[1]), AtomicLevels.spectroscopic), 0) - 1
            ts = map(ts) do t
                m = match(p1, t)
                if m != nothing
                    S = (parse(Int, m[1]) - 1)//2
                    map(eachmatch(r"([A-Z])([0-9]*)", m[2])) do mmm
                        mm = mmm.match
                        [Term(toL(mm[1]),S,p)
                         for j in 1:(length(mm)>1 ? parse(Int, mm[2:end]) : 1)]
                    end
                else
                    m = match(p2, t)
                    S = (parse(Int, m[1]) - 1)//2
                    Term(toL(m[2]), S, p)
                end
            end
            test_single_orbital_terms(orb, occ, vcat(vcat(ts...)...))
        end

        # Table taken from Cowan, p. 110
        # Numbers following term symbols indicate the amount of times
        # different terms with the same (L,S) occur.
        test_orbital_terms("s" => "2S")
        for o in ["s2", "p6", "d10", "f14"]
            test_orbital_terms(o => "1S")
        end
        test_orbital_terms("p" => "2P")
        test_orbital_terms("p2" => ["1(SD)", "3P"])
        test_orbital_terms("p3" => ["2(PD)", "4S"])
        test_orbital_terms("d" => "2D")
        test_orbital_terms("d2" => ["1(SDG)", "3(PF)"])
        test_orbital_terms("d3" => ["2(PD2FGH)", "4(PF)"])
        test_orbital_terms("d4" => ["1(S2D2FG2I)", "3(P2DF2GH)", "5D"])
        test_orbital_terms("d5" => ["2(SPD3F2G2HI)", "4(PDFG)", "6S"])
        test_orbital_terms("f" => "2F")
        test_orbital_terms("f2" => ["1(SDGI)", "3(PFH)"])
        test_orbital_terms("f3" => ["2(PD2F2G2H2IKL)", "4(SDFGI)"])
        test_orbital_terms("f4" => ["1(S2D4FG4H2I3KL2N)", "3(P3D2F4G3H4I2K2LM)", "5(SDFGI)"])
        test_orbital_terms("f5" => ["2(P4D5F7G6H7I5K5L3M2NO)", "4(SP2D3F4G4H3I3K2LM)", "6(PFH)"])
        test_orbital_terms("f6" => ["1(S4PD6F4G8H4I7K3L4M2N2Q)", "3(P6D5F9G7H9I6K6L3M3NO)", "5(SPD3F2G3H2I2KL)", "7F"])
        test_orbital_terms("f7" => ["2(S2P5D7F10G10H9I9K7L5M4N2OQ)", "4(S2P2D6F5G7H5I5K3L3MN)", "6(PDFGHI)", "8S"])
        # # Data below are from Xu2006
        # test_orbital_terms("g9" => [2(S8 P19 D35 F40 G52 H54 I56 K53 L53 M44 N40 O32 Q26 R19 T15 U9 V7 W4 X2 YZ) 4(S6 P16 D24 F34 G38 H40 I42 K39 L35 M32 N26 O20 Q16 R11 T7 U5 V3 WX) 6(S3P3D9F8G12H10I12K9L9M6N6O3Q3RT) 8(PDFGHIKL) 10(S)])
        # test_orbital_terms("h11" => [2(S36 P107 D173 F233 G283 H325 I353 K370 L376 M371 N357 O335 Q307 R275 T241 U207 V173 W142 X114 Y88 Z68 2150 2236 2325 2417 2511 267 274 282 29 30) 4 (S37 P89 D157 F199 G253 H277 I309 K313 L323 M308 N300 O271 Q251 R216 T190 U155 V131 W101 X81 Y59 Z45 2130 2222 2313 249 255 263 27 28) 6 (S12 P35 D55 F76 G90 H101 I109 K111 L109 M105 N97 O87 Q77 R65 T53 U43 V33 W24 X18 Y12 Z8 215 223 2324) 8(S4 P4 D12 F11 G17 H15 I19 K16 L18 M14 N14 O10 Q10 R6 T6 U3 V3WX) 10(PDFGHIKLMN) 12S])

        @test_throws ArgumentError terms(o"2p", 7)
    end
        
    @testset "Coupling" begin
        # Coupling two 2S terms should yields singlets and triplets of S,P,D
        @test couple_terms(Term(1,1//2,1), Term(1,1//2,1)) ==
            sort([Term(2,0,1), Term(2,1,1), Term(0,0,1), Term(1,1,1), Term(1,0,1), Term(0,1,1)])    
        
        @test couple_terms(Term(0,1//2,1), Term(1,1//2,-1)) ==
            sort([Term(1,0,-1), Term(1,1,-1)])
        
        function test_coupling(o1, o2)
            c1 = couple_terms(terms(o1), terms(o2))
            r = o1 + o2
            c2 = terms(r)
            @test c1 == c2
        end
        
        test_coupling(c"1s", c"2p")
        test_coupling(c"2s", c"2p")
        test_coupling(c"2p", c"3p")
        test_coupling(c"2p", c"3d")

        @test terms(c"1s") == [T"2S"]
        @test terms(c"1s2") == [T"1S"]
        @test terms(c"2p") == [T"2Po"]
        @test terms(c"[He]") == [T"1S"]
        @test terms(c"[Ne]") == [T"1S"]
        @test terms(c"[Ar]") == [T"1S"]
        @test terms(c"[Kr]") == [T"1S"]
        @test terms(c"[Xe]") == [T"1S"]
        @test terms(c"[Rn]") == [T"1S"]

        @test terms(c"1s ks") == [T"1S", T"3S"]
        @test terms(c"1s kp") == [T"1Po", T"3Po"]
    end
end
