@testset "CSFs" begin
    @testset "Construction" begin
        csf = CSF(rc"1s2 2p- 2p", [0, 1//2, 3//2], [0, 1//2, 2])
        @test csf isa CSF{RelativisticOrbital{Int},HalfInteger,HalfInteger}
        @test csf == csf
        @test csf != CSF(rc"1s2 2p- 2p", [0, 1//2, 3//2], [0, 1//2, 1])
    end

    @testset "CSF list generation" begin
        @testset "ls coupling" begin
            csf_1s2 = CSF(c"1s2", [IntermediateTerm(T"1S",0)],[T"1S"])
            @test csfs(c"1s2") == [csf_1s2]
            @test string(csf_1s2) == "1s²(₀¹S|¹S)+"

            csfs_1s_kp = [CSF(c"1s kp",
                              [IntermediateTerm(T"2S",1),IntermediateTerm(T"2Po",1)],
                              [T"2S", T"1Po"]),
                          CSF(c"1s kp",
                              [IntermediateTerm(T"2S",1),IntermediateTerm(T"2Po",1)],
                              [T"2S", T"3Po"])]
            @test csfs(c"1s kp") == csfs_1s_kp
        end
        @testset "jj coupling" begin
            # These are not real tests, they only make sure that the
            # generation does not crash.
            sort(csfs([rc"3s 3p2",rc"3s 3p- 3p"]))
            csfs(rc"[Kr] 5s2 5p-2 5p3 6s")
            csfs(rc"1s kp")

            @test string(CSF(rc"[Kr] 5s2 5p- 5p4 kd", [0//1, 1//2, 0//1, 5//2], [0//1, 1//2, 1//2, 3//1])) ==
                "[Kr]ᶜ 5s²(0|0) 5p⁻(1/2|1/2) 5p⁴(0|1/2) kd(5/2|3)-"
        end
    end
end
