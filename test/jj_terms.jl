@testset "JJ terms" begin
    @testset "jj coupling of equivalent electrons" begin
        # Table 4.5, Cowan 1981
        foreach([
            (o"1s",o"2p-") => [(0,2) => 0, 1 => 1//2],
            (o"2p",o"3d-") => [(0,4) => 0, (1,3) => 3//2, 2 => (0,2)],
            (o"3d",o"4f-") => [(0,6) => 0, (1,5) => 5//2, (2,4) => (0,2,4),
                              3 => (3//2,5//2,9//2)],
            (o"4f",o"5g-") => [(0,8) => 0, (1,7) => 7//2, (2,6) => (0,2,4,6),
                              (3,5) => (3//2,5//2,7//2,9//2,11/2,15//2),
                              4 => (0,2,2,4,4,5,6,8)],
            (o"5g",o"6h-") => [(0,10) => 0, (1,9) => 9//2, (2,8) => (0,2,4,6,8),
                              (3,7) => (3//2,5//2,7//2,9//2,9//2,11//2,13//2,15//2,17//2,21//2),
                              (4,6) => (0,0,2,2,3,4,4,4,5,6,6,6,7,8,8,9,10,12),
                              5 => (1//2,3//2,5//2,5//2,7//2,7//2,9//2,9//2,9//2,11//2,11//2,13//2,13//2,15//2,15//2,17//2,17//2,19//2,21//2,25//2)]
        ]) do (orbs,wsj)
            foreach(orbs) do orb
                foreach(wsj) do (ws,j)
                    j isa Number && (j = j:j)
                    foreach(ws) do w
                        # Should actually test without unique
                        @test unique(jj_terms(orb,w)) == unique(j)
                    end
                end
            end
        end
    end

    @testset "Coupling of jj terms" begin
        @test couple_terms(0, 1) == 1:1
        @test couple_terms(1//2,1//2) == 0:1
        @test couple_terms(1, 1) == 0:2
        @test couple_terms(1//1, 3//2) == 1//2:5//2
    end
end
