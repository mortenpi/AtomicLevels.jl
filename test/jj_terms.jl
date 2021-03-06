@testset "JJ terms" begin
    @testset "jj coupling of equivalent electrons" begin
        # Table 4.5, Cowan 1981
        foreach([
            (ro"1s",ro"2p-") => [(0,2) => 0, 1 => 1//2],
            (ro"2p",ro"3d-") => [(0,4) => 0, (1,3) => 3//2, 2 => (0,2)],
            (ro"3d",ro"4f-") => [(0,6) => 0, (1,5) => 5//2, (2,4) => (0,2,4),
                                3 => (3//2,5//2,9//2)],
            (ro"4f",ro"5g-") => [(0,8) => 0, (1,7) => 7//2, (2,6) => (0,2,4,6),
                                (3,5) => (3//2,5//2,7//2,9//2,11/2,15//2),
                                4 => (0,2,2,4,4,5,6,8)],
            (ro"5g",ro"6h-") => [(0,10) => 0, (1,9) => 9//2, (2,8) => (0,2,4,6,8),
                                (3,7) => (3//2,5//2,7//2,9//2,9//2,11//2,13//2,15//2,17//2,21//2),
                                (4,6) => (0,0,2,2,3,4,4,4,5,6,6,6,7,8,8,9,10,12),
                                5 => (1//2,3//2,5//2,5//2,7//2,7//2,9//2,9//2,9//2,11//2,11//2,13//2,13//2,15//2,15//2,17//2,17//2,19//2,21//2,25//2)]
        ]) do (orbs,wsj)
            foreach(orbs) do orb
                foreach(wsj) do (ws,j)
                    js = map(HalfInteger, isa(j, Number) ? [j] : j)
                    foreach(ws) do w
                        # Should actually test without unique
                        @test unique(intermediate_terms(orb,w)) == unique(js)
                    end
                end
            end
        end
    end
end
