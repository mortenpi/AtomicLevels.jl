@testset "Excited configurations" begin
    @test excited_configurations(c"1s2", os"2[s-d]") ==
        [c"1s2", c"1s 2s", c"2s2", c"2p-2", c"2p- 2p", c"2p2"]
end
