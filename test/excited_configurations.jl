@testset "Excited configurations" begin
    @test excited_configurations(c"1s2", os"2[s-p]"...) ==
        [c"1s2", c"1s 2s", c"2s2", c"2p2"]
    @test excited_configurations(c"1s2", os"k[s-p]"...) ==
        [c"1s2", c"1s ks", c"ks2", c"kp2"]

    @test excited_configurations(rc"1s2", ros"2[s-p]"...) ==
        [rc"1s2", rc"1s 2s", rc"2s2", rc"2p-2", rc"2p- 2p", rc"2p2"]
end
