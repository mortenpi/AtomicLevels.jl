using Random

@testset "Orbitals" begin
    @testset "Construction" begin
        @test o"1s" == Orbital(1,0)
        @test o"2p-" == Orbital(2,1,1//2)
        @test o"2p" == Orbital(2,1,3//2)
        @test o"2[1]" == Orbital(2,1,3//2)

        @test_throws ArgumentError AtomicLevels.orbital_from_string("sdkfl")
    end

    @testset "Order" begin
        @test sort(shuffle([o"1s", o"2s", o"2p-", o"2p", o"3s", o"3p-", o"3p"])) ==
            [o"1s", o"2s", o"2p-", o"2p", o"3s", o"3p-", o"3p"]
    end

    @testset "Parity" begin
        @test parity(o"1s") == 1
        @test parity(o"2p") == -1
        @test parity(o"3s") == 1
        @test parity(o"3p") == -1
        @test parity(o"3d") == 1
    end
end
