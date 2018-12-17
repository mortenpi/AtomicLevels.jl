using Random

@testset "Orbitals" begin
    @testset "Construction" begin
        @test o"1s" == Orbital(1,0)
        @test o"2p" == Orbital(2,1)
        @test o"2[1]" == Orbital(2,1)

        @test ro"1s" == RelativisticOrbital(1,0)
        @test ro"2p-" == RelativisticOrbital(2,1,1//2)
        @test ro"2p" == RelativisticOrbital(2,1,3//2)
        @test ro"2[1]" == RelativisticOrbital(2,1,3//2)

        @test o"kp" == Orbital(:k,1)
        @test o"ϵd" == Orbital(:ϵ,2)

        @test ro"kp" == RelativisticOrbital(:k,1,3//2)
        @test ro"ϵd-" == RelativisticOrbital(:ϵ,2,3//2)

        @test_throws ArgumentError AtomicLevels.orbital_from_string(Orbital, "2p-")
        @test_throws ArgumentError AtomicLevels.orbital_from_string(Orbital, "sdkfl")

        @test_throws ArgumentError Orbital(0, 0)
        @test_throws ArgumentError Orbital(1, 1)
        @test_throws ArgumentError RelativisticOrbital(0, 0)
        @test_throws ArgumentError RelativisticOrbital(1, 1)
        @test_throws ArgumentError RelativisticOrbital(1, 0, 3//2)
    end

    @testset "Order" begin
        @test sort(shuffle([o"1s", o"2s", o"2p", o"3s", o"3p"])) ==
            [o"1s", o"2s", o"2p", o"3s", o"3p"]
        @test sort([o"ls", o"kp", o"2p", o"1s"]) ==
            [o"1s", o"2p", o"kp", o"ls"]

        @test sort(shuffle([ro"1s", ro"2s", ro"2p-", ro"2p", ro"3s", ro"3p-", ro"3p"])) ==
            [ro"1s", ro"2s", ro"2p-", ro"2p", ro"3s", ro"3p-", ro"3p"]
        @test sort([ro"ls", ro"kp", ro"2p", ro"1s"]) ==
            [ro"1s", ro"2p", ro"kp", ro"ls"]
    end

    @testset "Range of orbitals" begin
        @test os"6[s-d] 5[d]" == [o"5d", o"6s", o"6p", o"6d"]
        @test os"k[s-g] l[s-g]" == [o"ks", o"kp", o"kd", o"kf", o"kg",
                                    o"ls", o"lp", o"ld", o"lf", o"lg"]

        @test ros"6[s-d] 5[d]" == [ro"5d-", ro"5d", ro"6s", ro"6p-", ro"6p", ro"6d-", ro"6d"]
        @test ros"k[s-g] l[s-g]" == [ro"ks", ro"kp-", ro"kp", ro"kd-", ro"kd", ro"kf-", ro"kf", ro"kg-", ro"kg",
                                     ro"ls", ro"lp-", ro"lp", ro"ld-", ro"ld", ro"lf-", ro"lf", ro"lg-", ro"lg"]
    end

    @testset "Flip j" begin
        @test AtomicLevels.flip_j(ro"1s") == ro"1s"
        @test AtomicLevels.flip_j(ro"2p-") == ro"2p"
        @test AtomicLevels.flip_j(ro"2p") == ro"2p-"
        @test AtomicLevels.flip_j(ro"3d-") == ro"3d"
        @test AtomicLevels.flip_j(ro"3d") == ro"3d-"

        @test AtomicLevels.flip_j(ro"kd-") == ro"kd"
        @test AtomicLevels.flip_j(ro"kd") == ro"kd-"
    end

    @testset "Degeneracy" begin
        @test degeneracy(o"1s") == 2
        @test degeneracy(o"2p") == 6
        @test degeneracy(o"3d") == 10

        @test degeneracy(o"kp") == 6

        @test degeneracy(ro"1s") == 2
        @test degeneracy(ro"2p-") == 2
        @test degeneracy(ro"2p") == 4
        @test degeneracy(ro"3d-") == 4
        @test degeneracy(ro"3d") == 6

        @test degeneracy(ro"kp-") == 2
        @test degeneracy(ro"kp") == 4
    end

    @testset "Parity" begin
        @test iseven(parity(o"1s"))
        @test isodd(parity(o"2p"))
        @test iseven(parity(o"3s"))
        @test isodd(parity(o"3p"))
        @test iseven(parity(o"3d"))

        @test isodd(parity(o"kp"))

        @test iseven(parity(ro"1s"))
        @test isodd(parity(ro"2p"))
        @test iseven(parity(ro"3s"))
        @test isodd(parity(ro"3p"))
        @test iseven(parity(ro"3d"))

        @test isodd(parity(ro"kp"))
    end
end
