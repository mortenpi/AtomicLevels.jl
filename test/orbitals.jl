using Random

@testset "Orbitals" begin
    @testset "Construction" begin
        @test o"1s" == Orbital(1,0)
        @test o"2p-" == Orbital(2,1,1//2)
        @test o"2p" == Orbital(2,1,3//2)
        @test o"2[1]" == Orbital(2,1,3//2)

        @test o"kp" == Orbital(:k,1,3//2)
        @test o"ϵd-" == Orbital(:ϵ,2,3//2)

        @test_throws ArgumentError AtomicLevels.orbital_from_string("sdkfl")

        @test_throws ArgumentError Orbital(0, 0)
        @test_throws ArgumentError Orbital(1, 1)
        @test_throws ArgumentError Orbital(1, 0, 3//2)
    end

    @testset "Order" begin
        @test sort(shuffle([o"1s", o"2s", o"2p-", o"2p", o"3s", o"3p-", o"3p"])) ==
            [o"1s", o"2s", o"2p-", o"2p", o"3s", o"3p-", o"3p"]
        @test sort([o"ls", o"kp", o"2p", o"1s"]) ==
            [o"1s", o"2p", o"kp", o"ls"]
    end

    @testset "Range of orbitals" begin
        @test os"6[s-d] 5[d]" == [o"5d-", o"5d", o"6s", o"6p-", o"6p", o"6d-", o"6d"]
        @test os"k[s-g] l[s-g]" == [o"ks", o"kp-", o"kp", o"kd-", o"kd", o"kf-", o"kf", o"kg-", o"kg",
                                    o"ls", o"lp-", o"lp", o"ld-", o"ld", o"lf-", o"lf", o"lg-", o"lg"]
    end

    @testset "Flip j" begin
        @test AtomicLevels.flip_j(o"1s") == o"1s"
        @test AtomicLevels.flip_j(o"2p-") == o"2p"
        @test AtomicLevels.flip_j(o"2p") == o"2p-"
        @test AtomicLevels.flip_j(o"3d-") == o"3d"
        @test AtomicLevels.flip_j(o"3d") == o"3d-"

        @test AtomicLevels.flip_j(o"kd-") == o"kd"
        @test AtomicLevels.flip_j(o"kd") == o"kd-"
    end

    @testset "Degeneracy" begin
        @test degeneracy(o"1s") == 2
        @test degeneracy(o"2p-") == 2
        @test degeneracy(o"2p") == 4
        @test degeneracy(o"3d-") == 4
        @test degeneracy(o"3d") == 6

        @test degeneracy(o"kp-") == 2
        @test degeneracy(o"kp") == 4

        @test non_rel_degeneracy(o"1s") == 2
        @test non_rel_degeneracy(o"2p-") == 6
        @test non_rel_degeneracy(o"2p") == 6
        @test non_rel_degeneracy(o"3d-") == 10
        @test non_rel_degeneracy(o"3d") == 10

        @test non_rel_degeneracy(o"kp-") == 6
        @test non_rel_degeneracy(o"kp") == 6
    end

    @testset "Parity" begin
        @test parity(o"1s") == 1
        @test parity(o"2p") == -1
        @test parity(o"3s") == 1
        @test parity(o"3p") == -1
        @test parity(o"3d") == 1

        @test parity(o"kp") == -1
    end
end
