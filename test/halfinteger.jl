# Largely borrowed from WignerSymbols.jl (https://github.com/Jutho/WignerSymbols.jl)
# Copyright (c) 2017: Jutho Haegeman; MIT "Expat" License

@testset "HalfInteger" begin
    @testset "HalfInteger type" begin
        @test convert(HalfInteger, 2) == HalfInteger(twoX = 4)
        @test convert(HalfInteger, 1//2) == HalfInteger(twoX = 1)
        @test convert(HalfInteger, 1.5) == HalfInteger(twoX = 3)
        @test_throws InexactError convert(HalfInteger, 1//3)
        @test_throws InexactError convert(HalfInteger, 0.6)
        @test convert(HalfInteger, 2) == 2
        @test convert(HalfInteger, 1//2) == 1//2
        @test convert(HalfInteger, 1.5) == 1.5
        @test_throws InexactError convert(Integer, HalfInteger(1//2))
        a = HalfInteger(2)
        b = HalfInteger(3//2)
        @test a + b == 2 + 3//2
        @test a - b == 2 - 3//2
        @test zero(a) == 0
        @test one(a) == 1
        @test a > b
        @test b < a
        @test b <= a
        @test a >= b
        @test a == a
        @test a != b

        @test 2 * a == HalfInteger(4)
        @test (-1) * b == HalfInteger(-3//2)

        @test HalfInteger(0) == HalfInteger(twoX = 0)
        @test HalfInteger(1) == HalfInteger(twoX = 2)
        @test HalfInteger(2) == HalfInteger(twoX = 4)
        @test HalfInteger(-30) == HalfInteger(twoX = -60)
        @test HalfInteger(0//2) == HalfInteger(twoX = 0)
        @test HalfInteger(1//2) == HalfInteger(twoX = 1)
        @test HalfInteger(-5//2) == HalfInteger(twoX = -5)

        @test hi"0" == HalfInteger(0)
        @test hi"1" == HalfInteger(1)
        @test hi"210938" == HalfInteger(210938)
        @test hi"-15" == HalfInteger(-15)
        @test hi"1/2" == HalfInteger(1//2)
        @test hi"-3/2" == HalfInteger(-3//2)

        @test string(hi"0") == "0"
        @test string(hi"1") == "1"
        @test string(hi"-1") == "-1"
        @test string(hi"1/2") == "1/2"
        @test string(hi"-3/2") == "-3/2"

        @test_throws ArgumentError parse(HalfInteger, "")
        @test_throws ArgumentError parse(HalfInteger, "-50/100")
        @test_throws ArgumentError parse(HalfInteger, "1/3")

        @test hash(a) == hash(2)
        @test hash(b) == hash(1.5)
        @test hash(b) == hash(3//2)
    end

    @testset "HalfIntegerRange" begin
        using AtomicLevels: HalfIntegerRange
        @test length(HalfIntegerRange(hi"0", hi"0")) == 1
        @test length(HalfIntegerRange(hi"0", hi"2")) == 3
        @test length(HalfIntegerRange(hi"-1/2", hi"1/2")) == 2
        @test size(HalfIntegerRange(hi"-1/2", hi"1/2")) == (2,)

        @test hi"5" : hi"7" == HalfIntegerRange(HalfInteger(5), HalfInteger(7))
        @test hi"-1/2" : hi"1/2" == HalfIntegerRange(hi"-1/2", hi"1/2")

        @test collect(hi"0" : hi"2") == [hi"0", hi"1", hi"2"]
        @test collect(hi"-3/2" : hi"1/2") == [hi"-3/2", hi"-1/2", hi"1/2"]

        @test hi"1/2" ∈ hi"-1/2" : hi"1/2"
        @test 1 ∈ hi"0" : hi"2"
        @test 1//2 ∈ hi"-1/2" : hi"7/2"
        @test !(hi"1/2" ∈ hi"0" : hi"1")
        @test !(1//2 ∈ hi"-1" : hi"7")

        r = hi"-3/2" : hi"3/2"
        @test r[1] == hi"-3/2"
        @test r[2] == hi"-1/2"
        @test r[3] == hi"1/2"
        @test r[4] == hi"3/2"
        @test_throws BoundsError r[0]
        @test_throws BoundsError r[5]
    end
end
