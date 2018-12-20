@testset "Levels & States" begin
    @testset "Construction" begin
        csfs_3s_3p = csfs(c"3s 3p")
        csf_1 = first(csfs_3s_3p)
        @test_throws ArgumentError Level(csf_1, hi"2")
        @test_throws ArgumentError State(Level(csf_1, hi"1"), hi"2")

        @test string(Level(csf_1, hi"1")) == "|3s(₁²S|²S) 3p(₁²Pᵒ|¹Pᵒ)-, J = 1⟩"
        @test string(State(Level(csf_1, hi"1"), hi"-1")) == "|3s(₁²S|²S) 3p(₁²Pᵒ|¹Pᵒ)-, J = 1, M_J = -1⟩"

        @test sort([State(Level(csf_1, hi"1"), M_J) for M_J ∈ reverse(hi"-1":hi"1")]) ==
            [State(Level(csf_1, hi"1"), M_J) for M_J ∈ hi"-1":hi"1"]

        @test states.(csfs_3s_3p) == [[[State(Level(csf_1, hi"1"), M_J) for M_J ∈ hi"-1":hi"1"]],
                                      [[State(Level(csfs_3s_3p[2], J), M_J) for M_J ∈ -J:J] for J ∈ hi"0":hi"2" ]]
    end
end
