@testitem "approx_zero" begin
    using BiochemicalVisualization: approx_zero
    @test approx_zero(0)
    @test !approx_zero(1)
    @test approx_zero(0.000001)
    @test approx_zero(-0.000001)
end