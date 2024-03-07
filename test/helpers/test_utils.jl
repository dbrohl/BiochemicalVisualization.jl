@testitem "hsv_to_rgb" begin
    using BiochemicalVisualization: hsv_to_rgb

    @test hsv_to_rgb(0, 1, 1)==(255, 0, 0)
    @test hsv_to_rgb(60, 1, 1)==(255, 255, 0)
    @test hsv_to_rgb(0, 0, 1)==(255, 255, 255)
    @test hsv_to_rgb(0, 1, 0)==(0, 0, 0)
end

@testitem "rgb_to_hex and hex_to_rgb" begin
    using BiochemicalVisualization: rgb_to_hex, hex_to_rgb

    @test rgb_to_hex((0, 0, 0)) == "000000"
    @test rgb_to_hex((255, 0, 0)) == "ff0000"
    @test rgb_to_hex((255, 0, 0), prefix="#") == "#ff0000"
    @test rgb_to_hex((255, 0, 0), prefix="0x") == "0xff0000"

    @test hex_to_rgb("000000")==(0, 0, 0)
    @test hex_to_rgb("ff0000")==(255, 0, 0)
    @test hex_to_rgb("#ff0000")==(255, 0, 0)
    @test hex_to_rgb("0xff0000")==(255, 0, 0)
end

@testitem "n_colors" begin
    using BiochemicalVisualization: n_colors
    colors = n_colors(30)
    @test allunique(colors)
end

@testitem "approx_zero" begin
    using BiochemicalVisualization: approx_zero
    @test approx_zero(0)
    @test !approx_zero(1)
    @test approx_zero(0.000001)
    @test approx_zero(-0.000001)
end