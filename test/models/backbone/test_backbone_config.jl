@testitem "BackboneConfig" begin
    using BiochemicalVisualization: BackboneConfig, complete_config, add_to_config!

    # check that types in struct are working
    @test_throws MethodError BackboneConfig()
    @test_throws MethodError BackboneConfig(0.2, 1.5,
                                            12, 
                                            Color.SECONDARY_STRUCTURE, 
                                            Color.SECONDARY_STRUCTURE, 
                                            Spline.CUBIC_B, 
                                            ControlPoints.MID_POINTS, 
                                            Frame.SECOND_SPLINE, 
                                            Filter.ANGLE)
    @test_throws MethodError BackboneConfig(nothing, 1.5,
                                            12, 
                                            Color.SECONDARY_STRUCTURE, 
                                            Color.SECONDARY_STRUCTURE, 
                                            Spline.CUBIC_B, 
                                            ControlPoints.MID_POINTS, 
                                            Frame.SECOND_SPLINE, 
                                            Filter.ANGLE)

    # test == operator
    a = BackboneConfig(0.2, 1.5,
        12, 
        BackboneType.CARTOON, 
        Color.SECONDARY_STRUCTURE, 
        Spline.CUBIC_B, 
        ControlPoints.MID_POINTS, 
        Frame.SECOND_SPLINE, 
        Filter.ANGLE)
    b = BackboneConfig(0.3, 1.5, 12, BackboneType.CARTOON, Color.SECONDARY_STRUCTURE, Spline.CUBIC_B, ControlPoints.MID_POINTS, Frame.SECOND_SPLINE, Filter.ANGLE)
    c = BackboneConfig(0.2, 1.5, 14, BackboneType.CARTOON, Color.SECONDARY_STRUCTURE, Spline.CUBIC_B, ControlPoints.MID_POINTS, Frame.SECOND_SPLINE, Filter.ANGLE)
    d = BackboneConfig(0.2, 1.5, 12, BackboneType.BACKBONE, Color.SECONDARY_STRUCTURE, Spline.CUBIC_B, ControlPoints.MID_POINTS, Frame.SECOND_SPLINE, Filter.ANGLE)
    e = BackboneConfig(0.2, 1.5, 12, BackboneType.CARTOON, Color.RESIDUE, Spline.CUBIC_B, ControlPoints.MID_POINTS, Frame.SECOND_SPLINE, Filter.ANGLE)
    f = BackboneConfig(0.2, 1.5, 12, BackboneType.CARTOON, Color.SECONDARY_STRUCTURE, Spline.CATMULL_ROM, ControlPoints.MID_POINTS, Frame.SECOND_SPLINE, Filter.ANGLE)
    g = BackboneConfig(0.2, 1.5, 12, BackboneType.CARTOON, Color.SECONDARY_STRUCTURE, Spline.CUBIC_B, ControlPoints.C_ALPHA, Frame.SECOND_SPLINE, Filter.ANGLE)
    h = BackboneConfig(0.2, 1.5, 12, BackboneType.CARTOON, Color.SECONDARY_STRUCTURE, Spline.CUBIC_B, ControlPoints.MID_POINTS, Frame.RMF, Filter.ANGLE)
    i = BackboneConfig(0.2, 1.5, 12, BackboneType.CARTOON, Color.SECONDARY_STRUCTURE, Spline.CUBIC_B, ControlPoints.MID_POINTS, Frame.SECOND_SPLINE, Filter.NONE)
    j = BackboneConfig(0.2, 1.0, 12, BackboneType.CARTOON, Color.SECONDARY_STRUCTURE, Spline.CUBIC_B, ControlPoints.MID_POINTS, Frame.SECOND_SPLINE, Filter.NONE)

    @test a==a
    for other in [b, c, d, e, f, g, h, i, j]
        @test a!=other
    end


    # test complete_config
    empty = PartialBackboneConfig()
    large = PartialBackboneConfig(tube_radius=100)
    full = BackboneConfig(0.2, 1.5,
                            12, 
                            BackboneType.CARTOON, 
                            Color.SECONDARY_STRUCTURE, 
                            Spline.CUBIC_B, 
                            ControlPoints.MID_POINTS, 
                            Frame.SECOND_SPLINE, 
                            Filter.ANGLE)

    @test complete_config(empty, full)==full
    @test complete_config(large, full)!=full
    @test complete_config(large, full).tube_radius==100
    @test complete_config(large, full).resolution_cross==12

    # test add_to_config!
    empty = PartialBackboneConfig()
    large = PartialBackboneConfig(tube_radius=100)
    full = PartialBackboneConfig(0.2, 1.5,
                            12, 
                            BackboneType.CARTOON, 
                            Color.SECONDARY_STRUCTURE, 
                            Spline.CUBIC_B, 
                            ControlPoints.MID_POINTS, 
                            Frame.SECOND_SPLINE, 
                            Filter.ANGLE)

    a = deepcopy(empty)
    add_to_config!(a, large)
    @test a==large

    a = deepcopy(empty)
    add_to_config!(a, full)
    @test a==full

    a = deepcopy(large)
    add_to_config!(a, full)
    @test a!=large && a!=full
    @test a.tube_radius==100
    @test a.resolution_cross==12

    a.tube_radius = 0.2
    @test a==full    


    a = deepcopy(full)
    add_to_config!(a, large)
    @test a==full
end