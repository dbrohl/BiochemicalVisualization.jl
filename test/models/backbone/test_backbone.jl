@testitem "backbone" begin
    using BiochemicalAlgorithms
    using BiochemicalVisualization: prepare_backbone_model, BackboneConfig
    function loadPDB(path)
        fdb = FragmentDB()
        pdb = load_pdb(path)
        normalize_names!(pdb, fdb);
        reconstruct_fragments!(pdb, fdb);
        add_secondary_structures!(pdb, path)
        return pdb
    end

    path = normpath(joinpath(@__DIR__, "..", "..", "data", "2mma.pdb"))
    pdb = loadPDB(path)

    configs = [
            #backbone types
            BackboneConfig(Float32(0.2), 12, BackboneType.BACKBONE, Color.UNIFORM, Spline.CATMULL_ROM, ControlPoints.C_ALPHA, Frame.RMF, Filter.ANGLE),
            BackboneConfig(Float32(0.2), 12, BackboneType.RIBBON, Color.UNIFORM, Spline.CATMULL_ROM, ControlPoints.C_ALPHA, Frame.RMF, Filter.ANGLE),
            BackboneConfig(Float32(0.2), 12, BackboneType.CARTOON, Color.UNIFORM, Spline.CATMULL_ROM, ControlPoints.C_ALPHA, Frame.RMF, Filter.ANGLE),
            
            #colors
            BackboneConfig(Float32(0.2), 12, BackboneType.BACKBONE, Color.CHAIN, Spline.CATMULL_ROM, ControlPoints.C_ALPHA, Frame.RMF, Filter.ANGLE),
            BackboneConfig(Float32(0.2), 12, BackboneType.BACKBONE, Color.RAINBOW, Spline.CATMULL_ROM, ControlPoints.C_ALPHA, Frame.RMF, Filter.ANGLE),
            BackboneConfig(Float32(0.2), 12, BackboneType.BACKBONE, Color.SECONDARY_STRUCTURE, Spline.CATMULL_ROM, ControlPoints.C_ALPHA, Frame.RMF, Filter.ANGLE),
            BackboneConfig(Float32(0.2), 12, BackboneType.BACKBONE, Color.RESIDUE, Spline.CATMULL_ROM, ControlPoints.C_ALPHA, Frame.RMF, Filter.ANGLE),

            # other config
            BackboneConfig(Float32(0.2), 12, BackboneType.BACKBONE, Color.UNIFORM, Spline.CATMULL_ROM, ControlPoints.MID_POINTS, Frame.RMF, Filter.ANGLE),
            BackboneConfig(Float32(0.2), 12, BackboneType.BACKBONE, Color.UNIFORM, Spline.CATMULL_ROM, ControlPoints.MID_POINTS, Frame.SECOND_SPLINE, Filter.ANGLE),
            BackboneConfig(Float32(0.2), 12, BackboneType.BACKBONE, Color.UNIFORM, Spline.CUBIC_B, ControlPoints.C_ALPHA, Frame.RMF, Filter.ANGLE),
            BackboneConfig(Float32(0.2), 12, BackboneType.BACKBONE, Color.UNIFORM, Spline.CUBIC_B, ControlPoints.MID_POINTS, Frame.RMF, Filter.ANGLE),
            BackboneConfig(Float32(0.2), 12, BackboneType.BACKBONE, Color.UNIFORM, Spline.CUBIC_B, ControlPoints.MID_POINTS, Frame.SECOND_SPLINE, Filter.ANGLE)
    ]
    prepare_backbone_model(pdb, configs[1])
    for (i, config) in enumerate(configs)
        println("Config $i: $config")
        stats = @timed prepare_backbone_model(pdb, config)
        println("$(stats.time)")
        @test stats.time<2
    end
end