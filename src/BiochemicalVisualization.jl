module BiochemicalVisualization

using BiochemicalAlgorithms

using Colors
import GeometryBasics
import Meshes
using JSServe
using LinearAlgebra
using MsgPack
using Statistics
using Rotations
using Dates
using StaticArrays
using EnumX

using KernelAbstractions


using BenchmarkTools # TODO not necessary for package
using Infiltrator


asset_path(parts...) = normpath(joinpath(@__DIR__, "..", "assets", parts...))

export asset_path

include("helpers/Logger.jl")

if Base.find_package("CUDA") !== nothing
    using CUDA
    using CUDA.CUDAKernels
    if(CUDA.functional())
        const backend = CUDABackend()
        CUDA.allowscalar(false)
        log_info(gpu, "Using CUDA")
        # CUDA.versioninfo()
    else
        const backend = CPU()
        log_info(gpu, "Using CPU")
    end
else
    const backend = CPU()
    log_info(gpu, "Using CPU")
end # TODO other backends

include("helpers/meshes/PlainMesh.jl")
include("helpers/meshes/PlainNonStdMesh.jl")
include("helpers/meshes/ColoredMesh.jl")
include("helpers/meshes/MeshHelpers.jl")
include("helpers/gpu_rotation_test.jl")
include("helpers/utils.jl")
include("helpers/point_filter.jl")
include("helpers/color_maps.jl")

include("core/representation.jl")




include("core/export.jl")

include("models/backbone/backbone_config.jl")
include("models/backbone/frame_construction.jl")
include("models/backbone/backbone.jl")
include("models/backbone/backbone_gpu.jl")
include("models/ball_and_stick.jl")
include("models/stick.jl")
include("models/van_der_waals.jl")

include("splines/SplineHelper.jl")
include("splines/CatmullRom.jl")
include("splines/CubicB.jl")


include("core/visualize.jl")


end # module BiochemicalVisualization
