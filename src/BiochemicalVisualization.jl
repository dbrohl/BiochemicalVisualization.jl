module BiochemicalVisualization

using BiochemicalAlgorithms

using Colors
import GeometryBasics
using Meshes
using JSServe
using LinearAlgebra
using MsgPack
using Statistics
using Rotations
using Dates
using StaticArrays

using KernelAbstractions


using BenchmarkTools # TODO not necessary for package


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

include("helpers/PlainMesh.jl")
include("helpers/PlainNonStdMesh.jl")
include("helpers/MeshHelpers.jl")

include("core/representation.jl")
include("core/visualize.jl")

include("splines/CatmullRom.jl")


include("core/export.jl")



include("models/backbone.jl")
include("models/backbone_gpu.jl")
include("models/ball_and_stick.jl")
include("models/stick.jl")
include("models/van_der_waals.jl")


end # module BiochemicalVisualization
