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


using BenchmarkTools # TODO not necessary for package


asset_path(parts...) = normpath(joinpath(@__DIR__, "..", "assets", parts...))

export asset_path

include("helpers/Logger.jl")
include("helpers/ColoredMesh.jl")
include("helpers/MeshHelpers.jl")

include("core/representation.jl")
include("core/visualize.jl")

include("splines/CatmullRom.jl")


include("core/export.jl")



include("models/backbone.jl")
include("models/ball_and_stick.jl")
include("models/stick.jl")
include("models/van_der_waals.jl")


end # module BiochemicalVisualization
