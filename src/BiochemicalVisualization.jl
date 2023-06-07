module BiochemicalVisualization

using BiochemicalAlgorithms

using Colors
import GeometryBasics
using Meshes
using JSServe
using LinearAlgebra
using MsgPack
using Statistics


asset_path(parts...) = normpath(joinpath(@__DIR__, "..", "assets", parts...))

export asset_path

include("core/representation.jl")
include("core/visualize.jl")
include("core/export.jl")

include("helpers/ColoredMesh.jl")

include("models/backbone.jl")
include("models/ball_and_stick.jl")
include("models/stick.jl")
include("models/van_der_waals.jl")


end # module BiochemicalVisualization
