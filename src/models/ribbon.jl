function prepare_ribbon_model(
    ac::AbstractAtomContainer{T}; 
    stick_radius=T(0.2), resolution=30) where {T<:Real}

    start_time = now()
    U = Float64
    if(T <: AbstractFloat)
        U = T
    end

    log_info(types, "Types: ", T, " ", U)

    chain_meshes::Vector{PlainMesh{U}} = []
    chain_colors = map(c->map(channel->Int(channel*255), (c.r, c.g, c.b)), collect(distinguishable_colors(nchains(ac)+1))[2:end])
    for (chain_num, chain) in enumerate(eachchain(ac))

        c_alphas = filter(x -> x.element==Elements.C && x.name=="CA", atoms(chain))
        @assert length(c_alphas)>=2 # TODO was sonst?

        from_protein_chain(chain)

        #push!(chain_meshes, spline_mesh)
    end
    #temp = merge_multiple_meshes(chain_meshes)
    #result = Representation(temp)
    #log_info(types, "Type of result: ", typeof(result))

    #log_info(time_info, "Generated backbone mesh in $((now()-start_time).value/1000) seconds. ($(length(result.vertices) ÷ 3) vertices)")

    #return result
end