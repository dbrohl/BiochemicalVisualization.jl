"""
Returns a tuple (R, G, B) with values in [0; 255]. 
Based on https://de.wikipedia.org/wiki/HSV-Farbraum#Umrechnung_HSV_in_RGB

* hue in [0;360]
* saturation in [0;1]
* value in [0;1]
"""
function hsv_to_rgb(hue, saturation, value)
    hi = floor(hue/60)
    f = (hue/60 - hi)
    p = value * (1-saturation)
    q = value * (1 - saturation*f)
    t = value * (1 - saturation*(1-f))

    value = Int(round(value*255))
    t = Int(round(t*255))
    p = Int(round(p*255))
    q = Int(round(q*255))

    if hi==0 || hi==6
        return (value, t, p)
    elseif hi==1
        return (q, value, p)
    elseif hi==2
        return (p, value, t)
    elseif hi==3
        return (p, q, value)
    elseif hi==4
        return (t, p, value)
    else
        return (value, p, q)
    end
end

function rgb_to_hex(rgb::NTuple{3, Int}; prefix="")
    result = prefix
    for channel=1:3
        hex = string(rgb[channel], base=16)
        if(length(hex)==1)
            hex = "0"*hex
        end
        result = result * hex
    end
    return result
end

"Converts a string (possibly containing a prefix like \"#\"), into an NTuple{3, Int}(r, g, b). "
function hex_to_rgb(hex::String)
    for i=length(hex)-1:-1:1
        if !('0'<=hex[i]<='9' || 'a'<=hex[i]<='f' || 'A'<=hex[i]<='F')
            hex = hex[i+1:end]
            break
        end
    end
    if length(hex)!=6
        throw(ArgumentError("hex string with length!=6: $hex"))
    end
    return (parse(Int, hex[1:2], base=16), parse(Int, hex[3:4], base=16), parse(Int, hex[5:6], base=16))
end


"""
Generates a rainbow color where pos in [0;1] determines the hue. 
"""
function rainbow(pos)
    return hsv_to_rgb(pos*360, 1, 1)
end

"""
Generates n different colors. 
"""
function n_colors(n)
    colors = Vector{NTuple{3, Int}}(undef, n)
    for (i, hue) in enumerate(collect(range(0, 360, n+1))[1:end-1])
        colors[i] = hsv_to_rgb(hue, 1, 1)
    end
    return colors
end

"""
Inplace cross()
"""
function cross!(dest::AbstractArray, a::AbstractVector, b::AbstractVector)
    if !(length(a) == length(b) == length(dest) == 3)
        throw(DimensionMismatch("cross product is only defined for vectors of length 3"))
    end
    a1, a2, a3 = a
    b1, b2, b3 = b
    dest .= (a2*b3-a3*b2, a3*b1-a1*b3, a1*b2-a2*b1)
end

""""
Normalizes a 3-element column of an array. 
"""
function normalize_col!(arr, i)
    dist = 0
    for j=1:3
        dist += arr[j, i]^2
    end
    dist = sqrt(dist)
    for j=1:3
        arr[j, i] /= dist
    end
end

function approx_zero(value)
    return abs(value)< 10^-5
end

"""
Handles logistics for generic models with more than one chain. Assumptions: 
* config!==nothing => config.color existiert
* single_chain_function is called with types (chain{T}, Union{config, Nothing}; fixed_color=Union{Nothing, NTuple{3, Int}})
* single_chain_function throws ErrorExceptions when something goes wrong
"""
function handle_multichain_model(ac::System{T}, config, fixed_color, single_chain_function) where T
    start_time = now()
    log_info(types, "Type: ", T)

    if config!==nothing && config.color==Color.CHAIN
        chain_colors = n_colors(nchains(ac))
    end

    chain_reps::Vector{Union{Missing, Representation{T}}} = fill(missing, (nchains(ac)))
    for (chain_idx, chain) in enumerate(BiochemicalAlgorithms.chains(ac))
        try
            color = nothing
            if config!==nothing && config.color==Color.UNIFORM && fixed_color!==nothing
                color = fixed_color
            end
            if config!==nothing && config.color==Color.CHAIN
                color = chain_colors[chain_idx]
            end

            chain_rep = single_chain_function(chain, config, fixed_color=color)

            chain_reps[chain_idx] = chain_rep
        catch e
            if e isa ErrorException #TODO
                log_warning("Skipped chain $(chain.name), because an error occured: $(e.msg)")
            else
                rethrow(e)
            end
        end
    end
    result_reps::Vector{Representation{T}} = filter(!ismissing, chain_reps)
    if(length(result_reps)==0)
        log_info(misc, "No chain representations were generated. ")
    end
    result = merge_representations(result_reps)
    log_info(types, "Type of result: ", typeof(result))

    log_info(time_info, "Generated representation in $((now()-start_time).value/1000) seconds. ")

    return result
end


function get_string_color(color::Color.T, atom, fixed_color::Union{NTuple{3, Int}, Nothing})
    if color==Color.UNIFORM || color==Color.CHAIN
        color = (fixed_color===nothing ? (255, 0, 0) : fixed_color)
        return rgb_to_hex(color, prefix="#")

    elseif color==Color.RAINBOW #rainbow per residue
        chain = parent_chain(atom)
        min_idx, max_idx = extrema(fragments_df(chain).idx)

        frag = parent_fragment(atom)
        return rgb_to_hex(rainbow((frag.idx-min_idx)/(max_idx-min_idx)), prefix="#")

    elseif color==Color.SECONDARY_STRUCTURE 
        frag = parent_fragment(atom)
        return rgb_to_hex(SS_COLORS[frag.properties[:SS]], prefix="#")

    elseif color==Color.RESIDUE
        frag = parent_fragment(atom)
        return rgb_to_hex(AA_COLORS[frag.name], prefix="#")

    elseif color==Color.ELEMENT
        return element_color_web(atom.element)
    else
        throw(ArgumentError("unknown color config $(color)"))
    end
end