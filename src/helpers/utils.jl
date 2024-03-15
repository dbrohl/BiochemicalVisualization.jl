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

function rgb_to_hex(rgb; prefix="")
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

function approx_zero(value)
    return abs(value)< 10^-5
end

function isless_tolerance(a, b, tol=10^-4)
    return (a-b) <= -tol
end

function islessorqual_tolerance(a, b, tol=10^-4)
    return (a - b) < tol
end

function isgreater_tolerance(a, b, tol=10^-4)
    return (a-b) >= tol
end

function isgreaterorequal_tolerance(a, b, tol=10^-4)
    return (a-b) > -tol
end

function iszero_tolerance(a, tol=10^-4)
    return abs(a) < tol
end

function isequal_tolerance(a, b, tol=10^-4)
    return abs(a - b) < tol
end

function isnotequal_tolerance(a, b, tol=10^-4)
    return abs(a - b) >= tol
end

