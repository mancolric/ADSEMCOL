include("../../src/src_NonHydrostaticWaterWaves/NonHydrostaticWaterWaves.jl")
include("../ResUbi.jl")

function LatexString(model::NHWW, vble::String)

    if vble=="h"
        return "h"
    elseif vble=="q1"
        return "q_1"
    elseif vble=="q2"
        return "q_2"
    elseif vble=="q3"
        return "q_3"
    elseif vble=="v1"
        return "v_1"
    elseif vble=="v2"
        return "v_2"
    elseif vble=="v3"
        return "v_3"
    elseif vble=="vr"
        return "v_r"
    elseif vble=="b"
        return "b"
    elseif vble=="eta"
        return "\\eta"
    elseif vble=="p"
        return "p"
    elseif vble=="P"
        return "h\\xi"
    elseif vble=="xi_h"
        return "\\xi/h"
    else 
        @warn "Variable $(vble) not supported"
        return vble
    end
    
end

#Default options for plot:
# LabelSize   = 10
# TickSize    = 8
PyPlot.matplotlib.rc("mathtext",fontset="cm")
PyPlot.matplotlib.rc("font",family="serif")
