include("../../src/src_NonHydrostaticWaterWaves/NonHydrostaticWaterWaves.jl")
include("ResUbi.jl")

function LatexString(vble::String)

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
        return "P"
    else 
        @warn "Variable $(vble) not supported"
        return vble
    end
    
end

function errL2L2_(t_n::Float64, errL2L2::Float64, t_np1::Float64, e_np1::Float64)

    errL2L2     = sqrt((errL2L2^2*t_n + e_np1^2*(t_np1-t_n))/t_np1)
    return errL2L2
    
end

#Default options for plot:
# LabelSize   = 10
# TickSize    = 8
PyPlot.matplotlib.rc("mathtext",fontset="cm")
PyPlot.matplotlib.rc("font",family="serif")
