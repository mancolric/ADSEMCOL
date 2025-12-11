include("../../src/src_ShallowWater/ShallowWater.jl")
include("../ResUbi.jl")

function LatexString(model::SWE, vble::String)

    if vble=="h"
        return "h"
    elseif vble=="q1"
        return "q_1"
    elseif vble=="q2"
        return "q_2"
    elseif vble=="v1"
        return "v_1"
    elseif vble=="v2"
        return "v_2"
    elseif vble=="vr"
        return "v_r"
    elseif vble=="v"
        return "|\\mathbf{v}|"
    elseif vble=="Fr"
        return "\\mathrm{Fr}"
    elseif vble=="b"
        return "b"
    elseif vble=="eta"
        return "\\eta"
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
