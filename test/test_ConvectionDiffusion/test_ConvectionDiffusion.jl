include("../ResUbi.jl")
include("../../src/src_ConvectionDiffusion/ConvectionDiffusion.jl")

function errL2L2_(t_n::Float64, errL2L2::Float64, t_np1::Float64, e_np1::Float64)

    errL2L2     = sqrt((errL2L2^2*t_n + e_np1^2*(t_np1-t_n))/t_np1)
    return errL2L2
    
end
