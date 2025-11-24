include("../../src/src_CompressibleFlow/CompressibleFlow.jl")
include("InputData/CompressibleFlow_SC.jl")
include("../ResUbi.jl")

function LatexString(model::GasModel, vble::String)

    if vble=="rho"
        return "\\rho"
    elseif vble=="rhovx"
        return "m_1"
    elseif vble=="rhovy"
        return "m_2"
    elseif vble=="rhoE"
        return "\\mathcal{E}"
    elseif vble=="Y_F"
        return "Y_F"
    elseif vble=="vx"
        return "v_x"
    elseif vble=="vy"
        return "v_y"
    elseif vble=="v"
        return "||\\mathbf{v}||"
    elseif vble=="rhoe"
        return "\\rho e"
    elseif vble=="e"
        return "e"
    elseif vble=="p"
        return "p"
    elseif vble=="a"
        return "a"
    elseif vble=="lambda_max"
        return "\\lambda_{max}"
    elseif vble=="M"
        return "M"
    elseif vble=="RT"
        return "RT"
    elseif vble=="T"
        return "T"
    elseif vble=="epsilon"
        return "\\epsilon"
    elseif vble=="nu"
        return "\\nu"
    elseif vble=="beta"
        return "\\beta"
    elseif vble=="kappa_rho_cv"
        return "D_T"
    elseif vble=="rhos"
        return "\\rho s"
    elseif vble=="gamma"
        return "\\gamma"
    elseif vble=="e_i"
        return "e_i"
    elseif length(vble)>=5 && vble[1:5]=="rhoY_"
        S   = vble[6:length(vble)]
        return "\\rho_{$(S)}"
    elseif length(vble)>=2 && vble[1:2]=="Y_"
        S   = vble[3:length(vble)]
        return "Y_{$(S)}"
    elseif length(vble)>=2 && vble[1:2]=="X_"
        S   = vble[3:length(vble)]
        return "X_{$(S)}"
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
