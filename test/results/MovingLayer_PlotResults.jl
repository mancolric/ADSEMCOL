include("PlotResults.jl")

#Initial condition:
function utheor_MovingLayer(t::Float64, x::Vector{<:AAF64})
    
    u   = @avxt @. 2.0 + tanh((x[1]-1.0*t)/deltax)
    return [u]
    
end

function CompareSpaceAdaptation_MovingLayer(StudyCase::String; nb::Int=1, 
    SaveFig::Bool=false, w::Float64=8.50, h::Float64=8.50)

    SCvv0   = NaN
    SCvv1   = NaN
#     nb      = NaN
    
    if StudyCase=="normal"
    
        #correct etaS, dynamic AMA_rDeltat:
        SCvv1       = [ 30063:30066, 
                        30067:30070, 
                        30071:30074 ]
#         nb          = 1
        
    elseif StudyCase=="normal_L2"
    
        #correct etaS, dynamic AMA_rDeltat:
        SCvv1       = [ 30075:30078, 
                        30079:30082, 
                        30083:30086 ]
#         nb          = 1
      
    end
    
    #------------------------------------------------------------
    
    hpvv1, errvv1, etavv1, 
        tCPUvv1, pvv1       = GetVbles(SCvv1, ["hp", "errL2L2", "etaL2L2", "tCPU", "FesOrder"], nb=nb)
    EOCvv1                  = ExpOrderConv(hpvv1, errvv1)

    PyPlotFigure(w=w, h=h, bottom=1.5)
    leg                     = Legend(pvv1, "p", format="%d")
    loglog(NaN,NaN, "-k", linewidth=0.5)
    loglog(NaN,NaN, "--k", linewidth=0.5)
    leg                     = vcat(leg, L"e_{ST}", latexstring("\\mathcal{E}_{ST}"))
    PlotXY(hpvv1, errvv1, linestyle="solid", marker="s")
    PlotXY(hpvv1, etavv1, linestyle="dashed", marker="s")
    xlabel(L"h/p")
    tick_params(axis="both", which="both", labelsize=TickSize)
    legend(leg)
    grid("on")
    if SaveFig
        savefig("$(FigUbi)MovingLayer_SpaceAdapt1.eps", dpi=800, pad_inches=0)
    end
    
    PyPlotFigure(w=w, h=h, bottom=1.5)
    leg                     = Legend(pvv1, "p", format="%d")
    PlotXY(tCPUvv1, errvv1, linestyle="solid", marker="s")
    xlabel(latexstring(GetString("tCPU")))
    ylabel(latexstring(GetString("errL2L2")), rotation=0, labelpad=7)
    tick_params(axis="both", which="both", labelsize=TickSize)
    legend(leg)
    grid("on")
    if SaveFig
        savefig("$(FigUbi)MovingLayer_SpaceAdapt2.eps", dpi=800, pad_inches=0)
    end
    
    display(EOCvv1)
    
    return
    
end

function CompareTimeAdaptation_MovingLayer(StudyCase::String; 
    SaveFig::Bool=false, w::Float64=8.50, h::Float64=8.50)

    SCvv0   = NaN
    SCvv1   = NaN
    nb      = NaN
    
    if StudyCase=="0"
    
        #No SpaceAdapt:
        SCvv0       = [ 30034:30038 ]
                    
        #SpaceAdapt:
        SCvv1       = [ 30039:30043 ]
                        
        nb          = 1
        
    elseif StudyCase=="normal"
    
        #correct etaS, dynamic AMA_rDeltat:
        SCvv0       = [ 30053:30057 ]
        SCvv1       = [ 30058:30062 ]
        nb          = 1
      
    end
    
    #------------------------------------------------------------
    
    Deltatvv0, errvv0, etavv0, 
        tCPUvv0                 = GetVbles(SCvv0, ["Deltat", "errL2L2", "etaL2L2", "tCPU"], 
                                    nb=nb)
    EOCvv0                      = ExpOrderConv(Deltatvv0, errvv0)
    Deltatvv1, errvv1, etavv1, 
        tCPUvv1                 = GetVbles(SCvv1, ["Deltat", "errL2L2", "etaL2L2", "tCPU"],
                                    nb=nb)
    EOCvv1                      = ExpOrderConv(Deltatvv1, errvv1)
    
#     figure()
#     PlotXY(Deltatvv0, errvv0, linestyle="solid", marker="o")
#     PlotXY(Deltatvv0, etavv0, linestyle="dashed", marker="o")
#     PlotXY(Deltatvv1, errvv1, linestyle="solid", marker="s")
#     PlotXY(Deltatvv1, etavv1, linestyle="dashed", marker="s")
#     grid("on")
    
    PyPlotFigure(w=w, h=h, bottom=1.5)
    loglog(Deltatvv0[1], errvv0[1], color="r", linewidth=0.5, linestyle="solid", marker="o", markersize=3.5)
    loglog(Deltatvv0[1], etavv0[1], color="r", linewidth=0.5, linestyle="dashed", marker="o", markersize=3.5)
    loglog(Deltatvv1[1], errvv1[1], color="b", linewidth=0.5, linestyle="solid", marker="s", markersize=3.5)
    loglog(Deltatvv1[1], etavv1[1], color="b", linewidth=0.5, linestyle="dashed", marker="s", markersize=3.5)
    xlabel(latexstring("\\tau"))
    leg     = vcat(L"Unif., e_{ST}", latexstring("Unif., \\eta_{ST}"), 
                    L"Adapt., e_{ST}", latexstring("Adapt., \\eta_{ST}"))
    legend(leg)
    grid("on")
    if SaveFig
        savefig("$(FigUbi)Test1_TimeAdapt1.png", dpi=800, pad_inches=0)
    end
            
#     figure()
#     PlotXY(tCPUvv0, errvv0, linestyle="solid", marker="o")
#     PlotXY(tCPUvv1, errvv1, linestyle="solid", marker="s")
#     grid("on")
    
    PyPlotFigure(w=w, h=h, bottom=1.5, left=1.5)
    loglog(tCPUvv0[1], errvv0[1], color="r", linewidth=0.5, linestyle="solid", marker="o", markersize=3.5)
    loglog(tCPUvv1[1], errvv1[1], color="b", linewidth=0.5, linestyle="solid", marker="s", markersize=3.5)
    grid("on")
    xlabel(latexstring(GetString("tCPU")))
    ylabel(latexstring(GetString("errL2L2")), rotation=0)
    leg     = vcat(L"Unif.", L"Adapt.")
    legend(leg)
    if SaveFig
        savefig("$(FigUbi)Test1_TimeAdapt2.png", dpi=800, pad_inches=0)
    end
    
    display(EOCvv0)
    display(EOCvv1)
    
    return
    
end

function TableResults_SpaceAdapt_MovingLayer(StudyCase::String)

    vbles       = [ "FesOrder", "TolS_max", 
                    "NDOF",     "hp", 
                    "errL2L2",  "EOC_hp_errL2L2",   "etaST",    "EOC_hp_etaST", 
                    "iST",      "CFLmax",           "tCPU"]
    formats     = [ "%d",       "%.2E",     
                    "%.0f",     "%.2E", 
                    "%.2E",     "%.2f",             "%.2E",     "%.2f", 
                    "%.2f",     "%.2E",             "%.1f"]
    
    #------------------------------------------------------------
    
    SCvv1   = NaN
    nb      = NaN
    
    if StudyCase=="normal"
    
        #correct etaS, dynamic AMA_rDeltat:
        SCvv1       = [ 30063:30066, 
                        30067:30070, 
                        30071:30074 ]
        nb          = 1

    elseif StudyCase=="normal_L2"
    
        #correct etaS, dynamic AMA_rDeltat:
        SCvv1       = [ 30075:30078, 
                        30079:30082, 
                        30083:30086 ]
        nb          = 1

    end
    
    #Save results for IMA and AMA:
    table           = string("\\begin{tabular}{", repeat("r", length(vbles)), "} \n", 
                        "\\hline \n", 
                        "\$p\$ & ", 
                        "\$\\Tol_S\$ & \$\\NDOF\$ & \$h/p\$ & \$e_{ST}\$ ", 
                        "& \$\\EOC\$ & \$\\Ecal_{ST}\$ & \$\\EOC\$ & \$i_{ST}\$ ", 
                        "& \$\\CFL^{\\max}\$ & \$t_{\\CPU} [s]\$ \\\\", 
                        "\\hline \n")
    for ii=1:length(SCvv1)
        SCv         = SCvv1[ii]
        FesOrder    = GetVbles(SCv[1], nb, ["FesOrder"])[1]
        table_SC    = TableVbles(SCv, vbles, formats, nb=nb)
        table       = string(table, table_SC, "\\hline \n")
    end
    table           = string(table, "\\end{tabular}")
    print(table)
    
    write("$(FigUbi)MovingLayer_SpaceAdapt.txt", table)
    
    return
    
end

function TableResults_MovingLayer(StudyCase::String)
       
    vbles       = [ "TolS_bool",    "TolT", 
                    "Nt",           "Deltat", 
                    "errL2L2",      "EOC_Deltat_errL2L2",   "etaST",    "EOC_Deltat_etaST", 
                    "iST",          "CFLmax"       ]
    formats     = [ "%.2E",         "%.2E",     
                    "%d",           "%.2E", 
                    "%.2E",         "%.2F",                 "%.2E",     "%.2F",     
                    "%.3F",         "%.2E"      ]
    
#     header      = "\$\\Tol_S\$ & \$\\Tol_T\$ & \$\\NDOF\$ & \$N_t\$ & \$\\tau\$ & \$e_{ST}\$ & \$\\EOC\$ & \$\\Ecal_{ST}\$ & \$\\EOC\$ & \$i_{ST}\$ & \$t_{\\CPU} [s]\$ \\\\" 
    header      = "\$\\Tol_S\$ & \$\\Tol_T\$ & \$N_t\$ & \$\\tau\$ & \$e_{ST}\$ & \$\\EOC\$ & \$\\Ecal_{ST}\$ & \$\\EOC\$ & \$i_{ST}\$ & \$\\CFL^{\\max}\$ \\\\" 
    
    #------------------------------------------------------------
    
    SCvv0   = NaN
    SCvv1   = NaN
    nb      = NaN
    
    if StudyCase=="0"
    
        #------------------------------------------------------------
        #New test1 (tanh in time, boundary layer in space), 
        #deltat=0.2, deltax=0.001, TolT=0.01*TolS:
        #
        #Now with factor AMA_TolS/etaS, AMA_OptimalSizeEq, AMA_SizeOrder=FesOrder, AMA_AnisoOrder=2
        
        SCvv0       = [ 30034:30038 ]
        SCvv1       = [ 30039:30043 ]
        nb          = 1
        
    elseif StudyCase=="0_CFL"
    
        #------------------------------------------------------------
        #New test1 (tanh in time, boundary layer in space), 
        #deltat=0.2, deltax=0.001, TolT=0.01*TolS:
        #
        #Now with factor AMA_TolS/etaS, AMA_OptimalSizeEq, AMA_SizeOrder=FesOrder, AMA_AnisoOrder=2
        
        SCvv0       = [ 30034:30038 ]
        SCvv1       = [ 30048:30052 ]
        nb          = 1
        
    elseif StudyCase=="normal"
    
        #correct etaS, dynamic AMA_rDeltat:
        SCvv0       = [ 30053:30057 ]
        SCvv1       = [ 30058:30062 ]
        nb          = 1
        
    elseif StudyCase=="normal_L2"
    
        #Since the space error is negligible, we take the same results than for "normal".
        
        #correct etaS, dynamic AMA_rDeltat:
        SCvv0       = [ 30053:30057 ]
        SCvv1       = [ 30058:30062 ]
        nb          = 1
        
    end
    
    #Save results for IMA and AMA:
    table           = string("\\begin{tabular}{", repeat("r", length(vbles)), "} \n", 
                        "\\hline \n", 
                        header,
                        "\\hline \n")
    for ii=1:length(SCvv0)
        SCv         = SCvv0[ii]
        FesOrder    = GetVbles(SCv[1], nb, ["FesOrder"])[1]
        table_SC    = TableVbles(SCv, vbles, formats, nb=nb)
        table       = string(table, table_SC, "\\hline \n")
        SCv         = SCvv1[ii]
        table_SC    = TableVbles(SCv, vbles, formats, nb=nb)
        table       = string(table, table_SC, "\\hline \n")
    end
    table           = string(table, "\\end{tabular}")
    print(table)
    
#     write("$(FigUbi)MovingLayer_table.txt", table)
    
    return
    
end

function NElems_MovingLayer(; SaveFig::Bool=false, w::Float64=8.50, h::Float64=8.50)

    @warn "Not used"
    
    SC      = 30039
    nb      = 1
    
    solver  = GetSolver(SC, nb)
    validv  = ones(Bool, length(solver.tv))
    if true
        @. validv  = solver.validv .== 1
    end
    
    PyPlotFigure(w=w, h=h, bottom=1.2, left=2.0)
    semilogy(solver.tv[validv], solver.nElemsv[validv], "k", linewidth=0.5)
    xlabel(L"t")
    ylabel(latexstring("\\#\\mathcal{T}^{h,n}"), rotation=0, labelpad=5.0)
    tick_params(axis="both", which="both", labelsize=TickSize)
    
    SC      = 30039
    nb      = 1
    
    solver  = GetSolver(SC, nb)
    validv  = ones(Bool, length(solver.tv))
    if true
        @. validv  = solver.validv .== 1
    end
    
    if SaveFig
        savefig("$(FigUbi)MovingLayer_NElems.png", dpi=800, pad_inches=0)
    end
    
    return
    
end

function Deltat_MovingLayer(; SaveFig::Bool=false, w::Float64=6.50, h::Float64=6.50)

#     SC0     = 30035
#     SC1     = 30040
#     nb      = 1
    
    SC0     = 30054
    SC1     = 30059
    nb      = 1
    
    #Without space adaptation:
    solver  = GetSolver(SC0, nb)
    validv  = ones(Bool, length(solver.tv))
    markv   = @. solver.validv == 0 && solver.tv>0.0
    Deltatv = vcat(diff(solver.tv), NaN)
    
    fig     = PyPlotFigure(w=w, h=h, bottom=1.0, left=1.2)
    semilogy(solver.tv, Deltatv, "m", linewidth=0.7, marker="x", markersize=4.0)
    
    #With space adaptation:
    solver          = GetSolver(SC1, nb)
    Deltatv         = vcat(diff(solver.tv), NaN)
    aux             = @. Deltatv <= 0.0
    Deltatv         = Deltatv[.!aux]
    solver.tv       = solver.tv[.!aux]
    
    semilogy(solver.tv, Deltatv, "c", linewidth=0.7, marker="+", markersize=4.0)
    
    legend(["Fixed mesh", 
            string("AMA., ", latexstring("\\mathrm{Tol}_S^{\\max}=", sprintf1("%.2E", solver.TolS_max)))],
            fontsize=TickSize, 
            loc="lower left")
    tick_params(axis="both", which="both", labelsize=TickSize)
    
#     xlim(5.0, 7.0)
    xlim(1.0, 3.0)
    ylim(0.9*0.01, 1.1*0.1)
    
    xlabel(latexstring("t"))
    ylabel(latexstring("\\tau"), rotation=0, labelpad=0.0)
    
    if SaveFig
        savefig("$(FigUbi)MovingLayer_Deltat.png", dpi=800, pad_inches=0)
    end
    
    return
    
end

function Errors_MovingLayer(; SaveFig::Bool=false, w::Float64=7.50, h::Float64=6.50)

#     SC      = 30040
#     SC      = 30049
    SC      = 30059
    nb      = 1
    solver  = GetSolver(SC, nb)
    
    #Erase t=0:
    aux                 = @. solver.tv == 0.0
    solver.tv[aux]      .= NaN
    solver.etaSv[aux]   .= NaN
    solver.etaTv[aux]   .= NaN
    
    fig     = PyPlotFigure(w=w, h=h, bottom=1.0, left=1.0)
    markv   = @. solver.validv == 0 && solver.tv>0.0
    
    #all steps:
    semilogy(solver.tv, solver.etaSv, "b", linewidth=0.7, marker="x", markersize=1.0)
    semilogy(solver.tv, solver.etaTv, "g", linewidth=0.7, marker="x", markersize=1.0)
    
    #tolerances:
    semilogy(solver.tv, 0.0.*solver.tv .+ solver.TolS_max, "--b", linewidth=0.7)
    semilogy(solver.tv, 0.0.*solver.tv .+ solver.TolT, "--g", linewidth=0.7)
    
    #rejected steps:
    semilogy(solver.tv[markv], solver.etaSv[markv], "r", linestyle="none", marker="s", markersize=1.5)
    semilogy(solver.tv[markv], solver.etaTv[markv], "r", linestyle="none", marker="s", markersize=1.5)

    legend([latexstring("\\mathcal{E}_{S},~\\mathrm{all~steps}"), 
            latexstring("\\mathcal{E}_{T},~\\mathrm{all~steps}"), 
            latexstring("\\mathrm{Tol}_{S}"), 
            latexstring("\\mathrm{Tol}_{T}"),
            latexstring("\\mathrm{Rejected~steps}")], 
            loc=(0.01,0.42), 
            fontsize=TickSize)
    
    tick_params(axis="both", which="both", labelsize=TickSize)
    
#     xlim(5.0, 7.0)
    xlim(1.0, 3.0)
    ylim(0.9e-8, 1.2*solver.TolT)
    xlabel(L"t")
    
    if SaveFig
        savefig("$(FigUbi)MovingLayer_Errors.png", dpi=800, pad_inches=0)
    end
    
    return
    
end

function Errors2_MovingLayer(SC::Int, nb::Int; SaveFig::Bool=false, w::Float64=7.50, h::Float64=6.50)

    solver      = GetSolver(SC, nb)
    filename    = GetFileName(SC,nb)
    
    #Erase t=0:
    aux                 = @. solver.tv == 0.0
    solver.tv[aux]      .= NaN
    solver.etaSv[aux]   .= NaN
    solver.etaTv[aux]   .= NaN
    solver.etaAv[aux]   .= NaN
    
    fig     = PyPlotFigure(w=w, h=h, bottom=1.0, left=1.0)
    markv   = @. solver.validv == 1 && solver.tv>0.0 && solver.etaAv>1e-12
    
    #all steps:
    semilogy(solver.tv[markv], solver.etaSv[markv], "b", linewidth=0.7, marker="x", markersize=1.0)
    semilogy(solver.tv[markv], solver.etaTv[markv], "g", linewidth=0.7, marker="x", markersize=1.0)
    semilogy(solver.tv[markv], solver.etaAv[markv], "r", linewidth=0.7, marker="x", markersize=1.0)
    errL2v  = load(filename, "errLqv")
    semilogy(unique(solver.tv[solver.validv.==1]), errL2v, "k", linewidth=0.7, marker="x", markersize=1.0)
    
    #tolerances:
    semilogy(solver.tv, 0.0.*solver.tv .+ solver.TolS_max, "--b", linewidth=0.7)
    semilogy(solver.tv, 0.0.*solver.tv .+ solver.TolT, "--g", linewidth=0.7)
    
    legend([latexstring("\\mathcal{E}_{S}"), 
            latexstring("\\mathcal{E}_{T}"),
            latexstring("\\mathcal{E}_{A}"), 
            latexstring("||u^{h\\tau}-u||_X"), 
            latexstring("\\mathrm{Tol}_{S}"), 
            latexstring("\\mathrm{Tol}_{T}")],
            loc=(0.01,0.42), 
            fontsize=TickSize)
    
    tick_params(axis="both", which="both", labelsize=TickSize)
    
#     xlim(5.0, 7.0)
#     xlim(1.0, 3.0)
    ylim(1e-11, 1e-3)
    xlabel(L"t")
    
    if SaveFig
        savefig("$(FigUbi)MovingLayer_Errors_SC$(SC)_$(nb).png", dpi=800, pad_inches=0)
    end
    
    return
    
end
