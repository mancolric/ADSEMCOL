include("PlotResults.jl")

#Initial condition:
function utheor_Test1(t::Float64, x::Vector{<:AAF64})
    
    r   = @avxt @. sqrt(x[1]*x[1]+x[2]*x[2])+1e-14
    r0  = sqrt(2)/2.0
    u   = @avxt @. cos(lambda1*x[1])*cos(lambda2*x[2]) - tanh((r-r0)/delta) + 2.5
    return [u]
    
end
    
function Contour_Test1(SC::Int, nb::Int; SaveFig::Bool=false, w::Float64=8.50, h::Float64=8.50)
    
    solver      = GetSolver(SC, nb)
    
    uP          = GetPolynomials(solver.u, solver.fes)
    delta       = 1e-4
    umin        = minimum(uP[1])-delta
    umax        = maximum(uP[1])+delta
    Nz          = 2000
    iv          = 1:Nz
    lambda      = 10.0/Nz
    zv          = linspace(umin, umax, Nz)
    
#     PyPlotFigure(w=w, h=h, top=0.7, bottom=0.7, left=1.0, right=0.5)
    PyPlotFigure(w=w, h=h, top=0.0, bottom=0.0, left=0.0, right=0.0)
    v_plot      = PlotContour(solver.u[1], solver.fes, cmap="jet", zv=zv)
    PlotMesh!(solver.mesh, color="k") 
    tick_params(axis="both", which="both", labelsize=TickSize)
#     title(latexstring("p=",solver.FesOrder), fontsize=LabelSize)
    text(0.1, 0.1, latexstring("p=",solver.FesOrder), color="w", fontsize=LabelSize)
    axis("off")
    if SaveFig
        savefig("$(FigUbi)SC$(SC).png", dpi=800, pad_inches=0)
    end
    
    zv          = vcat(linspace(umin, 0.9*umax-delta, 100), linspace(0.9*umax, umax, 1000))
#     PyPlotFigure(w=w, h=h, top=0.7, bottom=0.7, left=1.1, right=0.5)
    PyPlotFigure(w=w, h=h, top=0.0, bottom=0.0, left=0.0, right=0.0)
    PlotContour(solver.u[1], solver.fes, cmap="jet", zv=zv)
    PlotMesh!(solver.mesh, color="k") 
    xlim(0.97,1.0)
    ylim(0.97,1.0)
    axis("off")
    tick_params(axis="both", which="both", labelsize=TickSize)
#     title(latexstring("p=",solver.FesOrder), fontsize=LabelSize)
    text(0.973, 0.973, latexstring("p=",solver.FesOrder), color="w", fontsize=LabelSize)
    if SaveFig
        savefig("$(FigUbi)SC$(SC)_Zoom.png", dpi=800, pad_inches=0)
    end
    
    #Colorbar (Careful: if the colorbar is small, visualization in screen is different 
    #from saved figure).
    PyPlotFigure(w=1.2, h=h, bottom=0.2, top=0.2, left=0.9, right=0.0)
    PyPlotColorbar!(minimum(v_plot), maximum(v_plot), 4, tickformat="%.2f")
    if SaveFig
        savefig("$(FigUbi)SC$(SC)_Colorbar.png", dpi=800, pad_inches=0)
    end
        
    return
    
end

function CompareSpaceAdaptation_Test1(StudyCase::String; nb::Int=1,
    SaveFig::Bool=false, w::Float64=8.50, h::Float64=8.50)

    #Best case: "001_3"
    
    #------------------------------------------------------------
    
    SCvv0   = NaN
    SCvv1   = NaN
#     nb      = NaN
    
    if StudyCase=="InitialCondition"
    
        #delta=0.02 (initial condition)
        
        #NO AMA:
        SCvv0       = [ 20101:20110, 
                        20111:20120,
                        20121:20130 ]
                        
        SCvv1       = [ 20131:20140, 
                        20141:20150,
                        20151:20160 ]
                        
        nb          = 0
                
    elseif StudyCase=="InitialCondition2"
    
        #------------------------------------------------------------
        #delta=0.02 (initial condition, new reconstruction)
        
        #NO AMA:
        SCvv0       = [ 20161:20170, 
                        20171:20180,
                        20181:20190 ]
                        
        SCvv1       = [ 20191:20199, 
                        20201:20210,
                        20211:20220 ]
                    
        nb          = 0
        
    elseif StudyCase=="01"
    
        #------------------------------------------------------------
        #New test1 (tanh in time, boundary layer in space), 
        #deltat=0.2, deltax=0.01, TolT=0.01*TolS:
        
        #IMA:
        SCvv0       = [ 20245:20248, 
                        20249:20252,
                        20253:20256 ]
                    
        #AMA:
        SCvv1       = [ 20257:20260, 
                        20261:20264,
                        20265:20268 ]
                        
        nb          = 1
    
    elseif StudyCase=="001"
    
        #------------------------------------------------------------
        #New test1 (tanh in time, boundary layer in space), 
        #deltat=0.2, deltax=0.001, TolT=0.01*TolS:
        
        #IMA:
        SCvv0       = [ 20269:20272, 
                        20273:20276,
                        20277:20280 ]
                    
        #AMA:
        SCvv1       = [ 20281:20284, 
                        20285:20288,
                        20289:20292 ]
                        
        nb          = 1
        
    elseif StudyCase=="001_2"
    
        #------------------------------------------------------------
        #New test1 (tanh in time, boundary layer in space), 
        #deltat=0.2, deltax=0.001, TolT=0.01*TolS:
        #
        #Now with factor AMA_TolS/etaS, AMA_SizeOrder=FesOrder, AMA_AnisoOrder=3
        
        #IMA:
        SCvv0       = [ 20316:20319, 
                        20320:20323,
                        20324:20327 ]
                    
        #AMA:
        SCvv1       = [ 20328:20331, 
                        20332:20335,
                        20336:20339 ]
                        
        nb          = 1
        
    elseif StudyCase=="001_3"
    
        #------------------------------------------------------------
        #New test1 (tanh in time, boundary layer in space), 
        #deltat=0.2, deltax=0.001, TolT=0.01*TolS:
        #
        #Now with factor AMA_TolS/etaS, AMA_OptimalSizeEq, AMA_SizeOrder=FesOrder, AMA_AnisoOrder=2
        
        #IMA:
        SCvv0       = [ 20316:20319, 
                        20320:20323,
                        20324:20327 ] #Same as 001_2 (AMA_SizeOrder=FesOrder, AMA_AnisoOrder=0)
                    
        #AMA:
        SCvv1       = [ 20452:20455, 
                        20456:20459,
                        20460:20463 ]
#                         20488:20491 ]
                        
        nb          = 1
        
    elseif StudyCase=="001_4"
    
        #------------------------------------------------------------
        #New test1 (tanh in time, boundary layer in space), 
        #deltat=0.2, deltax=0.001, TolT=0.01*TolS:
        #
        #Now with factor AMA_TolS/etaS, AMA_OptimalSizeEq, AMA_SizeOrder=2, AMA_AnisoOrder=2
        
        #IMA:
        SCvv0       = [ 20464:20467, 
                        20468:20471,
                        20472:20474 ]
                    
        #AMA:
        SCvv1       = [ 20476:20479, 
                        20480:20483,
                        20484:20487 ]
                        
        nb          = 1
        
    elseif StudyCase=="001_CFL"
    
        #With correct CFL and parameters set as in article:
        SCvv0       = [ 20504:20507, 
                        20508:20511, 
                        20512:20515 ]
        SCvv1       = [ 20492:20495, 
                        20496:20499, 
                        20500:20503 ]
        nb          = 1
        
    elseif StudyCase=="normal"
    
        #Correct etaS, dynamic AMA_rDeltat:
        
        #With correct CFL and parameters set as in article:
        SCvv0       = [ 20528:20531, 
                        20532:20535, 
                        20536:20539 ]
        SCvv1       = [ 20516:20519, 
                        20520:20523, 
                        20524:20527 ]
#         nb          = 1
        
    elseif StudyCase=="normal_L2"
    
        #Correct etaS, dynamic AMA_rDeltat:
        
        #With correct CFL and parameters set as in article:
        SCvv0       = [ 20552:20555, 
                        20556:20559, 
                        20560:20563 ] 
#         SCvv1       = [ 20540:20543, 
#                         20544:20547, 
#                         20548:20551 ]
        SCvv1       = [ vcat(20541:20543,20575), 
                        vcat(20572,20544:20547), 
                        vcat(20573,20574,20548,20549) ]
                        
#         nb          = 1

    elseif StudyCase=="normal_L2_bubbles"
    
        #ignoring bubbles at reconstruction:
        SCvv0       = [ 20540:20543, 
                        20544:20547, 
                        20548:20551 ]
        #not ignoring them:
        SCvv1       = [ 20576:20579, 
                        20580:20583,
                        20584:20587 ]
        
        
    end
        
    #------------------------------------------------------------
    
    hpvv0, errvv0, etavv0, 
        tCPUvv0, pvv0       = GetVbles(SCvv0, ["hp", "errL2L2", "etaL2L2", "tCPU", "FesOrder"], 
                                nb=nb)
    EOCvv0                  = ExpOrderConv(hpvv0, errvv0)
    hpvv1, errvv1, etavv1, 
        tCPUvv1             = GetVbles(SCvv1, ["hp", "errL2L2", "etaL2L2", "tCPU"], nb=nb)
    EOCvv1                  = ExpOrderConv(hpvv1, errvv1)
    
#     @warn "Hiding IMA"
    PyPlotFigure(w=w, h=h, bottom=1.5)
    leg                     = Legend(pvv0, "p", format="%d")
    loglog(NaN,NaN, color="k", marker="o", markersize=3.5, linestyle="none")
    loglog(NaN,NaN, color="k", marker="s", markersize=3.5, linestyle="none")
    loglog(NaN,NaN, "-k", linewidth=0.5)
    loglog(NaN,NaN, "--k", linewidth=0.5)
    leg                     = vcat(leg, L"IMA", L"AMA", L"e_{ST}", latexstring("\\mathcal{E}_{ST}"))
    PlotXY(hpvv0, errvv0, linestyle="solid", marker="o")
    PlotXY(hpvv0, etavv0, linestyle="dashed", marker="o")
    PlotXY(hpvv1, errvv1, linestyle="solid", marker="s")
    PlotXY(hpvv1, etavv1, linestyle="dashed", marker="s")
    xlabel(L"h/p")
    tick_params(axis="both", which="both", labelsize=TickSize)
    legend(leg)
    grid("on")
    if SaveFig
        savefig("$(FigUbi)Test1_SpaceAdapt1.eps", dpi=800, pad_inches=0)
    end
    
    PyPlotFigure(w=w, h=h, bottom=1.5, left=1.5)
    leg                     = Legend(pvv0, "p", format="%d")
    loglog(NaN,NaN, color="k", marker="o", markersize=3.5, linestyle="none")
    loglog(NaN,NaN, color="k", marker="s", markersize=3.5, linestyle="none")
    
    leg                     = vcat(leg, L"IMA", L"AMA")
    PlotXY(tCPUvv0, errvv0, linestyle="solid", marker="o")
    PlotXY(tCPUvv1, errvv1, linestyle="solid", marker="s")
    grid("on")
    xlabel(latexstring(GetString("tCPU")))
    ylabel(latexstring(GetString("errL2L2")), rotation=0, labelpad=7)
    tick_params(axis="both", which="both", labelsize=TickSize)
    legend(leg)
    if SaveFig
        savefig("$(FigUbi)Test1_SpaceAdapt2.eps", dpi=800, pad_inches=0)
    end
    
    PyPlotFigure(w=w, h=h, bottom=1.5, left=1.5)
    leg                     = Legend(pvv0, "p", format="%d")
    PlotXY(tCPUvv1, errvv1, linestyle="solid", marker="s")
    grid("on")
    xlabel(latexstring(GetString("tCPU")), fontsize=LabelSize)
    ylabel(latexstring(GetString("errL2L2")), rotation=0, fontsize=LabelSize, labelpad=10.0)
    legend(leg, fontsize=TickSize)
    tick_params(axis="both", which="both", labelsize=TickSize)
    if SaveFig
        savefig("$(FigUbi)Test1_SpaceAdapt3.png", dpi=800, pad_inches=0)
    end
    
    display(EOCvv0)
    display(EOCvv1)
    
    return
    
end

function CompareTimeAdaptation_Test1(StudyCase::String; 
    SaveFig::Bool=false, w::Float64=8.50, h::Float64=8.50)

    SCvv0   = NaN
    SCvv1   = NaN
    nb      = NaN
    
    if StudyCase=="01"
    
        #------------------------------------------------------------
        #New test1 (tanh in time, boundary layer in space), 
        #deltat=0.01, deltax=0.1.
        #
        #In principle, one could be tempted to set TolS=0.01*TolT. However,
        #note that space adaptation implies also time adaptation, so the only
        #way to study time adaptation is to set SpaceAdapt=false.
        #
        #Furthermore, if there is no TimeAdapt, it is impossible to prescribe 
        #a value for TolT.
        
        #No TimeAdapt:
        SCvv0       = [ 20293:20295,
                        20293:20295 ]
                    
        #TimeAdapt:
        SCvv1       = [ 20297:20300,
                        20297:20300 ]
                        
        nb          = 1
        
    elseif StudyCase=="001"
    
        #------------------------------------------------------------
        #New test1 (tanh in time, boundary layer in space), 
        #deltat=0.001, deltax=0.1.
        #
        #In principle, one could be tempted to set TolS=0.01*TolT. However,
        #note that space adaptation implies also time adaptation, so the only
        #way to study time adaptation is to set SpaceAdapt=false.
        #
        #Furthermore, if there is no TimeAdapt, it is impossible to prescribe 
        #a value for TolT.
        
        #No TimeAdapt:
        SCvv0       = [ 20301:20303,
                        [20312,20313,20314,20301] ]
                    
        #TimeAdapt:
        SCvv1       = [ vcat(0,20305:20308) ]
                        
        #Best results:
        SCvv0       = [ [20312,20313,20314,20301] ]
        SCvv1       = [ vcat(20315,20305:20308) ]
        
        nb          = 1
        
    elseif StudyCase=="v01" #test 2, variable BC
    
        #------------------------------------------------------------
        #New test1 (tanh in time, boundary layer in space), 
        #deltat=0.01, deltax=0.1.
        #
        #In principle, one could be tempted to set TolS=0.01*TolT. However,
        #note that space adaptation implies also time adaptation, so the only
        #way to study time adaptation is to set SpaceAdapt=false.
        #
        #Furthermore, if there is no TimeAdapt, it is impossible to prescribe 
        #a value for TolT.
        
        #No TimeAdapt:
        SCvv0       = [ 30025:30028 ]
                    
        #TimeAdapt:
        SCvv1       = [ 30025:30028 ]
                               
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
    xlabel(L"N_t^{-1}")
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

function TableResults_SpaceAdapt_Test1(StudyCase::String)

    vbles       = [ "FesOrder", "SAType",           "TolS_max", 
                    "NDOF",     "hp", 
                    "errL2L2",  "EOC_hp_errL2L2",   "etaST",    "EOC_hp_etaST", 
                    "iST",      "CFLmax",           "tCPU"]
    formats     = [ "%d",       "%s",               "%.2E",     
                    "%.0f",     "%.2E", 
                    "%.2E",     "%.2f",             "%.2E",     "%.2f", 
                    "%.2f",     "%.2E",             "%.1f"]
    
    #------------------------------------------------------------
    
    SCvv0   = NaN
    SCvv1   = NaN
    nb      = NaN
    
    if StudyCase=="001_3"
    
        #------------------------------------------------------------
        #New test1 (tanh in time, boundary layer in space), 
        #deltat=0.2, deltax=0.001, TolT=0.01*TolS:
        #
        #Now with factor AMA_TolS/etaS, AMA_OptimalSizeEq, AMA_SizeOrder=FesOrder, AMA_AnisoOrder=2
        
        #IMA:
        SCvv0       = [ 20316:20319, 
                        20320:20323,
                        20324:20327 ] #Same as 001_2 (AMA_SizeOrder=FesOrder, AMA_AnisoOrder=0)
                    
        #AMA:
        SCvv1       = [ 20452:20455, 
                        20456:20459,
                        20460:20463 ]
                        
        nb          = 1
        
    elseif StudyCase=="001_CFL"
    
        #With correct CFL and parameters set as in article:
        SCvv0       = [ 20504:20507, 
                        20508:20511, 
                        20512:20515 ]
        SCvv1       = [ 20492:20495, 
                        20496:20499, 
                        20500:20503 ]
        nb          = 1
        
    elseif StudyCase=="normal_L2"
    
        #Correct etaS, dynamic AMA_rDeltat:
        
        #With correct CFL and parameters set as in article:
        SCvv0       = [ 20552:20555, 
                        20556:20559, 
                        20560:20563 ] 
        SCvv1       = [ 20540:20543, 
                        20544:20547, 
                        20548:20551 ]
        nb          = 1
        
    end
    
    #Save results for IMA and AMA:
    table           = string("\\begin{tabular}{", repeat("r", length(vbles)), "} \n", 
                        "\\hline \n", 
                        "\$p\$ & \\makecell[cr]{Space\\\\Adapt.} & ", 
                        "\$\\Tol_S\$ & \$\\NDOF\$ & \$h/p\$ & \$e_{ST}\$ ", 
                        "& \$\\EOC\$ & \$\\Ecal_{ST}\$ & \$\\EOC\$ & \$i_{ST}\$ ", 
                        "& \$\\CFL^{\\max}\$ & \$t_{\\CPU} [s]\$ \\\\", 
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
    
    write("$(FigUbi)Test1_SpaceAdapt.txt", table)
    
    return
    
end

function TableResults_TimeAdapt_Test1(StudyCase::String)

    vbles       = [ "TolT_bool", 
                    "Nt",       "Deltat", 
                    "errL2L2",  "EOC_Deltat_errL2L2",   "etaST",    "EOC_Deltat_etaST", 
                    "iST",      "CFLmax",               "tCPU"]
    formats     = [ "%.2E",     
                    "%.0f",     "%.2E", 
                    "%.2E",     "%.2f",                 "%.2E",     "%.2f", 
                    "%.2f",     "%.2E",                 "%.1f"]
    
    #------------------------------------------------------------
    
    SCvv0   = NaN
    SCvv1   = NaN
    nb      = NaN
    
    if StudyCase=="001"
    
        #------------------------------------------------------------
        #New test1 (tanh in time, boundary layer in space), 
        #deltat=0.2, deltax=0.001, TolT=0.01*TolS:
        #
        #Now with factor AMA_TolS/etaS, AMA_OptimalSizeEq, AMA_SizeOrder=FesOrder, AMA_AnisoOrder=2
        
        SCvv0       = [ [20312,20313,20314,20301] ]
        SCvv1       = [ vcat(20315,20305,20306,20307) ]
        nb          = 1
        
    elseif StudyCase=="normal" || StudyCase=="normal_L2"
    
        SCvv0       = [ 20564:20567 ]
        SCvv1       = [ 20568:20571 ]
        nb          = 1
        
    end
    
    #Save results for IMA and AMA:
    table           = string("\\begin{tabular}{", repeat("r", length(vbles)), "} \n", 
                        "\\hline \n", 
                        "\$\\Tol_T\$ & \$N_t\$ & \$\\tau\$ & \$e_{ST}\$ & \$\\EOC\$ & \$\\Ecal_{ST}\$ & \$\\EOC\$ & \$i_{ST}\$ & \$\\CFL^{\\max}\$ & \$t_{\\CPU} [s]\$ \\\\", 
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
    
    write("$(FigUbi)Test1_TimeAdapt.txt", table)
    
    return
    
end
