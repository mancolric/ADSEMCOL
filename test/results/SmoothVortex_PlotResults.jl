include("PlotResults.jl")

function utheor_SmoothVortex(t::Float64, x::Vector{<:VecOrMat{Float64}})
        
    #Problem data:
    gamma           = 1.4
    vortex_st       = 5.0   #vortex strength
    u_inf           = 1.0   #vortex horizontal velocity
    
    xm          = hcat(x[1][:], x[2][:])
    rho         = rho_SmoothVortex(xm, gamma, vortex_st, t, 0.0, u_inf, 0.0)
    rhovx       = rho .* vx_SmoothVortex(xm, gamma, vortex_st, t, 0.0, u_inf, 0.0)
    rhovy       = rho .* vy_SmoothVortex(xm, gamma, vortex_st, t, 0.0, u_inf, 0.0)
    rhoE        = @. rho * $e_SmoothVortex(xm, gamma, vortex_st, t, 0.0, u_inf, 0.0) + 
                            0.5*(rhovx^2 + rhovy^2)/rho
                            
    rho         = reshape(rho, size(x[1]))
    rhovx       = reshape(rhovx, size(x[1]))
    rhovy       = reshape(rhovy, size(x[1]))
    rhoE        = reshape(rhoE, size(x[1]))
    
    return [rho, rhovx, rhovy, rhoE]
    
end
    
function Lines_SmoothVortex(SC::Int, nb::Int; 
    SaveFig::Bool=false, w::Float64=9.50, h::Float64=6.50,
    PlotVars::Vector{String}=["rho", "v", "p", "M"], mFig::Int=2, nFig::Int=2)

    FileName    = GetFileName(SC, nb)
    solver      = GetSolver(SC, nb)
    
    #Theoretical solution:
    x1v     = linspace(-7.0, 17.0, 500)
    x2v     = zeros(size(x1v))
    xm      = [ x1v, x2v ]
    um      = utheor_SmoothVortex(solver.t, xm)
    
    #Interpolate solution:
    u_terp, = SolutionCompute(solver.u, solver.fes, [x1v, x2v])
    
    #Define gas:
    GasModel                = GasIdeal()
    GasModel.gamma          = 1.4
    GasModel.epsilon        = 0.0
    GasModel.nu             = 0.0
    GasModel.beta           = 0.0
    GasModel.kappa_rho_cv   = 0.0
#     GasModel.CSS            = 0.1
    
    #Loop plot variables:
    fig     = PyPlotSubPlots(mFig, nFig, w=w, h=h, left=0.9, right=0.4, bottom=1.1, top=1.0)
    for ii=1:length(PlotVars)
        PyPlot.subplot(mFig, nFig, ii)
        PyPlot.cla()
        plot(x1v, DepVars(GasModel, solver.t, Vector{<:AMF64}(xm), Vector{<:AMF64}(um), [PlotVars[ii]])[1][1][:], "--k", linewidth=0.5)
        plot(x1v, DepVars(GasModel, solver.t, Vector{<:AMF64}(xm), Vector{<:AMF64}(u_terp), [PlotVars[ii]])[1][1][:], "b", linewidth=0.5)
        xlabel(latexstring("x_1"), fontsize=10)
        title(latexstring(LatexString(PlotVars[ii]),
            "; t^n=", sprintf1("%.2e", solver.t)), 
            fontsize=10)
    end
    if SaveFig
        savefig("$(FigUbi)SC$(SC)_Line1.png", dpi=800, pad_inches=0)
    end
    
    return
    
end

function CompareSpaceAdaptation_SmoothVortex(StudyCase::String; nb::Int=1, 
    SaveFig::Bool=false, w::Float64=8.50, h::Float64=8.50)

    #------------------------------------------------------------
    
    SCvv0   = NaN
    SCvv1   = NaN
#     nb      = NaN
    
    #Choose normal5_L2:
    
    if StudyCase=="normal"

                        
        SCvv1       = [ 50001:50005, 
                        50006:50010,
                        50011:50015
                        ]
                        
        nb          = 1
        
    elseif StudyCase=="normal2" #SizeOrder=FesOrder, AnisoOrder=2

                        
        SCvv1       = [ 50020:50024, 
                        50025:50029,
                        50030:50034
                        ]
                        
        nb          = 1
        
    elseif StudyCase=="normal3" #SizeOrder=FesOrder, AnisoOrder=2, with Eq. loop

                        
        SCvv1       = [ 50035:50038, 
                        50040:50044,
                        50045:50049
                        ]
                        
        nb          = 1
        
    elseif StudyCase=="normal4" #SizeOrder=FesOrder, AnisoOrder=2, with Eq. loop only for initial condition

                        
        SCvv1       = [ 50050:50053, 
                        50055:50059,
                        50060:50064
                        ]
                        
        nb          = 1
        
    elseif StudyCase=="normal5" #with dynamic AMA_rDeltat and correct etaS:
    
        SCvv1       = [ 50066:50069, 
                        vcat(50078,50071:50073),
                        50074:50077
                        ]
                        
#         nb          = 0

     elseif StudyCase=="normal5_L2" #with dynamic AMA_rDeltat and correct etaS:
    
        SCvv1       = [ 50083:50086, 
                        50087:50090,
                        50091:50094
                        ]
                        
#         nb          = 0

    elseif StudyCase=="normal5_L2_SS" #with dynamic AMA_rDeltat and correct etaS:
    
        SCvv1       = [ vcat(50108:50112,50123,50124), 
                        vcat(50113:50117,50125),
                        vcat(50118:50122,50127)
                        ]
        
        #All with TolT=0.01TolS:
        SCvv1       = [ vcat(50129:50133,50123,50124), 
                        vcat(50134:50138,50125),
                        vcat(50139:50142,50122,50127)
                        ]
                        
#         nb          = 0
        
    end
        
    #------------------------------------------------------------
    
    hpvv1, errvv1, etavv1, 
        tCPUvv1, pvv1       = GetVbles(SCvv1, ["hp", "errL2L2", "etaL2L2", "tCPU", "FesOrder"], nb=nb)
    EOCvv1                  = ExpOrderConv(hpvv1, errvv1)
    
    PyPlotFigure(w=w, h=h, bottom=1.5, left=1.7)
    leg                     = Legend(pvv1, "p", format="%d")
    loglog(NaN,NaN, "-k", linewidth=0.5)
    loglog(NaN,NaN, "--k", linewidth=0.5)
    leg                     = vcat(leg, L"e_{ST}", latexstring("\\mathcal{E}_{ST}"))
    PlotXY(hpvv1, errvv1, linestyle="solid", marker="s")
    PlotXY(hpvv1, etavv1, linestyle="dashed", marker="s")
    xlabel(L"h/p", fontsize=LabelSize)
    tick_params(axis="both", which="both", labelsize=TickSize)
    legend(leg, fontsize=TickSize)
    grid("on")
    if SaveFig
        savefig("$(FigUbi)SmoothVortex_SpaceAdapt1.png", dpi=800, pad_inches=0)
    end
    
    PyPlotFigure(w=w, h=h, bottom=1.5, left=1.5)
    leg                     = Legend(pvv1, "p", format="%d")
    PlotXY(tCPUvv1, errvv1, linestyle="solid", marker="s")
    grid("on")
    xlabel(latexstring(GetString("tCPU")), fontsize=LabelSize)
    ylabel(latexstring(GetString("errL2L2")), rotation=0, fontsize=LabelSize, labelpad=10.0)
    legend(leg, fontsize=TickSize)
    tick_params(axis="both", which="both", labelsize=TickSize)
    if SaveFig
        savefig("$(FigUbi)SmoothVortex_SpaceAdapt2.png", dpi=800, pad_inches=0)
    end
    
    println("EOC")
    display(EOCvv1)
    println("i_eff")
    for ii=1:length(etavv1)
        println(etavv1[ii]./errvv1[ii])
    end
    
    return
    
end

function CompareTimeAdaptation_SmoothVortex(StudyCase::String; 
    SaveFig::Bool=false, w::Float64=8.50, h::Float64=8.50)

    SCvv0   = NaN
    SCvv1   = NaN
    nb      = NaN
    
    if StudyCase=="normal"
    
        #TimeAdapt:
        SCvv1       = [ 50016:50019 ]
                        
        nb          = 1
        
    end
    
    #------------------------------------------------------------
    
    Deltatvv1, errvv1, etavv1, 
        tCPUvv1, CFLvv1         = GetVbles(SCvv1, ["Deltat", "errL2L2", "etaL2L2", 
                                                    "tCPU", "CFLmax"], nb=nb)
    EOCvv1                      = ExpOrderConv(Deltatvv1, errvv1)
    
    PyPlotFigure(w=w, h=h, bottom=1.5)
    loglog(Deltatvv1[1], errvv1[1], color="b", linewidth=0.5, linestyle="solid", marker="s", markersize=3.5)
    loglog(Deltatvv1[1], etavv1[1], color="b", linewidth=0.5, linestyle="dashed", marker="s", markersize=3.5)
    xlabel(latexstring("\\tau"))
    leg     = vcat(L"e_{ST}", latexstring("\\mathcal{E}_{ST}"))
    legend(leg)
    tick_params(axis="both", which="both", labelsize=TickSize)
    grid("on")
    if SaveFig
        savefig("$(FigUbi)SmoothVortex_TimeAdapt1.png", dpi=800, pad_inches=0)
    end
               
    PyPlotFigure(w=w, h=h, bottom=1.5, left=1.5)
    loglog(tCPUvv1[1], errvv1[1], color="b", linewidth=0.5, linestyle="solid", marker="s", markersize=3.5)
    grid("on")
    xlabel(latexstring(GetString("tCPU")))
    ylabel(latexstring(GetString("errL2L2")), rotation=0)
    tick_params(axis="both", which="both", labelsize=TickSize)
    if SaveFig
        savefig("$(FigUbi)SmoothVortex_TimeAdapt2.png", dpi=800, pad_inches=0)
    end
    
    PyPlotFigure(w=w, h=h, bottom=1.5, left=1.5)
    loglog(CFLvv1[1], errvv1[1], color="b", linewidth=0.5, linestyle="solid", marker="s", markersize=3.5)
    grid("on")
    xlabel(latexstring(GetString("CFLmax")))
    ylabel(latexstring(GetString("errL2L2")), rotation=0)
    tick_params(axis="both", which="both", labelsize=TickSize)
    if SaveFig
        savefig("$(FigUbi)SmoothVortex_TimeAdapt3.png", dpi=800, pad_inches=0)
    end
    
    display(EOCvv1)
    
    return
    
end

function TableResults_SpaceAdapt_SmoothVortex(StudyCase::String)

    vbles       = [ "FesOrder", "TolS_max", 
                    "NDOF",     "hp", 
                    "errL2L2",  "EOC_hp_errL2L2",   "etaST",    "EOC_hp_etaST", 
                    "iST",      "CFLmax",           "tCPU"]
    formats     = [ "%d",       "%.2E",     
                    "%.0f",     "%.2E", 
                    "%.2E",     "%.2f",             "%.2E",     "%.2f", 
                    "%.2f",     "%.2E",             "%.1f"]
    
    header      = string("\$p\$ & \$\\Tol_S\$ & \$\\NDOF\$ & \$h/p\$ & ",
                        "\$e_{ST}\$ & \$\\EOC\$ & \$\\Ecal_{ST}\$ & \$\\EOC\$ & ",
                        "\$i_{ST}\$ & \$\\CFL^{\\max}\$ & \$t_{\\CPU} [s]\$ \\\\")
 
#     vbles       = [ "FesOrder", "TolS_max", 
#                     "NDOF",     "hp", 
#                     "errL2L2",  "EOC_hp_errL2L2",   "etaST",    "EOC_hp_etaST", 
#                     "iST",      "tCPU"]
#     formats     = [ "%d",       "%.2E",     
#                     "%.0f",     "%.2E", 
#                     "%.2E",     "%.2E",             "%.2E",     "%.2E", 
#                     "%.2E",     "%.1f"]
#     
#     header      = string("\$p\$ & \$\\Tol_S\$ & \$\\NDOF\$ & \$h/p\$ & ",
#                         "\$e_{ST}\$ & \$\\EOC\$ & \$\\Ecal_{ST}\$ & \$\\EOC\$ & ",
#                         "\$i_{ST}\$ & \$t_{\\CPU} [s]\$ \\\\")
                        
    #------------------------------------------------------------
    
    SCvv1   = NaN
    nb      = NaN
    
    if StudyCase=="normal4" #SizeOrder=FesOrder, AnisoOrder=2, with Eq. loop, correct CFL

                        
        SCvv1       = [ 50050:50054, 
                        50055:50059,
                        50060:50064
                        ]
                        
        nb          = 1
        
    elseif StudyCase=="normal5" #with dynamic AMA_rDeltat and correct etaS:
    
        SCvv1       = [ 50066:50069, 
                        vcat(50078,50071:50073),
                        50074:50077
                        ]
                        
        nb          = 1
        
    elseif StudyCase=="normal5_L2" #with dynamic AMA_rDeltat and correct etaS:
    
        SCvv1       = [ 50083:50086, 
                        50087:50090,
                        50091:50094
                        ]
                        
        nb          = 1
        
    end
    
    #Save results for IMA and AMA:
    table           = string("\\begin{tabular}{", 
                        repeat("r", length(vbles)), 
                        "} \n", 
                        "\\hline \n", 
                        header, 
                        "\\hline \n")
    for ii=1:length(SCvv1)
        SCv         = SCvv1[ii]
        table_SC    = TableVbles(SCv, vbles, formats, nb=nb)
        table       = string(table, table_SC, "\\hline \n")
    end
    table           = string(table, "\\end{tabular}")
    print(table)
    
    write("$(FigUbi)IsentropicVortex_SpaceAdapt.txt", table)
    
    return
    
end

function TableResults_SpaceAdapt_SS_SmoothVortex(StudyCase::String)

    vbles       = [ "FesOrder", 
                    "NDOF",     "hp", 
                    "errL2L2",  "EOC_hp_errL2L2",   "etaST",    "EOC_hp_etaST", 
                    "iST",      "CFLmax",           "tCPU"]
    formats     = [ "%d",       
                    "%.0f",     "%.2E", 
                    "%.2E",     "%.2f",             "%.2E",     "%.2f", 
                    "%.2f",     "%.2E",             "%.1f"]
    
    header      = string("\$p\$ & \$\\NDOF\$ & \$h/p\$ & ",
                        "\$e_{ST}\$ & \$\\EOC\$ & \$\\Ecal_{ST}\$ & \$\\EOC\$ & ",
                        "\$i_{ST}\$ & \$\\CFL^{\\max}\$ & \$t_{\\CPU} [s]\$ \\\\")
                        
    #------------------------------------------------------------
    
    SCvv1   = NaN
    nb      = NaN
    
    if StudyCase=="normal5_L2_SS" #NO AMA
    
        #All with TolT=0.01TolS:
        SCvv1       = [ vcat(50129:50133,50123,50124), 
                        vcat(50134:50138,50125),
                        vcat(50139:50142,50122,50127)
                        ]
        #Last 4 cases only:
        SCvv1       = [ vcat(50132:50133,50123,50124), 
                        vcat(50136:50138,50125),
                        vcat(50141:50142,50122,50127)
                        ]
                        
        nb          = 1
        
    end
    
    #Save results for IMA and AMA:
    table           = string("\\begin{tabular}{", 
                        repeat("r", length(vbles)), 
                        "} \n", 
                        "\\hline \n", 
                        header, 
                        "\\hline \n")
    for ii=1:length(SCvv1)
        SCv         = SCvv1[ii]
        table_SC    = TableVbles(SCv, vbles, formats, nb=nb)
        table       = string(table, table_SC, "\\hline \n")
    end
    table           = string(table, "\\end{tabular}")
    print(table)
    
    write("$(FigUbi)IsentropicVortex_SS_SpaceAdapt.txt", table)
    
    return
    
end

function TableResults_TimeAdapt_SmoothVortex(StudyCase::String)

    vbles       = [ "TolT", 
                    "Nt",       "Deltat", 
                    "errL2L2",  "EOC_Deltat_errL2L2",   
                    "etaST",    "EOC_Deltat_etaST", 
                    "iST",      "CFLmax", 
                    "tCPU"]
    formats     = [ "%.2E",       
                    "%.0f",     "%.2E", 
                    "%.2E",     "%.2f", 
                    "%.2E",     "%.2f", 
                    "%.2f",     "%.2E", 
                    "%.1f"]
    
    header      = string("\$\\Tol_T\$ & \$N_t\$ & \$\\tau\$ & ",
                        "\$e_{ST}\$ & \$\\EOC\$ & \$\\Ecal_{ST}\$ & \$\\EOC\$ & ",
                        "\$i_{ST}\$ & \$\\CFL^{\\max}\$ & \$t_{\\CPU} [s]\$ \\\\")
    
    #------------------------------------------------------------
    
    SCvv1   = NaN
    nb      = NaN
    
    if StudyCase=="normal"
    
        #TimeAdapt:
        SCvv1       = [ 50016:50019 ]
                        
        nb          = 1
        
    elseif StudyCase=="normal5"
    
        SCvv1       = [ 50079:50082 ]
                        
        nb          = 1
        
    elseif StudyCase=="normal5_L2"
    
        #Since space errors are negligible, we take the results from normal5.
        
        SCvv1       = [ 50079:50082 ]
                        
        nb          = 1
        
    end
    
    #Save results for IMA and AMA:
    table           = string("\\begin{tabular}{", 
                        repeat("r", length(vbles)), 
                        "} \n", 
                        "\\hline \n", 
                        header, 
                        "\\hline \n")
    for ii=1:length(SCvv1)
        SCv         = SCvv1[ii]
        table_SC    = TableVbles(SCv, vbles, formats, nb=nb)
        table       = string(table, table_SC, "\\hline \n")
    end
    table           = string(table, "\\end{tabular}")
    print(table)
    
    write("$(FigUbi)IsentropicVortex_TimeAdapt.txt", table)
    
    return
    
end

function NElems_SmoothVortex(; SaveFig::Bool=false, w::Float64=8.50, h::Float64=8.50)

#     SC      = 50027
    SC      = 50085
    nb      = 1
    
    solver  = GetSolver(SC, nb)
    validv  = ones(Bool, length(solver.tv))
    if true
        @. validv  = solver.validv .== 1
    end
    
    PyPlotFigure(w=w, h=h, bottom=1.2, left=1.5)
    semilogy(solver.tv[validv], solver.nElemsv[validv], "k", linewidth=0.5)
    xlabel(L"t^n", fontsize=LabelSize)
    ylabel(latexstring("\\#\\mathcal{T}^{h,n}"), rotation=0, labelpad=-10.0)
    ylim(4000, 5000)
#     yticks([])
#     yticks([175, 180, 185, 190, 195, 200, 205, 210], 
#             ["175", "", "", "", "", "", "", "210"])
    tick_params(axis="both", which="both", labelsize=TickSize)
    if SaveFig
        savefig("$(FigUbi)SmoothVortex_NElems.png", dpi=800, pad_inches=0)
    end
    
    return
    
end

function TableResults_ProfTimes_SmoothVortex(StudyCase::String)

    vbles       = [ "FesOrder",     "TolS_max", "NDOF",         "tmeanCPU",     "EOC_NDOF_tmeanCPU",
                    "tmeanAMA",     "tmeanJm",  "tmeanSCILU",   "tmeanLS",
                    "LS_mean" ]
                    
    formats     = [ "%d",       "%.2E",      "%.0f",     "%.3f",     "%.1f",
                    "%.3f",     "%.3f",     "%.3f",     "%.3f",
                    "%.1f"  ]
    
    header      = string("\$p\$ & \$\\Tol_S\$ & \$\\NDOF\$",
                        " & \\makecell[cr]{Total}", 
                        " & EOC",
                        " & \\makecell[cr]{Mesh \\\\ adapt.}", 
                        " & \\makecell[cr]{Jacobian \\\\ comput.}", 
                        " & \\makecell[cr]{Matrix \\\\ fact.}",
                        " & \\makecell[cr]{Solution of \\\\ linear system}",
                        " & \\makecell[cr]{Linear solver \\\\ iterations}",
                        "\\\\")
                
#     vbles       = [ "FesOrder",     "TolS_max", "NDOF",         "tmeanCPU",     "EOC_NDOF_tmeanCPU",
#                     "%tmeanAMA",    "%tmeanJm", "%tmeanSCILU",  "%tmeanLS",
#                     "LS_mean" ]
#                     
#     formats     = [ "%d",       "%.2E",      "%.0f",     "%.3f",     "%.1f",
#                     "%.1f",     "%.1f",     "%.1f",     "%.1f",
#                     "%.1f"  ]
#     
#     header      = string("\$p\$ & \$\\Tol_S\$ & \$\\NDOF\$", 
#                         " & \\makecell[cr]{Total \\\\ ~[s]}", 
#                         " & EOC",
#                         " & \\makecell[cr]{Mesh \\\\ adaptation}", 
#                         " & \\makecell[cr]{Jacobian \\\\ computation}", 
#                         " & \\makecell[cr]{Matrix \\\\ factorization}",
#                         " & \\makecell[cr]{Solution of \\\\ linear system}",
#                         " & \\makecell[cr]{Linear solver \\\\ iterations}",
#                         "\\\\")
                        
    #------------------------------------------------------------
    
    SCvv1   = NaN
    nb      = NaN
    
    if StudyCase=="normal4" #SizeOrder=FesOrder, AnisoOrder=2, with Eq. loop

                        
        SCvv1       = [ 50050:50054, 
                        50055:50059,
                        50060:50064
                        ]
                        
#         SCvv1       = [ vcat(50052, 50059, 50064) ]

        nb          = 1
        
    elseif StudyCase=="normal5_L2" #with dynamic AMA_rDeltat and correct etaS:
    
        SCvv1       = [ [50083,50084,50095,50085,50086], 
                        50087:50090,
                        50091:50094
                        ]
                        
        nb          = 1
        
    end
    
    #Save results for IMA and AMA:
    table           = string("\\begin{tabular}{", 
                        repeat("r", length(vbles)), 
                        "} \n", 
                        "\\hline \n", 
                        header, 
                        "\\hline \n")
    for ii=1:length(SCvv1)
        SCv         = SCvv1[ii]
        table_SC    = TableVbles(SCv, vbles, formats, nb=nb)
        table       = string(table, table_SC, "\\hline \n")
    end
    table           = string(table, "\\end{tabular}")
    print(table)
    
    write("$(FigUbi)IsentropicVortex_ProfTimes.txt", table)
    
    return
    
end

function TableResults_ProfTimes2_SmoothVortex(StudyCase::String)

    vbles       = [ "FesOrder",     "NDOF",     "tmeanCPU",     
                    "tmeanAMA",     "tmeanJm",  "tmeanSCILU",   "tmeanLS",
                    "LS_mean" ]
                    
    formats     = [ "%d",       "%.0f",     "%.3f",
                    "%.3f",     "%.3f",     "%.3f",     "%.3f",
                    "%.1f"  ]
    
    header      = string("\$p\$ & \$\\NDOF\$ & Total",
                        " & \\makecell[cr]{Mesh \\\\ adaptation}", 
                        " & \\makecell[cr]{Jacobian \\\\ computation}", 
                        " & \\makecell[cr]{Matrix \\\\ factorization}",
                        " & \\makecell[cr]{Solution of \\\\ linear system}",
                        " & \\makecell[cr]{Linear solver \\\\ iterations}",
                        "\\\\")
                        
    #------------------------------------------------------------
    
    SCvv1   = NaN
    nb      = NaN
    
    if StudyCase=="normal4" #SizeOrder=FesOrder, AnisoOrder=2, with Eq. loop
                        
        SCvv1       = [ vcat(50052, 50059, 50064) ]

        nb          = 1
        
    end
    
    #Save results for IMA and AMA:
    table           = string("\\begin{tabular}{", 
                        repeat("r", length(vbles)), 
                        "} \n", 
                        "\\hline \n", 
                        header, 
                        "\\hline \n")
    for ii=1:length(SCvv1)
        SCv         = SCvv1[ii]
        table_SC    = TableVbles(SCv, vbles, formats, nb=nb)
        table       = string(table, table_SC, "\\hline \n")
    end
    table           = string(table, "\\end{tabular}")
    print(table)
    
    write("$(FigUbi)IsentropicVortex_ProfTimes2.txt", table)
    
    return
    
end

function Contour_SmoothVortex(SC::Int, nb::Int; SaveFig::Bool=false, w::Float64=8.50, h::Float64=8.50, 
    PlotVars::Vector{String}=["rho", "v", "p", "M"], mFig::Int=2, nFig::Int=2)

    #Problem data:
    gamma           = 1.4
    
    #Define gas:
    GasModel                = GasIdeal()
    GasModel.gamma          = 1.4
    GasModel.epsilon        = 0.0
    GasModel.nu             = 0.0
    GasModel.beta           = 0.0
    GasModel.kappa_rho_cv   = 0.0
#     GasModel.CSS            = 0.1

    FileName    = GetFileName(SC, nb)
    solver      = GetSolver(SC, nb)
    
    for ii=1:length(PlotVars)
        
        PyPlotFigure(w=w, h=h, bottom=1.5, top=1.0)
        
        #Numerical solution:
        PlotContour(solver, GasModel, PlotVars[ii], delta=1e-5)
        PlotMesh!(SC, nb, color="w")
        title(latexstring(LatexString(PlotVars[ii]),"; t^n=", sprintf1("%.2e", solver.t)),
            fontsize=10)
        tick_params(axis="both", which="both", labelsize=TickSize)
        axis("off")
        
        if SaveFig
            savefig("$(FigUbi)SC$(SC)_Contour_$(PlotVars[ii]).png", dpi=800, pad_inches=0)
        end
    
    end
    
    return
    
end

function Video_SmoothVortex(SC::Int; SaveFig::Bool=false, w::Float64=8.50, h::Float64=w*14/24+0.75, 
    PlotVars::Vector{String}=["rho", "v", "p", "M"], mFig::Int=2, nFig::Int=2)

    #Problem data:
    gamma           = 1.4
    
    #Define gas:
    GasModel                = GasIdeal()
    GasModel.gamma          = 1.4
    GasModel.epsilon        = 0.0
    GasModel.nu             = 0.0
    GasModel.beta           = 0.0
    GasModel.kappa_rho_cv   = 0.0
#     GasModel.CSS            = 0.1

    fig     = PyPlotSubPlots(mFig, nFig, w=w, h=h, top=0.75, left=0.01, right=0.01, bottom=0.01)
    
    kSave   = 0
    for nb=0:49
    
        FileName    = GetFileName(SC, nb)
        solver      = GetSolver(SC, nb)
    
        for ii=1:length(PlotVars)
            
            PyPlot.subplot(mFig,nFig,ii)
            PyPlot.cla()
            
            PlotContour(solver, GasModel, PlotVars[ii], delta=1e-5)
            PlotMesh!(SC, nb, color="w")
   
            title(latexstring(LatexString(PlotVars[ii]),"; t^n=", sprintf1("%.2e", solver.t)),
                fontsize=10)
            tick_params(axis="both", which="both", labelsize=TickSize)
            axis("off")
            
        end
        
        if SaveFig
            savefig("$(FigUbi)SC$(SC)_$(kSave).png", dpi=800, pad_inches=0)
        end
        kSave       += 1
        pause(0.1)
        
    end
    
    return
    
end

function Errors_SmoothVortex(; SaveFig::Bool=false, w::Float64=7.50, h::Float64=6.50)

    SC      = 50061
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
            latexstring("\\mathrm{Tol}_{S}^{\\max}"), 
            latexstring("\\mathrm{Tol}_{T}"),
            latexstring("\\mathrm{Rejected~steps}")], 
            loc=(0.01,0.42), 
            fontsize=TickSize)
    
    tick_params(axis="both", which="both", labelsize=TickSize)
    
#     xlim(5.0, 7.0)
#     ylim(0.9e-8, 1.2*solver.TolT)
    xlabel(L"t")
    
    if SaveFig
        savefig("$(FigUbi)SmoothVortex_Errors.png", dpi=800, pad_inches=0)
    end
    
    return
    
end


    
