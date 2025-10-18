include("PlotResults.jl")

function Contour_KPP(SC::Int, nb::Int; SaveFig::Bool=false, w::Float64=8.50, h::Float64=8.50)

    solver      = GetSolver(SC, nb)
    
    PyPlotFigure(w=w, h=h, top=0.0, bottom=0.0, left=0.0, right=0.0)
#     uP          = GetPolynomials(solver.u, solver.fes)
#     v_plot      = PlotContour(uP[1], solver.fes.PSpace, cmap="jet")
    v_plot      = PlotContour(solver.u[1], solver.fes)
    axis("off")
    axis("equal")
    if SaveFig
        savefig("$(FigUbi)SC$(SC).png", dpi=800, pad_inches=0)
    end
    
    PyPlotFigure(w=1.4, h=h, bottom=0.2, top=0.2, left=1.1, right=0.0)
    PyPlotColorbar!(minimum(v_plot), maximum(v_plot), 4, tickformat="%.2f")
    if SaveFig
        savefig("$(FigUbi)SC$(SC)_Colorbar.png", dpi=800, pad_inches=0)
    end
    
    println("u min = ", minimum(v_plot))
    println("u max = ", maximum(v_plot))
    uP          = GetPolynomials(solver.u, solver.fes)
    println("u min = ", minimum(uP[1]))
    println("u max = ", maximum(uP[1]))
    #The deviation of the two first values is usually smaller because of the resolution
    #of the contour plot. Try a larger value of nSubtriangles to have similar deviations.
    
    return
    
end

function Mesh_KPP(SC::Int, nb::Int; SaveFig::Bool=false, w::Float64=8.50, h::Float64=8.50)

    solver      = GetSolver(SC, nb)
    
    PyPlotFigure(w=w, h=h, top=0.0, bottom=0.0, left=0.0, right=0.0)
    PlotContour(solver.u[1], solver.fes, cmap="jet")
    PlotMesh!(solver.mesh, color="w") 
    axis("off")
    if SaveFig
        savefig("$(FigUbi)SC$(SC)_Mesh.png", dpi=800, pad_inches=0)
    end
    
    PyPlotFigure(w=w, h=h, top=0.0, bottom=0.0, left=0.0, right=0.0)
    PlotContour(solver.u[1], solver.fes, cmap="jet")
    PlotMesh!(solver.mesh, color="w") 
    axis("off")
    xlim(-0.7,-0.3)
    ylim(0.6,1.0)
    if SaveFig
        savefig("$(FigUbi)SC$(SC)_MeshZoom.png", dpi=800, pad_inches=0)
    end
    
    PyPlotFigure(w=w, h=h, top=0.0, bottom=0.0, left=0.0, right=0.0)
    PlotContour(solver.u[1], solver.fes, cmap="jet")
    PlotMesh!(solver.mesh, color="w") 
    xlim(-0.5-0.001,-0.5+0.001)
    ylim(0.7395-0.001,0.7395+0.001)
    axis("off")
    if SaveFig
        savefig("$(FigUbi)SC$(SC)_MeshZoom2.png", dpi=800, pad_inches=0)
    end
    
    return
    
end

function CFL_KPP(; SaveFig::Bool=false, w::Float64=8.50, h::Float64=8.50)

    SC      = 40001
    nb      = 1
    
    solver  = GetSolver(SC, nb)
    validv  = ones(Bool, length(solver.tv))
    if true
        @. validv  = solver.validv .== 1
    end
    
    PyPlotFigure(w=w, h=h, bottom=1.5, left=1.5)
    semilogy(solver.tv[validv], solver.CFLv[validv], "k", linewidth=0.5)
    xlabel(L"t")
    ylabel(L"CFL", rotation=0)
    if SaveFig
        savefig("$(FigUbi)SC$(SC)_NElems.png", dpi=800, pad_inches=0)
    end
    
    return
    
end

function Errors_KPP(SC::Int, nb::Int; SaveFig::Bool=false, w::Float64=7.50, h::Float64=6.50)

#     SC      = 40010
#     nb      = 1
    solver  = GetSolver(SC, nb)
    
    #Erase t=0:
    aux                 = @. solver.tv == 0.0
    solver.tv[aux]      .= NaN
    solver.etaSv[aux]   .= NaN
    solver.etaTv[aux]   .= NaN
    
    fig     = PyPlotFigure(w=w+3.5, h=h, bottom=1.0, left=1.0, right=3.5)
    markv   = @. solver.validv == 0 && solver.tv>0.0
    validv  = @. solver.validv == 1
    
    #all steps:
    semilogy(solver.tv, solver.etaSv, "b", linewidth=0.7, marker="x", markersize=1.0)
    semilogy(solver.tv, solver.etaTv, "g", linewidth=0.7, marker="x", markersize=1.0)
    semilogy(solver.tv, solver.etaAv, "c", linewidth=0.7, marker="x", markersize=1.0)
    
    #=
    #accepted time steps:
    semilogy(solver.tv[validv], solver.etaSv[validv], "b", linewidth=0.7, marker="x", markersize=1.0)
    semilogy(solver.tv[validv], solver.etaTv[validv], "g", linewidth=0.7, marker="x", markersize=1.0)
    =#
    
    #tolerances:
    semilogy(solver.tv, 0.0.*solver.tv .+ solver.TolS_max, "--b", linewidth=0.7)
    semilogy(solver.tv, 0.0.*solver.tv .+ solver.TolT, "--g", linewidth=0.7)
    
    #rejected steps:
    semilogy(solver.tv[markv], solver.etaSv[markv], "r", linestyle="none", marker="s", markersize=1.5)
    semilogy(solver.tv[markv], solver.etaTv[markv], "r", linestyle="none", marker="s", markersize=1.5)
    
    legend([latexstring("\\mathcal{E}_{S},~\\mathrm{all~steps}"), 
            latexstring("\\mathcal{E}_{T},~\\mathrm{all~steps}"),
            latexstring("\\mathcal{E}_{A},~\\mathrm{all~steps}"), 
            latexstring("\\mathrm{Tol}_{S}"),
            latexstring("\\mathrm{Tol}_{T}"),
            latexstring("\\mathrm{Rejected~steps}")
            ],
            loc=(0.01,0.42), 
            bbox_to_anchor=(1.05, 0.3),
            fontsize=TickSize)
    
    #=
    legend([latexstring("\\mathcal{E}_{S},"), 
            latexstring("\\mathcal{E}_{T},"), 
            latexstring("\\mathrm{Tol}_{S}"), 
            latexstring("\\mathrm{Tol}_{T}"),
            ],
            loc=(0.01,0.42), 
            fontsize=TickSize)
    =#
    
    tick_params(axis="both", which="both", labelsize=TickSize)
    
    xlim(0.01,0.02)
    ylim(1e-9, 10*solver.TolT)
    xlabel(L"t")
    
#     xlim(0.0,0.005)
#     ylim(1e-9, 10*solver.TolT)
#     xlabel(L"t")
    
    if SaveFig
        savefig("$(FigUbi)KPP_Errors.png", dpi=800, pad_inches=0)
    end
    
    return
    
end

function ErrorsValid_KPP(SC::Int, nb::Int; SaveFig::Bool=false, w::Float64=7.50, h::Float64=6.50)

    solver  = GetSolver(SC, nb)
    
    #Erase t=0:
    aux                 = @. solver.tv == 0.0
    solver.tv[aux]      .= NaN
    solver.etaSv[aux]   .= NaN
    solver.etaTv[aux]   .= NaN
    
    fig     = PyPlotFigure(w=w+2.5, h=h, bottom=1.0, left=1.0, right=2.5)
    markv   = @. solver.validv == 0 && solver.tv>0.0
    validv  = @. solver.validv == 1 && solver.tv>0.0
    
    #accepted time steps:
    semilogy(solver.tv[validv], solver.etaSv[validv], "b", linewidth=0.7, marker="x", markersize=1.0)
    semilogy(solver.tv[validv], solver.etaTv[validv], "g", linewidth=0.7, marker="x", markersize=1.0)
#     semilogy(solver.tv[validv], solver.etaAv[validv], "r", linewidth=0.7, marker="x", markersize=1.0)
    
    #tolerances:
    semilogy(solver.tv, 0.0.*solver.tv .+ solver.TolS_max, "--b", linewidth=0.7)
    semilogy(solver.tv, 0.0.*solver.tv .+ solver.TolT, "--g", linewidth=0.7)
    
    legend([latexstring("\\mathcal{E}_{S},"), 
            latexstring("\\mathcal{E}_{T},"),
#             latexstring("\\mathcal{E}_{A},"), 
            latexstring("\\mathrm{Tol}_{S}"), 
            latexstring("\\mathrm{Tol}_{T}"),
            ],
            loc=(0.01,0.42), 
            bbox_to_anchor=(1.05, 0.3),
            fontsize=TickSize)
    
    tick_params(axis="both", which="both", labelsize=TickSize)
    
#     xlim(0.0005,0.0045)
    xlim(0.0, 0.005)
#     ylim(1e-9, 10*solver.TolT)
    ylim(1e-7, 10*solver.TolT)
    xlabel(L"t")
    
    if SaveFig
        savefig("$(FigUbi)KPP_Errors.png", dpi=800, pad_inches=0)
    end
    
    return
    
end

function Video_KPP(SC::Int; SaveFig::Bool=false, w::Float64=8.50, h::Float64=8.50)

#     fig     = PyPlotSubPlots(1, 2, w=w, h=h, top=1.0)
    fig     = PyPlotSubPlots(1, 2, w=w, h=h, left=0.1, right=0.1, bottom=0.1, top=0.7)
    
    kSave   = 0
    for nb=0:49
    
        FileName    = GetFileName(SC, nb)
        solver      = GetSolver(SC, nb)
    
        PyPlot.subplot(1,2,1)
        PyPlot.cla()
        PlotContour(solver.u[1], solver.fes)
#         PlotMesh!(SC, nb, color="w")
        tick_params(axis="both", which="both", labelsize=TickSize)
        axis("off")
#         title(latexstring("t^n=", sprintf1("%.2e", solver.t)),
#             fontsize=10)
        
        PyPlot.subplot(1,2,2)
        PyPlot.cla()
        PlotContour(solver.u[1], solver.fes)
        PlotMesh!(SC, nb, color="w")
        tick_params(axis="both", which="both", labelsize=TickSize)
        axis("off")
        suptitle(latexstring("t^n=", sprintf1("%.2e", solver.t)),
            fontsize=10)
        
        if SaveFig
            savefig("$(VideosUbi)KPP/SC$(SC)_$(kSave).png", dpi=800, pad_inches=0)
        end
        kSave       += 1
        
        pause(0.1)
        
    end
    
    return
    
end

function MinMax_KPP(SC::Int, nb::Int; 
    SaveFig::Bool=false, w::Float64=8.0, h::Float64=6.0)

    FileName    = GetFileName(SC, nb)
    tv          = load(FileName, "tv")
    uminv       = load(FileName, "uminv")
    umaxv       = load(FileName, "umaxv")
    Deltav      = @. max(abs.(umaxv.-3.5*pi), abs.(uminv.-0.25*pi))
    
    PyPlotFigure(w=w, h=h, bottom=1.0, left=1.0, right=0.2, top=0.2)
#     semilogy(tv, abs.(umaxv.-3.5*pi), "b")
#     semilogy(tv, abs.(uminv.-0.25*pi), "g")
    semilogy(tv, Deltav, "b")
    tick_params(axis="both", which="both", labelsize=TickSize)
    xlabel(L"t^n", fontsize=LabelSize)
    leg         = vcat(latexstring("|\\max ~ u^{h\\tau,n-}-1|"), 
                        latexstring("|\\min ~ u^{h\\tau,n-}|"))
#     ylim(1e-6,1e-2)
    legend(leg, fontsize=TickSize) 
    if SaveFig
        savefig("$(FigUbi)KPP_Oscillations.png", dpi=800, pad_inches=0)
    end
    
    Nt          = length(uminv)
    println("Initial condition: Delta u_min=", abs(uminv[1]-0.25*pi), ", Delta u_max=", abs(umaxv[1]-3.5*pi))
    println("Final solution: Delta u_min=", abs(uminv[Nt]-0.25*pi), ", Delta u_max=", abs(umaxv[Nt]-3.5*pi))
    
    return
    
end

function OscillationIndicator_KPP(StudyCase::String; 
    SaveFig::Bool=false, w::Float64=8.1, h::Float64=6.0)

    SCv         = []
    nb          = NaN
    
    if StudyCase=="TolS_1e-4"
        SCv         = 40020:40022
        nb          = 1
    elseif StudyCase=="TolS_1e-5"
        SCv         = 40013:40015 #40016 is also available
        nb          = 1
    end
    
    nSC         = length(SCv)
    colorv      = PyPlotColors("jet2", nSC)
    
    fig1        = PyPlotFigure(w=w, h=h, bottom=1.0, left=1.5, right=0.1, top=0.2)
    leg         = String[]
    for ii=1:nSC
    
        SC          = SCv[ii]
        FileName    = GetFileName(SC, nb)
        tv          = load(FileName, "tv")
        uminv       = load(FileName, "uminv")
        umaxv       = load(FileName, "umaxv")
        Deltav      = @. max(abs.(umaxv.-3.5*pi), abs.(uminv.-0.25*pi))
        solver      = GetSolver(SC,nb)
        TolT        = solver.TolT
        
        figure(fig1.number)
        semilogy(tv, Deltav, color=colorv[ii])
    
        leg         = vcat(leg, latexstring("\\mathrm{Tol}_T=", sprintf1("%.2E", TolT)))
        
    end
    
    figure(fig1.number)
    xlabel(L"t", fontsize=LabelSize)
    ylabel(latexstring("\\Delta_{ST}"), fontsize=LabelSize, rotation=0, labelpad=5.0)
#     ylim(1e-6,1e-2)
    legend(leg, fontsize=TickSize, labelspacing=0.1) 
    tick_params(axis="both", which="both", labelsize=TickSize)
    if SaveFig
        savefig("$(FigUbi)KPP_OscIndicator_$(StudyCase).png", dpi=800, pad_inches=0)
    end
    
    return
    
end

function TableResults_PrePost_1_KPP()

    vbles       = [ "TolT", 
                    "Nt",       "Nt_space_rejected", 
                    "NDOF",     "tCPU" ]
    formats     = [ "%.2E",     
                    "%d",       "%d", 
                    "%d",       "%.1f" ]
    
    #------------------------------------------------------------
    
    SCv0    = [ 40032, 40034, 40036 ]
    SCv1    = [ 40033, 40035, 40037 ]
    nb      = 4
    
    #Save results for pre- and post-time-step adaptaton approaches:
#     table           = string("\\begin{tabular}{", repeat("r", length(vbles)+2), "} \n", 
#                         "\\hline \n", 
#                         "\\makecell[cr]{AMA \\\\ type} ", 
#                         " & \$\\Tol_T\$ ",
#                         " & \\makecell[cr]{Accepted \\\\ time steps} ", 
#                         " & \\makecell[cr]{Rejected \\\\ time steps} ", 
#                         " & NDOF ",
#                         " & \$ t_{\\mathrm{CPU}}[s] \$ ", 
#                         " & \\makecell[cr]{\$ t_{\\mathrm{CPU}} \$ \\\\ reduction} ", 
#                         "\\\\ \\hline \n")
#     table           = string("\\begin{tabular}{", repeat("r", length(vbles)+2), "} \n", 
#                         "\\hline \n", 
#                         "\\makecell[cr]{AMA \\\\ type} ", 
#                         " & \$\\Tol_T\$ ",
#                         " & \\makecell[cr]{Accepted \\\\ time steps} ", 
#                         " & \\makecell[cr]{Rejected \\\\ time steps \\\\ (\$ \\Ecal^{S,(n+1)-}<\\Tol_S \$) } ", 
#                         " & NDOF ",
#                         " & \$ t_{\\mathrm{CPU}}[s] \$ ", 
#                         " & \\makecell[cr]{\$ t_{\\mathrm{CPU}} \$ \\\\ reduction} ", 
#                         "\\\\ \\hline \n")
    table           = string("\\begin{tabular}{", repeat("r", length(vbles)+2), "} \n", 
                        "\\hline \n", 
                        "\\makecell[cr]{AMA \\\\ type} ", 
                        " & \$ \\Tol_T\$ ",
                        " & \$ N_t \$ ",
                        " & \$ N_t^{\\mathrm{rej},S} \$ ",
                        " & NDOF ",
                        " & \$ t_{\\mathrm{CPU}}[s] \$ ", 
                        " & \\makecell[cr]{\$ t_{\\mathrm{CPU}} \$ \\\\ reduction} ", 
                        "\\\\ \\hline \n")
    
    #Plot table:
    for ii=1:length(SCv0)
        
        #Get just one row:
        table_SC1       = TableVbles([SCv1[ii]], vbles, formats, nb=nb, LastSeparator="\t")
        table_SC0       = TableVbles([SCv0[ii]], vbles, formats, nb=nb, LastSeparator="\t")
        
        #Get CPU time reduction:
        tCPU1           = GetVbles(SCv1[ii], nb, ["tCPU"])[1]
        tCPU0           = GetVbles(SCv0[ii], nb, ["tCPU"])[1]
        tCPU_red        = (tCPU1-tCPU0)/tCPU1
        table           = string(table, "Post \t & \t", 
                                        table_SC1, 
                                        "\t & - \t \\\\ \n", 
                                        "Pre  \t & \t", 
                                        table_SC0, 
                                        "\t & ", sprintf1("%2.1f", tCPU_red*100), "\\%", 
                                        "\t \\\\ \\hline \n")
    end
    table           = string(table, "\\end{tabular}")
    print(table)
    
    write("$(FigUbi)KPP_PrePost.txt", table)
    
    return
    
end

function TableResults_PrePost_2_KPP()

    vbles       = [ "TolT", 
                    "Nt",       "Nt_rejected",      "Nt_space_rejected", 
                    "NDOF",     "tCPU" ]
    formats     = [ "%.2E",     
                    "%d",       "%d",               "%d", 
                    "%d",       "%.1f" ]
    
    #------------------------------------------------------------
    
    SCv0    = [ 40032, 40034, 40036 ]
    SCv1    = [ 40033, 40035, 40037 ]
    nb      = 4
    
    #Save results for pre- and post-time-step adaptaton approaches:
    table           = string("\\begin{tabular}{", repeat("r", length(vbles)+2), "} \n", 
                        "\\hline \n", 
                        "\\makecell[cr]{AMA \\\\ type} ", 
                        " & \$ \\Tol_T\$ ",
                        " & \$ N_t \$ ",
                        " & \$ N_t^{\\mathrm{rej}} \$ ",
                        " & \$ N_t^{\\mathrm{rej},S} \$ ",
                        " & NDOF ",
                        " & \$ t_{\\mathrm{CPU}}[s] \$ ", 
                        " & \\makecell[cr]{\$ t_{\\mathrm{CPU}} \$ \\\\ reduction} ", 
                        "\\\\ \\hline \n")
    
    #Plot table:
    for ii=1:length(SCv0)
        
        #Get just one row:
        table_SC1       = TableVbles([SCv1[ii]], vbles, formats, nb=nb, LastSeparator="\t")
        table_SC0       = TableVbles([SCv0[ii]], vbles, formats, nb=nb, LastSeparator="\t")
        
        #Get CPU time reduction:
        tCPU1           = GetVbles(SCv1[ii], nb, ["tCPU"])[1]
        tCPU0           = GetVbles(SCv0[ii], nb, ["tCPU"])[1]
        tCPU_red        = (tCPU1-tCPU0)/tCPU1
        table           = string(table, "Post \t & \t", 
                                        table_SC1, 
                                        "\t & - \t \\\\ \n", 
                                        "Pre  \t & \t", 
                                        table_SC0, 
                                        "\t & ", sprintf1("%2.1f", tCPU_red*100), "\\%", 
                                        "\t \\\\ \\hline \n")
    end
    table           = string(table, "\\end{tabular}")
    print(table)
    
    write("$(FigUbi)KPP_PrePost.txt", table)
    
    return
    
end

function TableResults_PrePost_KPP()

    vbles       = [ "Nt",       "Nt_space_rejected", 
                    "NDOF",     "tCPU" ]
    formats     = [ "%d",       "%d", 
                    "%d",       "%.1f" ]
    
    #------------------------------------------------------------
    
    SCv0    = [ 40032, 40034, 40036 ]
    SCv1    = [ 40033, 40035, 40037 ]
    nb      = 10
    
    #Save results for pre- and post-time-step adaptation approaches:
    table           = string("\\begin{tabular}{", repeat("r", length(vbles)+3), "} \n", 
                        "\\hline \n", 
                        "\$ \\Tol_T\$ ",
                        " & \\makecell[cr]{AMA \\\\ type} ", 
                        " & \$ N_t \$ ",
                        " & \$ N_t^{\\mathrm{rej},S} \$ ",
                        " & NDOF ",
                        " & \$ t_{\\mathrm{CPU}}[s] \$ ", 
                        " & \\makecell[cr]{\$ t_{\\mathrm{CPU}} \$ \\\\ reduction} ", 
                        "\\\\ \\hline \n")
    
    #Plot table:
    for ii=1:length(SCv0)
        
        #Get just one row:
        table_SC1       = TableVbles([SCv1[ii]], vbles, formats, nb=nb, LastSeparator="\t")
        table_SC0       = TableVbles([SCv0[ii]], vbles, formats, nb=nb, LastSeparator="\t")
        
        #Get CPU time reduction:
        tCPU1           = GetVbles(SCv1[ii], nb, ["tCPU"])[1]
        tCPU0           = GetVbles(SCv0[ii], nb, ["tCPU"])[1]
        TolT1           = GetVbles(SCv1[ii], nb, ["TolT"])[1]
        TolT0           = GetVbles(SCv0[ii], nb, ["TolT"])[1]
        tCPU_red        = (tCPU1-tCPU0)/tCPU1
        table           = string(table, sprintf1("%.2E", TolT1), " \t & \t ", 
                                        "Post \t & \t", 
                                        table_SC1, 
                                        "\t & - \t \\\\ \n", 
                                        sprintf1("%.2E", TolT0), " \t & \t ", 
                                        "Pre  \t & \t", 
                                        table_SC0, 
                                        "\t & ", sprintf1("%2.1f", tCPU_red*100), "\\%", 
                                        "\t \\\\ \\hline \n")
    end
    table           = string(table, "\\end{tabular}")
    print(table)
    
    write("$(FigUbi)KPP_PrePost.txt", table)
    
    return
    
end

function SpaceErrors_KPP(SC::Int, nb::Int; 
    SaveFig::Bool=false, w::Float64=7.50, h::Float64=6.50, PlotLeg::Bool=false)

    solver  = GetSolver(SC, nb)
    
    #Erase t=0:
    aux                 = @. solver.tv == 0.0
    solver.tv[aux]      .= NaN
    solver.etaSv[aux]   .= NaN
    solver.etaTv[aux]   .= NaN
    
    #Deltat_i = t_i - t_(i-1)
    Deltatv         = vcat(0.0, diff(solver.tv)) 
    
    fig     = PyPlotFigure(w=w, h=h, bottom=1.0, left=1.0)
    
    #all steps:
    markv   = @. Deltatv!=0.0
    semilogy(solver.tv[markv], solver.etaSv[markv], "b", linewidth=0.7, marker="s", markersize=3.5)
    
    #Valid time levels:
    markv   = @. solver.validv == 1 && Deltatv!=0.0
    semilogy(solver.tv[markv], solver.etaSv[markv], "c", linewidth=0.7, marker="x", markersize=3.5)
    
    #tolerances:
    semilogy(solver.tv, 0.0.*solver.tv .+ solver.TolS_max, "--b", linewidth=0.7)
    
    #rejected steps:
    markv   = @. solver.validv == 0 && Deltatv!=0.0
    semilogy(solver.tv[markv], solver.etaSv[markv], "r", linestyle="none", marker="x", markersize=3.5)
    
    if PlotLeg
        legend([latexstring("\\mathcal{E}_{S},~\\mathrm{all~steps}"), 
                latexstring("\\mathcal{E}_{S},~\\mathrm{accepted~steps}"), 
                latexstring("\\mathrm{Tol}_{S}"), 
                latexstring("\\mathcal{E}_{S},~\\mathrm{rejected~steps}")],  
                loc="best", 
                fontsize=TickSize)
    end
    
    tick_params(axis="both", which="both", labelsize=TickSize)
    
    xlim(0.82, 0.83)
    ylim(8e-7, 2*solver.TolS_max)
    xlabel(L"t")
    
    if SaveFig
        savefig("$(FigUbi)KPP_SC$(SC)_Errors.png", dpi=800, pad_inches=0)
    end
    
    return
    
end
