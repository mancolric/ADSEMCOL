include("PlotResults.jl")
include("../test_CompressibleFlow/test_CompressibleFlow.jl")

function Contour_TriplePoint(SC::Int, nb::Int; SaveFig::Bool=false, w::Float64=8.50, h::Float64=w*3/7, 
    PlotVars::Vector{String}=["rho", "v", "e"], mFig::Int=1, nFig::Int=3)
    
    #Define gas:
    GasModel                = GasIdeal()
    GasModel.gamma          = 1.4
    GasModel.epsilon        = 0.0
    GasModel.nu             = 0.0
    GasModel.beta           = 0.0
    GasModel.kappa_rho_cv   = 0.0

    FileName    = GetFileName(SC, nb)
    solver      = GetSolver(SC, nb)
    
    for ii=1:length(PlotVars)
        
        
        #Numerical solution:
        PyPlotFigure(w=w, h=h, bottom=0.0, top=0.0, left=0.0, right=0.0)
        v_plot  = PlotContour(solver, GasModel, PlotVars[ii], delta=1e-5)
#         title(latexstring(LatexString(PlotVars[ii]),"; t^n=", sprintf1("%.2e", solver.t)),
#             fontsize=10)
        if PlotVars[ii]=="e"
            text(-0.5, 0.9, latexstring(LatexString(PlotVars[ii])), color="k", fontsize=10)
        else
            text(-0.5, 0.9, latexstring(LatexString(PlotVars[ii])), color="w", fontsize=10)
        end
        tick_params(axis="both", which="both", labelsize=TickSize)
        axis("off")
        axis("equal")
        if SaveFig
            savefig("$(FigUbi)SC$(SC)_Contour_$(PlotVars[ii]).png", dpi=800, pad_inches=0)
        end
    
        #Colorbar (Careful: if the colorbar is small, visualization in screen is different 
        #from saved figure).
        PyPlotFigure(w=1.2, h=h, bottom=0.2, top=0.2, left=0.9, right=0.0)
        PyPlotColorbar!(max(0.0,minimum(v_plot)-1e-10), maximum(v_plot)+1e-10, 4, tickformat="%.2f")
        if SaveFig
            savefig("$(FigUbi)SC$(SC)_Colorbar_$(PlotVars[ii]).png", dpi=800, pad_inches=0)
        end
        
        println(PlotVars[ii], ": min=", minimum(v_plot), ", max=", maximum(v_plot))
        
    end
    
    return
    
end

function Contour_TriplePoint2(SC::Int, nb::Int; SaveFig::Bool=false, w::Float64=8.50, h::Float64=w*3/7, 
    PlotVars::Vector{String}=["rho", "v", "e"], mFig::Int=1, nFig::Int=3)
    
    #Define gas:
    GasModel                = GasIdeal()
    GasModel.gamma          = 1.4
    GasModel.epsilon        = 0.0
    GasModel.nu             = 0.0
    GasModel.beta           = 0.0
    GasModel.kappa_rho_cv   = 0.0

    FileName    = GetFileName(SC, nb)
    solver      = GetSolver(SC, nb)
            
    for ii=1:length(PlotVars)
        
        figure()
        v_plot  = PlotNodes3D(solver, GasModel, PlotVars[ii])
    
        #Numerical solution:
        PyPlotFigure(w=w, h=h+0.7, bottom=0.0, top=0.7, left=0.0, right=0.0)
        title(latexstring(LatexString(PlotVars[ii]),"; t^n=", sprintf1("%.2e", solver.t)),
            fontsize=10)
        if PlotVars[ii]=="e"
            v_plot  = PlotContour(solver, GasModel, PlotVars[ii], delta=1e-5, z1=0.24)
        else
            v_plot  = PlotContour(solver, GasModel, PlotVars[ii], delta=1e-5)
        end
        tick_params(axis="both", which="both", labelsize=TickSize)
        axis("off")
        axis("equal")
        if SaveFig
            savefig("$(FigUbi)SC$(SC)_$(nb)_Contour_$(PlotVars[ii]).png", dpi=800, pad_inches=0)
        end
    
        #Colorbar (Careful: if the colorbar is small, visualization in screen is different 
        #from saved figure).
        PyPlotFigure(w=1.2, h=h+0.7, bottom=0.2, top=0.9, left=0.9, right=0.0)
        if PlotVars[ii]=="e"
            PyPlotColorbar!(0.24, maximum(v_plot)+1e-10, 4, tickformat="%.2f")
        else
            PyPlotColorbar!(max(0.0,minimum(v_plot)-1e-10), maximum(v_plot)+1e-10, 4, tickformat="%.2f")
        end
        if SaveFig
            savefig("$(FigUbi)SC$(SC)_$(nb)_Colorbar_$(PlotVars[ii]).png", dpi=800, pad_inches=0)
        end
#         title("", fontsize=10)
            
        println(PlotVars[ii], ": min=", minimum(v_plot), ", max=", maximum(v_plot))
        
    end
    
    return
    
end

function Mesh_TriplePoint(SC::Int, nb::Int; SaveFig::Bool=false, w::Float64=8.50, h::Float64=w*3/7)


    FileName    = GetFileName(SC, nb)
    solver      = GetSolver(SC, nb)
    
    PyPlotFigure(w=w, h=h, bottom=0.0, top=0.0, left=0.0, right=0.0)
    PlotContour(solver.u[1], solver.fes)
    PlotMesh!(SC, nb, color="w")
    tick_params(axis="both", which="both", labelsize=TickSize)
    axis("off")
    axis("equal")
    if SaveFig
        savefig("$(FigUbi)SC$(SC)_Mesh.png", dpi=800, pad_inches=0)
    end
    
#     xlimv   = [3.0, 5.5]
#     ylimv   = [-1.5, 1.0]
    xlimv   = [5.3, 5.8]
    ylimv   = [-1.25, -0.75]
    hZoom   = h
    wZoom   = hZoom*(xlimv[2]-xlimv[1])/(ylimv[2]-ylimv[1])
    #
    PyPlotFigure(w=wZoom, h=hZoom, bottom=0.0, top=0.0, left=0.0, right=0.0)
    PlotContour(solver.u[1], solver.fes)
    PlotMesh!(SC, nb, color="w")
    xlim(xlimv)
    ylim(ylimv)
    tick_params(axis="both", which="both", labelsize=TickSize)
    axis("off")
    if SaveFig
        savefig("$(FigUbi)SC$(SC)_MeshZoom.png", dpi=800, pad_inches=0)
    end
    
    return
    
end

function Video_TriplePoint_old(SC::Int; SaveFig::Bool=false, w::Float64=8.50, h::Float64=8.50)

    fig     = PyPlotSubPlots(1, 2, w=w, h=h, top=1.0)
    
    kSave   = 0
    for nb=vcat(0:38,1000)
    
        FileName    = GetFileName(SC, nb)
        solver      = GetSolver(SC, nb)
    
        PyPlot.subplot(1,2,1)
        PyPlot.cla()
        PlotContour(solver.u[1], solver.fes)
#         PlotMesh!(SC, nb, color="w")
        tick_params(axis="both", which="both", labelsize=TickSize)
        title(latexstring(LatexString("rho"),"; t^n=", sprintf1("%.2e", solver.t)),
            fontsize=10)
        
        PyPlot.subplot(1,2,2)
        PyPlot.cla()
        PlotContour(solver.u[1], solver.fes)
        PlotMesh!(SC, nb, color="w")
        tick_params(axis="both", which="both", labelsize=TickSize)
        title(latexstring(LatexString("rho"),"; t^n=", sprintf1("%.2e", solver.t)),
            fontsize=10)
        
        if SaveFig
            savefig("$(FigUbi)SC$(SC)_$(kSave).png", dpi=800, pad_inches=0)
        end
        kSave       += 1
        
        pause(1.0)
        
    end
    
    return
    
end

function Video_TriplePoint(SC::Int; SaveFig::Bool=false, w::Float64=8.50, h::Float64=w*3/7+0.75, 
    PlotVars::Vector{String}=["rho", "rho", "v", "e"], mFig::Int=2, nFig::Int=2)

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
    for nb=0:50
    
        FileName    = GetFileName(SC, nb)
        solver      = GetSolver(SC, nb)
    
        for ii=1:length(PlotVars)
            
            PyPlot.subplot(mFig,nFig,ii)
            PyPlot.cla()
            
            PlotContour(solver, GasModel, PlotVars[ii], delta=1e-5)
            if ii==1
                PlotMesh!(SC, nb, color="w")
            end
            
            title(latexstring(LatexString(PlotVars[ii]),"; t^n=", sprintf1("%.2e", solver.t)),
                fontsize=10)
            tick_params(axis="both", which="both", labelsize=TickSize)
            axis("off")
            
        end
        
        if SaveFig
            savefig("$(VideosUbi)/TriplePoint/SC$(SC)_$(kSave).png", dpi=800, pad_inches=0)
        end
        kSave       += 1
        pause(0.1)
            
    end
    
    return
    
end

function TableResults_ProfTimes2_TriplePoint(StudyCase::String)

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
    
    if StudyCase=="normal" #SizeOrder=FesOrder, AnisoOrder=2, with Eq. loop
                        
        SCvv1       = [ [80005, 80004] ]

        nb          = 1
        
    elseif StudyCase=="normal2" #epsilon=1e-4
                        
        SCvv1       = [ [80015, 80018] ]

        nb          = 50
        
    elseif StudyCase=="normal3" #epsilon=3e-5
                        
        SCvv1       = [ [80028, 80029] ]

        nb          = 7
        
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
    
#     write("$(FigUbi)TriplePoint_ProfTimes2.txt", table)
    
    return
    
end
