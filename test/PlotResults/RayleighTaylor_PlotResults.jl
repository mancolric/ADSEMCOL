include("PlotResults.jl")

function Contour_RayleighTaylor(SC::Int, nb::Int; SaveFig::Bool=false, w::Float64=8.50, h::Float64=8.50, 
    PlotVars::Vector{String}=["rho", "v", "p", "M"], mFig::Int=1, nFig::Int=3)
    
    #Define gas:
    GasModel                = GasIdeal()
    GasModel.gamma          = 5/3
    GasModel.epsilon        = 0.0
    GasModel.nu             = 0.0
    GasModel.beta           = 0.0
    GasModel.kappa_rho_cv   = 0.0

    FileName    = GetFileName(SC, nb)
    solver      = GetSolver(SC, nb)
    
    for ii=1:length(PlotVars)
        
        PyPlotFigure(w=w, h=h, bottom=1.5, top=1.0)
        
        #Numerical solution:
        PlotContour(solver, GasModel, PlotVars[ii], delta=1e-5)
        title(latexstring(LatexString(GasModel, PlotVars[ii]),"; t^n=", sprintf1("%.2e", solver.t)),
            fontsize=10)
        tick_params(axis="both", which="both", labelsize=TickSize)
        ylim(-0.25, 0.25)
        
        if SaveFig
            savefig("$(FigUbi)SC$(SC)_Contour_$(PlotVars[ii]).png", dpi=800, pad_inches=0)
        end
    
    end
    
    return
    
end

function Mesh_RayleighTaylor(SC::Int, nb::Int; SaveFig::Bool=false, w::Float64=8.50, h::Float64=8.50)


    FileName    = GetFileName(SC, nb)
    solver      = GetSolver(SC, nb)
    
    PyPlotFigure(w=w, h=h, bottom=1.5)
    PlotContour(solver.u[1], solver.fes)
    PlotMesh!(SC, nb, color="w")
    ylim(-0.25, 0.25)
    tick_params(axis="both", which="both", labelsize=TickSize)
    if SaveFig
        savefig("$(FigUbi)SC$(SC)_Mesh.png", dpi=800, pad_inches=0)
    end
    
    PyPlotFigure(w=w, h=h, bottom=1.5)
    PlotContour(solver.u[1], solver.fes)
    PlotMesh!(SC, nb, color="w")
    tick_params(axis="both", which="both", labelsize=TickSize)
    x0          = 0.07
    y0          = -0.07
    deltax      = 0.05
    xlim(x0 - deltax, x0 + deltax)
    ylim(y0 - deltax, y0 + deltax)
    if SaveFig
        savefig("$(FigUbi)SC$(SC)_MeshZoom.png", dpi=800, pad_inches=0)
    end
    
    return
    
end

function Video_RayleighTaylor(SC::Int; SaveFig::Bool=false, w::Float64=8.50, h::Float64=8.50)

    fig         = PyPlotSubPlots(1, 2, w=w, h=h, top=1.0)
    
    GasModel    = GasIdeal()
    kSave       = 0
    for nb=vcat(0:38,1000)
    
        FileName    = GetFileName(SC, nb)
        solver      = GetSolver(SC, nb)
    
        PyPlot.subplot(1,2,1)
        PyPlot.cla()
        PlotContour(solver.u[1], solver.fes)
#         PlotMesh!(SC, nb, color="w")
        tick_params(axis="both", which="both", labelsize=TickSize)
        title(latexstring(LatexString(GasModel, "rho"),"; t^n=", sprintf1("%.2e", solver.t)),
            fontsize=10)
        ylim(-0.25, 0.25)
        
        PyPlot.subplot(1,2,2)
        PyPlot.cla()
        PlotContour(solver.u[1], solver.fes)
        PlotMesh!(SC, nb, color="w")
        tick_params(axis="both", which="both", labelsize=TickSize)
        title(latexstring(LatexString(GasModel, "rho"),"; t^n=", sprintf1("%.2e", solver.t)),
            fontsize=10)
        ylim(-0.25, 0.25)
        
        if SaveFig
            savefig("$(FigUbi)SC$(SC)_$(kSave).png", dpi=800, pad_inches=0)
        end
        kSave       += 1
        
        pause(1.0)
        
    end
    
    return
    
end
