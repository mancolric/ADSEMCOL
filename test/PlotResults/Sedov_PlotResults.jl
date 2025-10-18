include("PlotResults.jl")
include("../test_CompressibleFlow/InputData/CompressibleFlow_SC.jl")

#Theoretical solution:
function utheor_Sedov(t::Float64, x::Vector{<:VecOrMat{Float64}})
    
    gamma           = 1.4
    alpha           = 0.9840418765432098
    rho0            = 1.0
    E0              = alpha*1.0^4    #alpha*Rf^4
    
    Rshock          = (E0*t^2/rho0/alpha)^(1/4)
    ushock          = 0.5*Rshock/(t+1e-14)
    rho1            = 1.0
    rho2            = (gamma+1.0)/(gamma-1.0)*rho1
    v2              = 2.0/(gamma+1.0)*ushock
    p2              = 2.0/(gamma+1.0)*rho1*ushock^2.0
    xi              = @mlv sqrt(x[1]^2.0+x[2]^2.0)/(Rshock+1e-14)
    rho             = @mlv rho2 * $Sedov2D_frho(xi)
    vr              = @mlv v2 * $Sedov2D_fv(xi)
    p               = @mlv p2 * $Sedov2D_fp(xi)
    rhovx           = @mlv rho * vr * x[1]/sqrt.(x[1]^2.0 + x[2]^2.0 .+ 1e-10)
    rhovy           = @mlv rho * vr * x[2]/sqrt.(x[1]^2.0 + x[2]^2.0 .+ 1e-10)
    rhoE            = @mlv p/(gamma-1.0) .+ 0.5*rho*vr^2.0 #Note that vr^2=vx^2+vy^2
    
    return [ rho, rhovx, rhovy, rhoE ]
end

function RadialPlot_Sedov(SC::Int, nb::Int; 
    SaveFig::Bool=false, w::Float64=9.50, h::Float64=6.50,
    PlotVars::Vector{String}=["rho", "v", "p"], mFig::Int=1, nFig::Int=3)

    #Problem data:
    gamma           = 1.4
    alpha           = 0.9840418765432098
    rho0            = 1.0
    E0              = alpha*1.0^4    #alpha*Rf^4
    
    FileName    = GetFileName(SC, nb)
    solver      = GetSolver(SC, nb)
    
    #Max R to plot:
    Rshock      = (E0*solver.t^2/rho0/alpha)^(1/4)
    Rplot       = min(1.2, 1.2*Rshock)
            
    #Theoretical solution:
    delta   = 1e-3
    x1v     = vcat(linspace(1e-10, Rshock-delta, 100), 
                linspace(Rshock-delta, Rshock+delta, 1000), 
                linspace(Rshock+delta, Rplot, 100))
    x2v     = zeros(size(x1v))
    xm      = [ x1v, x2v ]
    um      = utheor_Sedov(solver.t, xm)
    
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
    for ii=1:length(PlotVars)
        PyPlotFigure(w=w, h=h, bottom=1.0, top=0.7, left=1.0)
        
        #Theoretical solution:
        plot(x1v, DepVars(GasModel, solver.t, Vector{<:AMF64}(xm), Vector{<:AMF64}(um), [PlotVars[ii]])[1][1][:], "-xr", linewidth=1.0)
        
        #Numerical solution:
        splot_fun(x1,x2)    = @mlv sqrt(x1^2+x2^2)
        PlotNodes(splot_fun, solver, GasModel, PlotVars[ii])
        
        tick_params(axis="both", which="both", labelsize=TickSize)
        xlim([0.0,Rplot])
        xlabel(latexstring("r"), fontsize=10)
        title(latexstring(LatexString(PlotVars[ii]),
            "; t^n=", sprintf1("%.2e", solver.t)), 
            fontsize=10)
        legend(["Exact", "Num."], fontsize=TickSize)
           
        #Save figure:
        if SaveFig
            savefig("$(FigUbi)Sedov_RadialPlot_$(PlotVars[ii]).png", dpi=800, pad_inches=0)
        end
        
    end
    
    #Idem with zoom:
    for ii=1:length(PlotVars)
        PyPlotFigure(w=w, h=h, bottom=1.0, top=0.7, left=1.0)
        
        #Theoretical solution:
        plot(x1v, DepVars(GasModel, solver.t, Vector{<:AMF64}(xm), Vector{<:AMF64}(um), [PlotVars[ii]])[1][1][:], "-r", linewidth=1.0)
        
        #Numerical solution:
        splot_fun(x1,x2)    = @mlv sqrt(x1^2+x2^2)
        PlotNodes(splot_fun, solver, GasModel, PlotVars[ii])
        
        tick_params(axis="both", which="both", labelsize=TickSize)
        xlim([Rshock-delta, Rshock+delta])
        xlabel(latexstring("r"), fontsize=10)
        xticks([])
        xticks([Rshock-delta, Rshock, Rshock+delta], 
                [sprintf1("%.3f", Rshock-delta),  
                sprintf1("%.3f", Rshock),  
                sprintf1("%.3f", Rshock+delta)]) 
        title(latexstring(LatexString(PlotVars[ii]),
            "; t^n=", sprintf1("%.2e", solver.t)), 
            fontsize=10)
        legend(["Exact", "Num."], fontsize=TickSize)
        
        #Save figure:
        if SaveFig
            savefig("$(FigUbi)Sedov_RadialPlotZoom_$(PlotVars[ii]).png", dpi=800, pad_inches=0)
        end
        
    end
    
    
    return
    
end

function Contour_Sedov(SC::Int, nb::Int; SaveFig::Bool=false, w::Float64=8.50, h::Float64=8.50, 
    PlotVars::Vector{String}=["rho", "v", "p"], mFig::Int=1, nFig::Int=3, 
    PlotTitle::Bool=false)

    #Problem data:
    gamma           = 1.4
    alpha           = 0.9840418765432098
    rho0            = 1.0
    E0              = alpha*1.0^4    #alpha*Rf^4
    
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
        
        if PlotTitle
            PyPlotFigure(w=w, h=h+0.8, top=0.8, bottom=0.0, left=0.0, right=0.0)
        else
            PyPlotFigure(w=w, h=h, top=0.0, bottom=0.0, left=0.0, right=0.0)
        end
    
        #Numerical solution:
        PlotContour(solver, GasModel, PlotVars[ii], delta=1e-5)
        tick_params(axis="both", which="both", labelsize=TickSize)
        axis("off")
        if PlotTitle
            title(latexstring(LatexString(PlotVars[ii]),
                "; t^n=", sprintf1("%.2e", solver.t)), 
                fontsize=10)
        end
        if SaveFig
            savefig("$(FigUbi)SC$(SC)_$(nb)_Contour_$(PlotVars[ii]).png", dpi=800, pad_inches=0)
        end
    
    end
    
    return
    
end

function Mesh_Sedov(SC::Int, nb::Int; 
    PlotVar::String="p", 
    SaveFig::Bool=false, w::Float64=8.50, h::Float64=8.50,
    PlotTitle::Bool=false)

    #Define gas:
    GasModel                = GasIdeal()
    GasModel.gamma          = 1.4
    GasModel.epsilon        = 0.0
    GasModel.nu             = 0.0
    GasModel.beta           = 0.0
    GasModel.kappa_rho_cv   = 0.0

    FileName    = GetFileName(SC, nb)
    solver      = GetSolver(SC, nb)
    
#     PyPlotFigure(w=w, h=h, bottom=1.5)
    if PlotTitle
        PyPlotFigure(w=w, h=h+0.8, top=0.8, bottom=0.0, left=0.0, right=0.0)
    else
        PyPlotFigure(w=w, h=h, top=0.0, bottom=0.0, left=0.0, right=0.0)
    end
    PlotContour(solver, GasModel, PlotVar, delta=1e-5)
    PlotMesh!(SC, nb, color="w")
    tick_params(axis="both", which="both", labelsize=TickSize)
    axis("off")
    if PlotTitle
        title(latexstring(LatexString(PlotVar),
            "; t^n=", sprintf1("%.2e", solver.t)), 
            fontsize=10)
    end
    if SaveFig
        savefig("$(FigUbi)SC$(SC)_$(nb)_Mesh.png", dpi=800, pad_inches=0)
    end
    
#     PyPlotFigure(w=w, h=h, bottom=1.5)
    PyPlotFigure(w=w, h=h, top=0.0, bottom=0.0, left=0.0, right=0.0)
#     PlotContour(solver.u[1], solver.fes)
    PlotContour(solver, GasModel, PlotVar, delta=1e-5)
    PlotMesh!(SC, nb, color="w")
    tick_params(axis="both", which="both", labelsize=TickSize)
    deltax      = 5e-2
    xlim(1.0/sqrt(2.0) - deltax, 1.0/sqrt(2.0) + deltax)
    ylim(1.0/sqrt(2.0) - deltax, 1.0/sqrt(2.0) + deltax)
    axis("off")
    if SaveFig
        savefig("$(FigUbi)SC$(SC)_$(nb)_MeshZoom.png", dpi=800, pad_inches=0)
    end
    display([1.0/sqrt(2.0) - deltax, 1.0/sqrt(2.0) + deltax])
    
#     PyPlotFigure(w=w, h=h, bottom=1.5)
    PyPlotFigure(w=w, h=h, top=0.0, bottom=0.0, left=0.0, right=0.0)
#     PlotContour(solver.u[1], solver.fes)
    PlotContour(solver, GasModel, PlotVar, delta=1e-3)
    PlotMesh!(SC, nb, color="w")
    tick_params(axis="both", which="both", labelsize=TickSize)
    deltax      = 5e-4
    xlim(1.0/sqrt(2.0) - deltax, 1.0/sqrt(2.0) + deltax)
    ylim(1.0/sqrt(2.0) - deltax, 1.0/sqrt(2.0) + deltax)
    axis("off")
    if SaveFig
        savefig("$(FigUbi)SC$(SC)_$(nb)_MeshZoom2.png", dpi=800, pad_inches=0)
    end
    display([1.0/sqrt(2.0) - deltax, 1.0/sqrt(2.0) + deltax])
    
    return
    
end

function Video_Sedov(SC::Int; SaveFig::Bool=false, w::Float64=8.50, h::Float64=w)

    #Define gas:
    GasModel                = GasIdeal()
    GasModel.gamma          = 1.4
    GasModel.epsilon        = 0.0
    GasModel.nu             = 0.0
    GasModel.beta           = 0.0
    GasModel.kappa_rho_cv   = 0.0
    
    fig     = PyPlotSubPlots(1, 2, w=w, h=h, top=1.0)
    
    kSave   = 0
    for nb=vcat(0:49,1000)
    
        FileName    = GetFileName(SC, nb)
        solver      = GetSolver(SC, nb)
    
        #Position of the shock:
        gamma   = 1.4
        alpha   = 0.9840418765432098
        rho0    = 1.0
        E0      = alpha*1.0^4    #alpha*Rf^4
        Rshock  = (E0*solver.t^2/rho0/alpha)^(1/4)
        Rplot   = min(1.2, 1.2*Rshock)
        
        #Theoretical solution:
        delta   = 1e-3
        x1v     = vcat(linspace(1e-10, Rshock-delta, 100), 
                    linspace(Rshock-delta, Rshock+delta, 1000), 
                    linspace(Rshock+delta, Rplot, 100))
        x2v     = zeros(size(x1v))
        xm      = [ x1v, x2v ]
        um      = utheor_Sedov(solver.t, xm)
        
        PyPlot.subplot(1,2,1)
        PyPlot.cla()
#         PlotContour(solver.u[1], solver.fes)
        sv      = @. sqrt(solver.fes.PSpace.NodesCoords[:,1]^2 + solver.fes.PSpace.NodesCoords[:,2]^2)
        uP      = GetPolynomials(solver.u, solver.fes)
        plot(sv, uP[1], ".b", markersize=0.5)
        plot(xm[1], um[1], "k", linewidth=1.0)
#         xlim(0, Rplot)
        xlim(Rshock-1e-3, Rshock+1e-3)
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
