include("PlotResults.jl")
include("../test_CompressibleFlow/InputData/CompressibleFlow_SC.jl")

#Theoretical solution:
function utheor_ShuOsher(t::Float64, x::Vector{<:VecOrMat{Float64}})
    
    gamma           = 1.4
    x1dims          = size(x[1])
    
    udisc           = 3.549648
    xm              = hcat(x[1][:], x[2][:])
    rho             = ShuOsher_FV(xm, "rho", 3200, "5Diaz")
    vx              = ShuOsher_FV(xm, "u", 3200, "5Diaz")
    rhovx           = @. rho*vx
    rhovy           = @. 0.0*rhovx
    p               = ShuOsher_FV(xm, "p", 3200, "5Diaz")
    rhoE            = @. p/(gamma-1.0) + 0.5*rhovx^2/rho
        
    rho             = reshape(rho, x1dims)
    rhovx           = reshape(rhovx, x1dims)
    rhovy           = reshape(rhovy, x1dims)
    rhoE            = reshape(rhoE, x1dims)

    return [ rho, rhovx, rhovy, rhoE ]
end
        
function x1Plot_ShuOsher(SC::Int, nb::Int; 
    SaveFig::Bool=false, w::Float64=9.50, h::Float64=6.50,
    PlotVars::Vector{String}=["rho"], mFig::Int=1, nFig::Int=3)

    #Problem data:
    gamma           = 1.4
    xdisc           = -4.0
    ML              = 3.0
    rhoL            = 3.857143
    uL              = 2.629369
    pL              = 10.33333
    rhoEL           = pL/(gamma-1.0)+0.5*rhoL*uL^2
    pR              = 1.0
    
    FileName    = GetFileName(SC, nb)
    solver      = GetSolver(SC, nb)

    #Define gas:
    GasModel                = GasIdeal()
    GasModel.gamma          = 1.4
    GasModel.epsilon        = load(FileName, "epsilon")
    GasModel.nu             = load(FileName, "nu")
    GasModel.beta           = load(FileName, "beta")
    GasModel.kappa_rho_cv   = load(FileName, "kappa_rho_cv")
    
    #Theoretical solution:
    Nx      = 200
    x1v     = linspace(-5.0, 5.0, Nx)
    x2v     = zeros(Nx)
    xm      = [ x1v, x2v ]
    um      = utheor_ShuOsher(solver.t, xm)
    
    #No zoom - loop plot variables:
    for ii=1:length(PlotVars)
        
        PyPlotFigure(w=w+0.7, h=h, bottom=1.0, top=0.7, left=0.6)
        plot(x1v, DepVars(GasModel, solver.t, Vector{<:AMF64}(xm), Vector{<:AMF64}(um), [PlotVars[ii]])[1][1][:], "-xr", linewidth=0.5)
        
        splot_fun(x1,x2)    = @mlv x1
        PlotNodes(splot_fun, solver, GasModel, PlotVars[ii])
        
        legend(["Ref.", "Num."], fontsize=TickSize)
        tick_params(axis="both", which="both", labelsize=TickSize)
        xlabel(latexstring("x_1"), fontsize=10)
        title(latexstring(LatexString(GasModel, PlotVars[ii]),
            "; t^n=", sprintf1("%.2e", solver.t)), 
            fontsize=10)
        if SaveFig
            savefig("$(FigUbi)SC$(SC)_x1Plot_$(PlotVars[ii]).png", dpi=800, pad_inches=0)
        end
       
    end
    
    #Theoretical solution:
    Nx      = 200
    x1v     = linspace(0.5, 2.5, Nx)
    x2v     = zeros(Nx)
    xm      = [ x1v, x2v ]
    um      = utheor_ShuOsher(solver.t, xm)
    
    #With zoom - loop plot variables:
    for ii=1:length(PlotVars)

        PyPlotFigure(w=w, h=h, bottom=1.0, top=0.7, left=0.2)
        plot(x1v, DepVars(GasModel, solver.t, Vector{<:AMF64}(xm), Vector{<:AMF64}(um), [PlotVars[ii]])[1][1][:], "-xr", linewidth=0.5)
        
        splot_fun(x1,x2)    = @mlv x1
        PlotNodes(splot_fun, solver, GasModel, PlotVars[ii])
        
        legend(["Ref.", "Num."], fontsize=TickSize)
        tick_params(axis="both", which="both", labelsize=TickSize)
        xlim(0.5, 2.5)
        xlabel(latexstring("x_1"), fontsize=10)
        title(latexstring(LatexString(GasModel, PlotVars[ii]),
            "; t^n=", sprintf1("%.2e", solver.t)), 
            fontsize=10)
        if SaveFig
            savefig("$(FigUbi)SC$(SC)_x1Plot_$(PlotVars[ii])_Zoom.png", dpi=800, pad_inches=0)
        end
        
    end
    
    #Theoretical solution:
    Nx      = 100
    x1v     = linspace(2.38, 2.41, Nx)
    x2v     = zeros(Nx)
    xm      = [ x1v, x2v ]
    um      = utheor_ShuOsher(solver.t, xm)
    
    #With larger zoom - loop plot variables:
    for ii=1:length(PlotVars)

        PyPlotFigure(w=w, h=h, bottom=1.0, top=0.7, left=0.3)
        plot(x1v, DepVars(GasModel, solver.t, Vector{<:AMF64}(xm), Vector{<:AMF64}(um), [PlotVars[ii]])[1][1][:], "-r", linewidth=1.0)
        
        splot_fun(x1,x2)    = @mlv x1
        PlotNodes(splot_fun, solver, GasModel, PlotVars[ii])
        
        legend(["Ref.", "Num."], fontsize=TickSize)
        tick_params(axis="both", which="both", labelsize=TickSize)
        xlim(2.38, 2.41)
        xlabel(latexstring("x_1"), fontsize=10)
        yticks([1, 2, 3, 4], ["", "", "", ""])
        title(latexstring(LatexString(GasModel, PlotVars[ii]),
            "; t^n=", sprintf1("%.2e", solver.t)), 
            fontsize=10)
        if SaveFig
            savefig("$(FigUbi)SC$(SC)_x1Plot_$(PlotVars[ii])_Zoom2.png", dpi=800, pad_inches=0)
        end
        
    end
    
    return
    
end

function Mesh_ShuOsher(SC::Int, nb::Int; SaveFig::Bool=false, w::Float64=8.50, h::Float64=8.50)


    FileName    = GetFileName(SC, nb)
    solver      = GetSolver(SC, nb)
    
    PyPlotFigure(w=w, h=h, top=0.0, bottom=0.0, left=0.0, right=0.0)
    PlotContour(solver.u[1], solver.fes)
    PlotMesh!(SC, nb, color="w")
    tick_params(axis="both", which="both", labelsize=TickSize)
    axis("equal")
    axis("off")
    if SaveFig
        savefig("$(FigUbi)SC$(SC)_Mesh.png", dpi=800, pad_inches=0)
    end
    
    return
    
end
