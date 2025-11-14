include("PlotResults.jl")
include("../test_NonHydrostaticWaterWaves/SolitonExact.jl")

function Nodes_SolitonExact(SC::Int, nb::Int; 
    SaveFig::Bool=false, wFig::Float64=9.50, hFig::Float64=6.50,
    PlotVars::Vector{String}=["h", "v1", "v3", "p"], mFig::Int=2, nFig::Int=2)

    FileName    = GetFileName(SC, nb)
    solver      = GetSolver(SC, nb)
    
    #Load data:
    alpha               = 10.0
    g                   = 9.8
    h0                  = 1.0
    A                   = 0.2
    
    #Define model:
    model               = NHWW()
    model.g             = g
    model.c             = alpha*sqrt(g*h0)

    #Domain limits:
    c0                  = sqrt(g*(A+h0))
    l0                  = h0*sqrt((A+h0)/h0)
    x11                 = c0*solver.t-10*l0
    x12                 = c0*solver.t+10*l0
            
    #Exact solution:
    x1v                 = linspace(x11, x12, 500)
    x2v                 = zeros(size(x1v))
    xm                  = [ x1v, x2v ]
    eta, q1, q2, q3, p  = SolitonExactFun(solver.t, Matrix{Float64}(xm[1]), A=A, h0=h0, g=g)
    P                   = @. 0.5*eta*eta*(1.0+sqrt(1.0-12*p/(model.c*model.c)))
    b                   = @. 0.0*eta
    um                  = [eta, q1, q2, q3, P, b]
    vm                  = DepVars(model, solver.t, Vector{<:AMF64}(xm), um, PlotVars)
    
    #Loop plot variables:
    fig                 = PyPlotSubPlots(mFig, nFig, w=wFig, h=hFig, left=0.9, right=0.4, bottom=1.1, top=1.0)
    for ii=1:length(PlotVars)
        PyPlot.subplot(mFig, nFig, ii)
        PyPlot.cla()
        splot_fun(x1,x2)    = @mlv x1
        PlotNodes(splot_fun, solver, model, PlotVars[ii])
        plot(x1v, vm[ii][1], "r", linewidth=0.5)
        xlabel(latexstring("x_1"), fontsize=10)
        title(latexstring(LatexString(PlotVars[ii]),
                            "; t^n=", sprintf1("%.2e", solver.t)),
                fontsize=10)
        tick_params(axis="both", which="both", labelsize=TickSize)
        xlim([x11, x12])
    end
    if SaveFig
        savefig("$(FigUbi)SC$(SC)_Nodes.png", dpi=800, pad_inches=0)
    end
    
    return
    
end
