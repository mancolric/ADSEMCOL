include("PlotResults.jl")

function Contour_DetonationWave(SC::Int, nb::Int; SaveFig::Bool=false, w::Float64=30.0, h::Float64=8.50,
                                  PlotVars::Vector{String}=["Q_reac", "T", "p"], mFig::Int=1, nFig::Int=3)
    N = length(PlotVars)
    #Define gas:
    GasModel    = GasH2()

    FileName    = GetFileName(SC, nb)
    solver      = GetSolver(SC, nb)

    fig         = PyPlotSubPlots(N, 1, w=w, h=h, top=1.0)

    for ii=1:length(PlotVars)

        PyPlot.subplot(N,1,ii)
        PyPlot.cla()

        #Numerical solution:
        splot_fun(x1,x2)    = @mlv x1
        PlotNodes(splot_fun, solver, GasModel, PlotVars[ii])
        title(latexstring(LatexString(PlotVars[ii])),
              fontsize=10)
        tick_params(axis="both", which="both", labelsize=TickSize)

        if SaveFig
            savefig("$(FigUbi)SC$(SC)_Contour_$(PlotVars[ii]).png", dpi=800, pad_inches=0)
        end

    end

    return

end

function Line_DetonationWave(SC::Int, nb::Int; SaveFig::Bool=false, w::Float64=15.0, h::Float64=4.50,
                                PlotVar::String="Q_reac")

    #Define gas:
    GasModel    = GasH2()

    FileName    = GetFileName(SC, nb)
    solver      = GetSolver(SC, nb)

    fig         = PyPlotFigure(w=w, h=h, top=0.5, bottom=0.1)

    #Numerical solution:
    splot_fun(x1,x2)    = @mlv x1
    PlotNodes(splot_fun, solver, GasModel, PlotVar)
    xlim(0,0.4)
#     title(latexstring(LatexString(PlotVar)),
#             fontsize=10)
    text(0.35, -0.4e9, latexstring("t = $(nb) \\mu s"))
    #text(0.34, -0.4e9, latexstring("\\epsilon = 1·10^{-1}"))
    tick_params(axis="both", which="both", labelsize=TickSize)

    if SaveFig
        savefig("$(FigUbi)SC$(SC)_$(nb)_Line_$(PlotVar).png", dpi=800, pad_inches=0)
    end

    return

end

function Mesh_DetonationWave(SC::Int, nb::Int; PlotVar::String = "T", SaveFig::Bool=false, w::Float64=15.0, h::Float64=4.50)


    FileName    = GetFileName(SC, nb)
    solver      = GetSolver(SC, nb)

    PyPlotFigure(w=w, h=h, bottom=1.5)
#     PlotContour(solver.u[1], solver.fes)
#     PlotContour(solver, GasH2(), PlotVar, delta=1e-5)
    PlotMesh!(SC, nb, color="b")
    text(0.4, 0.0037, latexstring("t = $(nb) \\mu s"))
    text(0.4, 0.0029, latexstring("\\epsilon = 1·10^{-2}"))
    tick_params(axis="both", which="both", labelsize=TickSize)
    if SaveFig
        savefig("$(FigUbi)SC$(SC)_$(nb)_Mesh.png", dpi=800, pad_inches=0)
    end

    return

end

function Video_DetonationWave(SC::Int, var::String; SaveFig::Bool=false, w::Float64=30.0, h::Float64=8.50)

    fig      = PyPlotSubPlots(2, 1, w=w, h=h, top=1.0)
    GasModel = GasH2()

    kSave   = 0
    for nb=0:1000

        FileName    = GetFileName(SC, nb)
        solver      = GetSolver(SC, nb)

        PyPlot.subplot(2,1,1)
        PyPlot.cla()
        splot_fun(x1,x2)    = @mlv x1
        PlotNodes(splot_fun, solver, GasModel, var)
        #         PlotMesh!(SC, nb, color="w")
        tick_params(axis="both", which="both", labelsize=TickSize)
        title(latexstring(LatexString(var),"; t^n=", sprintf1("%.2e", solver.t)),
              fontsize=10)

        PyPlot.subplot(2,1,2)
        PyPlot.cla()
#         PlotContour(solver.u[1], solver.fes)
        PlotMesh!(SC, nb, color="b")
        tick_params(axis="both", which="both", labelsize=TickSize)
        title(latexstring(LatexString(var),"; t^n=", sprintf1("%.2e", solver.t)),
              fontsize=10)

        if SaveFig
            savefig("$(VideosUbi)SC$(SC)_$(kSave).png", dpi=800, pad_inches=0)
        end
        kSave       += 1

        pause(1.0)

    end

    return

end

function Line_DetonationWave2(SC::Int, var::String; SaveFig::Bool=false, w::Float64=30.0, h::Float64=8.50)

    GasModel    = GasH2()
    SC_final    = load(GetFileName(SC, "info"), "nb_SaveRes")

#     iv          = [157, 188, 219, 235]
#     iv          = [0, 31, 62, 94, 125]
    iv          = Int.(ceil.(linspace(0.0, SC_final*1.0, 8)))
    
    N           = length(iv)
    fig         = PyPlotSubPlots(N, 1, w=w, h=h, top=1.0)
    cont        = 1
    for i = iv

        FileName    = GetFileName(SC, i)
        solver      = GetSolver(SC, i)

        PyPlot.subplot(N,1,cont)
#         PyPlot.cla()
        splot_fun(x1,x2)    = @mlv x1
        PlotNodes(splot_fun, solver, GasModel, var)
        tick_params(axis="both", which="both", labelsize=TickSize)
        title(latexstring("t^n=", sprintf1("%.0f", solver.t*1e6), " \\mu s"),
              fontsize=10)
        cont += 1
    end

    if SaveFig
        savefig("$(FigUbi)SC$(SC)_$(t0).png", dpi=800, pad_inches=0)
    end

end

function Err_DetonationWave(SC::Int, nb::Int; SaveFig::Bool=false, w::Float64=30.0, h::Float64=8.50,
                                PlotVars::Vector{String}=["rho", "T", "p"], mFig::Int=1, nFig::Int=3)
    N = 3
    #Define gas:
    GasModel    = GasH2()

    FileName    = GetFileName(SC, nb)
    solver      = GetSolver(SC, nb)

    figure()
    PyPlot.cla()
    semilogy(solver.tv, solver.etaSv, ".-b")
    semilogy(solver.tv, solver.etaTv, ".-g")
    semilogy(solver.tv, solver.etaAv, ".-r")
    if true
        validv  = solver.validv .== 1
        semilogy(solver.tv[validv], solver.etaSv[validv], "sb")
        semilogy(solver.tv[validv], solver.etaTv[validv], "sg")
        semilogy(solver.tv[validv], solver.etaAv[validv], "sr")
    end
    legend(["space", "time", "algebraic"])
    xlabel(L"t")

    if SaveFig
        savefig("$(FigUbi)SC$(SC)_err.png", dpi=800, pad_inches=0)
    end

    return

end

function escalon(delta::Vector{Float64}; SaveFig::Bool=false, w::Float64=10.0, h::Float64=5.0)

    x = [i for i = -5.0:0.01:5.0]

    PyPlotFigure(w=w, h=h, bottom=0.6)

    for j = 1:length(delta)
        plot(x, SmoothHeaviside.(x, delta[j], 2.0, 5.0))
    end

    legend([latexstring("\\delta = $(delta[1])"),latexstring("\\delta = $(delta[2])"),latexstring("\\delta = $(delta[3])")])

    if SaveFig
        savefig("$(FigUbi)deltas.png", dpi=800, pad_inches=0)
    end

end
