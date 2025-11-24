include("../test_NonHydrostaticWaterWaves/SolitonRelaxed.jl")
include("PlotResults.jl")

function Nodes_SolitonRelaxed(SC::Int, nb::Int; 
    SaveFig::Bool=false, w::Float64=9.50, h::Float64=6.50,
    PlotVars::Vector{String}=["h", "v1", "v3", "p"], mFig::Int=2, nFig::Int=2)

    FileName    = GetFileName(SC, nb)
    solver      = GetSolver(SC, nb)
    SolitonId   = load(FileName, "SolitonId")
    
    #Load data:
    SolitonData         = readdlm("$(@__DIR__)/../test_NonHydrostaticWaterWaves/SolitonRelaxed$(SolitonId)_data.txt")
    SolitonCoeffs       = readdlm("$(@__DIR__)/../test_NonHydrostaticWaterWaves/SolitonRelaxed$(SolitonId)_coefficients.txt")
    alpha               = SolitonData[1]
    g                   = SolitonData[2]
    h0                  = SolitonData[3]
    A                   = SolitonData[4]
    
    #Define model:
    model               = NHWW()
    model.g             = g
    model.c             = alpha*sqrt(g*h0)

    #Domain limits:
    c0      = sqrt(g*(A+h0))
    l0      = h0*sqrt((A+h0)/h0)
    x11     = c0*solver.t-10*l0
    x12     = c0*solver.t+10*l0
            
    #Exact solution:
    x1v     = linspace(x11, x12, 500)
    x2v     = zeros(size(x1v))
    xm      = [ x1v, x2v ]
    um      = SolitonRelaxedFun(solver.t, Vector{<:AMF64}(xm), SolitonData, SolitonCoeffs)
    vm      = DepVars(model, solver.t, Vector{<:AMF64}(xm), um, PlotVars)
    
    #Loop plot variables:
#     fig     = PyPlotSubPlots(mFig, nFig, w=w, h=h, left=0.9, right=0.4, bottom=1.1, top=1.0)
    fig     = PyPlotSubPlots(mFig, nFig, w=w, h=h, left=0.9, right=0.4, bottom=1.1, top=1.2) #for suptitle
    for ii=1:length(PlotVars)
        PyPlot.subplot(mFig, nFig, ii)
        PyPlot.cla()
        splot_fun(x1,x2)    = @mlv x1
        PlotNodes(splot_fun, solver, model, PlotVars[ii])
        plot(x1v, vm[ii][1], "r", linewidth=0.5)
        xlabel(latexstring("x_1"), fontsize=10)
        title(latexstring(LatexString(model, PlotVars[ii]),
                            "; t^n=", sprintf1("%.2e", solver.t)),
                fontsize=10)
        tick_params(axis="both", which="both", labelsize=TickSize)
        xlim([x11, x12])
        suptitle(latexstring("\\mathrm{Relaxed~soliton}, \\alpha=", sprintf1("%d",alpha), ", Tol_S=Tol_T=", sprintf1("%.2e", solver.TolT)), fontsize=10)
    end
    if SaveFig
        savefig("$(FigUbi)SC$(SC)_Nodes.png", dpi=800, pad_inches=0)
    end
    
    return
    
end

function CompareSpaceAdaptation_SolitonRelaxed(StudyCase::String; nb::Int=1000, 
    SaveFig::Bool=false, w::Float64=8.50, h::Float64=8.50)

    #------------------------------------------------------------
    
    SCvv0   = NaN
    SCvv1   = NaN
#     nb      = NaN
    
    #Choose normal5_L2:
    if StudyCase=="normal"

                        
        SCvv1       = [ 110000:110004, 
                        110005:110009,
                        110010:110014
                        ]
        
    elseif StudyCase=="IRK"
                        
        SCvv1       = [ 110000:110004, 
                        110005:110009,
                        110040:110044,
                        110045:110049 
                        ]
        
    end
        
    #------------------------------------------------------------
    
    hpvv1, errvv1, etavv1, 
        tCPUvv1, alphavv1   = GetVbles(SCvv1, ["hp", "e_Lq", "etaST", "tCPU", "Soliton_alpha"], nb=nb)
    EOCvv1                  = ExpOrderConv(hpvv1, errvv1)
    
    PyPlotFigure(w=w, h=h, bottom=1.5, left=1.7)
    leg                     = Legend(alphavv1, latexstring("\\alpha"), format="%d")
#     loglog(NaN,NaN, "-k", linewidth=0.5)
#     loglog(NaN,NaN, "--k", linewidth=0.5)
#     leg                     = vcat(leg, L"e_{ST}", latexstring("\\mathcal{E}_{ST}"))
    PlotXY(hpvv1, errvv1, linestyle="solid", marker="s")
#     PlotXY(hpvv1, etavv1, linestyle="dashed", marker="s")
    xlabel(L"h/p", fontsize=LabelSize)
    ylabel(latexstring(GetString("errL2L2")), rotation=0, fontsize=LabelSize, labelpad=10.0)
    tick_params(axis="both", which="both", labelsize=TickSize)
    legend(leg, fontsize=TickSize)
    grid("on")
    if SaveFig
        savefig("$(FigUbi)SolitonRelaxed_SpaceAdapt1.png", dpi=800, pad_inches=0)
    end
    
    PyPlotFigure(w=w, h=h, bottom=1.5, left=1.5)
    leg                     = Legend(alphavv1, latexstring("\\alpha"), format="%d")
    PlotXY(tCPUvv1, errvv1, linestyle="solid", marker="s")
    grid("on")
    xlabel(latexstring(GetString("tCPU")), fontsize=LabelSize)
    ylabel(latexstring(GetString("errL2L2")), rotation=0, fontsize=LabelSize, labelpad=10.0)
    legend(leg, fontsize=TickSize)
    tick_params(axis="both", which="both", labelsize=TickSize)
    if SaveFig
        savefig("$(FigUbi)SolitonRelaxed_SpaceAdapt2.png", dpi=800, pad_inches=0)
    end
    
    println("EOC")
    display(EOCvv1)
    println("i_eff")
    for ii=1:length(etavv1)
        println(etavv1[ii]./errvv1[ii])
    end
    
    return
    
end

function CompareTimeAdaptation_SolitonRelaxed(StudyCase::String; 
    SaveFig::Bool=false, w::Float64=8.50, h::Float64=8.50)

    SCvv0   = NaN
    SCvv1   = NaN
    nb      = NaN
    
    if StudyCase=="normal"
    
        #TimeAdapt:
        SCvv1       = [ 110015:110019, 
                        110020:110024, 
                        110025:110026 ]
                        
        nb          = 1000
        
    elseif StudyCase=="IRK"
    
        #TimeAdapt:
        SCvv1       = [ 110015:110019, 
                        110020:110024, 
                        110030:110034, 
                        110035:110039 ]
                        
        nb          = 1000
        
    end
    
    #------------------------------------------------------------
    
    Deltatvv1, errvv1, etavv1, 
        tCPUvv1, CFLvv1, 
        alphavv1                = GetVbles(SCvv1, ["Deltat", "e_Lq", "etaST", 
                                                    "tCPU", "CFLmax", "Soliton_alpha"], nb=nb)
    EOCvv1                      = ExpOrderConv(Deltatvv1, errvv1)
    
    PyPlotFigure(w=w, h=h, bottom=1.5, left=1.7)
    leg                     = Legend(alphavv1, latexstring("\\alpha"), format="%d")
    PlotXY(Deltatvv1, errvv1, linestyle="solid", marker="s")
#     PlotXY(Deltatvv1, etavv1, linestyle="dashed", marker="s")
    xlabel(latexstring("\\tau"), fontsize=LabelSize)
    ylabel(latexstring(GetString("errL2L2")), rotation=0, fontsize=LabelSize, labelpad=10.0)
    tick_params(axis="both", which="both", labelsize=TickSize)
    legend(leg, fontsize=TickSize)
    grid("on")
    if SaveFig
        savefig("$(FigUbi)SolitonRelaxed_TimeAdapt1.png", dpi=800, pad_inches=0)
    end
            
    PyPlotFigure(w=w, h=h, bottom=1.5, left=1.5)
    leg                     = Legend(alphavv1, latexstring("\\alpha"), format="%d")
    PlotXY(tCPUvv1, errvv1, linestyle="solid", marker="s")
    grid("on")
    xlabel(latexstring(GetString("tCPU")), fontsize=LabelSize)
    ylabel(latexstring(GetString("errL2L2")), rotation=0, fontsize=LabelSize, labelpad=10.0)
    legend(leg, fontsize=TickSize)
    tick_params(axis="both", which="both", labelsize=TickSize)
    if SaveFig
        savefig("$(FigUbi)SolitonRelaxed_TimeAdapt2.png", dpi=800, pad_inches=0)
    end
    
    PyPlotFigure(w=w, h=h, bottom=1.5, left=1.5)
    leg                     = Legend(alphavv1, latexstring("\\alpha"), format="%d")
    PlotXY(CFLvv1, errvv1, linestyle="solid", marker="s")
    grid("on")
    xlabel(latexstring(GetString("CFLmax")), fontsize=LabelSize)
    ylabel(latexstring(GetString("errL2L2")), rotation=0, fontsize=LabelSize, labelpad=10.0)
    legend(leg, fontsize=TickSize)
    tick_params(axis="both", which="both", labelsize=TickSize)
    if SaveFig
        savefig("$(FigUbi)SolitonRelaxed_TimeAdapt3.png", dpi=800, pad_inches=0)
    end
    
    display(EOCvv1)
    
    return
    
end
