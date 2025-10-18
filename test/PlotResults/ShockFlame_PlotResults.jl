include("PlotResults.jl")

function GasFP_ShockFlame(SC::Int, nb::Int)

    #Extract data:
    FileName                = GetFileName(SC, nb)
    solver                  = GetSolver(SC, nb)
    gamma                   = 1.4
    q                       = load(FileName, "q")
    beta                    = load(FileName, "beta")
    ML                      = load(FileName, "ML")
    
    #Fuel state:
    p0                      = 1.0
    RT0                     = 1.0
    rho0                    = p0/RT0
    YF0                     = load(FileName, "YF0")
    Q0                      = load(FileName, "Q0")
    
    #Define gas:
    GasModel                = GasFP()
    GasModel.epsilon        = load(FileName, "epsilon")
    GasModel.gamma          = 1.4
    GasModel.nu             = load(FileName, "nu")
    GasModel.beta           = 0.0
    GasModel.kappa_rho_cv   = load(FileName, "kappa_rho_cv")
    GasModel.D              = load(FileName, "D")
    GasModel.B              = load(FileName, "B")
    GasModel.RTa            = load(FileName, "RTa")
    GasModel.hfF            = Q0
    GasModel.hfP            = 0.0
    
    return GasModel
    
end
    
function x1Plot_ShockFlame(SC::Int, nb::Int;
    PlotVars::Vector{String}=["Y_F", "v", "p", "RT", "rho", "M"],
    MarkSF::Bool=false,
    mFig::Int=3, nFig::Int=2, 
    SaveFig::Bool=false, w::Float64=9.50, h::Float64=6.50)
    
    #Extract data:
    if nb==1000
        nb                  = load("$(ResUbi)LIRKHyp_SC$(SC)_info.jld2", "nb_SaveRes")
    end
    GasModel    = GasFP_ShockFlame(SC,nb)
    FileName    = GetFileName(SC, nb)
    solver      = GetSolver(SC, nb)
    q           = load(FileName, "q")
    beta        = load(FileName, "beta")
    ML          = load(FileName, "ML")
    
    #Detect shock and flame positions:
    pref,YFref  = pYref_ShockFlame(SC)
    
    #No zoom - loop plot variables:
#     fig         = PyPlotSubPlots(mFig, nFig, w=w, h=h, top=1.0, right=0.1, left=2.0, bottom=0.7)
#     fig         = PyPlotSubPlots(mFig, nFig, w=w, h=h, top=1.0, right=0.1)
    fig         = PyPlotSubPlots(mFig, nFig, w=w, h=h, top=0.1, right=0.1, left=1.8, bottom=0.1, tTitle=0.8, bTitle=1.0)
    for ii=1:length(PlotVars)
        
        #Plot solution:
        subplot(mFig, nFig, ii)
        splot_fun(x1,x2)    = @mlv x1
        vbleString      = ""
        v_plot          = NaN
        if false && PlotVars[ii]=="vxRel" #unfinished - xsxdot not defined
            vxRel           = SolutionAtNodes_vxRel(solver, GasModel, xsdot)
            plot(solver.fes.PSpace.NodesCoords[:,1], vxRel, ".b", markersize=0.2)
            vbleString      = "v_x'"
            v_plot          = vxRel
        elseif false && PlotVars[ii]=="MRel"
            MRel            = SolutionAtNodes_MRel(solver, GasModel, xsdot)
            plot(solver.fes.PSpace.NodesCoords[:,1], MRel, ".b", markersize=0.2)
            vbleString      = "M'"
            v_plot          = MRel
        else    
            splot_fun(x1,x2)= @mlv x1       
            v_plot          = PlotNodes(splot_fun, solver, GasModel, PlotVars[ii])
            vbleString      = LatexString(PlotVars[ii])
        end
        tick_params(axis="both", which="both", labelsize=TickSize)
        xlabel(latexstring("x_1"), fontsize=10)
        ylabel(latexstring(LatexString(PlotVars[ii])), rotation=0, fontsize=10, labelpad=10)
            
        #Get positions of shock and flame:
        if MarkSF
            xs, xf      = xsxf_ShockFlame(solver, GasModel, pref, YFref)
            plot([xf, xf], [0.9*minimum(v_plot), 1.1*maximum(v_plot)], "--r")
            plot([xs, xs], [0.9*minimum(v_plot), 1.1*maximum(v_plot)], "--g")
        end
        
        #Plot M=1 line:
        if PlotVars[ii]=="M" || PlotVars[ii]=="MRel"
            xmin    = minimum(solver.fes.PSpace.NodesCoords[:,1])
            xmax    = maximum(solver.fes.PSpace.NodesCoords[:,1])
            plot([xmin, xmax], [1.0, 1.0], "--k")
        end
       
    end
    suptitle(latexstring("q=", q, ", \\beta=", beta, ", M=", ML, 
                        ", t=", sprintf1("%.2e", solver.t)), fontsize=10)
    
    if SaveFig
        savefig("$(FigUbi)SC$(SC)_$(nb)_x1Plot.png", dpi=800, pad_inches=0)
    end
        
    return
    
end

function Contour_ShockFlame(SC::Int, nb::Int; SaveFig::Bool=false, w::Float64=8.50, h::Float64=8.50, 
    PlotVars::Vector{String}=["Y_F", "v", "p", "RT", "rho", "M"], mFig::Int=3, nFig::Int=2, 
    PlotTitle::Bool=false)

    GasModel    = GasFP_ShockFlame(SC, nb)
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

function Mesh_ShockFlame(SC::Int, nb::Int; SaveFig::Bool=false, w::Float64=8.50, h::Float64=8.50)


    FileName    = GetFileName(SC, nb)
    solver      = GetSolver(SC, nb)
    
    PyPlotFigure(w=w, h=h, bottom=1.5)
    PlotContour(solver.u[1], solver.fes)
    PlotMesh!(SC, nb, color="w")
    tick_params(axis="both", which="both", labelsize=TickSize)
    if SaveFig
        savefig("$(FigUbi)SC$(SC)_Mesh.png", dpi=800, pad_inches=0)
    end
    
    return
    
end

#Video with velocity and Mach numbers relative to the shock.
function Video_ShockFlame(SC::Int; 
    PlotVars::Vector{String}=["Y_F", "v", "p", "RT", "rho", "M"],
    MarkSF::Bool=false, 
    mFig::Int=3, nFig::Int=2, 
    SaveFig::Bool=false, w::Float64=8.50, h::Float64=8.50,
    tPlot::Float64=30.0, n1::Int=0, n2=Int=1000)

    nbMax       = load("$(ResUbi)LIRKHyp_SC$(SC)_info.jld2", "nb_SaveRes")
    nbMax       = min(nbMax, n2)
    
    FileName    = GetFileName(SC, nbMax)
    q           = load(FileName, "q")
    beta        = load(FileName, "beta")
    ML          = load(FileName, "ML")
    GasModel    = GasFP_ShockFlame(SC,nbMax)
    
    tPlot_nFig  = tPlot/(nbMax-n1)
    
    #Track shock and flame positions:
    pref,YFref  = pYref_ShockFlame(SC)
    tv          = zeros(0)
    xsv         = zeros(0)
    xfv         = zeros(0)
    
#     fig         = PyPlotSubPlots(mFig, nFig, w=w, h=h, top=1.0)
    fig         = PyPlotSubPlots(mFig, nFig, w=w, h=h, top=0.1, right=0.1, left=1.8, bottom=0.1, tTitle=0.8, bTitle=1.0)
    kSave       = 0
    for nb=n1:nbMax
    
        t0          = time()
        
        FileName    = GetFileName(SC, nb)
        solver      = GetSolver(SC, nb)
        
        #Get positions of shock and flame
        xs, xf      = xsxf_ShockFlame(solver, GasModel, pref, YFref)
        push!(tv, solver.t)
        push!(xsv, xs)
        push!(xfv, xf)
        
        #Shock speed:
        xsdot       = 0.0
        if kSave>0
            xsdot   = (xsv[kSave+1]-xsv[kSave])/(tv[kSave+1]-tv[kSave])
        end
        
        for ii=1:length(PlotVars)
            subplot(mFig, nFig, ii)
            cla()
            vbleString      = ""
            v_plot          = NaN
            if PlotVars[ii]=="vxRel"
                vxRel           = SolutionAtNodes_vxRel(solver, GasModel, xsdot)
                plot(solver.fes.PSpace.NodesCoords[:,1], vxRel, ".b", markersize=0.2)
                vbleString      = "v_x'"
                v_plot          = vxRel
            elseif PlotVars[ii]=="MRel"
                MRel            = SolutionAtNodes_MRel(solver, GasModel, xsdot)
                plot(solver.fes.PSpace.NodesCoords[:,1], MRel, ".b", markersize=0.2)
                vbleString      = "M'"
                v_plot          = MRel
            else    
                splot_fun(x1,x2)= @mlv x1       
                v_plot          = PlotNodes(splot_fun, solver, GasModel, PlotVars[ii])
                vbleString      = LatexString(PlotVars[ii])
            end
            if MarkSF
                plot([xf, xf], [0.9*minimum(v_plot), 1.1*maximum(v_plot)], "--r")
                plot([xs, xs], [0.9*minimum(v_plot), 1.1*maximum(v_plot)], "--g")
            end
            if PlotVars[ii]=="M" || PlotVars[ii]=="MRel"
                xmin    = minimum(solver.fes.PSpace.NodesCoords[:,1])
                xmax    = maximum(solver.fes.PSpace.NodesCoords[:,1])
                plot([xmin, xmax], [1.0, 1.0], "--k")
            end
            
            tick_params(axis="both", which="both", labelsize=TickSize)
            xlabel(latexstring("x_1"), fontsize=10)
#             title(latexstring(vbleString,
#                 "; t^n=", sprintf1("%.2e", solver.t)), 
#                 fontsize=10)
            ylabel(latexstring(vbleString), rotation=0, fontsize=10, labelpad=10)

        end
        suptitle(latexstring("q=", q, ", \\beta=", beta, ", M=", ML, 
                    ", t=", sprintf1("%.2e", solver.t)), fontsize=10)
        
        if SaveFig
            savefig("$(FigUbi)SC$(SC)_$(kSave).png", dpi=800, pad_inches=0)
#             error("")
        end
        kSave       += 1
        
        t1          = time()
        
        pause(max(0.01, tPlot_nFig-(t1-t0)))
        
    end
    
    return
    
end

function Track_ShockFlame(SC::Int; SaveFig::Bool=false, w::Float64=10.0, h::Float64=10.0)

    nbMax       = load("$(ResUbi)LIRKHyp_SC$(SC)_info.jld2", "nb_SaveRes")
    
    FileName    = GetFileName(SC, nbMax)
    solver      = GetSolver(SC, nbMax)
    q           = load(FileName, "q")
    beta        = load(FileName, "beta")
    ML          = load(FileName, "ML")
    GasModel    = GasFP_ShockFlame(SC,nbMax)
    
    #Track shock and flame positions, as well as incoming (transmitted) relative Mach number:
    pref,YFref  = pYref_ShockFlame(SC)
    tv          = zeros(0)
    xsv         = zeros(0)
    xfv         = zeros(0)
    Mtv         = zeros(0)
    for nb=0:nbMax
    
        FileName    = GetFileName(SC, nb)
        solver      = GetSolver(SC, nb)
    
        #Get positions of shock and flame:
        xs, xf      = xsxf_ShockFlame(solver, GasModel, pref, YFref)
        push!(tv, solver.t)
        push!(xsv, xs)
        push!(xfv, xf)
        
        #Shock velocity:
        xsdot       = NaN
        if nb==0
            xsdot   = 0.0
        elseif nb==1
            xsdot   = (xsv[nb+1]-xsv[nb])/(tv[nb+1]-tv[nb])
        else
            xsdot   = (3*xsv[nb+1]-4*xsv[nb]+xsv[nb-1])/(3*tv[nb+1]-4*tv[nb]+tv[nb-1])
        end
        
        #Relative Mach number:
        if false && xf<xs #not working if flame does not pass shock (e.g. q=1, M=10)
            Mtv         = push!(Mtv, NaN)
        else
            Mrelv       = SolutionAtNodes_MRel(solver, GasModel, xsdot)
            x1v         = solver.mesh.NodesCoords[:,1]
            sortv       = sortperm(x1v)
            x1v         = x1v[sortv]
            Mrelv       = Mrelv[sortv]
            aux         = findfirst(x1v.>0.0)
            Mtv         = push!(Mtv, Mrelv[aux])
        end
        
    end
    
    #Plot xs(t) and xf(t):
    PyPlotFigure(w=w, h=h, top=1.0, bottom=1.0, left=1.5, right=0.2)
    plot(tv, xsv, "b")
    plot(tv, xfv, "r")
    tick_params(axis="both", which="both", labelsize=TickSize)
    xlabel(latexstring("t"), fontsize=10)
    ylabel(latexstring("x_1"), fontsize=10, rotation=0)    
    legend(["shock", "flame"], fontsize=8)
    title(latexstring("\\mathcal{Q}=", q, ", \\beta=", beta, ", \\mathcal{M}_i=", ML), fontsize=10)
    if SaveFig
        savefig("$(FigUbi)SC$(SC)_TrackShockFlame.png", dpi=800, pad_inches=0)
    end
    
    #Plot xf(t)-xs(t):
    if false
        PyPlotFigure(w=w, h=h, top=1.0, bottom=1.0, left=1.5, right=0.2)
        plot(tv, xfv-xsv, "g")
        tick_params(axis="both", which="both", labelsize=TickSize)
        xlabel(latexstring("t"), fontsize=10)
        ylabel(latexstring("x_f-x_s"), fontsize=10, rotation=0)    
        title(latexstring("\\mathcal{Q}=", q, ", \\beta=", beta, ", \\mathcal{M}_i=", ML), fontsize=10)
        if SaveFig
            savefig("$(FigUbi)SC$(SC)_TrackShockFlame2.png", dpi=800, pad_inches=0)
        end
    end
    
    #Plot Mrel(t):
    PyPlotFigure(w=w, h=h, top=1.0, bottom=1.0, left=1.5, right=0.2)
    plot(tv, Mtv, "b")
    tick_params(axis="both", which="both", labelsize=TickSize)
    xlabel(latexstring("t"), fontsize=10)
    ylabel(latexstring("\\mathcal{M}_t"), fontsize=10, rotation=0)    
    title(latexstring("\\mathcal{Q}=", q, ", \\beta=", beta, ", \\mathcal{M}_i=", ML), fontsize=10)
    if SaveFig
        savefig("$(FigUbi)SC$(SC)_Mt.png", dpi=800, pad_inches=0)
    end
    
    return
    
end

function xsxf_ShockFlame(solver::SolverDataSave, GasModel::GasFP,
    pref::Float64, YFref::Float64)

    #Pressure and YF at mesh nodes:
    x1v         = solver.mesh.NodesCoords[:,1]
    pv          = SolutionAtNodes(solver,GasModel,"p")
    YFv         = SolutionAtNodes(solver,GasModel,"Y_F")
        
    #Sort by x1 coordinate:
    sortv       = sortperm(x1v)
    x1v         = x1v[sortv]
    pv          = pv[sortv]
    YFv         = YFv[sortv]
    
    #Get shock and flame positions:
    xf          = NaN
    T           = findlast(YFv.>YFref)
    if T<length(x1v)
        xf      = x1v[T] + (x1v[T+1]-x1v[T])*(YFref-YFv[T])/(YFv[T+1]-YFv[T])
    end
    #
    xs          = NaN
    T           = findlast(pv.<pref)
    if T<length(x1v)
        xs      = x1v[T] + (x1v[T+1]-x1v[T])*(pref-pv[T])/(pv[T+1]-pv[T])
    end
        
    return xs, xf
    
end

function pYref_ShockFlame(SC::Int)

    solver      = GetSolver(SC,0)
    GasModel    = GasFP_ShockFlame(SC,0)
    
    pv          = SolutionAtNodes(solver,GasModel,"p")
    YFv         = SolutionAtNodes(solver,GasModel,"Y_F")
    
    pref        = 0.5*(maximum(pv)+minimum(pv))
    YFref       = 0.5*(maximum(YFv)+minimum(YFv))
   
    return pref, YFref
    
end

#Relative velocity v in x1 direction w.r.t. given velocity v0:
function SolutionAtNodes_vxRel(solver::SolverDataSave, GasModel::GasFP, v0::Float64)

    vx          = SolutionAtNodes(solver, GasModel, "vx")
    return @. vx-v0
    
end
function SolutionAtNodes_MRel(solver::SolverDataSave, GasModel::GasFP, v0::Float64)

    vx          = SolutionAtNodes(solver, GasModel, "vx")
    vy          = SolutionAtNodes(solver, GasModel, "vy")
    a           = SolutionAtNodes(solver, GasModel, "a")
    return @. sqrt((vx-v0)^2 + vy^2)/a
    
end

function MtPlot(q::Float64, beta::Float64; 
    SaveFig::Bool=false, w::Float64=10.0, h::Float64=10.0)

    Miv     = []
    Mt1v    = []
    Mt2v    = []
    gamma   = 1.4
    if q==5.0 && beta==25.0
        Miv     = [2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 10.0]
        Mt1v    = [2.6, 5.2, 5.4, 5.99, 6.65, 8.08, 15.71 ]
        Mt2v    = [2.6, 5.2, 5.4, 5.99, 6.65, 8.08, 15.72 ]
        #Case M1=2.5 is not totally stationary
    elseif q==5.0 && beta==50.0
        Miv     = [2.0, 2.5, 3.0, 3.5, 4.0]
        Mt1v    = [2.59, 3.42, NaN, 5.7, 5.5 ]
        Mt2v    = [2.59, 3.42, NaN, 6.3, 8.5 ]
    elseif q==1.0 && beta==25.0
        Miv     = [2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 10.0, 20.0]
        Mt1v    = [NaN, 3.17, 3.73, 4.30, 4.88, 6.07, 12.09, 24.13 ]
        Mt2v    = [NaN, 3.22, 3.76, 4.33, 4.91, 6.09, 12.10, 24.16 ]
    else
        error("Undefined case")
    end
    
    #Plot xs(t) and xf(t):
    PyPlotFigure(w=w, h=h, top=1.0, bottom=1.5, left=2.2, right=0.2)
    semilogx(Miv, Mt1v./Miv, "sb", markersize=3.0)
    semilogx(Miv, Mt2v./Miv, "sb", markersize=3.0)
    Miv_ex          = logspace(minimum(Miv), maximum(Miv), 100)
    Mtv_ex          = zeros(size(Miv_ex))
    for ii=1:length(Miv_ex)
        Mtv_ex[ii], = MtMrDeton(Miv_ex[ii], gamma, q) 
    end
    semilogx(Miv_ex, Mtv_ex./Miv_ex, color=:orange)
    for ii=1:length(Miv_ex)
        Mtv_ex[ii], = MtMrShock(Miv_ex[ii], gamma, q) 
    end
    semilogx(Miv_ex, Mtv_ex./Miv_ex, color=:cyan)
    
    tick_params(axis="both", which="both", labelsize=TickSize)
    xlabel(latexstring("\\mathcal{M}_i"), fontsize=10)
    ylabel(latexstring("\\mathcal{M}_t/\\mathcal{M}_i"), fontsize=10, rotation=0, labelpad=20.0)    
    title(latexstring("\\mathcal{Q}=", q, ", \\beta=", beta), fontsize=10)
    if SaveFig
        savefig("$(FigUbi)Mt_q$(q)_beta$(beta).png", dpi=800, pad_inches=0)
    end
    
end

function MCJ_(gamma::Float64, q::Float64)

    return sqrt((gamma+1)*q/2)+sqrt(1.0+(gamma+1)*q/2)
    
end

#Based on Napierialski et al (AMM2024):
function MtMrDeton(Mi::Float64, gamma::Float64, q::Float64)

    Rfun(M)     = ((gamma+1)*M^2)/((gamma-1)*M^2+2);
    Pfun(M)     = (2*gamma*M^2-gamma+1)/(gamma+1);
    
    #qq=Q, q=QQ, in Cesar's code
#     qq          = 0.0       
    qq          = q
    
    #Define residual function:
    Ri          = Rfun(Mi)
    Ai          = sqrt(Pfun(Mi)/Ri)
    function ResidualFun!(x::Vector{Float64}, f::Union{Vector{Float64},Nothing}, g::Union{Vector{Float64},Nothing})
        
#         println("Mt=",x[1])
        
        #Transmitted shock:
        Mt      = x[1]
        if Mt<0
            @warn("Negative Mach")
            return -1
        end
        kappa2  = ((Mt^2-1)^2 - 2*(gamma+1)*qq*Mt^2)/Mt^4
        if kappa2<0
            @warn("Negative kappa2")
            return -1
        end
        kappa   = sqrt(kappa2)
        Rt      = (gamma+1)*Mt^2 / ((gamma-kappa)*Mt^2 +1)
        
        #Reflected shock:
        Mr      = sqrt( (3 - gamma - 2*Mi^2 + 2*gamma*Mi^2 + (1+gamma)*(1+kappa)*Mt^2) / 
                        (2*(1+gamma*(2*Mi^2-1))) )
        Rr      = Rfun(Mr)
      
        #Residual:
        f[1]    = Mi*(1.0-1.0/Ri)-Ai*Mr*(1.0-1.0/Rr)-Mt*(1.0-1.0/Rt)/sqrt(1.0+q)    #Eq. (6)
        if !isnothing(g)
            faux    = zeros(1)
            ResidualFun!([Mt+1e-4], faux, nothing)
            J       = (faux[1]-f[1])/1e-4
            J       *= 10.0      #Safety factor to trick Anderson's method
            g[1]    = J\f[1]
        end
        
        return 1
        
    end
    
    #DEBUG:
    if false
        Nplot       = 100
        fplotv      = zeros(Nplot)
        gplotv      = zeros(Nplot)
        Mtplotv     = linspace(0.9*Mi, MCJ_(gamma, q)+Mi, Nplot);
        for ii=1:length(Mtplotv)
            faux    = zeros(1)
            gaux    = zeros(1)
            flag    = ResidualFun!([Mtplotv[ii]], faux, gaux)
            if flag<0
                fplotv[ii]  = NaN
                gplotv[ii]  = NaN
            else
                fplotv[ii]  = faux[1]
                gplotv[ii]  = gaux[1]
            end
        end
        figure()
        plot(Mtplotv, fplotv, "b")
        plot(Mtplotv, gplotv, "r")
        axis("equal")
        grid("on")
        display(MCJ_(gamma, q) + Mi)
        println("")
    #     return
    end
    
    #Initial conditions:
    Mcj         = MCJ_(gamma, q)
    Mt0         = Mi+Mcj
    
    #Call Anderson's method:
    NLS_output  = Anderson(FW_NLS((u,f,g)->ResidualFun!(u,f,g)), 
                    [Mt0], 
                    AbsTolX=0.0, RelTolX=0.0, 
                    AbsTolF=0e-14, RelTolF=0.0,
                    AbsTolG=1e-4, RelTolG=0.0,
                    Display="notify", MaxIter=100, memory=0)
    Mt          = NLS_output[1][1]
    kappa       = sqrt(((Mt^2-1)^2 - 2*(gamma+1)*qq*Mt^2)/Mt^4)
    Mr          = sqrt( (3 - gamma - 2*Mi^2 + 2*gamma*Mi^2 + (1+gamma)*(1+kappa)*Mt^2) / 
                        (2*(1+gamma*(2*Mi^2-1))) )
    if NLS_output[2].flag<0
        Mt      = NaN
        Mr      = NaN
    end
    
    return Mt, Mr
    
end

#Based on Napierialski et al (AMM2024):
function MtMrShock(Mi::Float64, gamma::Float64, q::Float64)

    Rfun(M)     = ((gamma+1)*M^2)/((gamma-1)*M^2+2);
    Pfun(M)     = (2*gamma*M^2-gamma+1)/(gamma+1);
    
    #Define residual function:
    Ri          = Rfun(Mi)
    Ai          = sqrt(Pfun(Mi)/Ri)
    function ResidualFun!(x::Vector{Float64}, f::Union{Vector{Float64},Nothing}, g::Union{Vector{Float64},Nothing})
    
        Mt      = x[1]
        Mr      = x[2]
        Rt      = Rfun(Mt)
        Rr      = Rfun(Mr)
        f[1]    = Mi*(1.0-1.0/Ri)-Ai*Mr*(1.0-1.0/Rr)-Mt*(1.0-1.0/Rt)/sqrt(1.0+q)            #Eq. (6)
        f[2]    = Mr^2*(1.0+gamma*(2.0*Mi^2-1)) - Mt^2*(gamma+1)-Mi^2*(gamma-1)+(gamma-1)   #Eq. (7)
        if !isnothing(g)
            J       = [ -(1-1/Rt)/sqrt(1+q)     -Ai*(1-1/Rr)    
                        -2*Mt*(gamma+1)         2*Mr*(1+gamma*(2*Mi^2-1)) ]
            g       .= J\f
#             g       .= f
        end
        
        return 1
        
    end
    
#     fv      = zeros(2)
#     ResidualFun!([2.58,1.26], fv, nothing)
#     display(fv)
#     return
    
    #Initial conditions:
    Mt0         = Mi
    Mr0         = sqrt( (Mt0^2*(gamma+1)+Mi^2*(gamma-1)+(gamma-1))/(1.0+gamma*(2.0*Mi^2-1)) )
    
    #Call Anderson's method:
    NLS_output  = Anderson(FW_NLS((u,f,g)->ResidualFun!(u,f,g)), 
                    [Mt0, Mr0], 
                    AbsTolX=0.0, RelTolX=0.0, 
                    AbsTolF=0e-14, RelTolF=0.0,
                    AbsTolG=1e-4, RelTolG=0.0,
                    Display="notify", MaxIter=100)
    Mt          = NLS_output[1][1]
    Mr          = NLS_output[1][2]
    if NLS_output[2].flag<0
        Mt      = NaN
        Mr      = NaN
    end
    
    return Mt, Mr
    
end

#NOTAS:
#   -¿Se podría definir de forma más rigurosa la posición de la llama y de la onda?
#   -¿Sería necesario hacerlo para reducir las oscilaciones?
#   -Faltaría modificar el código para que la onda esté quieta.
#   -Con beta=50.0 los resultados no parecen ser estacionarios.
#   -Hay casos en los que la llama se separa de la onda al principio y luego
#   se junta (utilizando el modelo de Arrhenius). Sería importante limitar el dominio
#   o el tiempo.
