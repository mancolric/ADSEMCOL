include("PlotResults.jl")

function GasFXP_ShockFlame(SC::Int, nb::Int)

    FileName                = GetFileName(SC, nb)
    
    #Define gas:
    GasModel                = GasFXP()
    GasModel.epsilon        = load(FileName, "epsilon")
    GasModel.gamma          = load(FileName, "gamma")
    GasModel.nu             = load(FileName, "nu")
    GasModel.beta           = 0.0
    GasModel.kappa_rho_cv   = load(FileName, "kappa_rho_cv")
    GasModel.hfF            = 0.0
    GasModel.hfX            = 0.0
    GasModel.hfP            = -load(FileName, "Q0")
    GasModel.D              = load(FileName, "D")
    GasModel.BR             = load(FileName, "BR")
    GasModel.BI             = load(FileName, "BI")
    GasModel.BB             = load(FileName, "BB")
    GasModel.RTI            = load(FileName, "RTI")
    GasModel.RTB            = load(FileName, "RTB")
    
    return GasModel
    
end
    
function x1Plot_ShockFlame(SC::Int, nb::Int;
    PlotVars::Vector{String}=["Y_F", "v", "p", "RT", "rho", "M"],
    MarkSF::Bool=false,
    mFig::Int=3, nFig::Int=2, 
    SaveFig::Bool=false, w::Float64=9.50, h::Float64=6.50)
    
    #Extract data:
    if nb==1000
        nb      = load("$(ResUbi)LIRKHyp_SC$(SC)_info.jld2", "nb_SaveRes")
    end
    GasModel    = GasFXP_ShockFlame(SC,nb)
    FileName    = GetFileName(SC, nb)
    solver      = GetSolver(SC, nb)
    q           = load(FileName, "q")
    ML          = load(FileName, "ML")
    
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
            vbleString      = LatexString(GasModel, PlotVars[ii])
        end
        tick_params(axis="both", which="both", labelsize=TickSize)
        xlabel(latexstring("x_1"), fontsize=10)
        ylabel(latexstring(LatexString(GasModel, PlotVars[ii])), rotation=0, fontsize=10, labelpad=10)
            
        #Get positions of shock and flame:
        if MarkSF
            xsv         = load(FileName, "xsv")
            xfv         = load(FileName, "xfv")
            xs          = xsv[length(xsv)]
            xf          = xfv[length(xfv)]
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
    suptitle(latexstring("\\mathcal{Q}=", q, ", M=", ML, 
                        ", t=", sprintf1("%.2e", solver.t)), fontsize=10)
    
    if SaveFig
        savefig("$(FigUbi)SC$(SC)_$(nb)_x1Plot.png", dpi=800, pad_inches=0)
    end
        
    return
    
end

function Contour_ShockFlame(SC::Int, nb::Int; SaveFig::Bool=false, w::Float64=8.50, h::Float64=8.50, 
    PlotVars::Vector{String}=["Y_F", "v", "p", "RT", "rho", "M"], mFig::Int=3, nFig::Int=2, 
    PlotTitle::Bool=false)

    GasModel    = GasFXP_ShockFlame(SC, nb)
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
            title(latexstring(LatexString(GasModel, PlotVars[ii]),
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
    ML          = load(FileName, "ML")
    GasModel    = GasFXP_ShockFlame(SC,nbMax)
    
    tPlot_nFig  = tPlot/(nbMax-n1)
    
    #Shock and flame locations:
    tv          = load(FileName, "tv")
    xsv         = load(FileName, "xsv")
    xfv         = load(FileName, "xfv")
    xsdotv      = xdotfun(tv, xsv)
    
#     fig         = PyPlotSubPlots(mFig, nFig, w=w, h=h, top=1.0)
    fig         = PyPlotSubPlots(mFig, nFig, w=w, h=h, top=0.1, right=0.1, left=1.8, bottom=0.1, tTitle=0.8, bTitle=1.0)
    kSave       = 0
    for nb=n1:nbMax
    
        t0          = time()
        
        FileName    = GetFileName(SC, nb)
        solver      = GetSolver(SC, nb)
        
        #Get positions of shock and flame
        Nt          = length(load(FileName, "tv"))
        xs          = xsv[Nt]
        xf          = xfv[Nt]
        xsdot       = xsdotv[Nt]
        
        #Plot:
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
                vbleString      = LatexString(GasModel, PlotVars[ii])
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
        suptitle(latexstring("\\mathcal{Q}=", q, ", M_i=", ML, 
                    ", t=", sprintf1("%.2e", solver.t)), fontsize=10)
        
        if SaveFig
            savefig("$(FigUbi)SC$(SC)_$(kSave).png", dpi=800, pad_inches=0)
#             error("")
        end
        kSave       += 1
        
        t1          = time()
        
        pause(min(max(0.01, tPlot_nFig-(t1-t0)),0.5))
        
    end
    
    return
    
end

function Track_ShockFlame(SC::Int; SaveFig::Bool=false, w::Float64=10.0, h::Float64=10.0, 
    Deltat::Float64=1.0, nbMax::Int=-1)

    if nbMax<0
        nbMax   = load("$(ResUbi)LIRKHyp_SC$(SC)_info.jld2", "nb_SaveRes")
    end
    
    FileName    = GetFileName(SC, nbMax)
    q           = load(FileName, "q")
    ML          = load(FileName, "ML")
    tv          = load(FileName, "tv")
    xsv         = load(FileName, "xsv")
    xfv         = load(FileName, "xfv")
    gamma       = load(FileName, "gamma")
#     Mtv         = load(FileName, "Mrelv")

    #Relative Mach number:
    u0              = ML*sqrt(gamma*(1.0+q))    #u0=uL, RT0=1, RTL=RT0*(1+q)
    a0              = sqrt(gamma*1.0)
    xsdotv          = xdotfun(tv, xsv, Deltat=Deltat)
    Mtv             = @. (u0-xsdotv)/a0
#     Mtv[xfv.<xsv]   .= NaN
    
    #Plot xs(t) and xf(t):
    PyPlotFigure(w=w, h=h, top=1.0, bottom=1.0, left=1.5, right=0.2)
    plot(tv, xsv, "b")
    plot(tv, xfv, "r")
    tick_params(axis="both", which="both", labelsize=TickSize)
    xlabel(latexstring("t"), fontsize=10)
    ylabel(latexstring("x_1"), fontsize=10, rotation=0)    
    legend(["shock", "flame"], fontsize=8)
    title(latexstring("\\mathcal{Q}=", q, ", \\mathcal{M}_i=", ML), fontsize=10)
    if SaveFig
        savefig("$(FigUbi)SC$(SC)_TrackShockFlame.png", dpi=800, pad_inches=0)
    end
    
    #Plot xsdot:
    PyPlotFigure(w=w, h=h, top=1.0, bottom=1.0, left=1.5, right=0.2)
    plot(tv, xsdotv, "g")
    tick_params(axis="both", which="both", labelsize=TickSize)
    xlabel(latexstring("t"), fontsize=10)
    #         ylabel(latexstring("x_f-x_s"), fontsize=10, rotation=0)    
    ylabel(latexstring("\\dot{x}_s"), fontsize=10, rotation=0)
    ylim(-0.1, 1.2*maximum(xsdotv[.!isnan.(xsdotv)]))
    title(latexstring("\\mathcal{Q}=", q, ", \\mathcal{M}_i=", ML), fontsize=10)
    grid("on")
    if SaveFig
        savefig("$(FigUbi)SC$(SC)_TrackShockFlame2.png", dpi=800, pad_inches=0)
    end
    
    #Plot Mrel(t):
    PyPlotFigure(w=w, h=h, top=1.0, bottom=1.0, left=1.5, right=0.2)
    plot(tv, Mtv, "b")
    tick_params(axis="both", which="both", labelsize=TickSize)
    xlabel(latexstring("t"), fontsize=10)
    ylabel(latexstring("\\mathcal{M}_t"), fontsize=10, rotation=0, labelpad=10.0)    
    title(latexstring("\\mathcal{Q}=", q, ", \\mathcal{M}_i=", ML), fontsize=10)
    ylim(-0.1, 1.2*maximum(Mtv[.!isnan.(Mtv)]))
    grid("on")
    if SaveFig
        savefig("$(FigUbi)SC$(SC)_Mt.png", dpi=800, pad_inches=0)
    end
    println("M_rel=", Mtv)
    
    return
    
end

#Relative velocity v in x1 direction w.r.t. given velocity v0:
function SolutionAtNodes_vxRel(solver::SolverDataSave, GasModel::GasFXP, v0::Float64)

    vx          = SolutionAtNodes(solver, GasModel, "vx")
    return @. vx-v0
    
end
function SolutionAtNodes_MRel(solver::SolverDataSave, GasModel::GasFXP, v0::Float64)

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
    if q==5.0 && beta==10.0
        Miv     = [ 2.0,    2.4,    2.8,    2.86,   2.92, 
                    2.99,   3.05,   3.11,   3.18,   3.26, 
                    3.33,   3.4, 
                    4.0,    4.9,    5.8,    6.9,    8.4, 
                    10.0    ]
        Mt1v    = [ 2.58,   3.26,   3.95,   4.05,   4.16, 
                    4.29,   4.41,   5.56,   5.60,   5.69, 
                    5.78,   5.86, 
                    6.65,   7.93,   9.27,   10.94,  13.24, 
                    15.71   ] 
        Mt2v    = copy(Mt1v)
    else
        error("Undefined case")
    end
    #OLD CASES:
    #=
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
    =#
    
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
    #
    tick_params(axis="both", which="both", labelsize=TickSize)
    xlabel(latexstring("\\mathcal{M}_i"), fontsize=10)
    ylabel(latexstring("\\mathcal{M}_t/\\mathcal{M}_i"), fontsize=10, rotation=0, labelpad=20.0)    
    title(latexstring("\\mathcal{Q}=", q, ", T_I/T_s^{CJ}=T_B/T_s^{CJ}=", beta), fontsize=10)
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
    NLS_output  = Anderson(FW_NLS((u,g)->ResidualFun!(u,[0.0],g)), 
                    [Mt0], 
                    AbsTolX=0.0, RelTolX=0.0, 
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
    NLS_output  = Anderson(FW_NLS((u,g)->ResidualFun!(u,[0.0, 0.0],g)), 
                    [Mt0, Mr0], 
                    AbsTolX=0.0, RelTolX=0.0, 
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

#Compute shock and flame velocities:
function xdotfun(tv::Vector{Float64}, xv::Vector{Float64};
    Deltat::Float64=1.0)
#     return xdotfun_nofilter(tv, xv)
#     return xdotfun_filter1(tv, xv)
#     return xdotfun_filter2(tv, xv)
#     return xdotfun_filter3(tv, xv)
    return xdotfun_filter4(tv, xv; Deltat=Deltat)
end

function xdotfun_nofilter(tv::Vector{Float64}, xv::Vector{Float64})

    N       = length(tv)
    if N==1
        return [ 0.0 ]
    end
    
    xdotv           = zeros(N)
    xdotv[1]        = (xv[2]-xv[1])/(tv[2]-tv[1])
    xdotv[2]        = (xv[2]-xv[1])/(tv[2]-tv[1])
    for ii=3:N
        xdotv[ii]   = (3*xv[ii]-4*xv[ii-1]+xv[ii-2])/(3*tv[ii]-4*tv[ii-1]+tv[ii-2])
#         xdotv[ii]   = (xv[ii]-xv[ii-1])/(tv[ii]-tv[ii-1])
    end
    
    return xdotv
    
end

function xdotfun_filter1(tv::Vector{Float64}, xv::Vector{Float64})

    N               = length(tv)
    if N==1
        return [ 0.0 ]
    end
    
    #Evaluate mean position at t^(n+1/2):
    t12v            = 0.5*(tv[1:N-1]+tv[2:N])
    x12v            = 0.5*(xv[1:N-1]+xv[2:N])
    
    #Velocities from linear interpolation of mean position:
    xdotv           = zeros(N)
    xdotv[1]        = NaN
    for ii=2:N-1
        xdotv[ii]   = (x12v[ii]-x12v[ii-1])/(t12v[ii]-t12v[ii-1])
    end
    xdotv[N]        = NaN
    
    return xdotv
    
end

function xdotfun_filter2(tv::Vector{Float64}, xv::Vector{Float64}; 
    DeltaN::Int=10, Nfilter::Int=1)
    
    N       = length(tv)
    if N==1
        return [ 0.0 ]
    end
    
    #Compute velocities with second order formula. The result is noisy:
    xdotv           = zeros(N)
    xdotv[1]        = (xv[2]-xv[1])/(tv[2]-tv[1])
    xdotv[2]        = (xv[2]-xv[1])/(tv[2]-tv[1])
    for ii=3:N
        xdotv[ii]   = (3*xv[ii]-4*xv[ii-1]+xv[ii-2])/(3*tv[ii]-4*tv[ii-1]+tv[ii-2])
    end
    xdot0v          = 1.0*xdotv
    
    #Filter velocities:
    for ifilter=1:Nfilter
        mass            = 0.0
        inertia1        = 0.0
        tbar            = 0.0
        xdothatv        = 0.0*xdotv
        for ii=2:1:N
        
            #Update length and barycenter and inertia:
            i0              = max(1, ii+1-DeltaN)
            L_new           = tv[ii]-tv[i0]
            tbar_new        = 0.5*(tv[i0]+tv[ii])
            inertia1_new    = inertia1 + mass*(tbar-tbar_new)
            mass_new        = mass
            
            #Update mass and momentum respect to center of time interval [t_(i-1), t_i] 
            #with new element:
            h               = tv[ii]-tv[ii-1]
            Delta_mass      = 0.5*h*xdotv[ii-1] + 0.5*h*xdotv[ii]
            mass_new        += Delta_mass
            inertia1_new    += h*h/12*(xdotv[ii]-xdotv[ii-1]) + 
                                Delta_mass*(0.5*(tv[ii]+tv[ii-1])-tbar_new)
            
            #Update mass and momentum respect to center of time interval [t_(i-1), t_i] 
            #with outgoing element:
            if ii>DeltaN
                h               = tv[ii-DeltaN+1]-tv[ii-DeltaN]
                Delta_mass      = 0.5*h*xdotv[ii-DeltaN+1] + 0.5*h*xdotv[ii-DeltaN]
                mass_new        -= Delta_mass
                inertia1_new    -= h*h/12*(xdotv[ii-DeltaN+1]-xdotv[ii-DeltaN]) + 
                                    Delta_mass*(0.5*(tv[ii-DeltaN+1]+tv[ii-DeltaN])-tbar_new)
            end
            
            #Seek linear xdothat such that mass and momentum are the same in
            #the interval [tv[1], tv[ii]]:
            #   0.5*L*v1 + 0.5*L*v2     = mass_new
            #   -L^2/12*v1 + L^2/12*v2  = inertia1_new
            #v1 and v2 are the values of xdothat at tv[1] and tv[ii]
            #   v1 + v2                 = mass_new/(0.5*L)
            #   -v1 + v2                = inertia1_new*12/L^2
            xdothatv[ii]    = 0.5*( mass_new/(0.5*L_new) + inertia1_new*12/L_new^2 )
        
            #Update mass and inertia1:
            mass            = mass_new
            inertia1        = inertia1_new
            tbar            = tbar_new
            
        end
        xdotv               .= xdothatv
    end
    
    #NaNs for first time steps:
    xdotv[1:DeltaN]     .= NaN
        
#     figure()
#     plot(tv, xdot0v, "x-g")
#     plot(tv, xdotv, "+-r")
#     ylim(-0.5, 2.5)
#     error("")
    
    return xdotv
    
end

function xdotfun_filter3(tv::Vector{Float64}, xv::Vector{Float64};
    DeltaN::Int=50)
    
    #Compute velocities with first order formula. Compute average of 
    #DeltaN last time steps:
    N               = length(xv)
    xdotv           = zeros(N)
    xdotv[1]        = NaN
    for ii=2:N
        jmax            = min(DeltaN,ii-1)
        for jj=1:jmax
            xdotv[ii]   += (xv[ii]-xv[ii-jj])/(tv[ii]-tv[ii-jj])
        end
        xdotv[ii]       /= jmax
    end
        
    return xdotv
    
end

function xdotfun_filter4(tv::Vector{Float64}, xv::Vector{Float64};
    Deltat::Float64=1.0)
    
    #Compute velocities with first order formula. Take a minimum time step
    #given by Deltat:
    N               = length(xv)
    xdotv           = zeros(N)
    xdotv[1]        = NaN
    for ii=2:N
        
        #Find maximum jj such that tv[ii]-tv[jj]>Deltat
        jj          = ii-1
        while true
            if tv[ii]-tv[jj]>Deltat
                break
            end
            if jj==1
                break
            end
            jj      -= 1
        end
            
        #Compute velocity:
        xdotv[ii]   = (xv[ii]-xv[jj])/(tv[ii]-tv[jj])
        
    end
        
    return xdotv
    
end

#NOTAS:
#   -¿Se podría definir de forma más rigurosa la posición de la llama y de la onda?
#   -¿Sería necesario hacerlo para reducir las oscilaciones?
#   -Faltaría modificar el código para que la onda esté quieta.
#   -Con beta=50.0 los resultados no parecen ser estacionarios.
#   -Hay casos en los que la llama se separa de la onda al principio y luego
#   se junta (utilizando el modelo de Arrhenius). Sería importante limitar el dominio
#   o el tiempo.

#=
function pYref_ShockFlame(SC::Int)

    solver      = GetSolver(SC,0)
    GasModel    = GasFXP_ShockFlame(SC,0)
    
    pv          = SolutionAtNodes(solver,GasModel,"p")
    YFv         = SolutionAtNodes(solver,GasModel,"Y_F")
    
    pref        = 0.5*(maximum(pv)+minimum(pv))
    YFref       = 0.5*(maximum(YFv)+minimum(YFv))
   
    return pref, YFref
    
end

function xsxf_ShockFlame(solver::SolverDataSave, GasModel::GasFXP,
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

=#
