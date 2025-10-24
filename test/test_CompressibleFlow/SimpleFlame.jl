include("test_CompressibleFlow.jl")

function SimpleFlame(MeshCase::Int, FesOrder::Int, tfv::Array{Float64,1};
    RKMethod::String="Ascher3", epsilon::Float64=0.0,
    Re::Float64=1e4, Pr::Float64=0.73, beta::Float64=10.0, q::Float64=5.0,
    Le::Float64=1.0, Da::Float64=1e3,
    PlotVars::Vector{String}=String[], PlotCode::Vector{String}=fill("nodes", length(PlotVars)), 
    SaveFig::Bool=false, wFig::Float64=9.50, hFig::Float64=6.50, 
    mFig::Int=max(1,length(PlotCode)), nFig::Int=Int(ceil(length(PlotCode)/mFig)), Nt_SaveFig::Int=1, cmap::String="jet",
    SC::Int=0, CSS::Float64=0.1, CDC::Float64=5.0, Deltat0::Float64=1e-4,
    AMA_MaxIter::Int=20, TolS::Float64=1e-3, TolT::Float64=1e-3, AMA_SizeOrder::Int=2)

    #---------------------------------------------------------------------
    #PROBLEM DATA:
    
    #Data:
    rho0            = 1.0
    RT0             = 1.0               #Initial temperature (T) (multiplied by Rg)
    RTf             = (1.0+q)*RT0       #Flame temperature
    YF0             = 0.02
    gamma           = 1.4
    Q0              = q*gamma/(gamma-1)*RT0/YF0
    
    #Define gas:
    GasModel                = GasFP()
    GasModel.gamma          = gamma
    GasModel.epsilon        = epsilon
    GasModel.nu             = 1.0/Re
    GasModel.beta           = 0.0
    GasModel.kappa_rho_cv   = GasModel.nu/Pr * GasModel.gamma
    GasModel.D              = GasModel.kappa_rho_cv/GasModel.gamma/Le
    GasModel.B              = Da*exp(beta)
    GasModel.RTa            = beta*RTf
    GasModel.hfF            = Q0
    GasModel.hfP            = 0.0
    GasModel.CSS            = CSS

    #Initial conditions:
    function u0fun(x::Vector{Matrix{Float64}})
    
        YF0_fun         = (x) -> YF0*SmoothHeaviside(0.7-x, 0.02)
        RT0_fun         = (x) -> RT0 + (RTf-RT0)*SmoothHeaviside(x-0.7, 0.02)
        
        rhoYF           = @mlv rho0*YF0_fun(x[1])
        rhoYP           = @mlv rho0*(1.0-YF0_fun(x[1]))
        rhovx           = @mlv 0.0*x[1]
        rhovy           = @mlv 0.0*x[1]
        rhoE            = @mlv rho0*(RT0_fun(x[1])/(GasModel.gamma-1.0) + Q0*YF0_fun(x[1]))
        
        return [rhoYF, rhoYP, rhovx, rhovy, rhoE]
        
    end
    
    #Boundary conditions:
    BC_walls        = SlipAdiabatic()
    function pOutlet(t::Float64, x::Vector{Matrix{Float64}})
        
        p0          = rho0*RT0
        p           = @mlv 0.0*x[1]+p0
        return [p]
        
    end
    BC_outlet       = SubsonicOutlet1(FWt11((t,x)->pOutlet(t,x)))
    
    #---------------------------------------------------------------------
    #PRE-PROCESS STAGE:
    
    #Load LIRKHyp solver structure with default data. Modify the default data if necessary:
    solver              = LIRKHyp_Start(GasModel)
    solver.ProblemName  = "SimpleFlame"
    solver.SC           = SC
    solver.MeshFile     = "$(MeshUbi)SodStandard/MeshCase$(MeshCase).geo"
    solver.nBounds      = 4
    solver.FesOrder     = FesOrder
    solver.RKMethod     = RKMethod
    solver.Deltat0      = Deltat0
    solver.tf           = tfv[length(tfv)]
    solver.AMA_MaxIter  = AMA_MaxIter
    solver.AMA_SizeOrder= AMA_SizeOrder
    solver.TolS_max     = TolS
    solver.TolS_min     = 0.01*TolS
    solver.TolT         = TolT
    
    #Set initial and boundary conditions:
    solver.u0fun        = FW11((x) -> u0fun(x)) 
    solver.BC           = [BCW(BC_walls), BCW(BC_walls), BCW(BC_walls), BCW(BC_outlet)]
    
    #-----------------------------------------------------------------------------
    #INITIAL CONDITION:
    
    #Compute initial condition:
    LIRKHyp_InitialCondition!(solver)

    #Open plots (variables, errors):
    figv                = Vector{Figure}(undef,3)
    if !isempty(PlotVars)
        figv[1]         = PyPlotSubPlots(mFig, nFig, w=wFig, h=hFig, left=0.9, right=0.4, bottom=1.1, top=1.0)
        figv[2]         = figure()
        figv[3]         = figure()
    end    
    t_lastFig           = 0.0
    NPlotSave           = 0
    NSave               = 0
    function PlotSol()
        if isempty(PlotVars)
            return
        end
        figure(figv[1].number)
        for ii=1:length(PlotVars)
            PyPlot.subplot(mFig, nFig, ii)
            PyPlot.cla()
#             plot3D(solver, PlotVars[ii])
            splot_fun(x1,x2)    = @mlv x1
            PlotNodes(splot_fun, solver, PlotVars[ii])
            xlabel(latexstring("x_1"), fontsize=10)
            title(latexstring(LatexString(PlotVars[ii]),
                "; t^n=", sprintf1("%.2e", solver.t)), 
                fontsize=10)
        end
        if SaveFig
            savefig("$(VideosUbi)SimpleFlame_SC$(SC)_$(NSave).png", dpi=400, pad_inches=0)
        end
        
        figure(figv[2].number)
        PyPlot.cla()
        plot(solver.mesh)
        axis("equal")
        if SaveFig
            savefig("$(VideosUbi)SimpleFlame_Mesh$(SC)_$(NSave).png", dpi=400, pad_inches=0)
        end
        
        figure(figv[3].number)
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
        if SaveFig && solver.t==tfv[length(tfv)]
            savefig("$(VideosUbi)SimpleFlame_Errors$(SC).png", dpi=400, pad_inches=0)
        end
        
        NSave   += 1
        
        sleep(0.001)
        return
        
    end
    PlotSol()
    
    #-----------------------------------------------------------------------------
    #MARCH IN TIME:
    
    while solver.t<tfv[length(tfv)]
    
        flag    = LIRKHyp_Step!(solver)
        
        if flag<=0
            break
        end
        
        NPlotSave   += 1
        if solver.t-t_lastFig>=1e-2 || solver.t==tfv[length(tfv)] || NPlotSave==Inf
            NPlotSave   = 0
            t_lastFig   = solver.t
            PlotSol()
        end
        
        
    end
    
    return solver
    
end
