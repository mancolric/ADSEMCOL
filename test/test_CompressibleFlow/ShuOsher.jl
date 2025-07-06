cd(@__DIR__)
include("test_CompressibleFlow.jl")

function ShuOsher(hp0::Float64, FesOrder::Int;
    tf::Float64=1.8, RKMethod::String="Ascher3",
    epsilon::Float64=1e-4, nu::Float64=0e-3, beta::Float64=0.0, kappa_rho_cv::Float64=0e-6, 
    delta0::Float64=max(epsilon,nu), WaveAmpl::Float64=0.2, 
    #
    TolS::Float64=1e-4, AMA_MaxIter::Int=200, AMA_SizeOrder::Int=FesOrder, AMA_AnisoOrder::Int=2, 
    #
    TolT::Float64=1e-4, Deltat0::Float64=1e-5,
    #
    PlotFig::Bool=false, wFig::Float64=9.50, hFig::Float64=6.50, 
    PlotVars::Vector{String}=String[], PlotCode::Vector{String}=fill("nodes", length(PlotVars)), 
    SaveFig::Bool=false, Nt_SaveFig::Int=5, Deltat_SaveFig::Float64=Inf,
    mFig::Int=max(1,length(PlotCode)), nFig::Int=Int(ceil(length(PlotCode)/mFig)), cmap::String="jet",
    #
    SaveRes::Bool=false, Nt_SaveRes::Int=typemax(Int), Deltat_SaveRes::Float64=0.01, 
    #
    CSS::Float64=0.1, 
    #
    SC::Int=0)

    #---------------------------------------------------------------------
    #PROBLEM DATA:
    
    #Problem data:
    gamma           = 1.4
    xdisc           = -4.0
    ML              = 3.0
    rhoL            = 3.857143
    uL              = 2.629369
    pL              = 10.33333
    rhoEL           = pL/(gamma-1.0)+0.5*rhoL*uL^2
    pR              = 1.0
    
    #Define gas:
    GasModel                = GasIdeal()
    GasModel.gamma          = gamma
    GasModel.epsilon        = epsilon
    GasModel.nu             = nu
    GasModel.beta           = beta
    GasModel.kappa_rho_cv   = kappa_rho_cv
    GasModel.CSS            = CSS

    #Initial conditions:
    function u0fun(x::Vector{Matrix{Float64}})
        
        rho         = @. SmoothHeaviside(x[1]-xdisc, delta0, rhoL, 1.0) + 
                        WaveAmpl*sin(5.0*x[1])*SmoothHeaviside(x[1]-xdisc, delta0, 0.0, 1.0)
        rhovx       = @. rho * SmoothHeaviside(x[1]-xdisc, delta0, uL, 0.0)
        rhovy       = @. 0.0*x[1]
        p           = @. SmoothHeaviside(x[1]-xdisc, delta0, pL, pR)
        rhoE        = @. p/(gamma-1.0) + 0.5*rhovx^2/rho
        
        return [rho, rhovx, rhovy, rhoE]
        
    end

    #Boundary conditions:
    BC_walls        = SlipAdiabatic()
    function pOutlet(t::Float64, x::Vector{Matrix{Float64}})
        
        p           = @mlv 0.0*x[1]+pR
        return [p]
        
    end
    BC_outlet       = SubsonicOutlet1(FWt11((t,x)->pOutlet(t,x)))
    function uInlet(t::Float64, x::Vector{Matrix{Float64}})
    
        return u0fun(x)
        
    end
    BC_inlet        = SupersonicInlet1(FWt11((t,x)->uInlet(t,x)))
    
    #---------------------------------------------------------------------
    #PRE-PROCESS STAGE:
    
    #Mesh:
    MeshFile                = "../temp/ShuOsher_SC$(SC).geo"
    NX                      = Int(ceil(10.0/(hp0*FesOrder)))
    NY                      = Int(ceil(1.0/(hp0*FesOrder)))
    TrMesh_Rectangle_Create!(MeshFile, -5.0, 5.0, NX, -0.5, 0.5, NY)
    
    #Load LIRKHyp solver structure with default data. Modify the default data if necessary:
    solver                  = LIRKHyp_Start(GasModel)
    solver.ProblemName      = "ShuOsher"
    solver.SC               = SC
    solver.MeshFile         = MeshFile
    solver.nBounds          = 4
    solver.FesOrder         = FesOrder
    solver.RKMethod         = RKMethod
    solver.Deltat0          = Deltat0
    solver.tf               = tf
    solver.AMA_MaxIter      = AMA_MaxIter
    solver.AMA_SizeOrder    = AMA_SizeOrder
    solver.AMA_AnisoOrder   = AMA_AnisoOrder
    solver.TolS_max         = TolS
    solver.TolS_min         = 0.01*TolS
    solver.TolT             = TolT
#     solver.CA_max           = Inf
    
    #Set initial and boundary conditions:
    solver.u0fun        = FW11((x) -> u0fun(x)) 
    solver.BC           = [BCW(BC_walls), BCW(BC_outlet), BCW(BC_walls), BCW(BC_inlet)]
    
    #-----------------------------------------------------------------------------
    #INITIAL CONDITION:
    
    #Compute initial condition:
    ConvFlag            = LIRKHyp_InitialCondition!(solver)

    #Function to plot solution:
    figv                = Vector{Figure}(undef,3)
    if PlotFig
        figv[1]         = PyPlotSubPlots(mFig, nFig, w=wFig, h=hFig, left=0.9, right=0.4, bottom=1.1, top=1.0)
        for ii=2:length(figv)
            figv[ii]    = figure()
        end
    end    
    t_lastFig           = 0.0
    ct_SaveFig          = 0
    nb_SaveFig          = 0
    function PlotSol()
    
        ct_SaveFig      += 1
        
        if PlotFig && ( solver.t-t_lastFig>=Deltat_SaveFig || 
                        ct_SaveFig==Nt_SaveFig || solver.t==tf || solver.t==0.0 )
        
            figure(figv[1].number)
    
            #Loop plot variables:
            for ii=1:length(PlotVars)
                PyPlot.subplot(mFig, nFig, ii)
                PyPlot.cla()

                splot_fun(x1,x2)    = @mlv x1
                PlotNodes(splot_fun, solver, PlotVars[ii])
                xlabel(latexstring("x_1"), fontsize=10)
                title(latexstring(LatexString(PlotVars[ii]),
                    "; t^n=", sprintf1("%.2e", solver.t)), 
                    fontsize=10)
            end
            if SaveFig
                savefig("$(VideosUbi)ShuOsher_SC$(SC)_$(nb_SaveFig).png", dpi=400, pad_inches=0)
            end
        
            figure(figv[2].number)
            PyPlot.cla()
            plot(solver.mesh)
            if SaveFig
                savefig("$(VideosUbi)ShuOsher_Mesh$(SC)_$(nb_SaveFig).png", dpi=400, pad_inches=0)
            end
            
            figure(figv[3].number)
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
            if SaveFig && solver.t==tf
                savefig("$(VideosUbi)ShuOsher_Errors$(SC)_$(nb_SaveFig).png", dpi=400, pad_inches=0)
            end
            
            t_lastFig           += Deltat_SaveFig
            ct_SaveFig          = 0
            nb_SaveFig          += 1
            
        end
        return
        
    end
    PlotSol()
    
    #Function to save intermediate results:
    t_lastRes           = -Deltat_SaveRes
    ct_SaveRes          = 0
    nb_SaveRes          = 0
    function SaveSol()
        
        ct_SaveRes      += 1 
        if SaveRes && ( solver.t-t_lastRes>=Deltat_SaveRes || 
                        ct_SaveRes==Nt_SaveRes || solver.t==tf || solver.t==0.0 )
            save("$(ResUbi)LIRKHyp_SC$(SC)_$(nb_SaveRes).jld2", "StudyCase", "ShuOsher",
                "ConvFlag", ConvFlag, "solver", save(solver), 
                "epsilon", epsilon, "nu", nu, "beta", beta, "kappa_rho_cv", kappa_rho_cv,
                "delta0", delta0, "WaveAmpl", WaveAmpl, "TolS", TolS, "TolT", TolT)
            t_lastRes   += Deltat_SaveRes
            ct_SaveRes  = 0
            nb_SaveRes  += 1
        end
        return
        
    end
    SaveSol()
    
    #-----------------------------------------------------------------------------
    #MARCH IN TIME:
    
    while solver.t<tf
    
        ConvFlag    = LIRKHyp_Step!(solver)
        
        if ConvFlag<=0
            break
        end
        
        PlotSol()
        SaveSol()
        
    end
    
    return solver
    
end
