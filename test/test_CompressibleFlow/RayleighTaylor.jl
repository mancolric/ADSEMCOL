cd(@__DIR__)
include("test_CompressibleFlow.jl")

#NOTE: This is a diffusion-dominated problem. Hence, delta = sqrt(epsilon*tf)

function RayleighTaylor(hp0::Float64, FesOrder::Int;
    tf::Float64=1.95, RKMethod::String="Ascher3",
    epsilon::Float64=1e-8, nu::Float64=0e-3, beta::Float64=0.0, kappa_rho_cv::Float64=0e-6, 
    delta0::Float64=1e-4,
    #
    TolS::Float64=1e-6, AMA_MaxIter::Int=200, AMA_SizeOrder::Int=FesOrder, AMA_AnisoOrder::Int=2, 
    #
    TolT::Float64=1e-3, Deltat0::Float64=1e-6,
    #
    PlotFig::Bool=false, wFig::Float64=4.50, hFig::Float64=6.50, 
    PlotVars::Vector{String}=String[], 
    SaveFig::Bool=false, Nt_SaveFig::Int=5, Deltat_SaveFig::Float64=Inf,
    mFig::Int=max(1,length(PlotVars)), nFig::Int=Int(ceil(length(PlotVars)/mFig)), cmap::String="jet",
    #
    SaveRes::Bool=false, Nt_SaveRes::Int=typemax(Int), Deltat_SaveRes::Float64=Inf, 
    #
    CSS::Float64=0.1, 
    #
    SC::Int=0)

    #---------------------------------------------------------------------
    #PROBLEM DATA:
    
    #Problem data:
    gamma           = 5/3
    
    #Define gas:
    GasModel                = GasIdeal()
    GasModel.gamma          = gamma
    GasModel.epsilon        = epsilon
    GasModel.nu             = nu
    GasModel.beta           = beta
    GasModel.kappa_rho_cv   = kappa_rho_cv
    GasModel.CSS            = CSS
    GasModel.g              = 1.0

    #Initial conditions:
    function u0fun(x::Vector{Matrix{Float64}})
        
        rho         = @. SmoothHeaviside(x[2], delta0, 1.0, 2.0)
        rhovx       = @mlv 0.0*x[1] 
        rhovy       = @mlv 0.03*(exp(-2*pi*x[2]^2)-exp(-pi/2))*cos(8*pi*x[1])
        p           = @. (2.0-2.0*x[2])*(x[2]>=0) + (2.0-1.0*x[2])*(x[2]<0)
        rhoE        = @mlv p/(gamma-1.0) + 0.5*rhovy^2/rho
        
        return [rho, rhovx, rhovy, rhoE]
        
    end
    
    #Boundary conditions:
    BC_walls        = SlipAdiabatic()
    
    #---------------------------------------------------------------------
    #PRE-PROCESS STAGE:
    
    #Mesh:
    MeshFile                = "../temp/RayleighTaylor_SC$(SC).geo"
    NX                      = Int(ceil(0.25/(hp0*FesOrder)))
    NY                      = Int(ceil(1.0/(hp0*FesOrder)))
    TrMesh_Rectangle_Create!(MeshFile, -0.125, 0.125, NX, -0.5, 0.5, NY)
    
    #Load LIRKHyp solver structure with default data. Modify the default data if necessary:
    solver                  = LIRKHyp_Start(GasModel)
    solver.ProblemName      = "RayleighTaylor"
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
    
    #Set initial and boundary conditions:
    solver.u0fun        = FW11((x) -> u0fun(x)) 
    solver.BC           = [BCW(BC_walls), BCW(BC_walls), BCW(BC_walls), BCW(BC_walls)]
    
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
                v_plot  = PlotContour(solver, solver.model, PlotVars[ii])
                title(latexstring(LatexString(PlotVars[ii]),
                    "; t^n=", sprintf1("%.2e", solver.t)), 
                    fontsize=10)
                println(PlotVars[ii], ", min=", minimum(v_plot), ", max=", maximum(v_plot))
            end
            if SaveFig
                savefig("$(VideosUbi)RayleighTaylor_SC$(SC)_$(nb_SaveFig).png", dpi=400, pad_inches=0)
            end
        
            figure(figv[2].number)
            PyPlot.cla()
            PlotContour(solver.u[1], solver.fes)
            PlotMesh!(solver.mesh, color="w")
            if SaveFig
                savefig("$(VideosUbi)RayleighTaylor_Mesh_SC$(SC)_$(nb_SaveFig).png", dpi=400, pad_inches=0)
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
                savefig("$(VideosUbi)RayleighTaylor_Errors$(SC)_$(nb_SaveFig).png", dpi=400, pad_inches=0)
            end
            
            t_lastFig           += Deltat_SaveFig
            ct_SaveFig          = 0
            nb_SaveFig          += 1
            
        end
        return
        
    end
    PlotSol()
    
    #Function to save intermediate results:
    t_lastRes           = 0.0
    ct_SaveRes          = 0
    nb_SaveRes          = 0
    function SaveSol()
        
        ct_SaveRes      += 1 
        if SaveRes && ( solver.t-t_lastRes>=Deltat_SaveRes || 
                        ct_SaveRes==Nt_SaveRes || solver.t==tf || solver.t==0.0 )
            save("$(ResUbi)LIRKHyp_SC$(SC)_$(nb_SaveRes).jld2", "StudyCase", "RayleighTaylor",
                "ConvFlag", ConvFlag, "solver", save(solver), 
                "epsilon", epsilon, "nu", nu, "beta", beta, "kappa_rho_cv", kappa_rho_cv,
                "delta0", delta0, "TolS", TolS, "TolT", TolT)
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
    
    #Save results:
    if SaveRes
        save("$(ResUbi)LIRKHyp_SC$(SC)_1000.jld2", "StudyCase", "RayleighTaylor",
            "ConvFlag", ConvFlag, "solver", save(solver) )
    end
    
    return solver
    
end
