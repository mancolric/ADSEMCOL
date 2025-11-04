include("test_CompressibleFlow.jl")

function TriplePoint(hp0::Float64, FesOrder::Int;
    tf::Float64=5.0, RKMethod::String="Ascher3",
    epsilon::Float64=1e-3, nu::Float64=0e-3, beta::Float64=0.0, kappa_rho_cv::Float64=0e-6, 
    delta0::Float64=epsilon,
    #
    TolS::Float64=1e-6, TolS0=0.1*TolS, AMA_MaxIter::Int=200, AMA_SizeOrder::Int=FesOrder, AMA_AnisoOrder::Int=2, 
    PostAdapt::Bool=false, 
    #
    TolT::Float64=1e-4, Deltat0::Float64=1e-5,
    #
    PlotFig::Bool=true, wFig::Float64=9.50, hFig::Float64=6.50, 
    PlotVars::Vector{String}=String["rho", "v", "e"], 
    SaveFig::Bool=false, Nt_SaveFig::Int=5, Deltat_SaveFig::Float64=Inf,
    mFig::Int=max(1,length(PlotVars)), nFig::Int=Int(ceil(length(PlotVars)/mFig)), cmap::String="jet",
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
        
        rho         = @. SmoothHeaviside(x[1], delta0, 1.0, 0.0)*1.0 + 
                            SmoothHeaviside(x[1], delta0, 0.0, 1.0)*
                                SmoothHeaviside(x[2], delta0, 1.0, 0.125)
        rhovx       = @mlv 0.0*x[1] 
        rhovy       = @mlv 0.0*x[1] 
        p           = @. SmoothHeaviside(x[1], delta0, 1.0, 0.1)
        rhoE        = @mlv p/(gamma-1.0)
        
        return [rho, rhovx, rhovy, rhoE]
        
    end
    
    #Boundary conditions:
    BC_walls        = SlipAdiabatic()
    
    #---------------------------------------------------------------------
    #PRE-PROCESS STAGE:
    
    #Mesh:
    MeshFile                = "$(@__DIR__)/../../temp/TriplePoint_SC$(SC).geo"
    NX                      = Int(ceil(7.0/(hp0*FesOrder)))
    NY                      = Int(ceil(3.0/(hp0*FesOrder)))
    TrMesh_Rectangle_Create!(MeshFile, -1.0, 6.0, NX, -1.5, 1.5, NY)
    
    #Load LIRKHyp solver structure with default data. Modify the default data if necessary:
    solver                  = LIRKHyp_Start(GasModel)
    solver.ProblemName      = "TriplePoint"
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
    solver.TolS_max         = TolS0
    solver.TolS_min         = 0.0*TolS
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
                println(PlotVars[ii], ": min=", minimum(v_plot), ", max=", maximum(v_plot))
            end
            if SaveFig
                savefig("$(VideosUbi)TriplePoint_SC$(SC)_$(nb_SaveFig).png", dpi=400, pad_inches=0)
            end
        
            figure(figv[2].number)
            PyPlot.cla()
            PlotContour(solver.u[1], solver.fes)
            PlotMesh!(solver.mesh, color="w")
            title(latexstring(LatexString("rho"),
                    "; t^n=", sprintf1("%.2e", solver.t)), 
                    fontsize=10)
            if SaveFig
                savefig("$(VideosUbi)TriplePoint_Mesh_SC$(SC)_$(nb_SaveFig).png", dpi=400, pad_inches=0)
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
                savefig("$(VideosUbi)TriplePoint_Errors$(SC)_$(nb_SaveFig).png", dpi=400, pad_inches=0)
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
            save("$(ResUbi)LIRKHyp_SC$(SC)_$(nb_SaveRes).jld2", "StudyCase", "TriplePoint",
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
    
    #Set new tolerance:
    solver.TolS_max     = TolS
    
    while solver.t<tf
    
#         ConvFlag            = LIRKHyp_Step!(solver)
        if PostAdapt
            ConvFlag            = LIRKHyp_Step_Post!(solver)
        else
#             ConvFlag            = LIRKHyp_Step_Pre!(solver)
            ConvFlag            = IRK_Step!(solver)
        end
        if ConvFlag<=0
            break
        end
        
        PlotSol()
        SaveSol()
        
    end
    
    #Save results:
    if SaveRes
        save("$(ResUbi)LIRKHyp_SC$(SC)_1000.jld2", "StudyCase", "TriplePoint",
            "ConvFlag", ConvFlag, "solver", save(solver) )
    end
    
    return solver
    
end
