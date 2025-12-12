include("test_ShallowWater.jl")

function CProperty2(hp0::Float64, FesOrder::Int;
    tf::Float64=0.48, RKMethod::String="BPR3",
    epsilon::Float64=0e-3, delta::Float64=1e-2, Deltah::Float64=0e-2,
    eta0::Float64=1e0, g::Float64=1.0, Deltaeta0::Float64=1e-2,  
    #
    TolS::Float64=1e-4, TolS0=0.01*TolS, AMA_MaxIter::Int=200, AMA_SizeOrder::Int=FesOrder, AMA_AnisoOrder::Int=2,
    #
    TolT::Float64=1e-4, Deltat0::Float64=1e-5, TimeAdapt::Bool=true,
    #
    CA::Float64=1e-2, 
    #
    PlotFig::Bool=false, wFig::Float64=9.50, hFig::Float64=6.50,
    PlotVars::Vector{String}=String[],
    SaveFig::Bool=false, Nt_SaveFig::Int=5, Deltat_SaveFig::Float64=Inf,
    mFig::Int=max(1,length(PlotVars)), nFig::Int=Int(ceil(length(PlotVars)/mFig)), cmap::String="jet",
    #
    SaveRes::Bool=false, Nt_SaveRes::Int=typemax(Int), Deltat_SaveRes::Float64=0.1,
    #
    CSS::Float64=0.1,
    #
    SC::Int=0)

    #---------------------------------------------------------------------
    #PROBLEM DATA:

    ##Define model:
    model               = SWE()
    model.epsilon       = epsilon
    model.g             = g
    model.CSS           = CSS
    av                  = [0.56, 0.71, 0*0.13, 0*0.29]
    model.b             = FW1( (x)-> @. av[1]*sin(2*pi*x[1])+av[2]*cos(4*pi*x[1]) + 
                                        av[3]*sin(3*pi*x[2])+av[4]*cos(5*pi*x[1]) - eta0 ) 
    function u0fun(x::Vector{Matrix{Float64}})

        b               = model.b(x)
        eta             = @tturbo @. eta0 + Deltaeta0*exp(-(x[1]/0.1)^2)
        q1              = @tturbo @. 0.0*x[1]
        q2              = @tturbo @. 0.0*x[1]
        
        return [eta, q1, q2, b]

    end

    #---------------------------------------------------------------------
    #PRE-PROCESS STAGE:

    #Mesh:
    MeshFile                = "$(@__DIR__)/../../temp/CProperty2$(SC).geo"
    NX                      = Int(ceil(2.0/(hp0*FesOrder)))
    NY                      = Int(ceil(1.0/(hp0*FesOrder)))
    TrMesh_Rectangle_Create!(MeshFile, -1.0, 1.0, NX, -0.5, 0.5, NY)

    #Load LIRKHyp solver structure with default data. Modify the default data if necessary:
    solver                  = LIRKHyp_Start(model)
    solver.ProblemName      = "CProperty2"
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
    solver.TolS_min         = 0.01*TolS0
    solver.TolT             = TolT
    solver.TimeAdapt        = TimeAdapt
    solver.CA               = CA

    #Set initial and boundary conditions:
    solver.u0fun        = FW11((x) -> u0fun(x))
    BC_walls            = SlipAdiabatic()
    function uInOut(t::Float64, x::Vector{Matrix{Float64}})
        eta     = @. eta0 + 0.0*x[1]
        q1      = @. 0.0*x[1]
        q2      = @. 0.0*x[1]
        return [eta, q1, q2]
    end
    BC_outer            = SubsonicInOutlet1(FWt11((t,x)->uInOut(t,x)))
#     BC_outer            = SupersonicOutlet1()
    solver.BC           = [BCW(BC_walls), BCW(BC_walls), BCW(BC_walls), BCW(BC_outer)]

    #-----------------------------------------------------------------------------
    #INITIAL CONDITION:

    #Compute initial condition:
    ConvFlag            = LIRKHyp_InitialCondition!(solver)
#     CheckJacobian(solver, Plot_df_du=false, Plot_df_dgradu=false, 
#         Plot_dQ_du=true, Plot_dQ_dgradu=false)
#     for ii = 2
#         BC_CheckJacobian(solver, ii, Plot_df_du=true)
#     end
#     return
    
    #Change TolT:
    if TolT==0.0
        TolT            = 0.01*solver.etaS
        solver.TolT     = TolT
    end
    
    #Function to plot solution:
    figv                = Vector{Figure}(undef,2)
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
                PlotMesh!(solver.mesh, color="k")
                title(latexstring(LatexString(model, PlotVars[ii]),
                    "; t^n=", sprintf1("%.2e", solver.t)),
                    fontsize=10)
                println(PlotVars[ii], ": min=", minimum(v_plot), ", max=", maximum(v_plot))
            end
            if SaveFig
                savefig("$(VideosUbiTFG)CProperty2$(SC)_$(nb_SaveFig).png", dpi=400, pad_inches=0)
            end
            figure(figv[2].number)
            #Loop plot variables:
            for ii=1:length(PlotVars)
                PyPlot.subplot(mFig, nFig, ii)
                PyPlot.subplots_adjust(hspace=0.8)
                PyPlot.cla()
                if PlotVars[ii]=="eta"
                    plot(solver.fes.PSpace.NodesCoords[:,1], 
                        solver.u[1][1:solver.fes.PSpace.nDof].-eta0, ".b", markersize=0.5)
                    title(latexstring("\\eta-\\eta_0, t^n=", 
                                        sprintf1("%.2e", solver.t)),
                        fontsize=10)
                else
                    splot_fun(x1,x2)    = @mlv x1
                    PlotNodes(splot_fun, solver, PlotVars[ii])
                    title(latexstring(LatexString(model, PlotVars[ii]),
                                    "; t^n=", sprintf1("%.2e", solver.t)),
                            fontsize=10)
                end
                xlabel(latexstring("x_1"), fontsize=10)
            end
            if SaveFig
                savefig("$(VideosUbiTFG)CProperty2Nodes$(SC)_$(nb_SaveFig).png", dpi=400, pad_inches=0)
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
            save("$(ResUbi)LIRKHyp_SC$(SC)_$(nb_SaveRes).jld2", "StudyCase", "CProperty2",
                "ConvFlag", ConvFlag, "solver", save(solver),
                 "TolS", TolS, "TolT", TolT)
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
    solver.TolS_min     = 0.01*TolS
    
    while solver.t<tf
    
        #Allocate equilibrium solution:
        u_eq                = 0.0*solver.uv
        u_eq_views          = GetViews(u_eq, solver.nVars, solver.fes.nDof)
        
        #Interpolate unit function:
        u1_terp             = zeros(solver.fes.nDof)
        u1_terp[1:solver.fes.PSpace.nDof]  .= 1.0
        
        #Compute approximation to eta_average:
#         eta_approx          = dot(u1_terp, solver.MII, solver.u[1])/
#                                 dot(u1_terp, solver.MII, u1_terp)
#         eta_approx          = 1.0
        eta_approx          = eta0
        
        #DEBUG:
        b1_proj2            = zeros(solver.fes.nDof)
        Integ2D             = TrInt(solver.mesh, solver.AMA_ProjOrder, N=solver.AMA_ProjN)
        L2Projection!(Integ2D, FW11((x) -> [ 0.0*x[1].+eta_approx ]), 
                        solver.fes, [b1_proj2])
        b1_proj2            *= 1.0
        
        #Project eta_approx onto fes:
        b1_proj             = zeros(solver.fes.nDof)
        Integ2D             = TrInt(solver.mesh, solver.AMA_ProjOrder, N=solver.AMA_ProjN)
        L2Projection!(Integ2D, FW11((x) -> [ 0.0*x[1].+1.0 ]), 
                        solver.fes, [b1_proj])
        b1_proj             *= eta_approx
        u1_proj             = zeros(solver.fes.nDof)
        Interpolate!(FW11((x) -> [ 0.0*x[1].+eta_approx ]), [ u1_proj ] , solver.fes)
        flag, nIter, etaA_I = LS_gmres!(solver.MII_LS, u1_proj, b1_proj, 
                                AbsTol=solver.TolA_min*eta_approx, Display="notify")
        
        #Correction factor corr_fact such that eta_eq = c_fact * eta_approx and 
        #mass of eta_eq = mass of solution
        corr_fact           = dot(u1_terp, solver.MII, solver.u[1])/
                                dot(u1_terp, solver.MII, u1_proj) 
        println("corr_fact=", corr_fact, ", eta_eq=", eta_approx*corr_fact)
        
        #Correct eta and b:
        u_eq_views[1]       .= corr_fact*u1_proj
        u_eq_views[4]       .= solver.u[4]
                
#         figure()
#         plot(u_eq_views[1]-solver.u[1], "g")
#         plot(u_eq_views[1][1:solver.fes.PSpace.nDof]-solver.u[1][1:solver.fes.PSpace.nDof], "b")
#         plot(u_eq_views[1][1:solver.fes.PSpace.nDofU]-solver.u[1][1:solver.fes.PSpace.nDofU], "g")
#         error("")
        
        #March in time:
        ConvFlag    = LIRKHyp_Step!(solver)
#         ConvFlag    = WB_LIRKHyp_Step!(solver, u_eq)
        
        if ConvFlag<=0
            break
        end

        PlotSol()
        SaveSol()

    end

    hmean               = 2.0*sqrt(solver.Omega/solver.mesh.nElems/TrElem_Area)
    Deltat_mean         = solver.tf/(solver.Nt-1)
    println("hmean=", sprintf1("%.2e", hmean), ", Deltat_mean=", sprintf1("%.2e", Deltat_mean))

    #Save results:
    if SaveRes
        save("$(ResUbi)LIRKHyp_SC$(SC)_1000.jld2", "StudyCase", "CProperty2",
            "ConvFlag", ConvFlag, "solver", save(solver), "hmean", hmean, "Deltat_mean", Deltat_mean)
    end

    return solver

end
