cd(@__DIR__)
include("test_CompressibleFlow.jl")

function ShockFlame(FesOrder::Int;
    tf::Float64=20.0, RKMethod::String="Ascher3",
    epsilon::Float64=1e-2, nu::Float64=0e-3, kappa_rho_cv::Float64=0e-6, 
    Pr::Float64=0.73, Le::Float64=1.0, YF0::Float64=0.02, 
    delta_flame::Float64=1.0, delta_shock::Float64=max(epsilon,nu),
    beta::Float64=10.0, q::Float64=5.0, ML::Float64=2.0, 
    Lx::Float64=120.0, 
    #
    TolS::Float64=1e-6, AMA_MaxIter::Int=200, AMA_SizeOrder::Int=FesOrder, AMA_AnisoOrder::Int=2, 
    #
    TolT::Float64=1e-5, Deltat0::Float64=1e-6,
    #
    PlotFig::Bool=false, wFig::Float64=9.50, hFig::Float64=6.50, 
    PlotVars::Vector{String}=String[], PlotCode::Vector{String}=fill("nodes", length(PlotVars)), 
    SaveFig::Bool=false, Nt_SaveFig::Int=5, Deltat_SaveFig::Float64=Inf,
    mFig::Int=max(1,length(PlotCode)), nFig::Int=Int(ceil(length(PlotCode)/mFig)), cmap::String="jet",
    #
    SaveRes::Bool=false, Nt_SaveRes::Int=typemax(Int), Deltat_SaveRes::Float64=0.1, 
    #
    CSS::Float64=0.1, 
    #
    SC::Int=0)

    #---------------------------------------------------------------------
    #PROBLEM DATA:
    
    #Fuel state:
    gamma           = 1.4
    p0              = 1.0
    RT0             = 1.0
    rho0            = p0/RT0
    Q0              = q*gamma/(gamma-1)*RT0/YF0
    
    #Burned gas state:
    RTL             = (1.0+q)*RT0
    pL              = p0
    rhoL            = pL/RTL
    YFL             = 0.0
    
    #Conditions at the right:
    uL              = ML*sqrt(gamma*RTL)
    rhoR            = rhoL / ( (gamma-1.0)/(gamma+1.0) + 2.0/(gamma+1.0)/ML^2 )
    uR              = rhoL*uL/rhoR
    pR              = pL+rhoL*uL^2-rhoR*uR^2
    YFR             = 0.0
    u0              = uL
    RTR             = pR/rhoR
    
    #Chapman--Jouget detonation:
    CJdet           = ChapmanJouget(gamma,q)
    rhoprime        = CJdet[5]#/CJdet[2]*rho0
    uprime          = CJdet[6]#/CJdet[3]*u0
    pprime          = CJdet[7]#/CJdet[4]*p0
    RTprime         = pprime/rhoprime
    betaprime       = beta*RT0/RTprime
    
    #Define gas:
    GasModel                = GasFP()
    GasModel.epsilon        = epsilon
    GasModel.gamma          = gamma
    GasModel.nu             = nu
    GasModel.beta           = 0.0
    GasModel.kappa_rho_cv   = nu/Pr * GasModel.gamma
    GasModel.D              = nu/Pr/Le
    GasModel.B              = (uprime * exp(betaprime))/(1.0*betaprime)
    GasModel.RTa            = beta*RT0
    GasModel.hfF            = Q0
    GasModel.hfP            = 0.0

    #Define domain:
    x2                  = Lx
    y2                  = 1.0*delta_flame
    x0_flame            = 10.0*delta_flame
    x0_shock            = 20.0*delta_flame
    MeshFile            = "../temp/ShockFlame_SC$(SC).geo"
    TrMesh_Rectangle_Create!(MeshFile, 0.0, x2, 5, 0.0, y2, 2)
    
    #Initial conditions:
    function u0fun(x::Vector{Matrix{Float64}})
        
        rhoYF           = @mlv SmoothHeaviside(x[1]-x0_flame, delta_flame, rho0*YF0, 0.0)
        rhoYP           = @mlv SmoothHeaviside(x[1]-x0_flame, delta_flame, rho0*(1.0-YF0), rhoL) + 
                                SmoothHeaviside(x[1]-x0_shock, delta_shock, 0.0, rhoR-rhoL)
        rhovx           = @mlv SmoothHeaviside(x[1]-x0_flame, delta_flame, rho0*u0, rhoL*uL) + 
                                SmoothHeaviside(x[1]-x0_shock, delta_shock, 0.0, rhoR*uR-rhoL*uL)
        rhovy           = @mlv 0.0*x[1]
        rhoE            = @mlv (SmoothHeaviside(x[1]-x0_flame, delta_flame, p0, pL) + 
                                SmoothHeaviside(x[1]-x0_shock, delta_shock, 0.0, pR-pL)) / (gamma-1.0) + 
                                0.5*rhovx^2/(rhoYF+rhoYP) + rhoYF*Q0
        return [rhoYF, rhoYP, rhovx, rhovy, rhoE]
        
    end

    #Boundary conditions:
    BC_walls        = SlipAdiabatic()
    function pOutlet(t::Float64, x::Vector{Matrix{Float64}})
        
        p           = @mlv 0.0*x[1]+pR
        return [p]
        
    end
    function uInlet(t::Float64, x::Vector{Matrix{Float64}})
        
        rhoYF       = @mlv 0.0*x[1] + rho0*YF0
        rhoYP       = @mlv 0.0*x[1] + rho0*(1.0-YF0)
        rhovx       = @mlv 0.0*x[1] + rho0*u0
        rhovy       = @mlv 0.0*x[1]
        rhoE        = @mlv p0/(gamma-1.0) + 0.5*rho0*u0^2 + rhoYF*GasModel.hfF + rhoYP*GasModel.hfP
        return [rhoYF, rhoYP, rhovx, rhovy, rhoE]
        
    end
#     BC_outlet       = SubsonicOutlet1(FWt11((t,x)->pOutlet(t,x)))
#     BC_outlet       = SupersonicOutlet1()
    BC_outlet       = DoNothing1()
    BC_inlet        = SupersonicInlet1(FWt11((t,x)->uInlet(t,x)))
    
    #---------------------------------------------------------------------
    #PRE-PROCESS STAGE:
    
    #Load LIRKHyp solver structure with default data. Modify the default data if necessary:
    solver                  = LIRKHyp_Start(GasModel)
    solver.ProblemName      = "ShockFlame"
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
                savefig("$(VideosUbi)ShockFlame_SC$(SC)_$(nb_SaveFig).png", dpi=400, pad_inches=0)
            end
        
            figure(figv[2].number)
            PyPlot.cla()
            plot(solver.mesh)
            if SaveFig
                savefig("$(VideosUbi)ShockFlame_Mesh$(SC)_$(nb_SaveFig).png", dpi=400, pad_inches=0)
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
                savefig("$(VideosUbi)ShockFlame_Errors$(SC)_$(nb_SaveFig).png", dpi=400, pad_inches=0)
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
            save("$(ResUbi)LIRKHyp_SC$(SC)_$(nb_SaveRes).jld2", "StudyCase", "ShockFlame",
                "ConvFlag", ConvFlag, "solver", save(solver), 
                "epsilon", epsilon, "nu", nu, "kappa_rho_cv", kappa_rho_cv,
                "Pr", Pr, "Le", Le, "YF0", YF0,
                "delta_flame", delta_flame, "delta_shock", delta_shock, 
                "beta", beta, "q", q, "ML", ML,
                "D", GasModel.D, "B", GasModel.B, "RTa", GasModel.RTa, "Q0", GasModel.hfF, 
                "x2", x2, "y2", y2,
                "TolS", TolS, "TolT", TolT)
            t_lastRes   += Deltat_SaveRes
            ct_SaveRes  = 0
            nb_SaveRes  += 1
        end
        if SaveRes && solver.t==tf
#             save("$(ResUbi)../../ShockFlame/LIRKHyp_SC$(SC)_info.jld2", "nb_SaveRes", nb_SaveRes-1)
            save("$(ResUbi)LIRKHyp_SC$(SC)_info.jld2", "nb_SaveRes", nb_SaveRes-1)
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

function ChapmanJouget(gamma::Float64, q::Float64)

    #We assume the state at the end p2=1, rho2=1, u2=sqrt(gamma), YF2=0.
    p2      = 1.0
    rho2    = 1.0
    u2      = sqrt(gamma*p2/rho2)
    
    #Then, we compute the state "0" by imposing conservation equations:
    #   rho0*u0                                                 = M
    #   p0+rho0*u0^2                                            = P
    #   (gamma/(gamma-1)*p0 + 0.5*rho0*u0^2 + rho0*YF0*Q0)*u0   = HU
    #Note that rho0*Q0*YF0 = q*gamma/(gamma-1)*p0, so the last equation is
    #   ((1+q)*gamma/(gamma-1)*p0 + 0.5*rho0*u0^2)*u0           = HU
    #If we write the equations in terms of u0:
    #A u0^2 + B*u0 + HU = 0
    M   = rho2*u2
    P   = p2+rho2*u2^2
    HU  = (gamma/(gamma-1.0)*p2 + 0.5*rho2*u2^2)*u2
    #
    A   = ((1.0+q)*gamma/(gamma-1.0)-0.5)*M
    B   = -(1.0+q)*gamma/(gamma-1.0)*P
    C   = HU
    #
    u0      = (-B+sqrt(B^2-4*A*C))/(2*A)
    rho0    = M/u0
    p0      = P-M*u0
    M0      = u0/sqrt(gamma*p0/rho0)
    
    #Rescale solution with entry conditions:
    rho2    *= 1.0/rho0
    p2      *= 1.0/p0
    u2      = sqrt(gamma*p2/rho2)
    #
    rho0    = 1.0
    p0      = 1.0
    u0      = M0*sqrt(gamma)
    
    #=
    #Check results:
    display(rho0*u0 - rho2*u2)
    display(p0+rho0*u0^2 - (p2+rho2*u2^2))
    display(((1+q)*gamma/(gamma-1)*p0 + 0.5*rho0*u0^2)*u0 - 
            (gamma/(gamma-1)*p2 + 0.5*rho2*u2^2)*u2)
    =#
    
    #State at region 1:
    rho1            = rho0 / ( (gamma-1.0)/(gamma+1.0) + 2.0/(gamma+1.0)/M0^2 )
    u1              = rho0*u0/rho1
    p1              = p0+rho0*u0^2-rho1*u1^2
    
    return M0, rho0, u0, p0, rho1, u1, p1, rho2, u2, p2
    
end
