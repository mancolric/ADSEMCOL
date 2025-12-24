include("test_CompressibleFlow.jl")
include("InputData/ReactiveGasAuxFunctions.jl")

function ShockFlame_FXP(; FesOrder::Int=5, 
    tf::Float64=20.0, RKMethod::String="Ascher3",
    epsilon::Float64=1e-2, nu::Float64=0e-3, kappa_rho_cv::Float64=0e-6, 
    Pr::Float64=0.73, Le::Float64=1.0, YF0::Float64=1.0, 
    delta_flame::Float64=1.0, delta_shock::Float64=max(epsilon,nu),
    q::Float64=5.0, ML::Float64=2.0, 
    RThatI::Float64=10.0, RThatB::Float64=10.0, #RThat=RT_a/RT_1, #1: means post-shock condition
    BIhat::Float64=2.5e-3, BBhat::Float64=4.4e5, 
    Lx::Float64=120.0, 
    #
    TolS::Float64=1e-6, AMA_MaxIter::Int=200, AMA_SizeOrder::Int=FesOrder, AMA_AnisoOrder::Int=2, 
    #
    TolT::Float64=1e-5, Deltat0::Float64=1e-6,
    #
    PlotFig::Bool=false, MarkSF::Bool=true, wFig::Float64=9.50, hFig::Float64=6.50, 
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
    
    #Chapman--Jouget state with same incoming conditions:
    Min, rhoin, uin, pin, rhoCJ, uCJ, pCJ, _ = ChapmanJouget(gamma,q)
    rhoCJ           *= rho0/rhoin
    pCJ             *= p0/pin
    RTCJ            = pCJ/rhoCJ
    uCJ             *= sqrt(p0/rho0) / sqrt(pin/rhoin)
    
    #Define gas:
    GasModel                = GasFXP()
    GasModel.epsilon        = epsilon
    GasModel.gamma          = gamma
    GasModel.nu             = nu
    GasModel.beta           = 0.0
    GasModel.kappa_rho_cv   = nu/Pr * GasModel.gamma
    GasModel.hfF            = 0.0
    GasModel.hfX            = 0.0
    GasModel.hfP            = -Q0
    GasModel.D              = nu/Pr/Le
    GasModel.BR             = 1.0   #tc = BR^{-1} 
    GasModel.BI             = BIhat * GasModel.BR 
    GasModel.BB             = BBhat * GasModel.BR/rhoCJ
    GasModel.RTI            = RThatI * RTCJ
    GasModel.RTB            = RThatB * RTCJ
    
    #Define domain:
    x2                  = Lx
    y2                  = 1.0*delta_flame
    x0_flame            = 10.0*delta_flame
    x0_shock            = 20.0*delta_flame
    MeshFile            = "$(@__DIR__)/../../temp/SC$(SC).geo"
    TrMesh_Rectangle_Create!(MeshFile, 0.0, x2, 5, 0.0, y2, 2)
    
    #Initial conditions:
    function u0fun(x::Vector{Matrix{Float64}})
        
        rhoYF           = @mlv SmoothHeaviside(x[1]-x0_flame, delta_flame, rho0*YF0, 0.0)
        rhoYX           = @mlv SmoothHeaviside(x[1]-x0_flame, delta_flame, rho0*(1.0-YF0), 0.0)
        rhoYP           = @mlv SmoothHeaviside(x[1]-x0_flame, delta_flame, 0.0, rhoL) + 
                                SmoothHeaviside(x[1]-x0_shock, delta_shock, 0.0, rhoR-rhoL)
        rhovx           = @mlv SmoothHeaviside(x[1]-x0_flame, delta_flame, rho0*u0, rhoL*uL) + 
                                SmoothHeaviside(x[1]-x0_shock, delta_shock, 0.0, rhoR*uR-rhoL*uL)
        rhovy           = @mlv 0.0*x[1]
        rhoE            = @mlv (SmoothHeaviside(x[1]-x0_flame, delta_flame, p0, pL) + 
                                SmoothHeaviside(x[1]-x0_shock, delta_shock, 0.0, pR-pL)) / (gamma-1.0) + 
                                0.5*rhovx^2/(rhoYF+rhoYX+rhoYP) - rhoYP*Q0
        return [rhoYF, rhoYX, rhoYP, rhovx, rhovy, rhoE]
        
    end

    #Boundary conditions:
    BC_walls        = SlipAdiabatic()
    function pOutlet(t::Float64, x::Vector{Matrix{Float64}})
        
        p           = @mlv 0.0*x[1]+pR
        return [p]
        
    end
    function uInlet(t::Float64, x::Vector{Matrix{Float64}})
        
        return u0fun(x)
        
    end
#     BC_outlet       = SubsonicOutlet1(FWt11((t,x)->pOutlet(t,x)))
#     BC_outlet       = SupersonicOutlet1()
    BC_outlet       = DoNothing1()
    BC_inlet        = SupersonicInlet1(FWt11((t,x)->uInlet(t,x)))
    
    #---------------------------------------------------------------------
    #PRE-PROCESS STAGE:
    
    #Load LIRKHyp solver structure with default data. Modify the default data if necessary:
    solver                  = LIRKHyp_Start(GasModel)
    solver.ProblemName      = "ShockFlame_FXP"
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
#     CheckJacobian(solver, Plot_dQ_du=true)
#     return
    
    #Reference values of p and YF:
    pref, YFref         = pYref_ShockFlame(solver)
    
    #Location of shock and flame and relative Mach number:
    xs, xf              = xsxf_ShockFlame(solver, pref, YFref)
    tv                  = [ 0.0 ]
    xsv                 = [ xs ]
    xfv                 = [ xf ]
    Mrel                = Mrel_ShockFlame(solver, tv, xsv, xfv)
    Mrelv               = [ Mrel ]
    
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
        
            println("Plotting")
            
            figure(figv[1].number)
    
            #Loop plot variables:
            for ii=1:length(PlotVars)
            
                PyPlot.subplot(mFig, nFig, ii)
                PyPlot.cla()

                splot_fun(x1,x2)    = @mlv x1
                v_plot              = PlotNodes(splot_fun, solver, PlotVars[ii])
                xlabel(latexstring("x_1"), fontsize=10)
                title(latexstring(LatexString(GasModel, PlotVars[ii]),
                    "; t^n=", sprintf1("%.2e", solver.t)), 
                    fontsize=10)
                    
                if MarkSF
                    plot([xf, xf], [0.9*minimum(v_plot), 1.1*maximum(v_plot)], "--r")
                    plot([xs, xs], [0.9*minimum(v_plot), 1.1*maximum(v_plot)], "--g")
                end
                
            end
            if SaveFig
                savefig("$(VideosUbi)ShockFlame_FXP_SC$(SC)_$(nb_SaveFig).png", dpi=400, pad_inches=0)
            end
        
            figure(figv[2].number)
            PyPlot.cla()
            plot(solver.mesh)
            if SaveFig
                savefig("$(VideosUbi)ShockFlame_FXP_Mesh$(SC)_$(nb_SaveFig).png", dpi=400, pad_inches=0)
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
                savefig("$(VideosUbi)ShockFlame_FXP_Errors$(SC)_$(nb_SaveFig).png", dpi=400, pad_inches=0)
            end
            
            t_lastFig           += Deltat_SaveFig
            ct_SaveFig          = 0
            nb_SaveFig          += 1
            
            println("Plot finished")
            
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
            save("$(ResUbi)LIRKHyp_SC$(SC)_$(nb_SaveRes).jld2", "StudyCase", "ShockFlame_FXP",
                "ConvFlag", ConvFlag, "solver", save(solver), 
                "epsilon", epsilon, "nu", nu, "kappa_rho_cv", kappa_rho_cv, "gamma", gamma, 
                "Pr", Pr, "Le", Le, "YF0", YF0,
                "delta_flame", delta_flame, "delta_shock", delta_shock, 
                "q", q, "ML", ML,
                "D", GasModel.D, "BI", GasModel.BI, "BR", GasModel.BR, "BB", GasModel.BB, 
                "RTI", GasModel.RTI, "RTB", GasModel.RTB, "Q0", -GasModel.hfP, 
                "rhoCJ", rhoCJ, "pCJ", pCJ, "RTCJ", RTCJ, "uCJ", uCJ, 
                "x2", x2, "y2", y2,
                "TolS", TolS, "TolT", TolT, 
                "tv", tv, "xsv", xsv, "xfv", xfv, "Mrelv", Mrelv)
            t_lastRes   += Deltat_SaveRes
            ct_SaveRes  = 0
            nb_SaveRes  += 1
        end
        if SaveRes && solver.t==tf
#             save("$(ResUbi)../../ShockFlame_FXP/LIRKHyp_SC$(SC)_info.jld2", "nb_SaveRes", nb_SaveRes-1)
            save("$(ResUbi)LIRKHyp_SC$(SC)_info.jld2", "nb_SaveRes", nb_SaveRes-1)
        end
        return
        
    end
    SaveSol()
    
    #-----------------------------------------------------------------------------
    #MARCH IN TIME:
    
    while solver.t<tf
    
        ConvFlag    = LIRKHyp_Step!(solver)
#         ConvFlag    = IRK_Step!(solver)
        
        if ConvFlag<=0
            break
        end
        
        #Location of shock and flame, and relative Mach number to the shock speed:
        xs, xf              = xsxf_ShockFlame(solver, pref, YFref)
        push!(tv, solver.t)
        push!(xsv, xs)
        push!(xfv, xf)
        Mrel                = Mrel_ShockFlame(solver, tv, xsv, xfv)
        push!(Mrelv, Mrel)
        
        PlotSol()
        SaveSol()
        
    end
    
    return solver
    
end

#--------------------------------------------------------------------------

function pYref_ShockFlame(solver::SolverData)

    GasModel    = solver.model
    
    pv          = SolutionAtNodes(solver,GasModel,"p")
    YFv         = SolutionAtNodes(solver,GasModel,"Y_F")
    
    pref        = 0.5*(maximum(pv)+minimum(pv))
    YFref       = 0.5*(maximum(YFv)+minimum(YFv))
   
    return pref, YFref
    
end

function xsxf_ShockFlame(solver::SolverData, pref::Float64, YFref::Float64)

    GasModel    = solver.model
    
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

function Mrel_ShockFlame(solver::SolverData, tv::Vector{Float64}, 
    xsv::Vector{Float64}, xfv::Vector{Float64})

    Nt          = length(tv)
    
    #Check if this is not the first time level:
    if Nt==1
        return NaN
    end
    
    #Check that flame has passed the shock:
    if xfv[Nt]<xsv[Nt]
        return NaN
    end
    
    #Compute shock speed:
    if Nt>2
        xsdot   = (3*xsv[Nt]-4*xsv[Nt-1]+xsv[Nt-2])/(3*tv[Nt]-4*tv[Nt-1]+tv[Nt-2])
    else
        xsdot   = (xsv[Nt]-xsv[Nt-1])/(tv[Nt]-tv[Nt-1])
    end
    
    #Compute relative Mach number at mesh nodes:
    GasModel    = solver.model
    vx          = SolutionAtNodes(solver, GasModel, "vx")
    vy          = SolutionAtNodes(solver, GasModel, "vy")
    a           = SolutionAtNodes(solver, GasModel, "a")
    Mrel        = @tturbo @. sqrt((vx-xsdot)^2 + vy^2)/a
    
    #Get Mrel at left boundary:
    x1v         = solver.mesh.NodesCoords[:,1]
    sortv       = sortperm(x1v)
    x1v         = x1v[sortv]
    Mrel        = Mrel[sortv]
    aux         = findfirst(x1v.>0.0)
    
    return Mrel[aux]
    
end
