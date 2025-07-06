cd(@__DIR__)
include("test_CompressibleFlow.jl")

function AcousticWave(MeshCase::Int, FesOrder::Int, tfv::Array{Float64,1};
    RKMethod::String="Ascher3", 
    epsilon::Float64=0.0, nu::Float64=0e-6, beta::Float64=0.0, kappa_rho_cv::Float64=0e-6,
    PlotVars::Vector{String}=String[], PlotCode::Vector{String}=fill("nodes", length(PlotVars)), 
    SaveFig::Bool=false, wFig::Float64=7.50, hFig::Float64=7.50, 
    mFig::Int=max(1,length(PlotCode)), nFig::Int=Int(ceil(length(PlotCode)/mFig)), Nt_SaveFig::Int=1, cmap::String="jet",
    SC::Int=0, CSS::Float64=0.1, CDC::Float64=5.0, CFLa::Float64=1.0, 
    AMA_MaxIter::Int=20, TolS::Float64=1e-3, TolT::Float64=1e-3, AMA_ProjOrder::Int=2)

    #---------------------------------------------------------------------
    #PROBLEM DATA:
    
    #Problem data:
    gamma           = 7/5
    vortex_st       = 5.0   #vortex strength
    u_inf           = 1.0   #vortex horizontal velocity
    
    #Boundary conditions:
    BC_all          = SlipAdiabatic()
    
    #Initial conditions:
    function u0fun(t::Float64, x::Vector{Matrix{Float64}})

        lambda      = 0.25
        omega       = @mlv 0.01*(1-cos(2*pi*(x[1]-lambda)/lambda))*(lambda<=x[1]<=2*lambda)
        rho         = @mlv 1.0 + omega
        rhovx       = @mlv rho*omega
        rhovy       = @mlv 0.0*x[1]
        p           = @mlv 1.0 + omega
        rhoE        = @mlv p/(gamma-1.0) + 0.5*(rhovx*rhovx+rhovy*rhovy)/rho
        return [rho, rhovx, rhovy, rhoE]
        
    end
    
    #---------------------------------------------------------------------
    #PRE-PROCESS STAGE:
    
    #Define gas:
    GasModel                = GasIdeal()
    GasModel.gamma          = gamma
    GasModel.epsilon        = epsilon
    GasModel.nu             = nu
    GasModel.beta           = beta
    GasModel.kappa_rho_cv   = kappa_rho_cv
    GasModel.CSS            = CSS
    
    #Load LIRKHyp solver structure with default data. Modify the default data if necessary:
    solver              = LIRKHyp_Start(GasModel)
    solver.ProblemName  = "AcousticWave"
    solver.SC           = SC
    solver.MeshFile     = "$(MeshUbi)SodStandard/MeshCase$(MeshCase).db"
    solver.nBounds      = 4
    solver.FesOrder     = FesOrder
    solver.RKMethod     = RKMethod
    solver.CFL          = CFLa
    solver.tf           = tfv[end]
    solver.SpaceAdapt   = false
    solver.AMA_MaxIter  = AMA_MaxIter
    solver.AMA_ProjOrder= AMA_ProjOrder
    solver.TolS_max     = TolS
    solver.TolS_min     = 0.1*TolS
    solver.TolT         = TolT
    
    #Set initial and boundary conditions:
    solver.u0fun        = FW11((x) -> u0fun(0.0,x)) 
    solver.BC           = fill(BCW(BC_all),4)
    
    #-----------------------------------------------------------------------------
    #INITIAL CONDITION:
    
    #Compute initial condition:
    LIRKHyp_InitialCondition!(solver)
    
    #Compute Lq error:
#     errLq,              = LqError(solver, FW11((x) -> utheor(solver.t, x)), q=2.0)
    hmean               = 2.0*sqrt(solver.Omega/solver.mesh.nElems/TrElem_Area)
    
    #Open plots (variables, errors):
    figv                = Vector{Figure}(undef,1)
    if !isempty(PlotVars)
        figv[1]         = PyPlotSubPlots(mFig, nFig, w=wFig, h=hFig, left=0.9, right=0.4, bottom=1.1, top=0.5)
    end
    NPlotSave           = 0
    NSave               = 0
    function PlotSol()
        for ii=1:length(PlotVars)
            PyPlot.subplot(mFig, nFig, ii)
            PyPlot.cla()
            plotNodes((x1,x2)-> x1, solver, PlotVars[ii])
        end
    end
    PlotSol()
        
    #Plot solution:
#     figure()
#     plot3D(solver, "rho", color="b")
#     return
    
    #-----------------------------------------------------------------------------
    #MARCH IN TIME:
    
    while solver.t<tfv[end]
    
        flag    = LIRKHyp_Step!(solver)
        if flag<0
            break
        end
        
        PlotSol()
        
    end
    
    return solver
    
    #Boundary conditions:
    PD.BC_Type      = [ "SA", "SubO_p", "SA", "auto" ]
#     PD.BC_Type      = [ "SA", "auto", "SA", "auto" ]
    PD.BC_Fun       = Vector{Vector{Function}}(undef,4)
    PD.BC_Fun[2]    = [ (t,xm) -> fill(1.0, size(xm,1)) ]   #Output pressure=rho*RT=1
#     PD.BC_Fun[2]    = [ (t,xm) -> fill(1.0, size(xm,1)),    #Input rho
#                         (t,xm) -> fill(1.0, size(xm,1)),    #Input rhovn
#                         (t,xm) -> fill(0.0, size(xm,1)),    #Input rhovt
#                         (t,xm) -> fill(1.0/(gamma-1.0)+0.5, size(xm,1)) ]   #Input rhoE
    PD.BC_Fun[4]    = [ (t,xm) -> fill(1.0, size(xm,1)),    #Input rho
                        (t,xm) -> fill(-1.0, size(xm,1)),    #Input rhovn
                        (t,xm) -> fill(0.0, size(xm,1)),    #Input rhovt
                        (t,xm) -> fill(1.0/(gamma-1.0)+0.5, size(xm,1)) ]   #Input rhoE
    
    
    #Time step:
    PD.CFLa_Param       = CFLa
    PD.RKMethod         = RKMethod

    #Plotting options:
    xplotv          = linspace(-4.9999,14.9999,500)
    yplotv          = linspace(0.0,0.0,500)
#     splot           = (xm) -> @views @. sqrt( (xm[:,1]-5.0)^2 + xm[:,2]^2 ) 
    splot           = (xm) -> xm[:,1]

    
    #Allocate structures for errors:
#     rho_errs    = VarErrors()
#     rhovx_errs  = VarErrors()
#     rhovy_errs  = VarErrors()
#     rhoE_errs   = VarErrors()
    
    #---------------------------------------------------------------------
    #SOLVER STAGE:
    
    #Initial condition:
    solver          = CLG12_Start(PD, GasModel)
    
    #Plot initial condition:
    figv            = Vector{Figure}(undef,3)
    if !isempty(PlotCode)
        if any(PlotCode.=="nodes")
            figv[1] = PyPlotSubPlots(mFig, nFig, w=wFig, h=hFig, left=0.9, right=0.4, bottom=1.1, top=0.5)
        else
            figv[1] = PyPlotSubPlots(mFig, nFig, w=wFig, h=hFig, left=0.9, right=0.4, bottom=0.1, top=0.5)
        end
        for ii=2:length(figv)
            figv[ii]    = figure()
        end
    end
    NPlotSave       = 0
    NSave           = 0
    function PlotRes!()
        if !isempty(PlotCode)
            
            #Open subplots for fig1 and plot variables:
            figure(figv[1].number)
            PlotSolution!(GasModel, solver, PlotCode, PlotVars, mFig, nFig, splot)
            if SaveFig
                savefig("$(VideosUbi)SmoothVortex_SC$(SC)_$(NSave).png", dpi=400, pad_inches=0)
            end   
            
            #=
            bmesh           = solver.sol.bmesh[4]
            xm              = bmesh.NodesCoords
            ParentNodes     = bmesh.ParentNodes
            PyPlot.subplot(4,1,1)
            PyPlot.cla()
            plot(xm[:,2], solver.sol.rho.uv[1][ParentNodes].-1.0, ".b")
            PyPlot.subplot(4,1,2)
            PyPlot.cla()
            plot(xm[:,2], solver.sol.rhovx.uv[1][ParentNodes].-1.0, ".b")
            PyPlot.subplot(4,1,3)
            PyPlot.cla()
            plot(xm[:,2], solver.sol.rhovy.uv[1][ParentNodes].-0.0, ".b")
            PyPlot.subplot(4,1,4)
            PyPlot.cla()
            rhoe            = solver.sol.rhoE.uv[1][ParentNodes] .- 0.5
            plot(xm[:,2], solver.sol.rhoE.uv[1][ParentNodes].-(1.0/0.4+0.5), ".b")
#             plot(xm[:,2], rhoe.-(1.0/0.4), ".b")
            title("t=$(solver.sol.t)")
            =#
            
            figure(figv[2].number)
            Nt          = length(solver.tv)
            semilogy(solver.tv, solver.etav_S, ".-b")
            semilogy(solver.tv, fill(PD.TolS_max, Nt), "--b")
            semilogy(solver.tv, solver.etav_T, ".-g")
            semilogy(solver.tv, fill(PD.TolT, Nt), "--g")
            semilogy(solver.tv, solver.etav_A, ".-r")
            semilogy(solver.tv, max.(1e-8, PD.C_A_max.*min.(PD.TolS_max, PD.TolT, solver.etav_S, solver.etav_T)), "--r")
            
            figure(figv[3].number)
            semilogy(solver.tv, solver.nElemsv, ".-b")
            
            NSave       += 1
        end
    end
    PlotRes!()
    
    #Errors:
    tsimv       = [ 0.0 ]
    errv        = [ ComputeL2MeanError(solver.sol, utheor_fun) ]
    
    #Loop until final time:
    PD.tf       = tfv[length(PD.tf)]
    while solver.sol.t<PD.tf
        
        #Perform step:
        flag    = CLG12_Step!(solver)
        if flag<0
            break
        end
        
        #Append errors:
        tsimv       = vcat(tsimv, solver.sol.t)
        errv        = vcat(errv, ComputeL2MeanError(solver.sol, utheor_fun))
        println("<e>=",errv[length(errv)])
        
        #Plot:
        NPlotSave       += 1
        if NPlotSave==1
            NPlotSave   = 0
            PlotRes!()
        end
        
    end
#     PlotRes!()
    
    #Profile info:
    ProfileVars_display(solver.ProfVars)
    
    GC.gc()
    
    return solver
    
end

function test_loop(MCv::AbstractVector{Int}, orderv::AbstractVector{Int})

    Norder  = length(orderv)
    colorv  = PyPlotColors("jet2", Norder)
    
    NMC         = length(MCv)
    hpv         = zeros(NMC)
    errv        = zeros(NMC)
    errv_est    = zeros(NMC)
    
    figure()
    leg         = String[]
    for ii=1:Norder
        for jj=1:NMC
            hpv[jj], errv[jj], errv_est[jj]   = SmoothVortex(MCv[jj], orderv[ii], [0.0])
        end
        loglog(hpv*orderv[ii], errv, color=colorv[ii], linestyle="solid")
        loglog(hpv*orderv[ii], errv_est, color=colorv[ii], linestyle="dashed")
        leg     = vcat(leg, "p=$(orderv[ii])", "p=$(orderv[ii])")
        EOCv    = @. log(abs($diff(errv))) / log(abs($diff(hpv)))
        println(EOCv)
    end
    legend(leg)
    grid("on")
    
    return
    
end
