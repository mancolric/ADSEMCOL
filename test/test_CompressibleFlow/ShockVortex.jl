include("test_CompressibleFlow.jl")

function ShockVortex(MeshCase::Int, FesOrder::Int, tfv::Array{Float64,1};
    RKMethod::String="Ascher3",
    epsilon::Float64=0.0, nu::Float64=0e-6, beta::Float64=0.0, kappa_rho_cv::Float64=0e-6, 
    delta0::Float64=1e-4,
    PlotVars::Vector{String}=String[], PlotCode::Vector{String}=fill("nodes", length(PlotVars)), 
    SaveFig::Bool=false, wFig::Float64=9.50, hFig::Float64=6.50, 
    Nt_SaveFig::Int=typemax(Int), cmap::String="jet",
    SC::Int=0, CSS::Float64=0.1, CDC::Float64=5.0, pInf::Float64=1e-9, Deltat0::Float64=1e-5,
    AMA_MaxIter::Int=40, TolS::Float64=1e-6, TolT::Float64=1e-3, 
    AMA_SizeOrder::Int=3, AMA_AnisoOrder::Int=2,
    SaveRes::Bool=false)

    #---------------------------------------------------------------------
    #PROBLEM DATA:
    
    #Problem data:
    gamma           = 1.4
    xdisc           = 0.5
    ML              = 1.5
    rhoL            = 1.0
    pL              = 1.0
    uL              = ML*sqrt(gamma*pL/rhoL)
    rhoEL           = pL/(gamma-1.0)+0.5*rhoL*uL^2
    rhoR            = rhoL / ( (gamma-1.0)/(gamma+1.0) + 2.0/(gamma+1.0)/ML^2 )
    uR              = rhoL*uL/rhoR
    pR              = pL+rhoL*uL^2-rhoR*uR^2
    rhoER           = pR/(gamma-1.0)+0.5*rhoR*uR^2
    
    #Vortex:
    e0              = pL/rhoL/(gamma-1.0)
    a               = 0.075
    b               = 0.175
    x0              = 0.25
    y0              = 0.5
    Mv              = 0.9
    vm              = sqrt(gamma*pL/rhoL)*Mv
        
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
        
        #Vortex conditions:
        #
        #Radius:
        rv          = @mlv sqrt((x[1]-x0)^2 + (x[2]-y0)^2)
        #Internal energy:
        ev          = zeros(size(rv))
        aux         = findall(b.<=rv)
        @mlv ev[aux]= e0
        aux         = findall(a.<=rv.<b)
        @mlv ev[aux]= e0 - 1.0/gamma*vm^2*a^2/(a^2-b^2)^2 * (0.5*(b^2-rv[aux]^2)-2*b^2*log(b/rv[aux])-0.5*b^4*(1.0/b^2-1.0/rv[aux]^2))
        aux         = findall(rv.<a)
        @mlv ev[aux]= e0 - 1.0/gamma*vm^2*a^2/(a^2-b^2)^2 * (0.5*(b^2-a^2)-2*b^2*log(b/a)-0.5*b^4*(1.0/b^2-1.0/a^2)) - 
                        1.0/gamma*vm^2 * 0.5*(1.0-(rv[aux]/a)^2)
        #
        #Velocity:
        vv          = zeros(size(rv))
        @mlv vv     = vm*rv/a*(rv<=a) + vm*a/(a^2-b^2)*(rv-b^2/rv)*(a<rv<=b)
        #
        rhovortex   = @mlv rhoL*(ev/e0)^(1.0/(gamma-1.0))
        pvortex     = @mlv pL*(ev/e0)^(gamma/(gamma-1.0))
        
        #Initial condition:
        rho         = @mlv SmoothHeaviside(x[1]-xdisc, delta0, rhoL, rhoR) + 
                            (rhovortex-rhoL)
        rhovx       = @mlv rho*(SmoothHeaviside(x[1]-xdisc, delta0, uL, uR) - 
                            vv*(x[2]-y0)/(rv+1e-14))
        rhovy       = @mlv rho*vv*(x[1]-x0)/(rv+1e-14)
        rhoE        = @mlv (SmoothHeaviside(x[1]-xdisc, delta0, pL, pR) + 
                            pvortex-pL)/(gamma-1.0) + 
                             + 0.5*(rhovx*rhovx+rhovy*rhovy)/rho
#         rhoE        = @mlv rho*ev
        
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
        rho         = @mlv 0.0*x[1] + rhoL
        rhovx       = @mlv 0.0*x[1] + rhoL*uL
        rhovy       = @mlv 0.0*x[1]
        rhoE        = @mlv 0.0*x[1] + pL/(gamma-1.0) + 0.5*rhoL*uL^2
        return [rho, rhovx, rhovy, rhoE]
    end
    BC_inlet        = SupersonicInlet1(FWt11((t,x)->uInlet(t,x)))
    
    #---------------------------------------------------------------------
    #PRE-PROCESS STAGE:
    
    #Load LIRKHyp solver structure with default data. Modify the default data if necessary:
    solver                  = LIRKHyp_Start(GasModel)
    solver.ProblemName      = "ShockVortex, MeshCase=$MeshCase"
    solver.SC               = SC
    solver.MeshFile         = "$(MeshUbi)ShockVortex/MeshCase$(MeshCase).geo"
    solver.nBounds          = 4
    solver.FesOrder         = FesOrder
    solver.RKMethod         = RKMethod
    solver.Deltat0          = Deltat0
    solver.tf               = tfv[length(tfv)]
    solver.AMA_MaxIter      = AMA_MaxIter
    solver.AMA_SizeOrder    = AMA_SizeOrder
    solver.AMA_AnisoOrder   = AMA_AnisoOrder
    solver.TolS_max         = TolS
    solver.TolS_min         = 0.01*TolS
    solver.TolT             = TolT
    solver.LS_iters_target  = 30
#     @warn "CA"
#     solver.CA               = 1e-5/min(solver.TolS_max, solver.TolT)
    solver.CA_max           = Inf
    
    #Set initial and boundary conditions:
    solver.u0fun        = FW11((x) -> u0fun(x)) 
    solver.BC           = [BCW(BC_walls), BCW(BC_outlet), BCW(BC_walls), BCW(BC_inlet)]
    
    #-----------------------------------------------------------------------------
    #INITIAL CONDITION:
    
    #Compute initial condition:
    LIRKHyp_InitialCondition!(solver)

    #Open plots (variables, errors):
    figv                = Vector{Figure}(undef,2*length(PlotVars)+1)
    if !isempty(PlotVars)
        for ii=1:length(figv)
            figv[ii]    = figure()
        end
    end    
    t_lastFig           = 0.0
    NPlotSave           = 0
    NSave               = 0
    function PlotSol()
    
        #Compute variables at nodes:
        uNodes  = SolutionAtNodes(solver.u, solver.fes)
        xout    = DepVars(GasModel, solver.t, Vector{<:AMF64}(uNodes), 
                    Vector{<:AMF64}(uNodes), PlotVars)
            
        #Loop plot variables:
        for ii=1:length(PlotVars)
            figure(figv[ii].number)
            PyPlot.cla()
            PlotContour(reshape(xout[ii][1],:), solver.fes.PSpace)
            PlotMesh!(solver.mesh, color="k") 
            axis("equal")
            title(latexstring(LatexString(PlotVars[ii]),
                "; t^n=", sprintf1("%.2e", solver.t)), 
                fontsize=10)
            if SaveFig
                savefig("$(VideosUbi)ShockVortex_SC$(SC)_Mesh_$(PlotVars[ii])_$(NSave).png", dpi=400, pad_inches=0)
            end
        end
        
        for ii=1:length(PlotVars)
            figure(figv[length(PlotVars)+ii].number)
            PyPlot.cla()
            PlotContour(reshape(xout[ii][1],:), solver.fes.PSpace)
            axis("equal")
            title(latexstring(LatexString(PlotVars[ii]),
                "; t^n=", sprintf1("%.2e", solver.t)), 
                fontsize=10)
            if SaveFig
                savefig("$(VideosUbi)ShockVortex_SC$(SC)_$(PlotVars[ii])_$(NSave).png", dpi=400, pad_inches=0)
            end
        end
        
        figure(figv[length(figv)].number)
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
            savefig("$(VideosUbi)ShockVortex_Errors$(SC)_$(NSave).png", dpi=400, pad_inches=0)
        end
        
        NSave   += 1
        
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
        if solver.t-t_lastFig>=2e-3 || solver.t==tfv[length(tfv)] || NPlotSave==Nt_SaveFig
            NPlotSave   = 0
            t_lastFig   = solver.t
            PlotSol()
        end
        
        
    end
    
    #Save results:
    if SaveRes
        save("$(ResUbi)LIRKHyp_SC$(SC).jld2", "StudyCase", "ShockVortex",
            "ConvFlag", 1, "solver", save(solver) )
    end
    
    return solver
    
end

function test_LS()

    A           = load("../temp/test.jld2", "Am_LS").A
    u           = load("../temp/test.jld2", "u")
    b           = load("../temp/test.jld2", "b")
    eta_err     = load("../temp/test.jld2", "Am_LS").eta_err
    
    t_ini       = time()
    Am_LS       = LinearSystem2(A, eta_err=eta_err)
    println(time()-t_ini)
    
    u1          = copy(u)
    LSOutput    = LS_gmres!(Am_LS, u1, b, AbsTol=1e-8, Display="iter",
                            MaxIter=200)
    println(time()-t_ini)
    
    #-------------------------------------------
    
    Am_LS       = LinearSystem3(A, eta_err=eta_err)
    println(time()-t_ini)
    LinearSystem!(Am_LS, A, eta_err=eta_err)
    println(time()-t_ini)
    
    u2          = copy(u)
    LSOutput    = LS_gmres!(Am_LS, u2, b, AbsTol=1e-8, Display="iter",
                            MaxIter=200)
    println(time()-t_ini)
    
    display(norm(u1-u2, Inf))
    
    return
    
end

function test_LS2()

    A           = load("../temp/test.jld2", "AP")
    
    t_ini       = time()
    ilu         = ilu0(A)
    println(time()-t_ini)
    
    ILU         = ILU0_alloc(A)
    println(time()-t_ini)
    ILU0!(ILU, A)
    println(time()-t_ini)
    
    return
    
end

function test_LS3()

    A           = load("../temp/test.jld2", "AP")
    
    ILU         = ILU0_alloc(A)
    ILU0_rows!(ILU, A)
    L,U         = LUFactors_rows(ILU)
    
    ILU2        = ILU0_alloc(A)
    ILU0_ilu0!(ILU2, A)
    L2,U2       = LUFactors_rows(ILU2)
    
    ILU3        = ILU0_alloc(A)
    ILU0_columns!(ILU3, A)
    L3,U3       = LUFactors_columns(ILU3)
    
    AT          = permutedims(A)
    ILU4        = ILU0_alloc(AT)
    ILU0_rows!(ILU4, AT)
    U4,L4       = LUFactors_rows(ILU4)
    L4          = permutedims(L4)
    U4          = permutedims(U4)
    
    println(norm(A-L3*U3))
    println(norm(A-L4*U4))

    return A, L, U, L2, U2, L3, U3, L4, U4
    
end
