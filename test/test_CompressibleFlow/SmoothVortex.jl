include("test_CompressibleFlow.jl")

function SmoothVortex(hp::Float64, FesOrder::Int, tfv::Array{Float64,1};
    RKMethod::String="Ascher3", 
    epsilon::Float64=0.0, nu::Float64=0e-6, beta::Float64=0.0, kappa_rho_cv::Float64=0e-6,
    PlotVars::Vector{String}=String[], PlotCode::Vector{String}=fill("nodes", length(PlotVars)), 
    SaveFig::Bool=false, wFig::Float64=9.50, hFig::Float64=6.50, 
    mFig::Int=max(1,length(PlotCode)), nFig::Int=Int(ceil(length(PlotCode)/mFig)), Nt_SaveFig::Int=typemax(Int), cmap::String="jet",
    SC::Int=0, CSS::Float64=0.1, CDC::Float64=5.0, CFLa::Float64=1.0, Deltat0::Float64=1e-3,
    AMA_MaxIter::Int=200, TolS::Float64=1e-3, TolT::Float64=1e-3, 
    AMA_SizeOrder::Int=FesOrder, AMA_AnisoOrder::Int=2, AMA_ProjN::Int=1, AMA_ProjOrder::Int=0,
    SpaceAdapt::Bool=true, TimeAdapt::Bool=true,
    SaveRes::Bool=false, Nt_SaveRes::Int=typemax(Int), Deltat_SaveRes::Float64=0.01)

    #---------------------------------------------------------------------
    #PROBLEM DATA:
    
    #Problem data:
    gamma           = 1.4
    vortex_st       = 5.0   #vortex strength
    u_inf           = 1.0   #vortex horizontal velocity
    
    #Theoretical solution:
    function utheor(t::Float64, x::Vector{<:VecOrMat{Float64}})
        
        xm          = hcat(x[1][:], x[2][:])
        rho         = rho_SmoothVortex(xm, gamma, vortex_st, t, 0.0, u_inf, 0.0)
        rhovx       = rho .* vx_SmoothVortex(xm, gamma, vortex_st, t, 0.0, u_inf, 0.0)
        rhovy       = rho .* vy_SmoothVortex(xm, gamma, vortex_st, t, 0.0, u_inf, 0.0)
        rhoE        = @. rho * $e_SmoothVortex(xm, gamma, vortex_st, t, 0.0, u_inf, 0.0) + 
                                0.5*(rhovx^2 + rhovy^2)/rho
                                
        rho         = reshape(rho, size(x[1]))
        rhovx       = reshape(rhovx, size(x[1]))
        rhovy       = reshape(rhovy, size(x[1]))
        rhoE        = reshape(rhoE, size(x[1]))
        
        return [rho, rhovx, rhovy, rhoE]
        
    end
    
    #Boundary conditions:
    BC_horiz        = SlipAdiabatic()
    function uLeft(t::Float64, x::Vector{Matrix{Float64}})
        
        rho         = @mlv 0.0*x[1]+1.0
        rhovx       = @mlv 0.0*x[1]+1.0
        rhovy       = @mlv 0.0*x[1]
        
        return [rho, rhovx, rhovy]
        
    end
    BC_left         = SubsonicInlet1(FWt11((t,x)->uLeft(t,x)))
    function uRight(t::Float64, x::Vector{Matrix{Float64}})
        
        p           = @mlv 0.0*x[1]+1.0
        return [p]
        
    end
    BC_right        = SubsonicOutlet1(FWt11((t,x)->uRight(t,x)))
    
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
    
    #Mesh:
    MeshFile                = "$(@__DIR__)/../../temp/SmoothVortex_SC$(SC).geo"
    NX                      = Int(ceil(24.0/(hp*FesOrder)))
    NY                      = Int(ceil(14.0/(hp*FesOrder)))
    TrMesh_Rectangle_Create!(MeshFile, -7.0, 17.0, NX, -7.0, 7.0, NY)
    
    #Load LIRKHyp solver structure with default data. Modify the default data if necessary:
    solver                  = LIRKHyp_Start(GasModel)
    solver.ProblemName      = "SmoothVortex"
    solver.SC               = SC
#     solver.MeshFile         = "$(MeshUbi)SmoothVortex/MeshCase$(MeshCase).geo"
    solver.MeshFile         = MeshFile
    solver.nBounds          = 4
    solver.FesOrder         = FesOrder
    solver.RKMethod         = RKMethod
    solver.Deltat0          = Deltat0
    solver.tf               = tfv[length(tfv)]
    solver.AMA_MaxIter      = AMA_MaxIter
    solver.AMA_SizeOrder    = AMA_SizeOrder
    solver.AMA_AnisoOrder   = AMA_AnisoOrder
    solver.AMA_ProjN        = AMA_ProjN
    solver.AMA_ProjOrder    = AMA_ProjOrder
    solver.TolS_max         = TolS
    solver.TolS_min         = 0.01*TolS
    solver.TolT             = TolT
    solver.SpaceAdapt       = SpaceAdapt
    solver.TimeAdapt        = TimeAdapt
    
    #Set initial and boundary conditions:
    solver.u0fun        = FW11((x) -> utheor(0.0, x)) 
    solver.BC           = [BCW(BC_horiz), BCW(BC_right), BCW(BC_horiz), BCW(BC_left)]
    
    #-----------------------------------------------------------------------------
    #INITIAL CONDITION:
    
    #Compute initial condition:
    ConvFlag            = LIRKHyp_InitialCondition!(solver)
#     CheckJacobian(solver, Plot_dQ_du=true, Plot_dQ_dgradu=true)
#     BC_CheckJacobian(solver, 4, Plot_df_du=true, Plot_df_dgradu=true)
#     return
    
    #Change TolT:
    if TolT==0.0
        TolT            = 0.01*solver.etaS
        solver.TolT     = TolT
    end
    
    #Compute Lq error:
    errLq,              = LqError(solver, FW11((x) -> utheor(solver.t, x)), q=2.0)
    hmean               = 2.0*sqrt(solver.Omega/solver.mesh.nElems/TrElem_Area)
    errL2L2             = errLq
    etaL2L2             = solver.etaS
    
    println("hmean=", sprintf1("%.2e", hmean), ", e_L2L2=", sprintf1("%.2e", errL2L2))
    
    #Open plots (variables, errors):
    figv                = Vector{Figure}(undef,4)
    if !isempty(PlotVars)
        figv[1]         = PyPlotSubPlots(mFig, nFig, w=wFig, h=hFig, left=0.9, right=0.4, bottom=1.1, top=1.0)
        figv[2]         = figure()
        figv[3]         = figure()
        figv[4]         = figure()
    end    
    t_lastFig           = 0.0
    NSave               = 0
    NPlotSave           = 0
    function PlotSol()
        
        if isempty(PlotVars)
            return
        end
        
        figure(figv[1].number)
        
        #Theoretical solution:
        x1v     = linspace(-7.0, 17.0, 500)
        x2v     = zeros(size(x1v))
        xm      = [ x1v, x2v ]
        um      = utheor(solver.t, xm)
        
        #Interpolate solution:
        u_terp, = SolutionCompute(solver.u, solver.fes, [x1v, x2v])
            
        #Loop plot variables:
        for ii=1:length(PlotVars)
            PyPlot.subplot(mFig, nFig, ii)
            PyPlot.cla()
#             plot3D(solver, PlotVars[ii])
#             plot(xm[1], DepVars(solver.model, um, [PlotVars[ii]])[1][1], "k", linewidth=0.5)
#             plotNodes((x1,x2)-> x1, solver, PlotVars[ii])
            plot(x1v, DepVars(solver.model, solver.t, Vector{<:AMF64}(um), Vector{<:AMF64}(um), [PlotVars[ii]])[1][1][:], "--k", linewidth=0.5)
            plot(x1v, DepVars(solver.model, solver.t, Vector{<:AMF64}(u_terp), Vector{<:AMF64}(u_terp), [PlotVars[ii]])[1][1][:], "b", linewidth=0.5)
            xlabel(latexstring("x_1"), fontsize=10)
            title(latexstring(LatexString(PlotVars[ii]),
                "; t^n=", sprintf1("%.2e", solver.t)), 
                fontsize=10)
        end
        if SaveFig
            savefig("$(VideosUbi)SmoothVortex_SC$(SC)_$(NSave).png", dpi=400, pad_inches=0)
        end
        
        figure(figv[2].number)
        PyPlot.cla()
        plot(solver.mesh)
        axis("equal")
        if SaveFig
            savefig("$(VideosUbi)SmoothVortex_Mesh$(SC)_$(NSave).png", dpi=400, pad_inches=0)
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
            savefig("$(VideosUbi)SmoothVortex_Errors$(SC)_$(NSave).png", dpi=400, pad_inches=0)
        end
        
        figure(figv[4].number)
        semilogy(solver.tv, solver.nElemsv, ".-b")
        if true
            validv  = solver.validv .== 1
            semilogy(solver.tv[validv], solver.nElemsv[validv], "sb")
        end
        xlabel(L"t")
        ylabel(L"N_E")
        if SaveFig && solver.t==tfv[length(tfv)]
            savefig("$(VideosUbi)SmoothVortex_NE$(SC)_$(NSave).png", dpi=400, pad_inches=0)
        end
        
        NSave   += 1
        
    end
    PlotSol()
    
    #Function to save intermediate results:
    t_lastRes           = -Deltat_SaveRes
    ct_SaveRes          = 0
    nb_SaveRes          = 0
    function SaveSol()
        
        ct_SaveRes      += 1 
        if SaveRes && ( solver.t-t_lastRes>=Deltat_SaveRes || 
                        ct_SaveRes==Nt_SaveRes || solver.t==tfv[length(tfv)] || solver.t==0.0 )
            save("$(ResUbi)LIRKHyp_SC$(SC)_$(nb_SaveRes).jld2", "StudyCase", "SmoothVortex",
                "ConvFlag", ConvFlag, "solver", save(solver), 
                "epsilon", epsilon, "nu", nu, "beta", beta, "kappa_rho_cv", kappa_rho_cv,
                "errL2L2", errL2L2, "etaL2L2", etaL2L2 )
            t_lastRes   += Deltat_SaveRes
            ct_SaveRes  = 0
            nb_SaveRes  += 1
        end
        return
        
    end
    SaveSol()
    
    #-----------------------------------------------------------------------------
    #MARCH IN TIME:
    
#     flag    = LIRKHyp_Step!(solver)
#     return solver
    
    while solver.t<tfv[length(tfv)]
    
        ConvFlag            = LIRKHyp_Step!(solver)
        
        if ConvFlag<=0
            break
        end
        
        NPlotSave   += 1
        if solver.t-t_lastFig>=0.1 || solver.t==tfv[length(tfv)] || NPlotSave==Nt_SaveFig
            NPlotSave   = 0
            t_lastFig   = solver.t
            PlotSol()
        end
        
        #Compute Lq error:
        errLq,              = LqError(solver, FW11((x) -> utheor(solver.t, x)), q=2.0)
        hmean               = 2.0*sqrt(solver.Omega/solver.mesh.nElems/TrElem_Area)
        errL2L2             = errL2L2_(solver.t-solver.Deltat, errL2L2, solver.t, errLq)
        etaL2L2             = errL2L2_(solver.t-solver.Deltat, etaL2L2, solver.t, solver.etaS+solver.etaT)
        
        println("hmean=", sprintf1("%.2e", hmean), 
                ", Deltat_mean=", sprintf1("%.2e", solver.t/solver.Nt), 
                ", e_L2L2=", sprintf1("%.2e", errL2L2), 
                ", etaL2L2=", sprintf1("%.2e", etaL2L2))
        
        SaveSol()
        
    end
    
#     PlotSol()

    #Save results:
#     if SaveRes
#         save("$(ResUbi)LIRKHyp_SC$(SC)_1.jld2", "StudyCase", "SmoothVortex",
#             "ConvFlag", ConvFlag, "solver", save(solver), 
#             "errL2L2", errL2L2, "etaL2L2", etaL2L2 )
#     end
    
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
