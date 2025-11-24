include("test_CompressibleFlow.jl")

function Sedov(hp0::Float64, FesOrder::Int;
    tf::Float64=1.0, RKMethod::String="Ascher3",
    epsilon::Float64=1e-4, epsilon0::Float64=1e-3, t_epsilon::Float64=0.1, 
    nu::Float64=0e-4, beta::Float64=0.0, kappa_rho_cv::Float64=0e-6, 
    delta0::Float64=1e-2, pInf::Float64=1e-9, 
    #
    TolS::Float64=1e-6, AMA_MaxIter::Int=200, AMA_SizeOrder::Int=FesOrder, AMA_AnisoOrder::Int=2, 
    #
    TolT::Float64=1e-3, Deltat0::Float64=1e-8,
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
    alpha           = 0.9840418765432098
    rho0            = 1.0
    E0              = alpha*1.0^4    #alpha*Rf^4
    
    #Define gas:
    GasModel                = GasIdeal()
    GasModel.gamma          = gamma
    GasModel.epsilon        = epsilon #This viscosity will be modified in the while t<tf loop
    GasModel.nu             = nu
    GasModel.beta           = beta
    GasModel.kappa_rho_cv   = kappa_rho_cv
    GasModel.CSS            = CSS

    #Initial conditions:
    function u0fun(x::Vector{Matrix{Float64}})
        
        rho         = @mlv 0.0*x[1] + rho0
        rhovx       = @mlv 0.0*x[1] 
        rhovy       = @mlv 0.0*x[1] 
        
        rhoEinf     = pInf/(gamma-1.0)
        Delta_rhoE  = (E0-rhoEinf*2.4^2)/(pi*delta0^2)
        rhoE        = @mlv rhoEinf + Delta_rhoE*exp(-(x[1]^2+x[2]^2)/delta0^2)
        
        return [rho, rhovx, rhovy, rhoE]
        
    end
    
    #Boundary conditions:
    BC_walls        = SlipAdiabatic()
    function pOutlet(t::Float64, x::Vector{Matrix{Float64}})
        
        p           = @mlv 0.0*x[1]+pInf
        return [p]
        
    end
    BC_outlet       = SubsonicOutlet1(FWt11((t,x)->pOutlet(t,x)))
    
    #Theoretical solution:
    function utheor(t::Float64, x::Vector{<:VecOrMat{Float64}})
        Rshock          = (E0*t^2/rho0/alpha)^(1/4)
        ushock          = 0.5*Rshock/(t+1e-14)
        rho1            = 1.0
        rho2            = (gamma+1.0)/(gamma-1.0)*rho1
        v2              = 2.0/(gamma+1.0)*ushock
        p2              = 2.0/(gamma+1.0)*rho1*ushock^2.0
        xi              = @mlv sqrt(x[1]^2.0+x[2]^2.0)/(Rshock+1e-14)
        rho             = @mlv rho2 * $Sedov2D_frho(xi)
        vr              = @mlv v2 * $Sedov2D_fv(xi)
        p               = @mlv p2 * $Sedov2D_fp(xi)
        rhovx           = @mlv rho * vr * x[1]/sqrt.(x[1]^2.0 + x[2]^2.0 .+ 1e-10)
        rhovy           = @mlv rho * vr * x[2]/sqrt.(x[1]^2.0 + x[2]^2.0 .+ 1e-10)
        rhoE            = @mlv p/(gamma-1.0) .+ 0.5*rho*vr^2.0 #Note that vr^2=vx^2+vy^2
        return [ rho, rhovx, rhovy, rhoE ]
    end
    
    #---------------------------------------------------------------------
    #PRE-PROCESS STAGE:
    
    #Mesh:
    MeshFile                = "$(@__DIR__)/../../temp/Sedov_SC$(SC).geo"
    NX                      = Int(ceil(2.4/(hp0*FesOrder)))
    NY                      = Int(ceil(2.4/(hp0*FesOrder)))
    TrMesh_Rectangle_Create!(MeshFile, -1.2, 1.2, NX, -1.2, 1.2, NY)
    
    #Load LIRKHyp solver structure with default data. Modify the default data if necessary:
    solver                  = LIRKHyp_Start(GasModel)
    solver.ProblemName      = "Sedov"
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
#     solver.BC           = [BCW(BC_walls), BCW(BC_outlet), BCW(BC_outlet), BCW(BC_walls)]
    solver.BC           = [BCW(BC_outlet), BCW(BC_outlet), BCW(BC_outlet), BCW(BC_outlet)]
    
    #-----------------------------------------------------------------------------
    #INITIAL CONDITION:
    
    #Compute initial condition:
    ConvFlag            = LIRKHyp_InitialCondition!(solver)

    #Function to plot solution:
    figv                = Vector{Figure}(undef,6)
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
    
            #Max R to plot:
            Rshock          = (E0*solver.t^2/rho0/alpha)^(1/4)
            Rplot           = max(10.0*delta0, 1.2*Rshock)
            
            #Theoretical solution:
            xm      = [ linspace(1e-10, Rplot, 100), zeros(100) ]
            um      = utheor(solver.t, xm)
        
            #Loop plot variables:
            for ii=1:length(PlotVars)
                PyPlot.subplot(mFig, nFig, ii)
                PyPlot.cla()
    #             plot3D(solver, PlotVars[ii])
                plot(xm[1], 
                    DepVars(solver.model, solver.t, 
                            Vector{<:AMF64}(xm), 
                            Vector{<:AMF64}(um), 
                            [PlotVars[ii]])[1][1],
                    "k", linewidth=0.5)
                splot_fun(x1,x2)    = @mlv sqrt(x1^2+x2^2)
                PlotNodes(splot_fun, solver, PlotVars[ii])
                xlim([0.0,Rplot])
                xlabel(latexstring("r"), fontsize=10)
                title(latexstring(LatexString(GasModel, PlotVars[ii]),
                    "; t^n=", sprintf1("%.2e", solver.t)), 
                    fontsize=10)
            end
            if SaveFig
                savefig("$(VideosUbi)Sedov_SC$(SC)_$(nb_SaveFig).png", dpi=400, pad_inches=0)
            end
        
            figure(figv[2].number)
            PyPlot.cla()
            PlotContour(solver, solver.model, "rho", delta=1e-5)
            plot(solver.mesh, color="w")
            xlim([0.0,Rplot])
            ylim([0.0,Rplot])
    #         axis("equal")
            if SaveFig
                savefig("$(VideosUbi)Sedov_MeshZoom$(SC)_$(nb_SaveFig).png", dpi=400, pad_inches=0)
            end
        
            figure(figv[3].number)
            PyPlot.cla()
            plot(solver.mesh)
            axis("equal")
            if SaveFig
                savefig("$(VideosUbi)Sedov_Mesh$(SC)_$(nb_SaveFig).png", dpi=400, pad_inches=0)
            end
            
            figure(figv[4].number)
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
                savefig("$(VideosUbi)Sedov_Errors$(SC)_$(nb_SaveFig).png", dpi=400, pad_inches=0)
            end
            
            figure(figv[5].number)
            PyPlot.cla()
            Deltatv     = [0.0; diff(solver.tv)];
            markv       = @. Deltatv > 0.0
            plot(solver.tv[markv], Deltatv[markv], ".-k")
            xlabel(latexstring("t"))
            ylabel(latexstring("\\tau"), rotation=0)
            
            figure(figv[6].number)
            PyPlot.cla()
            plot(solver.tv, solver.nElemsv, ".-k")
            xlabel(latexstring("t"))
            ylabel(latexstring("NE"), rotation=0)
            
            #=
            figure(figv[5].number)
            PyPlot.cla()
            PlotContour(solver.u[1], solver.fes, 
                zv=linspace(-1e-10, (gamma+1.0)/(gamma-1.0)+1e-10, 100))
            PlotMesh!(solver.mesh, color="k") 
            if SaveFig && solver.t==tf
                savefig("$(VideosUbi)Sedov_Contour$(SC)_$(nb_SaveFig).png", dpi=400, pad_inches=0)
            end
            axis("equal")
            =#
            
            #=
            figure(figv[6].number)
            PyPlot.cla()
            PlotContour(solver.u[1], solver.fes, zv=linspace(-1e-10, (gamma+1.0)/(gamma-1.0)+1e-10, 100))
            PlotMesh!(solver.mesh, color="k") 
            if SaveFig && solver.t==tf
                savefig("$(VideosUbi)Sedov_ContourZoom$(SC)_$(nb_SaveFig).png", dpi=400, pad_inches=0)
            end
            xlim([-Rplot,Rplot])
            ylim([-Rplot,Rplot])
            =#
            
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
            save("$(ResUbi)LIRKHyp_SC$(SC)_$(nb_SaveRes).jld2", "StudyCase", "Sedov",
                "ConvFlag", ConvFlag, "solver", save(solver), 
                "epsilon", epsilon, "epsilon0", epsilon0, "nu", nu, "beta", beta, "kappa_rho_cv", kappa_rho_cv,
                "delta0", delta0, "pInf", pInf, "TolS", TolS, "TolT", TolT)
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
    
        #Modify viscosity:
        GasModel.epsilon        = epsilon + epsilon0*exp(-solver.t/t_epsilon);
        
        #Solve:
        ConvFlag    = LIRKHyp_Step!(solver)
        
        if ConvFlag<=0
            break
        end
        
        PlotSol()
        SaveSol()
        
    end
    
    #Save results:
    if SaveRes
        save("$(ResUbi)LIRKHyp_SC$(SC)_1000.jld2", "StudyCase", "Sedov",
            "ConvFlag", ConvFlag, "solver", save(solver) )
    end
    
    return solver
    
end
