include("test_CompressibleFlow.jl")

function DetonationWave(FesOrder::Int;
    tf::Float64=20.0, RKMethod::String="Ascher3", #"BPR3"
    delta::Float64=1e-2,
    epsilon::Float64=1000*delta,
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

    #Define gas:
    GasModel                = GasH2()
    GasModel.epsilon        = epsilon
    GasModel.nu             = 0.0
    GasModel.beta           = 0.0
    GasModel.kappa_rho_cv   = 0.0
    GasModel.D              = 0.0

    #Define domain:
    MeshFile            = "$(@__DIR__)/../../temp/DetonationWave_SC$(SC).geo"
    #Definition of the initial rectangle: x1,x2,Nx,y1,y2,Ny:
    TrMesh_Rectangle_Create!(MeshFile, 0.0, 0.45, 100, 0.0, 0.0045, 2) 

    #Initial conditions:
    W        = GasModel.Wi
    W        = [W[6], W[2], W[1], W[10], W[3]]
    R        = GasModel.R

    # 0-0.015: Ar (6), H2O (3), OH (10)
    X1       = [1.0, 0.0, 0.0, 0.01/8, 2/8]
    T1       = 3500.0
    p1       = 5.5e5
    R1       = R*sum(X1./sum(X1.*W))
    rho1     = p1/(R1*T1)
    Y1       = [(X1[i]*W[i])/sum(X1.*W) for i = 1:5]
    Y1_Ar    = Y1[1]
    Y1_O2    = Y1[2]
    Y1_H2    = Y1[3]
    Y1_OH    = Y1[4]
    Y1_H2O   = Y1[5]
    e1       = calc_ei(GasModel, fill(T1,1,1))
    E1       = sum(Y1.*e1[[6,2,1,10,3]])

    # 0.015-0.025: Ar (1), O2 (2), H2 (1), OH (10)
    X2       = [1.0, 1/7, 2/7, 0.01/7, 0.0]
    T2       = 350.0
    p2       = 6.67e3
    R2       = R*sum(X2./sum(X2.*W))
    rho2     = p2/(R2*T2)
    Y2       = [(X2[i]*W[i])/sum(X2.*W) for i = 1:5]
    Y2_Ar    = Y2[1]
    Y2_O2    = Y2[2]
    Y2_H2    = Y2[3]
    Y2_OH    = Y2[4]
    Y2_H2O   = Y2[5]
    e2       = calc_ei(GasModel, fill(T2,1,1))
    E2       = sum(Y2.*e2[[6,2,1,10,3]])

    # 0.025-0.45: Ar (6), O2 (2), H2 (1)
    X3       = [1.0, 1/7, 2/7, 0.0, 0.0]
    T3       = 298.0
    p3       = 6.67e3
    R3       = R*sum(X3./sum(X3.*W))
    rho3     = p3/(R3*T3)
    Y3       = [(X3[i]*W[i])/sum(X3.*W) for i = 1:5]
    Y3_Ar    = Y3[1]
    Y3_O2    = Y3[2]
    Y3_H2    = Y3[3]
    Y3_OH    = Y3[4]
    Y3_H2O   = Y3[5]
    e3       = calc_ei(GasModel, fill(T3,1,1))
    E3       = sum(Y3.*e3[[6,2,1,10,3]])

    function u0fun(x::Vector{Matrix{Float64}})

        rhoY_Ar     = @mlv SmoothHeaviside(x[1] - 0.015, delta, rho1*Y1_Ar, rho2*Y2_Ar) + SmoothHeaviside(x[1] - 0.025, delta, 0.0, rho3*Y3_Ar - rho2*Y2_Ar)
        rhoY_O2     = @mlv SmoothHeaviside(x[1] - 0.015, delta, rho1*Y1_O2, rho2*Y2_O2) + SmoothHeaviside(x[1] - 0.025, delta, 0.0, rho3*Y3_O2 - rho2*Y2_O2)
        rhoY_H2     = @mlv SmoothHeaviside(x[1] - 0.015, delta, rho1*Y1_H2, rho2*Y2_H2) + SmoothHeaviside(x[1] - 0.025, delta, 0.0, rho3*Y3_H2 - rho2*Y2_H2)
        rhoY_OH     = @mlv SmoothHeaviside(x[1] - 0.015, delta, rho1*Y1_OH, rho2*Y2_OH) + SmoothHeaviside(x[1] - 0.025, delta, 0.0, rho3*Y3_OH - rho2*Y2_OH)
        rhoY_H20    = @mlv SmoothHeaviside(x[1] - 0.015, delta, rho1*Y1_H2O, rho2*Y2_H2O) + SmoothHeaviside(x[1] - 0.025, delta, 0.0, rho3*Y3_H2O - rho2*Y2_H2O)
        rhoY_N2     = x[1].*0.0
        rhoY_He     = x[1].*0.0
        rhoY_CO     = x[1].*0.0
        rhoY_CO2    = x[1].*0.0
        rhoY_H      = x[1].*0.0
        rhoY_H2O2   = x[1].*0.0
        rhoY_HO2    = x[1].*0.0
        rhoY_O      = x[1].*0.0

        rhovx       = x[1].*0.0
        rhovy       = x[1].*0.0

        rhoE        = @mlv SmoothHeaviside(x[1] - 0.015, delta, rho1*E1[1], rho2*E2[1]) + SmoothHeaviside(x[1] - 0.025, delta, 0.0, rho3*E3[1] - rho2*E2[1])

        return [rhoY_H2, rhoY_O2, rhoY_H20, rhoY_N2, rhoY_He, rhoY_Ar, rhoY_CO, rhoY_CO2, rhoY_H, rhoY_OH, rhoY_HO2, rhoY_H2O2, rhoY_O, rhovx, rhovy, rhoE]

    end

    #Boundary conditions:
    BC_walls        = SlipAdiabatic()

    #---------------------------------------------------------------------
    #PRE-PROCESS STAGE:

    #Load LIRKHyp solver structure with default data. Modify the default data if necessary:
    solver                  = LIRKHyp_Start(GasModel)
    solver.ProblemName      = "DetonationWave"
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
    solver.BC           = [BCW(BC_walls), BCW(BC_walls), BCW(BC_walls), BCW(BC_walls)] # Abajo,dcha,arriba,izq

    #-----------------------------------------------------------------------------
    #INITIAL CONDITION:

    #Compute initial condition:
    ConvFlag            = LIRKHyp_InitialCondition!(solver)
#     CheckJacobian(solver, Plot_dQ_du=true, Plot_df_du=true)

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
                savefig("$(VideosUbi)DetonationWave_SC$(SC)_$(nb_SaveFig).png", dpi=400, pad_inches=0)
            end

            figure(figv[2].number)
            PyPlot.cla()
            plot(solver.mesh)
            if SaveFig
                savefig("$(VideosUbi)DetonationWave_Mesh$(SC)_$(nb_SaveFig).png", dpi=400, pad_inches=0)
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
                savefig("$(VideosUbi)DetonationWave_Errors$(SC)_$(nb_SaveFig).png", dpi=400, pad_inches=0)
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
            save("$(ResUbi)LIRKHyp_SC$(SC)_$(nb_SaveRes).jld2", "StudyCase", "DetonationWave",
                 "ConvFlag", ConvFlag, "solver", save(solver),
                 "epsilon", epsilon,
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
