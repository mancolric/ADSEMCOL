include("test_CompressibleFlow.jl")

#DetonationWave(tf=1.0, PlotFig=true, PlotVars=["X_Ar", "X_H2O", "X_OH", "X_O2", "X_H2", "p", "T"], delta=1e-3);
#El cálculo del jacobiano aumenta mucho, pero el número de elementos no
#aumenta tanto. Qué sucede??

#Section 5.3 in Johsonn and Kercher, JCP2020:

function DetonationWave(; FesOrder::Int=5, 
    tf::Float64=20.0, RKMethod::String="Ascher3", #"BPR3"
    delta::Float64=1e-3,
    epsilon::Float64=1000*delta,
    #
    TolS::Float64=1e-6, AMA_MaxIter::Int=200, AMA_SizeOrder::Int=FesOrder, AMA_AnisoOrder::Int=2,
    #
    TolT::Float64=1e-5, Deltat0::Float64=1e-8,
    #
    PlotFig::Bool=false, wFig::Float64=9.50, hFig::Float64=6.50,
    PlotVars::Vector{String}=String[], PlotCode::Vector{String}=fill("nodes", length(PlotVars)),
    SaveFig::Bool=false, Nt_SaveFig::Int=3, Deltat_SaveFig::Float64=Inf,
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
    MeshFile                = "$(@__DIR__)/../../temp/DetonationWave_SC$(SC).geo"
    #Definition of the initial rectangle: x1,x2,Nx,y1,y2,Ny:
    TrMesh_Rectangle_Create!(MeshFile, 0.0, 0.45, 100, 0.0, 0.0045, 2) 

    #Initial conditions:
    Wi       = GasModel.Wi
    R        = GasModel.R
    #Remind that the species are:
    # 1-H2  2-O2  3-H2O  4-N2  5-He  6-Ar  7-CO  8-CO2  9-H  10-OH  11-HO2 
    #12-H2O2 13-O
    
    #State 1: 0-0.015: Ar (6), H2O (3), OH (10)
    X1       = zeros(13)
    X1[6]    = 8.0
    X1[3]    = 2.0
    X1[10]   = 0.01
    X1       /= sum(X1)
    #
    T1       = 3500.0
    p1       = 5.5e5
    #Rgas = Rbar/W, W=sum_i X_i W_i
    W1       = sum(X1.*Wi)
    R1       = R/W1
    #
    rho1     = p1/(R1*T1)
    #Y_i = m_i/m = n_i W_i / (n W) = X_i W_i/W
    Y1       = @. X1*Wi/W1;
    #Energy per unit mass:
    e1       = calc_ei(GasModel, fill(T1,1,1))
    E1       = sum(Y1.*e1)
    
    #State 2: 0.015-0.025: 
    X2       = zeros(13)
    X2[6]    = 7.0
    X2[2]    = 1.0
    X2[1]    = 2.0
    X2[10]   = 0.01
    X2       /= sum(X2)
    #
    T2       = 350.0
    p2       = 6.67e3
    #Rgas = Rbar/W, W=sum_i X_i W_i
    W2       = sum(X2.*Wi)
    R2       = R/W2
    #
    rho2     = p2/(R2*T2)
    #Y_i = m_i/m = n_i W_i / (n W) = X_i W_i/W
    Y2       = @. X2*Wi/W2;
    #Energy per unit mass:
    e2       = calc_ei(GasModel, fill(T2,1,1))
    E2       = sum(Y2.*e2)

    #State 3: 0.025-0.45:
    X3       = zeros(13)
    X3[6]    = 7.0
    X3[2]    = 1.0
    X3[1]    = 2.0
    X3       /= sum(X3)
    #
    T3       = 298.0
    p3       = 6.67e3
    #
    W3       = sum(X3.*Wi)
    R3       = R/W3
    #
    rho3     = p3/(R3*T3)
    #Y_i = m_i/m = n_i W_i / (n W) = X_i W_i/W
    Y3       = @. X3*Wi/W3;
    #Energy per unit mass:
    e3       = calc_ei(GasModel, fill(T3,1,1))
    E3       = sum(Y3.*e3)

    function u0fun(x::Vector{Matrix{Float64}})

        u0          = Vector{Matrix{Float64}}(undef, GasModel.nVars)
        for II=1:GasModel.nSpecies
            u0[II]  = @. SmoothHeaviside(x[1] - 0.015, delta, 
                                    rho1*Y1[II], rho2*Y2[II]) + 
                        + SmoothHeaviside(x[1] - 0.025, delta, 
                            0.0, rho3*Y3[II]-rho2*Y2[II])
        end
        u0[GasModel.nSpecies+1]     = x[1].*0.0
        u0[GasModel.nSpecies+2]     = x[1].*0.0
        u0[GasModel.nSpecies+3]     = @. SmoothHeaviside(x[1] - 0.015, delta, 
                                            rho1*E1[1], rho2*E2[1]) + 
                                        + SmoothHeaviside(x[1] - 0.025, delta, 
                                                0.0, rho3*E3[1] - rho2*E2[1])

        return u0

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
        figv[1]         = PyPlotSubPlots(mFig, nFig, w=wFig, h=hFig, left=1.5, right=0.4, bottom=1.1, top=1.0)
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
                xlabel(latexstring("x_1"), fontsize=8)
                tick_params(axis="both", which="minor", labelsize=8)
                title(latexstring(LatexString(GasModel, PlotVars[ii]),
                                  "; t^n=", sprintf1("%.2e", solver.t)),
                fontsize=8)
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
