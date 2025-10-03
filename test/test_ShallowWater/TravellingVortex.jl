cd(@__DIR__)
include("test_ShallowWater.jl")

function TravellingVortex(hp0::Float64, FesOrder::Int;
    tf::Float64=10.0, RKMethod::String="BPR3",
    epsilon::Float64=0e-3, h0::Float64=1.0, g::Float64=9.8,
    #
    TolS::Float64=1e-4, AMA_MaxIter::Int=200, AMA_SizeOrder::Int=FesOrder, AMA_AnisoOrder::Int=2,
    SpaceAdapt::Bool=true, 
    #
    TolT::Float64=1e-4, Deltat0::Float64=1e-5, TimeAdapt::Bool=true,
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

    function u0fun(x::Vector{Matrix{Float64}})

        r               = @tturbo @. sqrt(x[1]^2 + x[2]^2 + 1e-10)
        h               = @tturbo @. h0 - (0.5*exp(-(r^2 - 1)))/g
        v_phi           = @tturbo @. r*exp(-0.5*(r^2 - 1))
        q1              = @tturbo @. -h*v_phi*x[2]/r + h       
        q2              = @tturbo @. h*v_phi*x[1]/r
        b               = @tturbo @. 0.0*x[1]
        eta             = h
        
        return [eta, q1, q2, b]

    end

    function utheorfun(t::Float64, x::Vector{Matrix{Float64}})

        x0      =  [x[1].-1.0*t, x[2]]
        return u0fun(x0)

    end

    #Boundary conditions:
    BC_horiz        = SlipAdiabatic()
    function uLeft(t::Float64, x::Vector{Matrix{Float64}})
        q1              = @tturbo @. 0.0*x[1] + h0*1.0
        q2              = @tturbo @. 0.0*x[1]
        return [q1, q2]

    end
    BC_left         = SubsonicInlet1(FWt11((t,x)->uLeft(t,x)))
    function uRight(t::Float64, x::Vector{Matrix{Float64}})
        h               = @tturbo @. 0.0*x[1] + h0
        eta             = h
        return [eta]
    end
    BC_right        = SubsonicOutlet1(FWt11((t,x)->uRight(t,x)))
    
    #---------------------------------------------------------------------
    #PRE-PROCESS STAGE:

    #Mesh:
    MeshFile                = "../temp/SteadyVortex$(SC).geo"
    NX                      = Int(ceil(7.0/(hp0*FesOrder)))
    NY                      = Int(ceil(3.0/(hp0*FesOrder)))
    TrMesh_Rectangle_Create!(MeshFile, -7.0, 17.0, NX, -7.0, 7.0, NY)

    #Load LIRKHyp solver structure with default data. Modify the default data if necessary:
    solver                  = LIRKHyp_Start(model)
    solver.ProblemName      = "SteadyVortex" 
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
    solver.TolS_min         = 0.0*TolS
    solver.SpaceAdapt       = SpaceAdapt
    solver.TolT             = TolT
    solver.TimeAdapt        = TimeAdapt

    #Set initial and boundary conditions:
#     solver.u0fun        = FW11((x) -> u0fun(x))
#     BC_walls            = SlipAdiabatic()
#     solver.BC           = [BCW(BC_walls), BCW(BC_walls), BCW(BC_walls), BCW(BC_walls)]

    #Set initial and boundary conditions:
    solver.u0fun        = FW11((x) -> u0fun(x))
    solver.BC           = [BCW(BC_horiz), BCW(BC_right), BCW(BC_horiz), BCW(BC_left)]

    #-----------------------------------------------------------------------------
    #INITIAL CONDITION:

    #Compute initial condition:
    ConvFlag            = LIRKHyp_InitialCondition!(solver)
#     CheckJacobian(solver, Plot_df_du=false, Plot_df_dgradu=true, 
#         Plot_dQ_du=false, Plot_dQ_dgradu=false)
#     for ii = 4
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
                title(latexstring(LatexString(PlotVars[ii]),
                    "; t^n=", sprintf1("%.2e", solver.t)),
                    fontsize=10)
                println(PlotVars[ii], ": min=", minimum(v_plot), ", max=", maximum(v_plot))
            end
            if SaveFig
                savefig("$(VideosUbiTFG)SteadyVortex$(SC)_$(nb_SaveFig).png", dpi=400, pad_inches=0)
            end
            figure(figv[2].number)
            #Loop plot variables:
            for ii=1:length(PlotVars)
                PyPlot.subplot(mFig, nFig, ii)
                PyPlot.subplots_adjust(hspace=0.8)
                PyPlot.cla()

                splot_fun(x1,x2)    = @mlv x1
                PlotNodes(splot_fun, solver, PlotVars[ii])
                xlabel(latexstring("x_1"), fontsize=10)
                title(latexstring(LatexString(PlotVars[ii]),
                                  "; t^n=", sprintf1("%.2e", solver.t)),
                fontsize=10)
            end
            if SaveFig
                savefig("$(VideosUbiTFG)SteadyVortex_Pts$(SC)_$(nb_SaveFig).png", dpi=400, pad_inches=0)
            end


#=
            figure(figv[2].number)
            PyPlot.cla()
            PlotContour(solver.u[1], solver.fes)
            PlotMesh!(solver.mesh, color="w")
            if SaveFig
                savefig("$(VideosUbi)Sod_Problem_Mesh$(SC)_$(nb_SaveFig).png", dpi=400, pad_inches=0)
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
                savefig("$(VideosUbi)Sod_Problem_Errors$(SC)_$(nb_SaveFig).png", dpi=400, pad_inches=0)
            end=#

            t_lastFig           += Deltat_SaveFig
            ct_SaveFig          = 0
            nb_SaveFig          += 1

#             figure(figv[4].number)
#
#             #Loop plot variables:
#             for ii=1:length(PlotVars)
#                 PyPlot.subplot(mFig, nFig, ii)
#                 PyPlot.cla()
#
#                 splot_fun(x1,x2)    = @mlv x1
#                 PlotNodes(splot_fun, solver, PlotVars[ii])
#                 xlabel(latexstring("x_1"), fontsize=10)
#                 title(latexstring(LatexString(PlotVars[ii]),
#                                   "; t^n=", sprintf1("%.2e", solver.t)),
#                 fontsize=10)
#             end
#             if SaveFig
#                 savefig("$(VideosUbi)Sod_Problem_Pts$(SC)_$(nb_SaveFig).png", dpi=400, pad_inches=0)
#             end

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
            save("$(ResUbi)LIRKHyp_SC$(SC)_$(nb_SaveRes).jld2", "StudyCase", "SteadyVortex",
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

    while solver.t<tf

        ConvFlag    = LIRKHyp_Step!(solver)

        if ConvFlag<=0
            break
        end


        PlotSol()
        SaveSol()

    end

    errLq,              = LqError(solver, FW11((x) -> utheorfun(tf, x)), q=2.0)
    hmean               = 2.0*sqrt(solver.Omega/solver.mesh.nElems/TrElem_Area)
    Deltat_mean         = solver.tf/(solver.Nt-1)
    println("hmean=", sprintf1("%.2e", hmean), ", Deltat_mean=", sprintf1("%.2e", Deltat_mean),
            ", e_Lq=", sprintf1("%.2e", errLq))

    #Save results:
    if SaveRes
        save("$(ResUbi)LIRKHyp_SC$(SC)_1000.jld2", "StudyCase", "SteadyVortex",
            "ConvFlag", ConvFlag, "solver", save(solver), "hmean", hmean, "e_Lq", errLq, "Deltat_mean", Deltat_mean)
    end

    return solver

end
