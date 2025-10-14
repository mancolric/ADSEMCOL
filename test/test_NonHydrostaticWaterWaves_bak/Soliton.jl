cd(@__DIR__)
include("test_NonHydrostaticWaterWaves.jl")

function Soliton(hp0::Float64, FesOrder::Int;
    tf::Float64=10.0, RKMethod::String="BPR3",
    epsilon::Float64=0e-3, alpha::Float64=5.0, g::Float64=9.8, A::Float64=0.2, 
    h0::Float64=1.0, gamma::Float64=2.0, xend::Float64=100.0, 
    #
    TolS::Float64=1e-4, AMA_MaxIter::Int=200, AMA_SizeOrder::Int=FesOrder, AMA_AnisoOrder::Int=2,
    SpaceAdapt::Bool=true, 
    #
    TolT::Float64=1e-4, Deltat0::Float64=1e-3, TimeAdapt::Bool=true,
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
    model               = NHWW()
    model.epsilon       = epsilon
    model.gamma         = gamma
    model.g             = g
    c                   = alpha*sqrt(g*h0)
    model.c             = c
    model.h0            = h0
    model.CSS           = CSS

#     function utheorfun(t::Float64, x::Vector{Matrix{Float64}}) 
# 
#         h, q1, q2, q3, p        = SolitonExact(t, x[1], A=A, gamma=gamma, 
#                                                 h0=h0, x0=0.0, g=g) + 
#                                     SolitonExact(t, x[2], A=A, gamma=gamma, 
#                                                 h0=h0, x0=0.0, g=g, theta=pi/2)
#         b                       = @. sin(x[1])+cos(x[2])
#         eta                     = @. h+b
#         P                       = @. h*(p-c*c*log(h/h0))
#         return [eta, q1, q2, q3, P, b]
# 
#     end
    function utheorfun(t::Float64, x::Vector{Matrix{Float64}}) 

        h, q1, q2, q3, p        = SolitonExact(t, x[1], A=A, gamma=gamma, 
                                                h0=h0, x0=0.0, g=g)
        b                       = @. 0.0*x[1]
        eta                     = h
        P                       = @. h*(p-c*c*log(h/h0))
        return [eta, q1, q2, q3, P, b]

    end
    
    function u0fun(x::Vector{Matrix{Float64}})

        return utheorfun(0.0, x)

    end

    #Boundary conditions:
    BC_horiz        = SlipAdiabatic()
    function uLeft(t::Float64, x::Vector{Matrix{Float64}})
        q1              = @tturbo @. 0.0*x[1]
        q2              = @tturbo @. 0.0*x[1]
        q3              = @tturbo @. 0.0*x[1]
        p               = @tturbo @. 0.0*x[1]
        return [q1, q2, q3, p]

    end
    BC_left         = SubsonicInlet1(FWt11((t,x)->uLeft(t,x)))
    function uRight(t::Float64, x::Vector{Matrix{Float64}})
        h               = @tturbo @. 0.0*x[1] + h0
        eta             = h
        p               = @tturbo @. 0.0*x[1]
        return [eta, p]
    end
    BC_right        = SubsonicOutlet1(FWt11((t,x)->uRight(t,x)))
#     BC_right        = DoNothing1()
    
    #---------------------------------------------------------------------
    #PRE-PROCESS STAGE:

    #Mesh:
    MeshFile                = "../temp/Soliton$(SC).geo"
    NX                      = Int(ceil((xend+20.0)/(hp0*FesOrder)))
    NY                      = 2
#     TrMesh_Rectangle_Create!(MeshFile, -5.0, 5.0, NX, -5.0, 5.0, NY)
    TrMesh_Rectangle_Create!(MeshFile, -20.0, xend, NX, -5.0, 5.0, NY)

    #Load LIRKHyp solver structure with default data. Modify the default data if necessary:
    solver                  = LIRKHyp_Start(model)
    solver.ProblemName      = "Soliton" 
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
    ConvFlag            = LIRKHyp_InitialCondition!(solver, AMA_RefineFactor=0.9, 
                            DEq_MaxIter=0)
#     CheckJacobian(solver, Plot_df_du=true, Plot_df_dgradu=true)
#     CheckJacobian(solver, Plot_dQ_du=true, Plot_dQ_dgradu=true)
#     for ii = 4
#         BC_CheckJacobian(solver, ii, Plot_df_du=true, Plot_df_dgradu=true)
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
                savefig("$(VideosUbi)Soliton$(SC)_$(nb_SaveFig).png", dpi=400, pad_inches=0)
            end
            
            figure(figv[2].number)
            #Domain limits:
            c0      = sqrt(g*(A+h0))
            l0      = h0*sqrt((A+h0)/h0)
            x11     = c0*solver.t-10*l0
            x12     = c0*solver.t+10*l0
            #Interpolate exact solution:
            x1v     = linspace(x11, x12, 500)
            x2v     = zeros(size(x1v))
            xm      = [ x1v, x2v ]
            um      = utheorfun(solver.t, Vector{<:AMF64}(xm))
            vm      = DepVars(solver.model, solver.t, Vector{<:AMF64}(xm), um, PlotVars)
            #Loop plot variables:
            for ii=1:length(PlotVars)
                PyPlot.subplot(mFig, nFig, ii)
                PyPlot.subplots_adjust(hspace=0.8)
                PyPlot.cla()
                splot_fun(x1,x2)    = @mlv x1
                PlotNodes(splot_fun, solver, PlotVars[ii])
                plot(x1v, vm[ii][1], "r", linewidth=0.5)
                xlabel(latexstring("x_1"), fontsize=10)
                title(latexstring(LatexString(PlotVars[ii]),
                                  "; t^n=", sprintf1("%.2e", solver.t)),
                        fontsize=10)
                xlim([x11, x12])
                
            end
            if SaveFig
                savefig("$(VideosUbi)SolitonNodes$(SC)_$(nb_SaveFig).png", dpi=400, pad_inches=0)
            end

            #=
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
            end
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
            save("$(ResUbi)LIRKHyp_SC$(SC)_$(nb_SaveRes).jld2", "StudyCase", "Soliton",
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
        save("$(ResUbi)LIRKHyp_SC$(SC)_1000.jld2", "StudyCase", "Soliton",
            "ConvFlag", ConvFlag, "solver", save(solver), "hmean", hmean, "e_Lq", errLq, "Deltat_mean", Deltat_mean)
    end

    return solver

end

function SolitonExact(t::Float64, x::AMF64; A::Float64=0.2, gamma::Float64=2.0, 
    h0::Float64=1.0, x0::Float64=0.0, g::Float64=9.8, theta::Float64=0.0)

    #Soliton velocity and length:
    c0      = sqrt(g*(A+h0))
    l0      = h0*sqrt((A+h0)/(A*gamma/2))
    z       = @. (x - (x0 + c0*t))/l0
    phi     = @. sech(z)
    dphi    = @. -sinh(z)/cosh(z)^2
    d2phi   = @. - 1/cosh(z) + 2*sinh(z)^2/cosh(z)^3
    
    #Solution:
    h       = @. h0 + A*phi^2
    q       = @. c0*(h-h0)
    q1      = @. q*cos(theta)
    q2      = @. q*sin(theta)
    q3      = @. -A*c0*h0/l0*phi*dphi
    p       = @. A*c0^2*h0^2/(2*l0^2*h^2)*( (2*h0-h)*dphi^2 + h*phi*d2phi )
    eta     = h
    
    return [ eta, q1, q2, q3, p ]
    
end

function SolitonRelaxed_old(; A::Float64=0.2, gamma::Float64=2.0, 
    h0::Float64=1.0, href::Float64=h0, g::Float64=9.8, c::Float64=10*sqrt(g*h0), 
    Deltaxi::Float64=1e-6, xif::Float64=10.0, theta::Float64=0.0)

    #Initial and final values for xi:
    xi0         = 0.0       #symmetry is applied later
    
    #Soliton velocity and length:
    c0          = sqrt(g*(A+h0))
#     l0          = h0*sqrt((A+h0)/h0)
    l0          = h0*sqrt((A+h0)/(A*gamma/2))
    
    #Load coefficients:
    RKCoeffs    = RK_Coefficients("BPR3")
    Bm          = RKCoeffs.AE
    ABm         = RKCoeffs.AI-RKCoeffs.AE
    cv          = RKCoeffs.c
    ss          = RKCoeffs.stages
    
    #Initial condition for u=[h,w]:
    un0         = [ h0+A, 0.0 ]
    
    #Function to march in time:
    function SolitonFun(u::Vector{Float64}, ComputeJ::Bool)
        h           = u[1]
        w           = u[2]
        p           = (c0^2*(h-h0)*h0/h - 0.5*g*(h^2-h0^2))/h
        hprime      = -2*l0*c^2*h/(c0*h0)*(c^2 + c0^2*h0/h + 
                            - 2*c0*h0^2/h^2 + 0.5*g*h0^2/h + 0.5*g*h)^(-1) * w
        wprime      = -2*l0/(c0*h0)*p
        f           = [hprime, wprime]
#         f           = [ 1.0, 1.0 ]
        J           = zeros(0,0)
        if ComputeJ
            J       = zeros(2,2)
            delta   = 1e-6
            for jj=1:size(J,2)
                upert       = copy(u)
                upert[jj]   += delta
                fpert       = SolitonFun(upert, false)[1]
                J[:,jj]     = (fpert-f)/(delta)
            end
        end
        return f, J
    end
    
    #Loop:
    xin         = xi0
    un          = un0
    xiv         = [xi0]
    um          = permutedims(un)
    fm          = zeros(length(un), ss)
    Jum         = zeros(length(un), ss)
    while xin<=xif
    
        #First stage:
        unp1            = copy(un)
        xinp1           = xin
        fm[:,1],Jn      = SolitonFun(un, true)
        Jum[:,1]        = Jn*un
        LSm             = eye(2) - Deltaxi*ABm[2,2]*Jn
        LSm_F           = factorize(LSm)
        
        #Loop stages:
        for kk=2:ss
        
            xinp1       = xin + cv[kk]*Deltaxi
            zn          = 1.0*un
            for ll=1:kk-1
                zn      .+= Deltaxi*(Bm[kk,ll]*fm[:,ll] + ABm[kk,ll]*Jum[:,ll])
            end
            unp1        .= LSm_F\zn
            fm[:,kk]    = SolitonFun(unp1, false)[1]
            Jum[:,kk]   = Jn*unp1
            
        end
        
        #Update solution:
        xin         = xinp1
        un          .= unp1
        xiv         = vcat(xiv, xin)
        um          = vcat(um, permutedims(un))
        
        println("xi=", xin)
        
    end
    
    #Extract solution:
    aux     = length(xiv):-1:1
    xiv     = vcat(-xiv[aux], xiv)
    hv      = vcat(um[aux,1], um[:,1])
    wv      = vcat(-um[aux,2], um[:,2])
    pv      = @. (c0^2*(hv-h0)*h0/hv - 0.5*g*(hv^2-h0^2))/hv
    Pv      = @. hv*(pv-c^2*log(abs(hv./href)))
    qv      = @. c0*(hv-h0)
    q1v     = @. qv*cos(theta)
    q2v     = @. qv*sin(theta)
    q3v     = @. hv*wv
    etav    = hv
    
    #Exact solution:
    zexact  = SolitonExact(0.0, l0*reshape(xiv,:,1), A=A, gamma=gamma, 
                                h0=h0, g=g, theta=theta)
    
    figure()
    plot(xiv, etav)
    plot(xiv, zexact[1], "r")
    xlabel("xi")
    ylabel("eta")
    
    figure()
    plot(xiv, q1v./hv)
    plot(xiv, zexact[2]./zexact[1], "r")
    xlabel("xi")
    ylabel("w")
    
    figure()
    plot(xiv, wv)
    plot(xiv, zexact[4]./zexact[1], "r")
    xlabel("xi")
    ylabel("w")
    
    figure()
    plot(xiv, pv)
    plot(xiv, zexact[5], "r")
    xlabel("xi")
    ylabel("p")
    
    save("Soliton.jld2", "xiv", xiv, "etav", etav, "q1v", q1v, "q2v", q2v, 
        "q3v", q3v, "Pv", Pv)
        
    return xiv, etav, q1v, q2v, q3v, Pv
            
end    

function SolitonRelaxed(; A::Float64=0.2, gamma::Float64=2.0, 
    h0::Float64=1.0, g::Float64=9.8, c::Float64=10*sqrt(g*h0), 
    Deltaxi::Float64=1e-6, xif::Float64=10.0, theta::Float64=0.0)

    #Define model:
    model               = NHWW()
    model.g             = g
    model.c             = c
    model.gamma         = gamma
    model.h0            = h0
    
    #Initial and final values for xi:
    xi0         = 0.0       #symmetry is applied later
    
    #Soliton velocity and length:
    c0          = sqrt(g*(A+h0))
    l0          = h0*sqrt((A+h0)/(A*gamma/2))
    
    #Load coefficients:
    RKCoeffs    = RK_Coefficients("BPR3")
    Bm          = RKCoeffs.AE
    ABm         = RKCoeffs.AI-RKCoeffs.AE
    cv          = RKCoeffs.c
    ss          = RKCoeffs.stages
    
    #Initial condition for u=[eta, q1, q2, q3, P, b] (q2=0):
    un0         = [ h0, 0.0, 0.0, 0.0, h0*(0.0-c*c*log(h0/h0))+1e-5, 0.0 ]
#     un0         = SolitonRelaxedPeak(A=A, gamma=gamma, h0=h0, g=g, c=c, theta=theta)
    
    #Function to march in time:
    function SolitonFun(u::Vector{Float64}, ComputeJ::Bool)
        
        #Reshape solution:
        usol            = Vector{Matrix{Float64}}(undef, model.nVars)
        dusol_dx        = Matrix{Matrix{Float64}}(undef, model.nVars, 2)
        for II=1:model.nVars
            usol[II]    = reshape([u[II]],1,1)
        end
        alloc!(dusol_dx, (1,1))     #we set du_dx=0 because this term is only needed to compute db/dx in source term
        
        #Compute fluxes and Jacobian:
        fsol            = Matrix{Matrix{Float64}}(undef, model.nVars, 2)
        dfsol_dusol     = Array{Matrix{Float64},3}(undef, model.nVars, 2, model.nVars)
        alloc!(fsol, (1,1))
        alloc!(dfsol_dusol, (1,1))
        HyperbolicFlux!(model, usol, true, fsol, dfsol_dusol)
        
        #Compute source terms:
        Qsol            = Vector{Matrix{Float64}}(undef, model.nVars)
        dQsol_dusol     = Matrix{Matrix{Float64}}(undef, model.nVars, model.nVars)
        dQsol_dusol_dx  = Array{Matrix{Float64},3}(undef, model.nVars, model.nVars, 2)
        alloc!(Qsol, (1,1))
        Source!(model, [zeros(1,1), zeros(1,1)], Inf, usol, dusol_dx, false, 
                Qsol, dQsol_dusol, dQsol_dusol_dx)
                
        #Extract terms:
        df_du           = Matrix{Float64}(undef, model.nVars, model.nVars)
        Q               = Vector{Float64}(undef, model.nVars)
        for II=1:model.nVars, JJ=1:model.nVars
            df_du[II,JJ]    = dfsol_dusol[II,1,JJ][1,1]
        end
        for II=1:model.nVars
            Q[II]           = Qsol[II][1,1]
        end
            
        #Solve system:
        du_dxi          = l0*((df_du-c0*eye(model.nVars))\Q)
        J               = zeros(0,0)
        if ComputeJ
            J       = zeros(model.nVars, model.nVars)
            delta   = 1e-6
            for jj=1:size(J,2)
                upert       = copy(u)
                upert[jj]   += delta
                fpert       = SolitonFun(upert, false)[1]
                J[:,jj]     = (fpert-du_dxi)/(delta)
            end
        end
        return du_dxi, J
        
    end
    
    #Loop:
    xin         = xi0
    un          = un0
    xiv         = [xi0]
    um          = permutedims(un)
    fm          = zeros(length(un), ss)
    Jum         = zeros(length(un), ss)
    while xin<=xif
    
        #First stage:
        unp1            = copy(un)
        xinp1           = xin
        fm[:,1],Jn      = SolitonFun(un, true)
        Jum[:,1]        = Jn*un
        LSm             = eye(model.nVars) - Deltaxi*ABm[2,2]*Jn
        LSm_F           = factorize(LSm)
        
        #Loop stages:
        for kk=2:ss
        
            xinp1       = xin + cv[kk]*Deltaxi
            zn          = 1.0*un
            for ll=1:kk-1
                zn      .+= Deltaxi*(Bm[kk,ll]*fm[:,ll] + ABm[kk,ll]*Jum[:,ll])
            end
            unp1        .= LSm_F\zn
            fm[:,kk]    = SolitonFun(unp1, false)[1]
            Jum[:,kk]   = Jn*unp1
            
        end
        
        #Exit if h<0
        if unp1[1]<0.0
            break
        end
        
        #Update solution:
        xin         = xinp1
        un          .= unp1
        xiv         = vcat(xiv, xin)
        um          = vcat(um, permutedims(un))
        
        println("xi=", xin)
        
    end
    
    #Extract solution:
#     aux     = length(xiv):-1:1
    aux     = zeros(Int,0)
    xiv     = vcat(xiv[aux], xiv)
    hv      = vcat(um[aux,1], um[:,1])
    qv      = vcat(um[aux,2], um[:,2])
    uv      = @. qv/hv
    q1v     = @. qv*cos(theta)
    q2v     = @. qv*sin(theta)
    q3v     = vcat(um[aux,4], um[:,4])
    wv      = @. q3v/hv
    Pv      = vcat(um[aux,5], um[:,5])
    pv      = @. Pv/hv + c*c*log(hv/h0)
    etav    = hv    
    bv      = vcat(um[aux,6], um[:,6])
    
    #Compute peak:
    _,ipeak  = findmax(etav)
    xiv     = xiv .- xiv[ipeak]
    
    #Compute exact solution:
    zexact  = SolitonExact(0.0, reshape(l0*xiv,:,1), A=A, gamma=gamma, h0=h0, g=g, theta=theta)
    
    figure()
    plot(xiv, etav, "b")
    plot(xiv, zexact[1], "r")
    xlabel(latexstring("\\xi"))
    ylabel(latexstring("\\eta"), rotation=0)
#     plot(xiv[ipeak], zpeak[1], "ok")
    legend(["Relaxed", "Exact"])
    
    figure()
#     plot(xiv, q1v, "b")
#     plot(xiv, zexact[2], "r")
    plot(xiv, uv, "b")
    plot(xiv, zexact[2]./zexact[1], "r")
    xlabel(latexstring("\\xi"))
#     ylabel(latexstring("q_1"), rotation=0)
    ylabel(latexstring("u_1"), rotation=0)
#     plot(xiv[ipeak], zpeak[2], "ok")
    
    figure()
#     plot(xiv, q3v, "b")
#     plot(xiv, zexact[4], "r")
    plot(xiv, wv, "b")
    plot(xiv, zexact[4]./zexact[1], "r")
    xlabel(latexstring("\\xi"))
#     ylabel(latexstring("q_3"), rotation=0)
    ylabel(latexstring("u_3"), rotation=0)
#     plot(xiv[ipeak], zpeak[4], "ok")
    
#     figure()
#     plot(xiv, Pv)
#     xlabel(latexstring("\\xi"))
#     ylabel(latexstring("P"), rotation=0)
#     #plot(xiv[ipeak], zpeak[5], "ok")
    
    figure()
    plot(xiv, pv, "-xb")
    plot(xiv, zexact[5], "r")
    pv2     = @. c0^2*(zexact[1]-h0)*h0/zexact[1]^2 + 0.5*g*(h0^2-zexact[1]^2)/zexact[1]
    plot(xiv, pv2, "c")
    xlabel(latexstring("\\xi"))
    ylabel(latexstring("p"), rotation=0)
    legend(["Relaxed", "Exact", "Ec. (30) JCP"])
    
    figure()
    plot(xiv, Pv, "b")
    plot(xiv, zexact[1].*(zexact[5]-c*c*log.(zexact[1]/h0)), "r")
    xlabel(latexstring("\\xi"))
    ylabel(latexstring("P"), rotation=0)
    legend(["Relaxed", "Exact"])
    
#     figure()
#     plot(xiv, zexact[5]./pv, "r")
#     xlabel(latexstring("\\xi"))
#     ylabel(latexstring("p"), rotation=0)
    
#     save("Soliton.jld2", "xiv#=", xiv, "etav", etav, "q1v", q1v, "q2v", q2v, 
#         "q3v", q3v, "Pv", Pv)
        
    return xiv, etav, q1v, q2v, q3v, Pv
            
end    

