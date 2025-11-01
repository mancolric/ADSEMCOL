include("test_ConvectionDiffusion.jl")

#Exponential boundary layer in space, tanh in time.

function ExponentialBoundaryLayer(hp0::Float64, FesOrder::Int;
    deltax::Float64=2e-2, deltat::Float64=0.2, epsilon::Float64=deltax,
    tf::Float64=1.0, RKMethod::String="Ascher3",
    #
    SpaceAdapt::Bool=true, SpaceNorm::Float64=2.0, TolS::Float64=1e-5, 
    AMA_MaxIter::Int=500, AMA_SizeOrder::Int=FesOrder, 
    AMA_AnisoOrder::Int=2, AMA_ProjN::Int=1, AMA_ProjOrder::Int=0,
    #
    TimeAdapt::Bool=true, TolT::Float64=1e-7, Deltat0::Float64=1e-4, 
    #
    PlotFig::Bool=false, wFig::Float64=9.50, hFig::Float64=6.50, 
    SaveFig::Bool=false, Nt_SaveFig::Int=typemax(Int), Deltat_SaveFig::Float64=0.01,
    #
    SaveRes::Bool=true, Nt_SaveRes::Int=typemax(Int), Deltat_SaveRes::Float64=Inf, 
    #
    SC::Int=0)
    
    #---------------------------------------------------------------------
    #PROBLEM DATA:
    
    #NOTE: The argument of the functions can be AbstractArray, and hence the output 
    #is a type-unstable AbstractArray. However, when FunctionWrapper is used, the 
    #argument is going to be <:AMF64 and the output a Matrix{Float64}.
    
    #Theoretical solution:
    function tfun(t::Float64)
#         return 1.0
        return 0.5+0.5*tanh((t-0.5)/deltat)
    end
    function dtfun(t::Float64)
#         return 0.0
        return 0.5/cosh((t-0.5)/deltat)^2 * 1.0/deltat
    end

    #In the steady state
    #   u   = x - (exp((x-1)/delta)-exp(-1/delta))/(1.0-exp(-1/delta))
    #       = x + c1*exp((x-1)/delta) + c2 
    c1      = -1.0/(1.0-exp(-1.0/deltax))
    c2      = exp(-1.0/deltax)/(1.0-exp(-1.0/deltax))
    function utheor(t::Float64, x::Vector{<:AAF64})
    
        #=
        #Try this with hp0=0.1, FesOrder=3, AnisoOrder=2 to plot G_K and mathcal{G}_K:
        u   = @. 0.5*x[2]^6 + 2.0*x[2]^1*x[1]^5 - 1.5*x[2]^3*x[1]^3
        return [u]
        =#
        
        u   = @avxt @. tfun(t) * 
                    (x[1]+c1*exp((x[1]-1.0)/deltax)+c2) *
                    (x[2]+c1*exp((x[2]-1.0)/deltax)+c2) + 1.0
        return [u]
        
    end
    function dk_utheor(t::Float64, x::Vector{<:AAF64})
    
        u       = @avxt @. tfun(t) * 
                    (x[1]+c1*exp((x[1]-1.0)/deltax)+c2) *
                    (x[2]+c1*exp((x[2]-1.0)/deltax)+c2) + 1.0
        du_dt   = @avxt @. dtfun(t) * 
                    (x[1]+c1*exp((x[1]-1.0)/deltax)+c2) *
                    (x[2]+c1*exp((x[2]-1.0)/deltax)+c2)
        du_dx   = @avxt @. tfun(t) * 
                    (1.0+c1*exp((x[1]-1.0)/deltax)/deltax) *
                    (x[2]+c1*exp((x[2]-1.0)/deltax)+c2)
        du_dy   = @avxt @. tfun(t) * 
                    (x[1]+c1*exp((x[1]-1.0)/deltax)+c2) *
                    (1.0+c1*exp((x[2]-1.0)/deltax)/deltax)
        d2u_dx2 = @avxt @. tfun(t) * 
                    (c1*exp((x[1]-1.0)/deltax)/deltax^2) *
                    (x[2]+c1*exp((x[2]-1.0)/deltax)+c2)
        d2u_dy2 = @avxt @. tfun(t) * 
                    (x[1]+c1*exp((x[1]-1.0)/deltax)+c2) *
                    (c1*exp((x[2]-1.0)/deltax)/deltax^2)
        
        return (du_dt, du_dx, du_dy, d2u_dx2, d2u_dy2)
        
    end
    uplot_min   = 1.0
    uplot_max   = 2.0
    
    #Velocity:
    function afun(t::Float64, x::Vector{<:AAF64}, u::Vector{<:AAF64})
        
        a1  = @avxt @. 0.5*u[1]+0.0
        a2  = @avxt @. 0.5*u[1]+0.0

        return [a1, a2]
    end
    function dafun(t::Float64, x::Vector{<:AAF64}, u::Vector{<:AAF64})
        
        da1  = @avxt @. 0.0*u[1]+0.5
        da2  = @avxt @. 0.0*u[1]+0.5

        return [da1, da2]
    end
    
    #Viscosity:
    function epsilonfun(t::Float64, x::Vector{<:AAF64}, u::Vector{<:AAF64})
    
        epsilon0    = @avxt @. 0.0*x[1]+epsilon
        return [ epsilon0 ]
        
    end
    function depsilonfun(t::Float64, x::Vector{<:AAF64}, u::Vector{<:AAF64})
    
        depsilon0   = @avxt @. 0.0*x[1]+0.0
        return [ depsilon0 ]
        
    end
    
    #Reaction terms:
    function Qfun(t::Float64, x::Vector{<:AAF64}, u::Vector{<:AAF64})
    
        uth                                     = utheor(t,x)[1]
        du_dt, du_dx, du_dy, d2u_dx2, d2u_dy2   = dk_utheor(t,x)
        
        #Q = u_t + a*u_x - epsilon*u_xx
#         Q0          = @avxt @. du_dt + 1.0*(du_dx + du_dy) - epsilon*(d2u_dx2+d2u_dy2)
        Q0          = @avxt @. du_dt + uth*(du_dx + du_dy) - epsilon*(d2u_dx2+d2u_dy2)
        return [ Q0 ]
        
    end
    function dQfun(t::Float64, x::Vector{<:AAF64}, u::Vector{<:AAF64})
    
       return [ 0.0*u[1] ]
        
    end
    
    #Boundary conditions:
    function uDir(t::Float64, x::Vector{<:AAF64})
        return utheor(t,x)
    end
    function qNeum2(t::Float64, x::Vector{<:AAF64})
        #Condition at boundary 2:
        # q= -epsilon*du/dx
        dku         = dk_utheor(t,x)
        return [-epsilon*dku[2]]
    end
    function qNeum3(t::Float64, x::Vector{<:AAF64})
        #Condition at boundary 3:
        # q= -epsilon*du/dy
        dku         = dk_utheor(t,x)
        return [-epsilon*dku[3]]
    end
    BC_Dir          = Dirichlet(FWt11((t,x)->uDir(t,x)))
    BC_Neum2        = Neumann(FWt11((t,x)->qNeum2(t,x)))
    BC_Neum3        = Neumann(FWt11((t,x)->qNeum3(t,x)))
    
    #Structure with convection--diffusion problem data:
    ProblemData     = SCD(  FWt21((t,x,u)->afun(t,x,u)), 
                            FWt21((t,x,u)->epsilonfun(t,x,u)),
                            FWt21((t,x,u)->Qfun(t,x,u)), 
                            FWt21((t,x,u)->dafun(t,x,u)), 
                            FWt21((t,x,u)->depsilonfun(t,x,u)),
                            FWt21((t,x,u)->dQfun(t,x,u)) )
    
    #Mesh:
    MeshFile            = "$(@__DIR__)/../../temp/ExponentialBoundaryLayer_SC$(SC).geo"
    NX                  = Int(ceil(2.0/(hp0*FesOrder)))
    TrMesh_Rectangle_Create!(MeshFile, 0.0, 1.0, NX, 0.0, 1.0, NX)
    
    #----------------------------------------------------------
    #INITIAL CONDITION:
    
    t_ini               = time()
    
    #Create solver structure:
    solver              = LIRKHyp_Start(ProblemData)
    solver.ProblemName  = "ExponentialBoundaryLayer"
    solver.SC           = SC
    solver.MeshFile     = MeshFile
    solver.nBounds      = 4
    solver.FesOrder     = FesOrder
    solver.Deltat0      = Deltat0
    solver.tf           = tf
    solver.SpaceAdapt   = SpaceAdapt
    solver.TolS_max     = TolS
    solver.TolS_min     = 0.0*TolS
    solver.SpaceNorm    = SpaceNorm
    solver.TimeAdapt    = TimeAdapt
    solver.TolT         = TolT
    solver.AMA_MaxIter  = AMA_MaxIter
    solver.AMA_SizeOrder    = AMA_SizeOrder
    solver.AMA_AnisoOrder   = AMA_AnisoOrder
    solver.AMA_ProjN        = AMA_ProjN
    solver.AMA_ProjOrder    = AMA_ProjOrder
#     solver.AMA_rDeltat      = 1.0
    solver.TolA_min         = 1e-12
    solver.RKMethod         = RKMethod
    
    #Set initial and boundary conditions:
    solver.u0fun        = FW11((x) -> utheor(0.0,x)) 
#     solver.BC           = [ BCW(BC_Dir), BCW(BC_Neum2), BCW(BC_Dir), BCW(BC_Dir) ]
    solver.BC           = [ BCW(BC_Dir), BCW(BC_Dir), BCW(BC_Dir), BCW(BC_Dir) ]
    
    #Initial condition:
    ConvFlag            = LIRKHyp_InitialCondition!(solver)
    tv                  = [0.0]
    
    #Compute Lq error:
    errLq,              = LqError(solver, FW11((x) -> utheor(solver.t, x)), q=SpaceNorm)
    hmean               = 2.0*sqrt(solver.Omega/solver.mesh.nElems/TrElem_Area)
    errL2L2             = errLq
    etaL2L2             = solver.etaS
    uP                  = GetPolynomials(solver.u, solver.fes)
    umin                = minimum(uP[1])
    umax                = maximum(uP[1])
    println("h=", sprintf1("%.2e",hmean), 
            ", eta_S=", sprintf1("%.2e",solver.etaS),
            ", e_S=", sprintf1("%.2e",errLq),
            ", eta_L2L2=", sprintf1("%.2e",etaL2L2),
            ", e_L2L2=", sprintf1("%.2e",errL2L2),
            ", umin=", sprintf1("%.2e",umin),
            ", umax=", sprintf1("%.2e",umax),
            ", umax-umin=", sprintf1("%.2e",umax-umin),
            ", time=", time()-t_ini)
                
    #Function to plot solution:
    figv                = Vector{Figure}(undef,4)
    if PlotFig
        for ii=1:length(figv)
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
            PyPlot.cla()
            PlotContour(solver.u[1], solver.fes, zv=linspace(uplot_min, uplot_max, 100))
#             PlotContour(solver.u[1].-utheor(solver.t,[zeros(1), zeros(1)])[1], solver.fes)
#             PlotMesh!(solver.mesh, color="k")
            axis("equal")
            if SaveFig
                savefig("$(VideosUbi)SC$(SC)_Contour_$(nb_SaveFig).png", dpi=800, pad_inches=0)
            end
            
            figure(figv[2].number)
            PyPlot.cla()
            PlotContour(solver.u[1], solver.fes, zv=linspace(uplot_min, uplot_max, 100), cmap="jet")
            PlotMesh!(solver.mesh, color="w")
            axis("equal")
            if SaveFig
                savefig("$(VideosUbi)SC$(SC)_Mesh_$(nb_SaveFig).png", dpi=800, pad_inches=0)
            end
            
            figure(figv[3].number)
            PyPlot.cla()
            semilogy(solver.tv, solver.etaSv, ".-b")
            semilogy(solver.tv, solver.etaTv, ".-g")
            semilogy(solver.tv, solver.etaAv, ".-r")
#             semilogy(solver.tv, solver.nElemsv, ".-k")
        
            figure(figv[4].number)
            PyPlot.cla()
            plot(solver.tv, solver.nDofv, "k")
            
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
            save("$(ResUbi)LIRKHyp_SC$(SC)_$(nb_SaveRes).jld2", "StudyCase", "ExponentialBoundaryLayer",
                "deltax", deltax, "deltat", deltat, 
                "epsilon", epsilon, 
                "ConvFlag", ConvFlag, "solver", save(solver), 
                "tv", tv,
                "errL2L2", errL2L2, "etaL2L2", etaL2L2 )
            t_lastRes   += Deltat_SaveRes
            ct_SaveRes  = 0
            nb_SaveRes  += 1
        end
        return
        
    end
    SaveSol()
    
    #----------------------------------------------------------
    #MARCH IN TIME:
    
    nIter               = 0
    while solver.t<solver.tf
    
        ConvFlag            = LIRKHyp_Step!(solver)
#         ConvFlag            = IRK_Step!(solver)
        if ConvFlag<=0
            break
        end
        push!(tv, solver.t)
        
        #Compute Lq error:
        errLq,              = LqError(solver, FW11((x) -> utheor(solver.t, x)), q=SpaceNorm)
        hmean               = 2.0*sqrt(solver.Omega/solver.mesh.nElems/TrElem_Area)
        errL2L2             = errL2L2_(solver.t-solver.Deltat, errL2L2, solver.t, errLq)
        etaL2L2             = errL2L2_(solver.t-solver.Deltat, etaL2L2, solver.t, solver.etaS+solver.etaT)
        uP                  = GetPolynomials(solver.u, solver.fes)
        umin                = minimum(uP[1])
        umax                = maximum(uP[1])
    
        #Print info:
        println("h=", sprintf1("%.2e",hmean), 
                ", eta_S=", sprintf1("%.2e",solver.etaS),
                ", e_S=", sprintf1("%.2e",errLq),
                ", eta_L2L2=", sprintf1("%.2e",etaL2L2),
                ", e_L2L2=", sprintf1("%.2e",errL2L2),
                ", umin=", sprintf1("%.2e",umin),
                ", umax=", sprintf1("%.2e",umax),
                ", umax-umin=", sprintf1("%.2e",umax-umin),
                ", time=", time()-t_ini)
        
        #Save figures and/or results?
        PlotSol()
        SaveSol()
        
        if errLq>0.1
            display(solver.t-solver.Deltat)
            display(solver.t-solver.Deltat .+ (solver.RK.c*solver.Deltat))
            error("too large error")
        end
        
    end

    return solver
    
end

