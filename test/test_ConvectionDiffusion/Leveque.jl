include("test_ConvectionDiffusion.jl")

function Leveque(hp0::Float64, FesOrder::Int;
    delta::Float64=1e-3,
    tf::Float64=1.0,
    #
    SpaceAdapt::Bool=true, TolS::Float64=1e-4, 
    AMA_MaxIter::Int=200, AMA_ProjN::Int=1, AMA_ProjOrder::Int=0,
    SpaceNorm::Float64=2.0, 
    PostAdapt::Bool=false, 
    #
    TimeAdapt::Bool=true, TolT::Float64=1e-4, Deltat0::Float64=1e-4, 
    #
    CA::Float64=0.01, 
    #
    PlotFig::Bool=false, wFig::Float64=9.50, hFig::Float64=6.50, 
    SaveFig::Bool=false, Nt_SaveFig::Int=typemax(Int), Deltat_SaveFig::Float64=0.05,
    #
    SaveRes::Bool=true, Nt_SaveRes::Int=typemax(Int), Deltat_SaveRes::Float64=Inf, 
    #
    SC::Int=0)
    
    #---------------------------------------------------------------------
    #PROBLEM DATA:
    
    #NOTE: The argument of the functions can be AbstractArray, and hence the output 
    #is a type-unstable AbstractArray. However, when FunctionWrapper is used, the 
    #argument is going to be <:AMF64 and the output a Matrix{Float64}.
    
    #Velocity:
    function afun(t::Float64, x::Vector{<:AAF64}, u::Vector{<:AAF64})
        
        a1  = @avxt @. -2.0*pi*x[2]
        a2  = @avxt @. 2.0*pi*x[1]
        return [a1, a2]
    end
    function dafun(t::Float64, x::Vector{<:AAF64}, u::Vector{<:AAF64})
        
        da1     = @avxt @. 0.0*x[2]
        da2     = @avxt @. 0.0*x[1]
        return [da1, da2]
    end
    
    #Viscosity:
    function epsilonfun(t::Float64, x::Vector{<:AAF64}, u::Vector{<:AAF64})
    
        epsilon0    = @mlv 0.0*x[1]+0.0
        return [ epsilon0 ]
        
    end
    function depsilonfun(t::Float64, x::Vector{<:AAF64}, u::Vector{<:AAF64})
    
        depsilon0   = @avxt @. 0.0*x[1]+0.0
        return [ depsilon0 ]
        
    end
    
    #Reaction terms:
    function Qfun(t::Float64, x::Vector{<:AAF64}, u::Vector{<:AAF64})
    
        #Q = u_t + a*u_x - epsilon*u_xx
        Q0          = @mlv  0.0*x[1]
        return [ Q0 ]
        
    end
    function dQfun(t::Float64, x::Vector{<:AAF64}, u::Vector{<:AAF64})
    
       return [ 0.0*u[1] ]
        
    end
    
    #Initial condition:
    function u0fun(x::Vector{<:AAF64})
        
        rcyl    = @avxt @. sqrt((x[1]-0.5)^2 + x[2]^2)
        rhump   = @avxt @. sqrt((x[1]+0.5)^2 + x[2]^2)
        rcone   = @avxt @. sqrt(x[1]^2 + (x[2]+0.5)^2)
        
        ucyl    = @avxt @. SmoothHeaviside(rcyl-0.25, delta, 1.0, 0.0) * #1 if inside the circle
                            (1.0 - 
                                SmoothHeaviside(abs(x[2])-0.05, delta, 1.0, 0.0)*
                                SmoothHeaviside(x[1]-0.4, delta, 0.0, 1.0) )
                                
        uhump   = @avxt @. cos(1.5*pi*rhump)^3*(rhump<1/3)
        
        ucone   = @avxt @. 0.0*(1.0 - rcone/0.25)*(rcone<0.25)
        
        return [ ucyl + uhump + ucone ]
        
    end
    
    #Theoretical solution:
    function utheor(t::Float64, x::Vector{<:AAF64})
    
        #Rotation matrix:
        omega       = 2*pi
        Q           = [ cos(omega*t)    sin(omega*t)
                        -sin(omega*t)   cos(omega*t) ]
        x0          = [ Q[1,1]*x[1]+Q[1,2]*x[2], 
                        Q[2,1]*x[1]+Q[2,2]*x[2] ]
        return u0fun(x0)
    end
    
    #Boundary conditions:
    function uDir(t::Float64, x::Vector{<:AAF64})
        return [0.0*x[1].+0.0]
    end
    BC_Dir          = Dirichlet(FWt11((t,x)->uDir(t,x)))
    
    #Structure with convection--diffusion problem data:
    ProblemData     = SCD(  FWt21((t,x,u)->afun(t,x,u)), 
                            FWt21((t,x,u)->epsilonfun(t,x,u)),
                            FWt21((t,x,u)->Qfun(t,x,u)), 
                            FWt21((t,x,u)->dafun(t,x,u)), 
                            FWt21((t,x,u)->depsilonfun(t,x,u)),
                            FWt21((t,x,u)->dQfun(t,x,u)) )
    
    #Mesh:
    MeshFile            = "$(@__DIR__)/../../temp/Leveque_SC$(SC).geo"
    NX                  = Int(ceil(2.0/(hp0*FesOrder)))
    TrMesh_Rectangle_Create!(MeshFile, -1.0, 1.0, NX, -1.0, 1.0, NX)
    
    #----------------------------------------------------------
    #INITIAL CONDITION:
    
    t_ini               = time()
    
    #Create solver structure:
    solver              = LIRKHyp_Start(ProblemData)
    solver.ProblemName  = "Leveque"
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
    solver.AMA_SizeOrder    = FesOrder
    solver.AMA_AnisoOrder   = 2
    solver.AMA_ProjN        = AMA_ProjN
    solver.AMA_ProjOrder    = AMA_ProjOrder
    solver.CA               = CA
    solver.CA_max           = 5*CA
    
    #Set initial and boundary conditions:
    solver.u0fun        = FW11((x) -> u0fun(x)) 
    solver.BC           = fill(BCW(BC_Dir), 4)
    
    #Initial condition:
    ConvFlag            = LIRKHyp_InitialCondition!(solver)
    tv                  = [0.0]
    uP                  = GetPolynomials(solver.u, solver.fes)
    uminv               = [minimum(uP[1])]
    umaxv               = [maximum(uP[1])]
    
    #Compute Lq error:
    errLq,              = LqError(solver, FW11((x) -> utheor(solver.t, x)), q=SpaceNorm)
    errLinf,            = LqError(solver, FW11((x) -> utheor(solver.t, x)), q=Inf)
    hmean               = 2.0*sqrt(solver.Omega/solver.mesh.nElems/TrElem_Area)
    errL2L2             = errLq
    etaL2L2             = solver.etaS
    errLqv              = [errLq]
    errLinfv            = [errLinf]
    
    #Print info:
    println("h=", sprintf1("%.2e",hmean), 
            ", eta_L2L2=", sprintf1("%.2e",etaL2L2),
            ", e_L2L2=", sprintf1("%.2e",errL2L2),
            ", e_Linf=", sprintf1("%.2e",errLinf),
            ", u_min=", sprintf1("%.4e",minimum(uP[1])),
            ", u_max=", sprintf1("%.4e",maximum(uP[1])),
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
        
#             figure(figv[1].number)
#             PyPlot.cla()
#             PlotContour(solver.u[1], solver.fes)
#             axis("equal")
#             if SaveFig
#                 savefig("$(VideosUbi)SC$(SC)_Contour_$(nb_SaveFig).png", dpi=800, pad_inches=0)
#             end
#             
#             figure(figv[2].number)
#             PyPlot.cla()
#             PlotContour(solver.u[1], solver.fes, cmap="bwr")
#             PlotMesh!(solver.mesh, color="k")
#             axis("equal")
#             if SaveFig
#                 savefig("$(VideosUbi)SC$(SC)_Mesh_$(nb_SaveFig).png", dpi=800, pad_inches=0)
#             end
            
            figure(figv[3].number)
            PyPlot.cla()
#             semilogy(solver.tv, solver.etaSv, ".-b")
#             semilogy(solver.tv, solver.etaTv, ".-g")
#             semilogy(solver.tv, solver.etaAv, ".-r")
            plot(solver.tv, solver.nDofv, ".-k")
            
            figure(figv[4].number)
            PyPlot.cla()
            plot(tv, umaxv.-1.0, ".-b")
            plot(tv, uminv, ".-g")
        
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
            save("$(ResUbi)LIRKHyp_SC$(SC)_$(nb_SaveRes).jld2", "StudyCase", "Leveque",
                "delta", delta,
                "ConvFlag", ConvFlag, "solver", save(solver), 
                "tv", tv, "uminv", uminv, "umaxv", umaxv, 
                "errL2L2", errL2L2, "etaL2L2", etaL2L2, 
                "errLqv", errLqv, "errLinfv", errLinfv )
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
    
#         ConvFlag            = LIRKHyp_Step!(solver)
        if PostAdapt
            ConvFlag            = LIRKHyp_Step_Post!(solver)
        else
            ConvFlag            = LIRKHyp_Step_Pre!(solver)
        end
        if ConvFlag<=0
            break
        end
        push!(tv, solver.t)
        uP                  = GetPolynomials(solver.u, solver.fes)
        push!(uminv, minimum(uP[1]))
        push!(umaxv, maximum(uP[1]))
        
        #Compute Lq error:
        errLq,              = LqError(solver, FW11((x) -> utheor(solver.t, x)), q=SpaceNorm)
        errLinf,            = LqError(solver, FW11((x) -> utheor(solver.t, x)), q=Inf)
        hmean               = 2.0*sqrt(solver.Omega/solver.mesh.nElems/TrElem_Area)
        errL2L2             = errL2L2_(solver.t-solver.Deltat, errL2L2, solver.t, errLq)
        etaL2L2             = errL2L2_(solver.t-solver.Deltat, etaL2L2, solver.t, solver.etaS+solver.etaT)
        push!(errLqv, errLq)
        push!(errLinfv, errLinf)
        
        #Print info:
        println("h=", sprintf1("%.2e",hmean), 
                ", eta_L2L2=", sprintf1("%.2e",etaL2L2),
                ", e_L2L2=", sprintf1("%.2e",errL2L2),
                ", i_ST=", sprintf1("%.4f", etaL2L2/errL2L2), 
                ", e_Linf=", sprintf1("%.2e",errLinf),
                ", u_min=", sprintf1("%.4e",minimum(uP[1])),
                ", u_max=", sprintf1("%.4e",maximum(uP[1])),
                ", time=", time()-t_ini)
        
        #Save figures and/or results?
        PlotSol()
        SaveSol()
        
    end

    return solver
    
end
