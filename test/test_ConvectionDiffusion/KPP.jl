cd(@__DIR__)
include("test_ConvectionDiffusion.jl")

function KPP(hp0::Float64, FesOrder::Int;
    epsilon::Float64=1e-4, 
    Deltat0::Float64=1e-6, 
    TolT::Float64=1e-3, TolS::Float64=1e-5, AMA_MaxIter::Int=200, tf::Float64=1.0,
    PostAdapt::Bool=false, 
    #
    AMA_rDeltat::Float64=0.3,
    SaveFig::Bool=false, wFig::Float64=9.50, hFig::Float64=6.50, 
    Nt_SaveFig::Int=5, 
    PlotRes::Bool=false, 
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
        
        a1  = @avxt @. sin(u[1])/u[1]
        a2  = @avxt @. cos(u[1])/u[1]
        return [a1, a2]
    end
    function dafun(t::Float64, x::Vector{<:AAF64}, u::Vector{<:AAF64})
        
        da1     = @avxt @. cos(u[1])/u[1] - sin(u[1])/u[1]^2
        da2     = @avxt @. -sin(u[1])/u[1] - cos(u[1])/u[1]^2
        return [da1, da2]
    end
    
    #Viscosity:
    function epsilonfun(t::Float64, x::Vector{<:AAF64}, u::Vector{<:AAF64})
    
        epsilon0    = @mlv 0.0*x[1]+epsilon
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
        
        r           = @avxt @. sqrt(x[1]*x[1]+x[2]*x[2])
        u0          = @avxt @. SmoothHeaviside(r-1.0, epsilon, 3.5*pi, 0.25*pi)
        return [ u0 ]
        
    end
    
    #Boundary conditions:
    function uLeft(t::Float64, x::Vector{<:AAF64})
        return [0.0*x[1].+1.0]
    end
    function uRight(t::Float64, x::Vector{<:AAF64})
        return [0.0*x[1].+1.0]
    end
    function qTopBottom(t::Float64, x::Vector{<:AAF64})
        return [ zeros(size(x[1])) ]
    end
#     BC_left         = Dirichlet(FWt11((t,x)->uLeft(t,x)))
#     BC_right        = Dirichlet(FWt11((t,x)->uRight(t,x)))
    BC_horiz        = Neumann(FWt11((t,x)->qTopBottom(t,x)))
    
    #Structure with convection--diffusion problem data:
    ProblemData     = SCD(  FWt21((t,x,u)->afun(t,x,u)), 
                            FWt21((t,x,u)->epsilonfun(t,x,u)),
                            FWt21((t,x,u)->Qfun(t,x,u)), 
                            FWt21((t,x,u)->dafun(t,x,u)), 
                            FWt21((t,x,u)->depsilonfun(t,x,u)),
                            FWt21((t,x,u)->dQfun(t,x,u)) )
    
    #Mesh:
    MeshFile            = "../temp/KPP_SC$(SC).geo"
    NX                  = Int(ceil(4.0/(hp0*FesOrder)))
    TrMesh_Rectangle_Create!(MeshFile, -2.0, 2.0, NX, -2.5, 1.5, NX)
    
    #----------------------------------------------------------
    #INITIAL CONDITION:
    
    t_ini               = time()
    
    #Create solver structure:
    solver              = LIRKHyp_Start(ProblemData)
    solver.ProblemName  = "PolarBurgers"
    solver.SC           = SC
    solver.MeshFile     = MeshFile
    solver.nBounds      = 4
    solver.FesOrder     = FesOrder
    solver.Deltat0      = Deltat0
    solver.tf           = tf
    solver.SpaceAdapt   = true
    solver.TolS_max     = TolS
    solver.TolS_min     = 0.0*TolS
#     solver.SfS          = 1.0
#     @warn "SfS"
    solver.TimeAdapt    = true
    solver.TolT         = TolT
    solver.AMA_MaxIter  = AMA_MaxIter
    solver.AMA_SizeOrder    = FesOrder
    solver.AMA_AnisoOrder   = 2
    solver.AMA_rDeltat      = AMA_rDeltat
    
    #Set initial and boundary conditions:
    solver.u0fun        = FW11((x) -> u0fun(x)) 
#     solver.BC           = [BCW(BC_horiz), BCW(BC_right), BCW(BC_horiz), BCW(BC_left)]
    solver.BC           = [BCW(BC_horiz), BCW(BC_horiz), BCW(BC_horiz), BCW(BC_horiz)]
    
    #Initial condition:
    ConvFlag            = LIRKHyp_InitialCondition!(solver)
    tv                  = [0.0]
    uP                  = GetPolynomials(solver.u, solver.fes)
    uminv               = [minimum(uP[1])]
    umaxv               = [maximum(uP[1])]
    
    #Function to plot solution:
    figv                = Vector{Figure}(undef,3)
    if PlotRes
        for ii=1:length(figv)
            figv[ii]    = figure()
        end
    end
    t_lastFig           = 0.0
    NPlotSave           = 0
    NSave               = 0
    function PlotSol()
    
        if !PlotRes
            return
        end
        
        figure(figv[1].number)
        PyPlot.cla()
#         plotNodes((x1,x2)->x1, solver.u[1], solver.fes)
        PlotContour(solver.u[1], solver.fes)
        axis("equal")
        
        figure(figv[2].number)
        PyPlot.cla()
        PlotContour(solver.u[1], solver.fes, cmap="bwr")
        PlotMesh!(solver.mesh, color="k")
        axis("equal")
        
        figure(figv[3].number)
        PyPlot.cla()
#         semilogy(solver.tv, solver.etaSv, ".-b")
#         semilogy(solver.tv, solver.etaTv, ".-g")
#         semilogy(solver.tv, solver.etaAv, ".-r")
        plot(solver.tv, solver.nElemsv, ".-k")
        
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
                        ct_SaveRes==Nt_SaveRes || solver.t==tf || solver.t==0.0 )
            save("$(ResUbi)LIRKHyp_SC$(SC)_$(nb_SaveRes).jld2", "StudyCase", "KPP",
                "tv", tv, "uminv", uminv, "umaxv", umaxv, 
                "epsilon", epsilon, "solver", save(solver))
            t_lastRes   += Deltat_SaveRes
            ct_SaveRes  = 0
            nb_SaveRes  += 1
        end
        return
        
    end
    SaveSol()
    
    #Display info:
    hmean       = 2.0*sqrt(solver.Omega/solver.mesh.nElems/TrElem_Area)
    etaL2,      = etaS_reconstruction(solver, solver.u, q=2)
    etaLinf,    = etaS_reconstruction(solver, solver.u, q=Inf)
    println("h=", sprintf1("%.2e",hmean), 
        ", etaL2=", sprintf1("%.2e",etaL2),
        ", etaLinf=", sprintf1("%.2e",etaLinf),
        ", Delta u_min=", sprintf1("%.4e",0.25*pi-minimum(uP[1])),
        ", Delta u_max=", sprintf1("%.4e",maximum(uP[1])-3.5*pi),
        ", time=", time()-t_ini)
        
    #----------------------------------------------------------
    #MARCH IN TIME:
    
    t_ini               = time()
    
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
        
        NPlotSave   += 1
        if solver.t-t_lastFig>=2e-3 || solver.t==tf || NPlotSave==Nt_SaveFig
            NPlotSave   = 0
            t_lastFig   = solver.t
            PlotSol()
        end
        SaveSol()
        
        #Print info:
        println("h=", sprintf1("%.2e",hmean), 
                ", Delta u_min=", sprintf1("%.4e",0.25*pi-minimum(uP[1])),
                ", Delta u_max=", sprintf1("%.4e",maximum(uP[1])-3.5*pi),
                ", time=", time()-t_ini)
                
    end
    
    return
    
end
