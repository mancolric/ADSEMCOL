include("../src_LIRKHyp/LIRKHyp.jl")

#------------------------------------------------------------------
#STRUCTURES WITH CONSTITUTIVE MODELS AND BOUNDARY CONDITIONS:

abstract type ScalarConvectionDiffusion <: ConstModels end

mutable struct SCD <: ScalarConvectionDiffusion

    #Model's characteristic fields. These functions receive (t,[x1,x2],[u]) and return
    a               ::FWt21                     #Returns velocity [a1, a2]
    epsilon         ::FWt21                     #Returns viscosity coefficient [epsilon]
    Q               ::FWt21                     #Returns source [Q]
    
    #Functions to compute the jacobians:
    da_du           ::FWt21
    depsilon_du     ::FWt21
    dQ_du           ::FWt21
    
    #Stabilization variables:
    CSS             ::Float64   #Subgrid stabilization
    CW              ::Float64   #Boundary penalty (50.0-200.0 for IIPG)
    CDC             ::Float64   #Discontinuity capturing
    
    #Mandatory fields:
    nVars           ::Int                 
    
    SCD()           = new()
    
end

function SCD(a::FWt21, epsilon::FWt21, Q::FWt21, 
    da_du::FWt21, depsilon_du::FWt21, dQ_du::FWt21)

    PD              = SCD()
    PD.a            = a
    PD.epsilon      = epsilon
    PD.Q            = Q
    PD.da_du        = da_du
    PD.depsilon_du  = depsilon_du
    PD.dQ_du        = dQ_du
    PD.CSS          = 0.1
    PD.CW           = 50.0
    PD.CDC          = 0.0
    PD.nVars        = 1
#     PD.DepVars  = String[]
    return PD
    
end

mutable struct Dirichlet <: BoundConds
    uDir            ::FWt11     #must return Dirichlet condition at the boundary [u]
end
mutable struct Neumann <: BoundConds
    q               ::FWt11     #must return diffusive flux [q=-epsilon*du/dn]
end

#-------------------------------------------------------------------------------
#LOAD LIRHyp and auxiliary functions:

include("ConvectionDiffusion_fluxes.jl")
include("ConvectionDiffusion_BC.jl")

#If x and u are AbstractVectors, they will automatically be converted to 
#AbstractMatrices. However, the output will still be a matrix:
function DepVars(model::SCD, t::Float64, x::Vector{<:AMF64}, u::Vector{<:AMF64}, 
    vout::Vector{String})
    
    nout        = length(vout)
    xout        = Vector{Vector{Matrix{Float64}}}(undef,nout)
    for ivar in eachindex(vout)
        vble    = vout[ivar]
        if vble=="u"
            xout[ivar]  = [ copy(u[1]) ]
        elseif vble=="lambda_max"
            a           = model.a(t, x, u)
            da_du       = model.da_du(t, x, u)
            #ahat_i = d(a_i u)/du = da_i/du * u + a_i 
            lambda      = @avxt @. sqrt( (da_du[1]*u[1]+a[1])^2 + (da_du[2]*u[1]+a[2])^2 ) 
            xout[ivar]  = [ lambda ]
        else
            error("Variable $(vble) not supported")
        end
    end
    
    return xout
    
end

#-------------------------------------------------------------------------------
#MANDATORY FUNCTIONS:

#Compute normalization factors from solution. Mass matrix has already been computed.
function nFactsCompute!(solver::SolverData{<:SCD})

    solver.nFacts[1]    = 1.0
    
    #Normalize with L2-norm of variable:
#     solver.nFacts[1]    = sqrt(dot(solver.u[1], solver.MII, solver.u[1])/solver.Omega)
        
    return
    
end

#Function to evaluate flux and source terms at quadrature nodes:
function FluxSource!(model::SCD, _qp::TrIntVars, ComputeJ::Bool)

    t               = _qp.t
    x               = _qp.x
    u               = _qp.u
    du              = _qp.gradu
    duB             = _qp.graduB
    metric          = _qp.Integ2D.mesh.metric
    
    #Natural viscosity:
    a                       = model.a(t,x,u)
    da_du                   = Vector{Matrix{Float64}}(undef,2)
    lambda_max              = DepVars(model,t,x,u,["lambda_max"])[1][1]
    epsilon                 = model.epsilon(t,x,u)[1]
    depsilon_du             = zeros(0,0)
    if ComputeJ
        da_du               = model.da_du(t,x,u)
        depsilon_du         = model.depsilon_du(t,x,u)[1]
    end
    
    #Add artificial viscosity:
    if true
        delta               = 1e-3
        C1                  = 0.0
        C2                  = 1.0
        hp_min              = _hmin(_qp.Integ2D.mesh) ./_qp.FesOrder * ones(1, _qp.nqp)
        S                   = @tturbo @. SmoothHeaviside(log10(hp_min/delta), C2, 1.0, 0.0)
#         S                   = @tturbo @. SmoothHeaviside(hp_min-delta, C2*delta, 1.0, 0.0)
#         @tturbo @. epsilon  += C1*lambda_max*delta*S
        @tturbo @. epsilon  += C1*1.0*delta*1.0
    end
    if false
        #Impossible model!!
        #   epsilon= C1*a*delta * min((C2|grad u|/(u0/delta))^m, 1.0) 
        delta               = 1e-3
        C1                  = 2.0
        C2                  = 5.0
        m                   = 2.0
        u0                  = 1.0
        #|grad u|:
        du_norm             = @tturbo @. sqrt(du[1]*du[1]+du[2]*du[2])
        #Shock sensor:
        S                   = @tturbo @. min( (C2*du_norm/(u0/delta))^m, 1.0)
        #Viscosity:
        @tturbo @. epsilon  += C1*lambda_max*delta*S
        
    end
    
    #Convective and diffusive fluxes:
    ConvectiveFlux!(model, a, da_du, u, du, ComputeJ, 
        _qp.f, _qp.df_du, _qp.df_dgradu)
    DiffusiveFlux!(model, epsilon, depsilon_du, u, du, ComputeJ, 
        _qp.f, _qp.df_du, _qp.df_dgradu)
    
    #Evaluate subgrid stabilization flux:
    A_Elems             = areas(_qp.Integ2D.mesh)
    h_Elems             = @tturbo @. sqrt(A_Elems)
    hp                  = h_Elems./_qp.FesOrder * ones(1, _qp.nqp)
    epsilonSS           = @avxt @. model.CSS*lambda_max*hp
    if ComputeJ
        @avxt @. depsilon_du        = 0.0
#         @avxt @. depsilon_du        = model.CSS*hp * (a[1]*da_du[1]+a[2]*da_du[2])/anorm
    end
    DiffusiveFlux!(model, epsilonSS, depsilon_du, 
        u, duB, ComputeJ, 
        _qp.fB, _qp.dfB_du, _qp.dfB_dgraduB)
    
    #Evaluate source terms:
    _qp.Q[1]            .= model.Q(t,x,u)[1]
    if ComputeJ
        _qp.dQ_du[1]    .= model.dQ_du(t,x,u)[1]
    end
    
    #Deltat imposed by CFL=1 (do not use @avxt, it does not work well with $ symbol)
    hp_min              = _hmin(_qp.Integ2D.mesh)./_qp.FesOrder * ones(1, _qp.nqp)
    Deltat_CFL_a        = @. $minimum(hp_min/lambda_max)
    Deltat_CFL_epsilon  = @. $minimum(hp_min^2/epsilon) 
    _qp.Deltat_CFL      = min(Deltat_CFL_a, Deltat_CFL_epsilon)
    
    return
    
end

#bflux! functions are in ConvectionDiffusion_BC
