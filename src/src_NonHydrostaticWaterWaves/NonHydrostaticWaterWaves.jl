include("../src_LIRKHyp/LIRKHyp.jl")

#------------------------------------------------------------------
#STRUCTURES WITH CONSTITUTIVE MODELS AND BOUNDARY CONDITIONS:

Base.@kwdef mutable struct NHWW <: ConstModels

    #Model's characteristic fields:
    epsilon         ::Float64           = 0.0   #Rusanov-type viscosity
    gamma           ::Float64           = 2.0   #Parameter
    g               ::Float64           = 9.8   #Gravity force in -z direction
    c               ::Float64           = 5*sqrt(10*10)
    h0              ::Float64           = 1.0
    CSS             ::Float64           = 0.1   #Subgrid stabilization
    CW              ::Float64           = 50.0  #Boundary penalty (50.0-200.0 for IIPG)
    b               ::FW1               = FW1( (x)-> zeros(size(x[1])) )
    
    #Mandatory fields:
    nVars           ::Int               = 6     #Number of discretization variables: h, q1, q2, b

end

#SlipAdiabatic boundary. 5 conditions (1 Dirichlet):
#   q_n                         = 0
#   normal diff flux for h      = 0
#   normal diff flux for q_t    = 0
#   normal diff flux for q_3    = 0
#   normal diff flux for p      = 0
mutable struct SlipAdiabatic <: BoundConds

end

#Subsonic inlet. 5 conditions (3 Dirichlet):
#   normal diff flux for eta    = 0
#   q1                          = q1_BC
#   q2                          = q2_BC
#   q3                          = q3_BC
#   normal diff flux for p      = p_BC
mutable struct SubsonicInlet1 <: BoundConds
    fun             ::FWt11     #must return [q_1, q_2]
end

#Subsonic outlet. 5 conditions (2 Dirichlet):
#   eta                     = eta_BC
#   normal diff flux for q1 = 0
#   normal diff flux for q2 = 0
#   normal diff flux for q3 = 0
#   p                       = p_BC
mutable struct SubsonicOutlet1 <: BoundConds
    fun             ::FWt11     #must return [h]
end

#Supersonic outlet/ do nothing: 3 conditions:
#   normal diff dlux for eta    = 0
#   normal diff dlux for q1     = 0
#   normal diff dlux for q2     = 0
mutable struct SupersonicOutlet1 <: BoundConds

end
mutable struct DoNothing1 <: BoundConds

end

#-------------------------------------------------------------------------------
#LOAD AUXILIARY FUNCTIONS:

#Return index corresponding to dependent variable "var":
function DepVarIndex(model::ConstModels, var::String)
    return findfirst(model.DepVars.==var)
end

function DepVars(model::NHWW, t::Float64, x::Vector{<:AMF64},
    u::Vector{<:AMF64}, vout::Vector{String})

    eta         = u[1]
    q1          = u[2]
    q2          = u[3]
    q3          = u[4]
    P           = u[5]
    b           = u[6]
    h           = eta-b
    c           = model.c
    nout        = length(vout)
    xout        = Vector{Vector{Array{Float64,ndims(q1)}}}(undef,nout)
    for ivar in eachindex(vout)
        vble    = vout[ivar]
        if vble=="eta"
            xout[ivar]      = [eta]
        elseif vble=="q1"
            xout[ivar]      = [q1]
        elseif vble=="q2"
            xout[ivar]      = [q2]
        elseif vble=="q3"
            xout[ivar]      = [q3]
        elseif vble=="P"
            xout[ivar]      = [P]
        elseif vble=="p"
            xout[ivar]      = [@tturbo @. P/h + c*c*log(h/model.h0)]
        elseif vble=="b"
            xout[ivar]      = [b]
        elseif vble=="v1"
            xout[ivar]      = [@tturbo @. q1/h]
        elseif vble=="v2"
            xout[ivar]      = [@tturbo @. q2/h]
        elseif vble=="v3"
            xout[ivar]      = [@tturbo @. q3/h]
        elseif vble=="vr"
            xout[ivar]      = [@tturbo @. sqrt(q1*q1+q2*q2)/h]
        elseif vble=="epsilon"
            xout[ivar]      = [fill(model.epsilon, size(u[1]))]
        elseif vble=="h"
            xout[ivar]      = [h]
        elseif vble=="D_penalty"
            D_penalty       = max(model.epsilon)
            xout[ivar]      = [fill(D_penalty, size(u[1]))]
       else
           error("Variable $(vble) not supported")
        end
    end

    return xout

end

#-------------------------------------------------------------------------------
#MANDATORY FUNCTIONS:


include("NonHydrostaticWaterWaves_fluxes.jl")
include("NonHydrostaticWaterWaves_BC.jl")

# #Compute normalization factors from solution. Mass matrix has already been computed.
function nFactsCompute!(solver::SolverData{NHWW})
    
    #Normalization factors:
    h0                      = solver.model.h0
    solver.nFacts[1]        = h0
    solver.nFacts[2:4]      .= h0*sqrt(solver.model.g*h0)
    solver.nFacts[5]        = h0*solver.model.g*h0          #~ h u^2
#     solver.nFacts[5]        = h0*(solver.model.c)^2
    solver.nFacts[6]        = h0
#     solver.nFacts           .= 1.0
    
    return

end

#Function to evaluate flux and source terms at quadrature nodes:
function FluxSource!(model::NHWW, _qp::TrIntVars, ComputeJ::Bool)

    t               = _qp.t
    x               = _qp.x
    u               = _qp.u
    du              = _qp.gradu
    duB             = _qp.graduB
    f               = _qp.f
    Q               = _qp.Q
    df_du           = _qp.df_du
    dQ_du           = _qp.dQ_du
    dQ_du_dx        = _qp.dQ_dgradu

    #Get variables:
    g               = model.g
    c               = model.c
    h0              = model.h0
    eta             = u[1]
    q1              = u[2]
    q2              = u[3]
    q3              = u[4]
    P               = u[5]
    b               = u[6]
    h               = eta-b
    p               = @tturbo @. P/h + (c*c)*log(h/h0)
    v1              = @tturbo @. q1/h
    v2              = @tturbo @. q2/h
    v3              = @tturbo @. q3/h
    
    #Mesh size:
    A_Elems         = areas(_qp.Integ2D.mesh)
    h_Elems         = @tturbo @. sqrt(A_Elems)
    hp              = h_Elems./_qp.FesOrder * ones(1, _qp.nqp)
    hp_min          = _hmin(_qp.Integ2D.mesh)./_qp.FesOrder * ones(1, _qp.nqp)
    
    #Characteristic velocities (original system and relaxed one)
    #   lambdac is used for subgrid stabilization and tau_char 
    #   lambda0 for definition of CFL
    lambda0             = @tturbo @. sqrt(v1*v1 + v2*v2 + v3*v3) + sqrt(g*h + p)
    lambdac             = @tturbo @. sqrt(v1*v1 + v2*v2 + v3*v3) + sqrt(g*h + p + c*c)
    
    #CFL number:
    Deltat_CFL_lambda0  = @tturbo @. $minimum(hp_min/(lambda0+1e-12))
    Deltat_CFL_lambdac  = @tturbo @. $minimum(hp_min/(lambdac+1e-12))
    Deltat_CFL_epsilon  = @tturbo @. $minimum(hp_min*hp_min/(model.epsilon+1e-12))
    _qp.Deltat_CFL      = min(Deltat_CFL_lambda0, Deltat_CFL_epsilon)
    
    #Compute fluxes and source term:
    HyperbolicFlux!(model, u, ComputeJ, f, df_du)
    #
    epsilon_qp          = @tturbo @. model.epsilon + 0.0*u[1]
    epsilonFlux!(model, epsilon_qp, du, ComputeJ, f, _qp.df_dgradu)
    #   
    tau_char            = min(Deltat_CFL_lambdac, Deltat_CFL_epsilon)
    Source!(model, x, tau_char, u, du, ComputeJ, Q, dQ_du, dQ_du_dx)

    #Subgrid viscosity:
    epsilon_SS          = @tturbo @. model.CSS*lambdac*hp
    epsilonFlux!(model, epsilon_SS, duB, ComputeJ, _qp.fB, _qp.dfB_dgraduB)
    
    return

end

#bflux! functions are in NonHydrostaticWaterWaves_BC
