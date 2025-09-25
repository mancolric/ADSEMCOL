#TODO:
#   -Correct Brenner's model for multiple species.
#   -avxt macro.
#   -Output JET to text file.

# include("../../AuxiliaryFunctions-v2/basic.jl")
# import Base

# include("../../FEM/FEM.jl") #this includes basic
# include("../../AuxiliaryFunctions-v2/ODESolvers.jl")
# import Base
# using ILUZero
# using IterativeSolvers
# using Formatting

include("../src_LIRKHyp/LIRKHyp.jl")

#------------------------------------------------------------------
#STRUCTURES WITH CONSTITUTIVE MODELS AND BOUNDARY CONDITIONS:

Base.@kwdef mutable struct SWE <: ConstModels

    #Model's characteristic fields:
    epsilon         ::Float64           = 0.0   #Rusanov-type viscosity
    gamma           ::Float64           = 0.0   #Friction coefficient
    g               ::Float64           = 9.8   #Gravity force in -y direction
    CSS             ::Float64           = 0.1   #Subgrid stabilization
    CW              ::Float64           = 50.0  #Boundary penalty (50.0-200.0 for IIPG)
    b               ::FW1               = FW1( (x)-> zeros(size(x[1])) )
    gradb           ::FW11              = FW11((x)-> [ 0.0*x[1], 0.0*x[2]] )
    
    #Dependent variables:
#     DepVars         ::Vector{String}    = ["h","q1","q2","v1","v2","epsilon","D_penalty"]

    #Mandatory fields:
    nVars           ::Int               = 3     #Number of discretization variables: h, q1, q2

end

#Slip condition at a wall:
#   q_n             = 0     zero normal velocity
#   epsilon*deta/dn = 0     zero mass flux
#   epsilon*dq_t/dn = 0     zero tangential velocity flux
mutable struct SlipAdiabatic <: BoundConds

end

#Subsonic inlet. 3 conditions:
#   normal diff flux for eta    = 0
#   q1                          = q1_BC
#   q2                          = q2_BC
mutable struct SubsonicInlet1 <: BoundConds
    fun             ::FWt11     #must return [q_1, q_2]
end

#Subsonic outlet. 3 conditions:
#   eta                     = eta_BC
#   normal diff flux for q1 = 0
#   normal diff flux for q2 = 0
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

# ConstModels     = Union{GasIdeal}
# BoundConds      = Union{SlipAdiabatic, SubsonicInlet1, SubsonicOutlet1, DoNothing1}


#-------------------------------------------------------------------------------
#LOAD AUXILIARY FUNCTIONS:

#Return index corresponding to dependent variable "var":
function DepVarIndex(model::ConstModels, var::String)
    return findfirst(model.DepVars.==var)
end

function DepVars(model::SWE, t::Float64, x::Vector{<:AMF64},
    u::Vector{<:AMF64}, vout::Vector{String})

    eta         = u[1]
    q1          = u[2]
    q2          = u[3]
    b           = model.b(x)
    h           = eta-b
    nout        = length(vout)
    xout        = Vector{Vector{Array{Float64,ndims(q1)}}}(undef,nout)
    for ivar in eachindex(vout)
        vble    = vout[ivar]
        if vble=="h"
            xout[ivar]      = [h]
        elseif vble=="q1"
            xout[ivar]      = [q1]
        elseif vble=="q2"
            xout[ivar]      = [q2]
        elseif vble=="b"
            xout[ivar]      = [b]
        elseif vble=="v1"
            xout[ivar]      = [@tturbo @. q1/h]
        elseif vble=="v2"
            xout[ivar]      = [@tturbo @. q2/h]
        elseif vble=="epsilon"
            xout[ivar]      = [fill(model.epsilon, size(u[1]))]
        elseif vble=="eta"
            xout[ivar]      = [eta]
        elseif vble=="D_penalty"
            D_penalty       = model.epsilon
            xout[ivar]      = [fill(D_penalty, size(u[1]))]
       else
           error("Variable $(vble) not supported")
        end
    end

    return xout

end

#-------------------------------------------------------------------------------
#MANDATORY FUNCTIONS:


include("ShallowWater_fluxes.jl")
include("ShallowWater_BC.jl")

# #Compute normalization factors from solution. Mass matrix has already been computed.
function nFactsCompute!(solver::SolverData{SWE})


    #Normalization factors:
    solver.nFacts                   .= 1.0

    return

end

#Function to evaluate flux and source terms at quadrature nodes:
function FluxSource!(model::SWE, _qp::TrIntVars, ComputeJ::Bool)

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
    eta             = u[1]
    q1              = u[2]
    q2              = u[3]
    b               = model.b(x)
    h               = eta-b
    v1              = @tturbo @. q1/h
    v2              = @tturbo @. q2/h
    g               = model.g
    
    #Mesh size:
    A_Elems         = areas(_qp.Integ2D.mesh)
    h_Elems         = @tturbo @. sqrt(A_Elems)
    hp              = h_Elems./_qp.FesOrder * ones(1, _qp.nqp)
    hp_min          = _hmin(_qp.Integ2D.mesh)./_qp.FesOrder * ones(1, _qp.nqp)
    
    #Characteristic velocity:
    lambda          = @tturbo @. sqrt(v1*v1 + v2*v2) + sqrt(g*h)
    
    #Compute fluxes and source term:
    HyperbolicFlux!(model, x, u, ComputeJ, f, df_du)
    #
    epsilon_qp      = @tturbo @. model.epsilon + 0.0*u[1]
    epsilonFlux!(model, epsilon_qp, du, ComputeJ, f, _qp.df_dgradu)
    #
    Source!(model, x, u, du, ComputeJ, Q, dQ_du, dQ_du_dx)

    #Subgrid viscosity:
    epsilon_SS          = @tturbo @. model.CSS*lambda*hp
    epsilonFlux!(model, epsilon_SS, duB, ComputeJ, _qp.fB, _qp.dfB_dgraduB)

    #CFL number:
    Deltat_CFL_lambda   = @tturbo @. $minimum(hp_min/lambda)
    Deltat_CFL_epsilon  = @tturbo @. $minimum(hp_min*hp_min/model.epsilon)
    _qp.Deltat_CFL      = min(Deltat_CFL_lambda, Deltat_CFL_epsilon)

    return

end

#bflux! functions are in ShallowWater_BC
