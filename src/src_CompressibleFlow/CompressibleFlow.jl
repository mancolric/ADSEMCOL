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

abstract type GasModel <: ConstModels end

Base.@kwdef mutable struct GasIdeal <: GasModel

    #Model's characteristic fields:
    epsilon         ::Float64           = 0.0
    nu              ::Float64           = 0.0
    beta            ::Float64           = 0.0
    kappa_rho_cv    ::Float64           = 0.0
    gamma           ::Float64           = 1.4
    g               ::Float64           = 0.0   #Gravity force in -y direction
    nReactions      ::Int               = 0
    nSpecies        ::Int               = 1
    CSS             ::Float64           = 0.1   #Subgrid stabilization
    CW              ::Float64           = 50.0  #Boundary penalty (50.0-200.0 for IIPG)
    
    #Dependent variables:
    DepVars         ::Vector{String}    = [ "rho", "p", "vx", "vy", "gamma", 
                                        "e_i", "lambda_max",
                                        "epsilon", "nu", "beta", "kappa_rho_cv",
                                        "D_penalty" ]

    #Mandatory fields:
    nVars           ::Int               = 4

end

abstract type ReactiveGas <: GasModel end

#Note:
#   Le      = kappa_rho_cp / D = kappa_rho_cv / gamma / D ~ 1
#   Da      = B*L0/U0 / exp(Ta/Tf)  ~ 1000
#   Ta/Tf   ~ 10
#   Tf/T0   ~ 6
#Normally, we do not solve the original problem (0), but a similar one
#in which rho~1, U~1, L~1. In that case,
#   p = (p0/rho0/U0^2) * 1 * 1 = a0^2 / gamma / U0^2 = 1/gamma/M0^2 <~ 1
#Indeed, we can check that p0/rho/U0^2 ~ 1 for p0=1bar, rho0=1Kg/m3, U0=340m/s.
#
#On the other hand, the thermodynamic variables do not depend on T, but on R*T.
#Note that R0*T0/U0^2 ~ p0/rho0/U0^2 <~ 1. Indeed, we can check that for air 
#(R0=287J/Kg/K, T0=298K).
#If we make T~1 in the similar problem, we must have R0*T0/U0^2 = R*T/U^2 ~ 1
#so R~1.
Base.@kwdef mutable struct GasFP <: ReactiveGas
    
    #Compressible flow characteristic fields:
    epsilon         ::Float64           = 0.0
    nu              ::Float64           = 0.0
    beta            ::Float64           = 0.0
    kappa_rho_cv    ::Float64           = 0.0
    gamma           ::Float64           = 1.4
    g               ::Float64           = 0.0   #Gravity force in -y direction
    nReactions      ::Int               = 1
    nSpecies        ::Int               = 2
    CSS             ::Float64           = 0.1   #Subgrid stabilization
    CW              ::Float64           = 50.0  #Boundary penalty (50.0-200.0 for IIPG)
    
    #Reaction's characteristic fields:
    hfF             ::Float64   = 500.0 #Q0=hFF-hFP=gamma/(gamma-1)*(RTf-RT0)/YF0
    hfP             ::Float64   = 0.0
    D               ::Float64   = kappa_rho_cv / gamma / 1.0
    B               ::Float64   = 1000*exp(10)
    RTa             ::Float64   = 10*6
    
    #Dependent variables:
    DepVars         ::Vector{String}    = [ "rho", "p", "vx", "vy", "gamma", 
                                        "e_i", "h_i", "D_i", 
                                        "mdot_i", "dmdot_ij", "lambda_max",
                                        "epsilon", "nu", "beta", "kappa_rho_cv",
                                        "D_penalty" ]
    
    #Mandatory fields:
    nVars           ::Int               = 5

    
end

#Slip + adiabatic condition at a wall. nSpecies+3 conditions:
#   m_n             = 0     zero normal velocity
#   fmass_(k,j) n_j = 0     no mass flux, k=1,...,nSpecies
#   tau_tj n_j      = 0     no viscous force in tangent direction
#   q_j n_j         = 0     no heat flux
mutable struct SlipAdiabatic <: BoundConds

end

#Subsonic inlet. nSpecies+3 conditions:
#   rhoY_k          = rhoY_k_BC
#   rhovx           = rhovx_BC 
#   rhovy           = rhovy_BC
#   (rhoe is extrapolated)
#   q_j n_j         = 0
mutable struct SubsonicInlet1 <: BoundConds
    fun             ::FWt11     #must return [rhoY_1,...,rhoY_nSpecies, rhovx, rhovy]
end

#Supersonic inlet. nSpecies+3 conditions:
#   rhoY_k          = rhoY_k_BC
#   rhovx           = rhovx_BC 
#   rhovy           = rhovy_BC
#   rhoE            = rhoE_BC
mutable struct SupersonicInlet1 <: BoundConds
    fun             ::FWt11     #must return [rhoY_1,...,rhoY_nSpecies, rhovx, rhovy, rhoE]
end

#Subsonic outlet: nSpecies+3 conditions:
#   p               = p_BC
#   (rhoY_k, rhov_i are extrapolated)
#   fmass_(k,j) n_j = 0
#   tau_ij n_j      = 0
mutable struct SubsonicOutlet1 <: BoundConds
    fun             ::FWt11     #must return [p]
end

#Supersonic outlet: nSpecies+3 conditions:
#   fmass_(k,j) n_j = 0
#   tau_nj n_j      = 0
#   v_t             = 0
#   q_j n_j         = 0
mutable struct SupersonicOutlet1 <: BoundConds

end

#Do nothing: nSpecies+3 conditions:
#   fmass_(k,j) n_j = 0
#   tau_ij n_j      = 0
#   q_j n_j         = 0
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

function DepVars(model::GasIdeal, t::Float64, x::Vector{<:AMF64},
    u::Vector{<:AMF64}, vout::Vector{String})
    
    nSpecies    = model.nSpecies
    rho         = sum(u[1:nSpecies])
    rhovx       = u[nSpecies+1]
    rhovy       = u[nSpecies+2]
    rhoE        = u[nSpecies+3]
    rhoe        = @mlv    rhoE-0.5*(rhovx*rhovx+rhovy*rhovy)/rho
    gamma       = model.gamma
    p           = @. abs((gamma-1.0)*rhoe)
#     p           = @mlv    (gamma-1.0)*rhoe
    nout        = length(vout)
    xout        = Vector{Vector{Array{Float64,ndims(rhoE)}}}(undef,nout)
    for ivar in eachindex(vout)
        vble    = vout[ivar]
        if vble=="rho"
            xout[ivar]      = [rho]
        elseif vble=="rhovx"
            xout[ivar]      = [rhovx]
        elseif vble=="rhovy"
            xout[ivar]      = [rhovy]
        elseif vble=="rhoE"
            xout[ivar]      = [rhoE]
        elseif vble=="vx"
            xout[ivar]      = [@mlv rhovx./rho]
        elseif vble=="vy"
            xout[ivar]      = [@mlv rhovy./rho]
        elseif vble=="v"
            xout[ivar]      = [@mlv sqrt(rhovx*rhovx+rhovy*rhovy)/rho]
        elseif vble=="rhoe"
            xout[ivar]      = [rhoe]
        elseif vble=="e"
            xout[ivar]      = [@avxt @. rhoe/rho]
        elseif vble=="p"
            xout[ivar]      = [p]
        elseif vble=="a"
            xout[ivar]      = [@mlv sqrt(abs(gamma*p/rho))]
        elseif vble=="lambda_max"
            xout[ivar]      = [@mlv sqrt(abs(gamma*p/rho)) + sqrt(rhovx*rhovx+rhovy*rhovy)/abs(rho)]
        elseif vble=="M"
            xout[ivar]      = [@mlv sqrt( (rhovx*rhovx+rhovy*rhovy)/(rho*rho) 
                                        / abs(gamma*p/rho)) ]
        elseif vble=="RT"
            xout[ivar]      = [@mlv abs(p/rho)]
        elseif vble=="epsilon"
            xout[ivar]      = [fill(model.epsilon, size(u[1]))]
        elseif vble=="nu"
            xout[ivar]      = [fill(model.nu, size(u[1]))]
#             nu              = model.nu
#             xout[ivar]      = [ @avxt @. SmoothHeaviside(x[1]-0.55, 5e-2, nu, 0.0) ]
        elseif vble=="beta"
            xout[ivar]      = [fill(model.beta, size(u[1]))]
        elseif vble=="kappa_rho_cv"
            xout[ivar]      = [fill(model.kappa_rho_cv, size(u[1]))]
        elseif vble=="D_penalty"
            D_penalty       = max(model.epsilon, model.nu, model.beta, model.kappa_rho_cv)
            xout[ivar]      = [fill(D_penalty, size(u[1]))]
        elseif vble=="rhos"
            xout[ivar]      = [@mlv rho*(log(abs(p))-gamma*log(abs(rho)))]
        elseif vble=="gamma"
            xout[ivar]      = [fill(gamma,size(u[1]))]
        elseif vble=="e_i"
            #(Assume all species have same energy)
            xout[ivar]      = Vector{Matrix{Float64}}(undef,nSpecies)
            for ss=1:nSpecies
                xout[ivar][ss]  = @mlv rhoe/rho
            end
        elseif vble=="dp_du"
            #From 
            #   rhoE    = p/(gamma-1) + 0.5*((rhovx)^2 + (rhovy)^2)/rho
            #we obtain
            #   0       = 1/(gamma-1) * dp/drhoYk - 0.5*v^2
            #   0       = 1/(gamma-1) * dp/drhovx + vx
            #   0       = 1/(gamma-1) * dp/drhovy + vy
            #   1       = 1/(gamma-1) * dp/drhoE 
            dp_dout         = Vector{Matrix{Float64}}(undef, nSpecies+3)
            dp_dout[1]      = @mlv (gamma-1.0)*(0.5*(rhovx*rhovx+rhovy*rho)/(rho*rho))
            dp_dout[2]      = @mlv -(gamma-1.0)*rhovx/rho
            dp_dout[3]      = @mlv -(gamma-1.0)*rhovy/rho
            dp_dout[4]      = fill(gamma-1.0, size(u[1]))
            xout[ivar]      = dp_dout
        else
            error("Variable $(vble) not supported")
        end
    end
    
    return xout
    
end

function DepVars(model::GasFP, t::Float64, x::Vector{<:AMF64},
    u::Vector{<:AMF64}, vout::Vector{String})
    
    #Recall that F and P have the same cv, so 
    #   e   = eF*YF + eP*YP = hfF*YF + hfP*YP + 1/(gamma-1)*RT
    
    nSpecies    = model.nSpecies
    rho         = sum(u[1:nSpecies])
    rhovx       = u[nSpecies+1]
    rhovy       = u[nSpecies+2]
    rhoE        = u[nSpecies+3]
    rhoe        = @mlv rhoE-0.5*(rhovx*rhovx+rhovy*rhovy)/rho
    gamma       = model.gamma
    hfF         = model.hfF
    hfP         = model.hfP
    p           = @mlv (rhoe-u[1]*hfF-u[2]*hfP)*(gamma-1.0)
    RT          = @mlv p/rho
    nout        = length(vout)
    xout        = Vector{Vector{Array{Float64,ndims(rhoE)}}}(undef,nout)
    for ivar in eachindex(vout)
        vble    = vout[ivar]
        if vble=="rho"
            xout[ivar]      = [rho]
        elseif vble=="rhoY_F"
            xout[ivar]      = [u[1]]
        elseif vble=="rhoY_P"
            xout[ivar]      = [u[2]]
        elseif vble=="Y_F"
            xout[ivar]      = [@mlv u[1]/rho]
        elseif vble=="Y_P"
            xout[ivar]      = [@mlv u[2]/rho]
        elseif vble=="rhovx"
            xout[ivar]      = [rhovx]
        elseif vble=="rhovy"
            xout[ivar]      = [rhovy]
        elseif vble=="rhoE"
            xout[ivar]      = [rhoE]
        elseif vble=="vx"
            xout[ivar]      = [@mlv rhovx./rho]
        elseif vble=="vy"
            xout[ivar]      = [@mlv rhovy./rho]
        elseif vble=="v"
            xout[ivar]      = [@mlv sqrt(rhovx*rhovx+rhovy*rhovy)/rho]
        elseif vble=="rhoe"
            xout[ivar]      = [rhoe]
        elseif vble=="e"
            xout[ivar]      = [@avxt @. rhoe/rho]
        elseif vble=="p"
            xout[ivar]      = [p]
        elseif vble=="a"
            xout[ivar]      = [@mlv sqrt(abs(gamma*p/rho))]
        elseif vble=="lambda_max"
            xout[ivar]      = [@mlv sqrt(abs(gamma*p/rho)) + sqrt(rhovx*rhovx+rhovy*rhovy)/abs(rho)]
        elseif vble=="M"
            xout[ivar]      = [@mlv sqrt( (rhovx*rhovx+rhovy*rhovy)/(rho*rho) 
                                        / abs(gamma*p/rho)) ]
        elseif vble=="RT"
            xout[ivar]      = [@mlv abs(p/rho)]
        elseif vble=="epsilon"
            xout[ivar]      = [fill(model.epsilon, size(u[1]))]
        elseif vble=="nu"
            xout[ivar]      = [fill(model.nu, size(u[1]))]
        elseif vble=="beta"
            xout[ivar]      = [fill(model.beta, size(u[1]))]
        elseif vble=="kappa_rho_cv"
            xout[ivar]      = [fill(model.kappa_rho_cv, size(u[1]))]
        elseif vble=="D_penalty"
            D_penalty       = max(model.epsilon, model.nu, model.beta, model.kappa_rho_cv)
            xout[ivar]      = [fill(D_penalty, size(u[1]))]
        elseif vble=="rhos"
            xout[ivar]      = [@mlv rho*(log(abs(p))-gamma*log(abs(rho)))]
        elseif vble=="gamma"
            xout[ivar]      = [fill(gamma,size(u[1]))]
        elseif vble=="e_i"
            xout[ivar]      = Vector{Matrix{Float64}}(undef,nSpecies)
            xout[ivar][1]   = @mlv hfF + RT/(gamma-1.0)
            xout[ivar][2]   = @mlv hfP + RT/(gamma-1.0)
        elseif vble=="h_i"
            xout[ivar]      = Vector{Matrix{Float64}}(undef,nSpecies)
            xout[ivar][1]   = @mlv hfF + RT*(gamma/(gamma-1.0))
            xout[ivar][2]   = @mlv hfP + RT*(gamma/(gamma-1.0))
        elseif vble=="D_i"
            xout[ivar]      = Vector{Matrix{Float64}}(undef,nSpecies)
            xout[ivar][1]   = fill(model.D, size(u[1]))
            xout[ivar][2]   = fill(model.D, size(u[1]))
        elseif vble=="mdot_i"
            xout[ivar]      = Vector{Matrix{Float64}}(undef,nSpecies)
            xout[ivar][1]   = @mlv -model.B*u[1]*min(1e40, exp(-model.RTa/RT))
            xout[ivar][2]   = @mlv -xout[ivar][1]
        elseif vble=="dmdot_ij"
            dmdot_ij            = Matrix{Matrix{Float64}}(undef,nSpecies,nSpecies+3)
            alloc!(dmdot_ij, size(u[1]))
            # dT                = 1/rho/cv * (d(rho e) - e_k d(rho Y_k)) 
            # d(rho e)          = d(rho E) - rhovx/rho * d(rhovx)  
            #                       - rhovy/rho * d(rhovy) 
            #                       + 0.5*((rhovx)^2 + (rhovy)^2)/rho^2 * d(rho)
            #R/rho/cv           = (gamma-1.0)/rho
            #Term due to u[1]
            @mlv dmdot_ij[1,1]  = -model.B*exp(-model.RTa/RT)
            #Term due to RT:
            aux                 = @mlv -model.B*u[1]*exp(-model.RTa/RT)*
                                        model.RTa/RT^2 * (gamma-1.0)/rho
            @mlv dmdot_ij[1,1]  += aux * (0.5*(rhovx*rhovx+rhovy*rhovy)/(rho*rho) -
                                    (hfF + RT/(gamma-1.0)))
            @mlv dmdot_ij[1,2]  = aux * (0.5*(rhovx*rhovx+rhovy*rhovy)/(rho*rho) -
                                    (hfP + RT/(gamma-1.0)))
            @mlv dmdot_ij[1,3]  = - aux * rhovx/rho
            @mlv dmdot_ij[1,4]  = - aux * rhovy/rho
            @mlv dmdot_ij[1,5]  = aux
            for JJ=1:5
                @mlv dmdot_ij[2,JJ]     = - dmdot_ij[1,JJ]
            end
            xout[ivar]          = reshape(dmdot_ij, :)
        elseif vble=="dp_du"
            #From 
            #   rhoE    = p/(gamma-1) + 0.5*((rhovx)^2 + (rhovy)^2)/rho + rhoYF*hfF + rhoYP*hfP,
            #we obtain
            #   0       = 1/(gamma-1) * dp/drhoYk - 0.5*v^2 + hfK
            #   0       = 1/(gamma-1) * dp/drhovx + vx
            #   0       = 1/(gamma-1) * dp/drhovy + vy
            #   1       = 1/(gamma-1) * dp/drhoE 
            dp_dout         = Vector{Matrix{Float64}}(undef, nSpecies+3)
            dp_dout[1]      = @mlv (gamma-1.0)*(0.5*(rhovx*rhovx+rhovy*rho)/(rho*rho)-hfF)
            dp_dout[2]      = @mlv (gamma-1.0)*(0.5*(rhovx*rhovx+rhovy*rho)/(rho*rho)-hfP)
            dp_dout[3]      = @mlv -(gamma-1.0)*rhovx/rho
            dp_dout[4]      = @mlv -(gamma-1.0)*rhovy/rho
            dp_dout[5]      = fill(gamma-1.0, size(u[1]))
            xout[ivar]      = dp_dout
        else
            error("Variable $(vble) not supported")
        end
    end
    
    return xout
    
end

#-------------------------------------------------------------------------------
#MANDATORY FUNCTIONS:

include("../src_CompressibleFlow/CompressibleFlow_fluxes.jl")
include("../src_CompressibleFlow/CompressibleFlow_BC.jl")

#Compute normalization factors from solution. Mass matrix has already been computed.
function nFactsCompute!(solver::SolverData{<:GasModel})

    #L2-norm of density and total energy:
    nVars           = solver.nVars
    rho             = sum(solver.u[1:nVars-3])[:]
    rho_L2          = sqrt( dot(rho, solver.MII, rho)/solver.Omega )
    rhoE            = solver.u[nVars]
    rhoE_L2         = sqrt( dot(rhoE, solver.MII, rhoE)/solver.Omega )
    
    #Normalization factors:
    solver.nFacts[1:nVars-3]        .= rho_L2
    solver.nFacts[nVars-2:nVars-1]  .= sqrt(rho_L2*rhoE_L2)
    solver.nFacts[nVars]            = rhoE_L2
    
    return
    
end

#Function to evaluate flux and source terms at quadrature nodes:
function FluxSource!(model::GasModel, _qp::TrIntVars, ComputeJ::Bool)

    t               = _qp.t
    x               = _qp.x
    u               = _qp.u
    du              = _qp.gradu
    duB             = _qp.graduB
    
    #Compute dependent variables:
    udep            = DepVars(model, t, x, u, model.DepVars)
    
    #Terms due to pressures and convection:
    HyperbolicFlux!(model, u, udep, ComputeJ, _qp.f, _qp.df_du)
    
    #Viscosities:
    epsilon         = @mlv udep[$DepVarIndex(model,"epsilon")][1]
    nu              = @mlv udep[$DepVarIndex(model,"nu")][1]         
    beta            = @mlv udep[$DepVarIndex(model,"beta")][1]       
    kappa_rho_cv    = @mlv udep[$DepVarIndex(model,"kappa_rho_cv")][1]
    
    #Terms due to viscous flux:
    ViscousFlux!(model, nu, beta, u, udep, du, 
        ComputeJ, _qp.f, _qp.df_du, _qp.df_dgradu)
    HeatFlux!(model, kappa_rho_cv, u, udep, du, 
        ComputeJ, _qp.f, _qp.df_du, _qp.df_dgradu)
#     BrennerMassFlux!(model, epsilon, u, udep, du, 
#         ComputeJ, _qp.f, _qp.df_du, _qp.df_dgradu)
    MassDiffusionFlux!(model, u, udep, du, 
        ComputeJ, _qp.f, _qp.df_du, _qp.df_dgradu)
  
    #Monolithic diffusion:
    epsilonFlux!(model, epsilon, du, ComputeJ, _qp.f, _qp.df_dgradu)
    
    #Subgrid stabilization - monolithic diffusion:
    lambda          = udep[DepVarIndex(model,"lambda_max")][1]
#     h_Elems         = _hElems(_qp.Integ2D.mesh)
    A_Elems         = areas(_qp.Integ2D.mesh)
    h_Elems         = @tturbo @. sqrt(A_Elems)
    hp              = h_Elems./_qp.FesOrder * ones(1, _qp.nqp)
    tau             = @mlv model.CSS*lambda*hp
    epsilonFlux!(model, tau, duB, ComputeJ, _qp.fB, _qp.dfB_dgraduB)
    
    #Source terms:
    source!(model, t, x, u, udep, ComputeJ, _qp.Q, _qp.dQ_du)
    
    #CFL number:
    hp_min              = _hmin(_qp.Integ2D.mesh)./_qp.FesOrder * ones(1, _qp.nqp)
    D_max               = @mlv max(epsilon, nu, beta, kappa_rho_cv)
    Deltat_CFL_lambda   = @. $minimum(hp_min/lambda)
    Deltat_CFL_D        = @. $minimum(hp_min^2/D_max)
    _qp.Deltat_CFL      = min(Deltat_CFL_lambda, Deltat_CFL_D)
    
    return
    
end

#bflux! functions are in CompressibleFlow_BC
