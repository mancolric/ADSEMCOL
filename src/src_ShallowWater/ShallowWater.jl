include("../src_LIRKHyp/LIRKHyp.jl")

@warn "tau_char"

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
    
    #Dependent variables:
#     DepVars         ::Vector{String}    = ["h","q1","q2","v1","v2","epsilon","D_penalty"]

    #Mandatory fields:
    nVars           ::Int               = 4     #Number of discretization variables: h, q1, q2, b

end

#Slip condition at a wall:
#   q_n             = 0     zero normal velocity
#   epsilon*deta/dn = 0     zero diff mass flux
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
    b           = u[4]
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
        elseif vble=="v"
            xout[ivar]      = [@tturbo @. sqrt(q1*q1+q2*q2)/h]
        elseif vble=="Fr"
            xout[ivar]      = [@tturbo @. (sqrt(q1*q1+q2*q2))/h / sqrt(model.g*h)]
        elseif vble=="epsilon"
            xout[ivar]      = [fill(model.epsilon, size(u[1]))]
        elseif vble=="eta"
            xout[ivar]      = [eta]
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


include("ShallowWater_fluxes.jl")
include("ShallowWater_BC.jl")

# #Compute normalization factors from solution. Mass matrix has already been computed.
function nFactsCompute!(solver::SolverData{SWE})

    #L2-norm of h:
    nVars           = solver.nVars
    eta             = solver.u[1]
    eta_L2          = sqrt( dot(eta, solver.MII, eta)/solver.Omega )
    h               = solver.u[1]-solver.u[4]
    h_L2            = sqrt( dot(h, solver.MII, h)/solver.Omega )
    
    #Normalization factors:
#     solver.nFacts[1]        = eta_L2
#     solver.nFacts[2:3]      .= h_L2*sqrt(solver.model.g*h_L2)
#     solver.nFacts[4]        = eta_L2
    solver.nFacts[1]        = h_L2
    solver.nFacts[2:3]      .= h_L2*sqrt(solver.model.g*h_L2)
    solver.nFacts[4]        = h_L2
    
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
    b               = u[4]
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
    lambda              = @tturbo @. sqrt(v1*v1 + v2*v2) + sqrt(g*h)
    
    #CFL number:
    Deltat_CFL_lambda   = @tturbo @. $minimum(hp_min/(lambda+1e-12))
    Deltat_CFL_epsilon  = @tturbo @. $minimum(hp_min*hp_min/(model.epsilon+1e-12))
    _qp.Deltat_CFL      = min(Deltat_CFL_lambda, Deltat_CFL_epsilon)
    
    #Compute fluxes and source term:
    HyperbolicFlux!(model, u, ComputeJ, f, df_du)
    #
    epsilon_qp          = @tturbo @. model.epsilon + 0.0*u[1]
    epsilonFlux!(model, epsilon_qp, du, ComputeJ, f, _qp.df_dgradu, IIv=Vector{Int}(1:3))
    #   
    tau_char            = 10.0*_qp.Deltat_CFL
    Source!(model, x, tau_char, u, du, ComputeJ, Q, dQ_du, dQ_du_dx)

    #Subgrid viscosity:
    epsilon_SS          = @tturbo @. model.CSS*lambda*hp
    epsilonFlux!(model, epsilon_SS, duB, ComputeJ, _qp.fB, _qp.dfB_dgraduB, IIv=Vector{Int}(1:3))
    
    return

end

#bflux! functions are in ShallowWater_BC

#-------------------------------------

function WB_LIRKHyp_Step!(solver::SolverData, u_eq::Vector{Float64})

    println()
    t_start     = time()
    
    #flag=1: solver succeeded.
    ConvFlag    = 0
    
    #Get current solution:
    t_n         = solver.t
    u_n         = solver.uv
    etaS_n      = solver.etaS
    etaT_n      = solver.etaT
    etaA_n      = solver.etaA
    Deltat_nm1  = solver.Deltat #previous time step
    
    #Reset RK matrices f_RK and Ju_RK:
    @mlv solver.f_RK        = 0.0
    @mlv solver.Ju_RK       = 0.0
    
    #Target Deltat and estimated etaT_np1:
    Deltat_n            = 0.0
    etaT_np1_est        = 0.0   #this is estimated value. It is necessary to set TolA:
    if solver.t==0.0 #first step
        Deltat_n        = solver.Deltat0    #initial/prescribed Deltat 
    elseif !solver.TimeAdapt
        Deltat_n        = min(solver.Deltat0, 1.8*Deltat_nm1) 
    else
        Deltat_n        = Deltat_nm1 * min(
                            max(0.15, (solver.SfT*solver.TolT/etaT_n)^(1.0/solver.RK.order)),
                            2.0,
                            (solver.LS_iters_target/solver.LS_iters)^0.5)
    end
    
    #Clip Deltat if necessary:
    t_np1                   = t_n+Deltat_n
    if t_np1>solver.tf
        Deltat_n            = solver.tf-t_n
        t_np1               = solver.tf
        etaT_np1_est        = etaT_n * (Deltat_n/Deltat_nm1)^(solver.RK.order)
    end
    
    #Estimate next etaS and etaT:
    etaS_np1_est            = etaS_n * exp(solver.dlogetaS_dt*Deltat_n)
    etaT_np1_est            = etaT_n * (Deltat_n/Deltat_nm1)^(solver.RK.order)
    
    #Adapt mesh before time step:
    TolS_req                = max(solver.SfS*solver.TolS_max * exp(-solver.dlogetaS_dt*Deltat_n),
                                    0.1*solver.SfS*solver.TolS_max)
    t_ini                   = time()
    AMA_flag, AMA_iters     = AdaptMesh!(solver, (solver.u, solver.fes), 
                                solver.TolS_min, min(solver.TolS_max, TolS_req) )
#     AMA_flag, AMA_iters     = AdaptMesh!(solver, (solver.u, solver.fes), 
#                                 solver.TolS_min, min(solver.TolS_max, TolS_req), 
#                                 AMA_RefineFactor=0.5, DEq_MaxIter=10)
    solver.tAMA             += time()-t_ini
    if AMA_flag<0
        printstyled("AMA algorithm failed. Aborting\n", 
                        color=:magenta)
        ConvFlag            = -2+AMA_flag
        return ConvFlag
    end
    
    #Redefine some variables if mesh has changed:
    t_ini                   = time()
    if AMA_iters>0
        
        #Change decrease Deltat by a factor:
        #(eta_n remains unchanged)
        Deltat_n        *= solver.AMA_rDeltat
        t_np1           = t_n + Deltat_n
        etaT_np1_est    = etaT_n * (Deltat_n/Deltat_nm1)^(solver.RK.order)
        
        #Solution and space errors:
        u_n             = solver.uv
        etaS_n          = solver.etaS
        etaS_np1_est    = etaS_n * exp(solver.dlogetaS_dt*Deltat_n)
        
        #RK variables:
        RKAlloc!(solver)
        
    end
    if t_n==0.0
        etaT_np1_est    = solver.TolT
    end
    
    solver.tAlloc       += time()-t_ini
    
    #Set algebraic tolerance:
    TolA                = max(solver.TolA_min, solver.CA*min(etaS_np1_est, etaT_np1_est, etaS_n))
    
    #Well-balancing: compute F(ueq):
    f_eq                    = zeros(solver.nVars*solver.fes.nDof)
    Rhs!(solver, t_n, u_eq, false, f_eq, solver.Jm)
    
    #Flux and Jacobian for first stage:
    t_ini                   = time()
    println("Computing Jacobian")
    Deltat_CFL              = Rhs!(solver, t_n, u_n, true, 
                                view(solver.f_RK,:,1), solver.Jm)
    solver.tJm              += time()-t_ini
    printstyled("Jacobian computed in ", time()-t_ini, " seconds, nnz=", 
        solver.Jm.colptr[solver.Jm.n+1], " \n", color=:cyan)
    view(solver.f_RK,:,1)   .-= f_eq
    view(solver.Ju_RK,:,1)  .= solver.Jm*(u_n-u_eq)

    #---------------------------------------------------------
    #Stages 2,...,s:
    
    #Allocate variables:
    u_k                 = 1.0*u_n
    u_k_pert            = u_n-u_eq
    b_1                 = solver.Mm*u_k_pert
    etaS_np1            = 0.0
    etaS_elems_np1      = zeros(0)
    urec_np1            = Vector{VectorView{Float64}}(undef,solver.nVars)
    etaT_np1            = 0.0
    etaA_np1            = 0.0
    
    #Loop stages until algebraic and time tolerances are met:
    RepeatTA            = true
    LSFailed_iters      = 0
    DCFailed_iters      = 0
    while RepeatTA #before this loop, Deltat_n, t_np1, TolA must be computed
        
        #Loop stages (solution for stage kk is overwritten in solver.t and solver.u):
        LSFlag                  = 0
        solver.LS_iters         = 0
        for kk=2:solver.RK.stages
        
            #Set time:
            t_k                     = t_n + solver.RK.c[kk]*Deltat_n
            
            #Compute system matrix and set up linear system, if necessary.
            if kk==2 || solver.RK.AI[kk,kk]!=solver.RK.AI[kk-1,kk-1]
                println("Setting up linear system")
                t_ini                   = time()
#                 @mlv solver.Am.nzval    = solver.Mm.nzval - 
#                                             (Deltat_n*solver.RK.AI[kk,kk])*solver.Jm.nzval
                BLAS.axpby!(1.0, solver.Mm.nzval, 0.0, solver.Am.nzval)
                BLAS.axpby!(-Deltat_n*solver.RK.AI[kk,kk], solver.Jm.nzval, 1.0, solver.Am.nzval)
                LinearSystem!(solver.Am_LS)
                solver.tSCILU           += time()-t_ini
                printstyled("Linear system set up in ", time()-t_ini, " seconds \n", color=:cyan)
#                 PlotLSPermutations(solver.Am_LS, solver.MII, solver.nVars)
            end
            
            #Compute r.h.s.:
            t_ini                   = time()
            @mlv solver.bv          = b_1
            for ll=1:kk-1
                @mlv solver.bv      += Deltat_n*((solver.RK.AI[kk,ll]-solver.RK.AE[kk,ll])*
                                                $view(solver.Ju_RK,:,ll) + 
                                        solver.RK.AE[kk,ll]*$view(solver.f_RK,:,ll))
            end
            if any(isnan.(solver.bv))
                display(kk)
                display(norm(b_1))
                display(norm(solver.f_RK))
                display(norm(solver.Ju_RK))
                error("NaNs in solver.bv")
            end
            solver.tb               += time()-t_ini
            
            #Save linear system to discuss role of CA:
            t_CA        = 0.8*solver.tf
            if false && solver.t<t_CA && t_np1>=t_CA && kk==solver.RK.stages
                filename    = "$(ResUbi)LIRKHyp_SC$(solver.SC)_LinearSystem.jld2"
                save(filename, "A", solver.Am_LS, "u", u_k, "b", solver.bv, "etaST", TolA/solver.CA, 
                    "solver", save(solver))
            end
            
            #Solve:
            t_ini           = time()
#             BLAS.axpby!(1.0, u_k, 0.0, u_k_pert)
#             BLAS.axpby!(-1.0, u_eq, 1.0, u_k_pert)
            @. u_k_pert     = u_k - u_eq
            #Correct u_k_pert if the latter is very small:
            u_k_pert_norm,  = LqMean(solver.Integ2D, 
                                GetViews(u_k_pert, solver.nVars, solver.fes.nDof), 
                                solver.fes, solver.nFacts, q=solver.SpaceNorm)
#             if u_k_pert_norm<10*TolA+Inf
#                 @. u_k_pert = 0.0
#             end
#             u_k_pert_old    = copy(u_k_pert)
            LSOutput        = LS_gmres!(solver.Am_LS, u_k_pert, solver.bv, AbsTol=TolA, Display="final",
                                MaxIter=solver.LS_iters_max)
            solver.tLS      += time()-t_ini
#             BLAS.axpby!(1.0, u_eq, 0.0, u_k)
#             BLAS.axpby!(1.0, u_k_pert, 1.0, u_k)
            @. u_k          = u_eq + u_k_pert
#             display(norm(u_k_pert-u_k_pert_old,Inf))
            LSFlag          = LSOutput[1]
            LSIter          = LSOutput[2]
            etaA_np1        = LSOutput[3]
            if LSFlag<=0
#                 printstyled("Linear solver did not converge\n", color=:light_yellow)
                break #loop for kk=1:nStages
            else
#                 println("Stage $kk, linear solver converged in $(LSIter) iterations")
            end
            solver.LS_iters     += LSIter
            solver.LS_total     += LSIter
            
            #Compute flux and save derivatives:
            t_ini                   = time()
            Rhs!(solver, t_k, u_k, false, view(solver.f_RK,:,kk), solver.Jm)
            solver.tRhs             += time()-t_ini
#             println("Rhs vector = ", time()-t_ini)
            view(solver.f_RK,:,kk)  .-= f_eq
            view(solver.Ju_RK,:,kk) .= solver.Jm*u_k_pert
            
        end
        
        #If LS converged:
        if LSFlag>0
        
            #Mean number of LS iterations:
            solver.LS_iters     /= solver.RK.stages-1
            
            #Space errors:
            printstyled("Computing space error\n", color=:white)
            t_ini               = time()
            etaS_np1, 
                etaS_elems_np1,
                urec_np1        = etaS_reconstruction(solver, GetViews(u_k, solver.nVars, solver.fes.nDof), 
                                    q=solver.SpaceNorm)
            solver.dlogetaS_dt  = (log(etaS_np1)-log(etaS_n))/Deltat_n
            solver.tetaS        += time()-t_ini
            printstyled("Space error computed in ", time()-t_ini, " seconds\n", color=:cyan)
                        
            #Time errors:
            printstyled("Computing time error\n", color=:white)
            t_ini               = time()
            @mlv solver.bv      = 0.0
            ss                  = solver.RK.stages
            for ll=1:ss
                @mlv solver.bv  += Deltat_n*( ( (solver.RK.AI[ss,ll]-solver.RK.AE[ss,ll])-
                                              (solver.RK.bhatI[ll]-solver.RK.bhatE[ll]) ) *
                                                $view(solver.Ju_RK,:,ll) + 
                                            ( (solver.RK.AE[ss,ll]-solver.RK.bhatE[ll]) * 
                                                $view(solver.f_RK,:,ll) ) )
            end
            ehat                = 0.0*u_k   #ehat=u-uhat
            ehat_views          = GetViews(ehat, solver.nVars, solver.fes.nDof)
            for II=1:solver.nVars
                flag,           = LS_gmres!(solver.MII_LS, ehat_views[II], solver.b[II], 
                                    AbsTol=TolA*solver.nFacts[II], Display="notify")
                if flag<=0
#                     printstyled("Unable to compute solution for embedded RK. Aborting\n", 
#                         color=:magenta)
#                     return -2
                    printstyled("Unable to compute solution for embedded RK\n", 
                        color=:light_yellow)
                end
            end
            etaT_np1,           = LqMean(solver.Integ2D, ehat_views, solver.fes, solver.nFacts, 
                                    q=solver.SpaceNorm)
            solver.tetaT        += time()-t_ini
            printstyled("Time error computed in ", time()-t_ini, " seconds\n", color=:cyan)
                     
            #Save errors:
            push!(solver.tv, t_np1)
            push!(solver.etaSv, etaS_np1)
            push!(solver.etaTv, etaT_np1)
            push!(solver.etaAv, etaA_np1)
            push!(solver.nElemsv, solver.mesh.nElems)
            push!(solver.nDofv, solver.nVars*solver.fes.nDof)
            push!(solver.CFLv, Deltat_n/Deltat_CFL)
            
            #If algebraic error does not satisfy tolerance, solve again:
            #NOTE: We should avoid this as much as possible:
            TolA_required       = max(solver.TolA_min, 
                                    solver.CA_max*min(etaS_n, etaS_np1, etaT_np1))
            if etaA_np1>TolA_required
            
                printstyled("Algebraic error (", sprintf1("%.2e", etaA_np1), 
                    ") larger than required (", sprintf1("%.2e", TolA_required),
                    ") \n", color=:light_yellow)
                    
                #Correct target algebraic tolerance:
                TolA        = max(solver.TolA_min, solver.CA*min(etaS_np1, etaT_np1))
                
                #Loop again with same time step and same mesh:
                RepeatTA    = true
                
                #Save iteration:
                push!(solver.validv, false)
                
            elseif solver.TimeAdapt && etaT_np1>solver.TolT
            
                printstyled("Time error (", sprintf1("%.2e", etaT_np1), 
                    ") larger than required (", sprintf1("%.2e", solver.TolT),
                    ") \n", color=:light_yellow)
                
                #Correct time step:
                Deltat_ratio    = max(0.2, (solver.SfT*solver.TolT/etaT_np1)^(1.0/solver.RK.order))
                Deltat_n        *= Deltat_ratio
                t_np1           = t_n + Deltat_n
                
                #Next iteration:
                RepeatTA        = true
                etaS_np1_est    = etaS_n * exp(solver.dlogetaS_dt * Deltat_n)
                etaT_np1_est    = etaT_np1 * Deltat_ratio^(solver.RK.order)
                TolA            = max(solver.TolA_min, solver.CA*min(etaS_np1_est, etaT_np1_est, etaS_n))
                
                #Save iteration:
                push!(solver.validv, false)
                
            elseif solver.SpaceAdapt && solver.AMA_MaxIter>0 && etaS_np1>solver.TolS_max
            
                printstyled("Space error (", sprintf1("%.2e", etaS_np1), 
                    ") larger than required (", sprintf1("%.2e", solver.TolS_max),
                    ") \n", color=:light_yellow)
                
                #Correct time step:
                Deltat_ratio    = max(0.2, log(solver.SfS*solver.TolS_max/etaS_n)/solver.dlogetaS_dt)
                Deltat_n        *= Deltat_ratio
                t_np1           = t_n + Deltat_n
                
                #Next iteration:
                RepeatTA        = true
                etaS_np1_est    = etaS_n * exp(solver.dlogetaS_dt * Deltat_n)
                etaT_np1_est    = etaT_np1 * Deltat_ratio^(solver.RK.order)
                TolA            = max(solver.TolA_min, solver.CA*min(etaS_np1_est, etaT_np1_est, etaS_n))
                
                #Save iteration:
                push!(solver.validv, false)
            
            #If all the errors satisfy tolerance, terminate iterations:
            else
            
                RepeatTA        = false
                ConvFlag        = 1
                push!(solver.validv, true)
                
                #Update AMA_rDeltat if there has been AMA:
                if AMA_iters>0
                    #Ideal Deltat to have etaT=0.8TolT
                    Deltat_opt  = Deltat_n*(0.8*solver.TolT/etaT_np1)^(1.0/solver.RK.order)
                    #etaT_np1<TolT, so Deltat_opt is not going to be much smaller that Deltat_n:
                    Deltat_opt  = min(Deltat_opt, 1.2*Deltat_n)
                    #Update rDeltat:
                    solver.AMA_rDeltat  = min(1.0, Deltat_opt/Deltat_nm1)
                    printstyled("Update AMA_rDeltat=", sprintf1("%.2f", solver.AMA_rDeltat),", \n", color=:cyan)
                end
                
            end
            
        #LS did not converge. Decrease step size or abort.
        else
        
            #Save errors:
            push!(solver.tv, t_np1)
            push!(solver.etaSv, NaN)
            push!(solver.etaTv, NaN)
            push!(solver.etaAv, etaA_np1)
            push!(solver.nElemsv, solver.mesh.nElems)
            push!(solver.nDofv, solver.nVars*solver.fes.nDof)
            push!(solver.CFLv, Deltat_n/Deltat_CFL)
            push!(solver.validv, 0)
            
            printstyled("Linear solver did not converge. etaA=", sprintf1("%.2e", etaA_np1), 
                ", TolA=", sprintf1("%.2e", TolA), 
                ". Decreasing time step\n", color=:light_yellow)
            if LSFailed_iters==10
                printstyled("Linear solver unable to converge after ",
                    LSFailed_iters, " attempts. Aborting\n", color=:magenta)
                RepeatTA            = false
                ConvFlag            = -1
            else
            
                #Reduce time step.:
                LSFailed_iters      += 1
                rDeltat             = max(0.2, (0.8*TolA/etaA_np1)^(2))
                Deltat_n            *= rDeltat
                t_np1               = t_n + Deltat_n
                
                #Next iteration:
                RepeatTA            = true
                etaS_np1_est        = etaS_n * exp(solver.dlogetaS_dt * Deltat_n)
                etaT_np1_est        = etaT_n * (Deltat_n/Deltat_nm1)^(solver.RK.order)
                if t_n==0.0
                    etaT_np1_est    = solver.TolT
                end
                TolA                = max(solver.TolA_min, 
                                        solver.CA*min(etaS_np1_est, etaT_np1_est, etaS_n))
                                    
            end
            
        end
        
    end
    
    #Time-algebraic loop finished. If method has converged...
    if ConvFlag==1
        
        #Save new solution:
        solver.t            = t_np1
        solver.Nt           += 1
        solver.Deltat       = Deltat_n
        solver.CFL          = Deltat_n/Deltat_CFL
        solver.uv           .= u_k
        solver.etaS         = etaS_np1
        solver.etaS_elems   = etaS_elems_np1
        solver.urec         = urec_np1
        solver.etaT         = etaT_np1
        solver.etaA         = etaA_np1
        solver.etaST        = sqrt( (solver.etaST^2*t_n + 
                                 (etaS_np1 + etaT_np1)^2*Deltat_n)
                                /t_np1 )
        
        println("SC=", sprintf1("%d", solver.SC),
            ", t=", sprintf1("%.2e", solver.t), 
            ", Deltat=", sprintf1("%.2e", solver.Deltat), 
            ", CFL=", sprintf1("%.2e", solver.CFL), 
            ", TotalDof=", sprintf1("%d", solver.nVars*solver.fes.nDof),
            ", MasterDof=", sprintf1("%d", solver.Am_LS.Pl.nMasters),
            ", hp_min=", sprintf1("%.2e", minimum(_hmin(solver.mesh))/solver.FesOrder), 
            ", etaS=", sprintf1("%.2e", solver.etaS), 
            ", etaT=", sprintf1("%.2e", solver.etaT), 
            ", etaA=", sprintf1("%.2e", solver.etaA),
            ", etaST=", sprintf1("%.2e", solver.etaST), 
            ", tCPU=", sprintf1("%.2e", solver.tCPU+time()-t_start))
        
    else
    
        printstyled("LIRKHyp method did not converge. Aborting\n", color=:magenta)
        
    end
    
    #Clear memory:
    if mod(solver.Nt,500)==0
        GC.gc()
    end
    
    solver.tCPU         += time()-t_start 
    
    return ConvFlag
    
end
