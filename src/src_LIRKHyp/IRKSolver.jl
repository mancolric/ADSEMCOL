#-----------------------------------------------------------------
#Implicit RK solver:

function IRK_Step!(solver::SolverData; LS_RelTol::Float64=0.0)

#     return IRK_Step_Pre!(solver; LS_RelTol=LS_RelTol)
    return IRK_Step_Post!(solver; LS_RelTol=LS_RelTol)
    
end

function IRK_Step_Pre!(solver::SolverData; LS_RelTol::Float64=0.0)
    
    println()
    t_start     = time()
    
    #flag=1: solver succeeded.
    ConvFlag    = 0
    
    #Get current solution:
    t_n         = solver.t
    u_n         = solver.uv
    u_n_views   = solver.u
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
    if isnan(TolS_req)  #May happen for the first time step, if etaS=0 and dlogetaS_dt=Inf:
        TolS_req            = 0.1*solver.SfS*solver.TolS_max
    end
    t_ini                   = time()
    AMA_flag, AMA_iters     = AdaptMesh!(solver, (solver.u, solver.fes), 
                                solver.TolS_min, min(solver.TolS_max, TolS_req) )
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
        u_n_views       = solver.u
        etaS_n          = solver.etaS
        etaS_np1_est    = etaS_n * exp(solver.dlogetaS_dt*Deltat_n)
        
        #RK variables:
        RKAlloc!(solver)
        
    end
    if t_n==0.0
        etaT_np1_est    = 0.1*solver.TolT
    end
    
    solver.tAlloc       += time()-t_ini
    
    #Set algebraic tolerance:
    TolA            = max(solver.TolA_min, solver.CA*min(etaS_np1_est, etaT_np1_est, etaS_n))
    
    #Flux and Jacobian for first stage:
    t_ini                   = time()
    println("Computing Jacobian")
    flag, Deltat_CFL        = Rhs!(solver, t_n, u_n, true, 
                                view(solver.f_RK,:,1), solver.Jm)
    solver.tJm              += time()-t_ini
    printstyled("Jacobian computed in ", time()-t_ini, " seconds, nnz=", 
        solver.Jm.colptr[solver.Jm.n+1], " \n", color=:cyan)
    
    #---------------------------------------------------------
    #Stages 2,...,s:
    
    #Allocate variables:
    u_k                 = zeros(size(u_n))
    etaS_np1            = 0.0
    etaS_elems_np1      = zeros(0)
    urec_np1            = Vector{VectorView{Float64}}(undef,solver.nVars)
    etaT_np1            = 0.0
    etaA_np1            = Inf
    
    #Save M*u[1]:
    b1                  = MuProduct(solver.MII, u_n, solver.nVars)

    #Scaling vector:
    scalv               = zeros(solver.nVars*solver.fes.nDof)
    scalv_views         = GetViews(scalv, solver.nVars, solver.fes.nDof)
    for II=1:solver.nVars
        scalv_views[II] .= solver.nFacts[II]
    end
    
    #Loop stages until algebraic and time tolerances are met:
    RepeatTA            = true
    LSFailed_iters      = 0
    DCFailed_iters      = 0
    Am_F                = nothing
    while RepeatTA #before this loop, Deltat_n, t_np1, TolA must be computed
    
        #First stage:
        @mlv u_k        = u_n
        
        #Loop stages (solution for stage kk is overwritten in solver.t and solver.u):
        LSFlag          = 0
        solver.LS_iters = 0
        for kk=2:solver.RK.stages
        
            #Set time:
            t_k         = t_n + solver.RK.c[kk]*Deltat_n
            
            #Compute system matrix and set up linear system, if necessary.
            if kk==2 || solver.RK.AI[kk,kk]!=solver.RK.AI[kk-1,kk-1]
                println("Setting up linear system")
                t_ini                   = time()
                BLAS.axpby!(1.0, solver.Mm.nzval, 0.0, solver.Am.nzval)
                BLAS.axpby!(-Deltat_n*solver.RK.AI[kk,kk], solver.Jm.nzval, 1.0, solver.Am.nzval)
                LinearSystem!(solver.Am_LS.LS)                
                solver.tSCILU           += time()-t_ini
                printstyled("Linear system set up in ", time()-t_ini, " seconds \n", color=:cyan)
#                 PlotLSPermutations(solver.Am_LS, solver.MII, solver.nVars)
            end
            
            #Compute r.h.s.:
            t_ini                   = time()
            @mlv solver.bv          = 0.0
            solver.bv               .= b1
            for ll=1:kk-1
                @mlv solver.bv      += (Deltat_n*solver.RK.AI[kk,ll])*$view(solver.f_RK,:,ll)
            end
            solver.tb               += time()-t_ini
            
            #Solve nonlinear system
            #   Mu - b - a_ii*Deltat_n*f(t_k, u_k) = 0
            #The nonlinear system is preconditioned with Mhat - alpha*tau*Jhat, where
            #Jhat is possibly a more sparse approximation to J.
            #Also, the preconditioner may be just an approximate factorization or some iterations
            #of preconditioned GMRES:
            if any(uppercase(solver.LSType).==["SCILU0_GMRES", "BLOCKJACOBI", 
                "ILU0", "ILU0_GMRES"])
            
                #Define preconditioned residual:
                function QNResidualJ!(uhat::Vector{Float64}, gres::Vector{Float64})
                    
                    t_ini2      = time()
                    
                    #Unscale u:
                    u           = @tturbo @. scalv*uhat
                    
#                     println("Unscale: ", time()-t_ini2)
                    
                    #Compute flux and save derivatives:
    #                 t_ini                   = time()
                    flag,       = Rhs!(solver, t_k, u, false, 
                                    view(solver.f_RK,:,kk), solver.Jm)
    #                 solver.tRhs             += time()-t_ini
                    if flag<0
                        return flag
                    end
                    
#                     println("Rhs: ", time()-t_ini2) #BOTTLENECK
                    
                    #Compute residual f:
                    #   f = M*u - b - Deltat*a_kk*f(t_k,u_k)
                    fres        = MuProduct(solver.MII, u, solver.nVars) - 
                                        solver.bv - (Deltat_n*solver.RK.AI[kk,kk])*
                                        view(solver.f_RK,:,kk)
                                        
#                     println("fres: ", time()-t_ini2)
                    
                    #Compute preconditioned residual. Since this equation is linear and 
                    #f does not need to be computed at each iteration, we ask for a very
                    #small tolerance. In that case, having an accurate Newton step gres will 
                    #allow for a reduction in the number of Newton steps (and evaluations of f):
                    LSOutput    = LS_solve!(solver.Am_LS.LS, gres, fres, 
                                    RelTol=LS_RelTol, AbsTol=TolA, 
                                    Display="notify", MaxIter=solver.LS_iters_max)
                                    
                    #Returns flag, nIters, etaA
                    flag        = LSOutput[1]
                    
                    #Scale gres:
                    @tturbo @. gres     /= scalv
                    
#                     println("gres: ", time()-t_ini2) #BOTTLENECK
                    
                    return flag
                    
                end
            
                #Solve:
                t_ini       = time()
                LSOutput    = Anderson(FW_NLS((uhat,gres)->QNResidualJ!(uhat,gres)), 
                                u_k./scalv, 
                                UpdateType=2, 
                                AbsTolX=1.0*TolA, RelTolX=0.0, 
                                AbsTolG=0.0*TolA, RelTolG=0.0, 
                                NormFun=FW_NLS_norm((x)->norm(x)/sqrt(length(u_k))), 
                                memory=100, MaxIter=solver.LS_iters_max, 
                                history=true, Display="final")
                solver.tLS  += time()-t_ini
                u_k         .= LSOutput[1].*scalv
                LSFlag      = LSOutput[2].flag
                LSIter      = LSOutput[2].nIter
                etaA_np1    = LSOutput[2].pnorms[ length(LSOutput[2].pnorms) ]
                if LSFlag<0
                    #NLS did not converge because something went wrong with the residual function
                    etaA_np1        = Inf   #In this way, time step is decreased as maximum as possible
                    break #loop for kk=2:nStages
                elseif LSFlag==0
                    #NLS reached maximum of iterations
                    break #loop for kk=2:nStages
                else
                end
                solver.LS_iters     += LSIter
                solver.LS_total     += LSIter
                #Here, LS in reality is NLS
                
            else
            
                error("Unknown method $(solver.LSType)")
                
            end #krylov if
            
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
                @mlv solver.bv  += Deltat_n*(solver.RK.AI[ss,ll]-solver.RK.bhatI[ll])*
                                    $view(solver.f_RK,:,ll)
            end
            ehat                = 0.0*u_k   #ehat=u-uhat
            ehat_views          = GetViews(ehat, solver.nVars, solver.fes.nDof)
            for II=1:solver.nVars
                flag,           = LS_solve!(solver.MII_LS, ehat_views[II], solver.b[II], 
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
                TolA        = max(solver.TolA_min, solver.CA*min(etaS_n, etaS_np1, etaT_np1))
#                 TolA        = TolA_required
                
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
                if isnan(TolA) #happens if etaS_n=0.0 (initial condition = 0 and dlogetaS_dt=0)
                    TolA        = solver.TolA_min
                end
                
                #Save iteration:
                push!(solver.validv, false)
                
            elseif solver.SpaceAdapt && solver.AMA_MaxIter>0 && etaS_np1>solver.TolS_max
            
                printstyled("Space error (", sprintf1("%.2e", etaS_np1), 
                    ") larger than required (", sprintf1("%.2e", solver.TolS_max),
                    ") \n", color=:light_yellow)
                
                #Correct time step:
                Deltat_ratio    = max(0.2, log(solver.SfS*solver.TolS_max/etaS_n)/solver.dlogetaS_dt)
                if isnan(Deltat_ratio) #happens if etaS_0=0 (initial condition = 0 and dlogetaS_dt=0)
                    Deltat_ratio    = 0.2
                end
                Deltat_n        *= Deltat_ratio
                t_np1           = t_n + Deltat_n
                
                #Next iteration:
                RepeatTA        = true
                etaS_np1_est    = etaS_n * exp(solver.dlogetaS_dt * Deltat_n)
                etaT_np1_est    = etaT_np1 * Deltat_ratio^(solver.RK.order)
                TolA            = max(solver.TolA_min, solver.CA*min(etaS_np1_est, etaT_np1_est, etaS_n))
                if isnan(TolA) #happens if etaS_0=0 (initial condition = 0 and dlogetaS_dt=0)
                    TolA        = solver.TolA_min
                end
                
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
            
            printstyled("Nonlinear solver did not converge. etaA=", sprintf1("%.2e", etaA_np1), 
                ", TolA=", sprintf1("%.2e", TolA), 
                ". Decreasing time step\n", color=:light_yellow)
            if LSFailed_iters==10
                printstyled("Nonlinear solver unable to converge after ",
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
            ", MasterDof=", sprintf1("%d", solver.Am_LS.LS.Pl.nMasters),
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
    if mod(solver.Nt,50)==0
        GC.gc()
    end
    
    solver.tCPU         += time()-t_start 
    
    return ConvFlag
    
end

function IRK_UpdatedJ_Step_Pre!(solver::SolverData)
    
    println()
    t_start     = time()
    
    #flag=1: solver succeeded.
    ConvFlag    = 0
    
    #Get current solution:
    t_n         = solver.t
    u_n         = solver.uv
    u_n_views   = solver.u
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
    if isnan(TolS_req)  #May happen for the first time step, if etaS=0 and dlogetaS_dt=Inf:
        TolS_req            = 0.1*solver.SfS*solver.TolS_max
    end
    t_ini                   = time()
    AMA_flag, AMA_iters     = AdaptMesh!(solver, (solver.u, solver.fes), 
                                solver.TolS_min, min(solver.TolS_max, TolS_req) )
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
        u_n_views       = solver.u
        etaS_n          = solver.etaS
        etaS_np1_est    = etaS_n * exp(solver.dlogetaS_dt*Deltat_n)
        
        #RK variables:
        RKAlloc!(solver)
        
    end
    if t_n==0.0
        etaT_np1_est    = 0.1*solver.TolT
    end
    
    solver.tAlloc       += time()-t_ini
    
    #Set algebraic tolerance:
    TolA                = max(solver.TolA_min, solver.CA*min(etaS_np1_est, etaT_np1_est, etaS_n))
    
    #Flux for first stage:
    t_ini               = time()
    flag, Deltat_CFL    = Rhs!(solver, t_n, u_n, false, 
                                view(solver.f_RK,:,1), solver.Jm)
    solver.tRhs         += time()-t_ini
    
    #---------------------------------------------------------
    #Stages 2,...,s:
    
    #Allocate variables:
    u_k                 = zeros(size(u_n))
    etaS_np1            = 0.0
    etaS_elems_np1      = zeros(0)
    urec_np1            = Vector{VectorView{Float64}}(undef,solver.nVars)
    etaT_np1            = 0.0
    etaA_np1            = Inf
    
    #Save M*u[1]:
    b1                  = MuProduct(solver.MII, u_n, solver.nVars)

    #Scaling vector:
    scalv               = zeros(solver.nVars*solver.fes.nDof)
    scalv_views         = GetViews(scalv, solver.nVars, solver.fes.nDof)
    for II=1:solver.nVars
        scalv_views[II] .= solver.nFacts[II]
    end
    
    #Loop stages until algebraic and time tolerances are met:
    RepeatTA            = true
    LSFailed_iters      = 0
    DCFailed_iters      = 0
    Am_F                = nothing
    while RepeatTA #before this loop, Deltat_n, t_np1, TolA must be computed
    
        #First stage:
        @mlv u_k        = u_n
        
        #Loop stages (solution for stage kk is overwritten in solver.t and solver.u):
        LSFlag          = 0
        solver.LS_iters = 0
        for kk=2:solver.RK.stages
        
            #Set time:
            t_k         = t_n + solver.RK.c[kk]*Deltat_n
            
            #Compute r.h.s.:
            t_ini                   = time()
            @mlv solver.bv          = 0.0
            solver.bv               .= b1
            for ll=1:kk-1
                @mlv solver.bv      += (Deltat_n*solver.RK.AI[kk,ll])*$view(solver.f_RK,:,ll)
            end
            solver.tb               += time()-t_ini
            
            #Solve nonlinear system
            #   Mu - b - a_ii*Deltat_n*f(t_k, u_k) = 0
            if any(uppercase(solver.LSType).==["SCILU0", "BLOCKJACOBI", 
                "DOMINANTILU0", "DOMINANT"])
            
                #Define preconditioned residual:
                function QNResidualJ!(uhat::Vector{Float64}, gres::Vector{Float64})
                    
                    t_ini2      = time()
                    
                    #Unscale u:
                    u           = @tturbo @. scalv*uhat
                    
#                     println("Unscale: ", time()-t_ini2)
                    
                    #Compute flux and Jacobian. Save flux directly in f_RK:
                    flag,       = Rhs!(solver, t_k, u, true, 
                                    view(solver.f_RK,:,kk), solver.Jm)
                    if flag<0
                        return flag
                    end
                    
#                     println("Rhs: ", time()-t_ini2) #BOTTLENECK

                    #
                    #Compute system matrix and set up linear system, if necessary.
                    BLAS.axpby!(1.0, solver.Mm.nzval, 0.0, solver.Am.nzval)
                    BLAS.axpby!(-Deltat_n*solver.RK.AI[kk,kk], solver.Jm.nzval, 1.0, solver.Am.nzval)
                    LinearSystem!(solver.Am_LS.LS)                
#                     println("Linear system: ", time()-t_ini2)
                    
                    #Compute residual f:
                    #   f = M*u - b - Deltat*a_kk*f(t_k,u_k)
                    fres        = MuProduct(solver.MII, u, solver.nVars) - 
                                        solver.bv - (Deltat_n*solver.RK.AI[kk,kk])*
                                        view(solver.f_RK,:,kk)
                                        
#                     println("fres: ", time()-t_ini2)
                    
                    #Compute preconditioned residual. Since this equation is linear and 
                    #f does not need to be computed at each iteration, we ask for a very
                    #small tolerance. In that case, having an accurate Newton step gres will 
                    #allow for a reduction in the number of Newton steps (and evaluations of f):
                    flag        = NaN
                    if false && uppercase(solver.LSType)=="SCILU0"
                        @warn "DEBUG: Direct solver"
                        flag    = 1
                        gres    .= solver.Am\fres
                    else
                        LSOutput    = LS_solve!(solver.Am_LS.LS, gres, fres, 
                                    RelTol=0e-3, AbsTol=0.1*TolA, 
                                    Display="notify", MaxIter=50)
                        #Returns flag, nIters, etaA
                        flag        = LSOutput[1]
                    end
                    
                    
                    #Scale gres:
                    @tturbo @. gres     /= scalv
                    
#                     println("gres: ", time()-t_ini2) #BOTTLENECK
                    
                    return flag
                    
                end
            
                #Solve:
                t_ini       = time()
                NLS_mem     = 100
                if true && uppercase(solver.LSType)=="SCILU0" 
                    #J(u)\f(u) is solved exactly, so this is the NR method.
                    NLS_mem = 0
                end
                LSOutput    = Anderson(FW_NLS((uhat,gres)->QNResidualJ!(uhat,gres)), 
                                u_k./scalv, 
                                UpdateType=2, 
                                AbsTolX=1.0*TolA, RelTolX=0.0, 
                                AbsTolG=0.0*TolA, RelTolG=0.0, 
                                NormFun=FW_NLS_norm((x)->norm(x)/sqrt(length(u_k))), 
                                memory=NLS_mem, MaxIter=solver.LS_iters_max, 
                                history=true, Display="final")
                solver.tLS  += time()-t_ini
                u_k         .= LSOutput[1].*scalv
                LSFlag      = LSOutput[2].flag
                LSIter      = LSOutput[2].nIter
                etaA_np1    = LSOutput[2].pnorms[ length(LSOutput[2].pnorms) ]
                if LSFlag<0
                    #NLS did not converge because something went wrong with the residual function
                    etaA_np1        = Inf   #In this way, time step is decreased as maximum as possible
                    break #loop for kk=2:nStages
                elseif LSFlag==0
                    #NLS reached maximum of iterations
                    break #loop for kk=2:nStages
                else
                end
                solver.LS_iters     += LSIter
                solver.LS_total     += LSIter
                #Here, LS in reality is NLS
                
            else
            
                error("Unknown method $(solver.LSType)")
                
            end #krylov if
            
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
                @mlv solver.bv  += Deltat_n*(solver.RK.AI[ss,ll]-solver.RK.bhatI[ll])*
                                    $view(solver.f_RK,:,ll)
            end
            ehat                = 0.0*u_k   #ehat=u-uhat
            ehat_views          = GetViews(ehat, solver.nVars, solver.fes.nDof)
            for II=1:solver.nVars
                flag,           = LS_solve!(solver.MII_LS, ehat_views[II], solver.b[II], 
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
#                 TolA        = max(solver.TolA_min, solver.CA*min(etaS_np1, etaT_np1))
                TolA        = TolA_required
                
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
                if isnan(TolA) #happens if etaS_n=0.0 (initial condition = 0 and dlogetaS_dt=0)
                    TolA        = solver.TolA_min
                end
                
                #Save iteration:
                push!(solver.validv, false)
                
            elseif solver.SpaceAdapt && solver.AMA_MaxIter>0 && etaS_np1>solver.TolS_max
            
                printstyled("Space error (", sprintf1("%.2e", etaS_np1), 
                    ") larger than required (", sprintf1("%.2e", solver.TolS_max),
                    ") \n", color=:light_yellow)
                
                #Correct time step:
                Deltat_ratio    = max(0.2, log(solver.SfS*solver.TolS_max/etaS_n)/solver.dlogetaS_dt)
                if isnan(Deltat_ratio) #happens if etaS_0=0 (initial condition = 0 and dlogetaS_dt=0)
                    Deltat_ratio    = 0.2
                end
                Deltat_n        *= Deltat_ratio
                t_np1           = t_n + Deltat_n
                
                #Next iteration:
                RepeatTA        = true
                etaS_np1_est    = etaS_n * exp(solver.dlogetaS_dt * Deltat_n)
                etaT_np1_est    = etaT_np1 * Deltat_ratio^(solver.RK.order)
                TolA            = max(solver.TolA_min, solver.CA*min(etaS_np1_est, etaT_np1_est, etaS_n))
                if isnan(TolA) #happens if etaS_0=0 (initial condition = 0 and dlogetaS_dt=0)
                    TolA        = solver.TolA_min
                end
                
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
            
            printstyled("Nonlinear solver did not converge. etaA=", sprintf1("%.2e", etaA_np1), 
                ", TolA=", sprintf1("%.2e", TolA), 
                ". Decreasing time step\n", color=:light_yellow)
            if LSFailed_iters==10
                printstyled("Nonlinear solver unable to converge after ",
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
            ", MasterDof=", sprintf1("%d", solver.Am_LS.LS.Pl.nMasters),
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
    if mod(solver.Nt,50)==0
        GC.gc()
    end
    
    solver.tCPU         += time()-t_start 
    
    return ConvFlag
    
end

function IRK_Step_Post!(solver::SolverData; LS_RelTol::Float64=0.0)
    
    println()
    t_start     = time()
    
    #flag=1: solver succeeded.
    ConvFlag    = 0
    
    #Get current solution. 
    #CAREFUL:
    #   u_n, fes_n, *_n is the solution BEFORE mesh adaptation.
    #   solver.u, solver.fes is the solution AFTER.
    t_n         = solver.t
    u_n         = solver.uv
    fes_n       = solver.fes
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
        Deltat_n        = Deltat_nm1 * max(0.15, min(2.0, 
                            (solver.SfT*solver.TolT/etaT_n)^(1.0/solver.RK.order),
                            (solver.LS_iters_target/solver.LS_iters)^0.5))
    end
    
    #Clip Deltat if necessary:
    t_np1                   = t_n+Deltat_n
    if t_np1>solver.tf
        Deltat_n            = solver.tf-t_n
        t_np1               = solver.tf
        etaT_np1_est        = etaT_n * (Deltat_n/Deltat_nm1)^(solver.RK.order)
    end
    
    #Correct dlogetaS_dt if etaS=0.0:
    if !isfinite(solver.dlogetaS_dt)
        solver.dlogetaS_dt  = 0.0
    end
    
    #Estimate next etaS and etaT:
    etaS_np1_est            = etaS_n * exp(solver.dlogetaS_dt*Deltat_n)
    etaT_np1_est            = etaT_n * (Deltat_n/Deltat_nm1)^(solver.RK.order)
    
    #Set algebraic tolerance:
    TolA            = max(solver.TolA_min, solver.CA*min(etaS_np1_est, etaT_np1_est, etaS_n))
    
    #Declare some variables:
    Deltat_CFL      = 0.0
    u_k             = zeros(0)
    urec_np1        = Vector{VectorView{Float64}}(undef,solver.nVars)
    etaS_np1        = 0.0
    etaS_elems_np1  = zeros(0)
    etaT_np1        = 0.0
    etaA_np1        = 0.0
    urec_n          = Vector{VectorView{Float64}}(undef,solver.nVars)
    etaS_n          = 0.0
    etaS_elems_n    = zeros(0)
    
    #Loop until space tolerance is achieved:
    RepeatS         = true
    AMA_iters       = 0
    AMA_RefineFact  = 0.9
    AMA_CoarseFact  = 5.0
    AMA_TolS        = NaN
    AMA_Strategy    = -1 #Undefined
    AMA_factor      = NaN
    while RepeatS
    
        #Flux and Jacobian for first stage:
        t_ini                   = time()
        println("Computing Jacobian...")
        flag, Deltat_CFL        = Rhs!(solver, t_n, solver.uv, true, 
                                    view(solver.f_RK,:,1), solver.Jm)
        solver.tJm              += time()-t_ini
        printstyled("Jacobian computed in ", time()-t_ini, " seconds, nnz=", 
            solver.Jm.colptr[solver.Jm.n+1], " \n", color=:cyan)
    
        #---------------------------------------------------------
        #Stages 2,...,s:
        
        #Allocate variables:
        u_k                 = zeros(size(solver.uv))
#         etaS_np1            = 0.0
#         etaS_elems_np1      = zeros(0)
#         urec_np1            = Vector{VectorView{Float64}}(undef,solver.nVars)
#         etaT_np1            = 0.0
#         etaA_np1            = Inf
        
        #Save M*u[1]:
        b1                  = MuProduct(solver.MII, solver.uv, solver.nVars)

        #Scaling vector:
        scalv               = zeros(solver.nVars*solver.fes.nDof)
        scalv_views         = GetViews(scalv, solver.nVars, solver.fes.nDof)
        for II=1:solver.nVars
            scalv_views[II] .= solver.nFacts[II]
        end
        
        #Loop stages until algebraic and time tolerances are met:
        RepeatTA            = true
        LSFailed_iters      = 0
        DCFailed_iters      = 0
        while RepeatTA #before this loop, Deltat_n, t_np1, TolA must be computed
    
            #First stage:
            @mlv u_k        = solver.uv
        
            #Loop stages (solution for stage kk is overwritten in solver.t and solver.u):
            LSFlag          = 0
            solver.LS_iters = 0
            for kk=2:solver.RK.stages
            
                #Set time:
                t_k         = t_n + solver.RK.c[kk]*Deltat_n
                
                #Compute system matrix and set up linear system, if necessary.
                if kk==2 || solver.RK.AI[kk,kk]!=solver.RK.AI[kk-1,kk-1]
                    println("Setting up linear system...")
                    t_ini                   = time()
                    BLAS.axpby!(1.0, solver.Mm.nzval, 0.0, solver.Am.nzval)
                    BLAS.axpby!(-Deltat_n*solver.RK.AI[kk,kk], solver.Jm.nzval, 1.0, solver.Am.nzval)
                    LinearSystem!(solver.Am_LS.LS)                
                    solver.tSCILU           += time()-t_ini
                    printstyled("Linear system set up in ", time()-t_ini, " seconds \n", color=:cyan)
    #                 PlotLSPermutations(solver.Am_LS, solver.MII, solver.nVars)
                end
            
                #Compute r.h.s.:
                t_ini                   = time()
                @mlv solver.bv          = 0.0
                solver.bv               .= b1
                for ll=1:kk-1
                    @mlv solver.bv      += (Deltat_n*solver.RK.AI[kk,ll])*$view(solver.f_RK,:,ll)
                end
                solver.tb               += time()-t_ini
            
                #Solve nonlinear system
                #   Mu - b - a_ii*Deltat_n*f(t_k, u_k) = 0
                #The nonlinear system is preconditioned with Mhat - alpha*tau*Jhat, where
                #Mhat and Jhat are possibly more sparse approximations to M and J.
                #Also, the preconditioner may be just an approximate factorization or some iterations
                #of preconditioned GMRES:
            
                #Define preconditioned residual:
                function QNResidualJ!(uhat::Vector{Float64}, gres::Vector{Float64})
                    
                    t_ini2      = time()
                    
                    #Unscale u:
                    u           = @tturbo @. scalv*uhat
                    
#                     println("Unscale: ", time()-t_ini2)
                    
                    #Compute flux and save derivatives:
    #                 t_ini                   = time()
                    flag,       = Rhs!(solver, t_k, u, false, 
                                    view(solver.f_RK,:,kk), solver.Jm)
    #                 solver.tRhs             += time()-t_ini
                    if flag<0
                        return flag
                    end
                    
#                     println("Rhs: ", time()-t_ini2) #BOTTLENECK
                    
                    #Compute residual f:
                    #   f = M*u - b - Deltat*a_kk*f(t_k,u_k)
                    fres        = MuProduct(solver.MII, u, solver.nVars) - 
                                        solver.bv - (Deltat_n*solver.RK.AI[kk,kk])*
                                        view(solver.f_RK,:,kk)
                                        
#                     println("fres: ", time()-t_ini2)
                    
                    #Compute preconditioned residual. Since this equation is linear and 
                    #f does not need to be computed at each iteration, we ask for a very
                    #small tolerance. In that case, having an accurate Newton step gres will 
                    #allow for a reduction in the number of Newton steps (and evaluations of f):
                    LSOutput    = LS_solve!(solver.Am_LS.LS, gres, fres, 
                                    RelTol=LS_RelTol, AbsTol=TolA, 
                                    Display="notify", MaxIter=solver.LS_iters_max)
                                    
                    #Returns flag, nIters, etaA
                    flag        = LSOutput[1]
                    
                    #Scale gres:
                    @tturbo @. gres     /= scalv
                    
#                     println("gres: ", time()-t_ini2) #BOTTLENECK
                    
                    return flag
                    
                end
            
                #Solve:
                t_ini       = time()
                LSOutput    = Anderson(FW_NLS((uhat,gres)->QNResidualJ!(uhat,gres)), 
                                u_k./scalv, 
                                UpdateType=2, 
                                AbsTolX=1.0*TolA, RelTolX=0.0, 
                                AbsTolG=0.0*TolA, RelTolG=0.0, 
                                NormFun=FW_NLS_norm((x)->norm(x)/sqrt(length(u_k))), 
                                memory=100, MaxIter=solver.LS_iters_max, 
                                history=true, Display="final")
                solver.tLS  += time()-t_ini
                u_k         .= LSOutput[1].*scalv
                LSFlag      = LSOutput[2].flag
                LSIter      = LSOutput[2].nIter
                etaA_np1    = LSOutput[2].pnorms[ length(LSOutput[2].pnorms) ]
                if LSFlag<0
                    #NLS did not converge because something went wrong with the residual function
                    etaA_np1        = Inf   #In this way, time step is decreased as maximum as possible
                    break #loop for kk=2:nStages
                elseif LSFlag==0
                    #NLS reached maximum of iterations
                    break #loop for kk=2:nStages
                else
                end
                solver.LS_iters     += LSIter
                solver.LS_total     += LSIter
                #Here, LS in reality is NLS
                
            end #stages loop
        
            #If LS converged:
            if LSFlag>0
            
                #Mean number of LS iterations:
                solver.LS_iters     /= solver.RK.stages-1
                
                #Space errors at t^n-:
                if AMA_iters==0
                    etaS_n                          = solver.etaS
                    etaS_elems_n                    = solver.etaS_elems
                    urec_n                          = solver.urec
                else 
                    t_ini                           = time()
                    etaS_n, etaS_elems_n, 
                        urec_n                      = etaS_reconstruction(solver, 
                                                        GetViews(solver.uv, solver.nVars, solver.fes.nDof), 
                                                        q=solver.SpaceNorm)
                    solver.tAMA                     += time()-t_ini 
                    #(in the pre- approach, this reconstruction is performed in AdaptMesh! and the time goes to
                    #tAMA, not to tetaS.)
                end
        
                #Space errors at t^(n+1)-:
#                 printstyled("Computing space error\n", color=:white)
                t_ini               = time()
                etaS_np1, 
                    etaS_elems_np1,
                    urec_np1        = etaS_reconstruction(solver, GetViews(u_k, solver.nVars, solver.fes.nDof), 
                                        q=solver.SpaceNorm)
                solver.dlogetaS_dt  = (log(etaS_np1)-log(etaS_n))/Deltat_n
                if !isfinite(solver.dlogetaS_dt)
                    solver.dlogetaS_dt  = 0.0   #happens if etaS=0 at two consecutive time levels
                end
                solver.tetaS        += time()-t_ini
#                 printstyled("Space error computed in ", time()-t_ini, " seconds\n", color=:cyan)
                            
                #Deactivate mesh adaptation if:
#                 etaS_n_np1          = max(etaS_n, etaS_np1)
                etaS_n_np1          = etaS_np1
                if AMA_Strategy==-1 && solver.TolS_min <= etaS_n_np1 <= solver.TolS_max
                    RepeatS         = false
                elseif AMA_Strategy==0 && solver.TolS_min <= etaS_n_np1 <= solver.TolS_max
                    RepeatS         = false
                elseif AMA_Strategy==1 && etaS_n_np1 <= solver.TolS_max
                    RepeatS         = false
                end
                    
                #Time errors:
#                 printstyled("Computing time error\n", color=:white)
                t_ini               = time()
                @mlv solver.bv      = 0.0
                ss                  = solver.RK.stages
                for ll=1:ss
                    @mlv solver.bv  += Deltat_n*(solver.RK.AI[ss,ll]-solver.RK.bhatI[ll])*
                                        $view(solver.f_RK,:,ll)
                end
                ehat                = 0.0*u_k   #ehat=u-uhat
                ehat_views          = GetViews(ehat, solver.nVars, solver.fes.nDof)
                for II=1:solver.nVars
                    flag,           = LS_solve!(solver.MII_LS, ehat_views[II], solver.b[II], 
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
#                 printstyled("Time error computed in ", time()-t_ini, " seconds\n", color=:cyan)
                        
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
                
                    #Save iteration and print info:
                    push!(solver.validv, false)
                    printstyled("Algebraic error (", sprintf1("%.2e", etaA_np1), 
                        ") larger than required (", sprintf1("%.2e", TolA_required),
                        ") \n", color=:light_yellow)
                        
                    #Correct target algebraic tolerance:
                    TolA        = max(solver.TolA_min, solver.CA*min(etaS_n, etaS_np1, etaT_np1))
#                     TolA        = TolA_required
                    
                    #Loop again with same time step and same mesh:
                    RepeatTA    = true
                    RepeatS     = true
                
                elseif solver.TimeAdapt && etaT_np1>solver.TolT
                
                    #Save iteration and print info:
                    push!(solver.validv, false)
                    printstyled("Time error (", sprintf1("%.2e", etaT_np1), 
                        ") larger than required (", sprintf1("%.2e", solver.TolT),
                        ") \n", color=:light_yellow)
                    
                    #Correct time step:
                    Deltat_ratio    = max(0.2, (solver.SfT*solver.TolT/etaT_np1)^(1.0/solver.RK.order))
                    Deltat_n        *= Deltat_ratio
                    t_np1           = t_n + Deltat_n
                    
                    #Next iteration:
                    RepeatTA        = true
                    RepeatS         = true
                    etaS_np1_est    = etaS_n * exp(solver.dlogetaS_dt * Deltat_n)
                    etaT_np1_est    = etaT_np1 * Deltat_ratio^(solver.RK.order)
                    TolA            = max(solver.TolA_min, solver.CA*min(etaS_np1_est, etaT_np1_est, etaS_n))
                
                elseif solver.SpaceAdapt && solver.AMA_MaxIter>0 && RepeatS
                
                    #Save iteration and print info:
                    push!(solver.validv, false)
                    printstyled("Space error (", sprintf1("%.2e", etaS_np1), 
                        ") larger than required (", sprintf1("%.2e", solver.TolS_max),
                        ") \n", color=:light_yellow)
                        
                    #Initialize AMA_TolS with estimated error:
                    if AMA_iters==0 #first value
                        AMA_TolS        = etaS_n_np1
                        AMA_Strategy    = 0     #coarse mesh
                        AMA_factor      = AMA_CoarseFact
                    end
                    
                    #Change strategy if necessary:
                    if etaS_n_np1>solver.TolS_max
                        AMA_Strategy    = 1 #refine mesh from now on
                        AMA_factor      = AMA_RefineFact
                    end
                    
                    #Exit if maximum iteration number has been reached:
                    if AMA_iters==solver.AMA_MaxIter
                        printstyled("Reached maximum mesh adaptation iterations\n", color=:magenta)
                        return -3
                    end
                    
                    #Update iteration number:
                    AMA_iters           += 1
                    
                    #-----------------------------------------------------------------
                    #Optimal metric for solution at t^(n+1)-:
                    
                    t_ini                   = time()
                    
                    #Compute optimal anisotropy:
                    emeanv0, emeanv, 
                        err0, err1, 
                        metric_np1          = AMA_OptimalAniso!(urec_np1, solver.rec_fes,
                                                solver.nFacts, max(1, solver.AMA_AnisoOrder), 
                                                DeltaAR=10.0, q=solver.SpaceNorm)
                                    
                    #Isotropic mesh adaptation:
                    if solver.AMA_AnisoOrder==0
                        metric_np1.sf[1]    .= 1.0
                        metric_np1.sf[2]    .= 1.0
                        emeanv              .= emeanv0
                    end
                    
                    #Set AMA tolerance:
                    AMA_TolS                *= AMA_factor
                    
                    #Compute optimal sizes and impose metric eigenvalues:
                    A1v                         = AMA_OptimalSize(etaS_elems_np1, 
                                                    solver.Av, fill(solver.AMA_SizeOrder+1, solver.mesh.nElems), 
                                                    fill(1,solver.mesh.nElems), AMA_TolS, 
                                                    q=solver.SpaceNorm)
                    @mlv metric_np1.lambda_bar  = TrElem_Area/A1v
                    nElems_est                  = AMA_CorrectMetric!(solver.Av, metric_np1)
                    
                    solver.tAMA                 += time()-t_ini
                    
                    #=
                    #-----------------------------------------------------------------
                    #Optimal metric for solution at t^(n-1)-:
                    
                    t_ini                   = time()
                    
                    #Compute optimal anisotropy:
                    emeanv0, emeanv, 
                        err0, err1, 
                        metric_n            = AMA_OptimalAniso!(urec_n, solver.rec_fes,
                                                solver.nFacts, max(1, solver.AMA_AnisoOrder), 
                                                DeltaAR=10.0, q=solver.SpaceNorm)
                                    
                    #Isotropic mesh adaptation:
                    if solver.AMA_AnisoOrder==0
                        metric_n.sf[1]      .= 1.0
                        metric_n.sf[2]      .= 1.0
                        emeanv              .= emeanv0
                    end
                    
                    #Compute optimal sizes and impose metric eigenvalues:
                    A1v                         = AMA_OptimalSize(etaS_elems_n, 
                                                    solver.Av, fill(solver.AMA_SizeOrder+1, solver.mesh.nElems), 
                                                    fill(1,solver.mesh.nElems), AMA_TolS, 
                                                    q=solver.SpaceNorm)
                    @mlv metric_n.lambda_bar    = TrElem_Area/A1v
                    nElems_est                  = AMA_CorrectMetric!(solver.Av, metric_n)
                    
                    solver.tAMA                 += time()-t_ini
                    =#
                    
                    #----------------------------------------------------------
                    #Generate new mesh:
                    
                    t_ini                   = time()
                    
                    #Compute metric intersection:
#                     metric                  = AMA_MetricIntersection([metric_n, metric_np1])
                    metric                  = metric_np1
                    
                    #Generate metric tensor for Bamg:
                    MetricFile              = "$(@__DIR__)/../../temp/BamgSC$(solver.SC).mtr"
                    MeshFileDb  = "$(@__DIR__)/../../temp/BamgSC$(solver.SC).db"
                    Bamg_MetricFile!(solver.mesh, metric, MetricFile)
                    
                    #Call bamg to refine. The new mesh is first written in .db.new file and 
                    #then overwritten in the solver.MeshFile:
                    run(pipeline(`bamg -KeepBackVertices -ratio 2.7 -nbv 1000000 -hmin 1e-9 -hmax 1e9 -v 0 -M $(MetricFile) -b $(MeshFileDb) -o $(MeshFileDb).new`, 
                        stdout="$(@__DIR__)/../../temp/BamgSC$(solver.SC).out", 
                        stderr="$(@__DIR__)/../../temp/BamgSC$(solver.SC).err")) 
                    mv(string(MeshFileDb,".new"), MeshFileDb, force=true)
                    ReadMesh!(solver)
                    solver.tAMA                 += time()-t_ini
                    
                    t_ini           = time()
                    
                    #Allocate variables for time integration:
                    RKAlloc!(solver)
                    
                    solver.tAlloc   += time()-t_ini
                    
                    #------------------------------------------------------------
                    #New solution:
                    
                    #Project "initial condition" u_n. The solution is overwritten in solver.u.
                    #Note, however, that solver.u is a new pointer after reading the mesh. Hence, the 
                    #previous solution is still stored in u_n.
                    L2Projection!(solver, u_n, fes_n)
                    
                    #Correct time step:
                    Deltat_ratio    = max(0.2, log(solver.SfS*solver.TolS_max/etaS_n)/solver.dlogetaS_dt)
                    if isnan(Deltat_ratio) #happens if etaS_0=0 (initial condition = 0 and dlogetaS_dt=0)
                        Deltat_ratio    = 0.2
                    end
                    Deltat_n        *= Deltat_ratio
                    t_np1           = t_n + Deltat_n
                    
                    #Adjust time step:
                    Deltat_n        = Deltat_nm1*solver.AMA_rDeltat
                    t_np1           = min(t_n + Deltat_n, solver.tf)
                    Deltat_n        = t_np1-t_n
                    
                    #Next iteration:
                    RepeatTA        = false
                    RepeatS         = true
                    etaS_np1_est    = etaS_n * exp(solver.dlogetaS_dt * Deltat_n)
                    etaT_np1_est    = etaT_n * Deltat_n^(solver.RK.order) 
                    TolA            = max(solver.TolA_min, solver.CA*min(etaS_n, etaS_np1_est, etaT_np1_est))
            
                    #Print mesh info:
#                     printstyled("AMA_iters=", AMA_iters, ", nElems=", solver.mesh.nElems, 
#                             ", nNodes=", solver.mesh.nNodes, ", TotalDof=", solver.fes.nDof*solver.nVars, 
#                             ", TolS_min=", sprintf1("%.2e", solver.TolS_min),
#                             ", etaS=", sprintf1("%.2e", etaS_np1), 
#                             ", TolS=", sprintf1("%.2e", AMA_TolS[solver.RK.stages]), 
#                             ", TolS_max=", sprintf1("%.2e", solver.TolS_max),
#                             ", TolS/etaS=", sprintf1("%.2e", solver.AMA_TolS_etaS),
#                             "\n")
                    printstyled("AMA_iters=", AMA_iters, ", nElems=", solver.mesh.nElems, 
                            ", nNodes=", solver.mesh.nNodes, ", TotalDof=", solver.fes.nDof*solver.nVars, 
                            "\n")
                            
                #If all the errors satisfy tolerance, terminate iterations:
                else
                
                    RepeatTA        = false
                    RepeatS         = false
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
                push!(solver.validv, false)
                
                printstyled("Nonlinear solver did not converge. etaA=", sprintf1("%.2e", etaA_np1), 
                    ", TolA=", sprintf1("%.2e", TolA), 
                    ". Decreasing time step\n", color=:light_yellow)
                if LSFailed_iters==10
                    printstyled("Nonlinear solver unable to converge after ",
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
            
        end #Time-algebraic loop
        
    end #Space adaptation loop
    
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
            ", MasterDof=", sprintf1("%d", solver.Am_LS.LS.Pl.nMasters),
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
    if mod(solver.Nt,50)==0
        GC.gc()
    end
    
    solver.tCPU         += time()-t_start 
    
    return ConvFlag
    
end
