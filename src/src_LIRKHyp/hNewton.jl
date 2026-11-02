function hNewton_Step!(solver::SolverData, Ste::Float64, D_L::Float64)
    
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
    if isnan(TolS_req)  #Happens for the first time step, if etaS=0 and dlogetaS_dt=Inf:
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
    if isnan(TolA)
        TolA            = solver.TolA_min
    end
    
    #Mass and stiffness matrices:
    MII                 = solver.MII
    KII                 = StiffnessMatrix(solver.Integ2D, solver.fes)
    nDofP               = solver.fes.PSpace.nDof
    nDof                = solver.fes.nDof
    
    #Zero bubbles:
    u_n[nDofP+1:nDof]   .= 0.0
    
    #Kirchhoff variable and diffusion at tn:
    kirchh_n            = zeros(solver.fes.nDof)
    kirchh_n[1:nDofP]   .= KirchhoffFun.(u_n[1:nDofP], Ste, D_L)
    D_n                 = zeros(solver.fes.nDof)
    D_n[1:nDofP]        .= DiffFun.(u_n[1:nDofP], Ste, D_L)
    #(Assume bubbles are zero)
    
    #Get Dirichlet dofs:
    DirDof                  = zeros(Bool, nDof)
    NextVertex              = [2,3,1]
    for jb=1:solver.nBounds*3
        bmesh               = solver.mesh.bmesh[jb]
        ib                  = bmesh.boundary_id
        if solver.BC[ib].BC isa Dirichlet
            ParentElems         = bmesh.ParentElems
            ParentFace          = bmesh.ParentFace
            LocalDof            = vcat(ParentFace, NextVertex[ParentFace], 
                                    3 .+ (solver.FesOrder-1)*(ParentFace-1) .+ 
                                        (1:solver.FesOrder-1))
            GlobalDof           = solver.fes.PSpace.ElemsDof[ParentElems, LocalDof][:]
            DirDof[GlobalDof]   .= true
        end
    end
#     display(sum(DirDof))
    DirDof[nDofP+1:nDof]    .= true      #Bubble dof are also set to zero
    NDirDof                 = .!DirDof
    
    #Flux for first stage:
    t_ini               = time()
    flag, Deltat_CFL    = Rhs_Stefan!(solver, t_n, u_n, KII, kirchh_n, D_n, 
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
                    u           = zeros(nDof)
                    u[NDirDof]  .= @tturbo @. scalv[NDirDof]*uhat
                    #DirDof are set to zero
                    
#                     println("Unscale: ", time()-t_ini2)
                    
#                     figure()
#                     plot3D(solver.fes.PSpace.NodesCoords[:,1], 
#                             solver.fes.PSpace.NodesCoords[:,2],
#                             u[1:nDofP], ".b")
#                     aux     = DirDof[1:nDofP]
#                     plot3D(solver.fes.PSpace.NodesCoords[aux,1], 
#                             solver.fes.PSpace.NodesCoords[aux,2],
#                             u[1:nDofP][aux], ".r")
#                     error("")
                    
                    #Compute Kirchhoff variables:
                    kirchh              = zeros(solver.fes.nDof)
                    kirchh[1:nDofP]     .= KirchhoffFun.(u[1:nDofP], Ste, D_L)
                    #Assume bubbles are zero
    
                    #Compute diffsion coefficients:
                    Dv          = zeros(solver.fes.nDof)
                    Dv[1:nDofP] .= DiffFun.(u[1:nDofP], Ste, D_L) 
                    Dm          = spdiagm(Dv)
                    
                    #Compute flux. Save flux directly in f_RK:
                    flag,       = Rhs_Stefan!(solver, t_k, u, KII, kirchh, Dv,  
                                    view(solver.f_RK,:,kk), solver.Jm)
                    if flag<0
                        return flag
                    end
                    
#                     println("Rhs: ", time()-t_ini2) #BOTTLENECK
                    
                    #Compute system matrix and set up linear system, if necessary.
                    AII         = MII + (Deltat_n*solver.RK.AI[kk,kk])*(KII*Dm) 
                    AII_UU      = AII[NDirDof, NDirDof]
                    
#                     println("Linear system: ", time()-t_ini2)
                    
                    #Compute residual f:
                    #   f = M*u - b - Deltat*a_kk*f(t_k,u_k)
                    fres        = MuProduct(solver.MII, u, solver.nVars) - 
                                        solver.bv - (Deltat_n*solver.RK.AI[kk,kk])*
                                        view(solver.f_RK,:,kk)
                    fres_U      = fres[NDirDof] # -A_UD*u_D = 0
                                        
#                     println("fres: ", time()-t_ini2)
                    
                    #Compute preconditioned residual. Since this equation is linear and 
                    #f does not need to be computed at each iteration, we ask for a very
                    #small tolerance. In that case, having an accurate Newton step gres will 
                    #allow for a reduction in the number of Newton steps (and evaluations of f):
                    flag        = NaN
                    if true && uppercase(solver.LSType)=="SCILU0"
                        flag    = 1
                        gres    .= AII_UU\fres_U
                    else
                        LSOutput    = LS_solve!(solver.Am_LS.LS, gres, fres, 
                                    RelTol=0e-3, AbsTol=0.1*TolA, 
                                    Display="notify", MaxIter=50)
                        #Returns flag, nIters, etaA
                        flag        = LSOutput[1]
                    end
                    
                    
                    #Scale gres:
                    @tturbo @. gres     /= scalv[NDirDof]
                    
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
                                u_k[NDirDof]./scalv[NDirDof], 
                                UpdateType=2, 
                                AbsTolX=1.0*TolA, RelTolX=0.0, 
                                AbsTolG=0.0*TolA, RelTolG=0.0, 
                                NormFun=FW_NLS_norm((x)->norm(x)/sqrt(length(u_k))), 
                                memory=NLS_mem, MaxIter=solver.LS_iters_max, 
                                history=true, Display="final")
                solver.tLS  += time()-t_ini
                u_k[NDirDof].= LSOutput[1].*scalv[NDirDof]
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


#Function to compute the r.h.s. term.
#NOTE: Here "f" is the r.h.s. vector
#NOTE: Most part of the time is spent on allocating matrices (nElems, DofPerElem^2)
#TODO: Instruction f[II] is not type-stable
function Rhs_Stefan!(solver::SolverData, t::Float64, uv::Vector{Float64}, 
    KII::SparseMatrixCSC{Float64,Int}, kirchhv::Vector{Float64}, Dv::Vector{Float64}, 
    fv::GenVector{Float64}, Jm::SparseMatrixCSC{Float64,Int})
    
    t_ini       = time()

    ComputeJ    = false
    
    #Extract variables:
    fes         = solver.fes
    fes_dof     = fes.nDof
    nVars       = solver.nVars
    
    u           = GetViews(uv, nVars, fes_dof)
    f           = GetViews(fv, nVars, fes_dof)
    @mlv fv     = 0.0
    
    #Allocate Jacobian:
    J_iv        = zeros(Int,0)
    J_jv        = zeros(Int,0)
    J_sv        = zeros(0)
    
    #Initialize monitor:
    solver.monitor  = Vector{Matrix{Float64}}(undef,0)
    
    #flag=1 if there are no NaNs nor Infs in the fluxes and source terms,
    #flag=-1 otherwise
    flag        = 1

    #----------------------------------------------------------------
    #Variables for domain integral:
    
    Integ2D         = solver.Integ2D
    Jinvm           = Integ2D.Jinv
    wdetJ           = Integ2D.wdetJ
    
    #Compute variables at quadrature nodes:
    _qp             = TrIntVars()
    _qp.Integ2D     = Integ2D
    _qp.nElems      = Integ2D.mesh.nElems
    _qp.nqp         = Integ2D.QRule.nqp
    _qp.FesOrder    = solver.FesOrder
    _qp.nFacts      = solver.nFacts
    #
    _qp.t           = t
    _qp.x           = QuadNodesCoords(Integ2D)
    _qp.u           = SolutionCompute(u, fes, Integ2D.QRule.xi)
    uB              = GetBubbles(u, fes)
    _qp.uB          = SolutionCompute(uB, fes.BSpace, Integ2D.QRule.xi)
    _qp.gradu       = Solution_gradCompute(u, fes, Integ2D.QRule.xi, Integ2D.Jinv)
    _qp.graduB      = Solution_gradCompute(uB, fes.BSpace, Integ2D.QRule.xi, Integ2D.Jinv)
    #
    _qp.f, _qp.df_du, _qp.df_dgradu         = FluxAllocate(nVars, (_qp.nElems, _qp.nqp), ComputeJ)
    _qp.fB, _qp.dfB_du, _qp.dfB_dgraduB     = FluxAllocate(nVars, (_qp.nElems, _qp.nqp), ComputeJ)
    _qp.Q, _qp.dQ_du, _qp.dQ_dgradu         = bFluxAllocate(nVars, (_qp.nElems, _qp.nqp), ComputeJ)
    
    #Shape functions:
    Nm              = NCompute(fes, Integ2D.QRule.xi)
    gradNm          = gradNCompute(fes, Integ2D.QRule.xi)   #Vector{Matrix{Float64}} with Grad1 and Grad2
    flucNm          = flucNCompute(fes, Integ2D.QRule.xi)   #fluc = phi-interp(phi)
    gradflucNm      = gradflucNCompute(fes, Integ2D.QRule.xi)
    ElemsDof        = fes.ElemsDof
    
    #Allocate memory for a(u,psi_\alpha):
    flux_ElemsDof               = Vector{Matrix{Float64}}(undef,nVars)
    for II=1:nVars
        flux_ElemsDof[II]       = zeros(solver.mesh.nElems, fes.DofPerElem)
    end
    J_ElemsDof                  = Matrix{Matrix{Float64}}(undef,nVars,nVars)
    if ComputeJ
        for JJ=1:nVars, II=1:nVars
            J_ElemsDof[II,JJ]   = zeros(solver.mesh.nElems, solver.Jm_localterms)
        end
    end
    flux_Ik                     = zeros(solver.mesh.nElems, Integ2D.QRule.nqp)
    
#     println("Memory allocation = ", time()-t_ini)
    
    #----------------------------------------------------------------
    #Contribution of flux and source terms in the domain:
    
    #Evaluate flux and source term:
    t_ini           = time()
#     FluxSource!(solver.model, _qp, ComputeJ)
#     println("Flux and source terms = ", time()-t_ini)
    
    #Check NaN's:
#     for II=1:solver.nVars
#         if !all(isfinite.(_qp.f[II]))
#             printstyled("flux $(II) is not finite", color=:orange)
#             flag    = -1
#         end
#         if !all(isfinite.(_qp.fB[II]))
#             printstyled("bubble flux $(II) is not finite", color=:orange)
#             flag    = -1
#         end
#         if !all(isfinite.(_qp.Q[II]))
#             printstyled("source term $(II) is not finite", color=:orange)
#             flag    = -1
#         end
#     end
    
    t_ini           = time()
    
    #Fluxes. We assume f=f(x, t, u, grad u) and fB=fB(x, t, u, grad uB),
    #with u the total (large+fine scale) solution and uB the bubble terms.
#     Rhs_Flux!(flux_ElemsDof, J_ElemsDof, flux_Ik, Integ2D, 
#         _qp.f, _qp.df_du, _qp.df_dgradu,
#         gradNm, Nm, gradNm, ComputeJ, solver.Jm_ilocal, solver.Jm_jlocal)
#     Rhs_Flux!(flux_ElemsDof, J_ElemsDof, flux_Ik, Integ2D, 
#         _qp.fB, _qp.dfB_du, _qp.dfB_dgraduB, 
#         gradflucNm, Nm, gradflucNm, ComputeJ, solver.Jm_ilocal, solver.Jm_jlocal)
        
    #Source:
    flux_I      = flux_Ik   #Change nomenclature for pointer
    Rhs_Source!(flux_ElemsDof, J_ElemsDof, flux_I, Integ2D, 
                _qp.Q, _qp.dQ_du, _qp.dQ_dgradu, Nm, Nm, gradNm, 
                ComputeJ, solver.Jm_ilocal, solver.Jm_jlocal)
    
#     println("Flux and source terms * shape functions = ", time()-t_ini)
    
#     for II=1:4, ii=1:2
#         display(norm(_qp.f[II,ii]))
#     end
#     for II=1:4, ii=1:2
#         display(norm(_qp.fB[II,ii]))
#     end
#     for II=1:4
#         display(norm(_qp.Q[II]))
#     end
#     error("")
    
    #---------------------------------------------------------------------
    #Boundary terms:
    
    t_ini           = time()
    tB1             = 0.0
    tB2             = 0.0
    tB3             = 0.0
    tB4             = 0.0
    tB5             = 0.0
    tB6             = 0.0
    
    #Loop boundaries:
    for jb=1:solver.nBounds*3
        
        t_ini2              = time()
        Binteg2D            = solver.Binteg2D[jb]
        bmesh               = Binteg2D.bmesh
        Jinvm               = Binteg2D.Jinv
        wdetJ               = Binteg2D.wdetJ
        bound_id            = bmesh.boundary_id
        
        #Compute variables at quadrature nodes:
        _bqp                = TrBintVars()
        _bqp.Binteg2D       = Binteg2D
        _bqp.nElems         = bmesh.nElems
        _bqp.nqp            = Binteg2D.QRule.nqp
        _bqp.FesOrder       = solver.FesOrder
        _bqp.nFacts         = solver.nFacts
        _bqp.tb             = Binteg2D.t_qp
        _bqp.nb             = Binteg2D.n_qp
        #
        _bqp.t              = t
        _bqp.x              = QuadNodesCoords(Binteg2D)
        _bqp.u              = SolutionCompute(Binteg2D, u, fes)
        _bqp.gradu          = Solution_gradCompute(Binteg2D, u, fes)
        #
        _bqp.f, _bqp.df_du, 
            _bqp.df_dgradu  = bFluxAllocate(nVars, (_bqp.nElems, _bqp.nqp), ComputeJ)
        
        #Shape functions at quadrature nodes:
        Nm                  = NCompute(fes, Binteg2D.xi2D)      #Vector{MFloat}
        gradNm              = gradNCompute(fes, Binteg2D.xi2D)  #[Vector{MFloat}, Vector{MFloat}]
        tB1                 += time()-t_ini2
        
        t_ini2              = time()
        #Allocate memory for f psi_alpha:
        flux_BElemsDof              = Vector{Matrix{Float64}}(undef,nVars)
        for II=1:nVars
            flux_BElemsDof[II]      = zeros(bmesh.nElems, fes.DofPerElem)
        end
        J_BElemsDof                 = Matrix{Matrix{Float64}}(undef,nVars,nVars)
        if ComputeJ
            for JJ=1:nVars, II=1:nVars
                J_BElemsDof[II,JJ]  = zeros(bmesh.nElems, solver.Jm_localterms)
            end
        end
        bflux_I                     = zeros(bmesh.nElems, Binteg2D.QRule.nqp)
        tB2                         += time()-t_ini2
        
        t_ini2                      = time()
        #Evaluate flux at the boundary:
        bflux!(solver.model, solver.BC[bound_id].BC, _bqp, ComputeJ)
        tB3                         += time()-t_ini2
        
        #Check NaN's:
        t_ini2                      = time()
        for II=1:solver.nVars
            if !all(isfinite.(_bqp.f[II]))
                printstyled("flux $(II) for boundary $(bound_id) is not finite \n", color=:orange)
                flag                = -1
            end
        end
        tB4                         += time()-t_ini2
        
        t_ini2                      = time()
        #Compute contribution of boundary term:
        Rhs_bFlux!(flux_BElemsDof, J_BElemsDof, bflux_I, Binteg2D,
            _bqp.f, _bqp.df_du, _bqp.df_dgradu, 
            Nm, Nm, gradNm, ComputeJ, solver.Jm_ilocal, solver.Jm_jlocal)
        tB5                         += time()-t_ini2
        
        #Update global values for flux_ElemsDof and J_ElemsDof:
        t_ini2                      = time()
        @inbounds for II=1:nVars, iDof=1:fes.DofPerElem, iElem=1:bmesh.nElems
            parent_elem                             = bmesh.ParentElems[iElem]
            flux_ElemsDof[II][parent_elem, iDof]    += flux_BElemsDof[II][iElem, iDof]
        end
        if ComputeJ
            @inbounds for JJ=1:nVars, II=1:nVars, iJ=1:solver.Jm_localterms, iElem=1:bmesh.nElems
                parent_elem                         = bmesh.ParentElems[iElem]
                J_ElemsDof[II,JJ][parent_elem,iJ]   += J_BElemsDof[II,JJ][iElem, iJ]
            end
        end
        tB6                         += time()-t_ini2
        
#         display("")
#         for II=1:4
#             display(norm(_bqp.f[II]))
#         end
#         error("")
        
    end
    
#     println("Boundary terms 1 = ", tB1)
#     println("Boundary terms 2 = ", tB2)
#     println("Boundary terms 3 = ", tB3)
#     println("Boundary terms 4 = ", tB4)
#     println("Boundary terms 5 = ", tB5)
#     println("Boundary terms 6 = ", tB6)
#     println("Boundary terms = ", time()-t_ini)

#     display(flux_ElemsDof[4])
#     error("")
    
    #---------------------------------------------------------------------
    #Assemble:
    
    t_ini               = time()
    for II=1:nVars
        VectorAssemble!(ElemsDof, flux_ElemsDof[II], f[II])
    end
#     println("Assembly vector = ", time()-t_ini)
    
    #Add -K*kirchh:
    f[1]                .= f[1] - KII*kirchhv
    
    t_ini               = time()
    if ComputeJ
        
        #Zero jacobian:
        BLAS.scal!(0.0, Jm.nzval)
        
        #Update blocks:
        for JJ=1:nVars, II=1:nVars
            sparse!(Jm, solver.Jm_pinv[II,JJ], J_ElemsDof[II,JJ])
        end
        
    end
#     println("Assembly Jacobian = ", time()-t_ini)
    
#     error("")
    
    #CFL number:
    hp_elems        = _hmin(_qp.Integ2D.mesh)./_qp.FesOrder 
    D_ElemsNodes    = Dv[solver.fes.PSpace.ElemsDof]
    Dmax_elems      = maximum(D_ElemsNodes, dims=2)
    _qp.Deltat_CFL  = @. $minimum(hp_elems^2 / Dmax_elems)
    
    return flag, _qp.Deltat_CFL
    
end

function KirchhoffFun(h::Float64, Ste::Float64, D_L::Float64)
    
    if h<=1.0
        return h
    elseif h<1.0+1.0/Ste
        return 1.0
    else
        return 1.0 + D_L*(h-1.0-1.0/Ste)
    end
    
end

function DiffFun(h::Float64, Ste::Float64, D_L::Float64)
    
    if h<=1.0
        return 1.0
    elseif h<1.0+1.0/Ste
        return 0.0
    else
        return D_L
    end
    
end

#DEPRECATED:

#= 

any(uppercase(solver.LSType).==["BLOCKJACOBIKRYLOV", "SCILU0KRYLOV"])
            
                #Approximation of J*g:
                Jg              = zeros(length(u_k))
                epsilon_Krylov  = 1e-5
                f_evals         = 0
            
                #Define preconditioned residual:
                #   fres = M*u - b - a_ii*Deltat_n * f(u)
                #   gres = A^{-1} fres
                function QNResidualKrylov!(uhat::Vector{Float64}, gres::Vector{Float64})
                    
                    #Unscale u:
                    u           = @tturbo @. scalv*uhat
                    
                    #Compute flux and save derivatives:
    #                 t_ini                   = time()
                    flag,               = Rhs!(solver, t_k, u, false, 
                                            view(solver.f_RK,:,kk), solver.Jm)
    #                 solver.tRhs             += time()-t_ini
                    f_evals             += 1
                    if flag<0
                        return flag
                    end
                    
                    #Compute residual f:
                    #   fres = M*u - b - Deltat*a_kk*f(t_k,u_k)
                    fres                = MuProduct(solver.MII, u, solver.nVars) - 
                                            solver.bv - (Deltat_n*solver.RK.AI[kk,kk])*
                                            view(solver.f_RK,:,kk)
                    
                    #Instead of solving Ahat*gres=fres, where Ahat is the preconditioner, 
                    #we solve A(t_k,u_k)*gres=fres, with A the exact Jacobian of fres, 
                    #in a Krylov fashion. That is, we solve for gres in
                    #   fres2   = M*gres - Deltat*a_kk*(f(u_k+epsilon*gres)-f(u_k))/epsilon - fres = 0
                    #This equation is preconditioned with Ahat. That is, we solve
                    #   gres2(gres) = Ahat\fres2(gres) = 0
                    function QNResidualG!(ghat::Vector{Float64}, gres2::Vector{Float64})
                    
                        #Unscale g:
                        g           = @tturbo @. scalv*ghat
                    
                        #Compute (f(t_k, u_k + epsilon*g)-f(t_k, u_k))/epsilon:
                        flag,       = Rhs!(solver, t_k, u+epsilon_Krylov*g, false, 
                                        Jg, solver.Jm)
                        if flag<=0
                            return flag
                        end
                        BLAS.axpby!(-1.0, view(solver.f_RK,:,kk), 1.0, Jg)
                        BLAS.scal!(Deltat_n*solver.RK.AI[kk,kk]/epsilon_Krylov, Jg)
                    
                        #Residual:
                        fres2       = MuProduct(solver.MII, g, solver.nVars) - Jg - fres
                        
                        #Compute preconditioned residual gres2 = Ahat\fres2. This system is linear
                        #and does not require to evaluate f. We solve the system with 
                        #high accuracy in order to reduce the number of evaluations of f in 
                        #outer loops.
                        LSOutput    = LS_solve!(solver.Am_LS.LS, gres2, fres2, 
                                        RelTol=0.0, AbsTol=TolA, 
                                        Display="notify", MaxIter=solver.LS_iters_max)
                        flag        = LSOutput[1]
                        
                        #Scale residual:
                        @tturbo @. gres2    /= scalv
                        
                        return flag 
                        
                    end
                    #Solve for gres. This value provides the Newton step. It is not necessary
                    #to ask for a small tolerance since this is only the Newton step, and anyway
                    #f has to be computed at each iteration:
                    t_ini       = time()
                    LSOutput    = Anderson(FW_NLS((ghat,gres2)->QNResidualG!(ghat,gres2)), 
                                    gres./scalv, 
                                    AbsTolX=0.0*TolA, RelTolX=0.0, 
                                    AbsTolG=1.0*TolA, RelTolG=1e-2, 
                                    NormFun=FW_NLS_norm((x)->norm(x)/sqrt(length(u_k))), 
                                    memory=100, MaxIter=solver.LS_iters_max, 
                                    history=true, Display="notify")
                    solver.tLS  += time()-t_ini
                    gres        .= LSOutput[1].*scalv
                    flag        = LSOutput[2].flag
                    f_evals     += LSOutput[2].nIter
                    
                    #Scale gres:
                    @tturbo @. gres     /= scalv
                    
                    return flag
                    
                end
            
                #Solve:
                t_ini       = time()
                LSOutput    = Anderson(FW_NLS((uhat,gres)->QNResidualKrylov!(uhat,gres)), 
                                u_k./scalv, 
                                AbsTolX=1.0*TolA, RelTolX=0.0, 
                                AbsTolG=0.0*TolA, RelTolG=0.0, 
                                NormFun=FW_NLS_norm((x)->norm(x)/sqrt(length(u_k))), 
                                memory=0, MaxIter=solver.LS_iters_max, 
                                history=true, Display="final")
                solver.tLS  += time()-t_ini
                u_k         .= LSOutput[1].*scalv
                LSFlag      = LSOutput[2].flag
                etaA_np1    = LSOutput[2].pnorms[ length(LSOutput[2].pnorms) ]
                solver.LS_iters     += f_evals
                solver.LS_total     += f_evals
                if LSFlag<0
                    #NLS did not converge because something went wrong with the residual function
                    etaA_np1        = Inf   #In this way, time step is decreased as maximum as possible
                    break #loop for kk=2:nStages
                elseif LSFlag==0
                    #NLS reached maximum of iterations
                    break #loop for kk=2:nStages
                else
                    printstyled("NLSolver converged, f_evals=", sprintf1("%d", f_evals), " \n", color=:cyan)
                end
                semilogy(LSOutput[2].gnorms)
                #Here, LS in reality is NLS
                
            elseif
            
=#
