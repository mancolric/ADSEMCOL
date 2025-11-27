#Create structure with default options:
function LIRKHyp_Start(model::ConstModel) where {ConstModel<:ConstModels}

    #Create empty structure:
    solver              = SolverData{ConstModel}()
    solver.model        = model
    
    #Problem name and StudyCase
    solver.ProblemName  = ""
    solver.SC           = 0
    
    #Mesh:
    #Unassigned
    
    #Initial condition and boundary conditions:
    #Unassigned
    
    #Solution:
    solver.nVars        = model.nVars
    solver.FesOrder     = 5
    #fes
    #uv
    #u
    solver.SpaceAdapt   = true
    #nFacts
    solver.SpaceNorm    = 2
    
    #Reconstruction:
    #Unassigned
    
    #Time integration method:
    #t
    #Nt
    solver.RKMethod     = "Ascher3"
    solver.Deltat0      = 1e-6
    #Deltat
    solver.CFL          = NaN
    solver.tf           = -1.0
    solver.TimeAdapt    = true
    #RK
    #f_RK
    #Ju_RK
    #fNStab_RK
    
    #Matrices/vectors:
    #Unassigned
    
    #Monitor:
    solver.monitor          = Vector{Matrix{Float64}}(undef,0)
    
    #Tolerances:
    solver.TolS_min         = 1e-3
    solver.TolS_max         = 1e-2
    solver.TolT             = 1e-4
    solver.TolA_min         = 1e-12
    solver.CA               = 0.01
    solver.CA_max           = 0.05
    solver.AMA_MaxIter      = 20
    solver.AMA_SizeOrder    = 2
    solver.AMA_AnisoOrder   = 2
    solver.AMA_ProjN        = 1
    solver.AMA_ProjOrder    = 0         #This is a flag to set a default parameter later
    solver.SfS              = 0.8
    solver.SfT              = 0.8
    solver.AMA_rDeltat      = 0.5
    solver.LS_iters_target  = 60
    solver.LS_iters_max     = 100
    
    #Errors:
    #Unassigned
    solver.AMA_TolS_etaS    = 1.0
    
    #Vector with time levels, errors, number of elements, etc.:
    #Unassigned
    
    #Total CPU time:
    solver.tCPU             = 0.0
    solver.tAMA             = 0.0
    solver.tJm              = 0.0
    solver.tSCILU           = 0.0
    solver.tLS              = 0.0
    solver.tRhs             = 0.0
    solver.tAlloc           = 0.0
    solver.tb               = 0.0
    solver.tetaS            = 0.0
    solver.tetaT            = 0.0
    solver.LS_total         = 0
    
    return solver
    
end

#Read (linear) mesh and allocate necessary variables:
function ReadMesh!(solver::SolverData)

    #Read linear mesh:
    solver.mesh         = TrMesh_Read("$(@__DIR__)/../../temp/BamgSC$(solver.SC).db", 1, solver.nBounds)
    
    #Create finite element space, matrices, integrals, auxiliary vectors:
    solver.fes          = TrPBSpace(solver.mesh, solver.FesOrder)
    solver.Integ2D      = TrInt(solver.mesh, 2*(solver.FesOrder+2)+1)
    solver.Binteg2D     = Vector{TrBint}(undef, solver.nBounds*3)
    for ib=1:solver.nBounds*3
        solver.Binteg2D[ib]     = TrBint(solver.mesh, ib, 2*(solver.FesOrder+2))
    end
    solver.uv           = zeros(solver.nVars*solver.fes.nDof)
    solver.u            = GetViews(solver.uv, solver.nVars, solver.fes.nDof)
    solver.MII          = MassMatrix(solver.Integ2D, solver.fes)
    solver.bv           = zeros(solver.nVars*solver.fes.nDof)
    solver.b            = GetViews(solver.bv, solver.nVars, solver.fes.nDof)
    solver.Av           = sum(solver.Integ2D.wdetJ, dims=2)[:]
    solver.Omega        = sum(solver.Av)
    
    #Variables for space reconstruction:
    solver.rec_fes      = TrQSpace(solver.mesh, solver.FesOrder+2)
    
    #Condense and precondition mass matrix. Note that 
    solver.MII_LS           = LinearSystem1(solver.MII, solver, [1.0])
    LinearSystem!(solver.MII_LS)
    
    return
    
end

#Allocate variables for RK method:
function RKAlloc!(solver::SolverData)

    t_ini               = time()
    
    #f and Ju terms:
    solver.f_RK         = zeros(solver.nVars*solver.fes.nDof, solver.RK.stages)
    solver.Ju_RK        = zeros(solver.nVars*solver.fes.nDof, solver.RK.stages)
    
    #Allocate expanded mass matrix, Jacobian matrix, linear system matrix
    #and indices for fast assembly in Jacobian matrix:
    MatricesAllocate!(solver)
#     println("Matrices allocate = ", time()-t_ini)
    
    #Allocate linear system:
    solver.Am_LS        = LinearSystem1(solver.Am, solver, solver.nFacts)
#     println("allocate LS = ", time()-t_ini)
    
    return
    
end

#Possible flags:
#   1: method converged
#   -1: linear solver for (M-a_kk*Deltat*J)^failed
#   -2: linear solver for embedded RK failed
#   -(2+i): AMA failed with flag=`-i 
function LIRKHyp_InitialCondition!(solver::SolverData; 
    AMA_RefineFactor::Float64=0.5, DEq_MaxIter::Int=10)

    t_start     = time()
    
    #Check inputs:
    if solver.tf<0.0
        error("Final time tf not set")
    end
    if solver.AMA_SizeOrder>solver.FesOrder
        @warn "Changing AMA_SizeOrder to $(solver.FesOrder)"
        solver.AMA_SizeOrder    = solver.FesOrder
    end
    if solver.AMA_AnisoOrder>solver.FesOrder
        @warn "Changing AMA_AnisoOrder to $(solver.FesOrder)"
        solver.AMA_AnisoOrder   = solver.FesOrder
    end
    if solver.AMA_rDeltat>1.0
        @warn "AMA_rDeltat must be smaller than 1.0. Setting AMA_rDeltat=1.0"
    end
    if solver.Deltat0<=0.0
        error("Initial Deltat must be positive")
    end
    
    #Set some default parameters:
    if solver.AMA_SizeOrder==0
        solver.AMA_SizeOrder    = solver.FesOrder
    end
    if solver.AMA_ProjOrder==0
        solver.AMA_ProjOrder    = 2*(solver.FesOrder+2+1)
    end
    
    #Copy db file or create from geo. The input is solver.MeshFile.
    #The output (copied db file or new db file from geo) is solver.MeshFileSC.
    MeshFileGeo         = "$(@__DIR__)/../../temp/BamgSC$(solver.SC).geo"
    MeshFileDb          = "$(@__DIR__)/../../temp/BamgSC$(solver.SC).db"
    MeshFileExtension   = solver.MeshFile[findlast('.',solver.MeshFile)+1:length(solver.MeshFile)]
    if MeshFileExtension=="db"
        #We need the geo file, because we need to copy it:
        error("Mesh file has to be of geo type")
        #Copy db file to DbFile:
        cp(solver.MeshFile, MeshFileSC, force=true)
    elseif MeshFileExtension=="geo"
        #Generate db file:
        cp(solver.MeshFile, MeshFileGeo, force=true)
        run(pipeline(`bamg -AbsError -NoRescaling -NbJacobi 3 -NbSmooth 5 -hmax 1e9 -hmin 1e-9 -ratio 0 -nbv 1000000 -v 0 -g $(MeshFileGeo) -o $(MeshFileDb)`, 
            stdout="$(@__DIR__)/../../temp/BamgSC$(solver.SC).out", stderr="$(@__DIR__)/../../temp/BamgSC$(solver.SC).err"))
    else
        error("Unknown mesh file extension $(MeshFileExtension)")
    end
    #For safety, change pointer to mesh file to avoid overwritting:
    solver.MeshFile     = ""
    
    #Read mesh and allocate memory:
    ReadMesh!(solver)
    
    #Allocate variables:
    solver.tv           = zeros(0)
    solver.etaSv        = zeros(0)
    solver.etaTv        = zeros(0)
    solver.etaAv        = zeros(0)
    solver.nElemsv      = zeros(Int,0)
    solver.nDofv        = zeros(Int,0)
    solver.CFLv         = zeros(0)
    solver.validv       = zeros(Int,0)
    
    #--------------------------------------------------------
    #Initial condition:
    
    ConvFlag            = 1
    
    #Set time:
    solver.t            = 0.0
    solver.Nt           = 0
    
    #Interpolate initial condition (necessary for iterative solvers):
    Interpolate!(solver.u0fun, solver.u, solver.fes)
    for II=1:solver.nVars
        if !all(isfinite.(solver.u[II]))
            error("Initial condition for variable $(II) is not finite")
        end
    end
    
    #Normalization factors (this remains unchanged all the simulation):
    solver.nFacts       = zeros(solver.nVars)
    nFactsCompute!(solver)
    if !all(isfinite.(solver.nFacts))
        error("Normalization factors are not finite")
    end
    if any(solver.nFacts.<=0.0)
        error("Normalization factors must be positive")
    end
    
    #Project initial condition. Note: If eta_i is the Lq-error in variable "i"
    #(normalized with the domain size), we need that 
    #   eta^2 := 1/nVars * sum_i (eta_i/n_i)^2 <= TolA^2
    #A sufficient condition for this is 
    #   eta_i/n_i <= TolA
    #So, the real tolerance for variable "i" is 
    #   n_i*TolA
    L2Projection!(solver.Integ2D, solver.u0fun, solver.fes, solver.b)
    etaA                    = 0.0
    for ii=1:solver.nVars
        if !all(isfinite.(solver.bv[ii]))
            error("Initial condition for variable $(ii) is not finite")
        end
#         solver.u[ii]    = cholesky(Symmetric(solver.MII))\solver.b[ii]
        flag, nIter, etaA_I = LS_gmres!(solver.MII_LS, solver.u[ii], solver.b[ii], 
                                AbsTol=solver.TolA_min*solver.nFacts[ii], Display="notify")
        etaA                += (etaA_I/solver.nFacts[ii])^2
        if flag<=0
            printstyled("Unable to compute initial condition. Aborting\n", 
                color=:magenta)
            return -2
        end
    end
    etaA                    = sqrt(etaA/solver.nVars)
    
    #Estimate errors:
    solver.etaS, 
        solver.etaS_elems, 
        solver.urec         = etaS_reconstruction(solver, solver.u, 
                                q=solver.SpaceNorm)
#     solver.etaT             = Inf   #we set Inf so that etaT_est=Inf and TolA=CA*TolS
    solver.etaT             = max(1e-9, 0.1*solver.TolT)
    solver.Deltat           = solver.Deltat0
    solver.etaA             = etaA
    solver.dlogetaS_dt      = 0.0
    printstyled("Initial mesh: nElems=", solver.mesh.nElems, 
        ", nNodes=", solver.mesh.nNodes, ", TotalDof=", solver.fes.nDof*solver.nVars, 
        ", etaS=", sprintf1("%.2e", solver.etaS), 
        "\n")
        
    #Save errors:
    push!(solver.tv, solver.t)
    push!(solver.etaSv, solver.etaS)
    push!(solver.etaTv, solver.etaT)
    push!(solver.etaAv, solver.etaA)
    push!(solver.nElemsv, solver.mesh.nElems)
    push!(solver.nDofv, solver.nVars*solver.fes.nDof)
    push!(solver.CFLv, NaN)
    if solver.TolS_min <= solver.etaS <= solver.SfS*solver.TolS_max
        push!(solver.validv, 1)
    else
        push!(solver.validv, 0)
    end
    
    #Adapt mesh if necessary:
    AMA_flag, AMA_iters = AdaptMesh!(solver, solver.u0fun, solver.TolS_min, solver.SfS*solver.TolS_max,
                            AMA_RefineFactor=AMA_RefineFactor, DEq_MaxIter=DEq_MaxIter)
    if AMA_flag<0
        printstyled("AMA algorithm failed. Aborting\n", 
                        color=:magenta)
        ConvFlag        = -2+AMA_flag
    end
    
    #Initialize RK variables:
    solver.RK           = RK_Coefficients(solver.RKMethod)
    RKAlloc!(solver)
    
    println("t=", sprintf1("%.2e", 0.0),
            ", etaS=", sprintf1("%.2e", solver.etaS))
    
    #Save total error:
    solver.etaST        = solver.etaS
    
    #Do not count initial condition for total time
#     solver.tCPU         += time()-t_start
    
    GC.gc()
    
    return ConvFlag
    
end

function LIRKHyp_Step_Pre!(solver::SolverData)

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
    TolA            = max(solver.TolA_min, solver.CA*min(etaS_np1_est, etaT_np1_est, etaS_n))
    
    #Flux and Jacobian for first stage:
    t_ini                   = time()
    println("Computing Jacobian")
    Deltat_CFL              = Rhs!(solver, t_n, u_n, true, 
                                view(solver.f_RK,:,1), solver.Jm)
    solver.tJm              += time()-t_ini
    printstyled("Jacobian computed in ", time()-t_ini, " seconds, nnz=", 
        solver.Jm.colptr[solver.Jm.n+1], " \n", color=:cyan)
    view(solver.Ju_RK,:,1)  .= solver.Jm*u_n
    
    #---------------------------------------------------------
    #Stages 2,...,s:
    
    #Allocate variables:
    u_k                 = zeros(size(u_n))
    b_1                 = solver.Mm*u_n
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
    
        #First stage:
        @mlv u_k        = u_n
        
        #Loop stages (solution for stage kk is overwritten in solver.t and solver.u):
        LSFlag          = 0
        solver.LS_iters = 0
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
            t_ini       = time()
            LSOutput    = LS_gmres!(solver.Am_LS, u_k, solver.bv, AbsTol=TolA, Display="final",
                            MaxIter=solver.LS_iters_max)
            solver.tLS  += time()-t_ini
            LSFlag      = LSOutput[1]
            LSIter      = LSOutput[2]
            etaA_np1    = LSOutput[3]
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
            view(solver.Ju_RK,:,kk) .= solver.Jm*u_k
            
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

#We have a solution contained in solver. Adapt mesh so as to achieve desired space tolerance.
#u0 is the function we have to project (it is a given function or an old solution).
#
#There are two strategies:
#   0: If etaS<TolS_min, the mesh is coarsed until etaS>=TolS_min, then we switch to strategy 1.
#With strategy 0, TolS_max is set to sqrt(etaS, TolS_max) to prevent excessive coarsening.
#   1: If etaS>TolS_max, the mesh is refined until etaS<=TolS_max. TolS_min is disregarded.
#
#Returns flag and number of mesh adaptation iterations. Possible flags:
#   -1: Reached maximum of iterations
#   -2: Failed to compute L2 projection of initial condition because of search--locate
#   -3: Failed to compute L2 projection of initial condition because of linear system
#   -4: Failed to interpolate initial condition because of search--locate
function AdaptMesh!(solver::SolverData, 
    u0::FunType, TolS_min::Float64, TolS_max::Float64;
    AMA_CoarseFactor::Float64=5.0, AMA_RefineFactor::Float64=0.9, DEq_MaxIter::Int=0) where FunType<:Union{FW11, 
    Tuple{Vector{<:GenVector{Float64}}, TrFES}}
    
    t_ini           = time()
    
    #Leave if no AMA has to be done:
    if solver.SpaceAdapt==0 || solver.AMA_MaxIter==0 || TolS_min<=solver.etaS<=TolS_max
        return 1, 0
    end
    
    #AMA loop variables:
    AMA_iters       = 0
    AMA_TolS        = solver.etaS#*solver.AMA_TolS_etaS
    AMA_Strategy    = 0         #0:Coarse mesh, 1: Refine mesh
    AMA_factor      = AMA_CoarseFactor
    MeshFileDb      = "$(@__DIR__)/../../temp/BamgSC$(solver.SC).db"
    
    #Print info of initial mesh:
    printstyled(
                "AMA_iters=", AMA_iters, 
                ", nElems=", solver.mesh.nElems, 
                ", TotalDof=", solver.fes.nDof*solver.nVars, 
                ", TolS_min=", sprintf1("%.2e", TolS_min),
                ", etaS=", sprintf1("%.2e", solver.etaS), 
                ", TolS=", sprintf1("%.2e", NaN), 
                ", TolS_max=", sprintf1("%.2e", TolS_max),
                "\n")
    
    #Loop for AMA_TolS:
    xv      = []
    yv      = []
#     figv    = Vector{Figure}(undef,0)
#     push!(figv, figure())
#     PlotMesh!(solver.mesh)
    while true
    
        xv      = vcat(xv,solver.mesh.nElems)
        yv      = vcat(yv,solver.etaS)
    
        #Change strategy:
        if solver.etaS>TolS_max
            AMA_Strategy    = 1 #refine mesh from now on
            AMA_factor      = AMA_RefineFactor
        end
        
        #Set AMA tolerance:
        AMA_TolS            *= AMA_factor
        
        #Loop for equidistribution:
        DEq_iters           = 0
        while true
            
            t_ini2              = time()
            
            #Update number of iterations:
            AMA_iters           += 1
            DEq_iters           += 1
            
            #----------------------------------------------------
            #DEFINE OPTIMAL MESH:
            
            #Compute optimal anisotropy:
            emeanv0, emeanv, 
                err0, err1, 
                metric          = AMA_OptimalAniso!(solver.urec, solver.rec_fes,
                                    solver.nFacts, max(1, solver.AMA_AnisoOrder), 
                                    DeltaAR=10.0, q=solver.SpaceNorm)
                
            #Isotropic mesh adaptation:
            if solver.AMA_AnisoOrder==0
                metric.sf[1]    .= 1.0
                metric.sf[2]    .= 1.0
                emeanv          .= emeanv0
            end
            
#             println("OptimalAniso=", time()-t_ini2)
#             t_ini2                  = time()
            
            #Compute optimal sizes. Correct error estimation for current anisotropy with
            #error estimation for optimal anisotropy:
            A1v                     = AMA_OptimalSize(solver.etaS_elems, 
                                        solver.Av, fill(solver.AMA_SizeOrder+1, solver.mesh.nElems), 
                                        fill(1,solver.mesh.nElems), AMA_TolS, 
                                        q=solver.SpaceNorm)
            @mlv metric.lambda_bar  = TrElem_Area/A1v
            
            #Correct metric (impose maximum radius):
            nElems_est              = AMA_CorrectMetric!(solver.Av, metric)
            
            if false && solver.t+solver.Deltat>5e-4
            
                DEq_MaxIter     = 5
                
                @warn "Plotting meshes"
                fig     = PyPlotFigure(w=8.5, h=8.5, top=0.0, bottom=0.0, left=0.0, right=0.0)
                Rplot   = 0.025
                
                #Plot mesh
                PlotMesh!(solver.mesh, color="k")
                xlim(0.0, Rplot)
                ylim(0.0, Rplot)
                savefig(string(SaveDir(), "Mesh_$(AMA_iters).png"), pad_inches=0, dpi=800)
                
                #Plot ellipses:
                AMA_PlotEllipses!(solver.mesh, metric, corr_factor=1.0)
                xlim(0.0, Rplot)
                ylim(0.0, Rplot)
                savefig(string(SaveDir(), "Ellipses_$(AMA_iters).png"), pad_inches=0, dpi=800)
                
                close(fig)
                
                #Plot initial condition:
                if AMA_iters==1
                    PlotContour(solver, solver.model, "rho", delta=1e-5)
                    PlotMesh!(solver.mesh, color="k")
                    xlim(0.0, 0.025)
                    ylim(0.0, 0.025)
                    savefig(string(SaveDir(), "Remap_$(0).png"), pad_inches=0, dpi=800)
                end
                
            end
            
#             println("OptimalSize=", time()-t_ini2)
#             t_ini2                  = time()
            
            #----------------------------------------------------
            #GENERATE NEW MESH:
            
            #Generate metric tensor for Bamg:
            MetricFile              = "$(@__DIR__)/../../temp/BamgSC$(solver.SC).mtr"
            Bamg_MetricFile!(solver.mesh, metric, MetricFile)
            
            #The optimal (theoretical) mesh defined above may not be 
            #Read mesh until number
            #Call bamg to refine. The new mesh is first written in .db.new file and 
            #then overwritten in the solver.MeshFile:
            #previous ratio 2.7
            run(pipeline(`bamg -AbsError -NoRescaling -NbJacobi 3 -NbSmooth 5 -ratio 2.7 -nbv $(1000000) -v 0 -M $(MetricFile) -b $(MeshFileDb) -o $(MeshFileDb).new`, 
                stdout="$(@__DIR__)/../../temp/BamgSC$(solver.SC).out", 
                stderr="$(@__DIR__)/../../temp/BamgSC$(solver.SC).err"))  
            mv(string(MeshFileDb,".new"), MeshFileDb, force=true)
            
#             println("BAMG=", time()-t_ini2)
#             t_ini2                  = time()
            
            #-----------------------------------------------------
            #NEW SOLUTION:
            
            #Read mesh and allocate memory:
            nElems0                 = solver.mesh.nElems
            ReadMesh!(solver)
#             figure()
#             PlotMesh!(solver.mesh)
    #         println("nElems=", solver.mesh.nElems, ", expected=", nElems_est)
        
#             println("Read mesh=", time()-t_ini2)
            t_ini2                  = time()
            
            #Interpolate initial condition (necessary for iterative solvers):
            if FunType==FW11
                Interpolate!(u0, solver.u, solver.fes)
            else
                
                #Destructure u0:
                u0_uv, u0_fes   = u0
                
                #Element guess for each mesh vertex:
                vert_elems  = solver.mesh.ElemsNodesSp.rowval[ 
                                solver.mesh.ElemsNodesSp.colptr[1:solver.mesh.nVerts]]
                #For each dof, ElemGuess is that corresponding to the first vertex:
                ElemGuess   = zeros(Int,solver.fes.PSpace.nDof)
                for iDof=1:solver.fes.PSpace.DofPerElem, iElem=1:solver.mesh.nElems
                    dof             = solver.fes.PSpace.ElemsDof[iElem,iDof]
                    vert1           = solver.fes.PSpace.ElemsDof[iElem,1]
                    ElemGuess[dof]  = min(u0_fes.mesh.nElems, vert_elems[vert1])
                end

                #Evaluate finite element solution:
                flag    = SolutionCompute!(u0_uv, u0_fes, [view(solver.fes.PSpace.NodesCoords,:,1), 
                                    view(solver.fes.PSpace.NodesCoords,:,2)], 
                                    GetPolynomials(solver.u, solver.fes),
                                ElemGuess=ElemGuess)
                if flag<=0
                    printstyled("Unable to interpolate initial condition in AMA algorithm\n", 
                        color=:light_yellow)
                    return -4, AMA_iters
                end
                
            end
            
#             println("SearchLocate=", time()-t_ini2)
#             t_ini2                  = time()
            
            #Project initial condition. Note: If eta_i is the Lq-error in variable "i"
            #(normalized with the domain size), we need that 
            #   eta^2 := 1/nVars * sum_i (eta_i/n_i)^2 <= TolA^2
            #A sufficient condition for this is 
            #   eta_i/n_i <= TolA
            #So, the real tolerance for variable "i" is 
            #   n_i*TolA
            etaA    = 0.0
            if true
                Integ2D = TrInt(solver.mesh, solver.AMA_ProjOrder, N=solver.AMA_ProjN)
                flag    = L2Projection!(Integ2D, u0, solver.fes, solver.b)
                if flag<=0
                    printstyled("Unable to compute L2 projection in AMA algorithm\n", 
                        color=:light_yellow)
                    return -2, AMA_iters
                end
                for ii=1:solver.nVars
                    flag, nIter, etaA_I = LS_gmres!(solver.MII_LS, solver.u[ii], solver.b[ii], 
                                            AbsTol=solver.TolA_min*solver.nFacts[ii], Display="notify")
                    etaA                += (etaA_I/solver.nFacts[ii])^2
                    if flag<=0
                        printstyled("Unable to solve L2 projection in AMA algorithm\n", 
                            color=:light_yellow)
                        return -3, AMA_iters
                    end
                end
                etaA                    = sqrt(etaA/solver.nVars)
            end
            
#             println("L2 projection=", time()-t_ini2)
#             t_ini2                  = time()
            
            #Estimate and save errors:
            etaS0                   = solver.etaS
            solver.etaS, 
                solver.etaS_elems,
                solver.urec         = etaS_reconstruction(solver, solver.u, 
                                        q=solver.SpaceNorm)
            #(solver.etaT is not modified)
            solver.etaA             = etaA
                
            #Update ratio TolS/etaS:
            solver.AMA_TolS_etaS    = AMA_TolS/solver.etaS
            
            #Save errors:
            push!(solver.tv, solver.t)
            push!(solver.etaSv, solver.etaS)
            push!(solver.etaTv, solver.etaT)
            push!(solver.etaAv, solver.etaA)
            push!(solver.nElemsv, solver.mesh.nElems)
            push!(solver.nDofv, solver.nVars*solver.fes.nDof)
            push!(solver.CFLv, NaN)
            #validv will be updated later:
            
            #Print info:
            printstyled(
                "AMA_iters=", AMA_iters, 
                ", nElems=", solver.mesh.nElems, 
                ", TotalDof=", solver.fes.nDof*solver.nVars, 
                ", TolS_min=", sprintf1("%.2e", TolS_min),
                ", etaS=", sprintf1("%.2e", solver.etaS), 
                ", TolS=", sprintf1("%.2e", AMA_TolS), 
                ", TolS_max=", sprintf1("%.2e", TolS_max),
                "\n")
            
            if false && solver.t+solver.Deltat>5e-4
            
                @warn "Plotting solution"
                fig     = PyPlotFigure(w=8.5, h=8.5, top=0.0, bottom=0.0, left=0.0, right=0.0)

                #Plot solution
                PlotContour(solver, solver.model, "rho", delta=1e-5)
                PlotMesh!(solver.mesh, color="k")
                xlim(0.0, 0.025)
                ylim(0.0, 0.025)
                savefig(string(SaveDir(), "Remap_$(AMA_iters).png"), pad_inches=0, dpi=800)
                
                close(fig)
                
            end
            
            #Exit DEq loop:
            if DEq_MaxIter==0
                break
            elseif 0.95*nElems0 <= solver.mesh.nElems <= 1.05*nElems0 && DEq_iters>1
                printstyled("Equidistribution reached\n", color=:cyan)
                break
            elseif DEq_iters==DEq_MaxIter
                printstyled("Reached maximum of equidistribution iterations in mesh adaptation in $(time()-t_ini) seconds\n", 
                    color=:light_yellow)
                break
            elseif AMA_iters==solver.AMA_MaxIter
                break
            else
                push!(solver.validv, 0)
            end
            
        end #DEq loop
        
        #Exit AMA_TolS loop if errors are under tolerance:
        if AMA_Strategy==1 && solver.etaS<=TolS_max
            push!(solver.validv, 1)
            printstyled("Mesh adaptation succeeded in $(time()-t_ini) seconds\n", color=:cyan)
            return 1, AMA_iters
        elseif AMA_Strategy==0 && solver.etaS>=TolS_min #(AMA_Strategy=0 implies etaS<=TolS_max)
            printstyled("Mesh adaptation succeeded in $(time()-t_ini) seconds\n", color=:cyan)
            push!(solver.validv, 1)
            return 1, AMA_iters
        elseif AMA_Strategy==0 && AMA_iters>=min(10, solver.AMA_MaxIter)
            printstyled("Mesh is not coarse enough, but error is under tolerance in $(time()-t_ini) seconds\n", color=:cyan)
            push!(solver.validv, 1)
            return 1, AMA_iters
        elseif AMA_iters>=solver.AMA_MaxIter
            printstyled("Reached maximum of iterations in mesh adaptation in $(time()-t_ini) seconds\n", 
                color=:light_yellow)
            push!(solver.validv, 0)
            return -1, AMA_iters
        else
            push!(solver.validv, 0)
        end
        
    end #AMA_TolS loop
    
    return
    
end

function etaS_reconstruction(solver::SolverData, u::Vector{<:GenVector{Float64}}; 
    q::Number=2)

    #Compute reconstruction:
    t_ini       = time()
    uP          = GetPolynomials(u, solver.fes)
    urec,       = AMA_BC_PatchReconstruction(uP, solver.fes.PSpace, solver.rec_fes)
    
    #Error = urec - projection(urec). Review this for curvilinear meshes:
    if solver.mesh.order>1
        @warn "Review computation of e_h. Should it be e_h = u_h^rec - MQQ^{-1}(MQP*u_h)??"
    end
    
    #Error = urec - L2proj(urec):
    erec                    = Vector{Vector{Float64}}(undef,solver.nVars)
    @inbounds for II=1:solver.nVars
        erec[II]            = zeros(size(urec[II]))
        for iElem=1:solver.mesh.nElems, iDof=solver.fes.PSpace.DofPerElem+1:solver.rec_fes.DofPerElem
            dof             = (iElem-1)*solver.rec_fes.DofPerElem + iDof
            erec[II][dof]   = urec[II][dof]
        end
    end
    
    #Lq norm of the error:
    etaS, etaS_elems        = LqMean(solver.Integ2D, erec, solver.rec_fes, 
                                        solver.nFacts, q=q)
    
    return etaS, etaS_elems, urec
    
end

#-----------------------------------------------------------------
#Post-time-step adaptation:

function LIRKHyp_Step_Post_Dolejsi!(solver::SolverData)

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
    
    #Set algebraic tolerance:
    TolA            = max(solver.TolA_min, solver.CA*min(etaS_np1_est, etaT_np1_est))
    
    #Declare some variables:
    Deltat_CFL      = 0.0
    u_np1           = zeros(0)
    urec_np1        = Vector{VectorView{Float64}}(undef,solver.nVars)
    etaS_np1        = 0.0
    etaS_elems_np1  = zeros(0)
    etaT_np1        = 0.0
    etaA_np1        = 0.0
    urec_n          = Vector{VectorView{Float64}}(undef,solver.nVars)
    etaS_n          = 0.0
    etaS_elems_n    = zeros(0)
    
    #Vector of optimal metrics and adaptation tolerances. Here the mesh is refined only.
    AMA_etaS        = zeros(solver.RK.stages)
    AMA_metric      = Vector{SpectralDecomp{Vector{Float64}, Int}}(undef, solver.RK.stages)
    AMA_TolS        = fill(NaN, solver.RK.stages)
    AMA_valid       = zeros(Int, solver.RK.stages)
    
    #Initialize, for first stage, optimal metric = current metric:
#     for ss=1
#         AMA_metric[ss]      = deepcopy(solver.mesh.metric)
#         AMA_TolS[ss]        = NaN
#         AMA_valid[ss]       = 1
#     end
    
    #Loop until space tolerance is achieved:
    RepeatS         = true
#     LSFailed_iters  = 0
    AMA_iters       = 0
    AMA_RefineFact  = 0.9
    AMA_CoarseFact  = 5.0
    while RepeatS
    
        #Flux and Jacobian for first stage:
        t_ini                   = time()
        println("Computing Jacobian")
        Deltat_CFL              = Rhs!(solver, t_n, solver.uv, true, 
                                    view(solver.f_RK,:,1), solver.Jm)
        solver.tJm              += time()-t_ini
        printstyled("Jacobian computed in ", time()-t_ini, " seconds, nnz=", 
            solver.Jm.colptr[solver.Jm.n+1], " \n", color=:cyan)
        view(solver.Ju_RK,:,1)  .= solver.Jm*solver.uv
    
        #Optimal metric for first stage:
        if AMA_iters==0
        
            etaS_n                          = solver.etaS
            etaS_elems_n                    = solver.etaS_elems
            urec_n                          = solver.urec
            
        elseif true && solver.SpaceAdapt
                
            kk                              = 1
            
            #Estimate space error:
            t_ini                           = time()
            etaS_n, etaS_elems_n, 
                urec_n                      = etaS_reconstruction(solver, 
                                                GetViews(solver.uv, solver.nVars, solver.fes.nDof), 
                                                q=solver.SpaceNorm)
#             solver.tetaS                    += time()-t_ini
            solver.tAMA                     += time()-t_ini 
            #(in the pre- approach, this reconstruction is performed in AdaptMesh! and the time goes to
            #tAMA, not to tetaS.)
            
        end
    
        #---------------------------------------------------------
        #Stages 2,...,s:
        
        #Allocate variables:
        u_np1               = zeros(size(solver.uv))
        b_1                 = solver.Mm*solver.uv
        
        #Loop stages until algebraic and time tolerances are met:
        RepeatTA            = true
        LSFailed_iters      = 0
        while RepeatTA #before this loop, Deltat_n, t_np1, TolA must be computed
        
            #First stage:
            @mlv u_np1      = solver.uv
            
            #Loop stages (solution for stage kk is overwritten in solver.t and solver.u):
            LSFlag          = 0
            solver.LS_iters = 0
            for kk=2:solver.RK.stages
            
                #Set time:
                t_k                     = t_n + solver.RK.c[kk]*Deltat_n
                
                #Compute system matrix and set up linear system, if necessary.
                if kk==2 || solver.RK.AI[kk,kk]!=solver.RK.AI[kk-1,kk-1]
                    println("Setting up linear system")
                    t_ini                   = time()
                    @mlv solver.Am.nzval    = solver.Mm.nzval - 
                                                (Deltat_n*solver.RK.AI[kk,kk])*solver.Jm.nzval
                    LinearSystem!(solver.Am_LS)
                    solver.tSCILU           += time()-t_ini
                    printstyled("Linear system set up in ", time()-t_ini, " seconds \n", color=:cyan)
                    #PlotLSPermutations(solver.Am_LS, solver.MII, solver.nVars)
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
    #                 display(norm(solver.f_RK[:,kk-1]))
                    display(norm(solver.f_RK))
                    display(norm(solver.Ju_RK))
                    error("NaNs in rhs vector")
                end
                solver.tb               += time()-t_ini
                
                #Solve:
                t_ini       = time()
                LSOutput    = LS_gmres!(solver.Am_LS, u_np1, solver.bv, AbsTol=TolA, Display="final",
                                MaxIter=solver.LS_iters_max)
                solver.tLS  += time()-t_ini
                LSFlag      = LSOutput[1]
                LSIter      = LSOutput[2]
                etaA_np1    = LSOutput[3]
                if LSFlag<=0
                    break #loop for kk=1:nStages
                end
                solver.LS_iters     += LSIter
                solver.LS_total     += LSIter
                        
                #Compute flux and save derivatives:
                t_ini                   = time()
                Rhs!(solver, t_k, u_np1, false, view(solver.f_RK,:,kk), solver.Jm)
                solver.tRhs             += time()-t_ini
                view(solver.Ju_RK,:,kk) .= solver.Jm*u_np1
                
                #Only for AMA:
                if false && kk==solver.RK.stages && solver.SpaceAdapt
                    
                    #Estimate space error:
                    t_ini                           = time()
                    etaS_np1, etaS_elems_np1, 
                        urec_np1                    = etaS_reconstruction(solver, 
                                                        GetViews(u_np1, solver.nVars, solver.fes.nDof), 
                                                        q=solver.SpaceNorm)
                    solver.tetaS                    += time()-t_ini
                    
#                     println("etaS_k=", etaS_np1)
                    
                    t_ini                           = time()
                    
                end
                
            end
        
            #If LS converged:
            if LSFlag>0
            
                #Mean number of LS iterations:
                solver.LS_iters     /= solver.RK.stages-1
                
                #Space errors:
                t_ini                           = time()
                etaS_np1, etaS_elems_np1, 
                    urec_np1                    = etaS_reconstruction(solver, 
                                                    GetViews(u_np1, solver.nVars, solver.fes.nDof), 
                                                    q=solver.SpaceNorm)
                solver.dlogetaS_dt              = (log(etaS_np1)-log(etaS_n))/Deltat_n
                solver.tetaS                    += time()-t_ini
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
                ehat                = 0.0*u_np1   #ehat=u-uhat
                ehat_views          = GetViews(ehat, solver.nVars, solver.fes.nDof)
                for II=1:solver.nVars
                    flag,           = LS_gmres!(solver.MII_LS, ehat_views[II], solver.b[II], 
                                        AbsTol=TolA*solver.nFacts[II], Display="notify")
                    if flag<=0
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
                                        solver.CA_max*min(etaS_np1, etaT_np1))
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
                    TolA            = max(solver.TolA_min, solver.CA*min(etaS_np1_est, etaT_np1_est))
                    
                    #Save iteration:
                    push!(solver.validv, false)
                    
                elseif solver.SpaceAdapt && max(etaS_n, etaS_np1)>solver.TolS_max
#                 elseif solver.SpaceAdapt && etaS_np1>solver.TolS_max
                
                    printstyled("Space error (", sprintf1("%.2e", max(etaS_n, etaS_np1)), 
                        ") larger than required (", sprintf1("%.2e", solver.TolS_max),
                        ") \n", color=:light_yellow)
                    
                    #Save iteration:
                    push!(solver.validv, false)
                    
                    #-------------------------------------------------------------------------------
                    #Optimal metric for solution at t^(n+1)-:
                    
                    kk                              = solver.RK.stages
                    
                    t_ini                           = time()
                    
                    #Set AMA tolerance:
                    if AMA_iters==0 #first value
#                         AMA_TolS[kk]    = max(solver.SfS*solver.TolS_max, AMA_RefineFact*etaS_np1)
#                         AMA_TolS[kk]    = max(solver.TolS_max, AMA_RefineFact*etaS_np1)
#                         AMA_TolS[kk]    = max(solver.SfS*solver.TolS_max, etaS_np1)
#                         AMA_TolS[kk]    = solver.SfS*solver.TolS_max
                        AMA_TolS[kk]    = AMA_RefineFact*etaS_np1
                    elseif etaS_np1>solver.TolS_max 
                        AMA_TolS[kk]    *= AMA_RefineFact
                    else    #otherwise, leave TolS_kk as is 
                    end
                    
                    #Optimal anisotropy:
                    emeanv0, emeanv, 
                        err0, err1, 
                        AMA_metric[kk]          = AMA_OptimalAniso!(urec_np1, solver.rec_fes,
                                                    solver.nFacts, max(1, solver.AMA_AnisoOrder), q=solver.SpaceNorm)
                    if solver.AMA_AnisoOrder==0
                        AMA_metric[kk].sf[1]    .= 1.0
                        AMA_metric[kk].sf[2]    .= 1.0
                        emeanv                  .= emeanv0
                    end
            
                    #Optimal sizes:
                    A1v                             = AMA_OptimalSize(etaS_elems_np1, 
                                                        solver.Av, fill(solver.AMA_SizeOrder+1, solver.mesh.nElems), 
                                                        fill(1,solver.mesh.nElems), AMA_TolS[kk], 
                                                        q=solver.SpaceNorm)
                    @mlv AMA_metric[kk].lambda_bar  = TrElem_Area/A1v
            
                    #Correct metric (impose maximum radius):
                    nElems_est              = AMA_CorrectMetric!(solver.Av, AMA_metric[kk])
                    
                    solver.tAMA             += time()-t_ini
                    
                    #-------------------------------------------------------------------------------
                    #Optimal metric for solution at t^n+:
                    
                    kk                              = 1
                    
                    t_ini                           = time()
                    
                    #Set AMA tolerance:
                    if AMA_iters==0 #first value
#                         AMA_TolS[kk]    = max(solver.SfS*solver.TolS_max, AMA_RefineFact*etaS_n)
#                         AMA_TolS[kk]    = max(solver.SfS*solver.TolS_max, etaS_n)
#                         AMA_TolS[kk]    = max(solver.TolS_max, AMA_RefineFact*etaS_np1) #Same tolerance for both meshes
#                         AMA_TolS[kk]    = solver.TolS_max
                        AMA_TolS[kk]    = AMA_TolS[solver.RK.stages]
                    elseif etaS_n>solver.TolS_max 
                        AMA_TolS[kk]    *= AMA_RefineFact
                    else    #otherwise, leave TolS_kk as is
                    end
                    
                    #Optimal anisotropy:
                    emeanv0, emeanv, 
                        err0, err1, 
                        AMA_metric[kk]          = AMA_OptimalAniso!(urec_n, solver.rec_fes,
                                                    solver.nFacts, max(1, solver.AMA_AnisoOrder), 
                                                    DeltaAR=10.0, q=solver.SpaceNorm)
                    if solver.AMA_AnisoOrder==0
                        AMA_metric[kk].sf[1]    .= 1.0
                        AMA_metric[kk].sf[2]    .= 1.0
                        emeanv                  .= emeanv0
                    end
            
                    #Optimal sizes:
                    A1v                             = AMA_OptimalSize(etaS_elems_n, 
                                                        solver.Av, fill(solver.AMA_SizeOrder+1, solver.mesh.nElems), 
                                                        fill(1,solver.mesh.nElems), AMA_TolS[kk], 
                                                        q=solver.SpaceNorm)
                    @mlv AMA_metric[kk].lambda_bar  = TrElem_Area/A1v
            
                    #Correct metric (impose maximum radius):
                    nElems_est              = AMA_CorrectMetric!(solver.Av, AMA_metric[kk])
                    
                    solver.tAMA             += time()-t_ini
                    
                    #-------------------------------------------------------------------------------
                    #Generate new mesh:
                    
                    #Print mesh info:
                    printstyled("AMA_iters=", AMA_iters, ", nElems=", solver.mesh.nElems, 
                            ", nNodes=", solver.mesh.nNodes, ", TotalDof=", solver.fes.nDof*solver.nVars, 
                            ", TolS_min=", sprintf1("%.2e", solver.TolS_min),
                            ", etaS=", sprintf1("%.2e", etaS_np1), 
                            ", TolS=", sprintf1("%.2e", AMA_TolS[solver.RK.stages]), 
                            ", TolS_max=", sprintf1("%.2e", solver.TolS_max),
                            ", TolS/etaS=", sprintf1("%.2e", solver.AMA_TolS_etaS),
                            "\n")
                
                    #Exit if maximum iteration number has been reached:
                    if AMA_iters==solver.AMA_MaxIter
                        printstyled("Reached maximum mesh adaptation iterations\n", color=:magenta)
                        return -3
                    end
                
                    #Update correction factor (we take the last stage as reference):
                    if AMA_iters>0
                        solver.AMA_TolS_etaS    = etaS_np1/AMA_TolS[solver.RK.stages]
                    end
                    
                    t_ini       = time()
                    
                    #Compute mesh intersection and adapt (refine) mesh until etaS_[1] is 
                    #under tolerance:
                        
                    #Metric intersection:
                    metric      = AMA_MetricIntersection(AMA_metric[[1,solver.RK.stages]])
                    
                    #Call BAMG:
                    MetricFile  = "$(@__DIR__)/../../temp/BamgSC$(solver.SC).mtr"
                    MeshFileDb  = "$(@__DIR__)/../../temp/BamgSC$(solver.SC).db"
                    Bamg_MetricFile!(solver.mesh, metric, MetricFile)
                    run(pipeline(`bamg -AbsError -NoRescaling -NbJacobi 3 -NbSmooth 5 -ratio 2.7 -nbv $(1000000) -v 0 -M $(MetricFile) -b $(MeshFileDb) -o $(MeshFileDb).new`, 
                        stdout="$(@__DIR__)/../../temp/BamgSC$(solver.SC).out", 
                        stderr="$(@__DIR__)/../../temp/BamgSC$(solver.SC).err"))  
                    mv(string(MeshFileDb,".new"), MeshFileDb, force=true)
                    ReadMesh!(solver)
                    println("BAMG called")
                
                    #Project "initial condition" u_n. The solution is overwritten in solver.u.
                    #Note, however, that solver.u is a new pointer after reading the mesh. Hence, the 
                    #previous solution is still stored in u_n.
                    L2Projection!(solver, u_n, fes_n)
                    println("L2 projection done")
                    
                    solver.tAMA     += time()-t_ini
                    
                    t_ini           = time()
                    
                    #Allocate variables for time integration:
                    RKAlloc!(solver)
                    println("RK vars allocated")
                    
                    solver.tAlloc   += time()-t_ini
                    
                    #Update mesh iterations:
                    AMA_iters       += 1
                    
                    #Adjust time step:
                    Deltat_n        = Deltat_nm1*solver.AMA_rDeltat
                    t_np1           = min(t_n + Deltat_n, solver.tf)
                    Deltat_n        = t_np1-t_n
                    
                    #Next iteration:
                    RepeatTA        = false
                    RepeatS         = true
                    etaS_np1_est    = etaS_n * exp(solver.dlogetaS_dt * Deltat_n)
                    etaT_np1_est    = etaT_n * Deltat_n^(solver.RK.order) #does not change
                    TolA            = max(solver.TolA_min, solver.CA*min(etaS_np1_est, etaT_np1_est))
                
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
                
                printstyled("Linear solver did not converge. etaA=", sprintf1("%.2e", etaA_np1), 
                    ", TolA=", sprintf1("%.2e", TolA), 
                    ". Decreasing time step\n", color=:light_yellow)
                if LSFailed_iters==10
                    printstyled("Linear solver unable to converge after ",
                        LSFailed_iters, " attempts. Aborting\n", color=:magenta)
                    RepeatTA            = false
                    RepeatS             = false
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
                    TolA                = max(solver.TolA_min, 
                                            solver.CA*min(etaS_np1_est, etaT_np1_est))
                                        
                end
                
            end #if LSFlag>0
            
        end #RepeatTA
        
    end #RepeatS
    
    #Time-algebraic loop finished. If method has converged...
    if ConvFlag==1
        
        #Save new solution:
        solver.t            = t_np1
        solver.Nt           += 1
        solver.Deltat       = Deltat_n
        solver.CFL          = Deltat_n/Deltat_CFL
        solver.uv           .= u_np1
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
            ", tCPU=", sprintf1("%.2e", solver.tCPU + (time()-t_start)))
        
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

function LIRKHyp_Step_Post!(solver::SolverData)
    return LIRKHyp_Step_Post_Dolejsi!(solver)
end

function L2Projection!(solver::SolverData, u0::Vector{Float64}, u0_fes::TrPBSpace)

    u0_uv           = GetViews(u0, solver.nVars, u0_fes.nDof)
    
    #Element guess for each mesh vertex:
    vert_elems  = solver.mesh.ElemsNodesSp.rowval[ 
                    solver.mesh.ElemsNodesSp.colptr[1:solver.mesh.nVerts]]
                
    #For each dof, ElemGuess is that corresponding to the first vertex:
    ElemGuess   = zeros(Int,solver.fes.PSpace.nDof)
    for iDof=1:solver.fes.PSpace.DofPerElem, iElem=1:solver.mesh.nElems
        dof             = solver.fes.PSpace.ElemsDof[iElem,iDof]
        vert1           = solver.fes.PSpace.ElemsDof[iElem,1]
        ElemGuess[dof]  = min(u0_fes.mesh.nElems, vert_elems[vert1])
    end

    #Evaluate finite element solution:
    flag    = SolutionCompute!(u0_uv, u0_fes, [view(solver.fes.PSpace.NodesCoords,:,1), 
                        view(solver.fes.PSpace.NodesCoords,:,2)], 
                        GetPolynomials(solver.u, solver.fes),
                    ElemGuess=ElemGuess)
    if flag<=0
        printstyled("Unable to interpolate initial condition in AMA algorithm\n", 
            color=:light_yellow)
        return -1
    end
    
    #Project initial condition. Note: If eta_i is the Lq-error in variable "i"
    #(normalized with the domain size), we need that 
    #   eta^2 := 1/nVars * sum_i (eta_i/n_i)^2 <= TolA^2
    #A sufficient condition for this is 
    #   eta_i/n_i <= TolA
    #So, the real tolerance for variable "i" is 
    #   n_i*TolA
    etaA    = 0.0
    if true
        Integ2D = TrInt(solver.mesh, solver.AMA_ProjOrder, N=solver.AMA_ProjN)
        flag    = L2Projection!(Integ2D, (u0_uv, u0_fes), solver.fes, solver.b)
        if flag<=0
            printstyled("Unable to compute L2 projection in AMA algorithm\n", 
                color=:light_yellow)
            return -1
        end
        for ii=1:solver.nVars
            flag, nIter, etaA_I = LS_gmres!(solver.MII_LS, solver.u[ii], solver.b[ii], 
                                    AbsTol=solver.TolA_min*solver.nFacts[ii], Display="notify")
            etaA                += (etaA_I/solver.nFacts[ii])^2
            if flag<=0
                printstyled("Unable to solve L2 projection in AMA algorithm\n", 
                    color=:light_yellow)
                return -1
            end
        end
        etaA                    = sqrt(etaA/solver.nVars)
    end
   
    return 1
    
end


#-----------------------------------------------------------------
#Choose default step function:

function LIRKHyp_Step!(solver::SolverData)
    return LIRKHyp_Step_Pre!(solver)
end

#-----------------------------------------------------------------
#Implicit RK solver:

function IRK_Step!(solver::SolverData)
    
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
        u_n_views       = solver.u
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
    TolA            = max(solver.TolA_min, solver.CA*min(etaS_np1_est, etaT_np1_est, etaS_n))
    
    #Flux and Jacobian for first stage:
    t_ini                   = time()
    println("Computing Jacobian")
    Deltat_CFL              = Rhs!(solver, t_n, u_n, true, 
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
    etaA_np1            = 0.0
    
    #Save M*u[1]:
    b1                  = solver.Mm*u_n

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
                @mlv solver.Am.nzval    = solver.Mm.nzval - 
                                            (Deltat_n*solver.RK.AI[kk,kk])*solver.Jm.nzval
                LinearSystem!(solver.Am_LS)                
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
            if any(isnan.(solver.bv))
                display(kk)
                display(norm(b_1))
                display(norm(solver.f_RK))
                display(norm(solver.Ju_RK))
                error("NaNs in solver.bv")
            end
            solver.tb               += time()-t_ini
            
            #Define preconditioned residual:
            etaA_np1        = NaN
            function QNResidual1!(u::Vector{Float64}, gres::Vector{Float64})
                
                #Compute flux and save derivatives:
#                 t_ini                   = time()
                Rhs!(solver, t_k, u, false, view(solver.f_RK,:,kk), solver.Jm)
#                 solver.tRhs             += time()-t_ini
                
                #Compute residual f:
                #   f = D*(M*(Du) - b - Deltat*a_kk*f_k(Du))
                fres        = solver.Mm*u - solver.bv - (Deltat_n*solver.RK.AI[kk,kk])*
                                    view(solver.f_RK,:,kk)
                
                #Compute preconditioned residual:
                LSOutput    = LS_gmres!(solver.Am_LS, gres, fres, 
#                                 RelTol=1e-2, AbsTol=0.0, 
                                RelTol=0.0, AbsTol=TolA, 
                                Display="notify", MaxIter=solver.LS_iters_max)
                #Returns flag, nIters, etaA
                flag        = LSOutput[1]
                etaA_np1    = LSOutput[3]
                
                return flag
                
            end
            
            #Solve:
            t_ini       = time()
            LSOutput    = Anderson(FW_NLS((u,gres)->QNResidual1!(u,gres)), 
                            u_k, 
                            AbsTolX=1.0*sqrt(length(u_k))*TolA, RelTolX=0.0, 
                            AbsTolG=0.0*sqrt(length(u_k))*TolA, RelTolG=0.0, 
                            memory=100, MaxIter=solver.LS_iters_max, Display="final")
            solver.tLS  += time()-t_ini
            u_k         .= LSOutput[1]
            LSFlag      = LSOutput[2].flag
            if LSFlag<=0
                @warn "LSFlag"
                LSFlag      = 1
                etaA_np1    = 0.0 
            end
            LSIter      = LSOutput[2].nIter
            #etaA_np1 is computed in the last call to QNResidual1
            if LSFlag<=0
#                 printstyled("Linear solver did not converge\n", color=:light_yellow)
                break #loop for kk=1:nStages
            else
#                 println("Stage $kk, linear solver converged in $(LSIter) iterations")
            end
            solver.LS_iters     += LSIter
            solver.LS_total     += LSIter
            #Here, LS in reality is NLS
            
#             splot_fun(x1,x2)    = @mlv x1
#             solver.uv           .= u_k
#             for var in ["h", "v1", "v3", "p"]
#                 figure()
#                 PlotNodes(splot_fun, solver, var)
#             end
#             error("")
            
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
