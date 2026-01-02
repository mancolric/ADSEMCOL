#Project function for given FES:
function L2Projection!(Integ2D::TrInt, fun::FunType, FES::TrFES, 
    bv::Vector{<:GenVector{Float64}}) where FunType<:Union{FW11, Tuple{Vector{<:GenVector{Float64}}, <:TrFES}}

    t_ini           = time()
    
    flag            = 1
    
    nElems          = Integ2D.mesh.nElems
    nqp             = Integ2D.QRule.nqp
    xim             = Integ2D.QRule.xi
    
    #Eval function at quadrature nodes:
    x_qp            = QuadNodesCoords(Integ2D)
    fun_qp          = Vector{Matrix{Float64}}
    if FunType==FW11
    
        #Evaluate function:
        fun_qp      = fun(x_qp)
        
    else
    
        #Destructure donor u and fes:
        fun_uv, fun_fes     = fun
        
        #Find barycenters of new mesh in old mesh:
        xbarm               = barycenters(Integ2D.mesh)
        xbar                = [ view(xbarm,:,1), view(xbarm,:,2) ]
        ElemGuess_xbar      = min.(fun_fes.mesh.nElems, 1:nElems)
#         _, Elem0_xbar,      = SearchLocate(Integ2D.mesh, xbar; ElemGuess=ElemGuess_xbar)
        _, Elem0_xbar,      = SearchLocate(fun_fes.mesh, xbar; ElemGuess=ElemGuess_xbar)
        
#         println("ElemGuess_0 ", time()-t_ini)
        
        #DEBUG:
        if any(Elem0_xbar.<=0)
            aux         = findall(Elem0_xbar.<=0)
            figure()
            plot(xbar[1][aux], xbar[2][aux], "*r")
            PlotMesh!(fun_fes.mesh, color="b")
            PlotMesh!(Integ2D.mesh, color="k")
            save("../temp/test.jld2", "mesh0", fun_fes.mesh, "mesh1", Integ2D.mesh,
                "xbar", [xbar[1][aux], xbar[2][aux]], "ElemGuess", ElemGuess_xbar[aux]) 
            error("")
        end
        
        #Element guess for each qp:
#         ElemGuess   = min.(fun_fes.mesh.nElems, reshape(Vector(1:nElems)*ones(Int,1,nqp),:))
        ElemGuess   = reshape(Elem0_xbar * ones(Int, 1, nqp), :)
        
        #Evaluate finite element solution:
        u_qp, flag  = SolutionCompute(fun_uv, fun_fes, 
                        [reshape(x_qp[1],:), reshape(x_qp[2],:)], 
                        ElemGuess=ElemGuess)
        fun_qp      = Vector{Matrix{Float64}}(undef, length(u_qp))
        for II=1:length(u_qp)
            fun_qp[II]  = reshape(u_qp[II], nElems, nqp)
        end
        
#         println("SolutionCompute_0 ", time()-t_ini)
            
    end
    
    #Load shape functions:
    Nm              = NCompute(FES, xim)
    
    #Create matrices with nnz elements and assemble:
    for ii=1:length(bv)
        omega_qp    = @mlv    fun_qp[ii]*Integ2D.wdetJ
        sm          = omega_qp*Nm
        VectorAssemble!(FES.ElemsDof[:], sm[:], bv[ii])
    end

    return flag
    
end

#Compute Lq error from given theoretical solution:
function LqError(solver::SolverData, utheor::FW11; q::Real=2)

    return LqError(solver.Integ2D, solver.u, solver.fes, utheor, solver.nFacts, q=q)
    
end
    
#Compute Lq error from given theoretical solution:
function LqError(Integ2D::TrInt, u::Vector{<:GenVector{Float64}}, fes::TrFES,
    utheor::FW11, nFacts::Vector{Float64} ; q::Real=2)
    
    nElems          = Integ2D.mesh.nElems
    nqp             = Integ2D.QRule.nqp
    xim             = Integ2D.QRule.xi
    wdetJ           = Integ2D.wdetJ
    nVars           = length(u)
    Av              = sum(Integ2D.wdetJ, dims=2)[:]
    Omega           = sum(Av)
    
    #Compute theoretical solution at quadrature nodes:
    x_qp            = QuadNodesCoords(Integ2D)
    utheor_qp       = utheor(x_qp)   #Vector{Matrix{Float64}}
    
    #Compute numerical solution at quadrature nodes:
    unum_qp         = SolutionCompute(u, fes, xim)
    
#     figure()
#     PyPlot.using3D()
#     plot3D(x_qp[1][:], x_qp[2][:], abs.(unum_qp[1]-utheor_qp[1])[:], ".b")
    
    #Global error at each qp. Recall that the global error is defined by:
    #   eta^q       = 1/nVars * sum_i (eta_i/n_i)^q
    #where eta_i is the Lq error in the variable "i",
    #   eta_i^q     = int( |u_i-uh_i|^q dOmega ) / Omega
    #So,
    #   eta^q       = 1/nVars * sum_i int( (|u_i-uh_i|/n_i)^q dOmega ) / Omega
    #               = int( 1/nVars sum_i (|u_i-uh_i|/n_i)^q dOmega ) / Omega
    if q==Inf
        
        #Then eta   = max_i |u_i-uh_i|/n_i
        e_Lq_elems          = zeros(nElems)
        for ii=1:nVars
            e_Lq_elems      = max.(e_Lq_elems, 
                                maximum(abs.(unum_qp[ii]-utheor_qp[ii])/nFacts[ii], dims=2)[:])
        end
        e_Lq                = maximum(e_Lq_elems)
        
        return e_Lq, e_Lq_elems
        
    else
        #Compute 1/nVars * sum_i (|u_i-uh_i|/n_i)^q:
        e_qp                = zeros(nElems, nqp)
        for ii=1:nVars
            @mlv    e_qp    += (abs(unum_qp[ii]-utheor_qp[ii])/nFacts[ii])^q
        end
        @mlv    e_qp        /= nVars
        
        #Errors at each element:
        @mlv    e_qp        *= wdetJ        #This is e^q dOmega
        e_Lq_elems          = (sum(e_qp, dims=2)[:]./Av).^(1.0/q)
        e_Lq                = (sum(e_qp)/Omega)^(1.0/q)
        
        return e_Lq, e_Lq_elems
        
    end
    
    return
    
end

#Compute Lq norm of finite element solution:
function LqMean(Integ2D::TrInt, u::Vector{<:GenVector{Float64}}, fes::TrFES,
    nFacts::Vector{Float64}; q::Real=2)
    
    nElems          = Integ2D.mesh.nElems
    nqp             = Integ2D.QRule.nqp
    xim             = Integ2D.QRule.xi
    wdetJ           = Integ2D.wdetJ
    nVars           = length(u)
    
    #Compute numerical solution at quadrature nodes:
    unum_qp         = SolutionCompute(u, fes, xim)
    
    #Global norm at each qp. Recall that the global norm is defined by:
    #   eta^q       = 1/nVars * sum_i (eta_i/n_i)^q
    #where eta_i is the Lq norm in the variable "i",
    #   eta_i^q     = int( |u_i-uh_i|^q dOmega ) / Omega
    #So,
    #   eta^q       = 1/nVars * sum_i int( (|u_i-uh_i|/n_i)^q dOmega ) / Omega
    #               = int( 1/nVars sum_i (|u_i-uh_i|/n_i)^q dOmega ) / Omega
    if q==Inf
        
        #Then eta   = max_i |u_i-uh_i|/n_i
        u_Lq_elems          = zeros(nElems)
        for ii=1:nVars
            u_Lq_elems      = max.(u_Lq_elems, 
                                maximum(abs.(unum_qp[ii])/nFacts[ii], dims=2)[:])
        end
        u_Lq                = maximum(u_Lq_elems)
        
        return u_Lq, u_Lq_elems
        
    else
    
        #Compute 1/nVars * sum_i (|u_i-uh_i|/n_i)^q:
        u_qp                = zeros(nElems, nqp)
        for ii=1:nVars
            @mlv    u_qp    += (abs(unum_qp[ii])/nFacts[ii])^q
        end
        @mlv    u_qp        /= nVars
        
        #Areas of each element and total area:
        A_elems             = sum(wdetJ, dims=2)[:]
        A_total             = sum(A_elems)
        
        #Errors at each element:
        @mlv    u_qp        *= wdetJ        #This is e^q dOmega
        u_Lq_elems          = (sum(u_qp, dims=2)[:]./A_elems).^(1.0/q)
        u_Lq                = (sum(u_qp)/A_total)^(1.0/q)
        
        return u_Lq, u_Lq_elems
        
    end
    
    return
    
end

#Compute L2 norm of one variable. Here MII is the mass matrix of a single block.
function SolutionNorm(MII::SparseMatrixCSC{Float64,Int}, x::GenVector{Float64}, 
    Omega::Float64)

    return sqrt(dot(x,MII,x)/Omega)
    
end

#Compute combined L2 norm of solution. Here Mm is the mass matrix of a single block.
function SolutionNorm(MII::SparseMatrixCSC{Float64,Int}, x::Vector{Float64}, 
    Omega::Float64, nFacts::Vector{Float64})
    
    nVars       = length(nFacts)
    xI          = GetViews(x, nVars, size(MII,2))
    SolNorm     = 0.0
    for II=1:nVars
        SolNorm     += dot(xI[II], MII, xI[II])/(nFacts[II]*nFacts[II])
    end
    SolNorm     = sqrt(SolNorm/(nVars*Omega))
    return SolNorm
    
end

function LinearSystem(A::SparseMatrixCSC{Float64,Int}, solver::SolverData,
    nFacts::Vector{Float64}; LSType::String="SCILU0")
    
    if uppercase(LSType)=="SCILU0"
        return LinearSystem_Complete(A, solver, nFacts)
    elseif uppercase(LSType)=="BLOCKJACOBI" || uppercase(LSType)=="BLOCKJACOBIKRYLOV"
        return LinearSystem_BlockJacobi(A, solver, nFacts)
    elseif uppercase(LSType)=="DOMINANTILU0"
        return LinearSystem_DominantILU0(A, solver, nFacts)
    else
        error("Unknown option")
    end
    
end
    
function LinearSystem_Complete(A::SparseMatrixCSC{Float64,Int}, solver::SolverData,
    nFacts::Vector{Float64})

    t_ini           = time()
    
    #Create linear system:
    LS              = LS_SCILU0_GMRES()
    
    #Save pointer to original matrix:
    LS.A            = A
    
    #Nb of variables, master dof and slave dof:
    nVars           = length(nFacts)
    nMaster         = maximum(view(solver.fes.ElemsDof,:,1:3*solver.fes.order))
    nDof            = solver.fes.nDof
    nSlaves         = nDof-nMaster
    
    #Optimal ordering for the degrees of freedom:
    dof_RCM_inv                 = zeros(Int, nDof)
    dof_RCM_inv[nMaster+1:nDof] = 1:nSlaves
    #Master dof are permuted according to symmetric Cuthill-McKee:
    dof_RCM_inv[1:nMaster]      = nSlaves .+ invperm(symrcm(solver.MII[1:nMaster,1:nMaster]))
    
    #For a single variable, dof "idof" goes to new position jdof=dof_RCM_inv[idof].
    #On the other hand, dof jdof of vble II goes to position
    #   (jdof-1)*nVars + II:
    LS.pinv                 = zeros(Int, nVars*nDof)
    LS.scaleP               = zeros(Float64, nVars*nDof)
    @inbounds for II=1:nVars, ii=1:nDof
        ii_new                      = dof_RCM_inv[ii]
        iiII_new                    = (ii_new-1)*nVars+II
        LS.pinv[(II-1)*nDof+ii]     = iiII_new
        LS.scaleP[iiII_new]         = nFacts[II]
    end
    LS.p            = invperm(LS.pinv)
    LS.scaleP_m     = LS.scaleP[nSlaves*nVars+1:(nSlaves+nMaster)*nVars]
    
    #Allocate permuted matrix and ILU factorization:
    iv,jv,          = findnz(A)
    S_aux           = sparse(LS.pinv[iv], LS.pinv[jv], 1:length(iv))
    LS.ssPP_ss      = S_aux.nzval
    LS.APP          = SparseMatrixCSC{Float64,Int}(S_aux.m, S_aux.n, S_aux.colptr, 
                        S_aux.rowval, A.nzval[LS.ssPP_ss])
    LS.Pl           = SCILU0_alloc(LS.APP, nSlaves*nVars)
    
    return LS
    
end

function LinearSystem_BlockJacobi(A::SparseMatrixCSC{Float64,Int}, solver::SolverData,
    nFacts::Vector{Float64})

    t_ini           = time()
    
    #Create linear system:
    LS              = LS_BlockJacobi()
    
    #Save pointer to original matrix:
    LS.A            = A
    
    #Nb of variables, master dof and slave dof:
    nVars           = length(nFacts)
    nMaster         = 0                 #no master dofs, system is solved directly
    nDof            = solver.fes.nDof
    nSlaves         = nDof-nMaster
    
    #Due to the structure of the matrix, it is not necessary to permute.
    #Nevertheless, variables are permuted so as to cluster more elements close to the 
    #diagonal (mayble this improves performance):
    dof_RCM_inv             = 1:nDof
    
    #For a single variable, dof "idof" goes to new position jdof=dof_RCM_inv[idof].
    #On the other hand, dof jdof of vble II goes to position
    #   (jdof-1)*nVars + II:
    LS.pinv                 = zeros(Int, nVars*nDof)
    LS.scaleP               = zeros(Float64, nVars*nDof)
    @inbounds for II=1:nVars, ii=1:nDof
        ii_new                      = dof_RCM_inv[ii]
        iiII_new                    = (ii_new-1)*nVars+II
        LS.pinv[(II-1)*nDof+ii]     = iiII_new
        LS.scaleP[iiII_new]         = nFacts[II]
    end
    LS.p            = invperm(LS.pinv)
    
    #Allocate permuted matrix and ILU factorization:
    iv,jv,          = findnz(A)
    S_aux           = sparse(LS.pinv[iv], LS.pinv[jv], 1:length(iv))
    LS.ssPP_ss      = S_aux.nzval
    LS.APP          = SparseMatrixCSC{Float64,Int}(S_aux.m, S_aux.n, S_aux.colptr, 
                        S_aux.rowval, A.nzval[LS.ssPP_ss])
    LS.Pl           = SCILU0_alloc(LS.APP, nSlaves*nVars)

    return LS
    
end

function LinearSystem_DominantILU0(A::SparseMatrixCSC{Float64,Int}, solver::SolverData,
    nFacts::Vector{Float64})

    t_ini           = time()
    
    #Create linear system:
    LS              = LS_DominantILU0()
    
    #Save pointer to original matrix:
    LS.A            = A
    
    #Nb of variables, master dof and slave dof:
    nVars           = length(nFacts)
    nMaster         = solver.fes.nDof                 #no master dofs, system is solved directly
    nDof            = solver.fes.nDof
    nSlaves         = 0
    
    #Optimal ordering for the degrees of freedom:
    dof_RCM_inv             = invperm(symrcm(solver.Mm[1:nMaster,1:nMaster]))
    
    #For a single variable, dof "idof" goes to new position jdof=dof_RCM_inv[idof].
    #On the other hand, dof jdof of vble II goes to position
    #   (jdof-1)*nVars + II:
    LS.pinv                 = zeros(Int, nVars*nDof)
    LS.scaleP               = zeros(Float64, nVars*nDof)
    @inbounds for II=1:nVars, ii=1:nDof
        ii_new                      = dof_RCM_inv[ii]
        iiII_new                    = (ii_new-1)*nVars+II
        LS.pinv[(II-1)*nDof+ii]     = iiII_new
        LS.scaleP[iiII_new]         = nFacts[II]
    end
    LS.p            = invperm(LS.pinv)
    
    #Allocate permuted matrix and ILU factorization:
    iv,jv,          = findnz(A)
    S_aux           = sparse(LS.pinv[iv], LS.pinv[jv], 1:length(iv))
    LS.ssPP_ss      = S_aux.nzval
    LS.APP          = SparseMatrixCSC{Float64,Int}(S_aux.m, S_aux.n, S_aux.colptr, 
                        S_aux.rowval, A.nzval[LS.ssPP_ss])
    LS.Pl           = SCILU0_alloc(LS.APP, nSlaves*nVars)

    return LS
    
end

function LinearSystem!(LS::LS_SCILU0_GMRES)

#     t_ini           = time()
    
    #Save data:
#     LS.APP.nzval    .= LS.A.nzval[LS.ssPP_ss]
#     println("Permuting matrix: ", time()-t_ini)
    @inbounds for ii=1:length(LS.ssPP_ss)
        LS.APP.nzval[ii]    = LS.A.nzval[LS.ssPP_ss[ii]]
    end
#     println("Permuting matrix: ", time()-t_ini)
    
    #Update preconditioner:
    SCILU0!(LS.Pl, LS.APP)
#     println("Updating preconditioner: ", time()-t_ini)
    
    #Set norm of master dof vector for stopping condition:
    aux_norm        = sqrt(LS.Pl.nMasters)
    LS.NormFun      = FW_NLS_norm((x)->norm(x)/aux_norm)    #Impose L2 average of the master dof < Tol
 
    return LS
    
end

function LinearSystem!(LS::LS_BlockJacobi)

#     t_ini           = time()
    
    #Save data:
    @inbounds for ii=1:length(LS.ssPP_ss)
        LS.APP.nzval[ii]    = LS.A.nzval[LS.ssPP_ss[ii]]
    end
#     println("Permuting matrix: ", time()-t_ini)
    
    #Update ILU0 factorization (exact for block Jacobi):
    SCILU0!(LS.Pl, LS.APP)
#     println("Updating factorization: ", time()-t_ini)
 
    return LS
    
end

function LinearSystem!(LS::LS_DominantILU0)

#     t_ini           = time()
    
    #Save data:
    @inbounds for ii=1:length(LS.ssPP_ss)
        LS.APP.nzval[ii]    = LS.A.nzval[LS.ssPP_ss[ii]]
    end
#     println("Permuting matrix: ", time()-t_ini)
    
    #Update ILU0 factorization:
    SCILU0!(LS.Pl, LS.APP)
#     println("Updating factorization: ", time()-t_ini)
    
    return LS
    
end

#NOTE: A refers to the exact Jacobian matrix M - a_ii Deltat_n J

#Solve A*u=b. Employ permutations, SCILU0 preconditioner and call to GMRES:
function LS_solve!(LS::LS_SCILU0_GMRES, u::GenVector{Float64}, b::GenVector{Float64};
    RelTol::Float64=0.0, AbsTol::Float64=1e-8, MaxIter::Int=100, Display::String="notify")
    
    #Notation:
    #u, b: original vectors:
    #uhat, bhat: permuted vectors. uhat is also scaled by LS.scalP
    #bred: statically condensed permuted vector 
    
    #Exit early if there are NaNs or Infs in b:
    if !all(isfinite.(b))
        return -1, 0, Inf   #flag, nIters, etaA
    end
        
    #Extract variables:
    F               = LS.Pl
    scaleP          = LS.scaleP
    scaleP_m        = LS.scaleP_m
    
    #Permute vectors:
    bhat            = b[LS.p]
    uhat            = u[LS.p]
    
    #Reduce r.h.s.:
    bred            = ReduceRhs(F, bhat)
    bred_m          = bred[F.nSlaves+1:F.N]
    uhat_m          = uhat[F.nSlaves+1:F.N]
    
    #Exit early if there are no master dofs:
    if F.nMasters==0
        ldiv_slave!(uhat, F, bred)
        @mlv u[LS.p]    = uhat
        return 1, 0, 0.0   #flag, nIters, etaA
    end
    
    #DEBUG:
    t_mult          = 0.0
    t_ldiv          = 0.0
    #ILU0 (quasi-Newton) iteration. Compute residual f=A_mm x_m-b_m and preconditioned residual g=approx(A)\f:
    f               = zeros(length(uhat_m));
    function QN_ILU0!(x::Vector{Float64}, g::Vector{Float64})
        t_ini       = time()
        @avxt @. x  *= scaleP_m
        mul!(f, F.Ared_mm, x)              #f = 1.0*A*x
        axpby!(-1.0, bred_m, 1.0, f)
        t_mult      += time()-t_ini
        if !isnothing(g)
            t_ini       = time()
            ldiv_master!(g, F, f)
            @avxt @. g  /= scaleP_m     #g ~~ Delta u. We have to return the scaled g
            t_ldiv      += time()-t_ini
        end
        @avxt @. x      /= scaleP_m
        return 1
    end
    
    #Solve with GMRES:
    flag                    = 1
    solver_output           = NLS_gmres(FW_NLS((x,g)->QN_ILU0!(x,g)), uhat_m./scaleP_m,
                                memory=100, MaxIter=MaxIter, 
                                AbsTolX=0.0, RelTolX=0.0, 
                                AbsTolG=AbsTol, RelTolG=RelTol,
                                Display=Display, history=true,
                                NormFun=LS.NormFun)
#     solver_output           = NLS_gmres(FW_NLS((x,g)->QN_ILU0!(x,g)), uhat_m./scaleP_m,
#                                 memory=100, MaxIter=MaxIter, 
#                                 AbsTolX=AbsTol, RelTolX=RelTol,
#                                 AbsTolG=0.0, RelTolG=0.0, 
#                                 Display="iter", history=true,
#                                 NormFun=LS.NormFun)
    @mlv uhat[F.nSlaves+1:F.N]  = solver_output[1]*scaleP_m
    ldiv_slave!(uhat, F, bred)
    @mlv u[LS.p]                = uhat
    nIter                       = solver_output[2].nIter 
    resv                        = solver_output[2].gnorms
#     resv                        = solver_output[2].pnorms
    etaA                        = resv[length(resv)]
    if solver_output[2].flag==2
        #Exit due to convergence in residual g
    elseif solver_output[2].flag==1
        #Exit due to convergence in unknown x
    else
        if false
            figure()
            ch          = solver_output[2]
            semilogy(ch.pnorms, ".-b")
            semilogy(ch.gnorms, ".-r")
            error("")
        end
        flag        = -1
#         flag        = 1
    end
    
#     println("t_mult=", t_mult)
#     println("t_ldiv=", t_ldiv)
    
    return flag, nIter, etaA
    
end

#Solve BLOCK JACOBI(A)*u=b:
function LS_solve!(LS::LS_BlockJacobi, u::GenVector{Float64}, b::GenVector{Float64};
    RelTol::Float64=0.0, AbsTol::Float64=1e-8, MaxIter::Int=100, Display::String="notify")
    
    #Optional arguments are necessary for compatibility with the calling function
    
    #Exit early if there are NaNs or Infs in b:
    if !all(isfinite.(b))
        return -1, 0, Inf   #flag, nIters, etaA
    end
        
    #Extract variables:
    F               = LS.Pl
    scaleP          = LS.scaleP
    
    #Permute vectors:
    bhat            = b[LS.p]
    uhat            = u[LS.p]
    
    #Reduce r.h.s.:
    bred            = ReduceRhs(F, bhat)
    
    #Solve:
    ldiv_slave!(uhat, F, bred)
    @mlv u[LS.p]    = uhat
    return 1, 0, 0.0   #flag, nIters, etaA
    
end

#Solve ILU0(Dominant(A))*u=b:
function LS_solve_old!(LS::LS_DominantILU0, u::GenVector{Float64}, b::GenVector{Float64};
    RelTol::Float64=0.0, AbsTol::Float64=1e-8, MaxIter::Int=100, Display::String="notify")
    
    #Optional arguments are necessary for compatibility with the calling function
    
    #Exit early if there are NaNs or Infs in b:
    if !all(isfinite.(b))
        return -1, 0, Inf   #flag, nIters, etaA
    end
        
    #Extract variables:
    F               = LS.Pl
    scaleP          = LS.scaleP
    
    #Permute vectors:
    bhat            = b[LS.p]
    uhat            = u[LS.p]
    
    #Solve ILU0(Aapprox) uhat = bhat (no iterative solver):
    ldiv_master!(uhat, F, bhat)     #Overwrites uhat
    u[LS.p]         .= uhat
    
    return 1, 0, 0.0   #flag, nIters, etaA
    
end

function LS_solve!(LS::LS_DominantILU0, u::GenVector{Float64}, b::GenVector{Float64};
    RelTol::Float64=0.0, AbsTol::Float64=1e-8, MaxIter::Int=100, Display::String="notify")
    
    #Optional arguments are necessary for compatibility with the calling function
    
    #Exit early if there are NaNs or Infs in b:
    if !all(isfinite.(b))
        return -1, 0, Inf   #flag, nIters, etaA
    end
        
    #Extract variables:
    F               = LS.Pl
    scaleP          = LS.scaleP
    
    #Permute vectors:
    bhat            = b[LS.p]
    uhat            = u[LS.p]
    
    #Solve Aapprox uhat = bhat with ILU0+GMRES:
    f               = zeros(length(uhat));
    function QN_ILU0!(x::Vector{Float64}, g::Vector{Float64}) #x=u at some iteration
        t_ini       = time()
        @avxt @. x  *= scaleP
        mul!(f, LS.APP, x)              #f = 1.0*A*x
        axpby!(-1.0, bhat, 1.0, f)
#         t_mult      += time()-t_ini
        if !isnothing(g)
            t_ini       = time()
            ldiv_master!(g, F, f)
            @avxt @. g  /= scaleP     #g ~~ Delta u. We have to return the scaled g
#             t_ldiv      += time()-t_ini
        end
        @avxt @. x      /= scaleP
        return 1
    end
    
    #Call GMRES:
    flag            = 1
    solver_output   = NLS_gmres(FW_NLS((x,g)->QN_ILU0!(x,g)), 
                                    uhat./scaleP,
                                    memory=100, MaxIter=MaxIter, 
                                    AbsTolX=AbsTol, RelTolX=AbsTol, 
                                    AbsTolG=0.0, RelTolG=0.0,
                                    Display="iter", history=true,
                                    NormFun=FW_NLS_norm((x)->norm(x)/sqrt(length(uhat))))
    @mlv uhat       = solver_output[1]*scaleP
    @mlv u[LS.p]    = uhat
    nIter           = solver_output[2].nIter 
    resv            = solver_output[2].gnorms
    etaA            = resv[length(resv)]
    if solver_output[2].flag==2
        #Exit due to convergence in residual g
    elseif solver_output[2].flag==1
        #Exit due to convergence in unknown x
    else
        if false
            figure()
            ch          = solver_output[2]
            semilogy(ch.pnorms, ".-b")
            semilogy(ch.gnorms, ".-r")
            error("")
        end
        flag        = -1
    end
    
#     println("t_mult=", t_mult)
#     println("t_ldiv=", t_ldiv)
    
    return flag, nIter, etaA
    
end

#Allocate full mass matrix, jacobian matrix and linear system matrix:
function MatricesAllocate!(solver::SolverData)

    fes             = solver.fes
    fes_dof         = fes.nDof
    nVars           = solver.nVars
    
    #(i,j) values for nnz numbers for one variable:
    imat, jmat      = ElemsDofCompute(fes, fes, iv=solver.Jm_ilocal, jv=solver.Jm_jlocal)
    
    #Construct sparse matrix and get fast-assembly index vector:
    JII, JII_pinv   = SpAlloc(Float64, imat, jmat)
    JII_nnz         = JII.colptr[JII.n+1]-1
    MII_nzval       = zeros(0)
    if uppercase(solver.JPatt)=="COMPLETE"
        MII_nzval   = solver.MII.nzval
    else
        MII         = MassMatrix(solver.Integ2D, solver.fes, iv=solver.Jm_ilocal, jv=solver.Jm_jlocal)
        MII_nzval   = MII.nzval
    end
    
    #Allocate pinv vectors:
    Jm_pinv         = Matrix{Vector{Int}}(undef, nVars, nVars)
    
    #We have the info for one block. Now loop blocks and compute info for the full matrix.
    #Note that Jm is a matrix, whereas pinv is a matrix nVars*nVars of index vectors.
    Jm_colptr       = zeros(Int, nVars*JII.n+1)
    Jm_nnz          = nVars*nVars*JII_nnz
    Jm_rowval       = zeros(Int, Jm_nnz)
    Jm_nzval        = zeros(Float64, Jm_nnz)
    Mm_nzval        = zeros(Float64, Jm_nnz)
    Jm_colptr[1]    = 1
    Jm_ss           = 0
    ssJm_ssJII      = Vector{Vector{Int}}(undef, nVars)
    alloc!(ssJm_ssJII, (JII_nnz,))
    @inbounds for JJ=1:nVars
        for jj=1:JII.n
    
            #Loop rows:
            for II=1:nVars, ss=JII.colptr[jj]:JII.colptr[jj+1]-1
                
                #Update pointer to Jacobian matrix:
                Jm_ss                   += 1
                
                #Save row:
                Jm_rowval[Jm_ss]        = (II-1)*JII.m + JII.rowval[ss]
                
                #Save values of mass matrix:
                if II==JJ
                    Mm_nzval[Jm_ss]     = MII_nzval[ss]
                end
                
                #Save vector that maps from nnz position in JII to position in Jm:
                ssJm_ssJII[II][ss]      = Jm_ss
                
            end
            
            #Save number of nnzs in current column:
            Jm_colptr[(JJ-1)*JII.n+jj+1]    = Jm_colptr[(JJ-1)*JII.n+jj] + 
                                                nVars*(JII.colptr[jj+1]-JII.colptr[jj])
        
        end
        
        #Update pinv for JJ-th block:
        for II=1:nVars
            Jm_pinv[II,JJ]      = ssJm_ssJII[II][JII_pinv]
        end

    end

    
    #Construct matrices:
    solver.Mm       = SparseMatrixCSC{Float64,Int}(nVars*JII.m, nVars*JII.n, Jm_colptr, 
                        Jm_rowval, Mm_nzval)
    solver.Jm       = SparseMatrixCSC{Float64,Int}(nVars*JII.m, nVars*JII.n, Jm_colptr, 
                        Jm_rowval, Jm_nzval)
    solver.Am       = SparseMatrixCSC{Float64,Int}(nVars*JII.m, nVars*JII.n, Jm_colptr, 
                        Jm_rowval, zeros(Float64, Jm_nnz))
    solver.Jm_pinv  = Jm_pinv
    
    return
    
end

function FluxAllocate(nVars::Int, udims::Tuple{Int,Int}, ComputeJ::Bool)

    flux_qp         = Matrix{Matrix{Float64}}(undef,nVars,2) #f[I,j] = f_Ij
    for II=1:nVars, jj=1:2
        flux_qp[II,jj]                      = zeros(udims)
    end
    dflux_du_qp     = Array{Matrix{Float64},3}(undef,nVars,2,nVars) #df[I,j,J] = df_{Ij}/du_J
    dflux_dgradu_qp = Array{Matrix{Float64},4}(undef,nVars,2,nVars,2) #df[I,j,J,k] = df_{Ij}/(du_J/dx_k)
    if ComputeJ
        for II=1:nVars, jj=1:2, JJ=1:nVars
            dflux_du_qp[II,jj,JJ]           = zeros(udims)
        end
        for II=1:nVars, jj=1:2, JJ=1:nVars, kk=1:2
            dflux_dgradu_qp[II,jj,JJ,kk]    = zeros(udims)
        end
    end
    
    return flux_qp, dflux_du_qp, dflux_dgradu_qp
    
end 

function bFluxAllocate(nVars::Int, dims::Tuple{Int,Int}, ComputeJ::Bool)
    #bFluxAllocate(nVars::Int, dims::NTuple{N,Int}, ComputeJ::Bool) where N<:Any

    flux_qp         = Vector{Matrix{Float64}}(undef,nVars) #f[I] = fn_I
    for II=1:nVars
        flux_qp[II] = zeros(dims)
    end
    dflux_du_qp     = Array{Matrix{Float64},2}(undef,nVars,nVars) #df[I,J] = dfn_{I}/du_J
    dflux_dgradu_qp = Array{Matrix{Float64},3}(undef,nVars,nVars,2) #df[I,J,k] = dfn_{I}/(du_J/dx_k) 
    if ComputeJ
        for II=1:nVars, JJ=1:nVars
            dflux_du_qp[II,JJ]          = zeros(dims)
        end
        for II=1:nVars, JJ=1:nVars, kk=1:2
            dflux_dgradu_qp[II,JJ,kk]   = zeros(dims)
        end
    end
    
    return flux_qp, dflux_du_qp, dflux_dgradu_qp
    
end

function ViscosityAllocate(nVars::Int, dims::Tuple{Int,Int}, ComputeJ::Bool)

    epsilon         = zeros(dims)
    depsilon_du     = Vector{Matrix{Float64}}(undef,nVars)
    depsilon_dgradu = Matrix{Matrix{Float64}}(undef,nVars,2)
    if ComputeJ
        alloc!(depsilon_du, dims)
        alloc!(depsilon_dgradu, dims)
    end
    return epsilon, depsilon_du, depsilon_dgradu
    
end

function ProjectFlux(flux::Matrix{MFloat}, dflux_du::Array{MFloat,3}, 
    dflux_dgradu::Array{MFloat,4}, nb::Vector{MFloat}, ComputeJ::Bool) where MFloat<:Matrix{Float64}
    
    nVars                   = size(flux, 1)
    dims                    = size(flux[1,1])
    fn, dfn_du, dfn_dgradu  = bFluxAllocate(nVars, dims, ComputeJ)
    for II=1:nVars
        @mlv    fn[II]      += flux[II,1]*nb[1]+flux[II,2]*nb[2]
    end
    if ComputeJ
        for II=1:nVars, JJ=1:nVars
            @mlv    dfn_du[II,JJ]           += dflux_du[II,1,JJ]*nb[1] + 
                                                dflux_du[II,2,JJ]*nb[2]
        end
        for II=1:nVars, JJ=1:nVars, kk=1:2
            @mlv    dfn_dgradu[II,JJ,kk]    += dflux_dgradu[II,1,JJ,kk]*nb[1] + 
                                                dflux_dgradu[II,2,JJ,kk]*nb[2]
        end
    end
    return fn, dfn_du, dfn_dgradu
    
end

#Project a flux onto a finite element space and compute derivatives:
function Rhs_Flux!(
    flux_ElemsDof   ::Vector{Matrix{Float64}},
    J_ElemsDof      ::Matrix{Matrix{Float64}},
    flux_Ik         ::Matrix{Float64}, 
    Integ2D         ::TrInt,
    f               ::Matrix{Matrix{Float64}},
    df_du           ::Array{Matrix{Float64},3},
    df_dgradu       ::Array{Matrix{Float64},4},
    gradNm_alpha    ::Vector{Matrix{Float64}},
    Nm_beta         ::Matrix{Float64},
    gradNm_beta     ::Vector{Matrix{Float64}},
    ComputeJ        ::Bool, 
    Jm_ilocal       ::Vector{Int}, 
    Jm_jlocal       ::Vector{Int} )
    
    t_ini           = time()
    
    nVars           = size(f,1)
    Jinvm           = Integ2D.Jinv
    wdetJ           = Integ2D.wdetJ

    #Compute 
    #   int_Omega f_{Ij} dpsi_alpha/dx_j dOmega
    #   = int_Omega f_{Ij} dxi_k/dx_j dpsi_alpha/dxi_k dOmega
    #   = int_Omega fhat_Ik dpsi_alpha/dxi_k dOmega
    #
    for II=1:nVars, kk=1:2
        #sum for j and assemble:
        @avxt @. flux_Ik    = ( f[II,1]*Jinvm[kk,1] + 
                                f[II,2]*Jinvm[kk,2] )*wdetJ
        BLAS.gemm!('N', 'N', 1.0, flux_Ik, gradNm_alpha[kk], 1.0, flux_ElemsDof[II])
    end
    
    #Compute 
    #   int_Omega df_{Ij}/du_J dpsi_alpha/dx_j psi_beta dOmega
    #   = int_Omega df_{Ij}/du_J dxi_k/dx_j dpsi_alpha/dxi_k psi_beta dOmega
    #   = int_Omega dfhat_{Ik}/du_J UpsilonGN_{k}_{alpha,beta} dOmega
    if ComputeJ
    
        #Change notation for pointers:
        dflux_IkJ           = flux_Ik   
        
        #Matrices UpsilonGN:
        UpsilonGN           = Vector{Matrix{Float64}}(undef,2)
        for ii=1:2
            #Upsilon[ii][:,(II.JJ) = gradNm[ii][:,II]*Nm[:,JJ]
            UpsilonGN[ii]   = UpsilonCompute(gradNm_alpha[ii], Nm_beta, iv=Jm_ilocal, jv=Jm_jlocal)    
        end
        
#         println("Upsilon GN  = ", time()-t_ini)
        
        #Sum for j and assemble:
        for JJ=1:nVars, II=1:nVars, kk=1:2
            @avxt @. dflux_IkJ  = ( df_du[II,1,JJ]*Jinvm[kk,1] + 
                                    df_du[II,2,JJ]*Jinvm[kk,2])*wdetJ
#             println("df_du nabla psi 1 = ", time()-t_ini)
            BLAS.gemm!('N', 'N', 1.0, dflux_IkJ, UpsilonGN[kk], 1.0, J_ElemsDof[II,JJ])
#             println("df_du nabla psi 2 = ", time()-t_ini)
        end
        
    end
    
    #Compute 
    #   int_Omega df_{Ij}/d(du_J/dx_l) dpsi_alpha/dx_j dpsi_beta/dx_l dOmega
    #   = int_Omega df_{Ij}/d(du_J/dx_l) dxi_k/dx_j dxi_m/dx_l dpsi_alpha/dxi_k dpsi_beta/dxi_m dOmega
    #   = int_Omega dfhat_{Ik}/d(du_J/dx_m) UpsilonGG_{k,m}_{alpha,beta} dOmega
    if ComputeJ
    
        #Change notation for pointers:
        dflux_IkJm          = flux_Ik   
        
        #Matrices UpsilonGG:
        UpsilonGG           = Matrix{Matrix{Float64}}(undef,2,2)
        for ii=1:2, jj=1:2
            #Upsilon[ii][:,(II.JJ) = gradNm[ii][:,II]*Nm[:,JJ]
            UpsilonGG[ii,jj]= UpsilonCompute(gradNm_alpha[ii], gradNm_beta[jj], iv=Jm_ilocal, jv=Jm_jlocal)
        end
        
        #Sum for j,l  and assemble:
        for JJ=1:nVars, II=1:nVars, mm=1:2, kk=1:2  
            @avxt @. dflux_IkJm = ( df_dgradu[II,1,JJ,1]*Jinvm[kk,1]*Jinvm[mm,1] + 
                                    df_dgradu[II,1,JJ,2]*Jinvm[kk,1]*Jinvm[mm,2] +
                                    df_dgradu[II,2,JJ,1]*Jinvm[kk,2]*Jinvm[mm,1] +
                                    df_dgradu[II,2,JJ,2]*Jinvm[kk,2]*Jinvm[mm,2])*wdetJ
            BLAS.gemm!('N', 'N', 1.0, dflux_IkJm, UpsilonGG[kk,mm], 1.0, J_ElemsDof[II,JJ])
        end

    end
    
    return
    
end

#Project a source term onto a finite element space and compute derivatives:
function Rhs_Source!(
    flux_ElemsDof   ::Vector{Matrix{Float64}},
    J_ElemsDof      ::Matrix{Matrix{Float64}},
    flux_I          ::Matrix{Float64}, 
    Integ2D         ::TrInt,
    Q               ::Vector{Matrix{Float64}},
    dQ_du           ::Array{Matrix{Float64},2},
    dQ_dgradu       ::Array{Matrix{Float64},3},
    Nm_alpha        ::Matrix{Float64},
    Nm_beta         ::Matrix{Float64},
    gradNm_beta     ::Vector{Matrix{Float64}},
    ComputeJ        ::Bool, 
    Jm_ilocal       ::Vector{Int}, 
    Jm_jlocal       ::Vector{Int}              )
    
    nVars           = size(Q,1)
    Jinvm           = Integ2D.Jinv
    wdetJ           = Integ2D.wdetJ

    #Compute 
    #   int_Omega Q_{I} psi_alpha dOmega
    for II=1:nVars
        @avxt @. flux_I = Q[II]*wdetJ
        BLAS.gemm!('N', 'N', 1.0, flux_I, Nm_alpha, 1.0, flux_ElemsDof[II])
    end
    
    #Compute 
    #   int_Omega dQ_I/du_J psi_alpha psi_beta dOmega
    if ComputeJ
    
        #Change notation for pointers:
        dflux_IJ                = flux_I   
        
        #Matrices UpsilonNN:
        UpsilonNN               = UpsilonCompute(Nm_alpha, Nm_beta, iv=Jm_ilocal, jv=Jm_jlocal)    #Upsilon[ii][:,(II.JJ) = gradflucNm[ii][:,II]*Nm[:,JJ]
        
        #Assemble:
        for JJ=1:nVars, II=1:nVars
            @avxt @. dflux_IJ   = dQ_du[II,JJ]*wdetJ
            BLAS.gemm!('N', 'N', 1.0, dflux_IJ, UpsilonNN, 1.0, J_ElemsDof[II,JJ])
        end
        
    end
    
    #Compute 
    #   int_Omega dQ_I/d(du_J/dx_j) psi_alpha dpsi_beta/dx_j dOmega = 
    #   int_Omega dQ_I/d(du_J/dx_j)*dxi_k/dx_j psi_alpha dpsi_beta/dxi_k
    if ComputeJ
    
        #Change notation for pointers:
        dflux_IJk               = flux_I   
        
        #Matrices UpsilonNG:
        UpsilonNG               = Vector{Matrix{Float64}}(undef,2)
        for ii=1:2
            UpsilonNG[ii]       = UpsilonCompute(Nm_alpha, gradNm_beta[ii], iv=Jm_ilocal, jv=Jm_jlocal)
        end
        
        #Assemble:
        for JJ=1:nVars, II=1:nVars, kk=1:2
            @avxt @. dflux_IJk  = (dQ_dgradu[II,JJ,1]*Jinvm[kk,1] + 
                                    dQ_dgradu[II,JJ,2]*Jinvm[kk,2])*wdetJ
            BLAS.gemm!('N', 'N', 1.0, dflux_IJk, UpsilonNG[kk], 1.0, J_ElemsDof[II,JJ])
        end
        
    end
    
    return
    
end

#Project a boundary flux onto a finite element space and compute derivatives:
function Rhs_bFlux!(
    flux_BElemsDof  ::Vector{Matrix{Float64}},
    J_BElemsDof     ::Matrix{Matrix{Float64}},
    bflux_I         ::Matrix{Float64}, 
    Binteg2D        ::TrBint,
    f               ::Vector{Matrix{Float64}},
    df_du           ::Array{Matrix{Float64},2},
    df_dgradu       ::Array{Matrix{Float64},3},
    Nm_alpha        ::Matrix{Float64},
    Nm              ::Matrix{Float64},
    gradNm          ::Vector{Matrix{Float64}},
    ComputeJ        ::Bool, 
    Jm_ilocal       ::Vector{Int}, 
    Jm_jlocal       ::Vector{Int}
    )
    
    nVars           = length(f)
    Jinvm           = Binteg2D.Jinv
    wdetJ           = Binteg2D.wdetJ
    
    #Compute 
    #   -int_Gamma f_I psi_alpha dGamma
    for II=1:nVars
        @avxt @. bflux_I    = f[II]*wdetJ
        BLAS.gemm!('N', 'N', -1.0, bflux_I, Nm_alpha, 1.0, flux_BElemsDof[II])
    end

    #Compute 
    #   -int_Gamma df_I/du_J psi_alpha psi_beta dGamma
    if ComputeJ
    
        #Change notation:
        bflux_IJ        = bflux_I
        
        #Matrix UpsilonNN:
        UpsilonNN       = UpsilonCompute(Nm_alpha, Nm, iv=Jm_ilocal, jv=Jm_jlocal)
        
        #Assemble:
        for JJ=1:nVars, II=1:nVars
            @avxt @. bflux_IJ   = df_du[II,JJ]*wdetJ
            BLAS.gemm!('N', 'N', -1.0, bflux_IJ, UpsilonNN, 1.0, J_BElemsDof[II,JJ])
        end
    
    end
    
    #Compute 
    #   -int_Gamma df_I/d(du_J/dx_k) psi_alpha dpsi_beta/dx_k dGamma = 
    #   -int_Gamma df_I/d(du_J/dx_k) dxi_l/dx_k psi_alpha dpsi_beta/dxi_l dGamma
    if ComputeJ
    
        #Change notation:
        bflux_IJl       = bflux_I
        
        #Matrix UpsilonNG:
        UpsilonNG       = Vector{Matrix{Float64}}(undef,2)
        for kk=1:2
            UpsilonNG[kk]   = UpsilonCompute(Nm_alpha, gradNm[kk], iv=Jm_ilocal, jv=Jm_jlocal)
        end
        
        #Sum for k=1,2 and assemble:
        for JJ=1:nVars, II=1:nVars, ll=1:2
            @avxt @. bflux_IJl  = ( df_dgradu[II,JJ,1]*Jinvm[ll,1] + 
                                    df_dgradu[II,JJ,2]*Jinvm[ll,2] )*wdetJ
            BLAS.gemm!('N', 'N', -1.0, bflux_IJl, UpsilonNG[ll], 1.0, J_BElemsDof[II,JJ])
        end
    
    end
    
    return
    
end

#Function to compute the r.h.s. term.
#NOTE: Here "f" is the r.h.s. vector
#NOTE: Most part of the time is spent on allocating matrices (nElems, DofPerElem^2)
#TODO: Instruction f[II] is not type-stable
function Rhs!(solver::SolverData, t::Float64, uv::Vector{Float64}, 
    ComputeJ::Bool,
    fv::GenVector{Float64}, Jm::SparseMatrixCSC{Float64,Int})
    
    t_ini       = time()

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
    FluxSource!(solver.model, _qp, ComputeJ)
#     println("Flux and source terms = ", time()-t_ini)
    
    #Check NaN's:
    for II=1:solver.nVars
        if !all(isfinite.(_qp.f[II]))
#             error("flux $(II) is not finite")
            flag    = -1
        end
        if !all(isfinite.(_qp.fB[II]))
#             error("bubble flux $(II) is not finite")
            flag    = -1
        end
        if !all(isfinite.(_qp.Q[II]))
#             error("source term $(II) is not finite")
            flag    = -1
        end
    end
    
    t_ini           = time()
    
    #Fluxes. We assume f=f(x, t, u, grad u) and fB=fB(x, t, u, grad uB),
    #with u the total (large+fine scale) solution and uB the bubble terms.
    Rhs_Flux!(flux_ElemsDof, J_ElemsDof, flux_Ik, Integ2D, 
        _qp.f, _qp.df_du, _qp.df_dgradu,
        gradNm, Nm, gradNm, ComputeJ, solver.Jm_ilocal, solver.Jm_jlocal)
    Rhs_Flux!(flux_ElemsDof, J_ElemsDof, flux_Ik, Integ2D, 
        _qp.fB, _qp.dfB_du, _qp.dfB_dgraduB, 
        gradflucNm, Nm, gradflucNm, ComputeJ, solver.Jm_ilocal, solver.Jm_jlocal)
        
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
                error("flux $(II) for boundary $(bound_id) is not finite")
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
    
    return flag, _qp.Deltat_CFL
    
end

#Check Jacobian for ib-th boundary condition:
function BC_CheckJacobian(solver::SolverData, ib::Int;
    delta::Float64=1e-6,
    Plot_df_du::Bool=false, 
    Plot_df_dgradu::Bool=false, 
    errmax::Float64=0.0)
    
    #Extract variables:
    t                   = solver.t
    uv                  = solver.uv
    fes                 = solver.fes
    fes_dof             = fes.nDof
    nVars               = solver.nVars
    u                   = GetViews(uv, nVars, fes_dof)
    Binteg2D            = solver.Binteg2D[ib]
    bmesh               = Binteg2D.bmesh
    Jinvm               = Binteg2D.Jinv
    wdetJ               = Binteg2D.wdetJ
    
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
        _bqp.df_dgradu  = bFluxAllocate(nVars, (_bqp.nElems, _bqp.nqp), true)
    
    #Evaluate flux at the boundary:
    bflux!(solver.model, solver.BC[ib].BC, _bqp, true)
    for II=1:solver.nVars, JJ=1:solver.nVars
#         println("II=$II, JJ=$JJ, |df_du|=",norm(_bqp.df_du[II,JJ],Inf))
    end
    
    #Save initial pointer:
    _bqp0               = _bqp
    
    #Perturb u and compute derivatives:
    for JJ=1:solver.nVars
        
        if Plot_df_du
            figure()
        end
        
        _bqp            = deepcopy(_bqp0)
        @mlv _bqp.u[JJ] -= delta
        zero!(_bqp.f)
        bflux!(solver.model, solver.BC[ib].BC, _bqp, false)
        f1              = deepcopy(_bqp.f)
        #
        _bqp            = deepcopy(_bqp0)
        @mlv _bqp.u[JJ] += delta
        zero!(_bqp.f)
        bflux!(solver.model, solver.BC[ib].BC, _bqp, false)
        f2              = deepcopy(_bqp.f)
        for II=1:solver.nVars
            df_du_num       = @mlv (f2[II]-f1[II])/(2*delta)
            err_num         = norm(_bqp0.df_du[II,JJ]-df_du_num,Inf)
            println("error in |dfn_{$(II)}/du_$(JJ)|=", err_num)
            if Plot_df_du && err_num>=errmax
                subplot(solver.nVars, 1, II)
                semilogy(abs.(_bqp0.df_du[II,JJ]), "xb")
                semilogy(abs.(df_du_num), "+g")
                xticks([])
                title(latexstring("\\partial fn_{", II, "}/\\partial u_{", JJ,"}"))
            end
        end
    end
    
    #Perturb gradu and compute derivatives:
    for jj=1:2, JJ=1:solver.nVars
        
        if Plot_df_dgradu
            figure()
        end
        
        _bqp                    = deepcopy(_bqp0)
        @mlv _bqp.gradu[JJ,jj]  -= delta
        zero!(_bqp.f)
        bflux!(solver.model, solver.BC[ib].BC, _bqp, false)
        f1                      = deepcopy(_bqp.f)
        #
        _bqp                    = deepcopy(_bqp0)
        @mlv _bqp.gradu[JJ,jj]  += delta
        zero!(_bqp.f)
        bflux!(solver.model, solver.BC[ib].BC, _bqp, false)
        f2                      = deepcopy(_bqp.f)
        for II=1:solver.nVars
            df_du_num           = @mlv (f2[II]-f1[II])/(2*delta)
            err_num             = norm(_bqp0.df_dgradu[II,JJ,jj]-df_du_num,Inf)
            println("error in df_{$(II)}/d(du_$(JJ)/dx_$(jj))=", err_num)
            if Plot_df_dgradu && err_num>=errmax
                subplot(solver.nVars, 1, II)
                semilogy(abs.(_bqp0.df_dgradu[II,JJ,jj]), "xb")
                semilogy(abs.(df_du_num), "+g")
                xticks([])
                title(latexstring("\\partial fn_{", II, "}/\\partial (\\partial_", jj, " u_{", JJ,"})"))
            end
        end
    end
        
    return
    
end

#Check Jacobian for volume integrals:
function CheckJacobian(solver::SolverData;
    delta::Float64=1e-6, 
    Plot_df_du::Bool=false, 
    Plot_df_dgradu::Bool=false, 
    Plot_dfB_du::Bool=false,
    Plot_dfB_dgraduB::Bool=false,
    Plot_dQ_du::Bool=false, 
    Plot_dQ_dgradu::Bool=false,
    errmax::Float64=0e-8)

    #Extract variables:
    t                   = solver.t
    uv                  = solver.uv
    fes                 = solver.fes
    fes_dof             = fes.nDof
    nVars               = solver.nVars
    u                   = GetViews(uv, nVars, fes_dof)
    Integ2D             = solver.Integ2D
    Jinvm               = Integ2D.Jinv
    wdetJ               = Integ2D.wdetJ
    
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
    _qp.f, _qp.df_du, _qp.df_dgradu         = FluxAllocate(nVars, (_qp.nElems, _qp.nqp), true)
    _qp.fB, _qp.dfB_du, _qp.dfB_dgraduB     = FluxAllocate(nVars, (_qp.nElems, _qp.nqp), true)
    _qp.Q, _qp.dQ_du, _qp.dQ_dgradu         = bFluxAllocate(nVars, (_qp.nElems, _qp.nqp), true)
    
    #Evaluate flux and source term:
    t_ini           = time()
    FluxSource!(solver.model, _qp, true)
    
    #Save initial pointer:
    _qp0                = _qp
    
    #Perturb u and compute derivatives of f:
    println("\nDerivatives of f w.r.t. u")
    for JJ=1:solver.nVars
        
        if Plot_df_du
            figure()
        end
        
        _qp             = deepcopy(_qp0)
        @mlv _qp.u[JJ] -= delta
        zero!(_qp.f)
        FluxSource!(solver.model, _qp, false)
        f1              = deepcopy(_qp.f)
        #
        _qp            = deepcopy(_qp0)
        @mlv _qp.u[JJ] += delta
        zero!(_qp.f)
        FluxSource!(solver.model, _qp, false)
        f2              = deepcopy(_qp.f)
        for ii=1:2, II=1:solver.nVars
            df_du_num       = @mlv (f2[II,ii]-f1[II,ii])/(2*delta)
            err_num         = norm(_qp0.df_du[II,ii,JJ]-df_du_num,Inf)
            println("error in |df_{$(II)$(ii)}/du_$(JJ)|=", err_num)
            if Plot_df_du && err_num>=errmax
                subplot(solver.nVars, 2, (II-1)*2+ii)
                semilogy(abs.(_qp0.df_du[II,ii,JJ]), "xb")
                semilogy(abs.(df_du_num), "+g")
                xticks([])
                title(latexstring("\\partial f_{", II, ii, "}/\\partial u_{", JJ,"}"))
            end
        end
    end
    
    #Perturb gradu and compute derivatives of f:
    println("\nDerivatives of f w.r.t. du/dx")
    for jj=1:2, JJ=1:solver.nVars
        
        if Plot_df_dgradu
            figure()
        end
        
        _qp                     = deepcopy(_qp0)
        @mlv _qp.gradu[JJ,jj]   -= delta
        zero!(_qp.f)
        FluxSource!(solver.model, _qp, false)
        f1                      = deepcopy(_qp.f)
        #
        _qp                     = deepcopy(_qp0)
        @mlv _qp.gradu[JJ,jj]   += delta
        zero!(_qp.f)
        FluxSource!(solver.model, _qp, false)
        f2                      = deepcopy(_qp.f)
        for ii=1:2, II=1:solver.nVars
            df_du_num       = @mlv (f2[II,ii]-f1[II,ii])/(2*delta)
            err_num         = norm(_qp0.df_dgradu[II,ii,JJ,jj]-df_du_num,Inf)
            println("error in |df_{$(II)$(ii)}/d(du_$(JJ)/dx_$(jj))|=", err_num)
            if Plot_df_dgradu && err_num>=errmax
                subplot(solver.nVars, 2, (II-1)*2+ii)
                semilogy(abs.(_qp0.df_dgradu[II,ii,JJ,jj]), "xb")
                semilogy(abs.(df_du_num), "+g")
                xticks([])
                title(latexstring("\\partial f_{", II, ii, "}/\\partial (\\partial_", jj, " u_{", JJ,"})"))
            end
        end
    end
    
    #Perturb u and compute derivatives of f':
    println("\nDerivatives of f' w.r.t. u")
    for JJ=1:solver.nVars
        
        if Plot_dfB_du
            figure()
        end
        
        _qp             = deepcopy(_qp0)
        @mlv _qp.u[JJ] -= delta
        zero!(_qp.fB)
        FluxSource!(solver.model, _qp, false)
        f1              = deepcopy(_qp.fB)
        #
        _qp            = deepcopy(_qp0)
        @mlv _qp.u[JJ] += delta
        zero!(_qp.fB)
        FluxSource!(solver.model, _qp, false)
        f2              = deepcopy(_qp.fB)
        for ii=1:2, II=1:solver.nVars
            df_du_num       = @mlv (f2[II,ii]-f1[II,ii])/(2*delta)
            err_num         = norm(_qp0.dfB_du[II,ii,JJ]-df_du_num,Inf)
            println("error in |df'_{$(II)$(ii)}/du_$(JJ)|=", err_num)
            if Plot_dfB_du && err_num>=errmax
                subplot(solver.nVars, 2, (II-1)*2+ii)
                semilogy(abs.(_qp0.dfB_du[II,ii,JJ]), "xb")
                semilogy(abs.(df_du_num), "+g")
                xticks([])
                title(latexstring("\\partial f'_{", II, ii, "}/\\partial u_{", JJ,"}"))
            end
        end
    end
    
    #Perturb graduB and compute derivatives of fB:
    println("\nDerivatives of f' w.r.t. du'/dx")
    for jj=1:2, JJ=1:solver.nVars
        
        if Plot_dfB_dgraduB
            figure()
        end
        
        _qp                     = deepcopy(_qp0)
        @mlv _qp.graduB[JJ,jj]  -= delta
        zero!(_qp.fB)
        FluxSource!(solver.model, _qp, false)
        f1                      = deepcopy(_qp.fB)
        #
        _qp                     = deepcopy(_qp0)
        @mlv _qp.graduB[JJ,jj]  += delta
        zero!(_qp.fB)
        FluxSource!(solver.model, _qp, false)
        f2                      = deepcopy(_qp.fB)
        for ii=1:2, II=1:solver.nVars
            df_du_num       = @mlv (f2[II,ii]-f1[II,ii])/(2*delta)
            err_num         = norm(_qp0.dfB_dgraduB[II,ii,JJ,jj]-df_du_num,Inf)
            println("error in |df'_{$(II)$(ii)}/d(du'_$(JJ)/dx_$(jj))|=", err_num)
            if Plot_dfB_dgraduB && err_num>=errmax
                subplot(solver.nVars, 2, (II-1)*2+ii)
                semilogy(abs.(_qp0.dfB_dgraduB[II,ii,JJ,jj]), "xb")
                semilogy(abs.(df_du_num), "+g")
                xticks([])
                title(latexstring("\\partial f'_{", II, ii, "}/\\partial (\\partial_", jj, " u'_{", JJ,"})"))
            end
        end
    end
    
    #Perturb u and compute derivatives of Q:
    println("\nDerivatives of Q w.r.t. u")
    for JJ=1:solver.nVars
        
        if Plot_dQ_du
            figure()
        end
        
        _qp             = deepcopy(_qp0)
        @mlv _qp.u[JJ] -= delta
        zero!(_qp.Q)
        FluxSource!(solver.model, _qp, false)
        f1              = deepcopy(_qp.Q)
        #
        _qp            = deepcopy(_qp0)
        @mlv _qp.u[JJ] += delta
        zero!(_qp.Q)
        FluxSource!(solver.model, _qp, false)
        f2              = deepcopy(_qp.Q)
        for II=1:solver.nVars
            dQ_du_num       = @mlv (f2[II]-f1[II])/(2*delta)
            err_num         = norm(_qp0.dQ_du[II,JJ]-dQ_du_num,Inf)
            println("error in |dQ_{$(II)}/du_$(JJ)|=", err_num)
            if Plot_dQ_du && err_num>=errmax
                subplot(solver.nVars, 1, II)
                semilogy(abs.(_qp0.dQ_du[II,JJ]), "xb")
                semilogy(abs.(dQ_du_num), "+g")
                xticks([])
                title(latexstring("\\partial Q_{", II, "}/\\partial u_{", JJ,"}"))
            end
        end
    end
    
    #Perturb gradu and compute derivatives of f:
    println("\nDerivatives of Q w.r.t. du/dx")
    for jj=1:2, JJ=1:solver.nVars
        
        if Plot_dQ_dgradu
            figure()
        end
        
        _qp                     = deepcopy(_qp0)
        @mlv _qp.gradu[JJ,jj]   -= delta
        zero!(_qp.Q)
        FluxSource!(solver.model, _qp, false)
        f1                      = deepcopy(_qp.Q)
        #
        _qp                     = deepcopy(_qp0)
        @mlv _qp.gradu[JJ,jj]   += delta
        zero!(_qp.Q)
        FluxSource!(solver.model, _qp, false)
        f2                      = deepcopy(_qp.Q)
        for II=1:solver.nVars
            dQ_du_num       = @mlv (f2[II]-f1[II])/(2*delta)
            err_num         = norm(_qp0.dQ_dgradu[II,JJ,jj]-dQ_du_num,Inf)
            println("error in |dQ_{$(II)}/d(du_$(JJ)/dx_$(jj))|=", err_num)
            if Plot_dQ_dgradu && err_num>=errmax
                subplot(solver.nVars, 1, II)
                semilogy(abs.(_qp0.dQ_dgradu[II,JJ,jj]), "xb")
                semilogy(abs.(dQ_du_num), "+g")
                xticks([])
                title(latexstring("\\partial Q_{", II, "}/\\partial (\\partial_", jj, " u_{", JJ,"})"))
            end
        end
    end
    
    return
    
end

function PlotLSPermutations(LS::LS_SCILU0_GMRES, MII::SparseMatrixCSC{Float64,Int}, nVars::Int)
    
    nMasters    = LS.Pl.nMasters
    nSlaves     = LS.Pl.nSlaves
    nMasters_II = Int(nMasters/nVars)
    nSlaves_II  = Int(nSlaves/nVars)
    
    nDof_II                                         = size(MII,1)
    pinv_II                                         = zeros(Int, nDof_II)
    pinv_II[nMasters_II+1:nMasters_II+nSlaves_II]   = 1:nSlaves_II
    pinv_II[1:nMasters_II]                          = nSlaves_II .+ 
                                                        invperm(symrcm(MII[1:nMasters_II,
                                                                        1:nMasters_II]))
                                                                        
    pinv        = zeros(Int, nVars*nDof_II)
    for II=1:nVars, iDof=1:nDof_II
        #Degree iDof of variable II goes to position (pinv[iDof]-1)*nVars + II:
        pinv[(II-1)*nDof_II+iDof]   = (pinv_II[iDof]-1)*nVars + II
    end
    p           = invperm(pinv)
    A1          = LS.A[p,p]
    
    fig     = PyPlotSubPlots(1, 2, w=15.0, h=12.0, left=0.9, right=0.4, bottom=1.1, top=1.0)
    subplot(1,2,1)
    spy(LS.A)
    subplot(1,2,2)
    spy(LS.APP)
    plot([nSlaves+1-0.5, nSlaves+nMasters], [nSlaves+1-0.5, nSlaves+1-0.5], "--b")
    plot([nSlaves+1-0.5, nSlaves+1-0.5], [nSlaves+1-0.5, nSlaves+nMasters], "--b")
    xlim(0, nSlaves+nMasters)
    ylim(nSlaves+nMasters, 0)
    
    error("")
    
end

#Update L2L2 error:
function errL2L2_(t_n::Float64, errL2L2::Float64, t_np1::Float64, e_np1::Float64)

    errL2L2     = sqrt((errL2L2^2*t_n + e_np1^2*(t_np1-t_n))/t_np1)
    return errL2L2
    
end

#Compute mass*u product. Note that solver.Mm may not be the mass matrix, put only the 
#Jacobi part:
function MuProduct(MII::SparseMatrixCSC{Float64,Int}, uv::Vector{Float64}, nVars::Int)

    bv          = zeros(size(uv))
    b_views     = GetViews(bv, nVars, size(MII,1))
    u_views     = GetViews(uv, nVars, size(MII,1))
    
    for II=1:nVars
        b_views[II]     .= MII*u_views[II]
    end
    
    return bv
    
end
