#Create structure for given mesh with given integration order:
function TrInt(mesh::TrMesh, order::Int; N::Int=1)

    Integ           = TrInt()
    Integ.mesh      = mesh
    Integ.order     = order
    QRule           = WVRule(order, N)
    Integ.QRule     = QRule
    
    #Compute quadrature nodes:
    x_qp            = Vector{Matrix{Float64}}(undef,2)
    for ii=1:length(x_qp)
        x_qp[ii]    = zeros(mesh.nElems, QRule.nqp)
    end
    xm              = view(mesh.ElemsCoords,:, 1:2:2*mesh.NodesPerElem-1)
    ym              = view(mesh.ElemsCoords,:, 2:2:2*mesh.NodesPerElem)
    Nm              = NCompute(mesh, QRule.xi)
    BLAS.gemm!('N', 'T', 1.0, xm, Nm, 0.0, x_qp[1]) 
    BLAS.gemm!('N', 'T', 1.0, ym, Nm, 0.0, x_qp[2]) 
    Integ.x         = x_qp
    
    #Compute Jacobian inverse and detJ:
    gradNm          = gradNCompute(mesh, QRule.xi)
    Jinv11          = zeros(mesh.nElems, QRule.nqp)
    Jinv21          = zeros(mesh.nElems, QRule.nqp)
    Jinv12          = zeros(mesh.nElems, QRule.nqp)
    Jinv22          = zeros(mesh.nElems, QRule.nqp)
    BLAS.gemm!('N', 'T', 1.0, ym, gradNm[2], 0.0, Jinv11)    # gemm!(Jinv21, ym, transpose(view(gradNm[1], :, 1:mesh.NodesPerElem) ), -1.0, 0.0)
    BLAS.gemm!('N', 'T', -1.0, ym, gradNm[1], 0.0, Jinv21)   # gemm!(Jinv21, ym, transpose(view(gradNm[1], :, 1:mesh.NodesPerElem) ), -1.0, 0.0)
    BLAS.gemm!('N', 'T', -1.0, xm, gradNm[2], 0.0, Jinv12)   # gemm!(Jinv21, ym, transpose(view(gradNm[1], :, 1:mesh.NodesPerElem) ), -1.0, 0.0)
    BLAS.gemm!('N', 'T', 1.0, xm, gradNm[1], 0.0, Jinv22)    # gemm!(Jinv21, ym, transpose(view(gradNm[1], :, 1:mesh.NodesPerElem) ), -1.0, 0.0)
    wdetJ           = @mlv    Jinv11*Jinv22 - Jinv21*Jinv12  #this is actually detJ
    @mlv    Jinv11  /= wdetJ
    @mlv    Jinv21  /= wdetJ
    @mlv    Jinv12  /= wdetJ
    @mlv    Jinv22  /= wdetJ
    
    #Multiply detJ by weights of quadrature rule:
    wv                      = QRule.w
    @avx for iqp=1:QRule.nqp, iElem=1:mesh.nElems
        wdetJ[iElem, iqp]   *= wv[iqp]
    end
    
    #Save pointers:
    Integ.Jinv11    = Jinv11
    Integ.Jinv21    = Jinv21
    Integ.Jinv12    = Jinv12
    Integ.Jinv22    = Jinv22
    Integ.Jinv      = Matrix{Matrix{Float64}}(undef,2,2)
    Integ.Jinv[1,1] = Jinv11
    Integ.Jinv[2,1] = Jinv21
    Integ.Jinv[1,2] = Jinv12
    Integ.Jinv[2,2] = Jinv22
    Integ.wdetJ     = wdetJ
    
    return Integ
    
end

#Create structure for given mesh with given integration order:
function TrBint(mesh::TrMesh, ib::Int, order::Int)
    
    Integ           = TrBint()
    Integ.mesh      = mesh
    bmesh           = mesh.bmesh[ib]
    Integ.bmesh     = bmesh
    Integ.order     = order
    QRule           = GLRule(order)
    Integ.QRule     = QRule
    wv              = QRule.w
    Integ.xi2D      = TrBint_xi_trace(QRule.xi)
    
    #Quadrature nodes:
    x_qp            = Vector{Matrix{Float64}}(undef,2)
    for ii=1:length(x_qp)
        x_qp[ii]    = zeros(bmesh.nElems, QRule.nqp)
    end
    xm              = view(bmesh.ElemsCoords,:, 1:2:2*bmesh.NodesPerElem-1)
    ym              = view(bmesh.ElemsCoords,:, 2:2:2*bmesh.NodesPerElem)
    Nm              = NCompute(bmesh, QRule.xi)
    BLAS.gemm!('N', 'T', 1.0, xm, Nm, 0.0, x_qp[1])   
    BLAS.gemm!('N', 'T', 1.0, ym, Nm, 0.0, x_qp[2])   
    Integ.x         = x_qp
    
    #Compute tangent vector, normal vector and Jacobian:
    gradNm          = gradNCompute(bmesh, QRule.xi)
    x1m             = view(bmesh.ElemsCoords,:,1:2:2*bmesh.NodesPerElem-1)
    x2m             = view(bmesh.ElemsCoords,:,2:2:2*bmesh.NodesPerElem)
    dt1_ds          = x1m*transpose(gradNm[1])
    dt2_ds          = x2m*transpose(gradNm[1])
    wdetJ           = @mlv    sqrt(dt1_ds*dt1_ds+dt2_ds*dt2_ds)
    t_qp            = [ dt1_ds./wdetJ, dt2_ds./wdetJ ]
    n_qp            = [ copy(t_qp[2]), -t_qp[1] ]
    @avx for iqp=1:QRule.nqp, iElem=1:bmesh.nElems
        wdetJ[iElem, iqp]   *= wv[iqp]
    end
    Integ.Jinv      = JinvCompute(Integ)
    
    #Save pointers:
    Integ.wdetJ     = wdetJ
    Integ.t_qp      = t_qp
    Integ.n_qp      = n_qp

    return Integ
    
end

#Compute position of quadrature nodes:
function QuadNodesCoords(Integ2D::TrInt)

#     error("Deprecated. Use Integ2D.x")
    return Integ2D.x
    
    #Allocate:
    x_qp            = Vector{Matrix{Float64}}(undef,2)
    for ii=1:length(x_qp)
        x_qp[ii]    = zeros(Integ2D.mesh.nElems, Integ2D.QRule.nqp)
    end
    
    #Apply isoparametric transformation:
    mesh            = Integ2D.mesh
    xm              = view(mesh.ElemsCoords,:, 1:2:2*mesh.NodesPerElem-1)
    ym              = view(mesh.ElemsCoords,:, 2:2:2*mesh.NodesPerElem)
    Nm              = NCompute(mesh, Integ2D.QRule.xi)
    BLAS.gemm!('N', 'T', 1.0, xm, Nm, 0.0, x_qp[1]) 
    BLAS.gemm!('N', 'T', 1.0, ym, Nm, 0.0, x_qp[2]) 

    return x_qp
    
end
function QuadNodesCoords(Integ2D::TrBint)

#     error("Deprecated. Use Integ2D.x")
    return Integ2D.x
    
    bmesh           = Integ2D.bmesh
    
    #Allocate:
    x_qp            = Vector{Matrix{Float64}}(undef,2)
    for ii=1:length(x_qp)
        x_qp[ii]    = zeros(bmesh.nElems, Integ2D.QRule.nqp)
    end
    
    #Apply isoparametric transformation:
    xm              = view(bmesh.ElemsCoords,:, 1:2:2*bmesh.NodesPerElem-1)
    ym              = view(bmesh.ElemsCoords,:, 2:2:2*bmesh.NodesPerElem)
    Nm              = NCompute(bmesh, Integ2D.QRule.xi)
    BLAS.gemm!('N', 'T', 1.0, xm, Nm, 0.0, x_qp[1])   
    BLAS.gemm!('N', 'T', 1.0, ym, Nm, 0.0, x_qp[2])   
    
    return x_qp
    
end

#For a given vector "xi" in a boundary (1D), compute corresponding values of "xi" in the triangle,
#for each possible face:
function TrBint_xi_trace(xiv::AbstractVector{Float64})

    #Coordinates of nodes (1,2,3,1) in the reference triangle:
    xi1         = [ -sqrt(3)/2, sqrt(3)/2,  0.0, -sqrt(3)/2 ]
    xi2         = [ -1/2,       -1/2,       1.0, -1/2 ]
    
    #Lagrangian shape functions associated with starting and ending nodes of the face:
    N_start     = @mlv    0.5*(1.0-xiv)
    N_end       = @mlv    0.5*(1.0+xiv)
    
    #Create xim for each face:
    xim         = Vector{Matrix{Float64}}(undef,3)
    nqp         = length(xiv)
    for face=1:3
        xim[face]       = zeros(nqp,2)
        xim[face][:,1]  = xi1[face]*N_start + xi1[face+1]*N_end
        xim[face][:,2]  = xi2[face]*N_start + xi2[face+1]*N_end
    end
    
    return xim
    
end

#Return vector with elements and natural coordinates in the triangle 
#for the quadrature nodes in a boundary integral:
function QuadNodes_trace(Integ2D::TrBint)

    bmesh       = Integ2D.bmesh
    QRule       = Integ2D.QRule
    xiv         = QRule.xi    
    
    #Compute matrices:
    elemv       = bmesh.ParentElems*ones(Int,1,QRule.nqp)
    xi2D        = TrBint_xi_trace(xiv)
    xim         = zeros(bmesh.nElems, QRule.nqp, 2)
    @inbounds for iElem=1:bmesh.nElems
        xim[iElem,:,:]      .= xi2D[bmesh.ParentFaces[iElem]]
    end
    
    return elemv, xim
    
end

#Compute solution at the quadrature nodes of the boundary:
function SolutionCompute(Binteg::TrBint, uv::AbstractVector{Float64}, FES::FES_Union)
    
    return SolutionCompute(Binteg, [uv], FES)[1]
    
end
function SolutionGradCompute(Binteg::TrBint, uv::Vector{Float64}, FES::FES_Union)

    Du_ElemsQp  = SolutionGradCompute(Binteg, [uv], FES)
    
    return [ Du_ElemsQp[1,1], Du_ElemsQp[1,2] ]
    
end
function SolutionCompute(Binteg::TrBint, uv::Vector{<:AbstractVector{Float64}}, FES::FES_Union) 

    bmesh       = Binteg.bmesh
    nVars       = length(uv)
    
    #Matrices <ngp, DofPerElem> with shape functions:
    Nm          = NCompute(FES, Binteg.xi2D)
    
    #Allocate matrices u <nElems, nqp> with values of solution at each 
    #quadrature node:
    u_ElemsQp   = Vector{Matrix{Float64}}(undef, length(uv))
    for iVar=1:nVars
        u_ElemsQp[iVar] = zeros(bmesh.nElems, Binteg.QRule.nqp)
    end
    
    #Loop faces:
    for face=1:3
        
        #Select elements:
        face_belems     = bmesh.ParentFaces.==face
        face_elems      = bmesh.ParentElems[face_belems]
        
        for iVar=1:nVars
            
            #Matrix <nElems,DofPerElem> with values for each element and dof:
            u_ElemsDof                      = uv[iVar][FES.ElemsDof[face_elems,:]]
    
            #Compute umat*Nm^T, which gives a matrix <nElems, ngp>:
            u_ElemsQp[iVar][face_belems,:]  .= u_ElemsDof*transpose(Nm[face])
            
        end
        
    end
    
    return u_ElemsQp #Vector{Matrix{Float64}}
    
end
function SolutionGradCompute(Binteg::TrBint, uv::Vector{<:AbstractVector{Float64}}, FES::FES_Union)

    bmesh       = Binteg.bmesh
    nVars       = length(uv)
    
    #Matrices <ngp, DofPerElem> with shape functions:
    gradNm      = gradNCompute(FES, Binteg.xi2D)
    
    #Allocate output:
    Du_ElemsQp  = Matrix{Matrix{Float64}}(undef,nVars,2)
    for iVar=1:nVars, jj=1:2
        Du_ElemsQp[iVar,jj] = zeros(bmesh.nElems, Binteg.QRule.nqp)
    end 
    
    #Loop faces:
    for face=1:3
        
        #Select elements:
        face_belems     = bmesh.ParentFaces.==face
        face_elems      = bmesh.ParentElems[face_belems]
        
        for iVar=1:nVars
            #Matrix <nElems,DofPerElem> with values for each element and dof:
            u_ElemsDof      = uv[iVar][FES.ElemsDof[face_elems,:]]
        
            #Compute umat*Nm^T, which gives a matrix <nElems, ngp>:
            Du_ElemsQp[iVar,1][face_belems,:]    .= u_ElemsDof*transpose(gradNm[1][face])
            Du_ElemsQp[iVar,2][face_belems,:]    .= u_ElemsDof*transpose(gradNm[2][face])
        end
        
    end
    
    return Du_ElemsQp
    
end
function Solution_gradCompute(Binteg::TrBint, uv::Vector{<:AbstractVector{Float64}}, FES::FES_Union)

    Du      = SolutionGradCompute(Binteg, uv, FES)
    du      = Solution_gradCompute(Du, Binteg.Jinv)
    return du
    
end     

#Compute Jinv at the quadrature nodes of the boundary:
function JinvCompute(Binteg::TrBint)

    mesh            = Binteg.mesh
    bmesh           = Binteg.bmesh
    QRule           = Binteg.QRule
    
    #Load gradients:
    gradNm          = gradNCompute(mesh, Binteg.xi2D)
    
    #Allocate output:
    Jinv11          = zeros(bmesh.nElems, QRule.nqp)
    Jinv21          = zeros(bmesh.nElems, QRule.nqp)
    Jinv12          = zeros(bmesh.nElems, QRule.nqp)
    Jinv22          = zeros(bmesh.nElems, QRule.nqp)
    
    #Loop faces and compute Jinv:
    for face=1:3
    
        face_belems     = bmesh.ParentFaces.==face
        face_elems      = bmesh.ParentElems[face_belems]
        
        #Mesh nodes:
        x1              = mesh.ElemsCoords[face_elems, 1:2:2*mesh.NodesPerElem-1]
        x2              = mesh.ElemsCoords[face_elems, 2:2:2*mesh.NodesPerElem]
        
        #Compute Jacobian inverses:
        Jinv11[face_belems,:]   = x2*transpose(gradNm[2][face])
        Jinv21[face_belems,:]   = -x2*transpose(gradNm[1][face])
        Jinv12[face_belems,:]   = -x1*transpose(gradNm[2][face])
        Jinv22[face_belems,:]   = x1*transpose(gradNm[1][face])
        
    end
    detJ            = @mlv    Jinv11*Jinv22 - Jinv21*Jinv12 
    @mlv    Jinv11  /= detJ
    @mlv    Jinv21  /= detJ
    @mlv    Jinv12  /= detJ
    @mlv    Jinv22  /= detJ

    #Pack Jinv:
    Jinv            = Matrix{Matrix{Float64}}(undef,2,2)
    Jinv[1,1]       = Jinv11
    Jinv[2,1]       = Jinv21
    Jinv[1,2]       = Jinv12
    Jinv[2,2]       = Jinv22
    
    return Jinv
    
end

#----------------------------------------------------------
#INTEGRALS THAT MAY BE USEFUL AS EXAMPLES:

#Compute mass matrix for given FES: PSpace+BSpace:
function MassMatrix(Integ2D::TrInt, FES::TrFES)

    return MassMatrix(Integ2D, FES, FES)
    
end
function MassMatrix(Integ2D::TrInt, fes1::TrFES, fes2::TrFES)

    nElems          = Integ2D.mesh.nElems
    nqp             = Integ2D.QRule.nqp
    xim             = Integ2D.QRule.xi
    
    #Create matrices Upsilon_IJ=N_I*N_J:
    N1m             = NCompute(fes1, xim)
    N2m             = NCompute(fes2, xim)
    Upsilon         = UpsilonCompute(N1m, N2m)
    
    #Create matrix with nnz elements:
    smat            = Integ2D.wdetJ*Upsilon
    
    #Indices:
    imat,jmat       = ElemsDofCompute(fes1, fes2)
    
    #Construct sparse matrix:
    Mm              = sparse(imat[:], jmat[:], smat[:])

    return Mm
    
end

#Instead of returning the mass matrix M, return matrices Mhat, P1, P2 
#such that M = P1^T Mhat P2, where Mhat is block diagonal, and 
#P1 and P2 are the prolongation matrices from fes1 and fes2 to 
#the corresponding discontinuous spaces with identical basis functions.
function ProlMassMatrix(Integ2D::TrInt, fes1::TrFES, fes2::TrFES)

    nElems          = Integ2D.mesh.nElems
    nqp             = Integ2D.QRule.nqp
    xim             = Integ2D.QRule.xi
    
    #Create matrices Upsilon_IJ=N_I*N_J:
    N1m             = NCompute(fes1, xim)
    N2m             = NCompute(fes2, xim)
    Upsilon         = UpsilonCompute(N1m, N2m)
    
    #Create matrix with nnz elements:
    Mhat_smat       = Integ2D.wdetJ*Upsilon
    
    #Indices for Mhat:
    Mhat_imat       = zeros(Int, nElems, fes1.DofPerElem*fes2.DofPerElem)
    Mhat_jmat       = zeros(Int, nElems, fes1.DofPerElem*fes2.DofPerElem)
    @inbounds for iElem=1:nElems, iDof=1:fes1.DofPerElem, jDof=1:fes2.DofPerElem
        kDof                    = iDof+(jDof-1)*fes1.DofPerElem
        Mhat_imat[iElem, kDof]  = iDof+(iElem-1)*fes1.DofPerElem
        Mhat_jmat[iElem, kDof]  = jDof+(iElem-1)*fes2.DofPerElem 
    end 
    
    #Construct sparse matrix:
    Mm              = sparse(view(Mhat_imat,:), view(Mhat_jmat,:), 
                        view(Mhat_smat,:))
    
    return Mm, ProlongationMatrix(fes1), ProlongationMatrix(fes2)
    
end

#Project function for given FES:
function L2Projection(Integ2D::TrInt, 
    fun::FunctionWrapper{Matrix{Float64}, Tuple{Vector{Matrix{Float64}}}}, 
    FES::TrFES)

    nElems          = Integ2D.mesh.nElems
    nqp             = Integ2D.QRule.nqp
    xim             = Integ2D.QRule.xi
    
    #Compute function*|J|:
    x_qp            = QuadNodesCoords(Integ2D)
    fun_qp          = fun(x_qp)
    omega_qp        = @mlv fun_qp*Integ2D.wdetJ
    
    #Load shape functions:
    Nm              = NCompute(FES, xim)
    
    #Create matrices with nnz elements:
    sm              = omega_qp*Nm
    
    #Assemble:
    bv              = VectorAssemble(FES.ElemsDof[:], sm[:])

    return bv
    
end

#Compute stiffness matrix for given FES:
function StiffnessMatrix(Integ2D::TrInt, FES::TrFES)

    nElems          = Integ2D.mesh.nElems
    nqp             = Integ2D.QRule.nqp
    xim             = Integ2D.QRule.xi
    Jinv            = Integ2D.Jinv
    wdetJ           = Integ2D.wdetJ
    
    #Create matrices Upsilon_IJ=N_I*N_J:
    gradNm          = gradNCompute(FES, xim)
    UpsilonGG       = Matrix{Matrix{Float64}}(undef,2,2) #One UpsilonGG_(I,J)_(i,j)=dN_i/dxi_I*dN_j/dxi_J
    for JJ=1:2, II=1:2
        UpsilonGG[II,JJ]    = UpsilonCompute(gradNm[II], gradNm[JJ])
    end
    
    #We have to compute 
    #   K_ij    = int( dphi_i/dx_k * dphi_j/dx_k * dOmega )
    #           = int( dphi_i/dxi_K * dphi_j/dxi_L * dxi_K/dx_k * dxi_L/dx_k * dOmega ) 
    #           = int ( (dxi_K/dx_k*dxi_L/dx_k) * Upsilon_(KL)_(ij) )
    thetam          = zeros(nElems, nqp)
    smat            = zeros(nElems, FES.DofPerElem*FES.DofPerElem)
    for LL=1:2, KK=1:2
        @mlv    thetam  = (Jinv[KK,1]*Jinv[LL,1]+Jinv[KK,2]*Jinv[LL,2])*wdetJ
        BLAS.gemm!('N', 'N', 1.0, thetam, UpsilonGG[KK,LL], 1.0, smat)
    end
    
    #Indices:
    imat,jmat       = ElemsDofCompute(FES, FES)
    
    #Construct sparse matrix:
    Km              = sparse(imat[:], jmat[:], smat[:])
    
    return Km
    
end

#Impose flux through boundary:
function BoundaryFlux(Integ2D::TrBint, 
    fun::FunctionWrapper{Matrix{Float64}, Tuple{Vector{Matrix{Float64}}}}, 
    FES::TrFES)
   
    bmesh           = Integ2D.bmesh
    nElems          = bmesh.nElems
    nqp             = Integ2D.QRule.nqp
    xiv             = Integ2D.QRule.xi
    wdetJ           = Integ2D.wdetJ
    
    #Quadrature nodes:
    x_qp            = QuadNodesCoords(Integ2D)
    
    #Impose flux at quadrature nodes:
    flux_qp         = fun(x_qp).*wdetJ
    
    #Shape functions at the trace of the boundary:
    Nm              = NCompute(FES, Integ2D.xi2D)
    
    #Compute -flux*N_i:
    smat            = mult(Integ2D, -1.0, flux_qp, Nm)
    imat            = FES.ElemsDof[bmesh.ParentElems,:]
    
    return VectorAssemble(imat[:], smat[:], N=FES.nDof)
    
end

#Jacobian and r.h.s. of boundary term associated with Dirichlet boundary condition.
#Since this term is 
#   - int( f_j n_h psi_I ds ) 
#   - int( - epsilon * du/dx_j n_j psi_I ds) - int(sigma * (u-uDir) psi_I ds)
#We can write it as A*u_J + b, with
#   A_IJ    = int(epsilon*dpsi_J/dx_j n_j psi_I ds) - int(sigma psi_J psi_I ds) 
#           = int(epsilon * Jinv_kj * n_j * dpsi_J/dxi_k * psi_I * ds)
#               - int(sigma psi_J psi_I ds) 
#   b_I     = int(sigma uDir psi_I ds)
function DirichletCondition(Binteg::TrBint, fes::TrFES, 
    uDir::FunctionWrapper{Matrix{Float64}, Tuple{Vector{Matrix{Float64}}}}, 
    epsilon::Float64, sigma::Float64)
    
    bmesh           = Binteg.bmesh
    nElems          = bmesh.nElems
    nqp             = Binteg.QRule.nqp
    xiv             = Binteg.QRule.xi
    wdetJ           = Binteg.wdetJ
    n_qp            = Binteg.n_qp
    
    #Shape functions at the trace of the boundary:
    Nm              = NCompute(fes, Binteg.xi2D)
    gradNm          = gradNCompute(fes, Binteg.xi2D)
    UpsilonNN       = UpsilonCompute(Nm, Nm)
    UpsilonNG       = Vector{Vector{Matrix{Float64}}}(undef,2)    #UpsilonNG[kk][:,(II,JJ)] contains Nm[:,II]*gradNm[kk][:,JJ]
    for kk=1:2
        UpsilonNG[kk]   = UpsilonCompute(Nm, gradNm[kk])
    end
    
    #Jacobian inverse at quadrature nodes:
    Jinv            = JinvCompute(Binteg)
    
    #Flux at quadrature nodes (term epsilon, k=1)
    flux_qp         = @mlv    epsilon * (Jinv[1,1]*n_qp[1] + Jinv[1,2]*n_qp[2]) * wdetJ
    smat            = mult(Binteg, 1.0, flux_qp, UpsilonNG[1])
    #Flux at quadrature nodes (term epsilon, k=2)
    flux_qp         = @mlv    epsilon * (Jinv[2,1]*n_qp[1] + Jinv[2,2]*n_qp[2]) * wdetJ
    mult!(Binteg, 1.0, flux_qp, UpsilonNG[2], 1.0, smat)
    #Flux at quadrature nodes (penalty term):
    flux_qp         = @mlv    (- sigma) * wdetJ
    mult!(Binteg, 1.0, flux_qp, UpsilonNN, 1.0, smat)
    #Assemble:
    ElemsDof        = fes.ElemsDof[bmesh.ParentElems,:]
    imat, jmat      = ElemsDofCompute(fes, fes, elems=bmesh.ParentElems)
    Am              = sparse(imat[:], jmat[:], smat[:], fes.nDof, fes.nDof)
    
    #Rhs vector:
    uDir_qp         = uDir(QuadNodesCoords(Binteg))
    flux_qp         = @mlv    sigma*uDir_qp*wdetJ
    smat            = mult(Binteg, 1.0, flux_qp, Nm)
    bv              = VectorAssemble(ElemsDof[:], smat[:], N=fes.nDof)
    
    return Am, bv
    
end

#----------------------------------------------------------
#Lq ERRORS:

#Compute Lq error from given theoretical solution:
function LqError(Integ2D::TrInt, uh::Vector{Float64}, 
    FES::TrFES,
    u::FunctionWrapper{Matrix{Float64},Tuple{Vector{Matrix{Float64}}}};
    q::Real=2)
    
    nElems          = Integ2D.mesh.nElems
    nqp             = Integ2D.QRule.nqp
    xim             = Integ2D.QRule.xi
    wdetJ           = Integ2D.wdetJ
    
    #Compute theoretical solution at quadrature nodes:
    x_qp            = QuadNodesCoords(Integ2D)
    utheor_qp       = u(x_qp)
    
    #Compute numerical solution at quadrature nodes:
    unum_qp         = SolutionCompute(uh, FES, xim)
    
    #Error at each qp:
    e_qp            = @mlv    unum_qp-utheor_qp
    
    #Return maximum error if q==Inf:
    if q==Inf
        e_Lq_elems  = maximum(abs.(e_qp), dims=2)[:]
        e_Lq        = maximum(e_Lq_elems)
        return e_Lq, e_Lq_elems
    end
    
    #Areas of each element and total area:
    A_elems         = sum(wdetJ, dims=2)[:]
    A_total         = sum(A_elems)
    
    #Lq error at each element:
    eta_Lq_elems    = @. $sum(abs(e_qp)^q * wdetJ, dims=2)[:]
    e_Lq_elems      = @mlv    (eta_Lq_elems/A_elems)^(1.0/q)
    eta_Lq          = sum(eta_Lq_elems)
    e_Lq            = (eta_Lq/A_total)^(1.0/q)
    
    return e_Lq, e_Lq_elems
    
end

#Compute Lq average of FES:
function LqMean(Integ2D::TrInt, uh::Vector{Float64}, FES::TrFES; q::Real=2)
    
    nElems          = Integ2D.mesh.nElems
    nqp             = Integ2D.QRule.nqp
    xim             = Integ2D.QRule.xi
    wdetJ           = Integ2D.wdetJ
    
    #Compute numerical solution at quadrature nodes:
    u_qp            = SolutionCompute(uh, FES, xim)
    
    #Return maximum value if q==Inf:
    if q==Inf
        u_Lq_elems  = maximum(abs.(u_qp), dims=2)[:]
        u_Lq        = maximum(u_Lq_elems)
        return u_Lq, u_Lq_elems
    end
    
    #Areas of each element and total area:
    A_elems         = sum(wdetJ, dims=2)[:]
    A_total         = sum(A_elems)
    
    #Lq error at each element:
    eta_Lq_elems    = @. $sum(abs(u_qp)^q * wdetJ, dims=2)[:]
    u_Lq_elems      = @mlv    (eta_Lq_elems/A_elems)^(1.0/q)
    eta_Lq          = sum(eta_Lq_elems)
    u_Lq            = (eta_Lq/A_total)^(1.0/q)
    
    return u_Lq, u_Lq_elems
    
end

