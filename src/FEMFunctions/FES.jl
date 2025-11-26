#---------------------------------------------------------
#CONSTRUCTION AND ASSEMBLY OF MATRICES:

#Return matrices N and gradN for several xi's:
function NCompute(fes::Union{TrFES, TrMesh}, xim::Vector{<:AbstractMatrix{Float64}})

    nxi     = length(xim)
    Nm      = Vector{Matrix{Float64}}(undef,nxi)
    for ii=1:nxi
        Nm[ii]  = NCompute(fes, xim[ii])
    end
    return Nm
    
end
function gradNCompute(fes::Union{TrFES, TrMesh}, xim::Vector{<:AbstractMatrix{Float64}})

    nxi     = length(xim)
    gradN1  = Vector{Matrix{Float64}}(undef, nxi)
    gradN2  = Vector{Matrix{Float64}}(undef, nxi)
    for ii=1:nxi
        gradN1[ii], gradN2[ii]  = gradNCompute(fes, xim[ii])
    end
    return [gradN1, gradN2]
    
end

#Given two matrices A(N,m) B(N,n), return a matrix C(N,m*n) such that 
#C(:,ii+(jj-1)*m)=A(:,ii)*B(:,jj):
function UpsilonCompute(A::AbstractMatrix{Float64}, B::AbstractMatrix{Float64})

    N       = size(A,1)
    m       = size(A,2)
    n       = size(B,2)
    C       = zeros(typeof(A[1,1]),N,m*n)
    for jj=1:n, ii=1:m
        @mlv   C[:,ii+(jj-1)*m] = A[:,ii]*B[:,jj]
    end
    return C
    
end
#Idem for vector of matrices. Useful for computing boundary integrals:
function UpsilonCompute(A::Vector{Matrix{Float64}}, B::Vector{Matrix{Float64}})

    Upsilonm    = Vector{Matrix{Float64}}(undef, length(A))
    for ii=1:length(A)
        Upsilonm[ii]    = UpsilonCompute(A[ii], B[ii])
    end
    return Upsilonm
    
end

#Given two FES FES_I and FES_J, we can compute two matrices, imat and jmat, with
#the indices (ii,jj) associated with NComputeMerge.
#The following function returns two matrices imat' and jmat' such that 
#imat'[ii,jj]=imat(FES_I,FES_J), (idem jmat').
#
#The numbering of imat' and jmat' is local. And offset must be added to assemble
#matrices.
function ElemsDofCompute(FES1::FES_Union, FES2::FES_Union; elems::AbstractVector{Int}=1:FES1.mesh.nElems)

    ElemsDof1       = FES1.ElemsDof[elems, :]
    ElemsDof2       = FES2.ElemsDof[elems, :]
    nElems, nDof1   = size(ElemsDof1)
    nDof2           = size(ElemsDof2,2)
    imat            = zeros(Int, nElems, nDof1*nDof2)
    jmat            = zeros(Int, nElems, nDof1*nDof2)
    for jj=1:nDof2, ii=1:nDof1
        kk                  = ii+(jj-1)*nDof1
        @mlv   imat[:,kk]   = ElemsDof1[:,ii]
        @mlv   jmat[:,kk]   = ElemsDof2[:,jj]
    end
    return imat, jmat
    
end

#Function to assemble (i,j) in a vector of N components:
function VectorAssemble(iv::Vector{Int}, sv::Vector{Float64}; N::Int=maximum(iv))

    yv      = zeros(N)
    VectorAssemble!(iv,sv,yv)
    return yv
    
end
# function VectorAssemble!(iv::Vector{Int}, sv::Vector{Float64}, yv::Vector{Float64})
#     @inbounds for ii=1:length(iv)
#         yv[iv[ii]]  += sv[ii]
#     end
#     return
# end
function VectorAssemble!(iv::AbstractVecOrMat{Int}, sv::AbstractVecOrMat{Float64}, yv::AbstractVecOrMat{Float64})
    @inbounds for ii=1:length(iv)
        yv[iv[ii]]  += sv[ii]
    end
    return
end

#Function to multiply matrix <nElems,nqp> with several matrices <nqp,nDof>. Useful for 
#assembly at a boundary integral:
function mult(Binteg::TrBint, alpha::Float64, flux_qp::AbstractMatrix{Float64}, Upsilonm::Vector{Matrix{Float64}})

    smat        = zeros(Binteg.bmesh.nElems, size(Upsilonm[1],2))
    mult!(Binteg, alpha, flux_qp, Upsilonm, 1.0, smat)
    return smat
    
end
function mult!(Binteg::TrBint, alpha::Float64, flux_qp::AbstractMatrix{Float64}, Upsilonm::Vector{Matrix{Float64}}, 
    beta::Float64, smat::AbstractMatrix{Float64})
    
    for face=1:3
        face_belems         = (Binteg.bmesh.ParentFaces.==face)
        smat[face_belems,:] = alpha*(flux_qp[face_belems,:]*Upsilonm[face]) + beta*smat[face_belems,:]
    end
    return
    
end
function mult_test!(Binteg::TrBint, alpha::Float64, flux_qp::AbstractMatrix{Float64}, Upsilonm::Vector{Matrix{Float64}}, 
    beta::Float64, smat::AbstractMatrix{Float64})
    
    nElems              = Binteg.bmesh.nElems
    nqp                 = Binteg.QRule.nqp
    nDof                = size(smat,2)  #May also be nDof*nDof
    ParentFaces         = Binteg.bmesh.ParentFaces
    
#     BLAS.scal!(beta, smat)
#     @inbounds for iDof=1:nDof, iqp=1:nqp, iElem=1:nElems
#         ParentFace          = ParentFaces[iElem]
#         smat[iElem,iDof]    += alpha*flux_qp[iElem,iqp]*Upsilonm[ParentFace][iqp,iDof]
#     end
#     @inbounds for iElem=1:nElems
#         ParentFace              = ParentFaces[iElem]
#         Upsilonm_elem           = Upsilonm[ParentFace]
#         @inbounds for iDof=1:nDof, iqp=1:nqp
#             smat[iElem,iDof]    += flux_qp[iElem,iqp]*Upsilonm_elem[iqp,iDof]
#             smat[iElem,iDof]    += BLAS.dot(nqp, view(flux_qp,iElem,:), 1, 
#                                         view(Upsilonm_elem,:,iDof), 1)
#         end
#     end
#     for face=1:3
#         face_belems         = findall(Binteg.bmesh.ParentFaces.==face)
#         Upsilon             = Upsilonm[face]
#         @inbounds for iDof=1:nDof, iqp=1:nqp, iElem=face_belems
#             smat[iElem,iDof]    += flux_qp[iElem,iqp]*Upsilon[iqp,iDof]
#         end
#     end

#     smat        *= beta
#     smat        += alpha*(flux_qp*Upsilonm[1])
#     mul!(smat, flux_qp, Upsilonm[1], alpha, beta)
#     for face=1:3
#         face_belems         = (Binteg.bmesh.ParentFaces.==face)
#         smat_face           = smat[face_belems, :]
#         BLAS.gemm!('N', 'N', alpha, flux_qp[face_belems, :], Upsilonm[face], 
#                 beta, smat_face)
#         smat[face_belems,:] = smat_face
#     end
    for face=1:3
        BLAS.gemm!('N', 'N', alpha, flux_qp, Upsilonm[face], 
                beta, smat)
    end
    
    return
    
end

#Prolongation matrix from given space fes to corresponding discontinuous
#space with the same basis functions:
function ProlongationMatrix(fes::TrFES)

    iv          = zeros(Int, fes.mesh.nElems, fes.DofPerElem)
    jv          = zeros(Int, fes.mesh.nElems, fes.DofPerElem)
    @inbounds for iElem=1:fes.mesh.nElems, ii=1:fes.DofPerElem
        iv[iElem, ii]   = (iElem-1)*fes.DofPerElem + ii
        jv[iElem, ii]   = fes.ElemsDof[iElem, ii]
    end
    return sparse(view(iv,:), view(jv,:), ones(length(iv)))
    
end
        
#--------------------------------------------------------
#EXTRACT VIEWS FROM SOLUTIONS:

#Given a vector u of size nVars*nDof, extract the views corresponding to
#each variable:
function GetViews(uv::Vector{Float64}, nVars::Int, nDof::Int)

    u               = Vector{VectorView{Float64}}(undef,nVars)
    for ii=1:nVars
        u[ii]       = view(uv, (ii-1)*nDof+1:ii*nDof)
    end
    return u
    
end
function GetViews(uv::MatrixCol{Float64}, nVars::Int, nDof::Int)

    u               = Vector{MatrixColView{Float64}}(undef,nVars)
    for ii=1:nVars
        u[ii]       = view(uv, (ii-1)*nDof+1:ii*nDof)
    end
    return u
    
end

#Given a vector <nVars> of vectors <nDof>, extract views corresponding to the bubble terms:
function GetBubbles(u::Vector{<:GenVector{Float64}}, FES::TrPBSpace)

    nVars           = length(u)
    uB              = Vector{VectorView{Float64}}(undef,nVars)
    for ii=1:nVars
        uB[ii]      = view(u[ii], FES.PSpace.nDof+1:FES.nDof)
    end
    return uB

end

#Given a vector <nVars> of vectors <nDof>, extract views corresponding to the bubble terms:
function GetPolynomials(u::Vector{<:GenVector{Float64}}, FES::TrPBSpace)

    nVars           = length(u)
    uP              = Vector{VectorView{Float64}}(undef,nVars)
    for ii=1:nVars
        uP[ii]      = view(u[ii], 1:FES.PSpace.nDof)
    end
    return uP

end

#---------------------------------------------------------
#EVALUATE SOLUTION AT GIVEN ELEMENTS AND SAME QUADRATURE NODES:

#Given a finite element space, compute numerical solution at given natural coordinates xim 
# at each element.
#NOTE: Here the natural coordinates are the same for each element.
#NOTE: The output is a matrix <nElems, size(xim,1)>
#NOTE: We compute the solution at the elements given by the optional parameter 'elems'
function SolutionCompute(uv::AbstractVector{Float64}, FES::FES_Union, xim::Matrix{Float64};
    elems::AbstractVector{Bool}=ones(Bool, FES.mesh.nElems))

    return SolutionCompute([uv], FES, xim, elems=elems)[1]
    
end

function SolutionGradCompute(uv::AbstractVector{Float64}, FES::FES_Union, xim::Matrix{Float64};
    elems::AbstractVector{Bool}=ones(Bool, FES.mesh.nElems))

    Du_ElemsQp  = SolutionGradCompute([uv], FES, xim, elems=elems)[1]
    
    return [Du_ElemsQp[1,1], Du_ElemsQp[1,2]]
    
end

function SolutionCompute(uv::Vector{<:AbstractVector{Float64}}, FES::FES_Union, xim::Matrix{Float64};
    elems::AbstractVector{Bool}=ones(Bool, FES.mesh.nElems))

    nElems      = sum(elems)
    
    #Matrix <nqp, DofPerElem> with shape functions:
    Nm          = NCompute(FES, xim)
    
    #Allocate matrices u <nElems, nqp> with values of solution at each 
    #quadrature node:
    u_ElemsQp   = Vector{Matrix{Float64}}(undef, length(uv))
    for ii=1:length(uv)
        
        #Matrix <nElems,DofPerElem> with values for each element and dof:
        u_ElemsDof      = uv[ii][FES.ElemsDof[elems,:]]
    
        #Compute umat*Nm^T, which gives a matrix <nElems, ngp>:
        u_ElemsQp[ii]   = zeros(nElems, size(xim,1))
        BLAS.gemm!('N', 'T', 1.0, u_ElemsDof, Nm, 0.0, u_ElemsQp[ii])
    end
    
    return u_ElemsQp
    
end

function SolutionGradCompute(uv::Vector{<:AbstractVector{Float64}}, FES::FES_Union, xim::Matrix{Float64};
    elems::AbstractVector{Bool}=ones(Bool, FES.mesh.nElems))

    nElems      = sum(elems)
    
    #Matrices <nqp, DofPerElem> with gradients of shape functions:
    gradNm      = gradNCompute(FES, xim)
    
    #Allocate matrices u <nElems, nqp> with values of solution at each 
    #quadrature node:
    Du_ElemsQp  = Matrix{Matrix{Float64}}(undef, length(uv), 2)
    for ii=1:length(uv)
        
        #Matrix <nElems,DofPerElem> with values for each element and dof:
        u_ElemsDof          = uv[ii][FES.ElemsDof[elems,:]]
    
        #Compute umat*gradNm^T, which gives a matrix <nElems, ngp>:
        Du_ElemsQp[ii,1]    = zeros(nElems, size(xim,1))
        Du_ElemsQp[ii,2]    = zeros(nElems, size(xim,1))
        BLAS.gemm!('N', 'T', 1.0, u_ElemsDof, gradNm[1], 0.0, Du_ElemsQp[ii,1])
        BLAS.gemm!('N', 'T', 1.0, u_ElemsDof, gradNm[2], 0.0, Du_ElemsQp[ii,2])
        
    end
    
    return Du_ElemsQp
    
end

function Solution_gradCompute(Du::Matrix{MFloat}, Jinv::Matrix{MFloat}) where 
    MFloat<:AbstractVecOrMat{Float64}
   
    nVars       = size(Du,1)
    du          = Matrix{MFloat}(undef,nVars,2)
    for II=1:nVars, jj=1:2
        du[II,jj]           = zeros(size(Du[II,jj]))
    end
    Solution_gradCompute!(Du, Jinv, du)
    return du
    
end

function Solution_gradCompute!(Du::Matrix{MFloat}, Jinv::Matrix{MFloat}, 
    du::Matrix{MFloat}) where MFloat<:AbstractVecOrMat{Float64}
    
    for ii=1:size(Du,1)
        @mlv   du[ii,1]     = Du[ii,1]*Jinv[1,1]+Du[ii,2]*Jinv[2,1]
        @mlv   du[ii,2]     = Du[ii,1]*Jinv[1,2]+Du[ii,2]*Jinv[2,2]
    end
    return
    
end

function Solution_gradCompute(uv::Vector{<:AbstractVector{Float64}}, FES::FES_Union, 
    xim::Matrix{Float64}, Jinv::Matrix{MFloat}) where 
    MFloat<:AbstractVecOrMat{Float64}
   
    Du  = SolutionGradCompute(uv, FES, xim)
    du  = Solution_gradCompute(Du, Jinv)
    return du
    
end

#Compute high-order derivatives, assuming mesh is straight.
function Solution_HODerivs_Compute(uv::Vector{<:GenVector{Float64}}, 
    fes::TrFES, n::Int, xim::Matrix{Float64}; 
    elems::AbstractVector{Bool}=ones(Bool, fes.mesh.nElems))
    
    #Extract variables:
    nVars       = length(uv)
    nElems      = sum(elems)
    nqp         = size(xim,1)
    
    #Allocate output:
    vm          = Array{Matrix{Float64},3}(undef,nVars,n+1,n+1) #PKD coefficients for each variable ande derivative
    HODm        = Array{Matrix{Float64},3}(undef,nVars,n+1,n+1)
    
    #Create quadrature rule in the reference triangle:
    QRule       = WVRule(2*(fes.order+2)) 
    
    #Jacobian inverse at barycenter:
    Jinv        = JinvCompute(fes.mesh, [0.0 0.0], elems=elems)
    
    #Shape functions at quadrature nodes of QRule:
    Nm          = PKD_Nm(QRule.xi, fes.order)
    GradNm      = PKD_gradNm(QRule.xi, fes.order)
    phim        = NCompute(fes, QRule.xi)
    DofPerElem  = div((fes.order+1)*(fes.order+2),2)
    Nm_req      = PKD_Nm(xim, fes.order)
    
    #First, we write the solution in a PKD basis of order N. If 
    #   u   = u_j phi_j = v_j N_j
    #then
    #   v_i = int(N_i * phi_j) u_j
    Nphim       = transpose(Nm)*diagm(QRule.w)*phim
    for II=1:nVars
        um              = uv[II][fes.ElemsDof[elems,:]]
        vm[II,1,1]      = um*transpose(Nphim)
        HODm[II,1,1]    = vm[II,1,1]*transpose(Nm_req)
    end
    
    #Similarly, the derivative of v
    #   dv/dxi_k    = dN_j/dxi_k * v_j
    #can be written as a PKD expansion
    #   dv_dxi_k    = N_j (dv/dxi_k)_j
    #where the coefficients of the derivative are
    #   (dv_dxi_k)_i = int(N_i dN_j/dxi_k) v_j
    N_DN            = Vector{Matrix{Float64}}(undef,2)
    N_DN[1]         = transpose(Nm)*diagm(QRule.w)*GradNm[1]
    N_DN[2]         = transpose(Nm)*diagm(QRule.w)*GradNm[2]
    
    #For each element "e", the coefficients of the derivatives w.r.t. "x" 
    #are
    #   (dv/dx_l)_i = (dv/dxi_m)_i * Jinv_ml = (N_DN[m])_ij * v_j * Jinv_ml
    
    #Compute derivatives (1,0) to (n,0) (i.e., l=1):
    for II=1:nVars, ii=1:n, jj=0

        #Update coefficients of PKD expansion:
        vm[II,ii+1,jj+1]    = zeros(nElems, DofPerElem)
        for ll=1, mm=1:2, jDof=1:DofPerElem, iDof=1:DofPerElem
            @mlv vm[II,ii+1,jj+1][:,iDof]   += $view(Jinv[mm,ll],:,1) *
                                                N_DN[mm][iDof,jDof] *
                                                $view(vm[II,ii,jj+1],:,jDof)
        end
        
        #Evaluate PKD expansion at required points:
        HODm[II,ii+1,jj+1]  = vm[II,ii+1,jj+1]*transpose(Nm_req)

    end

    #Compute derivatives (ii,1) to (ii,ii):
    for II=1:nVars, ii=0:n-1, jj=1:n-ii

        #Update coefficients of PKD expansion:
        vm[II,ii+1,jj+1]    = zeros(nElems, DofPerElem)
        for ll=2, mm=1:2, jDof=1:DofPerElem, iDof=1:DofPerElem
            @mlv vm[II,ii+1,jj+1][:,iDof]   += $view(Jinv[mm,ll],:,1) *
                                                N_DN[mm][iDof,jDof] *
                                                $view(vm[II,ii+1,jj],:,jDof)
        end
        
        #Evaluate PKD expansion at required points:
        HODm[II,ii+1,jj+1]  = vm[II,ii+1,jj+1]*transpose(Nm_req)

    end
    
    return HODm
   
end

#Given a matrix HODm[II,ii+1,jj+1][nElems,1] with the high-order derivatives,
#reconstruct derivatives of order p from derivatives of order p-1. The results 
#are overwritten in HODm:
function HODerivs_Reconstruct!(HODm::Array{Matrix{Float64},3}, mesh::TrMesh, p::Int)

    nVars       = size(HODm,1)
    
    #Create P1 solution:
    P1          = TrPSpace(mesh, 1)
    uP1         = zeros(P1.nDof)
    
    #Note: weights for linear reconstruction are 1/A = M_bar/A_ref
    
    #Jacobian inverse at barycenter:
    Jinv        = JinvCompute(mesh, [0.0 0.0])
    
    #Loop variables: and derivatives of order p-1:
    for II=1:nVars
        
        #Zero derivatives of order p:
        for jj=0:p, ii=p-jj
            @mlv HODm[II,ii+1,jj+1]             = 0.0
        end
        
        #Loop derivatives of order p-1:
        for jj=0:p-1, ii=p-1-jj
        
            #Reconstruct (p-1)th order derivatives of u_I:
            uP1     .= LinearReconstruction(mesh, mesh.metric.lambda_bar,
                        view(HODm[II,ii+1,jj+1],:,1))
                        
            #Compute derivatives:
            DuP1    = SolutionGradCompute([uP1], P1, [0.0 0.0])
            duP1    = Solution_gradCompute(DuP1, Jinv)
            
            #Update high-order derivatives:
            @mlv $view(HODm[II,ii+2,jj+1],:,1)  += duP1[1,1]
            @mlv $view(HODm[II,ii+1,jj+2],:,1)  += duP1[1,2]
            
        end
        
        #Take averages of derivatives jj=1:p-1
        for jj=1:p-1, ii=p-jj
            @mlv $view(HODm[II,ii+1,jj+1],:,1)  *= 0.5
        end
        
    end
    
    return
    
end

#---------------------------------------------------------
#EVALUATE SOLUTION AT ARBITRARY ELEMENTS AND QUADRATURE NODES:

#Given a vector of elements, and a matrix with the values of the 
#natural coordinates, compute the numerical solution.
#NOTE: Here the natural coordinates may not be the same for each element.
#NOTE: The output is a matrix size(elems)
function SolutionCompute(uv::Vector{Float64}, FES::FES_Union, 
    elems::AbstractMatrix{Int}, xim::AbstractArray{Float64,3})
    
    @warn "not used??"
    
    #Reshape elems and xim:
    nElems, nqp = size(elems)
    elems       = reshape(elems, :)
    xim         = reshape(xim, :, 2)
    
    #Evaluate shape functions at quadrature nodes:
    Nm          = NCompute(FES, xim)
    
    #Loop d.o.f.:
    uDof        = zeros(length(elems))
    usol        = zeros(length(elems))
    @inbounds for iDof=1:FES.DofPerElem
        uDof            .= uv[FES.ElemsDof[elems, iDof]]
        @mlv    usol    += uDof*Nm[:,iDof]
    end
    
    #Reshape usol so as to have the same dimensions than elems:
    return reshape(usol, nElems, nqp)
    
end

function SolutionCompute(uv::Vector{<:GenVector{Float64}}, 
    FES::FES_Union, 
    elems::AbstractVector{Int}, xim::AbstractMatrix{Float64})
    
    nVars       = length(uv)
    nElems      = length(elems)
    usol        = Vector{Vector{Float64}}(undef,nVars)
    alloc!(usol, (nElems,))
    SolutionCompute!(uv, FES, elems, xim, usol)
    return usol
    
end

function SolutionCompute!(uv::Vector{<:GenVector{Float64}}, 
    FES::FES_Union, 
    elems::AbstractVector{Int}, xim::AbstractMatrix{Float64},
    usol::Vector{<:GenVector{Float64}})
    
    nVars       = length(uv)
    nElems      = length(elems)
    
    #Evaluate shape functions at quadrature nodes:
    Nm          = NCompute(FES, xim)
    
    #Loop d.o.f.:
    ElemsDof    = zeros(Int, nElems)
    u_ElemsDof  = zeros(nElems)
    zero!(usol)
    for iDof=1:FES.DofPerElem
        ElemsDof            .= FES.ElemsDof[elems, iDof]
        for II=1:nVars
            u_ElemsDof      .= uv[II][ElemsDof]
            @mlv usol[II]   += u_ElemsDof*$view(Nm,:,iDof)
        end
    end
    
    return
    
end

function SolutionCompute(uv::Vector{<:GenVector{Float64}}, 
    fes::TrFES, 
    xm::Vector{<:GenVector{Float64}};
    ElemGuess::AbstractVector{Int}=ones(Int, length(xm[1])),
    xiGuess::Matrix{Float64}=zeros(length(xm[1]),2))
    
    usol        = Vector{Vector{Float64}}(undef,length(uv))
    alloc!(usol, (length(ElemGuess),))
    flag        = SolutionCompute!(uv, fes, xm, usol, ElemGuess=ElemGuess, xiGuess=xiGuess)
    return usol, flag
    
end

function SolutionCompute!(uv::Vector{<:GenVector{Float64}}, 
    fes::TrFES, 
    xm::Vector{<:GenVector{Float64}},
    usol::Vector{<:GenVector{Float64}};
    ElemGuess::AbstractVector{Int}=ones(Int, length(xm[1])),
    xiGuess::Matrix{Float64}=zeros(length(xm[1]),2))
    
    t_ini       = time()
    
    #Search locate x_qp on donor mesh:
    SL_output   = SearchLocate(fes.mesh, [reshape(xm[1],:), reshape(xm[2],:)],
                    ElemGuess=min.(ElemGuess, fes.mesh.nElems), 
                    xiGuess=xiGuess)
#     println("SearchLocate ", time()-t_ini)
    flag        = SL_output[1]
    elem0v      = @mlv abs(SL_output[2])
    xi0m        = SL_output[3]
    
    #Evaluate finite element solution:
    SolutionCompute!(uv, fes, elem0v, xi0m, usol)
    
#     println("SolutionCompute! ", time()-t_ini)
    
    return flag
    
end


#Given a vector of elements, and a matrix with the values of the 
#natural coordinates, compute the gradients of the numerical solution w.r.t. natural coordinates.
#NOTE: Here the natural coordinates may not be the same for each element.
#NOTE: The output is [matrix(size(elems)), matrix(size(elems))]
function SolutionGradCompute(uv::Vector{Float64}, FES::FES_Union, 
    elems::AbstractMatrix{Int}, xim::AbstractArray{Float64,3})
    
    @warn "not used??"
    
    #Reshape elems and xim:
    nElems, nqp = size(elems)
    elems       = reshape(elems, :)
    xim         = reshape(xim, :, 2)
    
    #Evaluate shape functions at quadrature nodes:
    gradNm      = gradNCompute(FES, xim)
    
    #Loop d.o.f.:
    uDof        = zeros(length(elems))
    d1usol      = zeros(length(elems))
    d2usol      = zeros(length(elems))
    @inbounds for iDof=1:FES.DofPerElem
        uDof            .= uv[FES.ElemsDof[elems, iDof]]
        @mlv    d1usol  += uDof*gradNm[1][:,iDof]
        @mlv    d2usol  += uDof*gradNm[2][:,iDof]
    end
    
    #Reshape usol so as to have the same dimensions than elems:
    return [ reshape(d1usol, nElems, nqp), reshape(d2usol, nElems, nqp) ]
    
end

#Evaluate solution at nodes of PSpace:
function SolutionAtNodes(uv::Vector{v}, fes::TrPBSpace) where v<:GenVector{Float64}

    nVars   = length(uv)
    uvv     = Vector{ViewType(uv[1])}(undef,nVars)
    for II=1:nVars
        uvv[II]     = view(uv[II],1:fes.PSpace.nDof)
    end
    return uvv
    
end

#---------------------------------------------
#PLOT SOLUTIONS:

import PyPlot.plot3D

function plot3D(uv::AbstractVector{Float64}, fes::TrPBSpace; color="b", marker=".", markersize=0.5)

    @warn "Deprecated"
    
    plot3D(uv[1:fes.PSpace.nDof], fes.PSpace, marker=marker,
            markersize=markersize, color=color)
            
    return
    
end

function plot3D(uv::AbstractVector{Float64}, fes::TrPSpace; color="b", marker=".", markersize=0.5)

    @warn "Deprecated"
    
    plot3D(fes.NodesCoords[:,1], fes.NodesCoords[:,2],
            uv, linestyle="none", marker=marker,
            markersize=markersize, color=color)
            
    return
    
end

function plot3D(uv::AbstractVector{Float64}, fes::TrQSpace; 
    elems::AbstractVector{Int}=1:fes.mesh.nElems,
    color="b", marker=".", markersize=0.5)

    @warn "Deprecated"
    
    #Lagrangian nodes:
    xi_Lag  = TrElem_EquiNodes(fes.order)
    x_Lag   = LagNodes(fes.mesh, xi_Lag)
    
    #Solution at Lagrangian nodes:
    um      = SolutionCompute(uv, fes, xi_Lag)
    plot3D(x_Lag[1][elems,:][:], x_Lag[2][elems,:][:], um[elems,:][:], 
            linestyle="none", marker=marker, markersize=markersize, color=color)
            
    return
    
end

function plot3D(xm::Vector{Matrix{Float64}},
    fun::FunctionWrapper{Matrix{Float64},Tuple{Vector{Matrix{Float64}}}};
    color="g", marker=".", markersize=0.5)
    
    @warn "Deprecated"
    
    plot3D(xm[1][:], xm[2][:], fun(xm)[:], 
        linestyle="none", marker=marker, markersize=markersize, color=color)
    
    return
    
end

function plot3D(xm::Matrix{Float64},
    fun::FunctionWrapper{Matrix{Float64},Tuple{Vector{Matrix{Float64}}}};
    color="g", marker=".", markersize=0.5)

    @warn "Deprecated"
    
    xmv     = [xm[:,[1]], xm[:,[2]]]
    plot3D(xmv, fun, marker=marker, markersize=markersize, color=color)
    
    return
    
end

function plotNodes(sfun::Function, uv::GenVector{Float64}, fes::TrPBSpace;
    color="b", marker=".", markersize=0.5)
    
    sv  = sfun(view(fes.PSpace.NodesCoords,:,1), 
                view(fes.PSpace.NodesCoords,:,2))
    plot(sv, uv[1:fes.PSpace.nDof], linestyle="none", marker=marker, 
        markersize=markersize, color=color)
    
    return
    
end

function PlotContour(uv::GenVector{Float64}, fes::TrPSpace;
    cmap::String="jet", symmetry::String="none",
    zv::Vector{Float64}=linspace(minimum(uv)-1e-5,maximum(uv)+1e-5,100),
    elems::AbstractVector{Int}=1:fes.mesh.nElems, 
    nSubtriangles::Int=fes.order, PlotBubbles::Bool=false)
    
    #I do not know why this does not work:
    if typeof(fes)==TrPBSpace && !PlotBubbles
        return PlotContour(view(uv,1:fes.PSpace.nDof), fes.PSpace, cmap=cmap, symmetry=symmetry,
            zv=zv, elems=elems, nSubtriangles=nSubtriangles)
    end
    
    mesh                = fes.mesh
    nElems              = length(elems)
    elems_bool          = zeros(Bool, mesh.nElems)
    elems_bool[elems]   .= true
    
    #Interpolate solution onto equispaced mesh:
    xim             = TrElem_EquiNodes_sorted(nSubtriangles)
    xm              = view(mesh.ElemsCoords, elems, 1:2:2*mesh.NodesPerElem-1)
    ym              = view(mesh.ElemsCoords, elems, 2:2:2*mesh.NodesPerElem)
    Nm              = NCompute(mesh, xim)
    ElemsCoords1    = xm*transpose(Nm)
    ElemsCoords2    = ym*transpose(Nm)
    uEC             = SolutionCompute(uv, fes, xim, elems=elems_bool)
    ElemsDof        = zeros(Int, nElems, fes.DofPerElem)
    for ii=1:nElems, jj=1:fes.DofPerElem
        ElemsDof[ii,jj]     = (ii-1)*fes.DofPerElem + jj
    end
    
    #Function that assigns a number to each (i,j) node in the triangle:
    kfun(i,j)       = div((j-1)*j,2)+i
    
    #Local dof for each reference subtriangle:
    LocalDof            = zeros(Int, nSubtriangles^2, 3)
    #Upside triangles:
    kk                  = 0
    for jj=1:nSubtriangles, ii=1:jj
        kk              +=1 
        LocalDof[kk,1]  = kfun(ii,jj)
        LocalDof[kk,2]  = kfun(ii,jj+1)
        LocalDof[kk,3]  = kfun(ii+1,jj+1)
    end
    #Downside triangles:
    for jj=1:nSubtriangles, ii=1:jj-1
        kk              += 1
        LocalDof[kk,1]  = kfun(ii,jj)
        LocalDof[kk,2]  = kfun(ii+1,jj+1)
        LocalDof[kk,3]  = kfun(ii+1,jj)
    end
    LocalDofs           = kfun(nSubtriangles+1, nSubtriangles+1)
    
    #For each triangle, we have nSubtriangles^2 subtriangles:
    trianglesm      = zeros(Int, nElems*nSubtriangles^2, 3)
    iTr             = 0
    for iElem=1:nElems, iSubTr=1:nSubtriangles^2
        iTr                 += 1
        aux                 = (iElem-1)*LocalDofs
        trianglesm[iTr,1]   = LocalDof[iSubTr,1]+aux
        trianglesm[iTr,2]   = LocalDof[iSubTr,2]+aux
        trianglesm[iTr,3]   = LocalDof[iSubTr,3]+aux
    end
    
    #Plot:
    EC1             = reshape(transpose(ElemsCoords1),:)
    EC2             = reshape(transpose(ElemsCoords2),:)
    usol            = reshape(transpose(uEC),:)
    maskv           = zeros(Bool, nElems*nSubtriangles^2)
    tricontourf(EC1, EC2, usol, zv, triangles=trianglesm.-1, 
        mask=maskv, cmap=cmap)
    if symmetry=="x" || symmetry=="xy"
        tricontourf(-EC1, EC2, usol, zv, triangles=trianglesm.-1, cmap=cmap)
    end
    if symmetry=="y" || symmetry=="xy"
        tricontourf(EC1, -EC2, usol, zv, triangles=trianglesm.-1, cmap=cmap)
    end
    if symmetry=="xy"
        tricontourf(-EC1, -EC2, usol, zv, triangles=trianglesm.-1, cmap=cmap)
    end    
        
    return usol
    
end

#=
function PlotContour(uv::GenVector{Float64}, fes::TrPSpace;
    cmap::String="jet", symmetry::String="none",
    z1::Float64=minimum(uv)-1e-10, 
    z2::Float64=maximum(uv)+1e-10, 
    Nz::Int=100)

    PlotContour(uv, fes, cmap=cmap, symmetry=symmetry, zv=linspace(z1, z2, Nz)
    
    return
    
end
=#

#=
function PlotContour(uv::GenVector{Float64}, fes::TrPBSpace;
    cmap::String="jet", symmetry::String="none", 
    z1::Float64=minimum(view(uv,1:fes.PSpace.nDof))-1e-10, 
    z2::Float64=maximum(view(uv,1:fes.PSpace.nDof))+1e-10, 
    Nz::Int=100)
    
    @warn "Deprecated"
    
    PlotContour(view(uv,1:fes.PSpace.nDof), fes.PSpace, cmap=cmap, symmetry=symmetry,
        zv=linspace(z1,z2,Nz))
    return
    
end
=#

function PlotContour(uv::GenVector{Float64}, fes::TrPBSpace;
    cmap::String="jet", symmetry::String="none", 
    zv::Vector{Float64}=linspace(   minimum(view(uv,1:fes.PSpace.nDof))-1e-5, 
                                    maximum(view(uv,1:fes.PSpace.nDof))+1e-5, 100 ),
    elems::AbstractVector{Int}=1:fes.mesh.nElems, 
    nSubtriangles::Int=fes.order)
    
    return PlotContour(view(uv,1:fes.PSpace.nDof), fes.PSpace, cmap=cmap, symmetry=symmetry,
        zv=zv, elems=elems, nSubtriangles=nSubtriangles)
    
end

function Plot3D(uv::GenVector{Float64}, fes::TrPSpace;
    cmap::String="jet",
    z1::Float64=minimum(uv)-1e-10, 
    z2::Float64=maximum(uv)+1e-10, 
    Nz::Int=100, 
    elems::AbstractVector{Int}=1:fes.mesh.nElems)
    
    zv                  = linspace(z1, z2, Nz)
    mesh                = fes.mesh
    nElems              = length(elems)
    elems_bool          = zeros(Bool, mesh.nElems)
    elems_bool[elems]   .= true
    
    #Interpolate solution onto equispaced mesh:
    xim             = TrElem_EquiNodes_sorted(fes.order)
    xm              = view(mesh.ElemsCoords ,elems, 1:2:2*mesh.NodesPerElem-1)
    ym              = view(mesh.ElemsCoords, elems, 2:2:2*mesh.NodesPerElem)
    Nm              = NCompute(mesh, xim)
    ElemsCoords1    = xm*transpose(Nm)
    ElemsCoords2    = ym*transpose(Nm)
    uEC             = SolutionCompute(uv, fes, xim, elems=elems_bool)
    ElemsDof        = zeros(Int, nElems, fes.DofPerElem)
    for ii=1:nElems, jj=1:fes.DofPerElem
        ElemsDof[ii,jj]     = (ii-1)*fes.DofPerElem + jj
    end
    
    #Function that assigns a number to each (i,j) node in the triangle:
    kfun(i,j)       = div((j-1)*j,2)+i
    
    #Local dof for each reference subtriangle:
    LocalDof            = zeros(Int, fes.order^2, 3)
    #Upside triangles:
    kk                  = 0
    for jj=1:fes.order, ii=1:jj
        kk              +=1 
        LocalDof[kk,1]  = kfun(ii,jj)
        LocalDof[kk,2]  = kfun(ii,jj+1)
        LocalDof[kk,3]  = kfun(ii+1,jj+1)
    end
    #Downside triangles:
    for jj=1:fes.order, ii=1:jj-1
        kk              += 1
        LocalDof[kk,1]  = kfun(ii,jj)
        LocalDof[kk,2]  = kfun(ii+1,jj+1)
        LocalDof[kk,3]  = kfun(ii+1,jj)
    end
    
    #For each triangle, we have fes.order^2 subtriangles:
    trianglesm      = zeros(Int, nElems*fes.order^2, 3)
    iTr             = 0
    for iElem=1:nElems, iSubTr=1:fes.order^2
        iTr                 += 1
        aux                 = (iElem-1)*fes.DofPerElem
        trianglesm[iTr,1]   = LocalDof[iSubTr,1]+aux
        trianglesm[iTr,2]   = LocalDof[iSubTr,2]+aux
        trianglesm[iTr,3]   = LocalDof[iSubTr,3]+aux
    end
    
#     figure()
#     plot3D(ElemsCoords1[:], ElemsCoords2[:], usol[:], ".b")
#     error("")
    
    #Plot:
    EC1             = reshape(transpose(ElemsCoords1),:)
    EC2             = reshape(transpose(ElemsCoords2),:)
    usol            = reshape(transpose(uEC),:)
    maskv           = zeros(Bool, nElems*fes.order^2)
    plot_trisurf(EC1, EC2, usol, triangles=trianglesm.-1, 
        mask=maskv, cmap=cmap, antialiased=false) #antialiased to make surface opaque
        
    return
    
end

function Plot3D(uv::GenVector{Float64}, fes::TrPBSpace;
    cmap::String="jet",
    z1::Float64=minimum(view(uv,1:fes.PSpace.nDof)), 
    z2::Float64=maximum(view(uv,1:fes.PSpace.nDof)), 
    Nz::Int=100, 
    elems::AbstractVector{Int}=1:fes.mesh.nElems)
    
    Plot3D(view(uv,1:fes.PSpace.nDof), fes.PSpace, cmap=cmap,
        z1=z1, z2=z2, Nz=Nz, 
        elems=elems)
    return
    
end

#-----------------------------------------------
#INTERPOLATE FUNCTIONS:

function Interpolate!(fun::FW11, 
    u::Vector{<:GenVector{Float64}}, fes::TrPBSpace)
    
    xm      = fes.PSpace.NodesCoords
    u0      = fun([xm[:,1], xm[:,2]])
    for ii=1:length(u)
        u[ii][1:fes.PSpace.nDof]    .= u0[ii][:]
    end
    return
    
end
