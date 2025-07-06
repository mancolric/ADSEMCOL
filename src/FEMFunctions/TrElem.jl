const TrElem_Area=(3.0*sqrt(3.0)/4.0)::Float64         #Area of the reference equilateral triangle 
const LineElem_Area=2.0::Float64

#Return linear Lagrangian shape functions in equilateral triangle:
function TrElem_LinearLag(xi0m::AbstractMatrix{Float64})

    Lm      = zeros(size(xi0m,1), 3)
    
    #Note that the edges of the equilateral
    #triangle are given by
    #   xi2/(-0.5)      = 1
    #   xi1/d + xi2     = 1
    #   -xi1/d + xi2    = 1
    #with (d,0) and (-d,0) two points of the triangle, given by
    #d/1 = sqrt(3)/2 / (1+1/2)
    dinv            = 3.0/sqrt(3)
    sq32            = sqrt(3.0)/2.0
    @mlv    Lm[:,1] = (xi0m[:,1]*dinv + xi0m[:,2] - 1.0)/(-sq32*dinv - 0.5 - 1.0)
    @mlv    Lm[:,2] = (-xi0m[:,1]*dinv + xi0m[:,2] - 1.0)/(-sq32*dinv - 0.5 - 1.0)
    @mlv    Lm[:,3] = (-xi0m[:,2]*2.0 - 1.0)/(-3.0)
    
    return Lm
    
end

#Transform from equilateral triangle 
#   (-sqrt(3)/2, -1/2) -- (sqrt(3)/2, 1/2 -- (0,1)
#to triangle (-1,-1) -- (1,-1) -- (-1,1).
function TrElem_Equi_to_Iso2(xi0m::AbstractMatrix{Float64})

    xim     = zeros(size(xi0m))
    
    #Use Lagrangian interpolation:
    Lm      = TrElem_LinearLag(xi0m)
    
    @mlv    xim[:,1]    = -1.0*Lm[:,1] + 1.0*Lm[:,2] - 1.0*Lm[:,3]
    @mlv    xim[:,2]    = -1.0*Lm[:,1] - 1.0*Lm[:,2] + 1.0*Lm[:,3]
    
    return xim
    
end

#Return Jacobian of transformation:
function TrElem_Equi_to_Iso2_Jacobian()

    #Derivatives of Lagrange polynomials in the equilateral triangle:
    dinv    = 3.0/sqrt(3)
    sq32    = sqrt(3.0)/2.0
    d1L     = [ dinv/(-sq32*dinv - 0.5 - 1.0), 
                -dinv/(-sq32*dinv - 0.5 - 1.0), 
                0.0 ]
    d2L     = [ 1.0/(-sq32*dinv - 0.5 - 1.0), 
                1.0/(-sq32*dinv - 0.5 - 1.0), 
                2.0/3.0 ]
    dL_dxi  = hcat(d1L, d2L) #dL_k/dxi_j
    
    #X at the Lagrangian nodes:
    Xm      = [ -1.0 1.0 -1.0
                -1.0 -1.0 1.0 ]
    
    #Jacobian of transformation:
    #   dX_i/dxi_j  = d/dxi_j (X_i(xi_k) * L_k(xi)) = X_i(xi_k) * dL_k/dxi_j(xi)
    J       = Xm*dL_dxi
    
    return J
    
end

#Transform from isosceles triangle 
#   (0,0) -- (1,0) -- (0,1)
#to equilateral one
#   (-sqrt(3)/2, -1/2) -- (sqrt(3)/2, 1/2 -- (0,1)
function TrElem_Iso1_to_Equi(xi0m::AbstractMatrix{Float64})

    xim                 = zeros(size(xi0m))
    
    #Use Lagrangian interpolation:
    L1v                 = @mlv    1.0-xi0m[:,1]-xi0m[:,2]
    L2v                 = @mlv    xi0m[:,1]
    L3v                 = @mlv    xi0m[:,2]
    sq32                = sqrt(3.0)/2.0
    @mlv    xim[:,1]    = -sq32*L1v + sq32*L2v
    @mlv    xim[:,2]    = -0.5*L1v - 0.5*L2v + 1.0*L3v
    
    return xim
    
end

#Transform from isosceles triangle 
#   (-1,-1) -- (1,-1) -- (-1,1)
#to equilateral one
#   (-sqrt(3)/2, -1/2) -- (sqrt(3)/2, 1/2 -- (0,1)
function TrElem_Iso2_to_Equi(xi0m::AbstractMatrix{Float64})

    xim                 = zeros(size(xi0m))
    
    #Use Lagrangian interpolation:
    L1v                 = @mlv    (xi0m[:,1]+xi0m[:,2])/(-2.0)
    L2v                 = @mlv    (xi0m[:,1]+1.0)/(2.0)
    L3v                 = @mlv    (xi0m[:,2]+1.0)/(2.0)
    sq32                = sqrt(3.0)/2.0
    @mlv    xim[:,1]    = -sq32*L1v + sq32*L2v
    @mlv    xim[:,2]    = -0.5*L1v - 0.5*L2v + 1.0*L3v
    
    return xim
    
end

#Compute Chebyshev nodes of the second kind on the triangle:
function TrElem_ChebyshevNodes(r::Int)

    #Return barycenter if r==0:
    if r==0
        return [0.0 0.0]
    end

    #Compute number of nodes:
    nNodes          = (r+1)*(r+2) ÷ 2

    #Allocate output:
	xim             = zeros(nNodes,2)

    #Get Chebyshev nodes on the interval [0,1]:
    nuv             = @. 0.5 - 0.5*cos((0:r)*pi/r)

    #Auxiliary functions to compute coordinates:
    xifun(i,j)      = [ (1+2*nuv[i]-nuv[j]-nuv[r+3-i-j])/3.0,
                        (1-nuv[i]+2*nuv[j]-nuv[r+3-i-j])/3.0 ]

    #Vertices:
    xim[1,:]        = xifun(1,1)
    xim[2,:]        = xifun(r+1,1)
    xim[3,:]        = xifun(1,r+1)

    #Edge nodes:
    iNode           = 4
    for ii=2:r
        xim[iNode,:]    = xifun(ii,1)
        iNode           += 1
    end
    for jj=2:r
        xim[iNode,:]    = xifun(r+2-jj,jj)
        iNode           += 1
    end
    for jj=r:-1:2
        xim[iNode,:]    = xifun(1,jj)
        iNode           += 1
    end

    #Inner points:
    # for jj=2:r-1, ii=2:r+1-jj
    #     xim[iNode,:]    = xifun(ii,jj)
    #     iNode           += 1
    # end
    ij0     = 2
    while iNode<=nNodes
        for ii=ij0:r+3-2*ij0
            xim[iNode,:]    = xifun(ii,ij0)
            iNode           += 1
        end
        for jj=ij0+1:r+2-2*ij0
            xim[iNode,:]    = xifun(r+3-ij0-jj,jj)
            iNode           += 1
        end
        for jj=r+3-2*ij0:-1:ij0+1
            xim[iNode,:]    = xifun(ij0,jj)
            iNode           += 1
        end
        ij0                 += 1
    end

    #Correct machine-precision errors:
    @. xim[abs(xim)<1e-15]  = 0.0

    #Return:
    return TrElem_Iso1_to_Equi(xim)

end

#Set of equispaced nodes in the triangle of order "r":
function TrElem_EquiNodes(r::Int)

    #Return barycenter if r==0:
    if r==0
        return [0.0 0.0]
    end
    
    #Allocate:
    nNodes  = div((r+1)*(r+2),2)
    xim     = zeros(nNodes,2)
    
    #Vertices:
    xim[1,1]    = 0.0
    xim[1,2]    = 1.0
    xim[2,1]    = -sqrt(3)/2
    xim[2,2]    = -0.5
    xim[3,1]    = sqrt(3)/2
    xim[3,2]    = -0.5

    #Function that returns coordinate for node at row j, column i:
    function xiNode_(iNode::Int, jNode::Int)
        y_ij    = 1.0 - 1.5/r*(jNode-1)
        x0      = (1.0-y_ij)/sqrt(3)
        x_ij    = -x0 + 2*x0/(jNode-1)*(iNode-1)
        return x_ij, y_ij
    end
    
    #Edge nodes:
    kk                          = 3
    for jj=2:r
        kk                      += 1
        xim[kk,1], xim[kk,2]    = xiNode_(1, jj)
    end
    for ii=2:r
        kk                      += 1
        xim[kk,1], xim[kk,2]    = xiNode_(ii, r+1)
    end
    for jj=r:-1:2
        kk                      += 1
        xim[kk,1], xim[kk,2]    = xiNode_(jj, jj)
    end
    
    #Inner nodes:
    for jj=3:r, ii=2:jj-1
        kk                      += 1
        xim[kk,1], xim[kk,2]    = xiNode_(ii,jj)
    end
    
    return xim
    
end

function TrElem_EquiNodes_sorted(r::Int)

    #Return error if r==0:
    if r==0
        error("r must be positive")
    end
    
    #Allocate:
    nNodes  = div((r+1)*(r+2),2)
    xim     = zeros(nNodes,2)
    
    #Function that returns coordinate for node at row j, column i:
    function xiNode_(iNode::Int, jNode::Int)
        if jNode==1
            return 0.0, 1.0
        end
        y_ij    = 1.0 - 1.5/r*(jNode-1)
        x0      = (1.0-y_ij)/sqrt(3)
        x_ij    = -x0 + 2*x0/(jNode-1)*(iNode-1)
        return x_ij, y_ij
    end
    
    #Loop:
    kk      = 0
    for jj=1:r+1, ii=1:jj
        kk                      += 1
        xim[kk,1], xim[kk,2]    = xiNode_(ii,jj)
    end
    
    return xim
    
end

#-------------------------------------------------------------------------------------
#PKD basis:

function PKD_Nm(xi0m::AbstractMatrix{Float64}, r::Int)

    #Change coordinates from equilateral triangle to triangle (-1,-1) -- (1,-1) -- (-1,1):
    xim         = TrElem_Equi_to_Iso2(xi0m)
    xi          = view(xim,:,1)
    eta         = view(xim,:,2)

    #Allocate memory:
    Nxi         = size(xim,1)
    Ndof        = (r+1)*(r+2)÷2
    Nm          = zeros(Nxi, Ndof)

    #Compute shape functions. We use formula from Giraldo & Warburton (JCP, 2005):
    kk          = 0
    for rr=0:r, jj=0:rr, ii=rr-jj
        
        kk          += 1
        
        #Compute P_j^(2i+1,0)(eta):
        Pj          = zeros(Nxi)
        alpha       = 2*ii+1
        if jj==0
            @. Pj   = 1.0
        elseif jj==1
            @. Pj   = (alpha+1) + (alpha+2)*0.5*(eta-1.0)
        else
            Pjm2    = ones(Nxi)
            Pjm1    = @. (alpha+1) + (alpha+2)*0.5*(eta-1.0)
            for nn=2:jj
                a           = nn+alpha
                b           = nn
                c           = a+b
                @. Pj        = @. ((c-1)*(c*(c-2)*eta+(a-b)*(c-2*nn))*Pjm1 - 2*(a-1)*(b-1)*c*Pjm2) / 
                                    (2*nn*(c-nn)*(c-2))
                @. Pjm2     = Pjm1
                @. Pjm1     = Pj
            end
        end

        #Compute P_i^(0,0)((2xi+eta+1)/(1-eta))*(0.5*(1-eta))^i. Let z=(2xi+eta+1)/(1-eta)
        Pi      = zeros(Nxi)
        alpha   = 0
        if ii==0
            @. Pi       = 1.0
        elseif ii==1
            #Note that 0.5*(z-1)*(0.5*(1-eta))=0.25*(z*(1-eta)-1*(1-eta))=0.25*(2xi+eta+1-1+eta)
            @. Pi       = (alpha+1)*0.5*(1.0-eta) + (alpha+2)*0.25*(2*xi+2*eta)
        else
            Pim2    = @. 0.5*(1.0-eta)                                          #P0*0.5*(1-eta)
            Pim1    = @. (alpha+1)*0.5*(1.0-eta) + (alpha+2)*0.25*(2*xi+2*eta)  #P1*0.5*(1-eta)
            for nn=2:ii
                a           = nn+alpha
                b           = nn
                c           = a+b
                #Apply recurrence relation and multiply by 0.5*(1-eta):
                @. Pi        = @. ((c-1)*(c*(c-2)*(2*xi+eta+1)*0.5+(a-b)*(c-2*nn)*0.5*(1-eta))*Pim1 - 2*(a-1)*(b-1)*c*Pim2*0.5*(1-eta)) /
                                    (2*nn*(c-nn)*(c-2))
                @. Pim2     = Pim1 * 0.5*(1-eta)
                @. Pim1     = Pi
            end
        end

        #Save results:
        #NOTE: The constant in Giraldo and and Warbunton is chosen so that
        #the integral of phi_k*phi_k is 1 over a triangle of area 4.
        #We multiply this constant by C=sqrt(4/(sqrt(3)/4)) so that the integral in
        #the equilateral triangle is
        #   int_(equil) phi_k phi_k = C^2 |equil|/|GW| int_(GW) phiGW_k phiGW_k = 1
        #NOTE: The constant here is the constant from Giraldo and Warbunton multiplied by 2, 
        #so as the integral of phi_k*phi_k over our reference triangle equals 1:
        C2                  = 2.0/TrElem_Area
        @. Nm[:,kk]         = sqrt(0.5*(2*ii+1)*(ii+jj+1)*C2) * Pi * Pj

    end

    return Nm

end

function PKD_gradNm(xi0m::AbstractMatrix{Float64}, r::Int)

    #Change coordinates from equilateral triangle 
    #to triangle (-1,-1) -- (1,-1) -- (-1,1):
    xim         = TrElem_Equi_to_Iso2(xi0m)
    Nxi         = size(xim,1)
    xi          = view(xim,:,1)
    eta         = view(xim,:,2)

    #Allocate memory:
    Ndof        = (r+1)*(r+2)÷2
    dN_dxi1     = zeros(Nxi)        #Derivatives w.r.t. coordinates in isosceles triangle ("xi")
    dN_dxi2     = zeros(Nxi)
    grad1Nm     = zeros(Nxi, Ndof)  #Derivatives w.r.t. coordinates in equilateral triangle ("X")
    grad2Nm     = zeros(Nxi, Ndof)
    
    #Derivatives of coordinates in isosceles triangle  w.r.t. coordinates in
    #equilateral one:
    dxi_dX      = TrElem_Equi_to_Iso2_Jacobian()
    
    #Compute shape functions. We use formula from Giraldo & Warburton (JCP, 2005):
    kk          = 0
    C2          = 2.0/TrElem_Area
    for rr=0:r, jj=0:rr, ii=rr-jj
        
        kk          += 1
        
        #Compute P_j^(2i+1,0)(eta):
        Pj          = zeros(Nxi)
        dPj_deta    = zeros(Nxi)
        alpha       = 2*ii+1
        if jj==0
            @. Pj           = 1.0
            @. dPj_deta     = 0.0
        elseif jj==1
            @. Pj           = (alpha+1) + (alpha+2)*0.5*(eta-1.0)
            @. dPj_deta     = (alpha+2)*0.5
        else
            Pjm2            = ones(Nxi)
            dPjm2_deta      = zeros(Nxi)
            Pjm1            = @. (alpha+1) + (alpha+2)*0.5*(eta-1.0)
            dPjm1_deta      = fill((alpha+2)*0.5, Nxi)
            for nn=2:jj
                a               = nn+alpha
                b               = nn
                c               = a+b
                @. Pj           = ((c-1)*(c*(c-2)*eta+(a-b)*(c-2*nn))*Pjm1 - 2*(a-1)*(b-1)*c*Pjm2) / 
                                    (2*nn*(c-nn)*(c-2))
                @. dPj_deta     = ((c-1)*(c*(c-2))*Pjm1 + (c-1)*(c*(c-2)*eta+(a-b)*(c-2*nn))*dPjm1_deta - 2*(a-1)*(b-1)*c*dPjm2_deta) / 
                                    (2*nn*(c-nn)*(c-2))
                @. Pjm2         = Pjm1
                @. dPjm2_deta   = dPjm1_deta
                @. Pjm1         = Pj
                @. dPjm1_deta   = dPj_deta
            end
        end

        #Compute P_i^(0,0)((2xi+eta+1)/(1-eta))*(0.5*(1-eta))^i. Let z=(2xi+eta+1)/(1-eta)
        Pi          = zeros(Nxi)
        dPi_dxi     = zeros(Nxi)
        dPi_deta    = zeros(Nxi)
        alpha       = 0
        if ii==0
            @. Pi       = 1.0
            @. dPi_dxi  = 0.0
            @. dPi_deta = 0.0
        elseif ii==1
            #Note that 0.5*(z-1)*(0.5*(1-eta))=0.25*(z*(1-eta)-1*(1-eta))=0.25*(2xi+eta+1-1+eta)
            @. Pi       = (alpha+1)*0.5*(1.0-eta) + (alpha+2)*0.25*(2*xi+2*eta)
            @. dPi_dxi  = (alpha+2)*0.25*2
            @. dPi_deta = (alpha+1)*0.5*(-1.0) + (alpha+2)*0.25*2
        else
            Pim2            = @. 0.5*(1.0-eta)                                          #P0*0.5*(1-eta)
            dPim2_dxi       = zeros(Nxi)
            dPim2_deta      = fill(-0.5,Nxi)
            Pim1            = @. (alpha+1)*0.5*(1.0-eta) + (alpha+2)*0.25*(2*xi+2*eta)  #P1*0.5*(1-eta)
            dPim1_dxi       = fill((alpha+2)*0.25*2,Nxi)
            dPim1_deta      = fill((alpha+1)*0.5*(-1.0) + (alpha+2)*0.25*2, Nxi)
            for nn=2:ii
                a               = nn+alpha
                b               = nn
                c               = a+b
                #Apply recurrence relation and multiply by 0.5*(1-eta):
                @. Pi           = ((c-1)*(c*(c-2)*(2*xi+eta+1)*0.5+(a-b)*(c-2*nn)*0.5*(1-eta))*Pim1 - 2*(a-1)*(b-1)*c*Pim2*0.5*(1-eta)) /
                                    (2*nn*(c-nn)*(c-2))
                @. dPi_dxi      = ((c-1)*(c*(c-2)*(2)*0.5)*Pim1 + (c-1)*(c*(c-2)*(2*xi+eta+1)*0.5+(a-b)*(c-2*nn)*0.5*(1-eta))*dPim1_dxi - 2*(a-1)*(b-1)*c*dPim2_dxi*0.5*(1-eta)) /
                                    (2*nn*(c-nn)*(c-2))
                @. dPi_deta     = ((c-1)*(c*(c-2)*(1.0)*0.5+(a-b)*(c-2*nn)*0.5*(-1.0))*Pim1 + (c-1)*(c*(c-2)*(2*xi+eta+1)*0.5+(a-b)*(c-2*nn)*0.5*(1-eta))*dPim1_deta - 2*(a-1)*(b-1)*c*dPim2_deta*0.5*(1-eta) - 2*(a-1)*(b-1)*c*Pim2*0.5*(-1.0)) /
                                    (2*nn*(c-nn)*(c-2))                                    
                @. Pim2         = Pim1 * 0.5*(1-eta)
                @. dPim2_dxi    = dPim1_dxi * 0.5*(1-eta)
                @. dPim2_deta   = dPim1_deta * 0.5*(1-eta) + Pim1 * 0.5*(-1.0)
                @. Pim1         = Pi
                @. dPim1_dxi    = dPi_dxi
                @. dPim1_deta   = dPi_deta
            end
        end

        #Save results:
        #NOTE: The constant here is the constant from Giraldo and Warbunton multiplied by 2, 
        #so as the integral of phi_k*phi_k over our reference triangle equals 1. Also, the derivatives w.r.t. xi and eta are multiplied
        #by 2 because we work in the triangle (0,0) - (1,0) - (0,1).
        @. dN_dxi1              = sqrt(0.5*(2*ii+1)*(ii+jj+1)*C2) * (dPi_dxi * Pj)
        @. dN_dxi2              = sqrt(0.5*(2*ii+1)*(ii+jj+1)*C2) * (dPi_deta * Pj + Pi * dPj_deta)
        @. grad1Nm[:,kk]        = dN_dxi1*dxi_dX[1,1] + dN_dxi2*dxi_dX[2,1]
        @. grad2Nm[:,kk]        = dN_dxi1*dxi_dX[1,2] + dN_dxi2*dxi_dX[2,2]
        
    end

    return [ grad1Nm, grad2Nm ]

end

#Return index "k" associated with (i,j) pair, i+j<=order. 
#NOTE: k=1 for i=j=0:
#index[ii+1,jj+1] corresponds to (ii,jj) pair.
function PKD_index(order::Int)
    
    indexm      = zeros(Int, order+1, order+1)
    kk          = 0
    for pp=0:order, jj=0:pp, ii=pp-jj
        kk                  += 1
        indexm[ii+1,jj+1]   = kk
    end
    
    return indexm
    
end

#-----------------------------------------------------------------------------
#Lagrangian shape functions:

#Write L_k(x_i) as P_j(x_i) a_jk, with P the PKD polynomials (2D case).
#xim are the coordinates of the Lagrangian nodes:
function Lag_PKDCoefficients(order::Int, xim::AbstractMatrix{Float64})

    return PKD_Nm(xim,order)\eye(size(xim,1))
    
end

#Write L_k(x_i) as P_j(x_i) a_jk, with P the Legendre polynomials (1D case):
#xiv are the coordinates of the Lagrangian nodes:
function Lag_LegCoefficients(order::Int, xiv::AbstractVector{Float64})

    return PolyLeg(xiv,order)\eye(length(xiv))
    
end

#Create PSpace of given order associated with given mesh:
function TrPSpace(mesh::TrMesh, order::Int)

    if order<=0
        error("order must be positive")
    end
    
    PSpace              = TrPSpace()
    PSpace.mesh         = mesh
    PSpace.order        = order
    PSpace.DofPerElem   = div((order+1)*(order+2),2)
    PSpace.nDof         = mesh.nVerts + mesh.nEdges*(order-1) + 
                            mesh.nElems*(PSpace.DofPerElem-3*order)
    PSpace.nDofU        = mesh.nVerts + mesh.nEdges*(order-1)
    PSpace.xi           = TrElem_ChebyshevNodes(order)
    PSpace.Lag_PKD      = Lag_PKDCoefficients(order, PSpace.xi)
    
    #Set ElemsDof and NodesCoords at mesh edges:
    xiv_mesh            = ChebyNodes(-1.0, 1.0, mesh.order+1)
    xiv_mesh            = vcat(-1.0, 1.0, xiv_mesh[2:mesh.order]) #Careful: Lagrangian nodes are defined in this order
    Lag_Leg             = Lag_LegCoefficients(mesh.order, xiv_mesh)
    xiv_dof             = ChebyNodes(-1.0, 1.0, PSpace.order+1)[2:PSpace.order]
    Lm                  = PolyLeg(xiv_dof, mesh.order) * Lag_Leg
    PSpace.ElemsDof     = hcat(mesh.ElemsNodes[:,1:3],
                            zeros(mesh.nElems, PSpace.DofPerElem-3))
    PSpace.NodesCoords  = vcat(mesh.NodesCoords[1:mesh.nVerts,:], 
                            zeros(PSpace.nDof-mesh.nVerts, 2))
    @inbounds for iEdge=1:mesh.nEdges
    
        #Mesh nodes and elements:
        mesh_nodes                      = mesh.EdgesNodes[iEdge,:]
        x1v                             = view(mesh.NodesCoords,mesh_nodes,1)
        x2v                             = view(mesh.NodesCoords,mesh_nodes,2)
        elem1                           = mesh.EdgesElems[iEdge,1]
        face1                           = mesh.EdgesFaces[iEdge,1]
        elem2                           = mesh.EdgesElems[iEdge,2]
        face2                           = mesh.EdgesFaces[iEdge,2]
        
        #Insert dof:
        for iDof=1:PSpace.order-1
            
            dof                         = mesh.nVerts + (iEdge-1)*(PSpace.order-1) + iDof
            for iNode=1:mesh.order+1
                PSpace.NodesCoords[dof,1]   += Lm[iDof,iNode]*x1v[iNode]
                PSpace.NodesCoords[dof,2]   += Lm[iDof,iNode]*x2v[iNode]
            end
            PSpace.ElemsDof[elem1, 3+(face1-1)*(PSpace.order-1)+iDof]  = dof
            if elem2>0
                PSpace.ElemsDof[elem2, 3+face2*(PSpace.order-1)+1-iDof]= dof
            end
        
        end
        
    end
    
    #Set ElemsDof and NodesCoords inside elements:
    xim_mesh    = mesh.xi
    xim_dof     = TrElem_ChebyshevNodes(PSpace.order)[3*PSpace.order+1:PSpace.DofPerElem,:]
    aux0        = mesh.nVerts+mesh.nEdges*(PSpace.order-1)
    aux1        = PSpace.DofPerElem-3*PSpace.order
    Lag_PKD     = Lag_PKDCoefficients(mesh.order, xim_mesh)
    Lm          = PKD_Nm(xim_dof, mesh.order)*Lag_PKD
    @inbounds for iElem=1:mesh.nElems
        x1          = view(mesh.ElemsCoords,iElem,1:2:2*mesh.NodesPerElem-1)
        x2          = view(mesh.ElemsCoords,iElem,2:2:2*mesh.NodesPerElem)
        for iDof=1:aux1
            dof                                         = aux0 + (iElem-1)*aux1 + iDof
            for iNode=1:mesh.NodesPerElem
                PSpace.NodesCoords[dof,1]               += x1[iNode]*Lm[iDof,iNode]
                PSpace.NodesCoords[dof,2]               += x2[iNode]*Lm[iDof,iNode]
            end
            PSpace.ElemsDof[iElem,3*PSpace.order+iDof]  = dof 
        end
        
    end
    
    return PSpace
    
end

function NCompute(FES::TrPSpace, xim::AbstractMatrix{Float64})  #FES=finite element space

    return PKD_Nm(xim, FES.order)*FES.Lag_PKD #L_k(x_i) = P_j(x_i)*a_jk

end

#Compute fluctuaction w.r.t. interpolation on FES0:
function flucNCompute(FES::TrFES, xim::AbstractMatrix{Float64},
    FES0::TrPSpace)

    NP      = NCompute(FES, xim)
    
    #Interpolate shape functions of FES onto P_order. For that purpose, 
    #we write the first (phi_i) as 
    #   phi_k(x) = L_j(x) * a_jk.
    #The coefficients a_jk are the interpolation of phi_k at the Chebyshev
    #nodes of FES0:
    a_jk    = NCompute(FES, FES0.xi)
    intNP   = NCompute(FES0, xim)*a_jk
    
    return NP.-intNP
    
end

function gradNCompute(FES::TrPSpace, xim::AbstractMatrix{Float64})  #FES=finite element space

    gradN_PKD   = PKD_gradNm(xim, FES.order)
    return [gradN_PKD[1]*FES.Lag_PKD, gradN_PKD[2]*FES.Lag_PKD]
    
end

#Compute gradients fluctuaction w.r.t. interpolation on FES0:
function gradflucNCompute(FES::TrFES, xim::AbstractMatrix{Float64},
    FES0::TrPSpace)

    gradNP      = gradNCompute(FES, xim)
    a_jk        = NCompute(FES, FES0.xi)
    gradNP0     = gradNCompute(FES0, xim)
    
    return [ gradNP[1].-gradNP0[1]*a_jk, gradNP[2].-gradNP0[2]*a_jk ]
    
end

#-----------------------------------------------------------------------------
#Discontinuous space with PKD shape functions:

#Create QSpace of given order associated with given mesh:
function TrQSpace(mesh::TrMesh, order::Int)

    QSpace              = TrQSpace()
    QSpace.mesh         = mesh
    QSpace.order        = order
    QSpace.DofPerElem   = div((order+1)*(order+2),2)
    QSpace.ElemsDof     = zeros(Int, mesh.nElems, QSpace.DofPerElem)
    for iElem=1:mesh.nElems, iDof=1:QSpace.DofPerElem
        QSpace.ElemsDof[iElem,iDof]     = (iElem-1)*QSpace.DofPerElem + iDof
    end
    QSpace.nDof         = mesh.nElems*QSpace.DofPerElem
    QSpace.nDofU        = QSpace.nDof
    return QSpace
    
end

function NCompute(FES::TrQSpace, xim::AbstractMatrix{Float64})  #FES=finite element space

    return PKD_Nm(xim, FES.order)

end

function gradNCompute(FES::TrQSpace, xim::AbstractMatrix{Float64})  #FES=finite element space

    return PKD_gradNm(xim, FES.order)
    
end

#-----------------------------------------------------------------------------
#Bubble space:

#Elementary bubble function:
function TrElem_Bubble_Nm(xim::AbstractMatrix{Float64})

    #Note that the edges of the equilateral
    #triangle are given by
    #   xi2/(-0.5)      = 1
    #   xi1/d + xi2     = 1
    #   -xi1/d + xi2    = 1
    #with (d,0) and (-d,0) two points of the triangle, given by
    #d/1 = sqrt(3)/2 / (1+1/2)
    dinv            = 3.0/sqrt(3)
    sq32            = sqrt(3.0)/2.0
    bv              = @mlv    -1.0*(-2.0*xim[:,2]-1.0)*
                        (xim[:,1]*dinv+xim[:,2]-1.0)*
                        (-xim[:,1]*dinv+xim[:,2]-1.0)

    return bv

end
function TrElem_Bubble_gradNm(xim::AbstractMatrix{Float64})

    #Note that the edges of the equilateral
    #triangle are given by
    #   xi2/(-0.5)      = 1
    #   xi1/d + xi2     = 1
    #   -xi1/d + xi2    = 1
    #with (d,0) and (-d,0) two points of the triangle, given by
    #d/1 = sqrt(3)/2 / (1+1/2)
    dinv            = 3.0/sqrt(3)
    sq32            = sqrt(3.0)/2.0
    bv              = @mlv    -1.0*(-2.0*xim[:,2]-1.0)*
                        (xim[:,1]*dinv+xim[:,2]-1.0)*
                        (-xim[:,1]*dinv+xim[:,2]-1.0)
    grad1bv         = @mlv    -1.0*(-2.0*xim[:,2]-1.0)*
                        (dinv*(-xim[:,1]*dinv+xim[:,2]-1.0) - 
                         (xim[:,1]*dinv+xim[:,2]-1.0)*dinv)
    grad2bv         = @mlv    -1.0*(
                        (-2.0)*
                        (xim[:,1]*dinv+xim[:,2]-1.0)*
                        (-xim[:,1]*dinv+xim[:,2]-1.0) +
                        (-2.0*xim[:,2]-1.0)*
                        1.0*
                        (-xim[:,1]*dinv+xim[:,2]-1.0) +
                        (-2.0*xim[:,2]-1.0)*
                        (xim[:,1]*dinv+xim[:,2]-1.0)*
                        1.0)
                        
    return [grad1bv, grad2bv]

end

#Create BSpace of given order associated with given mesh:
function TrBSpace(mesh::TrMesh, order::Int)

    BSpace              = TrBSpace()
    BSpace.mesh         = mesh
    BSpace.order        = order
    if order==0
        error("Unable to construct bubble space of order 0")
    elseif order==1
        BSpace.DofPerElem   = 1
    else
        #Multiply "b" only by PKD functions of degree p-2 and p-1:
        BSpace.DofPerElem   = 2*order-1
    end
    BSpace.ElemsDof     = zeros(Int, mesh.nElems, BSpace.DofPerElem)
    for iElem=1:mesh.nElems, iDof=1:BSpace.DofPerElem
        BSpace.ElemsDof[iElem,iDof]     = (iElem-1)*BSpace.DofPerElem + iDof
    end
    BSpace.nDof         = mesh.nElems*BSpace.DofPerElem
    
    #We have to compute the fluctuation w.r.t. PSpace of order r:
    #For that purpose, let us write 
    #   NB_j(x_i)   = NP_k(x_i) a_kj. 
    #By definition, 
    #   a_kj        = NB_j(xi_k), 
    #with xi the Lagrange nodes.
    BSpace.xi           = TrElem_ChebyshevNodes(order)
    BSpace.Lag_PKD      = Lag_PKDCoefficients(order, BSpace.xi)
    BSpace.Bub_Lag      = BX_NCompute(BSpace, BSpace.xi) 
       
    return BSpace
    
end

#Compute bubble functions without fluctuation:
function BX_NCompute(FES::TrBSpace, xim::AbstractMatrix{Float64})  #FES=finite element space

    #PKD shape functions:
    Nm_PKD      = PKD_Nm(xim, FES.order-1)
    indexm      = PKD_index(FES.order-1)
    
    #Retain functions with order max(0,(p-2)) to (p-1):
    p1          = max(0,FES.order-2)
    index1      = indexm[p1+1,1]
    index2      = indexm[1,FES.order-1+1]
    
    #Bubble shape functions:
    bv          = TrElem_Bubble_Nm(xim)
    Nm          = zeros(size(xim,1), FES.DofPerElem)
    for kk=1:FES.DofPerElem
        @mlv    Nm[:,kk]    = Nm_PKD[:,index1-1+kk]*bv
    end
    
    return Nm

end

function BX_gradNCompute(FES::TrBSpace, xim::AbstractMatrix{Float64})  #FES=finite element space

    #PKD shape functions:
    Nm_PKD      = PKD_Nm(xim, FES.order-1)
    gradNm_PKD  = PKD_gradNm(xim, FES.order-1)
    indexm      = PKD_index(FES.order-1)
    
    #Retain functions with order max(0,(p-2)) to (p-1):
    p1          = max(0,FES.order-2)
    index1      = indexm[p1+1,1]
    index2      = indexm[1,FES.order-1+1]
    
    #Bubble shape functions:
    bv          = TrElem_Bubble_Nm(xim)
    gradbv      = TrElem_Bubble_gradNm(xim)
    
    #grad(b*N) = gradb*N + b*gradN:
    grad1Nm     = zeros(size(xim,1), FES.DofPerElem)
    grad2Nm     = zeros(size(xim,1), FES.DofPerElem)
    for kk=1:FES.DofPerElem
        @mlv    grad1Nm[:,kk]   = gradNm_PKD[1][:,index1-1+kk]*bv + 
                                    Nm_PKD[:,index1-1+kk]*gradbv[1]
        @mlv    grad2Nm[:,kk]   = gradNm_PKD[2][:,index1-1+kk]*bv + 
                                    Nm_PKD[:,index1-1+kk]*gradbv[2]
    end
    
    return [ grad1Nm, grad2Nm ]
    
end

#Compute fluctuations. For that purpose, let us write 
#   NB_j(x_i)   = NP_k(x_i) a_kj. 
#By definition, 
#   a_kj        = NB_j(xi_k), 
#with xi the Lagrange nodes.
function NCompute(FES::TrBSpace, xim::AbstractMatrix{Float64})  #FES=finite element space

    #Bubble functions without fluctuations:
    NBm         = BX_NCompute(FES, xim)
    
    #Lagrangian functions:
    NPm         = PKD_Nm(xim, FES.order)*FES.Lag_PKD
    
    return NBm - NPm*FES.Bub_Lag

end

function gradNCompute(FES::TrBSpace, xim::AbstractMatrix{Float64})  #FES=finite element space

    #Bubble functions without fluctuations:
    gradNBm     = BX_gradNCompute(FES, xim)
    
    #Lagrangian functions:
    gradN_PKD   = PKD_gradNm(xim, FES.order)
    gradNPm     = [gradN_PKD[1]*FES.Lag_PKD, gradN_PKD[2]*FES.Lag_PKD]
    
    return [ gradNBm[1]-gradNPm[1]*FES.Bub_Lag, gradNBm[2]-gradNPm[2]*FES.Bub_Lag ]
    
end

#---------------------------------------------------------------------------------
#Polynomial + bubble space:

function TrPBSpace(mesh::TrMesh, order::Int)

    PBSpace             = TrPBSpace()
    PBSpace.mesh        = mesh
    PBSpace.order       = order
    PBSpace.PSpace      = TrPSpace(mesh, order)
    PBSpace.BSpace      = TrBSpace(mesh, order)
    PBSpace.DofPerElem  = PBSpace.PSpace.DofPerElem + PBSpace.BSpace.DofPerElem
    PBSpace.nDof        = PBSpace.PSpace.nDof + PBSpace.BSpace.nDof
    PBSpace.nDofU       = PBSpace.PSpace.nDofU
    PBSpace.ElemsDof    = hcat(PBSpace.PSpace.ElemsDof, PBSpace.BSpace.ElemsDof.+PBSpace.PSpace.nDof)
    PBSpace.dofP        = 1:PBSpace.PSpace.DofPerElem
    PBSpace.dofB        = PBSpace.PSpace.DofPerElem+1:PBSpace.DofPerElem
    return PBSpace
    
end

function NCompute(FES::TrPBSpace, xim::AbstractMatrix{Float64})

    return hcat(NCompute(FES.PSpace, xim), NCompute(FES.BSpace, xim))
    
end

function flucNCompute(FES::TrPBSpace, xim::AbstractMatrix{Float64})

    NP      = zeros(size(xim,1),FES.PSpace.DofPerElem)
    NB      = NCompute(FES.BSpace, xim)
    
    return hcat(NP, NB)
    
end

function gradNCompute(FES::TrPBSpace, xim::AbstractMatrix{Float64})

    gradNP  = gradNCompute(FES.PSpace, xim)
    gradNB  = gradNCompute(FES.BSpace, xim)
    
    return [ hcat(gradNP[1],gradNB[1]), hcat(gradNP[2],gradNB[2]) ]
    
end

function gradflucNCompute(FES::TrPBSpace, xim::AbstractMatrix{Float64})

    gradNP  = [ zeros(size(xim,1),FES.PSpace.DofPerElem), 
                zeros(size(xim,1),FES.PSpace.DofPerElem) ]
    gradNB  = gradNCompute(FES.BSpace, xim)
    
    return [ hcat(gradNP[1],gradNB[1]), hcat(gradNP[2],gradNB[2]) ]
    
end

#---------------------------------------------------------------------------------
#Quadrature rules:

function WVRule(order::Int)

    #Order and number of quadrature nodes needed for Witherden--Vincent quadratures:
    order_nqp   = [ 1 1
                    2 3
                    4 6
                    5 7
                    6 12
                    7 15
                    8 16
                    9 19
                    10 25
                    11 28
                    12 33
                    13 37
                    14 42
                    15 49
                    16 55
                    17 60
                    18 67
                    19 73
                    20 79
                    22 97
                    24 112
                    26 130
                    28 153
                    30 172
                    32 198 
                    34 225
                    36 246 
                    38 273 
                    40 300 ]

    #Get number of quadrature points:
    aux     = findfirst(view(order_nqp,:,1).>=order)
    if isnothing(aux)
        error("Rule of order at least $(order) not found")
    end
    order   = order_nqp[aux,1] #overwrite target order with real one
    nqp     = order_nqp[aux,2]
    
    #Load rule:
    QRule       = TrQRule()
    QRule.order = order
    QRule.nqp   = nqp
    auxm        = zeros(nqp,3)
    auxm        .= readdlm("$(@__DIR__)/WVRules-triangle/$(order)-$(nqp).txt")
    QRule.xi    = TrElem_Iso2_to_Equi(auxm[:,1:2])
    QRule.w     = auxm[:,3]*TrElem_Area/2.0
    
    return QRule
    
end

#Witherden-Vincent rule with N*N subtriangles:
function WVRule(order::Int, N::Int)

    QRule1          = WVRule(order)
    
    QRule           = TrQRule()
    QRule.order     = QRule1.order
    QRule.nqp       = QRule1.nqp * N*N
    
    #Segment at y0 goes from -x0 to x0:
    x0fun(y0)       = (1.0-y0)/sqrt(3)
    
    #Loop subtriangles:
    QRule.xi        = zeros(QRule.nqp, 2)
    QRule.w         = zeros(QRule.nqp)
    kk              = 0
    for jj=1:N
    
        #Top and bottom segments:
        y_top   = 1.0 - (jj-1)*1.5/N;
        y_bot   = 1.0 - jj*1.5/N;
        x_top   = x0fun(y_top)
        x_bot   = x0fun(y_bot)
        
        #Upside triangles:
        for ii=1:jj
            
            #Divide x_bot in jj+1 nodes and get node ii+1/2:
            xbar    = -x_bot + 2*x_bot/jj*(ii-0.5)
            #
            ybar    = (2*y_bot+y_top)/3
            
            #Expand QRule:
            kk                      += 1
            kk_index                = (kk-1)*QRule1.nqp+1:kk*QRule1.nqp
            @. QRule.xi[kk_index,1] = $view(QRule1.xi,:,1)/N + xbar
            @. QRule.xi[kk_index,2] = $view(QRule1.xi,:,2)/N + ybar
            QRule.w[kk_index]       = QRule1.w/(N*N)
            
        end
        
        #Downside triangles:
        for ii=2:jj
        
            #Divide x_bot in jj+1 nodes and get nodes 2:jj:
            xbar    = -x_bot + 2*x_bot/jj*(ii-1)
            #
            ybar    = (y_bot+2*y_top)/3
            
            #Expand QRule:
            kk                      += 1
            kk_index                = (kk-1)*QRule1.nqp+1:kk*QRule1.nqp
            @. QRule.xi[kk_index,1] = $view(QRule1.xi,:,1)/N + xbar
            @. QRule.xi[kk_index,2] = -$view(QRule1.xi,:,2)/N + ybar
            QRule.w[kk_index]       = QRule1.w/(N*N)
            
        end
        
    end
    
    return QRule
    
end

function GLRule(order::Int)

    #order=2*nqp-1, so
    nqp         = Int(ceil((order+1)/2))
    
    #Load rule:
    QRule       = LineQRule()
    QRule.order = 2*nqp-1
    QRule.nqp   = nqp
    aux         = gausslegendre(nqp)
    QRule.xi    = aux[1]
    QRule.w     = aux[2]
    
    return QRule
    
end
