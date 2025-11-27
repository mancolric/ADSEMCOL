#READ MESH:
#We only read linear meshes. We have only structures ElemsNodes (defined for mesh and bmesh)
#and NodesCoords.
#
#Then, the meshes are post processed to create inner nodes if necessary.
function TrMesh_CreateRegular(FesOrder::Int, x1::Float64, x2::Float64, Nx::Int, 
    y1::Float64, y2::Float64, Ny::Int)

    #Allocate outputs:
    nBoundaries             = 4
    mesh                    = TrMesh()
    boundary_meshes         = Vector{TrBmesh}(undef, nBoundaries)
    
    #Set up some fields:
    mesh.dim                = 2
    mesh.VertsPerElem       = 3
    mesh.NodesPerElem       = (FesOrder+1)*(FesOrder+2)÷2
    mesh.order              = FesOrder
    mesh.nBounds            = nBoundaries
    mesh.bmesh              = boundary_meshes
    for ii=1:nBoundaries
        boundary_meshes[ii]                 = TrBmesh()
        boundary_meshes[ii].dim             = 2
        boundary_meshes[ii].VertsPerElem    = 2
        boundary_meshes[ii].NodesPerElem    = FesOrder+1
        boundary_meshes[ii].order           = FesOrder
        boundary_meshes[ii].nVerts          = 0
        boundary_meshes[ii].nNodes          = 0
        boundary_meshes[ii].nElems          = 0
        boundary_meshes[ii].ElemsNodes      = zeros(0,2)
    end

    #--------------------------------------------------------------
    #CREATE LINEAR MESH:

    #Mesh nodes:
    xv                  = linspace(x1,x2,Nx)
    yv                  = linspace(y1,y2,Ny)
    xm,ym               = ndgrid(xv,yv)
    mesh.NodesCoords    = [xm[:] ym[:]]
    mesh.nVerts         = Nx*Ny
    mesh.nNodes         = Nx*Ny

    #Inner elements:
    mesh.nElems         = 2*(Nx-1)*(Ny-1)
    mesh.ElemsNodes     = zeros(Int,mesh.nElems,3)
    iElem               = 1
    for jj=1:Ny-1, ii=1:Nx-1
        mesh.ElemsNodes[iElem,1]    = ii+(jj-1)*Nx
        mesh.ElemsNodes[iElem,2]    = ii+1+(jj-1)*Nx
        mesh.ElemsNodes[iElem,3]    = ii+jj*Nx
        mesh.ElemsNodes[iElem+1,1]  = ii+1+(jj-1)*Nx
        mesh.ElemsNodes[iElem+1,2]  = ii+1+jj*Nx
        mesh.ElemsNodes[iElem+1,3]  = ii+jj*Nx
        iElem                       += 2
    end

    #Boundaries:
    boundary_meshes[1].nElems           = Nx-1
    boundary_meshes[3].nElems           = Nx-1
    boundary_meshes[1].ElemsNodes       = zeros(Nx-1,2)
    boundary_meshes[3].ElemsNodes       = zeros(Nx-1,2)
    for ii=1:Nx-1, jj=1:2, 
        boundary_meshes[1].ElemsNodes[ii,jj]    = ii+jj-1
        boundary_meshes[3].ElemsNodes[ii,jj]    = (Nx-ii-jj+2) + (Ny-1)*Nx
    end
    boundary_meshes[2].nElems           = Ny-1
    boundary_meshes[4].nElems           = Ny-1
    boundary_meshes[2].ElemsNodes       = zeros(Ny-1,2)
    boundary_meshes[4].ElemsNodes       = zeros(Ny-1,2)
    for ii=1:Ny-1, jj=1:2
        boundary_meshes[2].ElemsNodes[ii,jj]    = Nx + (ii+jj-2)*Nx
        boundary_meshes[4].ElemsNodes[ii,jj]    = 1 + (Ny-ii-jj+1)*Nx
    end

    TrMesh_Process!(mesh)
    
    #Return:
    return mesh

end

#Create mesh with BAMG for rectangular domain:
function TrMesh_Rectangle_Create!(MeshFile::String,
    x1::Float64, x2::Float64, Nx::Int, 
    y1::Float64, y2::Float64, Ny::Int)
    
    #Check file extension:
    MeshFile_length     = length(MeshFile)
    if MeshFile[MeshFile_length-3:MeshFile_length]!=".geo"
        error("MeshFile extension must be .geo")
    end 
    
    #Check vertices:
    if Nx<2 || Ny<2
        error("Number of vertices at each edge must be >1")
    end     
    
    #Values of x and y at the boundary:
    M 			= 2
    N 			= 2
    xv 			= linspace(x1,x2,M)
    yv 			= linspace(y1,y2,N)
    hx 			= (x2-x1)/(Nx-1)
    hy 			= (y2-y1)/(Ny-1)
    nVerts  	= 2*(M-1)+2*(N-1)
    hv 			= zeros(nVerts)
    hmin 		= 0.85*min(hx,hy)
    hmax 		= 1.17*max(hx,hy)

    #Create auxiliary file with Bamg input:
    file_ptr 	= open(MeshFile, "w")
    write(file_ptr,"MeshVersionFormatted 0\n")
    write(file_ptr,"Dimension 2\n")
    write(file_ptr,"Vertices $(nVerts)\n")
    
    #Write vertices:
    iVert 		= 0
    for ii=1:M-1
        iVert 		+= 1
        write(file_ptr, "$(xv[ii]) $(y1) $iVert\n")
        hv[iVert]   = hx
    end
    for jj=1:N-1
        iVert 		+= 1
        write(file_ptr, "$(x2) $(yv[jj]) $iVert\n")
        hv[iVert] 	= hy
    end
    for ii=M:-1:2
        iVert 		+= 1
        write(file_ptr, "$(xv[ii]) $(y2) $iVert\n")
        hv[iVert] 	= hx
    end
    for jj=N:-1:2
        iVert 		+= 1
        write(file_ptr, "$(x1) $(yv[jj]) $iVert\n")
        hv[iVert] 	= hy
    end
    hhv 	= 0.5*vcat( hv[nVerts]+hv[1], hv[1:nVerts-1]+hv[2:nVerts])

    #Write edges:
    write(file_ptr, "Edges $(nVerts)\n")
    iVert 			= 0
    for ii=1:M-1
        iVert 		+= 1
        write(file_ptr, "$(iVert) $(iVert+1) 1\n")
    end
    for jj=1:N-1
        iVert 		+= 1
        write(file_ptr, "$(iVert) $(iVert+1) 2\n")
    end
    for ii=M:-1:2
        iVert 		+= 1
        write(file_ptr, "$(iVert) $(iVert+1) 3\n")
    end
    for jj=N:-1:3
        iVert 		+= 1
        write(file_ptr, "$(iVert) $(iVert+1) 4\n")
    end
    write(file_ptr, "$(nVerts) 1 4\n")

    #Set target h for vertices
    write(file_ptr, "hVertices\n")
    for hh in hhv
        write(file_ptr, "$hh ")
    end
    write(file_ptr,"\n")
    
    close(file_ptr)
    
    return
    
end 
		
		
function TrMesh_Read(file::String, FesOrder::Int, nBoundaries::Int)

    file_extension      = file[findlast('.', file)+1:length(file)]
    if file_extension=="db"
        return TrMesh_ReadFromBamg(file, FesOrder, nBoundaries)
    else
        error("Mesh file extension $(file_extension) not supported")
    end

end

"
    mesh, boundary_meshes = TrXN_ReadFromBamg!(file,FesOrder,nBoundaries)

Read mesh from output Bamg file (version 2.2).
"
function TrMesh_ReadFromBamg(file::String, FesOrder::Int, nBoundaries::Int)

    #Allocate outputs:
    mesh                = TrMesh()
    boundary_meshes     = Vector{TrBmesh}(undef,nBoundaries)
    mesh.bmesh          = boundary_meshes
    mesh.nBounds        = nBoundaries
    
    #Set up some fields:
    mesh.dim                = 2
    mesh.VertsPerElem       = 3
    mesh.NodesPerElem       = (FesOrder+1)*(FesOrder+2)÷2
    mesh.order              = FesOrder
    for ii=1:nBoundaries
        boundary_meshes[ii]                 = TrBmesh()
        boundary_meshes[ii].dim             = 2
        boundary_meshes[ii].VertsPerElem    = 2
        boundary_meshes[ii].NodesPerElem    = FesOrder+1
        boundary_meshes[ii].order           = FesOrder
        boundary_meshes[ii].nVerts          = 0
        boundary_meshes[ii].nNodes          = 0
        boundary_meshes[ii].nElems          = 0
        boundary_meshes[ii].ElemsNodes      = zeros(0,2)
    end

    #Open file:
    file_ptr            = open(file,"r")

    #Get nodes:
    readuntil(file_ptr, "Vertices")
    readline(file_ptr)
    mesh.nVerts         = parse(Int,readline(file_ptr))
    mesh.nNodes         = mesh.nVerts
    mesh.NodesCoords    = zeros(mesh.nNodes,mesh.dim)
    for ii=1:mesh.nNodes
        for jj=1:2
            mesh.NodesCoords[ii,jj]         = parse(Float64,readuntil(file_ptr," "))
        end
        readuntil(file_ptr,"\n")
    end

    #Get boundary elements:
    readuntil(file_ptr,"Edges")
    readline(file_ptr)
    nBElems     = parse(Int,readline(file_ptr))
    for ii=1:nBElems
        v1      = parse(Int,readuntil(file_ptr," "))
        v2      = parse(Int,readuntil(file_ptr," "))
        iBound  = parse(Int,readuntil(file_ptr,"\n"))
        boundary_meshes[iBound].ElemsNodes  = vcat(boundary_meshes[iBound].ElemsNodes, [v1 v2])
    end
    for ii=1:nBoundaries
        boundary_meshes[ii].nElems          = size(boundary_meshes[ii].ElemsNodes,1)
    end
    
    #Get elements:
    readuntil(file_ptr,"Triangles")
    readline(file_ptr)
    mesh.nElems         = parse(Int,readline(file_ptr))
    mesh.ElemsNodes     = zeros(Int,mesh.nElems,mesh.NodesPerElem)
    for ii=1:mesh.nElems
        for jj=1:3
            mesh.ElemsNodes[ii,jj]          = parse(Float64,readuntil(file_ptr," "))
        end
        readuntil(file_ptr,"\n")
    end

    #Close msh file:
    close(file_ptr)

    TrMesh_Process!(mesh)
    
    #Return:
    return mesh

end

function TrMesh_Process!(mesh::TrMesh)

    #Create ElemsNodesSp with current (linear) nodes:
    imat                = zeros(mesh.nElems, mesh.VertsPerElem)
    for iElem=1:mesh.nElems, iNode=1:mesh.VertsPerElem
        imat[iElem,iNode]   = iElem
    end
    smat                = ones(mesh.nElems,mesh.VertsPerElem)
    mesh.ElemsNodesSp   = sparse(imat[:],mesh.ElemsNodes[:,1:mesh.VertsPerElem][:],smat[:])
    
    #Loop elements and local edges:
    node_start              = [1,2,3]
    node_end                = [2,3,1]
    mesh.nEdges             = 0
    mesh.ElemsNeigh         = zeros(Int, mesh.nElems, mesh.VertsPerElem)
    mesh.ElemsNeighFaces    = zeros(Int, mesh.nElems, mesh.VertsPerElem)
    mesh.ElemsEdges         = zeros(Int, mesh.nElems, mesh.VertsPerElem)
    mesh.EdgesNodes         = zeros(Int, mesh.nElems*mesh.VertsPerElem, mesh.order+1)
    mesh.EdgesElems         = zeros(Int, mesh.nElems*mesh.VertsPerElem, 2)
    mesh.EdgesFaces         = zeros(Int, mesh.nElems*mesh.VertsPerElem, 2)
    @inbounds for iElem=1:mesh.nElems, iEdge=1:mesh.VertsPerElem
        
        #Has neighbour already been detected?
        if mesh.ElemsNeigh[iElem,iEdge]==0  #No
            
            #Select nodes:
            node1       = mesh.ElemsNodes[iElem,node_start[iEdge]]
            node2       = mesh.ElemsNodes[iElem,node_end[iEdge]]
        
            #Save edge:
            mesh.nEdges                     += 1
            mesh.ElemsEdges[iElem,iEdge]    = mesh.nEdges
            mesh.EdgesNodes[mesh.nEdges,1]  = node1
            mesh.EdgesNodes[mesh.nEdges,2]  = node2
            mesh.EdgesElems[mesh.nEdges,1]  = iElem
            mesh.EdgesFaces[mesh.nEdges,1]  = iEdge
            
            #Elements common to each node:
            elems1      = mesh.ElemsNodesSp.rowval[nzrange(mesh.ElemsNodesSp, node1)]
            elems2      = mesh.ElemsNodesSp.rowval[nzrange(mesh.ElemsNodesSp, node2)]
            
            #Disable current element:
            elems1[elems1.==iElem]  .= -1.0
        
            #Find common element:
            CommonElem  = intersect(elems1,elems2)
            if !isempty(CommonElem)
            
                #Save data of neighbour element for current iElem:
                ElemNeigh                           = CommonElem[1]
                mesh.ElemsNeigh[iElem,iEdge]        = ElemNeigh
                NeighFace                           = findfirst(mesh.ElemsNodes[ElemNeigh,:].==node2)
                mesh.ElemsNeighFaces[iElem,iEdge]   = NeighFace
            
                #Save data of current iElem for neighbour element:
                mesh.ElemsNeigh[ElemNeigh,NeighFace]        = iElem
                mesh.ElemsNeighFaces[ElemNeigh,NeighFace]   = iEdge
                
                #Save data of current edge for neighbour element:
                mesh.EdgesElems[mesh.nEdges,2]      = ElemNeigh
                mesh.EdgesFaces[mesh.nEdges,2]      = NeighFace
                
            else    
                #if there is no common element, current (iElem,iEdge) is a boundary edge. Do nothing
            end
            
        else
        
            #The current (iElem, iEdge) has already been processed.
            
        end
    
    end
    #Chunk matrices:
    mesh.EdgesNodes     = mesh.EdgesNodes[1:mesh.nEdges,:]
    mesh.EdgesElems     = mesh.EdgesElems[1:mesh.nEdges,:]
    mesh.EdgesFaces     = mesh.EdgesFaces[1:mesh.nEdges,:]
    
    #Expand ElemsNodes and NodesCoords (EdgesNodes is already expanded):
    mesh.NodesPerElem   = div((mesh.order+1)*(mesh.order+2),2)
    mesh.nNodes         = mesh.nVerts + (mesh.order-1)*mesh.nEdges + 
                            mesh.nElems*(mesh.NodesPerElem-mesh.VertsPerElem*mesh.order)
    mesh.ElemsNodes     = hcat(mesh.ElemsNodes[:,1:mesh.VertsPerElem],
                            zeros(mesh.nElems, mesh.NodesPerElem-mesh.VertsPerElem))
    mesh.NodesCoords    = vcat(mesh.NodesCoords[1:mesh.nVerts,:], 
                            zeros(mesh.nNodes-mesh.nVerts, mesh.dim))
                            
    #Create high-order nodes at edges:
    xiv         = ChebyNodes(-1.0, 1.0, mesh.order+1)[2:mesh.order]
    L1v         = @mlv 0.5*(1.0-xiv)
    L2v         = @mlv 0.5*(1.0+xiv)
    @inbounds for iEdge=1:mesh.nEdges, iNode=1:mesh.order-1
        node                            = mesh.nVerts + (iEdge-1)*(mesh.order-1) + iNode
        node1                           = mesh.EdgesNodes[iEdge,1]
        node2                           = mesh.EdgesNodes[iEdge,2]
        mesh.NodesCoords[node,1]        = L1v[iNode]*mesh.NodesCoords[node1,1]+L2v[iNode]*mesh.NodesCoords[node2,1]
        mesh.NodesCoords[node,2]        = L1v[iNode]*mesh.NodesCoords[node1,2]+L2v[iNode]*mesh.NodesCoords[node2,2]
        mesh.EdgesNodes[iEdge,2+iNode]  = node
        elem1                           = mesh.EdgesElems[iEdge,1]
        face1                           = mesh.EdgesFaces[iEdge,1]
        mesh.ElemsNodes[elem1, mesh.VertsPerElem+
            (face1-1)*(mesh.order-1)+iNode]         = node
        elem2                           = mesh.EdgesElems[iEdge,2]
        if elem2>0
            face2                       = mesh.EdgesFaces[iEdge,2]
            mesh.ElemsNodes[elem2, mesh.VertsPerElem+
                face2*(mesh.order-1)+1-iNode]       = node
        end
    end
    
    #Create high-order nodes inside elements:
    xim         = TrElem_ChebyshevNodes(mesh.order)[mesh.VertsPerElem*mesh.order+1:mesh.NodesPerElem,:]
    aux0        = mesh.nVerts+mesh.nEdges*(mesh.order-1)
    aux1        = mesh.NodesPerElem-3*mesh.order
    Lm          = TrElem_LinearLag(xim)
    @inbounds for iElem=1:mesh.nElems
        x1          = mesh.NodesCoords[mesh.ElemsNodes[iElem,1],:]
        x2          = mesh.NodesCoords[mesh.ElemsNodes[iElem,2],:]    
        x3          = mesh.NodesCoords[mesh.ElemsNodes[iElem,3],:]
        for iNode=1:aux1
            node    = aux0 + (iElem-1)*aux1 + iNode
            mesh.NodesCoords[node,:]                    = x1*Lm[iNode,1] + x2*Lm[iNode,2] + x3*Lm[iNode,3]
            mesh.ElemsNodes[iElem,3*mesh.order+iNode]   = node 
        end
        
    end
    
    #Compute ElemsNodesSp (again, now with high-order nodes):
    imat                = zeros(mesh.nElems, mesh.NodesPerElem)
    for iElem=1:mesh.nElems, iNode=1:mesh.NodesPerElem
        imat[iElem,iNode]   = iElem
    end
    smat                = ones(mesh.nElems,mesh.NodesPerElem)
    mesh.ElemsNodesSp   = sparse(imat[:],mesh.ElemsNodes[:],smat[:])
    
    #Compute ElemsCoords:
    ElemsCoordsCompute!(mesh)
    
    #Variables to compute isoparametric transformations:
    mesh.xi         = TrElem_ChebyshevNodes(mesh.order)
    mesh.Lag_PKD    = Lag_PKDCoefficients(mesh.order, mesh.xi)
    
    #--------------------------------------------------------------
    #Boundary meshes:
    
    #NOTE: The fields 
    #     dim         ::Int
    #     nElems      ::Int
    #     nNodes      ::Int
    #     nVerts      ::Int
    #     NodesPerElem::Int
    #     VertsPerElem::Int
    #     order       ::Int
    #     ElemsNodes  ::Matrix{Int}
    #have already been defined 
    
    #Generate new bmeshes with same ParentFace for each element:
    bmeshes     = Vector{TrBmesh}(undef, 3*mesh.nBounds)
    
    #Split original meshes into child meshes with same ParentFace for each element:
    for ib=1:mesh.nBounds
    
        #Get current bmesh:
        bmesh                           = mesh.bmesh[ib]
        
        #Loop bmesh elements. Find parent info:
        bmesh.ParentElems   = zeros(Int, bmesh.nElems)
        bmesh_ParentFaces   = zeros(Int, bmesh.nElems)
        bmesh.ParentEdges   = zeros(Int, bmesh.nElems)
        for iElem=1:bmesh.nElems
            
            #Find common element:
            node1                       = bmesh.ElemsNodes[iElem,1]
            node2                       = bmesh.ElemsNodes[iElem,2]
            elems1                      = mesh.ElemsNodesSp.rowval[nzrange(mesh.ElemsNodesSp, node1)]
            elems2                      = mesh.ElemsNodesSp.rowval[nzrange(mesh.ElemsNodesSp, node2)]
            elem                        = intersect(elems1,elems2)[1]
            bmesh.ParentElems[iElem]    = elem
            
            #Find face and edge:
            face                        = findfirst(view(mesh.ElemsNodes,elem,:).==node1)
            bmesh_ParentFaces[iElem]    = face
            bmesh.ParentEdges[iElem]    = mesh.ElemsEdges[elem,face]
            
        end
        
        #Loop faces:
        for face=1:3
            
            #Create mesh:
            jb                          = (ib-1)*3+face
            bmeshes[jb]                 = TrBmesh()
            bmeshes[jb].boundary_id     = ib
            bmeshes[jb].dim             = 2
            
            #Find elements:
            aux                         = bmesh_ParentFaces.==face
            bmeshes[jb].nElems          = sum(aux)
            bmeshes[jb].nNodes          = 0 #will be updated later
            bmeshes[jb].nVerts          = 0 #will be updated later
            bmeshes[jb].NodesPerElem    = bmesh.order+1
            bmeshes[jb].VertsPerElem    = 2
            bmeshes[jb].order           = bmesh.order
            bmeshes[jb].ElemsNodes      = bmesh.ElemsNodes[aux,:]
            
            #Save parent face, elements and edges:
            bmeshes[jb].ParentFace      = face
            bmeshes[jb].ParentElems     = bmesh.ParentElems[aux]
            bmeshes[jb].ParentEdges     = bmesh.ParentEdges[aux]
            
        end 
        
    end 
    
    #Get rest of fields for boundaries:
    mesh.bmesh                  = bmeshes
#     figure()
#     PlotMesh!(mesh, color="k")
    for ib=1:mesh.nBounds*3
    
        #Get mesh:
        bmesh                   = mesh.bmesh[ib]
        
        #Re-write ElemsNodes with high-order nodes:
        bmesh.ElemsNodes        = mesh.EdgesNodes[bmesh.ParentEdges,:]
        
        #Renumber nodes:
        bmesh.ParentNodes       = unique(bmesh.ElemsNodes)
        bmesh.nNodes            = length(bmesh.ParentNodes)
        ChildNodes              = sparsevec(bmesh.ParentNodes, 1:bmesh.nNodes, mesh.nNodes)
        bmesh.ElemsNodes        = ChildNodes[bmesh.ElemsNodes]
        bmesh.NodesCoords       = mesh.NodesCoords[bmesh.ParentNodes,:]
        if bmesh.nElems>0
            bmesh.nVerts        = maximum(view(bmesh.ElemsNodes,:,1:2))
        end 
        bmesh.nEdges            = bmesh.nElems
        ElemsCoordsCompute!(bmesh)
       
        #Variables to compute isoparametric transformations:
        xi_bound        = ChebyNodes(-1.0, 1.0, mesh.order+1)
        xi_bound        = vcat(-1.0, 1.0, xi_bound[2:mesh.order]) #Careful: Lagrangian nodes are defined in this order
        Lag_Leg         = Lag_LegCoefficients(mesh.order, xi_bound)
        for ii=1:mesh.nBounds
            bmesh.xi        = xi_bound
            bmesh.Lag_Leg   = Lag_Leg
        end
    
#         colorv  = ["b", "g", "r"]
#         PlotMesh!(bmesh, color=colorv[mod(ib,3)+1], linewidth=2.0)
        
    end

    #Metric tensor:
    mesh.metric     = ComputeMetric(mesh)
    
    return mesh
    
end

function ElemsCoordsCompute!(mesh::TrUMesh)

    mesh.ElemsCoords            = zeros(mesh.nElems,mesh.dim*mesh.NodesPerElem)
    for ii=1:mesh.NodesPerElem
        mesh.ElemsCoords[:,mesh.dim*(ii-1)+1:mesh.dim*ii]   =
                                mesh.NodesCoords[mesh.ElemsNodes[:,ii],:]
    end

end

#Compute parent nodes and faces of boundary mesh:
function ParentElemsFaces!(bmesh::TrBmesh, mesh::TrMesh)

    bmesh.ParentElems   = zeros(Int, bmesh.nElems)
    bmesh.ParentFaces   = zeros(Int, bmesh.nElems)
    for iElem=1:bmesh.nElems
        #Find common element:
        node1                       = bmesh.ParentNodes[bmesh.ElemsNodes[iElem,1]]
        node2                       = bmesh.ParentNodes[bmesh.ElemsNodes[iElem,2]]
        elems1                      = mesh.ElemsNodesSp.rowval[nzrange(mesh.ElemsNodesSp, node1)]
        elems2                      = mesh.ElemsNodesSp.rowval[nzrange(mesh.ElemsNodesSp, node2)]
        elem                        = intersect(elems1,elems2)[1]
        bmesh.ParentElems[iElem]    = elem
        #Find face:
        bmesh.ParentFaces[iElem]    = findfirst(view(mesh.ElemsNodes,elem,:).==node1)
    end
    
    return
    
end

function PlotMesh!(mesh::TrMesh; elems::AbstractVector{Int}=1:mesh.nElems, 
    color::String="b", linewidth::Float64=0.5, symmetry::String="none")

    triplot(mesh.NodesCoords[1:mesh.nVerts,1], mesh.NodesCoords[1:mesh.nVerts,2], 
        triangles=mesh.ElemsNodes[elems,1:3].-1, 
        color=color, linewidth=linewidth)
    if symmetry=="x" || symmetry=="xy"
        triplot(-mesh.NodesCoords[1:mesh.nVerts,1], mesh.NodesCoords[1:mesh.nVerts,2], 
            triangles=mesh.ElemsNodes[elems,1:3].-1, 
            color=color, linewidth=linewidth)
    end 
    if symmetry=="y" || symmetry=="xy"
        triplot(mesh.NodesCoords[1:mesh.nVerts,1], -mesh.NodesCoords[1:mesh.nVerts,2], 
            triangles=mesh.ElemsNodes[elems,1:3].-1, 
            color=color, linewidth=linewidth)
    end
    if symmetry=="xy"
        triplot(-mesh.NodesCoords[1:mesh.nVerts,1], -mesh.NodesCoords[1:mesh.nVerts,2], 
            triangles=mesh.ElemsNodes[elems,1:3].-1, 
            color=color, linewidth=linewidth)
    end     
    
    return
    
end

function PlotMesh!(mesh::TrBmesh; color::String="b", linewidth::Float64=0.5)

#     plot(mesh.NodesCoords[1:mesh.nVerts,1], mesh.NodesCoords[1:mesh.nVerts,2], 
#         color=color, linewidth=linewidth)
    xcoords     = mesh.ElemsCoords[:,[1,3]]
    ycoords     = mesh.ElemsCoords[:,[2,4]]
    plot(transpose(xcoords), transpose(ycoords), color=color, linewidth=linewidth)
    return
    
end

function PlotNodes!(mesh::TrUMesh; nodes::AbstractVector{Int}=1:mesh.nNodes, 
    color::String="k", markersize::Float64=0.5)

    plot(mesh.NodesCoords[nodes,1], mesh.NodesCoords[nodes,2], linestyle="none", 
        color=color, markersize=markersize, marker=".")
    
    return
    
end

function PlotEllipse!(mesh::TrMesh, elem::Int;
    color::String="r", linewidth::Float64=0.5)

    xbar        = barycenters(mesh, elems=[elem])
    
    #Extract metric tensor:
    V               = zeros(2,2)
    for ii=1:2, jj=1:2
        V[ii,jj]    = mesh.metric.V[ii,jj][elem]
    end
    
    #Compute h=1/sqrt(x^T M x), with x=[cos(theta), sin(theta)]:
    nTheta      = 361
    thetav      = linspace(0.0, 2.0*pi, nTheta)
    x           = [ cos.(thetav), sin.(thetav) ]
    #y_i = (V^T x)_i = V_ji * x_j:
    y           = Vector{Vector{Float64}}(undef,2)
    for ii=1:2
        y[ii]               = zeros(nTheta)
        for jj=1:2
            @mlv y[ii]      += mesh.metric.V[jj,ii][elem] * x[jj]
        end 
    end 
    #z = y_i M_ij y_j:
    z           = @mlv mesh.metric.lambda_bar[elem] * 
                    (y[1]*mesh.metric.sf[1][elem]*y[1] + 
                     y[2]*mesh.metric.sf[2][elem]*y[2])
    
    #h = 1/sqrt(z):
    h           = @mlv 1.0/sqrt(z)
    
    #Plot ellipse:
    @. $plot(xbar[1,1] + h*x[1], xbar[1,2] + h*x[2], color=color, linewidth=linewidth)
    
    return h
    
end 
        
#--------------------------------------------------------
#Lagrangian shape functions associated with the mesh:

function NCompute(mesh::TrMesh, xim::AbstractMatrix{Float64}) 

    return PKD_Nm(xim, mesh.order)*mesh.Lag_PKD #L_k(x_i) = P_j(x_i)*a_jk

end

function gradNCompute(mesh::TrMesh, xim::AbstractMatrix{Float64})  #FES=finite element space

    gradN_PKD   = PKD_gradNm(xim, mesh.order)
    return [gradN_PKD[1]*mesh.Lag_PKD, gradN_PKD[2]*mesh.Lag_PKD]
    
end

function NCompute(mesh::TrBmesh, xiv::AbstractVector{Float64}) 

    return PolyLeg(xiv, mesh.order)*mesh.Lag_Leg #L_k(x_i) = P_j(x_i)*a_jk

end

function gradNCompute(mesh::TrBmesh, xiv::AbstractVector{Float64})  #FES=finite element space

    gradN_Leg       = dPolyLeg(xiv, mesh.order)
    return [gradN_Leg*mesh.Lag_Leg]  #Only one component (dN/dxi)
    
end

#Compute Jinv for given elements and quadrature nodes:
function JinvCompute(mesh::TrMesh, xim::AbstractArray{Float64,2};
    elems::AbstractVector{Bool}=ones(Bool,mesh.nElems))

    #Reshape elems and xim:
    nElems          = sum(elems)
    nqp             = size(xim,1)
    
    #Load gradients:
    gradNm          = gradNCompute(mesh, xim)
    
    #Loop mesh nodes:
    x1Dof           = zeros(nElems,nqp)
    x2Dof           = zeros(nElems,nqp)
    Jinv11          = zeros(nElems,nqp)
    Jinv21          = zeros(nElems,nqp)
    Jinv12          = zeros(nElems,nqp)
    Jinv22          = zeros(nElems,nqp)
    for iNode=1:mesh.NodesPerElem
        x1Dof           = mesh.ElemsCoords[elems, 2*iNode-1]
        x2Dof           = mesh.ElemsCoords[elems, 2*iNode]
        @mlv    Jinv11  += x2Dof*gradNm[2][:,iNode]
        @mlv    Jinv21  -= x2Dof*gradNm[1][:,iNode]
        @mlv    Jinv12  -= x1Dof*gradNm[2][:,iNode]
        @mlv    Jinv22  += x1Dof*gradNm[1][:,iNode]
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

#Compute Jinv for given elements and quadrature nodes:
function JinvCompute(mesh::TrMesh, elems::AbstractMatrix{Int}, xim::AbstractArray{Float64,3})

    #Reshape elems and xim:
    nElems, nqp     = size(elems)
    elems           = reshape(elems, :)
    xim             = reshape(xim, :, 2)
    
    #Load gradients:
    gradNm          = gradNCompute(mesh, xim)
    
    #Loop mesh nodes:
    x1Dof           = zeros(length(elems))
    x2Dof           = zeros(length(elems))
    Jinv11          = zeros(length(elems))
    Jinv21          = zeros(length(elems))
    Jinv12          = zeros(length(elems))
    Jinv22          = zeros(length(elems))
    @inbounds for iNode=1:mesh.NodesPerElem
        x1Dof           = mesh.ElemsCoords[elems, 2*iNode-1]
        x2Dof           = mesh.ElemsCoords[elems, 2*iNode]
        @mlv    Jinv11  += x2Dof*gradNm[2][:,iNode]
        @mlv    Jinv21  -= x2Dof*gradNm[1][:,iNode]
        @mlv    Jinv12  -= x1Dof*gradNm[2][:,iNode]
        @mlv    Jinv22  += x1Dof*gradNm[1][:,iNode]
    end
    detJ            = @mlv    Jinv11*Jinv22 - Jinv21*Jinv12 
    @mlv    Jinv11  /= detJ
    @mlv    Jinv21  /= detJ
    @mlv    Jinv12  /= detJ
    @mlv    Jinv22  /= detJ

    #Pack Jinv:
    Jinv            = Matrix{Matrix{Float64}}(undef,2,2)
    Jinv[1,1]       = reshape(Jinv11, nElems, nqp)
    Jinv[2,1]       = reshape(Jinv21, nElems, nqp)
    Jinv[1,2]       = reshape(Jinv12, nElems, nqp)
    Jinv[2,2]       = reshape(Jinv22, nElems, nqp)
    
    return Jinv
    
end

#Compute coordinates of Lagrange nodes for each element:
function LagNodes(mesh::TrMesh, xi_Lag::Matrix{Float64})
    
    Nmeshm      = NCompute(mesh, xi_Lag)
    x1_Lag      = view(mesh.ElemsCoords,:,1:2:2*mesh.NodesPerElem-1) * transpose(Nmeshm)
    x2_Lag      = view(mesh.ElemsCoords,:,2:2:2*mesh.NodesPerElem) * transpose(Nmeshm)
    x_Lag       = [ x1_Lag, x2_Lag ]
    return x_Lag
    
end 

#--------------------------------------------------------
#Other functions:

import PyPlot.plot

function plot(mesh::TrMesh; color::String="b", linewidth::Float64=0.5, 
    elems::Union{Vector{Bool}, Vector{Int}}=ones(Bool, mesh.nElems))

    PyPlot.triplot(mesh.NodesCoords[1:mesh.nVerts,1], mesh.NodesCoords[1:mesh.nVerts,2], 
        triangles=mesh.ElemsNodes[elems,1:3].-1, color=color, linewidth=linewidth)
        
    return
    
end

#Matrix with barycenters:
function barycenters(mesh::TrMesh; 
    elems::Union{AbstractVector{Bool}, Vector{Int}}=ones(Bool, mesh.nElems))

    Nbar    = NCompute(mesh, [0.0 0.0])
    xbar    = zeros(sum(elems),2)
    for iNode=1:mesh.NodesPerElem
        x1Dof           = mesh.ElemsCoords[elems, 2*iNode-1]
        x2Dof           = mesh.ElemsCoords[elems, 2*iNode]
        @mlv xbar[:,1]  += Nbar[1,iNode]*x1Dof
        @mlv xbar[:,2]  += Nbar[1,iNode]*x2Dof
    end
    
    return xbar
    
end

#Vector with areas:
function areas(mesh::TrMesh)

    return areas(mesh.metric)
    
end
function areas(metric::SpectralDecomp{Vector{Float64},Int})

    return @mlv TrElem_Area/metric.lambda_bar
    
end

#hmean, hmin and hmax of each element. We define h in such a way that
#h=1 when the area is the same than that of the reference element. That is:
#   h^2 = K/Khat = E/Ehat = 1/sqrt(|M|) = 1/Mbar
function _hElems(mesh::TrMesh)

    return _hElems(mesh.metric)

end     
function _hmin(mesh::TrMesh)

    return _hmin(mesh.metric)
    
end     
function _hmax(mesh::TrMesh)

    return _hmax(mesh.metric)
    
end     
function _hElems(metric::SpectralDecomp{Vector{Float64},Int})

    return @mlv 1.0/sqrt(metric.lambda_bar)

end     
function _hmin(metric::SpectralDecomp{Vector{Float64},Int})

    return @mlv 1.0/sqrt(metric.lambda_bar*metric.sf[1])
    
end     
function _hmax(metric::SpectralDecomp{Vector{Float64},Int})

    return @mlv 1.0/sqrt(metric.lambda_bar*metric.sf[2])
    
end     

#Perform linear reconstruction of an element-wise variable:
function LinearReconstruction(mesh::TrMesh, wv::Vector{Float64}, 
    u::GenVector{Float64})

    v       = zeros(mesh.nVerts)
    sumwv   = zeros(mesh.nVerts)
    @inbounds for iVert=1:3, iElem=1:mesh.nElems
        TargetNode          = mesh.ElemsNodes[iElem,iVert]
        v[TargetNode]       += wv[iElem]*u[iElem]
        sumwv[TargetNode]   += wv[iElem]
    end
    @mlv v  /= sumwv
    
    return v
    
end 

#Adapt mesh with Bamg:
function Bamg_MetricFile!(mesh::TrMesh,
    M::SpectralDecomp{Vector{Float64},Int},
    MetricFile::String)

    nVerts          = mesh.nVerts
    
    #Metric at elements:
    m_elems         = GetMatrix(M)
    
    #Metric at vertices. NOTE: The linear reconstruction at the vertices
    #is a linear combination of the values at the elements with positive coefficients.
    #Hence, if M^1, M^N are the metric tensors at the elements common to a certain vertex,
    #the linear combination (alpha1 M^1 + ... + alphaN M^N) yields another positive definite tensor.
    m_verts         = zeros(mesh.nVerts, 3)
    wv              = M.lambda_bar     #wv=(Khat/K)
    m_verts[:,1]    = LinearReconstruction(mesh, wv, m_elems[1,1])
    m_verts[:,2]    = LinearReconstruction(mesh, wv, m_elems[2,1])
    m_verts[:,3]    = LinearReconstruction(mesh, wv, m_elems[2,2])
    
    #Careful!! m_elems considers that the reference element is of radius 1.
    #m_verts considers that the reference element is of length 1 (area sqrt(3)/4).
    #M1_bar/M2_bar = Khat1/Khat2
    @mlv m_verts    *= sqrt(3)/4/TrElem_Area
    
    #Create metric file:
    file_ptr        = open(MetricFile, "w")
    write(file_ptr, "$(size(m_verts,1)) 3\n")
    for ii=1:nVerts
        write(file_ptr, "$(m_verts[ii,1])\t$(m_verts[ii,2])\t$(m_verts[ii,3])\n")
    end
    close(file_ptr)

    return
    
end

#Search--locate algorithm (Allievi & Bermejo, 1997):
function SearchLocate(mesh::TrMesh, xm::Vector{<:GenVector{Float64}};
    ElemGuess::Vector{Int}=ones(Int, length(xm[1])),
    xiGuess::Matrix{Float64}=zeros(length(xm[1]),2),
    XiTol::Float64=1e-8, XTol::Float64=0.0, MaxIter::Int=2*mesh.nElems)
    
    #Allocate memory for outputs:
    nPoints             = length(xm[1])
    xi0m                = zeros(nPoints,2)
    elem0v              = zeros(Int,nPoints)

    #Check element guess:
    if maximum(ElemGuess)>mesh.nElems
        error("Maximum of ElemGuess must not be greater than the number of elements")
    elseif minimum(ElemGuess)<=0
        error("All element guesses must be positive")
    end

    #Check dimensions:
    if length(xm)!=2 || size(xiGuess,2)!=2
        error("xm and xiGuess must have two columns")
    end
    if length(ElemGuess)!=nPoints || size(xiGuess,1)!=nPoints
        error("xm, ElemGuess and xiGuess must have the same number of rows")
    end

    #Initialize outputs:
    xi0m                .= xiGuess
    elem0v              .= ElemGuess

    #Initialize vector with the points to be searched. The following arrays will
    #have fixed sizes, but we will only use the first nSearched rows:
    nSearched           = nPoints               #Nb of points being searched currently
    PointsSearched      = Array(1:nPoints)
    IsOnEdge            = zeros(Int,nPoints)
    xiv                 = zeros(2)
    fv                  = zeros(2)
    Jinvm               = zeros(2,2)
    pv                  = zeros(2)
    bElems_LastChance   = zeros(Bool,nPoints)

    #Bisectors of each edge:
    nm                  = [ 0           -1
                            sqrt(3)/2   0.5
                            -sqrt(3)/2  0.5 ]
    projv               = zeros(3)
    
#     figure()
    
    #Iterate until we locate all the points:
    nIter               = 0
    flag 				= 1
#     @warn ""
    while true

        #Evaluate shape functions:
        PSRange         = 1:nSearched
        PS_xim          = view(xi0m, PointsSearched[1:nSearched], :)
        Nm              = NCompute(mesh, PS_xim)
        GradNm          = gradNCompute(mesh, PS_xim)
        
        #Loop searched points:
        nSearched_next          = 0
        PointSearched           = 0
        @inbounds for iS=1:nSearched

            #Get current searched point and element:
            PointSearched       = PointsSearched[iS]
            ElemSearched        = elem0v[PointSearched]
            xiv[1]              = PS_xim[iS,1]
            xiv[2]              = PS_xim[iS,2]

            #Compute residual:
            fv[1]               = -xm[1][PointSearched]
            fv[2]               = -xm[2][PointSearched]
            for ii=1:mesh.NodesPerElem
                fv[1]           += Nm[iS,ii]*mesh.ElemsCoords[ElemSearched,2*ii-1]
                fv[2]           += Nm[iS,ii]*mesh.ElemsCoords[ElemSearched,2*ii]
            end
            fnorm           = sqrt(fv[1]^2 + fv[2]^2)
            
#             PyPlot.cla()
#             plot(mesh)
#             plot(mesh, color="r", elems=[ElemSearched])
#             plot(xm[1][PointSearched]+fv[1], xm[2][PointSearched]+fv[2], "*k")
            
            #Compute jacobian and determinant of isoparametric transformation:
            Jinvm           .= 0.0
            for ii=1:mesh.NodesPerElem
                Jinvm[1,1]  += mesh.ElemsCoords[ElemSearched,2*ii]*GradNm[2][iS,ii]
                Jinvm[2,1]  -= mesh.ElemsCoords[ElemSearched,2*ii]*GradNm[1][iS,ii]
                Jinvm[1,2]  -= mesh.ElemsCoords[ElemSearched,2*ii-1]*GradNm[2][iS,ii]
                Jinvm[2,2]  += mesh.ElemsCoords[ElemSearched,2*ii-1]*GradNm[1][iS,ii]
            end
            detJ            = Jinvm[1,1]*Jinvm[2,2]-Jinvm[2,1]*Jinvm[1,2]
            Jinvm[1,1]      /= detJ
            Jinvm[2,1]      /= detJ
            Jinvm[1,2]      /= detJ
            Jinvm[2,2]      /= detJ

            #Compute Newton's update:
            pv              .= 0.0
            for ii=1:2, jj=1:2
                pv[ii]      -= Jinvm[ii,jj]*fv[jj]
            end
            pnorm           = sqrt(pv[1]^2 + pv[2]^2)
            
            #If convergence is achieved:
            if fnorm<=XTol || pnorm<=XiTol
            
#                 display("converged")
                
            else
                
                #Update xi (this modifies xi0m):
                PS_xim[iS,1]    += pv[1]
                PS_xim[iS,2]    += pv[2]
                
                #Project xi onto the bisectors of each edge:
                projv[1]        = nm[1,1]*PS_xim[iS,1]+nm[1,2]*PS_xim[iS,2]
                projv[2]        = nm[2,1]*PS_xim[iS,1]+nm[2,2]*PS_xim[iS,2]
                projv[3]        = nm[3,1]*PS_xim[iS,1]+nm[3,2]*PS_xim[iS,2]
                
                #If all the projections are below 0.5, the point is inside the triangle.
                proj_max, edge  = findmax(projv)
                if proj_max<0.5
                
                    #Keep iterating:
                    nSearched_next                  += 1
                    PointsSearched[nSearched_next]  = PointSearched
                
                else
                
                    #Point leaves the triangle by edge "edge".
                    
                    #The middle point is located at b=0.5*n. Therefore, a=0.5/proj_max*xi
                    #is located at the edge (because it has the same normal projection than b).
                    #Project (a-b) onto tangent vector (-n2, n1):
                    edge_distance   = (0.5*PS_xim[iS,1]/proj_max-0.5*nm[edge,1])*(-nm[edge,2]) + 
                                        (0.5*PS_xim[iS,2]/proj_max-0.5*nm[edge,2])*nm[edge,1]
                    
                    #Change search element:
                    elem_neigh      = mesh.ElemsNeigh[ElemSearched, edge]
                    if elem_neigh==0 && bElems_LastChance[PointSearched]
                        #We are already at the boundary and xi is still outside.
                        #Do not iterate more here.
                        elem0v[PointSearched]   = -ElemSearched
                    elseif elem_neigh==0 #&& !bElems_LastChance[PointSearched]
                    
                        #Give one more chance to iterate from the current edge.
                        bElems_LastChance[PointSearched]    = true
                        PS_xim[iS,1]    = 0.5*nm[edge,1] + edge_distance*(-nm[edge,2])
                        PS_xim[iS,2]    = 0.5*nm[edge,2] + edge_distance*nm[edge,1]
                        
                        #Iterate again:
                        nSearched_next                  += 1
                        PointsSearched[nSearched_next]  = PointSearched
                        
                    else
                        #Set new neighbor:
                        elem0v[PointSearched]   = elem_neigh
                        edge_neigh              = mesh.ElemsNeighFaces[ElemSearched, edge]
                        
                        #Position in new neighbor. We start from the middle point of the
                        #edge and sum edge_distance in the direction opposite to the tangent vector:
                        PS_xim[iS,1]    = 0.5*nm[edge_neigh,1] - edge_distance*(-nm[edge_neigh,2])
                        PS_xim[iS,2]    = 0.5*nm[edge_neigh,2] - edge_distance*nm[edge_neigh,1]
                        
                        #Iterate again:
                        nSearched_next                  += 1
                        PointsSearched[nSearched_next]  = PointSearched
                        
                    end
                
                end
                
            end #if converged
                
        end #loop searched points

        #Update nSearched and nIter:
        nIter       += 1
        nSearched   = nSearched_next

        #Break loop:
        if nSearched==0
            break
        end
        if nIter==MaxIter
            printstyled("Reached maximum of iterations in search--locate. ", 
            "Unable to find ", nSearched, " points \n", color=:light_yellow)
			flag 	= -1
            break
        end

    end

    #DEBUG:
#     if length(xm[1])==100
#         Nx      = 100
#         Nm      = NCompute(mesh, xi0m)
#         xm2     = [ zeros(size(xm[1])), zeros(size(xm[1])) ]
#         for ix=1:Nx, iNode=1:mesh.NodesPerElem
#             xm2[1][ix]   += Nm[ix, iNode]*mesh.ElemsCoords[elem0v[ix], 2*iNode-1]
#             xm2[2][ix]   += Nm[ix, iNode]*mesh.ElemsCoords[elem0v[ix], 2*iNode]
#         end     
#         figure()
#         err     = @mlv sqrt((xm[1]-xm2[1])^2 + (xm[2]-xm2[2])^2) 
#         semilogy(err, ".b")
#         
#     end 
        
    #Return output:
    return flag, elem0v, xi0m, nIter

end

#Find patch of a given element:
function ElemsPatch_Compute(mesh::TrMesh, iElem::Int, order::Int)

    #Get patch of elements:
    ElemsPatch          = [iElem]
    nElemsPatch         = 1
    PatchNodes          = unique(mesh.ElemsNodes[ElemsPatch,1:3][:])
    for ii=1:order
        
        #Loop PatchNodes and add elements:
        for node in PatchNodes
            ElemsPatch      = vcat(ElemsPatch, mesh.ElemsNodesSp.rowval[
                                mesh.ElemsNodesSp.colptr[node]:mesh.ElemsNodesSp.colptr[node+1]-1])
        end
        ElemsPatch          = unique(ElemsPatch)
        nElemsPatch         = length(ElemsPatch)
        PatchNodes          = unique(mesh.ElemsNodes[ElemsPatch,1:3][:])
    end
    
    return nElemsPatch, ElemsPatch, PatchNodes
    
end
