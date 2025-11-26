include("../AuxiliaryFunctions/basic.jl")
import FastGaussQuadrature.gausslegendre

#NOTE: Julia does not allow circular references inside structures.
#Hence, we cannot save the corresponding "mesh" in TrBmesh.
#
#Likewise, we cannot associate a finite element space with a given mesh and vice versa.

#-------------------------------------------------------------
#Meshes:

#Boundary triangular mesh:
mutable struct TrBmesh

    #Mesh data:
    boundary_id ::Int
    dim         ::Int
    nElems      ::Int
    nNodes      ::Int
    nVerts      ::Int
    nEdges      ::Int
    NodesPerElem::Int
    VertsPerElem::Int
    order       ::Int

    #Mesh matrices:
    ElemsNodes  ::Array{Int,2}
    NodesCoords ::Array{Float64,2}
    ElemsCoords ::Array{Float64,2}
#     ElemsNodesSp::SparseMatrixCSC{Int,Int}
#     ElemsNeigh  ::Array{Int,2}

    #Parent mesh:
#     ChildNodes  ::SparseVector{Int,Int}
    ParentEdges ::Vector{Int}
    ParentNodes ::Vector{Int}
    ParentElems ::Vector{Int}
#     ParentFaces ::Vector{Int}
    ParentFace  ::Int
    
    #Matrix to compute Lagrangian shape functions associated with isoparametric
    #transformations:
    Lag_Leg     ::Matrix{Float64}
    
    #Natural coordinates of the nodes:
    xi          ::Vector{Float64}
    
    #Create mesh:
    TrBmesh()   = new()

end

#Triangular mesh:
mutable struct TrMesh

    #Mesh data:
    dim             ::Int
    nElems          ::Int
    nNodes          ::Int
    nVerts          ::Int
    nEdges          ::Int
    NodesPerElem    ::Int
    VertsPerElem    ::Int
    order           ::Int

    #Mesh matrices:
    ElemsNodes      ::Matrix{Int}
    NodesCoords     ::Array{Float64,2}
    ElemsCoords     ::Array{Float64,2}
    ElemsNodesSp    ::SparseMatrixCSC{Int,Int}
    #
    ElemsNeigh      ::Matrix{Int}   #<nElems,3>
    ElemsNeighFaces ::Matrix{Int}   #<nElems,3>
    #
    ElemsEdges      ::Matrix{Int}   #<nElems,3>
    EdgesNodes      ::Matrix{Int}   #<nElems,order+1>
    EdgesElems      ::Matrix{Int}   #<nElems,2>
    EdgesFaces      ::Matrix{Int}   #<nElems,2>
    
    #Boundaries:
    nBounds         ::Int
    bmesh           ::Vector{TrBmesh}
    
    #Matrix to compute Lagrangian shape functions associated with isoparametric
    #transformations:
    Lag_PKD         ::Matrix{Float64}
    
    #Natural coordinates of the nodes:
    xi              ::Matrix{Float64}

    #Metric tensor at the barycenter, assuming straight mesh:
    metric          ::SpectralDecomp{Vector{Float64}, Int}
#     Av              ::Vector{Float64}
#     ATotal          ::Float
    
    #Create mesh:
    TrMesh()        = new()

end

TrUMesh     = Union{TrMesh, TrBmesh}    #Union of triangular meshes

#-------------------------------------------------------------
#Finite element spaces:

#Lagrangian continuous finite element space, with constant polynomial degree at each element:
mutable struct TrPSpace
    mesh            ::TrMesh
    order           ::Int
    DofPerElem      ::Int
    ElemsDof        ::Matrix{Int}
    nDof            ::Int
    nDofU           ::Int       #Dof nDofU+1:nDof can be eliminated by static condensation
    xi              ::Matrix{Float64}
    Lag_PKD         ::Matrix{Float64}
    NodesCoords     ::Matrix{Float64}
    TrPSpace()      = new()
end

#PKD discontinuous finite element space, with constant polynomial degree at each element:
mutable struct TrQSpace
    mesh            ::TrMesh
    order           ::Int
    DofPerElem      ::Int
    ElemsDof        ::Matrix{Int}
    nDof            ::Int
    nDofU           ::Int
    TrQSpace()      = new()
end

#Bubble space. Note that B_r = b*P_(r-1)
mutable struct TrBSpace
    mesh            ::TrMesh
    order           ::Int
    DofPerElem      ::Int
    ElemsDof        ::Matrix{Int}
    nDof            ::Int
    xi              ::Matrix{Float64}
    Lag_PKD         ::Matrix{Float64}
    Bub_Lag         ::Matrix{Float64}
    TrBSpace()      = new()
end

#Polynomial + bubble space. Note that B_r = b*P_(r-1)
mutable struct TrPBSpace
    mesh            ::TrMesh
    order           ::Int
    DofPerElem      ::Int
    ElemsDof        ::Matrix{Int}
    nDof            ::Int
    nDofU           ::Int
    PSpace          ::TrPSpace
    BSpace          ::TrBSpace
    dofP            ::Vector{Int}   #indicates which dof's from 1:DofPerElem correspond to the PSpace
    dofB            ::Vector{Int}
    TrPBSpace()     = new()
end

#Create unions of FES:
TrFES       = Union{TrPSpace, TrQSpace, TrBSpace, TrPBSpace}
FES_Union   = Union{TrFES}

#=
#In order to create vectors of FES, all these must be packed into the same 
#type FES. FES is a parametric type, depending on the specific FES_Type:
#DEPRECATED
mutable struct FES{FES_Type<:FES_Union}    
    FES                                             ::FES_Type
    FES(fes::FES_Type) where {FES_Type<:FES_Union}  = new{FES_Type}(fes)
end
=#

#-------------------------------------------------------------
#Integration:

#Quadrature rule in the triangle:
mutable struct TrQRule 
    order           ::Int
    nqp             ::Int
    w               ::Vector{Float64}
    xi              ::Matrix{Float64}
    TrQRule()       = new()
end

#Quadrature rule in a straight 1D line:
mutable struct LineQRule 
    order           ::Int
    nqp             ::Int
    w               ::Vector{Float64}
    xi              ::Vector{Float64}
    LineQRule()     = new()
end

#Integral over a triangular mesh:
mutable struct TrInt
    mesh            ::TrMesh
    order           ::Int           #tentative order
    QRule           ::TrQRule        #QRule.order is the real order  
    x               ::Vector{Matrix{Float64}}
    wdetJ           ::Matrix{Float64}
    Jinv11          ::Matrix{Float64}
    Jinv21          ::Matrix{Float64}
    Jinv12          ::Matrix{Float64}
    Jinv22          ::Matrix{Float64}
    Jinv            ::Matrix{Matrix{Float64}}   #Pack pointers
    TrInt()         = new()
end

#Integral over a boundary of a triangular mesh:
mutable struct TrBint
    mesh            ::TrMesh        #Parent mesh
    bmesh           ::TrBmesh       #Boundary mesh
    order           ::Int           #tentative order
    QRule           ::LineQRule     #QRule.order is the real order  
    x               ::Vector{Matrix{Float64}}
    xi2D            ::Vector{Matrix{Float64}}   #Trace of natural coordinates at each face
    wdetJ           ::Matrix{Float64}
    Jinv            ::Matrix{Matrix{Float64}}
    t_qp            ::Vector{Matrix{Float64}}
    n_qp            ::Vector{Matrix{Float64}}
    TrBint()        = new()
end

#-------------------------------------------------------------
#Load files:

include("FES.jl")
include("TrElem.jl")
include("TrMesh.jl")
include("TrInt.jl")
include("AMA.jl")
