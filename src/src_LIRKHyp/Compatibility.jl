#Boundary triangular mesh:
mutable struct TrBmesh_v1

    #Mesh data:
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

const JLD2_Bmesh_v1 = JLD2.ReconstructedMutable{:TrBmesh, (:dim, :nElems, :nNodes, :nVerts, :nEdges, :NodesPerElem, :VertsPerElem, :order, :ElemsNodes, :NodesCoords, :ElemsCoords, :ParentEdges, :ParentNodes, :ParentElems, :ParentFaces, :Lag_Leg, :xi), Tuple{Int64, Int64, Int64, Int64, Int64, Int64, Int64, Int64, Vararg{Any, 9}}}

const TrBmesh_v2    = TrBmesh

import Base.convert
function convert(TrBmesh_v2, bmesh0::TrBmesh_v1)

    bmesh                   = TrBmesh()
    
    #Add new field:
    bmesh.boundary_id       = 0
    
    #All the fields in v1 are defined also in v2, so just copy them all:
    for ii=1:nfields(solver_red)
        field_name          = fieldname(TrBmesh_v1, ii)
        setfield!(bmesh, field_name, getfield(bmesh0, field_name))
    end
    
    return bmesh
    
end
function convert(TrBmesh_v2, bmesh::JLD2_Bmesh_v1)

    @warn "Unfinished"
    
    bmesh                   = TrBmesh()
    
    return bmesh
    
end


