#-------------------------------------------------------------------
#FEM library contains the basic finite element routines.

#LIRKHyp library solves a problem of the form
#   du_I/dt + df_Ij(u,x,t)/dx_j = Q_I(u,x,t)
#
#In this file, we declare all the necessary structures.
#
#In "LIRKHypSolver", we have the functions to start the LIRKHyp solver,
#set default options, project initial condition and march in time.
#
#In "Matrices", we create some basic vectors / matrices, e.g., for the 
#L2 projection of the initial condition.
#
#In PlotUtils, we create some functions to easily plot the solution.
#
#The core of the method are the functions to compute f_Ij, Q_I. These 
#are the constitutive laws. 

#For structures with default parameters, we use kwdef. 
#For structures with undefined parameters, we define the constructor new().

include("../FEMFunctions/FEMFunctions.jl") #this includes basic
include("../AuxiliaryFunctions/ODESolvers.jl")
include("../AuxiliaryFunctions/NonlinearSolvers.jl")
include("../AuxiliaryFunctions/SCILU0.jl")
using Format
using CuthillMcKee

#FGMRES:
#https://copyprogramming.com/howto/how-to-implement-flexible-gmres-in-matlab#high-quality-flexible-gmres-fgmres-implementation

abstract type ConstModels end
abstract type BoundConds end

#Wrapper for boundary conditions:
mutable struct BCW{BCType<:BoundConds}
    BC                  ::BCType
end
BCW(bc::BCType) where BCType<:BoundConds = BCW{BCType}(bc)

#Wrapper for linear systems/preconditioners:
abstract type LinearSystems end

mutable struct LSW{LSType<:LinearSystems}
    LS                  ::LSType
end
LSW(ls::LSType) where LSType<:LinearSystems = LSW{LSType}(ls)

#Static condensation + ILU0 + GMRES. The solution of a linear system involves 
#1) static condensation, 2) ILU0+GMRES for master dofs, 3) computation of slave dofs
#from static condensation. It is not possible to separate preconditioner from linear 
#solver; everything is coupled.
mutable struct LS_SCILU0_GMRES <: LinearSystems

    A                   ::SparseMatrixCSC{Float64,Int}  #Original matrix
    APP                 ::SparseMatrixCSC{Float64,Int}  #Permuted matrix
    p                   ::Vector{Int}
    pinv                ::Vector{Int}
    ssPP_ss             ::Vector{Int}                   #Vector to compute APP quickly
    Pl                  ::SCILU0Fact{Float64, Int}
    scaleP              ::Vector{Float64}               #Vector to scale solution with NFactors
    scaleP_m            ::Vector{Float64}               #Vector to scale master dof solution with NFactors
    
    #nMaster, nSlaves, etc is in Pl
        
    NormFun             ::FW_NLS_norm                   #Function to compute the norm
    
    LS_SCILU0_GMRES()   = new()
    
end

#Block Jacobi preconditioner:
mutable struct LS_BlockJacobi <: LinearSystems

    A                   ::SparseMatrixCSC{Float64,Int}  #Original matrix
    APP                 ::SparseMatrixCSC{Float64,Int}  #Permuted matrix
    p                   ::Vector{Int}
    pinv                ::Vector{Int}
    ssPP_ss             ::Vector{Int}                   #Vector to compute APP quickly
    Pl                  ::SCILU0Fact{Float64, Int}      #ILU0 (exact) factorization of block Jacobi
    scaleP              ::Vector{Float64}               #Vector to scale solution with NFactors
    
    #nMaster, nSlaves, etc is in Pl
    
    LS_BlockJacobi()    = new()
    
end

#ILU0 factorization of row-dominant terms:
mutable struct LS_DominantILU0 <: LinearSystems

    A                   ::SparseMatrixCSC{Float64,Int}  #Original matrix
    APP                 ::SparseMatrixCSC{Float64,Int}  #Permuted matrix
    p                   ::Vector{Int}
    pinv                ::Vector{Int}
    ssPP_ss             ::Vector{Int}                   #Vector to compute APP quickly
    Pl                  ::SCILU0Fact{Float64, Int}      #ILU0 factorization of approximate Jacobian
    scaleP              ::Vector{Float64}               #Vector to scale solution with NFactors
    
    #nMaster, nSlaves, etc is in Pl
    
    LS_DominantILU0()   = new()
    
end

#Row-dominant terms:
mutable struct LS_Dominant <: LinearSystems

    A                   ::SparseMatrixCSC{Float64,Int}  #Original matrix
    APP                 ::SparseMatrixCSC{Float64,Int}  #Permuted matrix
    p                   ::Vector{Int}
    pinv                ::Vector{Int}
    ssPP_ss             ::Vector{Int}                   #Vector to compute APP quickly
    Pl                  ::SCILU0Fact{Float64, Int}      #ILU0 factorization of approximate Jacobian
    scaleP              ::Vector{Float64}               #Vector to scale solution with NFactors
    
    #nMaster, nSlaves, etc is in Pl
    
    LS_Dominant()       = new()
    
end

#Variables at quadrature nodes:
mutable struct TrIntVars

    Integ2D         ::TrInt
    nElems          ::Int
    nqp             ::Int
    FesOrder        ::Int
    nFacts          ::Vector{Float64}
    
    t               ::Float64
    x               ::Vector{Matrix{Float64}}
    u               ::Vector{Matrix{Float64}}
    uB              ::Vector{Matrix{Float64}}
    gradu           ::Matrix{Matrix{Float64}}
    graduB          ::Matrix{Matrix{Float64}}
    
#     fB              ::Matrix{Matrix{Float64}}
#     dfB_du          ::Array{Matrix{Float64},3}
#     dfB_duB         ::Array{Matrix{Float64},3}
#     dfB_dgradu      ::Array{Matrix{Float64},4}
#     dfB_dgraduB     ::Array{Matrix{Float64},4}
    
    f               ::Matrix{Matrix{Float64}}
    df_du           ::Array{Matrix{Float64},3}
    df_dgradu       ::Array{Matrix{Float64},4}
    
    fB              ::Matrix{Matrix{Float64}}
    dfB_du          ::Array{Matrix{Float64},3}
    dfB_dgraduB     ::Array{Matrix{Float64},4}
    
    Q               ::Vector{Matrix{Float64}}
    dQ_du           ::Matrix{Matrix{Float64}}
    dQ_dgradu       ::Array{Matrix{Float64},3}
    
    Deltat_CFL      ::Float64
    
    TrIntVars()     = new{}()
    
end
#
mutable struct TrBintVars

    Binteg2D        ::TrBint
    nElems          ::Int
    nqp             ::Int
    FesOrder        ::Int
    nFacts          ::Vector{Float64}
    
    tb              ::Vector{Matrix{Float64}}
    nb              ::Vector{Matrix{Float64}}
    
    t               ::Float64
    x               ::Vector{Matrix{Float64}}
    u               ::Vector{Matrix{Float64}}
#     uB              ::Vector{Matrix{Float64}}
    gradu           ::Matrix{Matrix{Float64}}
#     graduB          ::Matrix{Matrix{Float64}}
#     epsilon         ::Matrix{Float64}
    
    f               ::Vector{Matrix{Float64}}
    df_du           ::Array{Matrix{Float64},2}
#     df_duB          ::Array{Matrix{Float64},2}
    df_dgradu       ::Array{Matrix{Float64},3}
#     df_dgraduB      ::Array{Matrix{Float64},3}
    
    TrBintVars()    = new{}()
    
end

#Solver data:
mutable struct SolverData{ConstModel<:ConstModels}
    
    #Constitutive model:
    model           :: ConstModel
    
    #Problem name and ID:
    ProblemName     :: String
    SC              :: Int
    
    #Mesh:
    MeshFile        :: String       #Mesh file provided by user
    nBounds         :: Int
    mesh            :: TrMesh                   
    Omega           :: Float64
    Av              :: Vector{Float64}
    
    #Initial condition and boundary conditions:
    u0fun           :: FW11
    BC              :: Vector{BCW}
    
    #Solution:
    nVars           :: Int                      
    FesOrder        :: Int                      
    fes             :: TrPBSpace       
    uv              :: Vector{Float64}
    u               :: Vector{VectorView{Float64}} 
    SpaceAdapt      :: Bool
    nFacts          :: Vector{Float64}
    SpaceNorm       :: Float64
    
    #Reconstruction:
    rec_fes         :: TrQSpace                         #fes for the reconstruction
    urec            :: Vector{VectorView{Float64}}      #reconstruction
#     MQQ             :: SparseMatrixCSC{Float64,Int}     #mass matrix for reconstruction space
#     MQP             :: SparseMatrixCSC{Float64,Int}     #mass matrix from PSpace to reconstruction space
    
    #Discontinuity-capturing variables:
#     epsilon_fes     :: TrPSpace
#     epsilonv        :: Vector{Float64} 
#     MRR             :: SparseMatrixCSC{Float64,Int}     #mass matrix to compute residual
#     MRP             :: SparseMatrixCSC{Float64,Int}     #idem
#     res_fes         :: TrPBSpace
    
    #Time integration method:
    t               :: Float64
    Nt              :: Int
    RKMethod        :: String 
    Deltat0         :: Float64
    Deltat          :: Float64  
    CFL             :: Float64  
    tf              :: Float64  
    TimeAdapt       :: Bool   
    RK              :: IMEX_RK      #Structure with RK coefficients
    f_RK            :: Matrix{Float64}
    Ju_RK           :: Matrix{Float64}
    fNStab_RK       :: Matrix{Float64}
    
    #Linear solver type. Available options:
    #SCILU0:            All Jacobian terms are computed, system is solved via static condensation + ILU0 + GMRES
    #BlockJacobi:       Jacobian approximated by BlockJacobi preconditioner
    #BlockJacobiKrylov: BlockJacobi preconditioner, Krylov approximation for Jacobian 
    LSType          :: String       
    #Separate, if possible, preconditioner and Krylov approximation:
    JPatt           :: String
    KrylovApprox    :: Bool
    #(i,j) pairs of Jacobian terms to be computed at each element
    Jm_ilocal       :: Vector{Int}  
    Jm_jlocal       :: Vector{Int}
    Jm_localterms   :: Int
    
    
    #Matrices/vectors:
    Integ2D         :: TrInt
    Binteg2D        :: Vector{TrBint}
#     elems_RCM       :: Vector{Int}                      #Optimal permutations for the edge dofs
#     dof_RCM         :: Vector{Int}                      #Optimal permutations for degrees of freedom
    MII             :: SparseMatrixCSC{Float64,Int}     #Mass matrix for each variable 
    MII_LS          :: LS_SCILU0_GMRES
    bv              :: Vector{Float64}
    b               :: Vector{VectorView{Float64}}
    Mm              :: SparseMatrixCSC{Float64,Int}     #JPatt part (complete or block Jacobi) of the mass matrix
    Jm              :: SparseMatrixCSC{Float64,Int}
    Jm_pinv         :: Matrix{Vector{Int}}              #Matrix nVars*nVars with vectors for fast assembly
    Am              :: SparseMatrixCSC{Float64,Int}
    Am_LS           :: LSW                              #Wrapper to linear system
    
    #Monitor variables:
    monitor         :: Vector{Matrix{Float64}}
    
    #Tolerances and safety factors:
    TolS_min        :: Float64
    TolS_max        :: Float64
    TolT            :: Float64
    TolA_min        :: Float64
    CA              :: Float64
    CA_max          :: Float64
    AMA_MaxIter     :: Int
    AMA_SizeOrder   :: Int          #Assume err~h^(SizeOrder+1) for size adaptation
    AMA_AnisoOrder  :: Int          #Use derivatives of order (AnisoOrder+1) for anisotropy
    AMA_ProjN       :: Int          #Number of subtriangles for remeshing
    AMA_ProjOrder   :: Int          #Order of quadrature rule for remeshing
#     AMA_RefineFact  :: Float64      
#     AMA_CoarseFact  :: Float64 
    SfS             :: Float64      #Safety factor for space adaptation
    SfT             :: Float64      #Safety factor for time adaptation
    AMA_rDeltat     :: Float64      #Change in Deltat after AMA
    LS_iters_target :: Int
    LS_iters_max    :: Int
    
    #Errors:
    etaS            :: Float64      #Error at end of time step
    etaS_elems      :: Vector{Float64}
    etaT            :: Float64  
    etaA            :: Float64  
    dlogetaS_dt     :: Float64  
    etaST           :: Float64
    LS_iters        :: Float64      #Mean number of LS iters per stage
    AMA_TolS_etaS   :: Float64      #Ratio between input TolS and produced etaS
        
    #Vector with time levels, errors, number of elements, etc.:
    tv              :: Vector{Float64}
    etaSv           :: Vector{Float64}
    etaTv           :: Vector{Float64}
    etaAv           :: Vector{Float64}
    nElemsv         :: Vector{Int}
    nDofv           :: Vector{Int}
#     tCPUv           :: Vector{Float64}
    CFLv            :: Vector{Float64}
    validv          :: Vector{Int}          #0: unaccepted time step, 1: accepted, -1: not decided
    
    #Profiling vars:
    tCPU            :: Float64  #Total time
    tAMA            :: Float64  #Mesh adaptation
    tJm             :: Float64  #Time spent on computing the Jacobian
    tSCILU          :: Float64  #ILU factorization + static condensation
    tLS             :: Float64  #Solution of linear systems at RK stages
    tRhs            :: Float64  #Time spent on computing r.h.s. vector (excluding first stage, where Jacobian is computed)
    tAlloc          :: Float64  #Memory allocation
    tb              :: Float64  #Compute b 
    tetaS           :: Float64  #Computation of etaS
    tetaT           :: Float64  #Computation of etaT
    LS_total        :: Int      #Total number of LS iterations at RK stages
    
    #Function to create empty structure:
    SolverData{ConstModel}() where {ConstModel<:ConstModels} = new{ConstModel}() 
    
end

#Reduced solver data structure, useful to save results:
mutable struct SolverDataSave{ConstModel<:ConstModels}

    #Problem name and ID:
    ProblemName     :: String
    SC              :: Int
    
    #Constitutive model:
#     model           :: ConstModel     #do not save. It yields errors if ConstModel is changed
    
    #Mesh:
    MeshFile        :: String       #Mesh file provided by user
    nBounds         :: Int
    mesh            :: TrMesh                   
    Omega           :: Float64
    Av              :: Vector{Float64}
    
    #Solution:
    nVars           :: Int                      
    FesOrder        :: Int                      
    fes             :: TrPBSpace       
    uv              :: Vector{Float64}
    u               :: Vector{VectorView{Float64}} 
    SpaceAdapt      :: Bool
    nFacts          :: Vector{Float64}
    SpaceNorm       :: Float64
    
    #Reconstruction:
    rec_fes         :: TrQSpace                         #fes for the reconstruction
    urec            :: Vector{VectorView{Float64}}      #reconstruction
    
    #Time integration method:
    t               :: Float64
    Nt              :: Int
    RKMethod        :: String 
    Deltat0         :: Float64
    Deltat          :: Float64  
    CFL             :: Float64  
    tf              :: Float64  
    TimeAdapt       :: Bool   
    
    #Linear solver:
#     LSType          :: String
    
    #Monitor variables at quadrature nodes:
    monitor         :: Vector{Matrix{Float64}}
    
    #Tolerances and safety factors:
    TolS_min        :: Float64
    TolS_max        :: Float64
    TolT            :: Float64
    TolA_min        :: Float64
    CA              :: Float64
    CA_max          :: Float64
    AMA_MaxIter     :: Int
    AMA_SizeOrder   :: Int
    AMA_AnisoOrder  :: Int          #Use derivatives of order (AnisoOrder+1) for anisotropy
    AMA_ProjN       :: Int          #Number of subtriangles for remeshing
    AMA_ProjOrder   :: Int          #Order of quadrature rule for remeshing
#     AMA_RefineFact  :: Float64      
#     AMA_CoarseFact  :: Float64      
    SfS             :: Float64      #Safety factor for space adaptation
    SfT             :: Float64      #Safety factor for time adaptation
    AMA_rDeltat     :: Float64      #Change in Deltat after AMA
    LS_iters_target :: Int
    LS_iters_max    :: Int
    
    #Errors:
    etaS            :: Float64      #Error at end of time step
    etaS_elems      :: Vector{Float64}
    etaT            :: Float64  
    etaA            :: Float64  
    dlogetaS_dt     :: Float64  
    etaST           :: Float64
    LS_iters        :: Float64      #Mean number of LS iters per stage
    
    #Vector with time levels, errors, number of elements, etc.:
    tv              :: Vector{Float64}
    etaSv           :: Vector{Float64}
    etaTv           :: Vector{Float64}
    etaAv           :: Vector{Float64}
    nElemsv         :: Vector{Int}
    nDofv           :: Vector{Int}
#     tCPUv           :: Vector{Float64}
    CFLv            :: Vector{Float64}
    validv          :: Vector{Int}          #0: unaccepted time step, 1: accepted, -1: not decided
    
    #Profiling vars:
    tCPU            :: Float64  #Total time
    tAMA            :: Float64  #Mesh adaptation
    tJm             :: Float64  #Time spent on computing the Jacobian
    tSCILU          :: Float64  #ILU factorization + static condensation
    tLS             :: Float64  #Solution of linear systems at RK stages
    tRhs            :: Float64  #Time spent on computing r.h.s. vector (excluding first stage, where Jacobian is computed)
    tAlloc          :: Float64  #Memory allocation
    tb              :: Float64  #Compute b 
    tetaS           :: Float64  #Computation of etaS
    tetaT           :: Float64  #Computation of etaT
    LS_total        :: Int      #Total number of LS iterations at RK stages
    
    #Function to create empty structure:
    SolverDataSave{ConstModel}() where {ConstModel<:ConstModels} = new{ConstModel}() 
    
end
import JLD2:save
function save(solver::SolverData{ConstModel}) where ConstModel<:ConstModels
    
    #Create empty structure:
    solver_red      = SolverDataSave{ConstModel}()
    
    #Copy pointers:
    for ii=1:nfields(solver_red)
        field_name          = fieldname(SolverDataSave{ConstModel}, ii)
#         println(field_name)
        setfield!(solver_red, field_name, getfield(solver, field_name))
    end
    
    return solver_red
    
end

function DisplayTimes(solver::SD) where SD<:Union{SolverData,SolverDataSave}

    t_total     = 0.0
    for ii=1:nfields(solver)
        field_name      = String(fieldname(SD, ii))
        if field_name[1]=='t' && all(field_name.!=["t", "tv", "tf"])
            field_value     = getfield(solver, Symbol(field_name))
            println(field_name, "=", field_value, 
                    ", ", sprintf1("%.2f", 100*field_value/solver.tCPU), "%")
            t_total     += field_value
        end
    end
    t_total             -= solver.tCPU
    t_others            = solver.tCPU - t_total
    println("tNotCounted", "=", t_others, 
                    ", ", sprintf1("%.2f", 100*t_others/solver.tCPU), "%")
    return
    
end
    
include("LIRKHypSolver.jl")
include("Matrices.jl")
include("PlotUtils.jl")
include("Compatibility.jl")
