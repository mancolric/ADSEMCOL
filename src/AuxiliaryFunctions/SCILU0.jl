#STATIC CONDENSATION + ILU0 PRECONDITIONER.
#
#The degrees of freedom that are to be eliminated by static condensation (slave dofs) 
#must be placed in first order.
#
#The main idea is that we first pivot only in the slave degrees of freedom, and we arrive at 
#   A = L0 Ared
#where L0 is a lower triangular matrix with all the row operations performed until this point,
#and Ared is a not-yet-upper triangular matrix. At this point, due to the block matrix structure,
#the ILU decomposition is exact, that is, A = L0 Ared is not an approximation.
#
#Hence, if the system is Ax=b, we can solve Ared x = L0\b =:bhat. But the slave-master block 
#of Ared is zero, so this is equivalent to solving Ahat_mm x_m = bhat_m. To solve the latter,
#we keep doing the ILU factorization of Ared (so in the end we need the ILU factorization of A).
#However, when solving Ared_mm x_m = bhat_m, we only use the blocks (m,m) in the matrices L,U,
#as well as we need to 

include("basic.jl")

mutable struct SCILU0Fact{nb<:Number, int<:Integer} <: Factorization{nb}

    #LU matrices:
    LU      ::SparseMatrixCSC{nb,int}
    N       ::int
    
    #pointer to diagonal value of each column:
    diag    ::Vector{int}
    
    #pointer to first master element at each column:
    nSlaves ::Int
    nMasters::Int
    master  ::Vector{int}
    
    #Allocate master*master block of reduced matrix. Upper block is stored in U.
    Ared_mm ::SparseMatrixCSC{nb,int}
    
    #auxiliary vector to perform operations:
    aux     ::Vector{nb}
    
    #CSR matrix to point the upper diagonal terms for each row:
    u_ptr   ::Vector{int}
    u_cols  ::Vector{int}
    u_ss    ::Vector{int}
    
    #
    SCILU0Fact{nb,int}() where {nb<:Number,int<:Integer} = new{nb,int}()
    
end

#---------------------------------------------------------------------------
#STRUCTURE ALLOCATION:

function SCILU0_alloc(A::SparseMatrixCSC{nb,int}, nSlaves::Int) where {nb<:Number,int<:Integer}

    ILU             = SCILU0Fact{nb,int}()
    
    #Size:
    M,N             = size(A)
    if M!=N
        error("Matrix must be square")
    end
    ILU.N           = N
    ILU.nSlaves     = nSlaves
    ILU.nMasters    = N-nSlaves
    
    #Allocate nnz:
    ILU.LU          = deepcopy(A)
    
    #Pointer to diagonal terms and first master term:
    ILU.diag                            = zeros(int, N)
    ILU.master                          = zeros(int, N)
    for jj=1:N
        SearchDiag                      = true
        SearchMaster                    = true
        ILU.master[jj]                  = ILU.LU.colptr[jj+1] #Assume there is no master
        for ss=A.colptr[jj]:A.colptr[jj+1]-1
            row                         = A.rowval[ss]
            if SearchDiag && row==jj
                ILU.diag[jj]            = ss
                SearchDiag              = false
            end
            if SearchMaster && row>nSlaves
                ILU.master[jj]          = ss
                SearchMaster            = false
            end
            if !SearchMaster && !SearchDiag
                break
            end
        end
        
        #Check diagonal element has been found:
        if ILU.diag[jj]         == 0
            error("Diagonal element not found in column $(jj)")
        end
        
    end
    
    #Allocate vector to do operations:
    ILU.aux         = zeros(nb, N)
    
    #Pointers to upper diagonal terms for each row
    #
    #Count number of (upper diagonal) nonzeros per row:
    u_aux           = zeros(int, N)
    for jj=1:N, ss=A.colptr[jj]:ILU.diag[jj]-1
        row         = A.rowval[ss]
        u_aux[row]  += 1
    end
    #Save accumulated number of nonzeros per row:
    u_ptr           = ones(int, N+1)
    for jj=1:N
        u_ptr[jj+1] = u_ptr[jj]+u_aux[jj]
    end
    #Save columns and ss-pointer for each row:
    u_aux           .= 0    #local pointer
    u_cols          = zeros(int, u_ptr[N+1])
    u_ss            = zeros(int, u_ptr[N+1])
    for jj=1:N, ss=A.colptr[jj]:ILU.diag[jj]-1
        row         = A.rowval[ss]
        ptr         = u_ptr[row]+u_aux[row]
        u_cols[ptr] = jj
        u_ss[ptr]   = ss
        u_aux[row]  += 1
    end
    ILU.u_ptr       = u_ptr
    ILU.u_cols      = u_cols
    ILU.u_ss        = u_ss
    
    #Allocate reduced matrix:
    Ared_colptr     = zeros(int, ILU.nMasters+1)
    #Save nb of cum nnz values per column:
    Ared_colptr[1]  = 1
    for jj=1:ILU.nMasters
        #master elements in column nSlaves+jj of matrix A:
        Ared_colptr[jj+1]   = Ared_colptr[jj] + ILU.LU.colptr[nSlaves+jj+1] - 
                                ILU.master[nSlaves+jj]
    end
    #Save rows of each nnz element in Ared_mm:
    Ared_nnz        = Ared_colptr[ILU.nMasters+1]-1
    Ared_rowval     = zeros(int, Ared_nnz)
    ssred           = 0
    for jj=1:ILU.nMasters, ss=ILU.master[nSlaves+jj]:ILU.LU.colptr[nSlaves+jj+1]-1
        row                 = ILU.LU.rowval[ss]
        ssred               += 1
        Ared_rowval[ssred]  = row-nSlaves
    end
    Ared_nzval      = zeros(nb, Ared_nnz)
    ILU.Ared_mm     = SparseMatrixCSC{nb,int}(ILU.nMasters, ILU.nMasters, Ared_colptr, 
                        Ared_rowval, Ared_nzval)
    
    return ILU
    
end

#---------------------------------------------------------------------------
#ILU FACTORIZATION ALGORITHMS (they are the same irrespective of the static condensation):

function SCILU0(A::SparseMatrixCSC{<:Number, <:Integer}, nSlaves::Int)

    ILU           = SCILU0_alloc(A, nSlaves)
    SCILU0!(ILU, A)
    return ILU
    
end

function SCILU0!(ILU::SCILU0Fact{nb,int}, A::SparseMatrixCSC{nb,int}) where {nb<:Number, int<:Integer}

#     return SCILU0_columns!(ILU, A)
    return SCILU0_rows!(ILU, A)
    
end

function SCILU0_rows!(ILU::SCILU0Fact{nb,int}, A::SparseMatrixCSC{nb,int}) where {nb<:Number, int<:Integer}
    
    t_ini           = time()
    
    #Extract variables:
    N               = ILU.N
    LU              = ILU.LU
    colptr          = LU.colptr
    rows            = LU.rowval
    nzval           = LU.nzval
    nzval           .= A.nzval
    diagptr         = ILU.diag
    auxv            = ILU.aux
    @mlv auxv       = 0.0
    u_ptr           = ILU.u_ptr
    u_cols          = ILU.u_cols
    u_ss            = ILU.u_ss
    master          = ILU.master
    Ared_nzval      = ILU.Ared_mm.nzval
    nSlaves         = ILU.nSlaves
    
#     println("Alloc = ", time()-t_ini)
    
    #Loop rows:
    t_pivot         = 0.0
    t_lfacts        = 0.0
    t_columns       = 0.0
    t_erase         = 0.0
    t_reduced       = 0.0
    @inbounds for ii=1:N
        
        #L is currently stored under the diagonal and in columns 1:i-1. The diagonal terms are 1.
        #U is stored in the rest of the matrix.
        #Initially, L=I and U=A.
        
        t_ini2      = time()
        
        #Pivot inverse:
        s_ii        = diagptr[ii]
        a_ii        = nzval[s_ii]
        aii_inv     = 1.0/a_ii
        
        t_pivot     += time()-t_ini2
        t_ini2      = time()
        
        #Save l_ji factors:
        for ss=s_ii+1:colptr[ii+1]-1
            l_ji        = nzval[ss]*aii_inv
            nzval[ss]   = l_ji
            row         = rows[ss]
            auxv[row]   = l_ji
        end
        
        t_lfacts    += time()-t_ini2
        t_ini2      = time()
        
        #Loop columns i+1:N with a_ij!=0:
        for tt=u_ptr[ii]:u_ptr[ii+1]-1 
        
            jj          = u_cols[tt]
            u_ij        = nzval[u_ss[tt]]
            
            #   F'_i    = F_i
            #   F'_j    = F_j - a_ij/a_ii * F_i
            #which means:
            #   F_i     = F'_i
            #   F_j     = F_j' + a_ij/a_ii * F'_i 
            #That is, L^(i+1) A^(i+1) = A^i 
            
            #Perform row operations in U:
            for ss=u_ss[tt]+1:colptr[jj+1]-1
                row         = rows[ss]
                nzval[ss]   -= u_ij*auxv[row]
            end
            
        end
        
        t_columns   += time()-t_ini2
        t_ini2      = time()
        
        #Erase l_ji factors from auxiliary vector:
        for ss=s_ii+1:colptr[ii+1]-1
            row         = rows[ss]
            auxv[row]   = 0.0
        end
        
        t_erase     += time()-t_ini2
        t_ini2      = time()
        
        #Save reduced matrix:
        if ii==nSlaves
            ss_red                  = 0
            for jj=nSlaves+1:N, ss=master[jj]:colptr[jj+1]-1
                ss_red              += 1
                Ared_nzval[ss_red]  = nzval[ss]
            end
            
        end
            
        t_reduced     += time()-t_ini2
        
    end
    
#     println("Loop = ", time()-t_ini)
#     println("t_pivot = ", t_pivot)
#     println("t_lfacts = ", t_lfacts)
#     println("t_columns = ", t_columns)
#     println("t_erase = ", t_erase)
#     println("t_reduced = ", t_reduced)
            
    return
    
end

#---------------------------------------------------------------------------
#GET L,U MATRICES (NOT OPTIMIZED):

#Get sparse matrices for L and U factors:
function LUFactors(F::SCILU0Fact{<:Number,<:Integer})

#     return LUFactors_columns(F)
    return LUFactors_rows(F)
    
end

function LUFactors_rows(F::SCILU0Fact{<:Number,<:Integer}) 

    iv, jv, sv  = findnz(F.LU)
    aux         = findall(iv.<=jv)
    U           = sparse(iv[aux], jv[aux], sv[aux])
    aux         = findall(iv.==jv)
    sv[aux]     .= 1.0
    aux         = findall(iv.>=jv)
    L           = sparse(iv[aux], jv[aux], sv[aux])
    
    return L,U
    
end

#---------------------------------------------------------------------------
#SOLUTION OF LU SYSTEMS:

#Reduce r.h.s., that is, compute L0\b:
function ReduceRhs(F::SCILU0Fact{nb, int}, b::Vector{nb}) where {nb<:Number, int<:Integer}

    colptr  = F.LU.colptr
    rowval  = F.LU.rowval
    nzval   = F.LU.nzval
    diag    = F.diag
    N       = F.N
    auxv    = F.aux
    nSlaves = F.nSlaves
    bred    = copy(b)
    
    @inbounds for jj=1:nSlaves
        bred_jj         = bred[jj]
        for ss=diag[jj]+1:colptr[jj+1]-1
            row         = rowval[ss]
            bred[row]   -= nzval[ss]*bred_jj
        end
    end
    
    return bred
    
end

#Standard ILU solution U\(L\b)):
import LinearAlgebra: ldiv!
function ldiv!(x::AbstractVector{nb}, F::SCILU0Fact{nb, int}, 
    b::AbstractVector{nb}) where {nb<:Number, int<:Integer}
    
#     return ldiv_columns!(x, F, b)
    return ldiv_rows!(x, F, b)
    
end

function ldiv_columns!(x::AbstractVector{nb}, F::SCILU0Fact{nb, int}, 
    b::AbstractVector{nb}) where {nb<:Number, int<:Integer}
    
    colptr  = F.LU.colptr
    rows    = F.LU.rowval
    nzval   = F.LU.nzval
    diag    = F.diag
    N       = F.N
    auxv    = F.aux
    
    #Solve Ux = L\b. Both r.h.s. and Ux are stored in aux:
    @mlv auxv           = b
    @inbounds for jj=1:N
    
        #Solve L_jj*(Ux)_j = aux_j:
        s_jj            = diag[jj]
        Ux_jj           = auxv[jj]/nzval[s_jj]
        auxv[jj]        = Ux_jj     #overwrite aux[j] with (Ux)_j
        
        #Modify r.h.s. and save in aux:
        for ss=s_jj+1:colptr[jj+1]-1
            row         = rows[ss]
            auxv[row]   -= nzval[ss]*Ux_jj
        end
        
    end
    
    #Solve x=U\(Ux). R.h.s. is stored in aux and x is stored in x.
    @inbounds for jj=N:-1:1
    
        #Solve x_j = aux_j:
        x_jj            = auxv[jj]
        x[jj]           = x_jj
        
        #Modify r.h.s.:
        for ss=colptr[jj]:diag[jj]-1
            row         = rows[ss]
            auxv[row]   -= nzval[ss]*x_jj
        end
        
    end
    
    return
    
end

function ldiv_rows!(x::AbstractVector{nb}, F::SCILU0Fact{nb, int}, 
    b::AbstractVector{nb}) where {nb<:Number, int<:Integer}
    
    colptr  = F.LU.colptr
    rows    = F.LU.rowval
    nzval   = F.LU.nzval
    diag    = F.diag
    N       = F.N
    auxv    = F.aux
    
    #Solve Ux = L\b. Both r.h.s. and Ux are stored in aux:
    @mlv auxv           = b
    @inbounds for jj=1:N
    
        #Solve L_jj*(Ux)_j = aux_j. Since L_jj=1. (Ux)_j=aux_j
        Ux_jj           = auxv[jj]
        
        #Modify r.h.s. and save in aux:
        for ss=diag[jj]+1:colptr[jj+1]-1
            row         = rows[ss]
            auxv[row]   -= nzval[ss]*Ux_jj
        end
        
    end
    
    #Solve x=U\(Ux). R.h.s. is stored in aux and x is stored in x.
    @inbounds for jj=N:-1:1
    
        #Solve x_j = aux_j:
        s_jj            = diag[jj]
        x_jj            = auxv[jj]/nzval[s_jj]
        x[jj]           = x_jj
        
        #Modify r.h.s.:
        for ss=colptr[jj]:s_jj-1
            row         = rows[ss]
            auxv[row]   -= nzval[ss]*x_jj
        end
        
    end
    
    return
    
end

#Solve only for master dof, i.e.,
#   A_red x_red         = b_red
#is approximated by
#   L_red U_red x_red   = b_red
function ldiv_master!(x::AbstractVector{nb}, F::SCILU0Fact{nb, int}, 
    b::AbstractVector{nb}) where {nb<:Number, int<:Integer}
    
    colptr  = F.LU.colptr
    rows    = F.LU.rowval
    nzval   = F.LU.nzval
    diag    = F.diag
    N       = F.N
    auxv    = F.aux
    nSlaves = F.nSlaves
    nMasters= F.nMasters
    master  = F.master
    
    #Solve Ux = L\b. Both r.h.s. and Ux are stored in aux:
#     @mlv @views auxv[1:nMasters]     = b
    @inbounds for jj=1:nMasters
        auxv[jj]        = b[jj]
    end
    @inbounds for jj=1:nMasters
    
        #Solve L_jj*(Ux)_j = aux_j. Since L_jj=1. (Ux)_j=aux_j
        Ux_jj           = auxv[jj]
        
        #Modify r.h.s. and save in aux:
        for ss=diag[nSlaves+jj]+1:colptr[nSlaves+jj+1]-1
            row         = rows[ss]-nSlaves
            auxv[row]   -= nzval[ss]*Ux_jj
        end
        
    end
    
    #Solve x=U\(Ux). R.h.s. is stored in aux and x is stored in x.
    @inbounds for jj=nMasters:-1:1
    
        #Solve x_j = aux_j:
        s_jj            = diag[nSlaves+jj]
        x_jj            = auxv[jj]/nzval[s_jj]
        x[jj]           = x_jj
        
        #Modify r.h.s.:
        for ss=master[nSlaves+jj]:s_jj-1
            row         = rows[ss]-nSlaves
            auxv[row]   -= nzval[ss]*x_jj
        end
        
    end
    
    return

end

#Once master dof have been solved, compute slave dof. The values of x_master are stored
#in the slice nSlaves+1:N of x:
function ldiv_slave!(x::AbstractVector{nb}, F::SCILU0Fact{nb, int}, 
    bred::AbstractVector{nb}) where {nb<:Number, int<:Integer}
    
    colptr  = F.LU.colptr
    rows    = F.LU.rowval
    nzval   = F.LU.nzval
    diag    = F.diag
    N       = F.N
    auxv    = F.aux
    nSlaves = F.nSlaves
    nMasters= F.nMasters
    master  = F.master
    
    #Compute r.h.s. of slave dof:
    @mlv auxv       = bred
    @inbounds for jj=N:-1:nSlaves+1, ss=colptr[jj]:master[jj]-1
        row         = rows[ss]
        auxv[row]   -= nzval[ss]*x[jj]
    end
    
    #Solve Ux = aux:
    @inbounds for jj=nSlaves:-1:1
    
        #Solve x_j = aux_j:
        s_jj            = diag[jj]
        x_jj            = auxv[jj]/nzval[s_jj]
        x[jj]           = x_jj
        
        #Modify r.h.s.:
        for ss=colptr[jj]:s_jj-1
            row         = rows[ss]
            auxv[row]   -= nzval[ss]*x_jj
        end
        
    end
    
    return

end

#If the error "x" is measured by a weighted L1 norm of the form
#   ||x|| = sum_s a_s |x_s| + sum_m a_m |x_m|, with a_i>0,
#after reducing the number of d.o.f. we have to express this norm in terms of |x_m| only.
#For that purpose, we assume that the slaves dof are solved exactly, that is, given the master
#dofs u_m,
#   A_ss u_s + A_sm u_m = b_s
#is satisfied exactly. Since
#   A_ss u^ex_s + A_sm u^ex_m = b_s
#we have
#   A_ss x_s + A_sm x_m = 0_s
#The latter is approximated by
#   x_s     = - D_ss^(-1) A_sm x_m,
#with D_ss = diag(A_ss)
#If M_ss and M_mm are diagonal matrices with the weights a_s and a_m, then
#   M_ss x_s    = - M_ss D_ss^(-1) A_sm M_mm^(-1) M_mm x_m,
#Hence,
#   |M_ss x_s|      <= |M_ss D_ss^(-1) A_sm M_mm^(-1)| |M_mm x_m|
#   sum_s a_s |x_s| <= |M_ss D_ss^(-1) A_sm M_mm^(-1)| sum_m a_m |x_m|
#where |·| is the L1 norm. So, we have to compute the norm of Ahat_sm = M_ss D_ss^(-1) A_sm M_mm^(-1). 
#Note that 
#   (Ahat_sm)_ij    = (a_s)_i / (A_ss)_ii * (A_sm)_ij / (a_m)_j
#and
#   |Ahat_sm|       = max_j sum_i |(Ahat_sm)_ij|
#Then
#   ||x||           = sum_s a_s |x_s| + sum_m a_m |x_m| <= (|Ahat_sm|+1) sum_m a_m |x_m|
#Here we compute the factor |Ahat_sm|+1:
function L1Norm_CorrFactor(F::SCILU0Fact{nb, int}, normv::Vector{nb}) where {nb<:Number, int<:Integer}

    error("Unfinished")
    
    a_s                 = normv[1:F.nSlaves]
    a_m                 = normv[F.nSlaves+1:F.N]

    #Scaling vectors:
#     row_scal            = 

end
