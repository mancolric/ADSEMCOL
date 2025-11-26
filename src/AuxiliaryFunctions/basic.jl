using LinearAlgebra
using SparseArrays
using PyPlot;   pygui(true);
using LoopVectorization
import LoopVectorization:@tturbo
using DelimitedFiles
using FunctionWrappers
import FunctionWrappers:FunctionWrapper
using SparseDiffTools
using JLD2
# using Formatting
using Format

#Set number of BLAS threads according to input from command line:
BLAS.set_num_threads(Threads.nthreads())

#Multithreaded vectorized loop:
macro mlv(ex)
    
    return :(@__dot__($(esc(ex))))

#     ex1     = :(@__dot__($(esc(ex))))
#     return :(@tturbo($(esc(ex1))))
    
#     :(@tturbo @__dot__ $ex)
#     esc(:(LoopVectorization.@tturbo Base.@. $ex))
    
end

#=
To profile, use
    using Profile
    using PProf
    Profile.clear()
    @profile myfun()
    pprof()
See also 
    https://github.com/JuliaPerf/PProf.jl
for memory allocation.
=#

#Automatically convert vectors to (m,1) matrices and vice versa.
#Both the matrix and the vector share the same pointer.
Matrix{T}(A::Vector{T}) where T = reshape(A,:,1)
Vector{T}(A::Matrix{T}) where T = reshape(A,:)
function Vector{<:AbstractMatrix{T}}(v::Vector{<:AbstractVector{T}}) where T
    n   = length(v)
    w   = Vector{Matrix{T}}(undef,n)
    for ii=1:n
        w[ii]   = reshape(v[ii],:,1)
    end
    return w
end

#Define generalized vector, that is, any array/subarray of dimension 1:
VectorView{T}       = SubArray{T,1,Vector{T},Tuple{UnitRange{Int64}},true} where T<:Any
MatrixCol{T}        = SubArray{T,1,Matrix{T},Tuple{Base.Slice{Base.OneTo{Int64}}, Int64},true} where T<:Any
MatrixColView{T}    = SubArray{T,1,Matrix{T},Tuple{UnitRange{Int64},Int64},true} where T<:Any
GenVector0{T}       = Union{VectorView{T}, MatrixColView{T}, 
                        MatrixCol{T}} where T<:Any
GenVector{T}        = Union{Vector{T}, GenVector0{T}} where T<:Any

function ViewType(v::Vector{T}) where T<:Any
    return VectorView{T}
end
function ViewType(v::VectorView{T}) where T<:Any
    return VectorView{T}
end

function Vector{T}(v::GenVector0{T}) where T
    return one(T)*v
end
function Matrix{T}(v::GenVector0{T}) where T
    return reshape(Vector{T}(v),:,1)
end

#Shorcuts for abstract arrays:
const AVF64         = AbstractVector{Float64}
const AMF64         = AbstractMatrix{Float64}
const AAF64         = AbstractArray{Float64}

#Convert from GenVector to Vector:
# Vector{T}(v::GenVector{T}) where T = one(T)*v
# function Vector{Vector{T}}(v::Vector{<:GenVector{T}}) where T
#     n   = length(v)
#     w   = Vector{Vector{T}}(undef,n)
#     for ii=1:n
#         w[ii]   = Vector{T}(v[ii])
#     end
#     return w
# end
# function Vector{Matrix{T}}(v::Vector{<:GenVector{T}}) where T
#     n   = length(v)
#     w   = Vector{Matrix{T}}(undef,n)
#     for ii=1:n
#         w[ii]   = reshape(Vector{T}(v[ii]),:,1)
#     end
#     return w
# end

#Alias for function that receives vector of matrices, [x[1], ..., x[M]], and returns matrix M:
const FW1           = FunctionWrapper{Matrix{Float64}, 
                        Tuple{Vector{<:AMF64}}}
#Alias for function that receives vector of matrices, [x[1], ..., x[M]], and returns 
#vector of matrices [f[1], ..., f[N]]:
const FW11          = FunctionWrapper{Vector{Matrix{Float64}}, 
                        Tuple{Vector{<:AMF64}}}
#Alias for function that receives Float and vector of matrices, [x[1], ..., x[M]], and returns 
#vector of matrices [f[1], ..., f[N]]:
const FWt11         = FunctionWrapper{Vector{Matrix{Float64}}, 
                        Tuple{Float64, Vector{<:AMF64}}}
                        
#Alias for function that receives Float and two vector of matrices, [x[1], ..., x[M]], [u[1], ..., u[N]] and returns 
#vector of matrices [f[1], ..., f[P]]:
const FWt21         = FunctionWrapper{
                        Vector{Matrix{Float64}}, 
                        Tuple{  Float64, 
                                Vector{<:AMF64},
                                Vector{<:AMF64} } }

"
    x   = linspace(a,b,n)

Generates a vector with 'n' equally spaced points between 'a' and 'b'.

"
function linspace(a::Float64,b::Float64,n::Int)

    return Array(range(a, stop=b, length=n))

end

"
    x   = logspace(a,b,n)

Generates a vector with 'n' logarithmically equally spaced points between 'a' and 'b'.

"
function logspace(a::Float64,b::Float64,n::Int)

    return exp.(linspace(log(a),log(b),n))

end

"

    f   = Heaviside(x::Float64)

Returns Heaviside ffunction
"
function Heaviside(x::Float64)

    if x<0.0
        return 0.0
    elseif x==0
        return 0.5
    else
        return 1.0
    end

end

function SmoothHeaviside(x::Float64, delta::Float64)

    return 0.5+0.5*tanh(x/delta)

end

function SmoothHeaviside(x::Float64, delta::Float64, 
    a::Float64, b::Float64)

    return (a+b)/2 + (b-a)/2*tanh(x/delta)

end

function SmoothSin(x::Float64, delta::Float64, 
    a::Float64, b::Float64)

    if x<-delta/2
        return a
    elseif x<delta/2
        return (a+b)/2 + (b-a)/2*sin(pi*x/delta)
    else    
        return b
    end

end

#Cumulative factorial [0!, 1!, ..., N!]
function cumfact(N::Int)

    return vcat(1,cumprod(1:N))

end

"

    xv, wv = GaussLegQuad(a::Float64, b::Float64, ngp::Int)

Returns nodes and weights of Gauss--Legendre quadrature in the interval [a,b].

"
function GaussLegQuad(a::Float64, b::Float64, ngp::Int)

    xiv, wv     = gausslegendre(ngp)
    xv          = @. 0.5*(a+b) + 0.5*(b-a)*xiv
    @. wv       *= 0.5*(b-a)
    return      xv, wv
    
end

"
    imat,jmat = ndgrid(iv,jv)

Given vectors 'iv' and 'jv', of lengths 'M' and 'N', returns matrices 'imat'
and 'jmat' of size 'M*N' so that imat[i,j]=iv[i], jmat[i,j]=jv[j]
"
function ndgrid(iv::AbstractVector{T}, jv::AbstractVector{T}) where {T<:Real}

    imat    = zeros(T, length(iv), length(jv))
    jmat    = zeros(T, length(iv), length(jv))
    for jj=1:length(jv), ii=1:length(iv)
        imat[ii,jj]     = iv[ii]
        jmat[ii,jj]     = jv[jj]
    end
    return imat, jmat
    
end

"
    A   = eye(N)

Returns identity matrix of dimension N*N
"
function eye(N::Int)

    return diagm(ones(N))

end

#Perform c = alpha*a*b + beta*c:
@inline function ew_mult!(alpha::Float64, a::AbstractVector{Float64}, b::AbstractVector{Float64}, beta::Float64, c::AbstractVector{Float64})

    @mlv    c   = alpha*a*b + beta*c
    
    return
    
end

#Perform y = alpha*x1*x2*x3 + beta*y:
@inline function ew_mult!(alpha::Float64, x1::AbstractVector{Float64}, x2::AbstractVector{Float64}, x3::AbstractVector{Float64}, beta::Float64, y::AbstractVector{Float64})

    @mlv    y   = alpha*x1*x2*x3 + beta*y
    
    return
    
end

#Chebyshev nodes of the second kind in a line:
function ChebyNodes(a::Float64, b::Float64, N::Int)

    if N<1
        error("N must be positive")
    elseif N==1
        return [0.0]
    else
        return @. 0.5*(a+b)+0.5*(a-b)*cos((0:N-1)*pi/(N-1))
    end
    
end

"
    Tm  = PolyLeg(xv, N; a=-1, b=1)

Return a matrix such that Tm[ii,jj+1] is the value of the Legendre polynomial
T_jj (traslated to the interval [a,b]) at xv[ii]. 
"
function PolyLeg(xv::Vector{Float64}, N::Int; a::Float64=-1.0, b::Float64=1.0)

    #Center x:
    xhatv   = @. (xv-0.5*(a+b))/(0.5*(b-a))

    #Allocate matrix:
    Nx      = length(xv)
    Tm      = zeros(Nx, N+1)

    #Apply recurrence:
    #   T0      = 1
    #   T1      = x
    #   T_(n+1) = ((2n+1)x*T_n - n*T_(n-1))/(n+1)
    Tm[:,1]     .= 1.0
    if N>0
        Tm[:,2]     .= xhatv
        for nn=1:N-1
            @views @. Tm[:,nn+2]    = ((2.0*nn+1)*xhatv*Tm[:,nn+1] - nn*Tm[:,nn])/(nn+1)
        end
    end

    return Tm

end

"
    dTm  = dPolyLeg(xv, N; a=-1, b=1)

Return a matrix such that dTm[ii,jj+1] is the value of the derivative of the Legendre polynomial
dT_jj (traslated to the interval [a,b]) at xv[ii]. 
"
function dPolyLeg(xv::Vector{Float64}, N::Int; a::Float64=-1.0, b::Float64=1.0)

    #Center x:
    xhatv   = @. (xv-0.5*(a+b))/(0.5*(b-a))

    #Allocate matrix:
    Nx      = length(xv)
    dTm     = zeros(Nx, N+1)
    Tm      = PolyLeg(xhatv, N)

    #Apply recurrence:
    #   dT_(n) = (n)*T_(n-1) + x*dT_(n-1)
    dTm[:,1]    .= 0.0
    for nn=1:N
        dTm[:,nn+1]     = @views @. nn*Tm[:,nn] + xhatv*dTm[:,nn]
    end

    #d/dx = dxhat/dx*d/dxhat
    return dTm*2.0/(b-a)

end

"
    colorv  = PyPlotColors(map::String, nColors::Int)

Generates array with nColors in the specified map"
function PyPlotColors(map::String,nColors::Int)

    if nColors==1
        if map=="jet2"
            cmap 		= PyPlot.get_cmap("jet")
            return [ cmap(0.1) ]
        else
            cmap        = PyPlot.get_cmap(map)
            return [ cmap(0.0) ]
        end
    end

	if map=="jet2"
		cmap 		= PyPlot.get_cmap("jet")
		iv          = linspace(0.1,0.9,nColors)
    	return [ cmap(ii) for ii in iv ]
	end
	
    cmap            = PyPlot.get_cmap(map)
    iv              = linspace(0.0,1.0,nColors)
    return [ cmap(ii) for ii in iv ]

end

#Open figure with given total size and given margins:
function PyPlotFigure(; w::Float64=7.50, h::Float64=7.50, left::Float64=1.3, right::Float64=0.4, bottom::Float64=1.2, top::Float64=0.25)

    #Open figure:
    fig1        = PyPlot.figure(figsize=(w/2.54,h/2.54))
    if true

        fwidth      = w
        fheight     = h

        #---- define margins -> size in cm / figure dimension ----

        left_margin  = left / fwidth
        right_margin = right / fwidth
        bottom_margin = bottom / fheight
        top_margin = top / fheight

        #---- create axes ----

        # dimensions are calculated relative to the figure size

        x = left_margin    # horiz. position of bottom-left corner
        y = bottom_margin  # vert. position of bottom-left corner
        w = 1 - (left_margin + right_margin) # width of axes
        h = 1 - (bottom_margin + top_margin) # height of axes

        ax = fig1.add_axes([x, y, w, h])
    end

    return fig1

end

function PyPlotSubPlots(m::Int, n::Int; w::Float64=7.50, h::Float64=7.50, left::Float64=1.3, right::Float64=0.4, bottom::Float64=1.2, top::Float64=0.25, 
    tTitle::Float64=0.0, bTitle::Float64=0.0)

    #Total figure dimensions:
    fwidth      = w*n
    fheight     = h*m + tTitle + bTitle
    
    #Create figure with subplots:
    fig1, ax1   = subplots(m, n, figsize=(fwidth/2.54,fheight/2.54))

    #---- create axes ----
    #left is the position of the first left axis normalized with the figure width,
    #right is the position of the last right axis normalized with the figure width,
    #idem with top and bottom
    #wspace is the horizontal padding between subplots, normalized with the average axes width.
    #hspace is the vertical padding between subplots, normalized with the average axes height.
    #
    #The axis height, say h', satisfies
    #   (h'+t+b)    = h,
    #   h'          = h - t - b
    #   hspace      = (t+b)/h'  = (t+b)/(h-t-b) 
    #   hspace      =(top_margin+bottom_margin)/(1.0-mFig*(top_margin+bottom_margin))
    subplots_adjust(
        left=left/fwidth, 
        right=1.0-right/fwidth, 
        bottom=(bTitle+bottom)/fheight, 
        top=1.0-(top+tTitle)/fheight, 
        hspace=(top+bottom)/(h-top-bottom),
        wspace=(left+right)/(w-left-right))

    return fig1

end

function PyPlotSubPlots_old(m::Int, n::Int; w::Float64=7.50, h::Float64=7.50, left::Float64=1.3, right::Float64=0.4, bottom::Float64=1.2, top::Float64=0.25)

    #Create figure with subplots:
    fig1, ax1   = subplots(m, n, figsize=(w*n/2.54,h*m/2.54))

    #-------------------------------------
    #Adjust subplots:

    fwidth      = w*n
    fheight     = h*m

    #---- define margins -> size in cm / figure dimension ----

    left_margin  = left / fwidth
    right_margin = right / fwidth
    bottom_margin = bottom / fheight
    top_margin = top / fheight

    #---- create axes ----
    #left is the position (normalized with the figure width) of the first left axis,
    #right is the position (idem) of the last right axis, and so on.
    #hspace is the space between succesive top and bottom axes, normalized with the axis height. The axis height, say h', satisfies
    #   mFig * (h'+t+b) = h,
    #   h'      = h/mFig - t - b
    #   hspace  = (t+b)/h'  = (t+b)/(h/mFig-t-b) 
    #   hspace=(top_margin+bottom_margin)/(1.0-mFig*(top_margin+bottom_margin))
    subplots_adjust(left=left_margin, right=1.0-right_margin, bottom=bottom_margin, top=1.0-top_margin,
        hspace=(top_margin+bottom_margin) / (1.0/m - top_margin - bottom_margin), 
        wspace=(left_margin+right_margin) / (1.0/n - left_margin - right_margin))

    return fig1

end

#Plot colorbar in active figure:
function PyPlotColorbar!(y1::Float64, y2::Float64, nTicks::Int;
    ycontour::Vector{Float64}=linspace(y1,y2,100), 
    cmap::String="jet", 
    labelsize::Number=10, ticksize::Number=8, tickformat::String=".2E")

    xm,ym       = ndgrid([0.0,1.0],ycontour)
    zm          = copy(ym)

    fig1        = gcf()
    plot1       = PyPlot.contourf(xm,ym,zm,ycontour,cmap="jet")
    PyPlot.tick_params(
        axis="x",          # changes apply to the x-axis
        which="both",      # both major and minor ticks are affected
        bottom=false,      # ticks along the bottom edge are off
        top=false,         # ticks along the top edge are off
        labelbottom=false) # labels along the bottom edge are off
    yticksv_1   = linspace(y1, y2, nTicks)
    yticksv_2   = Vector{String}(undef, nTicks)
    for ii=1:nTicks
        yticksv_2[ii]   = sprintf1(tickformat, yticksv_1[ii])
    end
    PyPlot.yticks(yticksv_1, yticksv_2, fontfamily="serif", 
        fontsize=labelsize)
    PyPlot.tick_params(axis="both",labelsize=ticksize)
    PyPlot.tight_layout()

end

#Print Profile in file:
function ProfilePrint(filename::String)

    #format=:free, :count
    file    = open(filename, "w")
    Profile.print(IOContext(file,:displaysize => (24, 2000)),format=:tree)#,sortedby=:count)
    close(file)

end

"
    Ainv    = BlockDiagInv_GaussJordan(A,BlockSize)

Returns the inverse 'Ainv' of a block diagonal matrix 'A' with blocks of size
'BlockSize' with Gauss-Jordan algorithm.

Use this algorithm only for matrices with a great number of small blocks, and with
blocks that can be factorized without pivoting.
"
function BlockDiagInv_GaussJordan(A::SparseMatrixCSC{Float64,Int},BlockSize::Int)

    t_ini       = time()
    
    #Check that matrix is square:
    N           = size(A,1)
    if N!=size(A,2)
        error("Matrix must be square")
    end

    #Get nb of blocks:
    nBlocks     = size(A,1)÷BlockSize  #write \div and press tab to write ÷
    if N%BlockSize!=0
        error("Dimension of A must be a multiple of BlockSize")
    end

    #Copy input matrix:
    Ainv        = deepcopy(A)

#     println("BDI - deepcopy, ", time()-t_ini)
    
    #Apply operations in A to get diagonal matrix. The zeros in the resulting
    #matrix are replaced by the numbers we obtain by applying the same operations in
    #the identity matrix:
    @inbounds for ii=1:BlockSize

        aux         = vcat(1:ii-1,ii+1:BlockSize)
#         println("BDI - loop0, ", time()-t_ini)
        
        #Modify term (ii,ii):
        @views @mlv Ainv.nzval[ii+(ii-1)*BlockSize:BlockSize^2:end]       =
            1.0 ./ Ainv.nzval[ii+(ii-1)*BlockSize:BlockSize^2:end]

#         println("BDI - loop1, ", time()-t_ini)
        
        #Modify column ii:
        for kk in aux
            @views @mlv Ainv.nzval[kk+(ii-1)*BlockSize:BlockSize^2:end]   =
                - Ainv.nzval[kk+(ii-1)*BlockSize:BlockSize^2:end] *
                Ainv.nzval[ii+(ii-1)*BlockSize:BlockSize^2:end]
        end

#         println("BDI - loop2, ", time()-t_ini)
        
        #Modify rest of columns, except row ii:
        for kk in aux, jj in aux

#             @views @mlv Ainv.nzval[kk+(jj-1)*BlockSize:BlockSize^2:end]   =
#                 Ainv.nzval[kk+(jj-1)*BlockSize:BlockSize^2:end] +
#                 Ainv.nzval[kk+(ii-1)*BlockSize:BlockSize^2:end] *
#                 Ainv.nzval[ii+(jj-1)*BlockSize:BlockSize^2:end]
#             @views @. Ainv.nzval[kk+(jj-1)*BlockSize:BlockSize^2:kk+(jj-1)*BlockSize+(nBlocks-1)*BlockSize^2]   +=
#                 Ainv.nzval[kk+(ii-1)*BlockSize:BlockSize^2:kk+(ii-1)*BlockSize+(nBlocks-1)*BlockSize^2] *
#                 Ainv.nzval[ii+(jj-1)*BlockSize:BlockSize^2:ii+(jj-1)*BlockSize+(nBlocks-1)*BlockSize^2]
            for ll=1:nBlocks
                Ainv.nzval[kk+(jj-1)*BlockSize+(ll-1)*BlockSize*BlockSize]     += 
                    Ainv.nzval[kk+(ii-1)*BlockSize+(ll-1)*BlockSize*BlockSize] *
                    Ainv.nzval[ii+(jj-1)*BlockSize+(ll-1)*BlockSize*BlockSize]
            end
            
        end

#         println("BDI - loop3, ", time()-t_ini)
        
        #Modify row ii:
        for jj in aux
            @views @mlv Ainv.nzval[ii+(jj-1)*BlockSize: BlockSize^2:end]  =
                Ainv.nzval[ii+(jj-1)*BlockSize:BlockSize^2:end] *
                Ainv.nzval[ii+(ii-1)*BlockSize:BlockSize^2:end]
        end

#         println("BDI - loop4, ", time()-t_ini)
        
    end

#     println("BDI - loop, ", time()-t_ini)
    
    #Return:
    return Ainv

end

"
    ColorInfo   = SpColors(A)

Returns tuple (colorv, nColors, ColorCols) with information of the coloring of
matrix 'A'. 'colorv' indicates the color each column of A belongs to, 'nColors'
indicates the total number of colors and 'ColorCols[ii]' is an array that indicates
the columns that belong to color 'ii'.
"
function SpColors(A::SparseMatrixCSC)

    colorv      = matrix_colors(A)
    nColors     = maximum(colorv)
    ColorCols   = Vector{Vector{Int}}(undef,nColors)
    for ii=1:nColors

        #Count number of columns that belong to color "ii":
        nCols           = 0
        for jj=1:size(A,2)
            if colorv[jj]==ii
                nCols   += 1
            end
        end

        #Save columns that belong to color "ii":
        ColorCols[ii]   = zeros(Int,nCols)
        kk              = 1
        for jj=1:size(A,2)
            if colorv[jj]==ii
                ColorCols[ii][kk]   = jj
                kk                  += 1
            end
        end

    end

    #Return output:
    return colorv, nColors, ColorCols

end

"
    SpJacobEst!(J,f_fun!,x0,f0,ColorInfo=nothing,deltaxv=fill(1e-4,size(x0,1)))

Computes sparse jacobian estimation with finite differences. Overrides 'J'. 'J'
must have been previously allocated. 'f_fun!' is of the form f_fun!(f,x)"
function SpJacobEst!(J::SparseMatrixCSC{Float64,Int}, f_fun!::Function, x0::Array{Float64,1}, 
    f0::Union{Array{Float64,1},Nothing};
    ColorInfo=nothing, deltaxv::Array{Float64,1}=fill(1e-8,size(x0,1)))

    #If ColorInfo is not specified, compute it:
    if ColorInfo==nothing
        ColorInfo   = SpColors(J)
    end

    #Extract color info:
    colorv          = ColorInfo[1]
    nColors         = ColorInfo[2]
    ColorCols       = ColorInfo[3]

    #Loop colors:
    fpert           = zeros(size(J,1))
    for iColor=1:nColors

        #Related columns:
        cols        = ColorCols[iColor]

        #Perturb x:
        xpert       = copy(x0)
        xpert[cols] .= xpert[cols].+deltaxv[cols]

        #Evaluate perturbed function and derivatives:
        f_fun!(fpert,xpert)
        Deltaf      = fpert.-f0

        #Sweep jacobian matrix:
        for jj=cols
            for ss=J.colptr[jj]:(J.colptr[jj+1]-1)
                ii              = J.rowval[ss]
                J.nzval[ss]     = Deltaf[ii]/deltaxv[jj]
            end
        end

    end

    return
    
end

function SpJacobEst!(J::SparseMatrixCSC{Float64,Int}, f_fun!::Function, x0::Array{Float64,1};
    ColorInfo=nothing, deltaxv::Array{Float64,1}=fill(1e-8,size(x0,1)))

    #If ColorInfo is not specified, compute it:
    if ColorInfo==nothing
        ColorInfo   = SpColors(J)
    end

    #Extract color info:
    colorv          = ColorInfo[1]
    nColors         = ColorInfo[2]
    ColorCols       = ColorInfo[3]

    #Loop colors:
    fpert1          = zeros(size(J,1))
    fpert2          = zeros(size(J,1))
    for iColor=1:nColors

        #Related columns:
        cols        = ColorCols[iColor]

        #Perturb x:
        xpert1      = copy(x0)
        xpert1[cols].= xpert1[cols].-deltaxv[cols]
        xpert2      = copy(x0)
        xpert2[cols].= xpert2[cols].+deltaxv[cols]

        #Evaluate perturbed function and derivatives:
        f_fun!(fpert1,xpert1)
        f_fun!(fpert2,xpert2)
        Deltaf      = fpert2.-fpert1

        #Sweep jacobian matrix:
        for jj=cols
            for ss=J.colptr[jj]:(J.colptr[jj+1]-1)
                ii              = J.rowval[ss]
                J.nzval[ss]     = Deltaf[ii]/(2.0*deltaxv[jj])
            end
        end

    end

    return
    
end

function free!(v::Array{<:AbstractArray{T}}) where T<:Number

#     return
    
    for ii=1:length(v)
        v[ii]   = zeros(T, ntuple(_->0, length(size(v[ii]))))
    end
    
    #DO NOT USE GC.gc(). It is too slow if free is used many times.

    return
    
end


"
    sparse!(Sm,iv,jv,sv).

Writes into previously allocated sparse matrix. 'iv', 'jv', 'sv' can be matrices.
"
function sparse!(Sm::SparseMatrixCSC{N,Int}, iv::AbstractArray{Int},
    jv::AbstractArray{Int}, sv::AbstractArray{N}) where N<:Number

    #Nb of nnz in iv, jv, sv:
    nnz     = length(iv)
    @assert length(jv)==length(sv)==nnz

    #Check dimensions:
    if maximum(iv)>size(Sm,1) ||  maximum(jv)>size(Sm,2)
        error("(i,j) must be smaller than matrix size")
    end
    
    #Zero Sm:
    Sm.nzval    .= zero(N)
    #Loop nnz:
    @inbounds for kk=1:nnz

        #Extract values:
        row     = iv[kk]
        col     = jv[kk]
        val     = sv[kk]

        #Search row in current column:
        ss      = Sm.colptr[col]
        while true
            if Sm.rowval[ss]==row
                break
            end
            ss  += 1
            if ss== Sm.colptr[col+1]
                error("Element ($row,$col) is not allocated in previous matrix")
            end
        end

        #Write in element ss of Sm:
        Sm.nzval[ss]    += val

    end
end

#Direct construction of sparse matrix from AbstractArray:
import SparseArrays.sparse
function sparse(iv::Array{Int,T}, jv::Array{Int,T}, 
    sv::Array{N,T}) where {N<:Number, T<:Any}
    #do not use AbstractArray; otherwise, this function will call itself
    
    return sparse(view(iv,:), view(jv,:), view(sv,:))
    
end

#Set all arrays in Array{<:AbstractArray}} to zero:
function zero!(v::Array{<:AbstractArray{T}}) where T<:Number
    
    for ii=1:length(v)
        v0              = zero(T)
        @mlv    v[ii]   = v0
    end
    
    return
    
end

function zero!(v::Array{<:AbstractArray{T}}) where T<:Number
    
    for ii=1:length(v)
        v0              = zero(T)
        @mlv    v[ii]   = v0
    end
    
    return
    
end

function alloc!(v::Array{MNumber}, dims::NTuple{n,Int}) where {n<:Any, 
    N<:Number, MNumber<:AbstractArray{N}}

    for ii=1:length(v)
        v[ii]       = zeros(N, dims)
    end
    return
    
end

import PyPlot.plot3D

function plot3D(S::SparseMatrixCSC{Float64,Int}; opts=".b")
    
    i,j,s   = findnz(S)
    plot3D(i,j,s, opts)
    
    return
    
end

import PyPlot.plot

function plot(S::SparseMatrixCSC{Float64,Int}; opts=".b")
    
    plot(S.nzval, opts)
    
    return
    
end

function PrintToFile(filename::String, v::Any)

    open(filename,"w") do io
        println(io, v)
    end

end

function AppendToFile(filename::String, v::Any)

    open(filename,"a") do io
        println(io, v)
    end

end

#Upper bound for L2 norm of sparse matrix:
function opnorm2(S::SparseMatrixCSC{N,Int}) where N<:Number
    return sqrt(opnorm(S,1)*opnorm(S,Inf))
end

#Spectral decomposition of a symmetric matrix V*lambda*V^T of dimension d:
mutable struct SpectralDecomp{T<:Union{Float64, Vector{Float64}}, d<:Int}

    d           ::Int           #dimension
    lambda_bar  ::T              #"geometric mean" Mbar=(prod_i lambda_i)^(1/d)=det(M)^(1/d)
    sf          ::Vector{T}     #stretching factors, s1>s2>..., s1*s2*...=1
    V           ::Matrix{T}     #matrix with eigenvectors [v1,v2,...]
   
    SpectralDecomp{T,d}() where{T<:Union{Float64, Vector{Float64}}, d<:Int} = 
        new{T,d}()
    
    #Recall that
    #   A   = V*diag(lambda_i)*V^T = lambda_bar V diag(sf_i) V^T
    #where
    #   lambda_bar  = (prod_i lambda_i)^(1/d)
    #is the algebraic mean of the eigenvalues,
    #   sf_i        = lambda_i/(prod_i lambda_i)^(1/d)  
    #are the stretching factors.
    #lambda_i, lambda_bar and A_ij have the same dimensionss, s_i is dimensionless.
    #
    #Note that
    #   |A| = lambda_bar^d = prod_i lambda_i
    
end

function SpectralDecomp(N::Int, d::Int)

    A       = SpectralDecomp{Vector{Float64}, Int}()
    A.d     = d
    A.lambda_bar    = zeros(N)
    A.sf            = Vector{Vector{Float64}}(undef,d)
    alloc!(A.sf, (N,))
    A.V             = Matrix{Vector{Float64}}(undef,d,d)
    alloc!(A.V, (N,))
    return A
    
end

function GetMatrix(A::SpectralDecomp{Vector{Float64},Int}; 
    lambda_bar::Union{Float64, Vector{Float64}}=A.lambda_bar)

    Am  = Matrix{Vector{Float64}}(undef, A.d, A.d)
    alloc!(Am, size(A.sf[1]))
    for ii=1:A.d, jj=1:A.d
        for kk=1:A.d
            @mlv Am[ii,jj]  += A.V[ii,kk]*A.sf[kk]*A.V[jj,kk]
        end
        @mlv Am[ii,jj] *= lambda_bar
    end
    return Am
    
end

"
    msg     = Error2String(err::ErrorException)

Converts error exception to string.
https://discourse.julialang.org/t/is-it-possible-to-get-an-exception-message-as-a-string/3201

This is better than simply calling string(err)
"
@inline function Error2String(err::Any)    ::String
    bt      = catch_backtrace()
    msg     = sprint(showerror, err, bt)
    return  msg
end

#Auxiliary vectors with undefined length, useful for loops:
mutable struct IncreasingVector{T<:Number}

    k       ::Int
    N       ::Int
    v       ::Vector{T}
    DeltaN  ::Int
    IncreasingVector{T}() where T<:Number = new{T}()
    
end

function IncreasingVector(T::DataType; N::Int=1000, DeltaN::Int=N)

    v           = IncreasingVector{T}()
    v.k         = 0
    v.N         = N
    v.v         = zeros(T,N)
    v.DeltaN    = DeltaN
    return v
    
end

function expand!(v::IncreasingVector{T}, vi::T) where T<:Number

    v.k         += 1
    if v.k>v.N
        display(v.v)
        display(v.N)
        display(v.DeltaN)
        error("")
        #Realloc memory:
        v.v     = vcat(v.v, zeros(T,v.DeltaN))
        v.N     += v.DeltaN
    end
    v.v[v.k]    = vi
    return
    
end

function zero!(v::IncreasingVector{T}) where T<:Number
    v.k         = 0
    @mlv v.v    = zero(T)
    return
end

function chunk(v::IncreasingVector{<:Number})
    return v.v[1:v.k]
end

function SecondDegreeEquation(a::Float64, b::Float64, c::Float64)
    x1  = (-b+sqrt(b^2-4*a*c))/(2*a)
    x2  = (-b-sqrt(b^2-4*a*c))/(2*a)
    return min(x1,x2), max(x1,x2)
end

function SpBlockDiag(A::Matrix{T}, nBlocks::Int) where T<:Any

    m,n     = size(A)
    iv      = zeros(Int,m,n,nBlocks)
    jv      = zeros(Int,m,n,nBlocks)
    sv      = zeros(T,m,n,nBlocks)
    @inbounds for II=1:nBlocks
        i0      = (II-1)*m
        j0      = (II-1)*n
        for jj=1:n, ii=1:m
            iv[ii,jj,II]    = i0+ii
            jv[ii,jj,II]    = j0+jj
            sv[ii,jj,II]    = A[ii,jj]
        end
    end
    return sparse(reshape(iv,:), reshape(jv,:), reshape(sv,:))
    
end

#Create (empty) sparse matrix S and return index vector pinv such that 
#   S.nzval[pinv[kk]]   += smat[kk]
#where smat is an array with the values corresponding to the rows 
#and columns in iv, jv.
function SpAlloc(T::DataType, imat::Array{Int,N}, jmat::Array{Int,N}; m=maximum(imat), n=maximum(jmat)) where N

    #Reshape:
    iv              = reshape(imat,:)
    jv              = reshape(jmat,:)
    
    #Sort by columns and then by rows: k=(j-1)*nDof+i
    kv              = @tturbo @. (jv-1)*m + iv
    p               = sortperm(kv)
    kv_sorted       = kv[p]
    iv_sorted       = iv[p]
    jv_sorted       = jv[p]
    
    #Currently, kv is sorted, but not numbered consecutively. #Renumber from 1 to nnz,
    #with nnz the number of nnz in the sparse matrix:
    ssv             = zeros(Int, length(kv))
    last_j          = 0
    last_k          = 0
    nnz             = 0
    colptr          = zeros(Int, n+1)
    rowval          = zeros(Int, length(kv))
    @inbounds for ss=1:length(kv_sorted)
        if jv_sorted[ss] != last_j
            #We are in a different column. 
            last_j          = jv_sorted[ss]
            #Note that "nnz" has not been updated yet:
            colptr[last_j]  = nnz+1
        end
        if kv_sorted[ss] != last_k
            #We are in a different (i,j) pair:
            last_k          = kv_sorted[ss]
            nnz             += 1 
            rowval[nnz]     = iv_sorted[ss]
        end
        ssv[ss]             = nnz
    end
    colptr[n+1]             = nnz+1
    rowval                  = rowval[1:nnz]
    
    #Note that S.nzval[ssv[ii]] += sv[p[ii]], that is, component 
    #p[ii] in sv updates component ssv[ii] in S.nzval.
    #
    #We seek vector pinv such that component jj in sv updates 
    #component pinv[jj] in S.nzval. That is,
    #   S.nzval[pinv[jj]]       += sv[jj]
    #Note that, if we make jj=p[ii], then
    #   S.nzval[pinv[p[ii]]]    += sv[p[ii]]
    #On the other hand, we know that
    #   S.nzval[ssv[ii]]        += sv[p[ii]]
    #so
    #   pinv[p]                 = ssv
    #(It is also possible to demonstrate this by "drawing" the vectors)
    pinv            = zeros(Int, length(kv))
    pinv[p]         .= ssv
    
    #Create sparse matrix:
    nzval           = zeros(T, nnz)
    S               = SparseMatrixCSC{T,Int}(m, n, colptr, rowval, nzval)
    
    return S, pinv

end

#Update sparse matrix S according to
#   S.nzval                 = S.nzval * alpha
#   S.nzval[pinv[ii]]       += smat[ii]
#DO NOT USE this funtion to update a matrix in a block-wise function. The sentence 
#S.nzval *= alpha will be executed too many times.
function sparse!(alpha::N, S::SparseMatrixCSC{N,Int}, #pinv::Vector{Int}, 
    pinv::Array{Int,T1}, 
    smat::Array{N,T2}) where{T1<:Any, T2<:Any, N<:Number}

    S_nzval                 = S.nzval
    @avxt @. S_nzval        *= alpha
    @inbounds for ii=1:length(smat)
        S_nzval[pinv[ii]]   += smat[ii]
    end
    
    return
    
end

#Update sparse matrix S according to
#   S.nzval[pinv[ii]]       += smat[ii]
function sparse!(S::SparseMatrixCSC{N,Int}, 
    pinv::Array{Int,T1}, 
    smat::Array{N,T2}) where{T1<:Any, T2<:Any, N<:Number}

    S_nzval                 = S.nzval
    @inbounds for ii=1:length(smat)
        S_nzval[pinv[ii]]   += smat[ii]
    end
    
    return
    
end

#Keep sprintf1 from Formatting:
function sprintf1(format::String, x::Any)

    return cfmt(format, x)
    
end
