include("basic.jl")

#flags:
#   -3: unfinished
#   -2: maximum of line searches reached
#   -1: fun! threw an exception
#   0: maximum of iterations reached
#   1: step p under tolerance
#   2: f under tolerance
#   3: g under tolerance

Base.@kwdef mutable struct NLS_ConvHistory

    nIter       ::Int               = 0
    fEvals      ::Int               = 0
    pnorms      ::Vector{Float64}   = zeros(0)
    fnorms      ::Vector{Float64}   = zeros(0)
    gnorms      ::Vector{Float64}   = zeros(0)
    flag        ::Int               = -3
    
end

#Function wrapper for NonlinearSolvers.jl
const FW_NLS        = FunctionWrapper{ 
                        Int, 
                        Tuple{  Vector{Float64}, 
                                Vector{Float64}, 
                                Union{Nothing,Vector{Float64}} } }
                            
const FW_NLS_norm   = FunctionWrapper{
                        Float64,
                        Tuple{  Vector{Float64} } }
                        
#Function to display info:
function Anderson_Display(Display::String, ch::NLS_ConvHistory, nSigma::Int, xnorm::Float64, pnorm::Float64, 
    f0norm::Float64, fnorm::Float64, g0norm::Float64, gnorm::Float64, CPUtime::Float64)

    lcDisplay   = lowercase(Display)
    if lcDisplay=="none"
        return
    end
    
    nIter       = ch.nIter
    fEvals      = ch.fEvals
    iter_info   = string("nIter=", nIter, ", fEvals=", fEvals, ", nSigma=", nSigma, 
                    ", |x|=", sprintf1("%.2e", xnorm), ", |p|=", sprintf1("%.2e", pnorm), 
                    ", |f0|=", sprintf1("%.2e", f0norm), ", |f|=", sprintf1("%.2e", fnorm),
                    ", |g0|=", sprintf1("%.2e", g0norm), ", |g|=", sprintf1("%.2e", gnorm), 
                    ", CPU time=", sprintf1("%.2e", CPUtime), " sec.", 
                    "\n") 
    if ch.flag==-3
        #Iteration in progress:
        if any(lcDisplay.==["iter"])
            printstyled(iter_info, color=:white) 
        end
    elseif ch.flag==-2
        #Max of line searches
        if any(lcDisplay.==["notify", "final", "iter"])
            printstyled("Maximum number of line searches reached; ", iter_info, color=:yellow)
        end
    elseif ch.flag==-1
        #fun! evaluation failed:
        if any(lcDisplay.==["notify", "final", "iter"])
            printstyled("Function evaluation failed; ", iter_info, color=:yellow)
        end
    elseif ch.flag==0
        #Max of iterations reached
        if any(lcDisplay.==["notify", "final", "iter"])
            printstyled("Maximum number of iterations reached; ", iter_info, color=:yellow)
        end
    elseif ch.flag==1
        if any(lcDisplay.==["final", "iter"])
            printstyled("Step p under tolerance; ", iter_info, color=:cyan)
        end
    elseif ch.flag==2
        if any(lcDisplay.==["final", "iter"])
            printstyled("Residual f under tolerance; ", iter_info, color=:cyan)
        end
    elseif ch.flag==3
        if any(lcDisplay.==["final", "iter"])
            printstyled("Residual g under tolerance; ", iter_info, color=:cyan)
        end
    end
    
    return
    
end

function NLS_ch_update!(ch::NLS_ConvHistory, history::Bool, pnorm::Float64, 
    fnorm::Float64, gnorm::Float64)

    if history
        ch.pnorms   = vcat(ch.pnorms, pnorm)
        ch.fnorms   = vcat(ch.fnorms, fnorm)
        ch.gnorms   = vcat(ch.gnorms, gnorm)
    end
    
    return
    
end
    
"
    Anderson(fun!, x0)
    
Anderson method (written iin multisecant form) to solve f(x)=0. fun! is of the form
    flag=fun!(x::Vector{Float64}, f::Vector{Float64}, g::Union{Nothing,Vector{Float64}})
i.e., receives vector x, and overwrites residual f, term g=H*f (when Hf!=nothing) , where H
is an approximation to the Jacobian inverse. The initial condition ffor g is saved iin g.
Also, returns flag>0 iff fun! found no errors, flag<=0 otherwise.

In particular, we apply Anderson's method to the quasi-Newton iteration
    x   = x - g(x), 
with g(x)=H(x)*f(x). We implement this iin multisecant form so as to be able to ddo line searches 
iin f(x).
    
-fun!: the ffunction above

-x0: Initial guess.

-options:

    --memory[100]: maximum of solutions considered

    --RelTolX[0e-6], AbsTolX[1e-6]: the method stops when norm(Delta x) < norm(x)*RelTolX + AbsTolX

    --RelTolF[0e-6], AbsTolF[0e-6]: the method stops when norm(f) < norm(f0)*RelTolF + AbsTolF
    
    --RelTolG[0e-6], AbsTolG[0e-6]: the method stops when norm(g) < norm(g0)*RelTolG + AbsTolG
    
    --MaxIter[length(x0)]: maximum of iterations.

    --Restart[MaxIter+1]: number of iterations after which the method is restarted.

    --Display['notify']: verbosity. Other options are 'none', 'iter', 'final'.

    --LineSearch[0]: maximum of line searches. If 0, no line searches are performed.

    --history[false]: log convergence history
    
    --sigma_min[1e-10]: drop tolerance after SVD decomposition
    
    --NormFun: FunctionWrapper that receives vector x and returns the norm of x.

Output:

-x: solution.

-ch: convergence history. See NLS_ConvHistory

"
function Anderson(fun!::FW_NLS, x0::GenVector{Float64};
    memory::Int=100, 
    RelTolX::Float64=0e-6, AbsTolX::Float64=1e-6,
    RelTolF::Float64=0e-6, AbsTolF::Float64=0e-6, 
    RelTolG::Float64=0e-6, AbsTolG::Float64=0e-6, 
    MaxIter::Int=length(x0), Restart::Int=MaxIter+1, 
    history::Bool=false, Display::String="notify",
    sigma_min::Float64=1e-10,
    NormFun::FW_NLS_norm=FW_NLS_norm((x)->norm(x)))
    
    #-----------------------------------------------------------------------------
    #ALLOCATE VARIABLES:
    
    t_start         = time()
    
    #Initialize variables:
    Nx              = size(x0,1)
    Iter2Restart    = Restart
    ch              = NLS_ConvHistory()
    
    #Iteration 0:
    xn              = copy(x0)
    fn              = zeros(Nx)
    gn              = zeros(Nx)
    FunFlag         = fun!(xn, fn, gn)
    pn              = @mlv -gn         #Tentative Delta x
    
    #Compute norms:
    xnorm           = NormFun(xn)
    pnorm           = NormFun(pn)
    fnorm           = NormFun(fn)
    gnorm           = NormFun(gn)
    fnorm0          = fnorm
    gnorm0          = gnorm
    
    #Update convergence history:
    ch.fEvals       += 1
    function _NLS_ch_update!()
        return NLS_ch_update!(ch, history, pnorm, fnorm, gnorm)
    end
    _NLS_ch_update!()
    
    #Display function:
    nSigma          = 0     #number of nonzero singular values     
    function _NLS_Display()
        return Anderson_Display(Display, ch, nSigma, xnorm, pnorm, fnorm0, fnorm, gnorm0, gnorm, time()-t_start)
    end
    
    #Leave immediatly if function evaluation failed:
    if FunFlag<=0
        ch.flag     = -1
        _NLS_Display()
        return xn, ch
    end
    
    #Initialize matrices X (contain Deltax's) and G (contain Deltag's):
    Xm              = zeros(Nx, memory)
#     Fm              = zeros(Nx, memory)
    Gm              = zeros(Nx, memory)
    XGcols          = 0
    
    #NOTE: At the end of each iteration k, the variables xn, fn, gn, pn, the corresponding norms
    #X and G have been computed.
    
    #Allocate next iteration:
    xnp1            = zeros(Nx)
    fnp1            = zeros(Nx)
    gnp1            = zeros(Nx)
  
    #Iterate:
    while true

        #Display iter info:
        _NLS_Display()
        
        #Exit:
        if pnorm <= RelTolX*xnorm+AbsTolX
            ch.flag     = 1
            break
        elseif fnorm <= RelTolF*fnorm0+AbsTolF
            ch.flag     = 2
            break
        elseif gnorm <= RelTolG*gnorm0+AbsTolG
            ch.flag     = 3
            break
        elseif ch.nIter==MaxIter
            ch.flag     = 0
            break
        end
        
        #Next iter:
        ch.nIter        += 1
        Iter2Restart    -= 1
        
        #Find next x:
        @mlv xnp1       = xn + pn   #Tentative step before line search
        @mlv gnp1       = gn        #Initial condition for gnp1:
        FunFlag         = fun!(xnp1, fnp1, gnp1)
        ch.fEvals       += 1
        if FunFlag<=0
            ch.flag     = -1
            break
        end
        #(Line search should be done here)
        
        #Update matrices X and F:
        if Iter2Restart==0 || memory==0
        
            #Reset X ans F:
            @mlv Xm         = 0.0
#             @mlv Fm         = 0.0
            @mlv Gm         = 0.0
            XGcols          = 0
            
            Iter2Restart    = Restart
            
        else 
        
            #Move columns:
            for jj=min(XGcols,memory-1):-1:1
                BLAS.copy!(view(Xm,:,jj+1), view(Xm,:,jj))
#                 BLAS.copy!(view(Fm,:,jj+1), view(Fm,:,jj))
                BLAS.copy!(view(Gm,:,jj+1), view(Gm,:,jj))
            end
            
            #New columns:
            @mlv Gm[:,1]    = gnp1-gn
            Deltag_norm_inv = 1.0/(NormFun(view(Gm,:,1))+1e-10)
            @mlv Gm[:,1]    *= Deltag_norm_inv
            @mlv Xm[:,1]    = (xnp1-xn)*Deltag_norm_inv
#             @mlv Fm[:,1]    = (fnp1-fn)*Deltag_norm_inv
            XGcols          = min(memory, XGcols+1)
            
            #Normalize new columns:
            BLAS.scal!(Deltag_norm_inv, view(Xm,:,1))
#             BLAS.scal!(Deltag_norm_inv, view(Fm,:,1))
            BLAS.scal!(Deltag_norm_inv, view(Gm,:,1))
            
        end
        
        #Update variables:
        @mlv xn         = xnp1
        @mlv fn         = fnp1
        @mlv gn         = gnp1
        
        #Next step. Apply multisecant formula:
        #   pn  = - H g^n = - [I + (X-G) (G^T G)^{-1} G^T] g^n 
        #       = -g^n - X*gamma + G*gamma,
        #with gamma = (G^T G)^{-1} G^T g^n 
        #If G=U Sigma V^T, then
        #   gamma = V Sigma^{-1} U^T g^n
        #Note that Anderson acceleration would read
        #   x_star      = xn - X*gamma
        #   g_star      = gn - G*gamma
        #   x_(n+1)     = x_star - g_star = xn - X*gamma - (gn - G*gamma)
        #hence
        #   p_n         = - gn - X*gamma + G*gamma
        #
        #Similarly, instead of minimizing the residual "g", we could minimize the residual
        #f. In that case, the multisecant update reads
        #   H = H_0 - (X-G) (F^T F)^{-1} F^T
        #where G = H_0^{1} F. If H_0 is not constant, the matrix G is not H_0^{-1} F and the formula
        #is not correct.
        Gm_SVD              = svd(view(Gm,:,1:XGcols))
        Sinv                = Gm_SVD.S
        deg_bool            = Sinv.<sigma_min
        Sinv[deg_bool]      .= 0.0
        Sinv[.!deg_bool]    .= 1.0./Sinv[.!deg_bool]
        gammav              = Gm_SVD.V*(Sinv.*(transpose(Gm_SVD.U)*gn))
        BLAS.axpby!(-1.0, gn, 0.0, pn)
        BLAS.gemv!('N', -1.0, view(Xm,:,1:XGcols), gammav, 1.0, pn)
        BLAS.gemv!('N', 1.0, view(Gm,:,1:XGcols), gammav, 1.0, pn)
        nSigma              = sum(.!deg_bool)
        
#         Anderson_Display(Display, ch, nSigma, norm(xn), norm(-Xm[:,1:XGcols]*gammav), 
#                 fnorm0, NaN, gnorm0, norm(gn-Gm[:,1:XGcols]*gammav))
                
        #Update norms:
        pnorm           = NormFun(pn)
        fnorm           = NormFun(fn)
        gnorm           = NormFun(gn)
        _NLS_ch_update!()
        
    end
    
    _NLS_Display()
    
    return xn, ch

end

"
    NLS_gmres(fun!, x0)
    
GMRES method to solve a linear equation of the form f(x)=0. fun! is of the form
    flag=fun!(x::Vector{Float64}, f::Vector{Float64}, g::Union{Nothing,Vector{Float64}})
i.e., receives vector x, and overwrites residual f, term g=H*f (when Hf!=nothing) , where H
is an approximation to the Jacobian inverse. The initial condition ffor g is saved iin g. Also, returns flag>0 iff fun! found no errors, flag<=0 
otherwise.

In particular, we apply an acceleration method to the quasi-Newton iteration
    x   = x - g(x), 
with g(x)=H(x)*f(x).

-fun!: the ffunction above

-x0: Initial guess.

-options:

    --memory[100]: maximum of solutions considered

    --RelTolX[0e-6], AbsTolX[1e-6]: the method stops when norm(Delta x) < norm(x)*RelTolX + AbsTolX

    --RelTolF[0e-6], AbsTolF[0e-6]: the method stops when norm(f) < norm(f0)*RelTolF + AbsTolF
    
    --RelTolG[0e-6], AbsTolG[0e-6]: the method stops when norm(g) < norm(g0)*RelTolG + AbsTolG
    
    --MaxIter[length(x0)]: maximum of iterations.

    --Restart[MaxIter+1]: number of iterations after which the method is restarted.

    --Display['notify']: verbosity. Other options are 'none', 'iter', 'final'.

    --LineSearch[0]: maximum of line searches. If 0, no line searches are performed.

    --history[false]: log convergence history
    
    --sigma_min[1e-10]: drop tolerance after SVD decomposition
    
    --NormFun: FunctionWrapper that receives vector x and returns the norm of x.

Output:

-x: solution.

-ch: convergence history. See NLS_ConvHistory

"
function NLS_gmres(fun!::FW_NLS, x0::GenVector{Float64};
    memory::Int=100, 
    RelTolX::Float64=0e-6, AbsTolX::Float64=1e-6,
    RelTolF::Float64=0e-6, AbsTolF::Float64=0e-6, 
    RelTolG::Float64=0e-6, AbsTolG::Float64=0e-6, 
    MaxIter::Int=length(x0), Restart::Int=MaxIter+1, 
    history::Bool=false, Display::String="notify",
    sigma_min::Float64=1e-10,
    NormFun::FW_NLS_norm=FW_NLS_norm((x)->norm(x)))
    
    #-----------------------------------------------------------------------------
    #ALLOCATE VARIABLES:
    
    t_start         = time()
    
    #PROFILE:
    t_startup       = 0.0
    t_residual      = 0.0
    t_SVD           = 0.0
    t_Matrices      = 0.0
    
    #Initialize variables:
    Nx              = size(x0,1)
    Iter2Restart    = Restart
    ch              = NLS_ConvHistory()
    
    #Iteration 0:
    t_ini           = time()
    xn              = copy(x0)
    fn              = zeros(Nx)
    gn              = zeros(Nx)
    FunFlag         = fun!(xn, fn, gn)
    pn              = fill(Inf, Nx)
    t_residual      += time()-t_ini
    
    #Compute norms:
    xnorm           = NormFun(xn)
    pnorm           = NormFun(pn)
    fnorm           = NormFun(fn)
    gnorm           = NormFun(gn)
    fnorm0          = fnorm
    gnorm0          = gnorm
    
    #Update convergence history:
    ch.fEvals       += 1
    function _NLS_ch_update!()
        return NLS_ch_update!(ch, history, pnorm, fnorm, gnorm)
    end
    _NLS_ch_update!()
    
    #Display function:
    nSigma          = 0     #number of nonzero singular values     
    function _NLS_Display()
        return Anderson_Display(Display, ch, nSigma, xnorm, pnorm, fnorm0, fnorm, gnorm0, gnorm, time()-t_start)
    end
    
    #Leave immediatly if function evaluation failed:
    if FunFlag<=0
        ch.flag     = -1
        _NLS_Display()
        return xn, ch
    end
    
    #Initialize matrices X (contain Deltax's) and G (contain Deltag's):
    t_ini           = time()
    Xm              = zeros(Nx, memory)
    Fm              = zeros(Nx, memory)
    Gm              = zeros(Nx, memory)
    XGcols          = 0
    
    #NOTE: At the end of each iteration k, the variables xn, fn, gn, pn, the corresponding norms
    #X and G have been computed.
    
    #Allocate next iteration:
    xnp1            = zeros(Nx)
    fnp1            = zeros(Nx)
    gnp1            = zeros(Nx)
    t_startup       += time()-t_ini
    
    #Iterate:
    while true

        #Display iter info:
        _NLS_Display()
        
        #Exit:
        if pnorm <= RelTolX*xnorm+AbsTolX
            ch.flag     = 1
            break
        elseif fnorm <= RelTolF*fnorm0+AbsTolF
            ch.flag     = 2
            break
        elseif gnorm <= RelTolG*gnorm0+AbsTolG
            ch.flag     = 3
#             println("t_residual= ", t_residual)
#             println("t_SVD= ", t_SVD)
#             println("t_Matrices= ", t_Matrices)
            break
        elseif ch.nIter==MaxIter
            ch.flag     = 0
            break
        end
        
        #Next iter:
        ch.nIter        += 1
        Iter2Restart    -= 1
        
        #Test fixed point iteration. Save results temporarily in _np1:
        t_ini           = time()
        @mlv xnp1       = xn - gn
        @mlv gnp1       = gn    #Initial condition for g_np1
        FunFlag         = fun!(xnp1, fnp1, gnp1)
        ch.fEvals       += 1
        if FunFlag<=0
            ch.flag     = -1
            break
        end
        t_residual      += time()-t_ini
        
        #Update matrices X and F:
        t_ini           = time()
        if Iter2Restart==0
        
            #Reset X ans F:
            @mlv Xm         = 0.0
            @mlv Fm         = 0.0
            @mlv Gm         = 0.0
            XGcols          = 0
            
            Iter2Restart    = Restart
            
        else 
        
            for jj=min(XGcols,memory-1):-1:1
                BLAS.copy!(view(Xm,:,jj+1), view(Xm,:,jj))
                BLAS.copy!(view(Fm,:,jj+1), view(Fm,:,jj))
                BLAS.copy!(view(Gm,:,jj+1), view(Gm,:,jj))
            end
            
        end
        @mlv Gm[:,1]    = gnp1-gn
        Deltag_norm_inv = 1.0/(NormFun(view(Gm,:,1))+1e-20)
        @mlv Gm[:,1]    *= Deltag_norm_inv
        @mlv Xm[:,1]    = (xnp1-xn)*Deltag_norm_inv
        @mlv Fm[:,1]    = (fnp1-fn)*Deltag_norm_inv
        XGcols          = min(memory, XGcols+1)
        t_Matrices      += time()-t_ini
        
        #Choose x_np1 such that f_np1 or g_np1 (last option seems best) are minimized:
        t_ini               = time()
        Gm_SVD              = svd(view(Gm,:,1:XGcols))
        Sinv                = Gm_SVD.S
        deg_bool            = Sinv.<sigma_min
        Sinv[deg_bool]      .= 0.0
        Sinv[.!deg_bool]    .= 1.0./Sinv[.!deg_bool]
        gammav              = Gm_SVD.V*(Sinv.*(transpose(Gm_SVD.U)*gn))
        BLAS.gemv!('N', -1.0, view(Xm,:,1:XGcols), gammav, 0.0, pn)
        BLAS.axpby!(1.0, pn, 1.0, xn)
        BLAS.gemv!('N', -1.0, view(Fm,:,1:XGcols), gammav, 1.0, fn)
        BLAS.gemv!('N', -1.0, view(Gm,:,1:XGcols), gammav, 1.0, gn)
        nSigma              = sum(.!deg_bool)
        t_SVD           += time()-t_ini
        
        #Update norms:
        pnorm           = NormFun(pn)
        fnorm           = NormFun(fn)
        gnorm           = NormFun(gn)
        _NLS_ch_update!()
        
    end
    
    _NLS_Display()
    
    return xn, ch

end
