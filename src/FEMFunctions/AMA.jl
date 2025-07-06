#-----------------------------------------------------------
#FORMULAS FOR OPTIMAL SIZES:

function AMA_OptimalSize(emeanv::Vector{Float64}, A0v::Vector{Float64}, 
    ratev::Union{N, Vector{N}},
    DPEv::Union{Int, Vector{Int}}, TolS::Float64; q::Real=2) where N<:Number
   
   return AMA_OptimalSizeEq(emeanv, A0v, ratev, DPEv, TolS; q=q)
#    return AMA_OptimalSize2(emeanv, A0v, ratev, DPEv, TolS; q=q)
   
end

#Optimal size with continous mesh model. Here "ratev" is the rate of emean w.r.t. h:
function AMA_OptimalSizeBraack(emeanv::Vector{Float64}, A0v::Vector{Float64}, 
    ratev::Union{N, Vector{N}},
    DPEv::Union{Int, Vector{Int}}, TolS::Float64; q::Real=2,
    Amin_A0::Float64=0.5, Amax_A0::Float64=Inf) where N<:Number

    if any((DPEv.-DPEv[1]).!=0) || any((ratev.-ratev[1]).!=0.0)
        error("This routine is only valid for h-refinement")
    end
    
    #------------------------------------------------------------
    #Optimal size:

    #Continuous mesh model:
    #
    #minimize
    #   N       = int_\Omega d / K d\Omega 
    #subject to
    #   int_\Omega e^q d\Omega < TolS^q |Omega| = Eta
    #with e(K)=e0(K0)*(K/K0)^(rate/2)
    #
    #We define eta= e^q K = B K^alpha K
    #for some constant B and alpha=q*rate/2.
    #
    #Using Lagrange multipliers, we arrive at
    #   B_K K^(alpha+1) = eta_K = const = eta_*
    #We compute eta_* by imposing that 
    #   Eta     = sum_K eta_* K0 / K = sum_K eta_* (eta_K/eta_*)^(1/(alpha+1)) 
    #           = (eta_*)^(alpha/(alpha+1)) * sum_K eta_K^(1/(alpha+1))
    #from which
    #   eta_*   = [ Eta/sum_K eta_K^(1/(alpha+1)) ]^[(alpha+1)/alpha]
    #Then
    #   K       = K0 * (eta_* / eta_K)^(1/(alpha+1))
    if q==Inf
        #For the Linf norm,
        #   eta_K alpha_K / d_K = e_K^q K alpha_K / d_K
        #                       = lambda 
        #                       = Eta/(sum_l d_l/alpha_l)
        #                       = TolS^q |\Omega| / (sum_l d_l/alpha_l)
        #So,
        #   e_K^q K alpha_K / d_K = TolS^q |\Omega| / (sum_l d_l/alpha_l)
        #and hence, in the limit q->Inf
        #   e_K = TolS
        #and
        #   K   = K0 * (e_K/e_K^0)^(2/rate)
        A1v         = @avxt @. exp( log(A0v) + 2.0/ratev*( log(TolS)-log(emeanv) ) )
        A1v         = @. max(Amin_A0*A0v, min(A1v, Amax_A0*A0v)) 
        return A1v
    end
    
    #Otherwise:
    logEta          = q*log(TolS) + log(sum(A0v))
    logeta0         = @avxt @. q*log(emeanv) + log(A0v)
    alphav          = @mlv q*ratev/2
#     logeta          = @avxt @. (alpha+1)/alpha*(logEta - log(sum(exp(logeta0/(alpha+1)))))
    logeta          = @avxt @. (alphav+1)/alphav*(logEta - prod(logeta0/(alphav+1)))
    
    #Formula based on convergence rate:
    A1v             = @avxt @. exp( log(A0v) + (logeta - logeta0)/(alphav+1) )
    
    return A1v

end

#Optimal size with continous mesh model and equidistribution principle. 
#Here "ratev" is the rate of emean w.r.t. h:
function AMA_OptimalSizeEq(emeanv::Vector{Float64}, A0v::Vector{Float64}, 
    ratev::Union{N, Vector{N}},
    DPEv::Union{Int, Vector{Int}}, TolS::Float64; q::Real=2) where N<:Number

    #------------------------------------------------------------
    #Optimal size:

    #Continuous mesh model:
    #
    #minimize
    #   N       = int_\Omega d / K d\Omega 
    #subject to
    #   int_\Omega e^q d\Omega < TolS^q |Omega| = Eta
    #with e(K)=e0(K0)*(K/K0)^(rate/2)
    #
    #We define eta= e^q K = B K^alpha K
    #for some constant B and alpha=q*rate/2.
    #
    #Using Lagrange multipliers, we arrive at
    #   eta_K alpha_K / d_K = lambda = const
    #We compute lambda by imposing that 
    #   Eta     = sum_K eta_K = sum_K lambda*d_K/alpha_K = lambda * sum_K d_K/alpha_K
    #so
    #   lambda  = Eta/(sum_K d_K/alpha_K)
    #   eta_K   = lambda * d_K / alpha_K [=Eta/NElems for h-adaptation]
    #and
    #   K       = K0*(eta_K/eta0_K)^(1/(alpha+1))
    
    if q==Inf
        #For the Linf norm,
        #   eta_K alpha_K / d_K = e_K^q K alpha_K / d_K
        #                       = lambda 
        #                       = Eta/(sum_l d_l/alpha_l)
        #                       = TolS^q |\Omega| / (sum_l d_l/alpha_l)
        #So,
        #   e_K^q K alpha_K / d_K = TolS^q |\Omega| / (sum_l d_l/alpha_l)
        #and hence, in the limit q->Inf
        #   e_K = TolS
        #and
        #   K   = K0 * (e_K/e_K^0)^(2/rate)
        A1v         = @mlv exp( log(A0v) + 2.0/ratev*( log(TolS)-log(emeanv) ) )
        return A1v
    end
    
    #Otherwise:
    logEta          = q*log(TolS) + log(sum(A0v))
    alphav          = @mlv q*ratev/2
    loglambda       = logEta - log(sum(DPEv./alphav))
    logeta          = @mlv loglambda + log(DPEv) - log(alphav)
    logeta0         = @mlv q*log(emeanv) + log(A0v)
    
    #Formula based on convergence rate:
    A1v             = @mlv exp( log(A0v) + (logeta - logeta0)/(alphav+1) )
    
    return A1v

end

#Formula based on thresholds. Here ratev is irrelevant for h-refinement:
function AMA_OptimalSize1(emeanv::Vector{Float64}, A0v::Vector{Float64}, 
    ratev::Vector{N}, DPEv::Vector{Int}, TolS::Float64; q::Real=2, rmax::Float64=0.25,
    cmax::Float64=2.5) where N<:Number

    #------------------------------------------------------------
    #Optimal size:

    #Continuous mesh model:
    #
    #minimize
    #   N       = int_\Omega d / K d\Omega 
    #subject to
    #   int_\Omega e^q d\Omega < TolS^q |Omega| = Eta
    #with e(K)=e0(K0)*(K/K0)^(rate/2)
    #
    #We define eta= e^q K = B K^alpha K
    #for some constant B and alpha=q*rate/2.
    #
    #Using Lagrange multipliers, we arrive at
    #   eta_K alpha_K / d_K = lambda = const
    #We compute lambda by imposing that 
    #   Eta     = sum_K eta_K = sum_K lambda*d_K/alpha_K = lambda * sum_K d_K/alpha_K
    #so
    #   lambda  = Eta/(sum_K d_K/alpha_K)
    #   eta_K   = lambda * d_K / alpha_K [=Eta/NElems for h-adaptation]
    #and
    #   K       = K0*(eta_K/eta0_K)^(1/(alpha+1))
    
    #Otherwise:
    logEta          = q*log(TolS) + log(sum(A0v))
    alphav          = @mlv q*ratev/2
    loglambda       = logEta - log(sum(DPEv./alphav))
    logeta          = @mlv loglambda + log(DPEv) - log(alphav)
    logeta0         = @mlv q*log(emeanv) + log(A0v)
    if q==Inf
        @mlv logeta     = log(TolS)
        @mlv logeta0    = log(emeanv)
    end
    
    #Thresholds:
    logeta0_max     = maximum(logeta0)
    logeta0_min     = minimum(logeta0)
    A1v             = @mlv 1.0*A0v
    #Refine:
    aux             = @mlv logeta0 > logeta
    xi              = @mlv (logeta0[aux]-logeta[aux])/(logeta0_max-logeta[aux]+1e-14)
    @mlv A1v[aux]   *= 0.5*(1.0+rmax) + 0.5*(1-rmax)*cos(pi*xi)
    #Coarse:
    aux             = @mlv logeta0 < logeta
    xi              = @mlv (logeta[aux]-logeta0[aux])/(logeta[aux]-logeta0_min+1e-14)
    @mlv A1v[aux]   *= 0.5*(1.0+cmax) - 0.5*(cmax-1)*cos(pi*xi)
    
#     figure()
#     plot(logeta0-logeta, A1v./A0v, ".b")
#     error("")
    
    return A1v

end

#Formula based on equidistribution and power law with truncation:
function AMA_OptimalSize2_old(emean0v::Vector{Float64}, A0v::Vector{Float64}, 
    ratev::Vector{N}, DPEv::Vector{Int}, TolS::Float64; q::Real=2, rmax::Float64=0.25,
    cmax::Float64=2.5, eta_ratio::Float64=10.0) where N<:Number

    #------------------------------------------------------------
    #Optimal size:

    #Continuous mesh model:
    #
    #minimize
    #   N       = int_\Omega d / K d\Omega 
    #subject to
    #   int_\Omega e^q d\Omega < TolS^q |Omega| = Eta
    #with e(K)=e0(K0)*(K/K0)^(rate/2)
    #
    #We define eta= e^q K = B K^alpha K
    #for some constant B and alpha=q*rate/2.
    #
    #Using Lagrange multipliers, we arrive at
    #   eta_K alpha_K / d_K = lambda = const
    #We compute lambda by imposing that 
    #   Eta     = sum_K eta_K = sum_K lambda*d_K/alpha_K = lambda * sum_K d_K/alpha_K
    #so
    #   lambda  = Eta/(sum_K d_K/alpha_K)
    #   eta_K   = lambda * d_K / alpha_K [=Eta/NElems for h-adaptation]
    #and
    #   K       = K0*(eta_K/eta0_K)^(1/(alpha+1))
    
    if q!=Inf
        error("")
    end
    
    #Target tolerances:
    logeta0         = @mlv log(emean0v)
    logTolS         = log(TolS)
    logDelta        = log(eta_ratio)
    logeta          = @mlv max(logeta0-logDelta, min(logeta0+logDelta, logTolS))
    A1v             = @mlv A0v * exp(2.0/ratev*(logeta-logeta0))
    
#     figure()
#     plot(logeta0-logeta, A1v./A0v, ".b")
#     error("")

    return A1v

end
function AMA_OptimalSize2(emean0v::Vector{Float64}, A0v::Vector{Float64}, 
    ratev::Union{N, Vector{N}}, DPEv::Union{Int, Vector{Int}}, TolS::Float64; 
    q::Real=2, eta_ratio::Float64=Inf, A_ratio::Float64=5.0) where N<:Number

    #------------------------------------------------------------
    #Optimal size:

    #Continuous mesh model:
    #
    #minimize
    #   N       = int_\Omega d / K d\Omega 
    #subject to
    #   int_\Omega e^q d\Omega < TolS^q |Omega|
    #with e(K)=e0(K0)*(K/K0)^(rate/2)
    #
    #Using Lagrange multipliers, we arrive at
    #   e_K^q K rate_K / d_K = const
    #or, equivalently,
    #   e_K K^(1/q) (rate_K/d_K)^(1/q)  = lambda
    #
    #We compute lambda by imposing that 
    #   TolS^q |Omega| = sum_K e_K^q K = sum_K lambda^q * d_K/alpha_K = lambda^q * sum_K d_K/alpha_K
    #so
    #   lambda          = TolS |Omega|^(1/q)/(sum_K d_K/alpha_K)^(1/q)
    #   e_K K^(1/q)     = lambda * (d_K/rate_K)^(1/q) 
    #                   = TolS |Omega|^(1/q) m_K^(1/q)
    #                   := lambda_star
    #with 
    #   m_K             = d_K/rate_K / sum(d_J/rate_J)
    #Note that m_k = (1/NElems) for h-adaptation.
    #Dividing by e_K0 * K0^(1/q):
    #   (K/K0)^(rate/2+1/q) = lambda_star/(eta0 K0^(1/q))
    #
    #Note that
    #   lambdastar          = lambda * (d_K/rate_K)^(1/q)
    #   lambdastar0         = lambda0 * (d_K/rate_K)^(1/q)
    #   lambda              = e_K K^(1/q) (rate_K/d_K)^(1/q)    known from equidistribution
    #   lambda0             = e0_K K0^(1/q) (rate_K/d_K)^(1/q)  known from initial data
    #   lambda/lambda0      = (K/K0)^(rate/2+1/q) 
    #   lambdastar/lambdastar0  = lambda/lambda0
    logmK           = @avxt @. log(DPEv/ratev)-log($sum(DPEv/ratev))  #log(1/N) in h-adaptation
    Omega           = sum(A0v)
    loglambdastar0  = @avxt @. log(emean0v) + 1.0/q*log(A0v)
    loglambdastar   = @avxt @. max(loglambdastar0-log(eta_ratio), 
                                    min(log(TolS) + 1.0/q*(log(Omega) + logmK), 
                                        loglambdastar0+log(eta_ratio)))
#     A1v             = @mlv exp( log(A0v) + (loglambdastar-loglambdastar0)/(ratev/2+1.0/q) )
    A1v             = @mlv max(A0v/A_ratio, 
                            min( exp( log(A0v) + (loglambdastar-loglambdastar0)/(ratev/2+1.0/q) ), 
                                A_ratio*A0v ))
    
    return A1v

end

#Optimal size with continous mesh model, equidistribution principle and a target number of elements,
#Here "ratev" is the rate of emean w.r.t. h:
function AMA_OptimalSize_NE(emean0v::Vector{Float64}, A0v::Vector{Float64}, 
    rate::Real,
    NElems::Real; q::Real=2)

    #------------------------------------------------------------
    #Optimal size:

    #Continuous mesh model:
    #
    #minimize
    #   int_\Omega e^q d\Omega
    #subject to
    #   NE  = int_\Omega d / K d\Omega 
    #with e ~ K^(rate/2) 
    #
    #Using Lagrange multipliers, we arrive at
    #   e_K K^(1/q) = lambda = const
    #If we divide by lambda0_K := e0_K K0^(1/q), 
    #   (K/K0)^alpha = lambda/lambda0_K
    #with alpha = (rate/2+1/q)
    #
    #We compute lambda by imposing that 
    #   NE      ~= sum_K K0/K = sum_K (lambda0_K/lambda)^(1/alpha) = 1/lambda^(1/alpha) * sum_K lambda0^(1/alpha)
    #so
    #   lambda^(1/alpha) = sum_K lambda0_K^(1/alpha) / NE
    #and
    #   K       = K0 * lambda^(1/alpha)/lambda0_K^(1/alpha)
    
    alpha           = rate/2.0 + 1.0/q
    #mu = lambda^(1/alpha)
    logmu0          = @avxt @. 1.0/alpha * ( log(emean0v) + 1.0/q*log(A0v) )
    logmu           = log( sum(exp.(logmu0))/NElems )
    
    #Formula based on convergence rate:
    A1v             = @mlv exp( log(A0v) + (logmu - logmu0) )
    
    return A1v

end

#A0v is the vector of current areas, and A1v is the vector of target areas.
#Returns estimated number of elements:
function AMA_CorrectSize!(A0v::Vector{Float64}, A1v::Vector{Float64}; 
    Amax::Float64=sum(A0v))

    #Impose maximum area:
    @mlv A1v        = min(A1v, Amax)
    
    #Estimated number of elements NE = int (dOmega/A) ~= sum_K A0_K/A_K:
    nElems          = length(A0v)
    NE_est          = sum(A0v/A1v)
    NE_min          = 0.1*nElems
    NE_max          = 10.0*nElems
    
#     return NE_est
    
    #Correct Av so that the number of elements is at most multiplied by 10:
    if NE_est>NE_max
        display("eeeh")
        A1v         *= NE_est/NE_max
        NE_est      = NE_max
    elseif NE_est<NE_min
        display("heee = $(NE_est/NE_min)")
        A1v         *= NE_est/NE_min
        NE_est      = NE_min
    end
    
    return NE_est
        
end

#-----------------------------------------------------------
#FORMULAS FOR BOUNDING POLYNOMIALS:

#Auxiliary function to compute powers 0,...,order of cos(theta) and sin(theta): 
function csthetam(thetav::Vector{Float64}, order::Int)

    nTheta          = length(thetav)
    cthetam         = zeros(nTheta, order+1)
    sthetam         = zeros(nTheta, order+1)
    cthetav         = @mlv cos(thetav)
    sthetav         = @mlv sin(thetav)
    cthetam[:,1]    .= 1.0
    sthetam[:,1]    .= 1.0
    @inbounds for jj=1:order
        @mlv cthetam[:,jj+1]    = $view(cthetam,:,jj)*cthetav
        @mlv sthetam[:,jj+1]    = $view(sthetam,:,jj)*sthetav
    end
    
    return cthetam, sthetam
    
end

#Given coefficients a_i through matrix A, compute 
#P = | sum(|i|=p) a_i x^i | (i is a multi-index)
#and approximate it by (x^T G x)^(p/2)
function AMA_BoundPolynomial(Pm::Matrix{Float64}, order::Int; nTheta::Int=181)

    nElems      = size(Pm,1)
    Pv          = zeros(nElems) #polynomial at arbitrary point
    P1v         = zeros(nElems) #maximum value of P
    P2v         = zeros(nElems) #minimum value of P
    phiv        = zeros(nElems) #angle of maximum value
    
    #Loop theta:
    thetav              = linspace(0.0, Float64(pi), nTheta)
    cthetam, sthetam    = csthetam(thetav, 2*order)
    aux_max             = zeros(Bool, nElems)
    delta               = eps(Float64)
    @inbounds for iTheta=1:nTheta
    
        #Compute Pv:
        @mlv Pv                     = 0.0
        for jj=0:order, ii=order-jj
            BLAS.axpby!(cthetam[iTheta,ii+1]*sthetam[iTheta,jj+1], view(Pm,:,jj+1), 1.0, Pv)
        end
        @mlv Pv                     = max(delta, abs(Pv))
        
        #Update Pmaxv:
        @mlv aux_max                = Pv > P1v
        @mlv phiv[aux_max]          = thetav[iTheta]
        @mlv P1v[aux_max]           = Pv[aux_max]
        
    end
    
    #Compute P in direction orthogonal to thetamax:
    theta2v             = @mlv phiv + pi/2
    c2m, s2m            = csthetam(theta2v, order)
    @inbounds for jj=0:order, ii=order-jj
        ew_mult!(1.0, view(Pm,:,jj+1), view(c2m,:,ii+1), view(s2m,:,jj+1), 1.0, P2v) 
    end
    @mlv P2v            = max(delta, P2v)

    #Compute spectral decomposition of G=Q*Lambda*Q^T:
    #   lambda_i^(order/2) = P_i
    #   sf_i        = lambda_i /lambda_bar = (P_i/P_bar)^(2/order)
    G                   = SpectralDecomp(nElems, 2)
    Pbarv               = @mlv sqrt(P1v*P2v)
    @mlv G.lambda_bar   = Pbarv^(2/order)
    @mlv G.sf[1]        = (P1v/Pbarv)^(2/order)
    @mlv G.sf[2]        = (P2v/Pbarv)^(2/order)
    @mlv G.V[1,1]       = cos(phiv)
    @mlv G.V[2,1]       = sin(phiv)
    @mlv G.V[1,2]       = -sin(phiv)
    @mlv G.V[2,2]       = cos(phiv)
    
    return G
    
end

#Compute etaK for given area, anisotropy and orientation. G
#contains the leading coefficients of the error, such that 
#   eta = int_E |e|^q dOmega = int_E (x^T G x)^(q*order/2) dOmega = emean^q A:
function AMA_Polynomial_errors(G::SpectralDecomp{Vector{Float64},Int}, 
    M::SpectralDecomp{Vector{Float64},Int}, order::Int; q::Number=2)

    #We have to compute I_m = int_E (x^T G x)^(m/2), for m=0,...,
    #E is the ellipse described by x^T Q D Q^T x <= 1.
    #We make Q^T x = mu * D^{-1/2} * t,
    #with 0<=mu<=1 and t=[cos(theta),sin(theta)]
    #(this is the classical change written in matrix form). Then
    #x  = mu * Q * D^{-1/2} * t.
    #The jacobian of the transformation is mu*|D|^{-1/2}.
    #So,
    #I_m    = 2 int_0^(pi) dtheta int_0^1 dmu mu mu^m (t^T B t)^(m/2) |D|^{-1/2}
    #       = 2 * |D|^{-1/2}/(m+2) * int_0^(pi) dtheta (t^T B t)^(m/2)
    #with B=D^{-1/2} Q^T G Q D^{-1/2}

    #Write G in matrix form (G_bar is left out to improve accuracy):
    Gm                  = GetMatrix(G, lambda_bar=1.0)  #this is V^T S V
    
    #D^{-1/2}:
    lambda_m12          = Vector{Vector{Float64}}(undef,2)
    for ii=1:2
        lambda_m12[ii]  = @mlv 1.0/sqrt(M.sf[ii])
    end
                            
    #Compute B:
    nElems              = size(Gm[1,1],1)
    Bm                  = Matrix{Vector{Float64}}(undef,2,2)
    alloc!(Bm, (nElems,))
    for ll=1:2, kk=1:2, jj=1:2, ii=1:2
        @mlv Bm[ii,jj]  += lambda_m12[ii]*M.V[kk,ii]*Gm[kk,ll]*M.V[ll,jj]*lambda_m12[jj]
    end
        
    #Evaluate (t^T B t)^(q*order/2) at quadrature nodes:
    n_qp                = Int(min(100, ceil(q*order+1)))
    theta_qp, w_qp      = GaussLegQuad(0.0, 1.0*pi, n_qp) #q*order=2*n_qp-1
    x_qp                = Vector{Vector{Float64}}(undef,2)
    x_qp[1]             = @mlv cos(theta_qp)
    x_qp[2]             = @mlv sin(theta_qp)
    tBtm                = zeros(nElems, n_qp)
    for ii=1:2, jj=1:2
        BLAS.gemm!('N', 'T', 1.0, Bm[ii,jj], x_qp[ii].*x_qp[jj], 
            1.0, tBtm)
    end
    #Correct underflow errors:
    @mlv tBtm       = max(tBtm, eps(Float64))
    
    #Get maximum for infinity norm:
    if q==Inf
        tBt_max     = maximum(tBtm, dims=2)[:]
        emeanv      = @mlv (tBt_max*G.lambda_bar/M.lambda_bar)^(order/2)
        etav        = @mlv emeanv*Inf
        return emeanv, etav
    end
    
    #Otherwise,
    #e^q:
    @mlv tBtm   = tBtm^(q*order/2)
    
    #Compute integral:
    etav        = zeros(nElems)
    BLAS.gemv!('N', 1.0, tBtm, w_qp, 0.0, etav)
    @mlv etav   *= 2.0/M.lambda_bar*(G.lambda_bar/M.lambda_bar)^(q*order/2)/(q*order+2)
    
    #Areas and mean error:
    Av                  = @mlv TrElem_Area/M.lambda_bar
    emeanv              = @mlv (etav/Av)^(1.0/q)
    
    return emeanv, etav

end

#Get metric tensor of current mesh:
function ComputeMetric(mesh::TrMesh)

    #Compute jacobian inverse at barycenters:
    Jinv    = JinvCompute(mesh, [0.0 0.0])
    
    #Compute metric tensor in matrix form: (M_ij = Jinv_ki * Jinv_kj)
    Mm      = Matrix{Vector{Float64}}(undef,2,2)
    alloc!(Mm, (mesh.nElems,))
    for kk=1:2, jj=1:2, ii=1:2
        @mlv Mm[ii,jj]  += $view(Jinv[kk,ii],:,1)*$view(Jinv[kk,jj],:,1)
    end
    
    #Write metric tensor in spectral form:
    Mm2     = hcat(Mm[1,1], Mm[2,1]+Mm[1,2], Mm[2,2])
    metric  = AMA_BoundPolynomial(Mm2, 2)
    
    return metric
    
end

#Correct metric:
function AMA_CorrectMetric!(A0v::Vector{Float64},
    metric::SpectralDecomp{Vector{Float64},Int};
    Amax::Float64=sum(A0v))
    
#     A1v                     = areas(metric)
#     return sum(A0v./A1v)
    
    #Get maximum and minimum radius. Recall that 1/h^2 = Mbar * sf:
    hmaxv           = _hmax(metric)
    hminv           = _hmin(metric)
    A1v             = areas(metric)
    L0              = sqrt(sum(A0v))
    
    #Impose maximum anisotropy for large elements
    Rmax                    = @mlv 2.0*hmaxv
    aux                     = @mlv hmaxv>0.9*L0
    @mlv Rmax[hmaxv>0.9*L0] = 0.9*L0
    if any(aux) && false
        println("A_max=", maximum(A1v[aux]), "A_min=", minimum(A1v[aux]))
        println("hmax_max=", maximum(hmaxv[aux]), "hmax_min=", minimum(hmaxv[aux]))
        println("hmin_max=", maximum(hminv[aux]), "hmin_min=", minimum(hminv[aux]))
    end
    
    #Impose bound for hmax. Correct hmin so area is the same:
#     Amax            = @mlv 100.0*A0v
#     Rmax            = @mlv sqrt(Amax/TrElem_Area)

    #Correct sizes:
    @mlv hmaxv      = min(Rmax, hmaxv)  
    @mlv hminv      = A1v/(TrElem_Area*hmaxv)
    #Impose bound for hmin:
    @mlv hminv      = min(Rmax, hminv)
    #Re-compute triangle areas = 3*sqrt(3)/4pi * |E|:
    A1v             = @mlv TrElem_Area*hmaxv*hminv
    
    #Estimated number of elements NE = int (dOmega/A) ~= sum_K A0_K/A_K:
    nElems          = length(A0v)
    NE_est          = sum(A0v./A1v)
    NE_min          = 0.1*nElems
    NE_max          = 10.0*nElems
    
#     display(norm(A0v./A1v.-1.0,Inf))
#     display(nElems)
#     display(NE_est)
    
    #Correct Av so that the number of elements is at most multiplied by 10:
    if NE_est>NE_max
#         display("eeeh")
        A1v         *= NE_est/NE_max
        NE_est      = NE_max
    elseif NE_est<NE_min
#         display("heee = $(NE_est/NE_min)")
        A1v         *= NE_est/NE_min
        NE_est      = NE_min
    end
    
    #Save new metric:
    @mlv metric.lambda_bar  = TrElem_Area/A1v
    @mlv metric.sf[1]       = hmaxv/hminv
    @mlv metric.sf[2]       = 1.0/metric.sf[1]
    
    return NE_est
        
end

#-----------------------------------------------------------
#FORMULAS FOR OPTIMAL ANISOTROPIES:

#Interpolation error and optimal anisotropy:
function AMA_OptimalAniso!(uv::Vector{<:GenVector{Float64}}, 
    fes::TrFES, nFacts::Vector{Float64}, ProjOrder::Int; 
    q::Number=2, DeltaAR::Float64=10.0)

#     t0          = time()
    
    mesh        = fes.mesh
    rate        = ProjOrder+1
    nVars       = length(uv)
    
    #DEBUG: Find closest element to given point:
    xbar        = barycenters(mesh)
    x0          = [0.7/sqrt(2), 0.7/sqrt(2)]
    aux         = @mlv $findmin( (xbar[:,1]-x0[1])^2.0 + (xbar[:,2]-x0[2])^2.0 )
    ElemSel     = Int(aux[2])

    #Matrix with high order derivatives at each element:
    HODm        = Solution_HODerivs_Compute(uv, fes, rate, [0.0 0.0])
    #Less memory size than a mass/stiffness matrix.
    
    #If ProjOrder>=fes.order, reconstruct derivatives:
    for pp=fes.order+1:ProjOrder+1
        HODerivs_Reconstruct!(HODm, fes.mesh, pp)
    end
        
    #Compute polynomial
    #   e^2 = 1/nVars * sum_I (e_I/n_I)^2,
    #where
    #   e_I = 1/nVars sum_(i+j=rate) 1/(i!j!) a_ij x_1^i x_2^j
    #represents the interpolation error in the variable I,
    #and a_ij are the high order derivatives.
    #The resulting polynomial is writen as 
    #   e^2 = sum(i+j=2*rate) b_ij x_1^i x_2^j
    e2m             = zeros(mesh.nElems, 2*rate+1)
    factv           = cumfact(2*rate)
    for II=1:nVars, j1=0:rate, i1=rate-j1, j2=0:rate, i2=rate-j2
        @mlv e2m[:,j1+j2+1]     += $view(HODm[II,i1+1,j1+1],:,1) * 
                                    $view(HODm[II,i2+1,j2+1],:,1) / 
                                    (factv[i1+1]*factv[j1+1]*factv[i2+1]*factv[j2+1]*
                                    nFacts[II]*nFacts[II])
    end
    @mlv e2m        /= nVars
    
    #Compute matrix G that bounds e^2, i.e., e^2 < (x^T G x)^(2*rate/2). 
    #Output is return in a SpectralDecomp type.
    G               = AMA_BoundPolynomial(e2m, 2*rate)
    
    #Errors for current anisotropy:
    emeanv0,                = AMA_Polynomial_errors(G, fes.mesh.metric, rate, q=q) 
     
    #Optimal metric and corresponding errors:
    metric                  = deepcopy(G)
    @mlv metric.lambda_bar  = fes.mesh.metric.lambda_bar
    emeanv,                 = AMA_Polynomial_errors(G, metric, rate, q=q)
    
#     figure()
#     plot(metric.sf[1], ".b")
    
    #Limit anisotropy. Aspect ratio=sqrt(sf[1]/sf[2])=sf[1].
    AR0                     = fes.mesh.metric.sf[1]
    @mlv metric.sf[1]       = min(DeltaAR*AR0, metric.sf[1])
    @mlv metric.sf[2]       = 1.0/metric.sf[1]
    
#     plot(metric.sf[1], ".g")
#     error("")
    
    #Total errors:
    Av                      = @mlv TrElem_Area/fes.mesh.metric.lambda_bar
    Omega                   = sum(Av)
    err0                    = (sum(@mlv emeanv0^q * Av)/Omega)^(1.0/q)
    err1                    = (sum(@mlv emeanv^q * Av)/Omega)^(1.0/q)
    if q==Inf
        err0                = maximum(emeanv0)
        err1                = maximum(emeanv)
    end
    
    #Plot error curves at ElemSel:
    if false
    
        nTheta              = 360
        thetav              = linspace(0.0, 2.0*Float64(pi), nTheta)
        cthetam, sthetam    = csthetam(thetav, 2*rate)
        fig                 = PyPlotFigure(w=11.50, h=4.50, top=0.3, bottom=0.3, left=0.3, right=0.8)
        
        #=
        #Plot P=e^2:
        P2v                 = zeros(nTheta)
        for II=1:nVars
            Pv              = zeros(nTheta)
            for jj=0:rate, ii=rate-jj
                @mlv Pv     += HODm[II,ii+1,jj+1][ElemSel,1]*
                                cthetam[:,ii+1]*sthetam[:,jj+1]/
                                factv[ii+1]/factv[jj+1]/nFacts[II]
            end
            @mlv P2v        += Pv*Pv
        end
        @. $plot(P2v*cos(thetav), P2v*sin(thetav), "xr")
        =#
        
        #DATA:
        # u   = @. 0.5*x[2]^6 + 2.0*x[2]^1*x[1]^5 - 1.5*x[2]^3*x[1]^3
        # hp0=0.1, FesOrder=3, AnisoOrder=2, [0,1]^2
        # x0          = [0.7/sqrt(2), 0.7/sqrt(2)]
        
        #Plot P=e^2:
        Pv                  = zeros(nTheta)
        order               = 2*rate
        for jj=0:order, ii=order-jj
            @mlv Pv         += e2m[ElemSel,jj+1]*cthetam[:,ii+1]*sthetam[:,jj+1]
        end
        @. $plot(Pv*cos(thetav), Pv*sin(thetav), "b")
        
        _,iPv_max           = findmax(Pv)
        println("theta1=", thetav[iPv_max]*180/pi)
        
        #=
        #Plot P=(x^T G x)^(2*rate/2):
        Pv                  = zeros(nTheta)
        xm                  = [ cos.(thetav), sin.(thetav) ]
        Gm                  = GetMatrix(G)
        for jj=1:2, ii=1:2
            @mlv Pv         += Gm[ii,jj][ElemSel]*xm[ii]*xm[jj]
        end
        @mlv Pv             = Pv^rate
        @. $plot(Pv*cos(thetav), Pv*sin(thetav), "g")
        =#
        
        #theta1:
        Gmax    = (G.lambda_bar[ElemSel]*G.sf[1][ElemSel])^(ProjOrder+1)
        Gmax_x  = Gmax*G.V[1,1][ElemSel]
        Gmax_y  = Gmax*G.V[2,1][ElemSel]
        plot([0.0*Gmax_x, 1.1*Gmax_x], [0.0, 1.1*Gmax_y], "--r")
        text(1.15*Gmax_x, 1.0*Gmax_y, latexstring("\\theta=\\theta_1"), color="r")
        plot(Gmax_x, Gmax_y, "xr")
        plot(-Gmax_x, -Gmax_y, "xr")
        
        #theta2:
        Gperp   = (G.lambda_bar[ElemSel]*G.sf[2][ElemSel])^(ProjOrder+1)
        Gperp_x = Gperp*G.V[1,2][ElemSel]
        Gperp_y = Gperp*G.V[2,2][ElemSel]
        plot([0.0*Gperp_x, 1.6*Gperp_x], [0.0, 1.6*Gperp_y], "--k")
        text(5.0*Gperp_x, 1.7*Gperp_y, latexstring("\\theta=\\theta_2"), color="k")
        plot(Gperp_x, Gperp_y, "xk")
        plot(-Gperp_x, -Gperp_y, "xk")
        
#         legend([latexstring("G_K"), latexstring("\\mathcal{G}_K")])
        legend([latexstring("\\theta \\mapsto G_K(\\theta) \\mathbf{e}(\\theta)")])
        tick_params(axis="both", which="both", labelsize=8)
        axis("off")
        axis("equal")
        
#         savefig("/home/manu/Documents/Proyectos/shared/AMA/66acb481d4fb9dd2497f2599/v2/aniso.eps",
#             pad_inches=0, dpi=800)
        savefig("/home/manu/Documents/Proyectos/shared/AMA/6707f6b998313315eabaa3a7/aniso.eps",
            pad_inches=0, dpi=800)
        
        error("")
        
    end
    
    return emeanv0, emeanv, err0, err1, metric
    
end

#Proposed formula, based on Chakraborty et al. (CMA2022), but without the optimization process:
function AMA_OptimalAniso_CD!(erec::Vector{<:GenVector{Float64}}, 
    fes::TrFES, nFacts::Vector{Float64}; DeltaAR::Float64=10.0)

    mesh        = fes.mesh
    nVars       = length(erec)
    order       = fes.order
    
    #DEBUG: Find closest element to given point:
    xbar        = barycenters(mesh)
    x0          = [1.0, 0.5]
    aux         = @mlv $findmin( (xbar[:,1]-x0[1])^2.0 + (xbar[:,2]-x0[2])^2.0 )
    ElemSel     = Int(aux[2])
    
    #Matrix with high order derivatives at each element:
    HODm        = Solution_HODerivs_Compute(erec, fes, order, [0.0 0.0])
    #Less memory size than a mass/stiffness matrix.
    
    #Compute polynomial
    #   e^2 = 1/nVars * sum_I (e_I/n_I)^2,
    #where
    #   e_I = 1/nVars sum_(0<=i+j<=order) 1/(i!j!) a_ij x_1^i x_2^j
    #represents the interpolation error in the variable I,
    #and a_ij are the high order derivatives.
    #The resulting polynomial is writen as 
    #   e^2 = sum(0<=i+j<=2*order) b_ij x_1^i x_2^j
    e2m             = Vector{Matrix{Float64}}(undef, 2*order+1) #e2m[ii+1] saves coefficients of term of degree ii:
    factv           = cumfact(2*order)
    for pp=0:2*order
        e2m[pp+1]   = zeros(mesh.nElems, pp+1)
        for II=1:nVars, p1=max(0,pp-order):min(pp,order), p2=pp-p1, j1=0:p1, i1=p1-j1, j2=0:p2, i2=p2-j2
            @avxt @. e2m[pp+1][:, j1+j2+1]  += $view(HODm[II,i1+1,j1+1],:,1) * 
                                                $view(HODm[II,i2+1,j2+1],:,1) / 
                                                (factv[i1+1]*factv[j1+1]*factv[i2+1]*factv[j2+1]*
                                                nFacts[II]*nFacts[II])
        end
        @avxt @. e2m[pp+1]      /= nVars
    end
    
    #Loop even polynomial degrees. Bound polynomials with matrix G. The optimal metric is 
    #proportional to G:
    G               = Vector{ SpectralDecomp{Vector{Float64},Int} }(undef, 2*order+1)
    wv              = zeros(2*order+1)      #weights for optimization
    for pp=2:2:2*order
        #Bound polynomial:
        G[pp+1]                 = AMA_BoundPolynomial(e2m[pp+1], pp)
        #Compute error for this term:
        metric                  = deepcopy(G[pp+1])
        @mlv metric.lambda_bar  = fes.mesh.metric.lambda_bar
        emeanv,                 = AMA_Polynomial_errors(G[pp+1], metric, pp, q=2)
        Av                      = TrElem_Area./fes.mesh.metric.lambda_bar
        Omega                   = sum(Av)
        wv[pp]                  = (sum(@. emeanv^2 * Av)/Omega)^(1/2)
    end
    
    #Plot error curves at ElemSel:
    if true
    
        nTheta              = 360
        thetav              = linspace(0.0, 2.0*Float64(pi), nTheta)
        cthetam, sthetam    = csthetam(thetav, 2*order)
#         fig                 = PyPlotFigure(w=11.50, h=4.50, top=0.3, bottom=0.3, left=0.3, right=0.8)
        colorv              = PyPlotColors("jet2", 2*order+1)
        
        #Loop terms of degree p in Taylor expansion:
        for pp=2:2:2*order
        
            figure()
            
            #Compute terms from high order derivatives:
            Pv                  = zeros(nTheta)
            for II=1:nVars, p1=max(0,pp-order):min(pp,order), p2=pp-p1
                P1v             = zeros(nTheta)
                for j1=0:p1, i1=p1-j1
                    @mlv P1v    += HODm[II,i1+1,j1+1][ElemSel,1]*
                                    cthetam[:,i1+1]*sthetam[:,j1+1]/
                                    factv[i1+1]/factv[j1+1]/nFacts[II]
                end
                P2v             = zeros(nTheta)
                for j2=0:p2, i2=p2-j2
                    @mlv P2v    += HODm[II,i2+1,j2+1][ElemSel,1]*
                                    cthetam[:,i2+1]*sthetam[:,j2+1]/
                                    factv[i2+1]/factv[j2+1]/nFacts[II]
                end
                @mlv Pv         += P1v*P2v
            end
            @. $plot(Pv*cos(thetav), Pv*sin(thetav), color=colorv[pp+1], marker="x", linestyle="none")
        
            #Plot terms from matrix e2m:
            Pv              = zeros(nTheta)
            for jj=0:pp, ii=pp-jj
                @mlv Pv     += e2m[pp+1][ElemSel,jj+1]*cthetam[:,ii+1]*sthetam[:,jj+1]
            end
            @. $plot(Pv*cos(thetav), Pv*sin(thetav), color=colorv[pp+1], marker="+", linestyle="none")
            
            #Plot P=(x^T G x)^(2*rate/2):
            Pv                  = zeros(nTheta)
            xm                  = [ cos.(thetav), sin.(thetav) ]
            Gm                  = GetMatrix(G[pp+1])
            for jj=1:2, ii=1:2
                @mlv Pv         += Gm[ii,jj][ElemSel]*xm[ii]*xm[jj]
            end
            @mlv Pv             = Pv^(pp/2)
            @. $plot(Pv*cos(thetav), Pv*sin(thetav), color="k")
        
            axis("equal")
            title("p=$(pp)")
            
        end
        
#         _,iPv_max           = findmax(Pv)
#         println("theta1=", thetav[iPv_max]*180/pi)
        
        
        
        error("")
        
        #theta1:
        Gmax    = (G.lambda_bar[ElemSel]*G.sf[1][ElemSel])^(ProjOrder+1)
        Gmax_x  = Gmax*G.V[1,1][ElemSel]
        Gmax_y  = Gmax*G.V[2,1][ElemSel]
        plot([0.0*Gmax_x, 1.1*Gmax_x], [0.0, 1.1*Gmax_y], "--r")
        text(1.15*Gmax_x, 1.0*Gmax_y, latexstring("\\theta=\\theta_1"), color="r")
        plot(Gmax_x, Gmax_y, "xr")
        plot(-Gmax_x, -Gmax_y, "xr")
        
        #theta2:
        Gperp   = (G.lambda_bar[ElemSel]*G.sf[2][ElemSel])^(ProjOrder+1)
        Gperp_x = Gperp*G.V[1,2][ElemSel]
        Gperp_y = Gperp*G.V[2,2][ElemSel]
        plot([0.0*Gperp_x, 1.6*Gperp_x], [0.0, 1.6*Gperp_y], "--k")
        text(5.0*Gperp_x, 1.7*Gperp_y, latexstring("\\theta=\\theta_2"), color="k")
        plot(Gperp_x, Gperp_y, "xk")
        plot(-Gperp_x, -Gperp_y, "xk")
        
#         legend([latexstring("G_K"), latexstring("\\mathcal{G}_K")])
        legend([latexstring("\\theta \\mapsto G_K(\\theta) \\mathbf{e}(\\theta)")])
        tick_params(axis="both", which="both", labelsize=8)
        axis("off")
        axis("equal")
        
        error("")
        
    end
    
end

#-----------------------------------------------------------
#RECONSTRUCTION OPERATORS:

#Dolejsi's patch reconstruction. Valid for linear meshes:
function AMA_Dolejsi_PatchReconstruction(u0::Vector{<:GenVector{Float64}},
    fes0::TrPSpace, fes1::TrQSpace)

    t_ini       = time()
    
    mesh        = fes0.mesh
    if mesh !== fes1.mesh
        error("Mesh must be the same for both spaces")
    end
    if mesh.order>1
        @warn "Patch reconstruction is valid only for straight meshes"
    end
    
    #DEBUG: Find closest element to given point:
    xbar        = barycenters(mesh)
    aux         = findmin( (view(xbar,:,1).-0.5).^2.0 + 
                           (view(xbar,:,2).-0.5).^2.0   )
    ElemSel     = Int(aux[2])

    #Chebyshev nodes at the reference triangle:
    xi_Lag      = TrElem_ChebyshevNodes(fes1.order)
    nLag        = size(xi_Lag,1)
    x_Lag       = LagNodes(mesh, xi_Lag)
    #Lagrange polynomials in terms of PKD basis:
    Lag_PKD     = Lag_PKDCoefficients(fes1.order, xi_Lag)
    
    #Express solution in PKD basis (coefficients of degree >fes0.order are null):
    nVars       = length(u0)
    u0_PKD      = Vector{Matrix{Float64}}(undef,nVars)
    for II=1:nVars
        u0_PKD[II]  = zeros(mesh.nElems, fes0.DofPerElem)
        BLAS.gemm!('N', 'T', 1.0, u0[II][fes0.ElemsDof], fes0.Lag_PKD, 
            0.0, u0_PKD[II])
    end
    
    #Allocate mass matrices (CSC format):
    mass11_iv       = zeros(Int, fes1.DofPerElem, fes1.DofPerElem, mesh.nElems)
    mass11_jv       = zeros(Int, mesh.nElems*fes1.DofPerElem+1)
    mass11_sv       = zeros(fes1.DofPerElem, fes1.DofPerElem, mesh.nElems)
    bv_rec          = zeros(nVars*fes1.nDof)
    b_rec           = GetViews(bv_rec, nVars, fes1.nDof)
#     println("Reconstruction, t_alloc = ", time()-t_ini)
    
    #Loop elements:
#     @warn "inbounds"
    Av              = areas(mesh)
    Jinv            = JinvCompute(mesh, [0.0 0.0])
    xi_Lag_iElem    = zeros(nLag, 2)    #Coordinates of Lagrangian nodes w.r.t. central element
    Rm              = zeros(fes1.DofPerElem, fes1.DofPerElem)   #Restriction matrix
    @inbounds for iElem=1:mesh.nElems
    
        #Get patch:
        PatchOrder      = 1
        if any(view(mesh.ElemsNeigh,iElem,:).==0)
            PatchOrder  = 2
        end
        nElemsPatch, ElemsPatch,    = ElemsPatch_Compute(mesh, iElem, PatchOrder)
        
        #Loop patch. At each element in the patch, we express the PKD functions on iElem
        #as a linear combination of the PKD basis functions on the patch element.
        for jElem in ElemsPatch
            
            #Natural coordinates of jElem's Lagrange nodes w.r.t. iElem. Use Newton--Raphson 
            #formula with initial guess xi0_i = 0
            #   xi_i    = xi0_i - Jinv_ij*(x_j(xi0)-x_j)
            #If mesh is not straight, we must iterate until convergence.
            @mlv xi_Lag_iElem           = 0.0
            for ii=1:2, jj=1:2
#                 @mlv xi_Lag_iElem[:,ii] -= Jinv[ii,jj][iElem,1]*
#                                             (xbar[iElem,jj]-x_Lag[jj][jElem,:])
                @mlv $view(xi_Lag_iElem,:,ii)   -= Jinv[ii,jj][iElem,1]*
                                                    (xbar[iElem,jj]-$view(x_Lag[jj],jElem,:))
            end
            
            #Evaluate PKD basis at Lagrange nodes:
            Nm          = NCompute(fes1, xi_Lag_iElem)
            
            #At jElem, let phi_i be a basis function for the patch reconstruction, 
            #L_j a Lagrangian basis function (associated with node x^j),
            #and psi_k a basis (PKD) function. Then,
            #   phi_i(x)= phi_i(x^j) * L_j(x) = phi_i(x^j) * psi_k(x) * Lag_PKD_(k,j)
            #           = R_ik * psi_k(x)
            @mlv Rm         = 0.0
            for kk=1:fes1.DofPerElem, ii=1:fes1.DofPerElem, jj=1:nLag
                Rm[ii,kk]   += Nm[jj,ii] * Lag_PKD[kk,jj]
            end
            
            #If phi_i = R_ij psi_j, then
            #   (phi_i, phi_j)_jElem = R_ik R_jl (psi_k, psi_l)_jElem = R_ik R_jk |K_jElem|: 
            for jj=1:fes1.DofPerElem, ii=1:fes1.DofPerElem, kk=1:fes1.DofPerElem
                mass11_sv[ii,jj,iElem]      += Rm[ii,kk]*Rm[jj,kk]*Av[jElem]
            end
            
            #Similarly,
            #   (phi_i, u0)_jElem   = R_ik (psi_k, u0)_jElem = R_ik b_k,
            #with b_k the PKD coefficient. Note that b_k=0 when k is a shape function
            #of degree > fes0.order
            for II=1:nVars, ii=1:fes1.DofPerElem, kk=1:fes0.DofPerElem
                index_ii                = (iElem-1)*fes1.DofPerElem + ii
                b_rec[II][index_ii]     += Rm[ii,kk]*u0_PKD[II][jElem,kk]*Av[iElem]
            end
            
        end

    end
    
    #Indices of sparse matrix:
    @inbounds for iElem=1:mesh.nElems, jj=1:fes1.DofPerElem, ii=1:fes1.DofPerElem
        mass11_iv[ii,jj,iElem]  = (iElem-1)*fes1.DofPerElem + ii
    end
    mass11_jv[1]            = 1
    @inbounds for jj=1:fes1.nDof
        mass11_jv[jj+1]     = mass11_jv[jj] + fes1.DofPerElem
    end
    #
#     println("Reconstruction, t_loop = ", time()-t_ini)
    
    #Solve systems:
    M11_rec         = SparseMatrixCSC{Float64,Int}(fes1.nDof, fes1.nDof, 
                        mass11_jv, reshape(mass11_iv,:), reshape(mass11_sv,:))
    M11_rec_inv     = BlockDiagInv_GaussJordan(M11_rec, fes1.DofPerElem)
    urecv           = zeros(nVars*fes1.nDof)
    urec            = GetViews(urecv, nVars, fes1.nDof)
    for II=1:nVars
        urec[II]    .= M11_rec_inv*b_rec[II]
    end
    
    
#     println("Reconstruction, t_LS = ", time()-t_ini)
    
#     figure()
#     plot3D(urec[1], fes1)
    
    return urec

end

#Bermejo--Carpio (2009) reconstruction:
function AMA_BC_PatchReconstruction(u0::Vector{<:GenVector{Float64}},
    fes0::TrFES, fes1::TrQSpace)
    
    flag        = 1
    
    t_ini       = time()
    
    mesh        = fes0.mesh
    if mesh !== fes1.mesh
        error("Mesh must be the same for both spaces")
    end
    
    #DEBUG: Find closest element to given point:
    xbar        = barycenters(mesh)
    aux         = findmin( (view(xbar,:,1).-0.5).^2.0 + 
                           (view(xbar,:,2).-0.5).^2.0   )
    ElemSel     = Int(aux[2])

    #Chebyshev nodes of the expanded triangles Khat:
    xi_Lag      = TrElem_ChebyshevNodes(fes1.order).*(fes1.order/fes0.order)
    nLag        = size(xi_Lag,1)
    x_Lag       = LagNodes(mesh, xi_Lag)

    if false
        #First order patch:
        nElemsPatch, ElemsPatch,    = ElemsPatch_Compute(mesh, ElemSel, 4)
        
        fig     = PyPlotFigure(w=4.50, h=4.50, top=0.3, bottom=0.3, left=0.3, right=0.3)
        PlotMesh!(mesh, elems=ElemsPatch, color="gray")
        PlotMesh!(mesh, elems=[ElemSel], color="b")
        plot(x_Lag[1][ElemSel,:], x_Lag[2][ElemSel,:], "or", markersize=1.0)
        tick_params(axis="both", which="both", labelsize=TickSize)
        axis("off")
        xlim(0.2,1.0)
        ylim(0.1,0.9)
        savefig("/home/manu/Documents/Proyectos/shared/AMA/66acb481d4fb9dd2497f2599/v2/reconstruction.eps",
            pad_inches=0, dpi=800)
        error("")
    end
    
    #Search--locate nodes:
    ElemGuess   = reshape(Vector(1:mesh.nElems)*ones(Int,1,nLag),:)
    SL_output   = SearchLocate(mesh, [reshape(x_Lag[1],:), reshape(x_Lag[2],:)],
                    ElemGuess=ElemGuess)
    flag        = SL_output[1]
    elem0v      = SL_output[2]
    xi0m        = SL_output[3]
    if flag<0
        @warn "Search--locate algorithm failed"
    end
    
    #=
    #Find nodes outside the domain and extrapolate solution from K:
    aux         = findall(elem0v.<0)
    elem0v[aux] .= ElemGuess[aux]
    fv          = zeros(2)
    Jinvm       = zeros(2,2)
    Nm          = NCompute(mesh, [0.0 0.0])
    GradNm      = gradNCompute(mesh, [0.0 0.0])
    for iS in aux
        ElemSearched    = elem0v[iS]
        #Compute residual:
        fv[1]           = -x_Lag[1][iS]
        fv[2]           = -x_Lag[2][iS]
        for ii=1:mesh.NodesPerElem
            fv[1]       += Nm[1,ii]*mesh.ElemsCoords[ElemSearched,2*ii-1]
            fv[2]       += Nm[1,ii]*mesh.ElemsCoords[ElemSearched,2*ii]
        end
        #Compute jacobian and determinant of isoparametric transformation:
        Jinvm           .= 0.0
        for ii=1:mesh.NodesPerElem
            Jinvm[1,1]  += mesh.ElemsCoords[ElemSearched,2*ii]*GradNm[2][1,ii]
            Jinvm[2,1]  -= mesh.ElemsCoords[ElemSearched,2*ii]*GradNm[1][1,ii]
            Jinvm[1,2]  -= mesh.ElemsCoords[ElemSearched,2*ii-1]*GradNm[2][1,ii]
            Jinvm[2,2]  += mesh.ElemsCoords[ElemSearched,2*ii-1]*GradNm[1][1,ii]
        end
        detJ            = Jinvm[1,1]*Jinvm[2,2]-Jinvm[2,1]*Jinvm[1,2]
        Jinvm[1,1]      /= detJ
        Jinvm[2,1]      /= detJ
        Jinvm[1,2]      /= detJ
        Jinvm[2,2]      /= detJ
        
        #Compute Newton's update:
        for ii=1:2
            xi0m[iS,ii]         = 0.0
            for jj=1:2
                xi0m[iS,ii]     -= Jinvm[ii,jj]*fv[jj]
            end
        end
        
    end
    =#
    
    #Evaluate donor solution at nodes. For the nodes outside the domain, we extrapolate
    #the solution from the closest element:
    u_Lag       = SolutionCompute(u0, fes0, abs.(elem0v), xi0m)
#     u_Lag       = SolutionCompute(u0, fes0, elem0v, xi0m)
    
    #The solution "u" can be written in a PKD basis as 
    #   u(x) = sum_j a_j P_j(x)
    #The coefficients a_j are obtained from u from 
    #   u(x_i)  = sum_j a_j P_j(x_i)
    #   a_i     = P^{-1}_ij u_j
    Lag_PKD     = inv(NCompute(fes1, xi_Lag))
    
#     QRule       = WVRule(20)
#     QRule.xi    *= fes1.order/fes0.order
#     Nm_qp       = NCompute(fes1, QRule.xi)
#     x_qp        = LagNodes(mesh, QRule.xi)
    
    #Change basis:
    nVars       = length(u0)
    urecv       = zeros(nVars*fes1.nDof)
    urec        = GetViews(urecv, nVars, fes1.nDof)
    a_II        = zeros(mesh.nElems, fes1.DofPerElem) #PKD coefficients
    for II=1:nVars
        u_II    = reshape(u_Lag[II], mesh.nElems, nLag)
        #a_i (nElems, nDof) = u_j (nElems, nDof) * Lag_PKD^T
        BLAS.gemm!('N', 'T', 1.0, u_II, Lag_PKD, 0.0, a_II)
        #Assemble a_II onto u_rec:
        VectorAssemble!(fes1.ElemsDof, a_II, urec[II])
    end

    return urec, flag
    
end
function AMA_BC_PatchReconstruction(fes0::TrFES, fes1::TrQSpace)

    @warn "Deprecated"
    
    flag        = 1
    
    t_ini       = time()
    
    mesh        = fes0.mesh
    if mesh !== fes1.mesh
        error("Mesh must be the same for both spaces")
    end
    
    #DEBUG: Find closest element to given point:
    xbar        = barycenters(mesh)
    aux         = findmin( (view(xbar,:,1).-3.5).^2.0 + 
                           (view(xbar,:,2).-0.0).^2.0   )
    ElemSel     = Int(aux[2])

    #Chebyshev nodes of the expanded triangles Khat:
    xi_Lag      = TrElem_ChebyshevNodes(fes1.order).*(fes1.order/fes0.order)
    nLag        = size(xi_Lag,1)
    x_Lag       = LagNodes(mesh, xi_Lag)
    
    #Search--locate nodes:
    ElemGuess           = min.(mesh.nElems, 
                            reshape(Vector(1:mesh.nElems)*ones(Int,1,nLag),:))
    flag, elem0v, xi0m  = SearchLocate(fes0.mesh, [reshape(x_Lag[1],:), reshape(x_Lag[2],:)])
    elem0v              = @mlv abs(elem0v)
    if flag<0
        @warn "Search--locate algorithm failed"
    end
    
    #Evaluate shape functions of fes0 at interpolation nodes:
    Nm          = NCompute(fes0, xi0m)
    
    #Create sparse matrix that maps original solution to reconstructed solution at the 
    #interpolation nodes:
    iv          = zeros(Int, mesh.nElems*nLag, fes0.DofPerElem)
    jv          = zeros(Int, mesh.nElems*nLag, fes0.DofPerElem)
    sv          = zeros(mesh.nElems*nLag, fes0.DofPerElem)
    iv_aux      = reshape( fes1.ElemsDof, : )
    for iDof=1:fes0.DofPerElem
        @views @mlv iv[:,iDof]  = iv_aux
        @views @mlv jv[:,iDof]  = fes0.ElemsDof[elem0v, iDof]
        @views @mlv sv[:,iDof]  = Nm[:,iDof]
    end
    A1          = sparse(reshape(iv,:), reshape(jv,:), reshape(sv,:))
        
    #The solution "u" can be written in a PKD basis as 
    #   u(x) = sum_j a_j P_j(x)
    #The coefficients a_j are obtained from u from 
    #   u(x_i)  = sum_j a_j P_j(x_i)
    #   a_i     = P^{-1}_ij u_j
    Lag_PKD     = inv(NCompute(fes1, xi_Lag))
    A2          = SpBlockDiag(Lag_PKD, mesh.nElems)
    
    return A2*A1, flag
    
end

#Conservative reconstruction:
function AMA_CD_PatchReconstruction(u0::Vector{<:GenVector{Float64}}, 
    fes0::TrPSpace, fes1::TrQSpace)
    
    flag        = 1
    
    t_ini       = time()
    
    nVars       = length(u0)
    mesh        = fes0.mesh
    if mesh !== fes1.mesh
        error("Mesh must be the same for both spaces")
    end
    if fes0.order+3!=fes1.order
        error("Reconstruction order must be solution order + 3")
    end
    
    #DEBUG: Find closest element to given point:
    xbar        = barycenters(mesh)
    aux         = findmin( (view(xbar,:,1).-3.5).^2.0 + 
                           (view(xbar,:,2).-0.0).^2.0   )
    ElemSel     = Int(aux[2])

    #Chebyshev nodes of the expanded triangles Khat. We retain only the 3*(p+3) nodes 
    #at \partial Khat:
    xi_Lag      = TrElem_ChebyshevNodes(fes1.order)[1:3*fes1.order,:].*
                    (fes1.order/fes0.order)
    nLag        = size(xi_Lag,1)
    x_Lag       = LagNodes(mesh, xi_Lag)
    
#     figure()
#     PlotMesh!(mesh, elems=[ElemSel])
#     plot(x_Lag[1][ElemSel,:], x_Lag[2][ElemSel,:], ".k")
#     error("")
    
    #Search--locate nodes:
    ElemGuess   = reshape(Vector(1:mesh.nElems)*ones(Int,1,nLag),:)
    SL_output   = SearchLocate(mesh, [reshape(x_Lag[1],:), reshape(x_Lag[2],:)],
                    ElemGuess=ElemGuess)
    flag        = SL_output[1]
    elem0v      = SL_output[2]
    xi0m        = SL_output[3]
    if flag<0
        @warn "Search--locate algorithm failed"
    end
    
    #Interpolate solution of fes0 at interpolation nodes (including points that 
    #lie outside the domain):
    u_SL        = SolutionCompute(u0, fes0, abs.(elem0v), xi0m)
    u_Lag       = Vector{Matrix{Float64}}(undef,nVars)
    for II=1:nVars
        u_Lag[II]   = reshape(u_SL[II], mesh.nElems, nLag)
    end
    
    #Evaluate contribution of lower-order terms:
    u_l         = SolutionCompute(u0, fes0, xi_Lag)
    
    #The relation between Lagrangian and PKD coefficients is
    #   u = L_k u_k = P_j A_jk u_k = P_j a_j,
    #where a_j  = A_jk u_k
    #and A_jk = fes0.Lag_PKD
    #So, for the first (p+1)(p+2)/2 d.o.f., we impose
    #   a_j          = A_jk u_k,     j=1,...,(p+1)(p+2)/2
    #whereas for the interpolation nodes we impose
    #   P_j(x_i) a_j = u_Lag(x_i)-u_l(x_i),       i=1,...,number of interp. nodes inside the domain
    #   a_j          = Pinv_jk (u_Lag(x_k)-u_l(x_l))
    
    #Compute reconstruction:
    Pinv        = inv(NCompute(fes1, xi_Lag)[:,fes0.DofPerElem+1:fes1.DofPerElem])
    Lag_PKD     = fes0.Lag_PKD
    urecv       = zeros(nVars*fes1.nDof)
    urec        = GetViews(urecv, nVars, fes1.nDof)
    am          = zeros(mesh.nElems, fes1.DofPerElem)
    for II=1:nVars
        
        #Low order coefficients:
        u_II                        = u0[II][fes0.ElemsDof]
        am[:,1:fes0.DofPerElem]     = u_II*transpose(Lag_PKD)
        
        #Right hand side vector, as a matrix (nElems, nLag):
        rm                          = u_Lag[II]-u_l[II]
        
        #High-order coefficients:
        am[:,fes0.DofPerElem+1:fes1.DofPerElem]     = rm*transpose(Pinv)
        
        #Save coefficients into urec:
        urec[II]                    .= transpose(am)[:]
        
    end
    
#     figure()
#     plot3D(u0[1], fes0)
#     figure()
#     plot3D(urec[1], fes1)
#     error("")
    
    return urec, flag

end

#-----------------------------------------------------------
#OTHERS:

#Compute intersection of several metrics:
function AMA_MetricIntersection(metrics::Vector{SpectralDecomp{Vector{Float64},Int}};
    nTheta::Int=181)

    nMeshes     = length(metrics)
    nElems      = size(metrics[1].lambda_bar,1)
    Pv          = zeros(nElems) #polynomial at arbitrary point
    P1v         = zeros(nElems) #maximum value of P
    P2v         = zeros(nElems) #minimum value of P
    phiv        = zeros(nElems) #angle of maximum value
    
    #Metric tensor as matrices:
    Mm          = Vector{Matrix{Vector{Float64}}}(undef, nMeshes)
    for iMesh=1:nMeshes
        Mm[iMesh]       = GetMatrix(metrics[iMesh]) #Mm[iMesh][ii,jj] is the vector with M_ij at each element
    end
    
    #Loop meshes and theta:
    thetav              = linspace(0.0, Float64(pi), nTheta)
    xm                  = [ cos.(thetav), sin.(thetav) ]
    aux_max             = zeros(Bool, nElems)
    delta               = eps(Float64)
    @inbounds for iTheta=1:nTheta
    
        #Compute Pv:
        for iMesh=1:nMeshes
        
            #P=x^T M x:
            @avxt @. Pv             = 0.0
            for jj=1:2, ii=1:2
                @avxt @. Pv         += xm[ii][iTheta]*Mm[iMesh][ii,jj]*xm[jj][iTheta]
            end
            @avxt @. Pv             = max(delta, Pv)
        
            #Update Pmaxv:
            @avxt @. aux_max        = Pv > P1v
            @avxt @. phiv[aux_max]  = thetav[iTheta]
            @avxt @. P1v[aux_max]   = Pv[aux_max]
            
        end
        
    end
    
    #Compute P in direction orthogonal to thetamax:
    theta2v             = @avxt @. phiv + pi/2
    xm                  = [ cos.(theta2v), sin.(theta2v) ]
    @inbounds for iMesh=1:nMeshes
        
        #P=x^T M x:
        @avxt @. Pv             = 0.0
        for jj=1:2, ii=1:2
            @avxt @. Pv         += xm[ii]*Mm[iMesh][ii,jj]*xm[jj] 
        end
        @avxt @. Pv             = max(delta, Pv)
        
        #Update P2:
        @avxt @. P2v            = max(P2v, Pv)
        
    end
    
    #Compute spectral decomposition of metric tensor:
    #   lambda_i^(order/2) = P_i
    #   sf_i        = lambda_i /lambda_bar = (P_i/P_bar)^(2/order)
    M                   = SpectralDecomp(nElems, 2)
    Pbarv               = @mlv sqrt(P1v*P2v)
    @mlv M.lambda_bar   = Pbarv^(2/2)
    @mlv M.sf[1]        = (P1v/Pbarv)^(2/2)
    @mlv M.sf[2]        = (P2v/Pbarv)^(2/2)
    @mlv M.V[1,1]       = cos(phiv)
    @mlv M.V[2,1]       = sin(phiv)
    @mlv M.V[1,2]       = -sin(phiv)
    @mlv M.V[2,2]       = cos(phiv)
    
    return M
    
end

#Plot target ellipses:
function AMA_PlotEllipses!(mesh::TrMesh, metric::SpectralDecomp{Vector{Float64},Int};
    color="b", linewidth=1.1, corr_factor::Float64=1.0)

    #Write metric as matrix:
    metricm     = GetMatrix(metric)
    
    #Get barycenters:
    xbar        = barycenters(mesh)
    
    #Radial unitary vector:
    thetav      = linspace(0.0, 2.0*pi, 361)
    ethetam     = [ cos.(thetav), sin.(thetav) ]
    
    #The radius in theta direction satisfies R(theta)^2 e_theta^T M e_theta = Rhat^2 = 1
    h_theta     = zeros(length(thetav))
    for iElem=1:mesh.nElems
        @. h_theta              = 0.0
        for ii=1:2, jj=1:2
            @. h_theta          += metricm[ii,jj][iElem] * ethetam[ii] * ethetam[jj]
        end
        @. h_theta              = corr_factor*sqrt(1.0/h_theta)
#         @. $plot(xbar[iElem,1]+h_theta*ethetam[1], xbar[iElem,2]+h_theta*ethetam[2], color=color)
        @. $plot(xbar[iElem,1]+h_theta*ethetam[1], xbar[iElem,2]+h_theta*ethetam[2])
    end

end
