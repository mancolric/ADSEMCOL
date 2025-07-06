function ConvectiveFlux!(model::SCD, a::Vector{MFloat}, da_du::Vector{MFloat},
    u::Vector{MFloat}, du::Matrix{MFloat}, 
    ComputeJ::Bool, flux::Matrix{MFloat}, dflux_du::Array{MFloat,3},  
        dflux_dgradu::Array{MFloat,4}) where MFloat<:Matrix{Float64}

    @mlv flux[1,1]  += a[1]*u[1] 
    @mlv flux[1,2]  += a[2]*u[1] 
    if ComputeJ
        @mlv dflux_du[1,1,1]        += da_du[1]*u[1] + a[1]
        @mlv dflux_du[1,2,1]        += da_du[2]*u[1] + a[2]
    end
    
    return
    
end

function DiffusiveFlux!(model::SCD, epsilon::MFloat, depsilon_du::MFloat,
    u::Vector{MFloat}, du::Matrix{MFloat}, 
    ComputeJ::Bool, flux::Matrix{MFloat}, dflux_du::Array{MFloat,3},  
        dflux_dgradu::Array{MFloat,4}) where MFloat<:Matrix{Float64}

    @mlv flux[1,1]  -= epsilon*du[1,1]
    @mlv flux[1,2]  -= epsilon*du[1,2]
    if ComputeJ #&& false
        @mlv dflux_du[1,1,1]        -= depsilon_du*du[1,1]
        @mlv dflux_du[1,2,1]        -= depsilon_du*du[1,2]
        @mlv dflux_dgradu[1,1,1,1]  -= epsilon
        @mlv dflux_dgradu[1,2,1,2]  -= epsilon
    end
    
    return
    
end
