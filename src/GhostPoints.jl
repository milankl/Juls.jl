""" Extends the matrices u,v,η,sst with a halo of ghost points for boundary conditions."""
function add_halo(  P::Parameter,
                    C::Constants,
                    G::Grid,
                    u::Array{T,2},
                    v::Array{T,2},
                    η::Array{T,2},
                    sst::Array{T,2}) where {T<:AbstractFloat}

    @unpack nx,ny,nux,nuy,nvx,nvy,halo,haloη,halosstx,halossty = G

    # Add zeros to satisfy kinematic boundary conditions
    u = cat(zeros(T,halo,nuy),u,zeros(T,halo,nuy),dims=1)
    u = cat(zeros(T,nux+2*halo,halo),u,zeros(T,nux+2*halo,halo),dims=2)

    v = cat(zeros(T,halo,nvy),v,zeros(T,halo,nvy),dims=1)
    v = cat(zeros(T,nvx+2*halo,halo),v,zeros(T,nvx+2*halo,halo),dims=2)

    η = cat(zeros(T,haloη,ny),η,zeros(T,haloη,ny),dims=1)
    η = cat(zeros(T,nx+2*haloη,haloη),η,zeros(T,nx+2*haloη,haloη),dims=2)

    sst = cat(zeros(T,halosstx,ny),sst,zeros(T,halosstx,ny),dims=1)
    sst = cat(zeros(T,nx+2*halosstx,halossty),sst,zeros(T,nx+2*halosstx,halossty),dims=2)

    ghost_points!(P,C,u,v,η)
    ghost_points_sst!(P,G,sst)
    return u,v,η,sst
end

""" Copy ghost points for u from inside to the halo in the nonperiodic case. """
function ghost_points_u_nonperiodic!(C::Constants,u::AbstractMatrix)

    @unpack one_minus_α = C

    # tangential boundary condition
    @views @inbounds u[:,1] .= one_minus_α*u[:,4]
    @views @inbounds u[:,2] .= one_minus_α*u[:,3]
    @views @inbounds u[:,end-1] .= one_minus_α*u[:,end-2]
    @views @inbounds u[:,end] .= one_minus_α*u[:,end-3]
end

""" Copy ghost points for u from inside to the halo in the periodic case. """
function ghost_points_u_periodic!(C::Constants,u::AbstractMatrix)

    @unpack one_minus_α = C

    # periodic bc
    @views @inbounds u[1,:] .= u[end-3,:]
    @views @inbounds u[2,:] .= u[end-2,:]
    @views @inbounds u[end-1,:] .= u[3,:]
    @views @inbounds u[end,:] .= u[4,:]

    # tangential bc
    @views @inbounds u[:,1] .= one_minus_α*u[:,4]
    @views @inbounds u[:,2] .= one_minus_α*u[:,3]
    @views @inbounds u[:,end-1] .= one_minus_α*u[:,end-2]
    @views @inbounds u[:,end] .= one_minus_α*u[:,end-3]
end

""" Copy ghost points for v from inside to the halo in the nonperiodic case. """
function ghost_points_v_nonperiodic!(C::Constants,v::AbstractMatrix)

    @unpack one_minus_α = C

    # tangential boundary condition
    @views @inbounds v[1,:] .= one_minus_α*v[4,:]
    @views @inbounds v[2,:] .= one_minus_α*v[3,:]
    @views @inbounds v[end-1,:] .= one_minus_α*v[end-2,:]
    @views @inbounds v[end,:] .= one_minus_α*v[end-3,:]
end

""" Copy ghost points for v from inside to the halo in the periodic case. """
function ghost_points_v_periodic!(v::AbstractMatrix)

    # tangential boundary condition
    @views @inbounds v[1,:] .= v[end-3,:]
    @views @inbounds v[2,:] .= v[end-2,:]
    @views @inbounds v[end-1,:] .= v[3,:]
    @views @inbounds v[end,:] .= v[4,:]
end

""" Copy ghost points for η from inside to the halo in the nonperiodic case. """
function ghost_points_η_nonperiodic!(η::AbstractMatrix)

    # assume no gradients of η across solid boundaries
    # the 4 corner points are copied twice, but it's faster!
    @views @inbounds η[1,:] .= η[2,:]
    @views @inbounds η[end,:] .= η[end-1,:]

    @views @inbounds η[:,1] .= η[:,2]
    @views @inbounds η[:,end] .= η[:,end-1]
end

""" Copy ghost points for η from inside to the halo in the periodic case. """
function ghost_points_η_periodic!(η::AbstractMatrix)

    # corner points are copied twice, but it's faster!
    @views @inbounds η[1,:] .= η[end-1,:]
    @views @inbounds η[end,:] .= η[2,:]

    @views @inbounds η[:,1] .= η[:,2]
    @views @inbounds η[:,end] .= η[:,end-1]
end

""" Copy ghost points for η from inside to the halo in the nonperiodic case. """
function ghost_points_sst_nonperiodic!(G::Grid,sst::AbstractMatrix)

    @unpack halosstx,halossty = G

    # assume no gradients of η across solid boundaries
    # the 4 corner points are copied twice, but it's faster!
    for i ∈ 1:halosstx
        @views @inbounds sst[i,:] .= sst[halosstx+1,:]
        @views @inbounds sst[end-i+1,:] .= sst[end-halosstx,:]
    end

    for j ∈ 1:halossty
        @views @inbounds sst[:,j] .= sst[:,halossty+1]
        @views @inbounds sst[:,end-j+1] .= sst[:,end-halossty]
    end
end

""" Copy ghost points for η from inside to the halo in the periodic case. """
function ghost_points_sst_periodic!(G::Grid,sst::AbstractMatrix)

    @unpack halosstx,halossty = G

    # corner points are copied twice, but it's faster!
    for i ∈ 1:halosstx
        @views @inbounds sst[i,:] .= sst[end-2*halosstx+i,:]
        @views @inbounds sst[end-halosstx+i,:] .= sst[halosstx+i,:]
    end

    for j ∈ 1:halossty
        @views @inbounds sst[:,j] .= sst[:,halossty+1]
        @views @inbounds sst[:,end-j+1] .= sst[:,end-halossty]
    end
end

"""Decide on boundary condition P.bc which ghost point function to execute."""
function ghost_points!( P::Parameter,
                        C::Constants,
                        u::AbstractMatrix,
                        v::AbstractMatrix,
                        η::AbstractMatrix)

    if P.bc == "periodic"
        ghost_points_u_periodic!(C,u)
        ghost_points_v_periodic!(v)
        ghost_points_η_periodic!(η)
    else
        ghost_points_u_nonperiodic!(C,u)
        ghost_points_v_nonperiodic!(C,v)
        ghost_points_η_nonperiodic!(C,η)
    end
end

"""Decide on boundary condition P.bc which ghost point function to execute."""
function ghost_points_sst!( P::Parameter,
                            G::Grid,
                            sst::AbstractMatrix)

    if P.bc == "periodic"
        ghost_points_sst_periodic!(G,sst)
    else
        ghost_points_sst_nonperiodic!(G,sst)
    end
end