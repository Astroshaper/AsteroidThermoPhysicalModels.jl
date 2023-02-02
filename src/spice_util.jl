
"""
Obtain a dataframe of ephemerides 

# Arguments
- `targ`   : Target body name
- `ets`    : AbstractVetor of observer epoch
- `ref`    : Reference frame of output position vector
- `abcorr` : Aberration correction flag
- `obs`    : Observing body name
"""
function spkpos_df(targ, ets::AbstractVector, ref, abcorr, obs)
    Base.depwarn("`spkpos_df` is deprecated. Please use `SPICE.spkpos` and `DataFrame` instead.", :spkpos_df)
    df = DataFrame(et=Float64[], x=Float64[], y=Float64[], z=Float64[], lt=Float64[])
    for et in ets
        pos, lt = SPICE.spkpos(targ, et, ref, abcorr, obs)  # pos [km], lt [s]
        pos = SPICE.convrt.(pos, "km", "m")
        push!(df, (et, pos..., lt))
    end
    df
end


"""
Obtain a vector of position 

# Arguments
- `targ`   : Target body name
- `ets`    : AbstractVetor of observer epoch
- `ref`    : Reference frame of output position vector
- `abcorr` : Aberration correction flag
- `obs`    : Observing body name
"""
function spkpos(targ, ets::AbstractVector, ref, abcorr, obs)
    Base.depwarn("`spkpos` is deprecated. Please use `[SPICE.spkpos($targ, et, $ref, $abcorr, $obs)[1]*1000 for et in $ets]` instead.", :spkpos)
    positions = Vector{Float64}[]
    for et in ets
        pos, lt = SPICE.spkpos(targ, et, ref, abcorr, obs)  # pos [km], lt [s]
        pos = SPICE.convrt.(pos, "km", "m")
        push!(positions, collect(pos))
    end
    positions
end


"""
# Arguments
- `from` : Name of the frame to transform from
- `to`   : Name of the frame to transform to
- `ets`  : Epoch of the rotation matrix

# Return
- `rotate` : A rotation matrix
"""
function pxform(from, to, ets)
    Base.depwarn("`pxform` is deprecated. Please use `[SPICE.pxform($from, $to, et) for et in $ets]` instead.", :pxform)
    return [SPICE.pxform(from, to, et) for et in ets]
end
