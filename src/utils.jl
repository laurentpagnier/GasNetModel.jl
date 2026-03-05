function convert(x::Real, to=:mmBtu; btu_p_kwh=3412.14, MJ_p_kg=55.0)
    # assumes that the input x in model units
    if to == :mmBtu 
        return x * MJ_p_kg / 3600 * btu_p_kwh / 1E3
    elseif to == :MWh
        return x * MJ_p_kg / 3600 
    end
end
