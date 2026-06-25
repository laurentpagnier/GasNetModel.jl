function convert(x::Real, to=:mmBtu; btu_p_kwh=3412.14, MJ_p_kg=55.0)
    # assumes that the input x in model units
    if to == :mmBtu 
        return x * MJ_p_kg / 3600 * btu_p_kwh / 1E3
    elseif to == :MWh
        return x * MJ_p_kg / 3600 
    end
end


function GasInfo(;folder="../data/", pipe_fn=nothing, node_fn = nothing, param_fn = nothing, compressor_fn = nothing, 
    json_fn = nothing)
    if isnothing(json_fn)
        # parsing assumes a given structure
        pipes = DataFrame(CSV.File(isnothing(pipe_fn) ? joinpath(folder, "pipes.csv") : pipe_fn, skipto=3, header=2))
        params = DataFrame(CSV.File(isnothing(param_fn) ? joinpath(folder, "params.csv") : param_fn, skipto=3, header=2))
        nodes = DataFrame(CSV.File(isnothing(node_fn) ? joinpath(folder, "nodes.csv") : params, skipto=3, header=2))
        compressors = DataFrame(CSV.File(isnothing(compressor_fn) ? joinpath(folder, "compressors.csv") : compressor_fn, skipto=3, header=2))
        #node_ic = DataFrame(CSV.File(joinpath(folder, "initial_conditions/node_ic.csv"), skipto=3, header=2))
        #pipe_ic = DataFrame(CSV.File(joinpath(folder, "initial_conditions/pipe_ic.csv"), skipto=3, header=2))
        #comp_ic = DataFrame(CSV.File(joinpath(folder, "initial_conditions/compressor_ic.csv"), skipto=3, header=2))
        #nodes = leftjoin(nodes, node_ic, on = :number)
        #pipes = leftjoin(pipes, pipe_ic, on = :number)
        # TODO compressors are not yet implemented
        #compressors = leftjoin(compressors, comp_ic, on = :number)
        return GasInfo(pipes, nodes, compressors, params)
    else
        return parse_json(json_fn)
    end
end


function parse_json(fn)
    temp = JSON.parsefile(fn)
    nodes = DataFrame()
    for k in keys(temp["nodes"])
        row = Dict()
        for v in keys(temp["nodes"][k])
            row[v] = temp["nodes"][k][v]
        end
        push!(nodes, row, cols=:union)
    end
    
    pipes = DataFrame()
    for k in keys(temp["pipes"])
        row = Dict()
        for v in keys(temp["pipes"][k])
            row[v] = temp["pipes"][k][v]
        end
        push!(pipes, row, cols=:union)
    end
    
    compressors = DataFrame()
    #=
    for k in keys(temp["compressors"])
        row = Dict()
        for v in keys(temp["compressors"][k])
            row[v] = temp["compressors"][k][v]
        end
        push!(compressors, row, cols=:union)
    end
    =#
    
    params = DataFrame( Dict( k => temp["params"][k] for k in keys(temp["params"])) )
    GasInfo(sort!(pipes, :number), sort!(nodes, :number), compressors, params)
end
