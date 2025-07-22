export GasInfo, GasSystem, q_nodal

struct GasInfo
    pipes::DataFrame
    nodes::DataFrame
    compressors::DataFrame
    params::DataFrame
end

function GasInfo(folder="../data/")
    # parsing assumes a given structure
    pipes = DataFrame(CSV.File(joinpath(folder, "network_data/pipes.csv"), skipto=3, header=2))
    params = DataFrame(CSV.File(joinpath(folder, "network_data/params.csv"), skipto=3, header=2))
    nodes = DataFrame(CSV.File(joinpath(folder, "network_data/nodes.csv"), skipto=3, header=2))
    compressors = DataFrame(CSV.File(joinpath(folder, "network_data/compressors.csv"), skipto=3, header=2))
    node_ic = DataFrame(CSV.File(joinpath(folder, "initial_conditions/node_ic.csv"), skipto=3, header=2))
    pipe_ic = DataFrame(CSV.File(joinpath(folder, "initial_conditions/pipe_ic.csv"), skipto=3, header=2))
    comp_ic = DataFrame(CSV.File(joinpath(folder, "initial_conditions/compressor_ic.csv"), skipto=3, header=2))
    nodes = leftjoin(nodes, node_ic, on = :number)
    pipes = leftjoin(pipes, pipe_ic, on = :number)
    # TODO compressors are not yet implemented
    #compressors = leftjoin(compressors, comp_ic, on = :number)
    GasInfo(pipes, nodes, compressors, params)
end


function get_speed_of_sound(ginfo)
    γ = ginfo.params.specific_heat_capacity_ratio[1]
    T = ginfo.params.temperature[1]
    Z = ginfo.params.compressibility_factor[1]
    R = 8.31446261815324 / 0.01604 # molar R / (methane) molar mass
    a = sqrt(γ * Z * R * T)
end


#for more details equations and diagrams, refer to https://arxiv.org/abs/2310.18507

@mtkmodel Pipe begin
    @description "Model a gas pipeline"

    @structural_parameters begin
        L = 1
        dx = 1
        N = ceil(Int, L/dx)
        p_o = 1
        p_i = 1
        ϕ0 = 1
        a = 450
        β = 0.001
    end

    @variables begin
        #staggered grid
        ρ(t)[1:N] = collect(range(p_o, p_i, N)) / a^2
        ϕ(t)[1:(N + 1)] = ϕ0
    end
    
    @equations begin
        [D(ϕ[i]) ~ -(a^2/dx)*(ρ[i] - ρ[i-1]) - 2β*(ϕ[i]  * abs(ϕ[i])) / (ρ[i] + ρ[i-1]) for i=2:N] ∪ # i=1 and i=N+1 are done later
        [D(ρ[i]) ~ -(1 / dx) * (ϕ[i+1] - ϕ[i]) for i=1:N]
    end
end

@register_symbolic q_nodal(i,t)

@mtkmodel Substation begin
    @description "Model a node"

    @structural_parameters begin
        q0 = 1
        p = 1
        a = 450
        i = 1
    end

    @variables begin
        q(t)
        ρ(t) = p / a^2
    end
    
    @equations begin
        q ~ q_nodal(i,t)
    end 
end


function create_equations(sub, pipe, ginfo, dx)
    n_out = ginfo.pipes.start_node 
    n_in = ginfo.pipes.end_node
    diam = ginfo.pipes.diameter
    volumes = zeros(length(sub)) # volume of the node cell
    flows = Num.(zeros(length(sub))) # mass flow entering or exiting the cell
    eqs = Equation[]
    a = get_speed_of_sound(ginfo)
    βs = [p.friction_factor/ p.diameter/2 for p in eachrow(ginfo.pipes)]
    Ss  = [pi*p.diameter^2/4 for p in eachrow(ginfo.pipes)]
    # compute the time-derivative of first and last pipe fluxes using nodal densities
    # and define quantities for the mass conservation eq
    for (k, (n_i, n_o, d, β, S)) in enumerate(zip(n_in, n_out, diam, βs, Ss))
        volumes[n_o] += S * dx/2 
        flows[n_o] = flows[n_o] - S*pipe[k].ϕ[1]
        ϕ, ρ, ρn = pipe[k].ϕ[1], pipe[k].ρ[1], sub[n_o].ρ 
        eq = D(ϕ) ~ -(a^2/dx)*(ρ  - ρn) - 2β*(ϕ  * abs(ϕ)) / (ρ + ρn) 
        push!(eqs, eq)
        
        volumes[n_i] += S * dx/2 
        flows[n_i] = flows[n_i] + S*pipe[k].ϕ[end]
        ϕ, ρ, ρn = pipe[k].ϕ[end], pipe[k].ρ[end], sub[n_i].ρ 
        eq = D(ϕ) ~ -(a^2/dx)*(ρn - ρ) - 2β*(ϕ  * abs(ϕ)) / (ρ + ρn) 
        push!(eqs, eq)
    end
    # mass conservation at the node cell
    eqs ∪ [D(n.ρ) ~ (n.q + f)/ v for (n, f, v) in zip(sub, flows, volumes)]
end

@mtkmodel GasSystem begin
    @description "Connect the different component of which the system consists."
    
    @structural_parameters begin
        ginfo
        dx = 10_000
    end
    
    @components begin
        sub = [Substation(
                    q0 = ginfo.nodes.initial_nodal_flow[i],
                    p = ginfo.nodes.initial_nodal_pressure[i],
                    a = get_speed_of_sound(ginfo),
                    i = i,
                ) for i=1:size(ginfo.nodes,1)] 
        pipe = [Pipe(
                    L = ginfo.pipes.length[i],
                    dx = dx,
                    β = ginfo.pipes.friction_factor[i]/ginfo.pipes.diameter[i]/2,
                    a = get_speed_of_sound(ginfo),
                    p_i = ginfo.pipes.initial_pipe_pressure_in[i],
                    p_o = ginfo.pipes.initial_pipe_pressure_out[i],
                    ϕ0 = 4*ginfo.pipes.initial_pipe_flow[i] / ginfo.pipes.diameter[i]^2 / pi ,
                ) for  i=1:size(ginfo.pipes,1)]
    end

    @equations begin
        # conservation of mass in the nodal cell
        # and flux at the boundary
        create_equations(sub, pipe, ginfo, dx)
    end
end
