export GasInfo, GasSystem, q_nodal, create_gas_sys

@register_symbolic q_nodal(i,t)

struct GasInfo
    pipes::DataFrame
    nodes::DataFrame
    compressors::DataFrame
    params::DataFrame
end


function compute_linepack(ginfo, sys, sol, dx = 10_000)
    # linepack is kg of natural gas
    linepack = zeros(length(sol.t))
    global temp, id = sol, sys # need to be global to be called in Meta parse
    for (i, p) in enumerate(eachrow(ginfo.pipes))
        id1 = p.start_node 
        id2 = p.end_node 
        eval(Meta.parse("dens = temp[[id.sub_$id1.ρ; id.pipe_$i.ρ; id.sub_$id2.ρ]]"))
        section = pi * p.diameter^2 / 4
        volume = section * dx
        l = map(x-> x[1]/2 + sum(x[2:end-1]) + x[end]/2 , dens) * volume
        linepack += l
    end
    linepack
end


function get_speed_of_sound(ginfo)
    γ = ginfo.params.specific_heat_capacity_ratio[1]
    T = ginfo.params.temperature[1]
    Z = ginfo.params.compressibility_factor[1]
    R = 8.31446261815324 / 0.01604 # molar R / (methane) molar mass
    a = sqrt(γ * Z * R * T)
end


#for more details equations and diagrams, refer to https://arxiv.org/abs/2310.18507

function make_pipe(name; L = 2000, dx = 1000, p_o = 75E5, p_i = 75E5, ϕ0 = 0, a = 450, β = 0.001, grav = 0)
    N = ceil(Int, L/dx)
    @variables ρ(t)[1:N] = collect(range(p_o, p_i, N)) / a^2 ϕ(t)[1:(N + 1)] = fill(ϕ0, N+1)
    eqs = [
        [D(ϕ[i]) ~ -(a^2/dx)*(ρ[i] - ρ[i-1]) - 2*β*(ϕ[i]  * abs(ϕ[i])) / (ρ[i] + ρ[i-1]) - 0.5*(ρ[i] + ρ[i-1])*grav for i=2:N]...
        [D(ρ[i]) ~ -(1 / dx) * (ϕ[i+1] - ϕ[i]) for i=1:N]...
    ]
    System(eqs, t; name)
end

function make_substation(name; q0 = 0, p = 75E5, a = 450, i = 1)
    @variables ρ(t) = p / a^2
    @parameters i = i
    #eqs = q ~ q_nodal(i,t)
    System(Equation[], t, [ρ], [i]; name)
end


function create_gas_sys(ginfo; name, dx = 10_000)
    a = GasNetModel.get_speed_of_sound(ginfo)
    subs = System[]
    for (i, n) in enumerate(eachrow(ginfo.nodes))
        push!(subs, make_substation(Symbol(:sub_, i), a=a, i=i))
    end
    
    eqs = Equation[]
    volumes = zeros(length(subs)) # volume of the node cell
    flows = Num.(zeros(length(subs))) # mass flow entering or exiting the cell
    pipes = System[]
    for (k, p) in  enumerate(eachrow(ginfo.pipes))
        n_out = p.start_node 
        n_in = p.end_node
        diam = p.diameter
        L = p.length
        β = p.friction_factor / diam / 2 
        S = pi * p.diameter^2 / 4
        grav = 9.81 * (ginfo.nodes.height[n_in] - ginfo.nodes.height[n_out]) / L
        p_i = ginfo.nodes.default_nodal_pressure[n_in]
        p_o = ginfo.nodes.default_nodal_pressure[n_out]
        
        push!(pipes, make_pipe(Symbol(:pipe_, k), L=L, a=a, β=β, grav=grav, p_i=p_i, p_o=p_o, dx=10_000))
        volumes[n_out] += S * dx/2 
        volumes[n_in] += S * dx/2 
        flows[n_in] = flows[n_in] + S*pipes[k].ϕ[end]
        flows[n_out] = flows[n_out] - S*pipes[k].ϕ[1]
        ϕ, ρ, ρn = pipes[k].ϕ[1], pipes[k].ρ[1], subs[n_out].ρ
        eq = D(ϕ) ~ -a^2/dx * (ρ  - ρn) - 2*β*(ϕ  * abs(ϕ)) / (ρ + ρn) - 0.5*grav*(ρ + ρn)
        push!(eqs, eq)
        ϕ, ρ, ρn = pipes[k].ϕ[end], pipes[k].ρ[end], subs[n_in].ρ 
        push!(eqs, D(ϕ) ~ -a^2/dx * (ρn - ρ) - 2*β*(ϕ  * abs(ϕ)) / (ρ + ρn) - 0.5*grav*(ρ + ρn))
    end

    eqs = eqs ∪ [D(n.ρ) ~ (q_nodal(n.i,t) + f)/ v for (n, f, v) in zip(subs, flows, volumes)]
    sys = System(eqs, t; systems = [subs; pipes], name = :gas_sys)
end




