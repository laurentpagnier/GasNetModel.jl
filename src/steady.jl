function compute_steady_state(prob, ginfo, slack, p; p_sys = 75E5, max_iter=100, tol = 0.0001)
    a = GasNetModel.get_speed_of_sound(ginfo)
    global rho = p / a^2
    global rho_sys = p_sys / a^2
    global cthulhu = deepcopy(prob) # hopefully this var name isn't used elsewhere or Cthulhu's wrath will be unleashed
   
    # zeros flow in the pipes
    for i=1:size(ginfo.nodes,1)
        eval(Meta.parse("cthulhu[sys.sub_$slack.ρ] = rho"))
    end
    
    for i=1:size(ginfo.pipes,1)
        eval(Meta.parse("temp = length(cthulhu[sys.pipe_$i.ϕ])"))
        eval(Meta.parse("cthulhu[sys.pipe_$i.ϕ] = zeros(temp)"))
        eval(Meta.parse("cthulhu[sys.pipe_$i.ρ] = rho_sys * ones(temp)"))
    end
    eval(Meta.parse("cthulhu[sys.sub_$(slack).ρ] = 123456789"))
    u = copy(cthulhu.u0)
    id_slack = findfirst(u .== 123456789)
    u[id_slack] = rho
    id = setdiff(1:length(u), id_slack)
    for i=1:max_iter
        u[id] -= pinv(Matrix(prob.f.jac(u,[],0.0)[id,id])) * prob.f(u,[],0.0)[id]
        if sum(abs2, prob.f(u,[],0.0)) / length(u) < tol
            break
        end
        if i == max_iter
            @warn "Reached max iter"
        end
    end
    u
end
