function compute_steady_state(prob, sys, ginfo, slack, p; p_sys = 75E5, max_iter=100, tol = 0.0001)
    a = GasNetModel.get_speed_of_sound(ginfo)
    rho = p / a^2
    rho_sys = p_sys / a^2

    prob[getproperty(sys, Symbol(:sub_, slack)).ρ] = 123456789
    for i = setdiff(1:size(ginfo.nodes,1), slack)
        prob[getproperty(sys, Symbol(:sub_,i)).ρ] = rho_sys
    end
        
    for i=1:size(ginfo.pipes,1)
        prob[getproperty(sys, Symbol(:pipe_,i)).ρ] = rho_sys * ones(length(getproperty(sys, Symbol(:pipe_,1)).ρ))
        prob[getproperty(sys, Symbol(:pipe_,i)).ϕ] = zeros(length(getproperty(sys, Symbol(:pipe_,1)).ϕ))

    end
    
    u = copy(prob.u0)
    id_slack = findfirst(u .== 123456789)
    u[id_slack] = rho
    id = setdiff(1:length(u), id_slack)
    
    for i=1:max_iter
        u[id] -= pinv(Matrix(prob.f.jac(u,prob.p,0.0)[id,id])) * prob.f(u,prob.p,0.0)[id]
        if sum(abs2, prob.f(u,prob.p,0.0)) / length(u) < tol
            break
        end
        if i == max_iter
            @warn "Reached max iter"
        end
    end
    prob = ModelingToolkit.remake(prob; u0 = u)
end
