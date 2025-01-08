function [Y, t, eta_i, eta_f, param] = nnamm_phase_trajectory(I, h, ...
    system_param, solver_param, flag_init)
[t, eta_i, param] = nnamm_init_system(I, h, system_param, solver_param, ...
    flag_init);

if flag_init
    Y = [];
    eta_f = eta_i;
    return
end

Y = nnamm_solve_ode(solver_param.ode_solver, eta_i, numel(t), h, param);
eta_f = Y(end, :);
end
