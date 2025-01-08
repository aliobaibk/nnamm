function Y = nnamm_solve_ode(solver, eta, nts, h, param)
switch solver
    case 'heun-sto'
        Y = heun_sto(eta, nts, h, param);
end
end


%% ------------ Local functions ------------------------------------------------
function Y = heun_sto(eta, nts, h, param)
Y = [eta; zeros(nts - 1, size(eta, 2))].';

q_mu = param.q.mu.';

q_sigma = sparse(size(eta, 2), nts);
q_sigma(param.ode(:, 5), :) = param.q.sigma.' .* repmat(h * (param.a .* ...
    param.A), 1, nts);


% - Solution
for k_t = 1 : (nts - 1)
    k_1 = nnamm_velocity_vector(Y(:, k_t), q_mu(:, k_t), param);
    k_2 = nnamm_velocity_vector(Y(:, k_t) + h * k_1 + q_sigma(:, k_t), q_mu(:, ...
        k_t), param);
    Y(:, k_t + 1) = Y(:, k_t) + (h / 2) * (k_1 + k_2) + q_sigma(:, k_t);
end
Y = Y.';
end
