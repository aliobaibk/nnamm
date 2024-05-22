function [t, eta, param] = nnamm_init_system(I, h, system_param, ...
    solver_param, flag_init)
t = (I(1) : h : I(2)).';

param = init_param(system_param, t);

eta = init_first_state(solver_param, param, flag_init);
end


%% ------------ Local functions ------------------------------------------------
function param = init_param(i_param, t)
nodes = i_param.nodes;

% ode(:, i) is the indices of the i-th ODE for all neuron-astrocyte mass models
param.ode = reshape(1 : (14 * numel(nodes)), numel(nodes), 14);

param.Omega_Pyr = i_param.omega_Pyr * i_param.Omega_Pyr;
param.Omega_Ast = i_param.Omega_Ast;

param.A = vertcat(nodes.A);
param.B = vertcat(nodes.B);
param.a = vertcat(nodes.a);
param.b = vertcat(nodes.b);

param.C_Pyr_to_ExIn = vertcat(nodes.C_Pyr_to_ExIn);
param.C_ExIn_to_Pyr = vertcat(nodes.C_ExIn_to_Pyr);
param.C_Pyr_to_InIn = vertcat(nodes.C_Pyr_to_InIn);
param.C_InIn_to_Pyr = vertcat(nodes.C_InIn_to_Pyr);
param.C_Pyr_to_Pyr = vertcat(nodes.C_Pyr_to_Pyr);

param.SN_max = repmat(vertcat(nodes.nu_max), [1, 3]);
param.SN_rate = repmat(vertcat(nodes.r), [1, 3]);
param.SN_mode = repmat(vertcat(nodes.v_0), [1, 3]);
param.SN_off = zeros(size(param.SN_max));

w_r = vertcat(nodes.w_r);
w_d = vertcat(nodes.w_d);
W_Ast_orig = vertcat(nodes.W_Ast);
param.W_Ast = i_param.omega_Glu * (W_Ast_orig .* w_r);
param.W_Pyr = vertcat(nodes.W_Pyr) .* w_r;
param.w_s = w_r + w_d;
param.w_p = w_r .* w_d;
param.tau_Glu_Ast = vertcat(nodes.tau_Glu_Ast);

z_r = vertcat(nodes.z_r);
z_d = vertcat(nodes.z_d);
Z_Ast_orig = vertcat(nodes.Z_Ast);
param.Z_Ast = i_param.omega_GABA * (Z_Ast_orig .* z_r);
param.Z_InIn = vertcat(nodes.Z_InIn) .* z_r;
param.z_s = z_r + z_d;
param.z_p = z_r .* z_d;
param.tau_GABA_Ast = vertcat(nodes.tau_GABA_Ast);

param.SE_max = horzcat( ...
    vertcat(nodes.m_Glu_Pyr), ...
    vertcat(nodes.m_Glu_InIn), ...
    vertcat(nodes.m_GABA_Pyr), ...
    vertcat(nodes.V_Glu_e_to_Ast), ...
    vertcat(nodes.V_Glu_e_to_Pyr), ...
    vertcat(nodes.m_Glu_Ast));
param.SE_rate = horzcat( ...
    vertcat(nodes.r_Glu_Pyr), ...
    vertcat(nodes.r_Glu_InIn), ...
    vertcat(nodes.r_GABA_Pyr), ...
    vertcat(nodes.r_Glu_e_to_Ast), ...
    vertcat(nodes.r_Glu_e_to_Pyr), ...
    vertcat(nodes.r_Glu_Ast));
param.SE_mode = horzcat( ...
    vertcat(nodes.theta_Glu_Pyr), ...
    vertcat(nodes.theta_Glu_InIn), ...
    vertcat(nodes.theta_GABA_Pyr), ...
    vertcat(nodes.theta_Glu_e_to_Ast), ...
    vertcat(nodes.theta_Glu_e_to_Pyr), ...
    vertcat(nodes.theta_Glu_Ast));
param.SE_off = horzcat( ...
    vertcat(nodes.delta_Glu_Pyr), ...
    vertcat(nodes.delta_Glu_InIn), ...
    vertcat(nodes.delta_GABA_Pyr), ...
    vertcat(nodes.delta_Glu_e_to_Ast), ...
    vertcat(nodes.delta_Glu_e_to_Pyr), ...
    vertcat(nodes.delta_Glu_Ast));

param.MM_v_max = horzcat(vertcat(nodes.V_GABA_e_to_Ast), ...
    vertcat(nodes.V_GABA_e_to_InIn));
param.MM_K = horzcat(vertcat(nodes.K_GABA_e_to_Ast), ...
    vertcat(nodes.K_GABA_e_to_InIn));

param.q = init_q({nodes.q}, t);

% Parameters non trivial to recover
param.omega_Pyr = i_param.omega_Pyr;
param.omega_Glu = i_param.omega_Glu;
param.W_Ast_orig = W_Ast_orig;
param.w_r = w_r;
param.omega_GABA = i_param.omega_GABA;
param.Z_Ast_orig = Z_Ast_orig;
param.z_r = z_r;
end


function q = init_q(q_param, t)
q = repmat(struct('mu', zeros(numel(t), 1), ...
    'sigma', zeros(numel(t), 1)), numel(q_param), 1);
for k = 1 : numel(q_param)
    switch q_param{k}.method
        case 'gaussian'
            q(k) = init_gaussian_q(q_param{k}, t);
    end
end
q = struct('mu', cat(2, q.mu), 'sigma', cat(2, q.sigma));
end


function q = init_gaussian_q(q_param, t)
% Multiplication trick for q.mu to ensure correct output dimensions
q = struct('mu', q_param.mean(:) .* ones(numel(t), 1), ...
    'sigma', q_param.std * randn(numel(t), 1));
end


function eta = init_first_state(solver_param, param, flag_init)
if flag_init
    eta = zeros(1, param.ode(end));
    return
end
switch solver_param.init_method
    case 'load'
        eta = solver_param.eta;
    case 'zeros'
        eta = zeros(1, param.ode(end));
end
eta = eta(end, :);
end
