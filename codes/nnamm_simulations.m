function nnamm_simulations()
% Simulations for Pub[astrocyte-neuron network dialogue mechanisms][PCBI-2025]
% NOTES:
% - Be aware of the custom saving for the simulations, which was implemented to
%   speed up post-processing from slow drives.
% - When resampling the simulated signals, noise components are left untouched.
% - Monitor the machine's RAM during simulations. If memory is insufficient, you
%   can temporarily save data to disk and resume simulations later. Use the
%   saved final states from the previous integration as the initial states for
%   the new integration. Follow a strategy similar to that on lines 239–250.

proj_dir = './data+codes_pcbi2025_aliobaibk';
if ~isfolder(proj_dir)
    error(['The "proj_dir" variable is not set to an existing directory.\n', ...
        'Ensure "proj_dir" points to the directory where the repository was ', ...
        'cloned.\nFor example, if you cloned the repository into ', ...
        '"./data+codes_pcbi2025_aliobaibk",\n', ...
        'use: proj_dir = ''./data+codes_pcbi2025_aliobaibk''.'], '')
end

addpath(genpath(fullfile(proj_dir, 'codes', 'util')))
addpath(genpath(fullfile(proj_dir, 'codes', 'functions')))
addpath(genpath(fullfile(proj_dir, 'codes', 'nnamm')))

fprintf(['\n', ...
    '*** Dialogue mechanisms between astrocytic and neuronal networks ***\n\n'])


% - Outputs
o_opt.flag_resume = true; % keeps existing files (without checking) if true.
o_opt.simul_dir = fullfile(proj_dir, 'data', 'simulations');

o_opt.parameters = fullfile(o_opt.simul_dir, 'parameters.mat');
o_opt.grid = fullfile(o_opt.simul_dir, 'grid.mat');
o_opt.names = fullfile(o_opt.simul_dir, 'names.xls');


% - Connectomes (files need to be consistent, no checking is done)
c_opt.roi_name = 'laus18-3';
c_opt.net_name = 'lobes';
c_opt.hier_file = fullfile(proj_dir, 'data', 'connectomes', ...
    sprintf('hierarchy.%s.%s.xls', c_opt.roi_name, c_opt.net_name));

c_opt.neu_mat = fullfile(proj_dir, 'data', 'connectomes', ...
    'HCP-average_laus18-3_sc.csv');
c_opt.ast_mat = fullfile(proj_dir, 'data', 'connectomes', ...
    'ICBM152A09C_laus18-3_mid-fs-311558V_weights.csv');


% - Simulations
s_opt.name = 'papaya';
s_opt.n_rep = 10;
s_opt.I_init = [0; 360 + 10; 360]; % start; end; discard from start
s_opt.I = [0; 120; 0];
s_opt.ode_solver = 'heun-sto';
s_opt.h = 2^(-12); % 2^(-14);
s_opt.fs_save = 2^8; % no resampling if empty


% - Parameters
setup_investigation(o_opt)


% - Main
main(c_opt, s_opt, o_opt)
end



%% ------------ Local functions ------------------------------------------------
function setup_investigation(o_opt)
if (~o_opt.flag_resume || any(~isfile({o_opt.parameters, o_opt.grid, ...
        o_opt.names})))
    create_dir(o_opt.simul_dir, false, false)
end


% - For defining baseline
v_Glu_b = -0.34;
q_v_GABA_b = -60;
omega_Glu_b = 0.01;
omega_GABA_b = 0.01;


% - For defining exploration grid
v_Glu = [-0.3, 0.15];
q_v_GABA = [-50, 50];
n_omega_Glu = 35;
n_omega_GABA = 35;


% - Main parameterization
fprintf('>>> Parameterization through bifurcation analysis...\n')
P = [];
if (~o_opt.flag_resume || ~isfile(o_opt.parameters))
    P.omega_Pyr = 7.5;

    P.A = 3.25;
    P.B = 22;
    P.a = 100;
    P.b = 50;
    P.nu_max = 5;
    P.r = 0.56;
    P.v_0 = 6;
    P.C_Pyr_to_ExIn = 135;
    P.C_ExIn_to_Pyr = 108;
    P.C_Pyr_to_InIn = 33.75;
    P.C_InIn_to_Pyr = 33.75;
    P.C_Pyr_to_Pyr = 0;

    P.w_r = 90;
    P.w_d = 33;

    P.V_Glu_e_to_Ast = 4.5;
    P.r_Glu_e_to_Ast = 0.5;
    P.theta_Glu_e_to_Ast = 9;
    P.delta_Glu_e_to_Ast = 0;

    P.m_Glu_Ast = P.V_Glu_e_to_Ast;
    P.r_Glu_Ast = P.r_Glu_e_to_Ast;
    P.theta_Glu_Ast = P.theta_Glu_e_to_Ast;
    P.delta_Glu_Ast = P.delta_Glu_e_to_Ast;

    P.V_Glu_e_to_Pyr = 0.5;
    P.r_Glu_e_to_Pyr = 0.5;
    P.theta_Glu_e_to_Pyr = 9;
    P.delta_Glu_e_to_Pyr = 0;

    P.tau_Glu_Ast = 1 / 9;

    P.m_Glu_Pyr = 0.8;
    P.r_Glu_Pyr = 0.5;
    P.theta_Glu_Pyr = 10;
    P.delta_Glu_Pyr = -0.4;

    P.mu_Glu_InIn_by_Pyr = 0.5;

    P.m_Glu_InIn = P.mu_Glu_InIn_by_Pyr * P.m_Glu_Pyr;
    P.r_Glu_InIn = P.r_Glu_Pyr;
    P.theta_Glu_InIn = P.theta_Glu_Pyr;
    P.delta_Glu_InIn = P.mu_Glu_InIn_by_Pyr * P.delta_Glu_Pyr;

    P.z_r = 90;
    P.z_d = 33;

    P.V_GABA_e_to_Ast = 2;
    P.K_GABA_e_to_Ast = 8;

    P.V_GABA_e_to_InIn = 5;
    P.K_GABA_e_to_InIn = 24;

    P.tau_GABA_Ast = 1 / 9;

    P.m_GABA_Pyr = 4.2;
    P.r_GABA_Pyr = 0.25;
    P.theta_GABA_Pyr = 20;
    P.delta_GABA_Pyr = -2.1;

    P.q = struct('method', 'gaussian', 'mean', 240, 'std', 10);

    v_GABA_b = q_v_GABA_b * P.A / P.a;
    P.baseline = namm_setup_baseline(P, v_GABA_b, omega_GABA_b, v_Glu_b, ...
        omega_Glu_b);

    P.W_Ast = P.baseline.W;
    P.W_Pyr = P.baseline.W;

    P.Z_Ast = P.baseline.Z;
    P.Z_InIn = P.baseline.Z;

    P = orderfields(P);
    save(o_opt.parameters, '-struct', 'P', '-v7.3')
end


% - Exploration
fprintf('>>> Exploration parameter space...\n')
G = [];
if (~o_opt.flag_resume || ~isfile(o_opt.grid))
    if isempty(P)
        P = load(o_opt.parameters);
    end

    v_GABA = q_v_GABA * P.A / P.a;
    G = namm_setup_exploration_grid(P, v_GABA, n_omega_GABA, v_Glu, n_omega_Glu);

    mk_out_dir = @(K) strjoin(arrayfun(@int2str, K, 'UniformOutput', false), '-');
    G.names = G.loop_ind.Variables;
    G.names = arrayfun(@(a) mk_out_dir(G.names(a, :)), 1 : size(G.names, 1), ...
        'UniformOutput', false);

    save(o_opt.grid, '-struct', 'G', '-v7.3')
end


% - Output names
if (~o_opt.flag_resume || ~isfile(o_opt.names))
    if isempty(G)
        G = load(o_opt.grid);
    end

    P_e = cellfun(@(c) G.(c)(:), G.var(:), 'UniformOutput', false);
    N = table(G.names(:), P_e{:}, 'VariableNames', [{'name'}; G.var(:)]);

    writetable(N, o_opt.names, 'AutoFitWidth', true, 'PreserveFormat', true, ...
        'UseExcel', false, 'WriteMode', 'replacefile')
end
end


function main(c_opt, s_opt, o_opt)
sol_p = struct('ode_solver', s_opt.ode_solver, 'init_method', 'load', 'eta', 0);


% - Connectivity matrices
sys_p.Omega_Pyr = readmatrix(c_opt.neu_mat);
sys_p.Omega_Ast = readmatrix(c_opt.ast_mat);


% - Labels (assumed to be consistent with connectivity matrices)
in_file = fullfile(o_opt.simul_dir, 'labels.xls');
if (~o_opt.flag_resume || ~isfile(in_file))
    L = readtable(c_opt.hier_file, 'VariableNamingRule', 'preserve');
    L = L(:, {c_opt.roi_name, c_opt.net_name});
    writetable(L, in_file)
else
    L = readtable(in_file, 'VariableNamingRule', 'preserve');
end
L = L.(c_opt.roi_name);


% - Default parameters
sys_p = default_param(sys_p, 'global');
sys_p.nodes = duplicate_object(default_param(struct(), 'nodal'), numel(L));


% - Investigations
sys_p = add_to_sys_param(sys_p, load(o_opt.parameters), []);

G = load(o_opt.grid);
P_e = cell2struct(cellfun(@(c) G.(c)(:), G.var, 'UniformOutput', false), ...
    G.var, 1);
N = getfield(readtable(o_opt.names, 'VariableNamingRule', 'preserve'), 'name');


% - Preliminary simulations to obtain initial states
fprintf('>>> Calibration simulations...\n')
init_equilibrium(G, s_opt.I_init, s_opt.h, s_opt.fs_save, sys_p, sol_p, P_e, ...
    N, L, o_opt.flag_resume, s_opt.name, o_opt.simul_dir);


% - Phase trajectories estimation
fprintf('>>> Main simulations...\n')
for k_r = 1 : s_opt.n_rep
    run_simul(s_opt.I, s_opt.h, s_opt.fs_save, sys_p, sol_p, P_e, N, L, ...
        o_opt.simul_dir, true, o_opt.flag_resume, s_opt.name, k_r, o_opt.simul_dir)
end

fprintf('\n')
end


function init_equilibrium(G, I, h, fs_save, sys_p, sol_p, P_e, N, L, ...
    flag_resume, simul_name, o_dir)
sol_p.eta = get_init_state(G, numel(N), sys_p, sol_p);
run_simul(I, h, fs_save, sys_p, sol_p, P_e, N, L, [], false, flag_resume, ...
    simul_name, 0, o_dir)
end


function eta = get_init_state(G, n_simul, sys_p, sol_p)
% - Initial state (nodal)
sys_p_nodal = sys_p;
sys_p_nodal.Omega_Pyr = zeros(n_simul, n_simul);
sys_p_nodal.Omega_Ast = zeros(n_simul, n_simul);
sys_p_nodal.nodes = sys_p_nodal.nodes(ones(n_simul, 1));

[~, ~, eta_nodal, ~, ode_nodal] = nnamm_phase_trajectory([0, 0], 1, ...
    sys_p_nodal, sol_p, true);
ode_nodal = ode_nodal.ode;

eta_nodal = zeros(size(eta_nodal));
eta_nodal(ode_nodal(:, 1)) = G.E(:, 1);
eta_nodal(ode_nodal(:, 2)) = G.E(:, 2);
eta_nodal(ode_nodal(:, 3)) = G.E(:, 3);
eta_nodal(ode_nodal(:, 7)) = G.Glu_upt_r;
eta_nodal(ode_nodal(:, 9)) = G.Glu_e;
eta_nodal(ode_nodal(:, 10)) = G.Glu_Ast;
eta_nodal(ode_nodal(:, 11)) = G.GABA_upt_r;
eta_nodal(ode_nodal(:, 13)) = G.GABA_e;
eta_nodal(ode_nodal(:, 14)) = G.GABA_Ast;


% - Initial state (network)
[~, ~, eta] = nnamm_phase_trajectory([0, 0], 1, sys_p, sol_p, true);
n_nodes = numel(sys_p.nodes);
eta = repmat({eta}, n_simul, 1);
for k_s = 1 : n_simul
    eta{k_s}(:) = flatten(repmat(eta_nodal(ode_nodal(k_s, :)), n_nodes, 1));
end
end


function run_simul(I, h, fs, sys_param, solver_param, explore_param, N, L, ...
    init_dir_root, flag_compress_outputs, flag_resume, simul_name, rep, o_dir)
ode = 1 : 14;
if flag_compress_outputs
    % Keeps: {E_Pyr, E_ExIn_and_Pyr, E_InIn, J_Glu, Glu_e, Glu_Ast, J_GABA, GABA_e,
    % GABA_Ast}
    ode = [1, 2, 3, 7, 9, 10, 11, 13, 14];
end

if (~isempty(init_dir_root) && (rep == 0))
    error('Repetition 0 is reserved for calibration simulations...')
end

name = sprintf('%s-r%d', simul_name, rep);
fprintf('>>> > Processing [%s] ', name)


% - Preliminary outputs
out_dir_root = fullfile(o_dir, name);
create_dir(out_dir_root, false, false)

out_file = fullfile(out_dir_root, sprintf('%s.asv', mfilename));
copyfile(sprintf('%s.m', mfilename('fullpath')), out_file);

out_file = fullfile(out_dir_root, 'labels.orig.txt');
writecell(L, out_file)


% - Simulations
loop_ind_e_p = ones(1, numel(fieldnames(explore_param)));
n_simul = numel(N);
init_dir = fullfile(init_dir_root, [simul_name, '-r0']);
t_start = tic;
n_iter = 0;
n_bytes = 0;
for k_s = 1 : n_simul
    n = N{k_s};
    update_progress(n_simul)


    % - Outputs
    out_dir = fullfile(out_dir_root, n, 'activity');
    create_dir(out_dir, false, false)

    out_file = fullfile(out_dir, 'simul.mat');

    if (flag_resume && isfile(out_file))
        continue
    end


    % - System parameters
    sys_p = add_to_sys_param(sys_param, explore_param, k_s * loop_ind_e_p);
    sol_p = solver_param;
    if ~isempty(init_dir_root)
        sol_p.eta = load_init_state(init_dir, n);
    end
    if iscell(sol_p.eta)
        sol_p.eta = sol_p.eta{k_s};
    end


    % - Phase trajectory
    [Y, t, ~, eta_final, param] = nnamm_phase_trajectory(I([1, 2]), h, sys_p, ...
        sol_p, false);


    % - Saving
    mask = find((t >= sum(I([1, 3]))) & (t <= I(2)));

    Y = Y(mask, :);
    t = t(mask) - t(mask(1));
    if ~isempty(fs)
        [Y, t] = resample(Y, t, fs, Dimension = 1);
    end
    param.q.mu = sparse(param.q.mu(mask, :));
    param.q.sigma = sparse(param.q.sigma(mask, :));

    S = struct('eta_final', eta_final, 't', t, 'labels', {L}, 'ode', ode);
    S.q = param.q;
    S.param = rmfield(param, {'q'});
    for k_ode = 1 : numel(ode)
        S.(sprintf('Y_%d', ode(k_ode))) = Y(:, param.ode(:, ode(k_ode)));
    end
    save(out_file, '-struct', 'S', '-v7.3')


    % - Cleaning
    Y = []; %#ok
    t = []; %#ok
    eta_final = []; %#ok
    param = []; %#ok
    S = []; %#ok
end
t_end = toc(t_start);
fprintf(' Done in %.2f seconds.\n', t_end)


% - Local functions
    function update_progress(n)
        n_iter = n_iter + 1;
        fprintf(repmat('\b', 1, n_bytes))
        n_bytes = fprintf('(%d%%)...', round(n_iter / n * 100));
    end
end


function eta = load_init_state(in_dir_root, simul_name)
% Assumes that all ODEs are present (i.e., no output compression)
in_file = fullfile(in_dir_root, simul_name, 'activity', 'simul.mat');
ode = 1 : size(getfield(load_from_mat(in_file, 'param'), 'ode'), 2);
Y = arrayfun(@(a) load_from_mat(in_file, sprintf('Y_%d', a)), ode, ...
    'UniformOutput', false);
Y = cat(2, Y{:});
eta = Y(randi(size(Y, 1)), :);
end


function P = default_param(P, type)
if strcmp(type, 'global')
    P.omega_Pyr = 1;
    P.omega_Glu = 1;
    P.omega_GABA = 1;

elseif strcmp(type, 'nodal')
    % Neurons
    P.q = struct('method', 'gaussian', 'mean', 12, 'std', 1);

    P.A = 3.25;
    P.B = 22;
    P.a = 100;
    P.b = 50;
    P.nu_max = 5;
    P.r = 0.56;
    P.v_0 = 6;
    P.C_Pyr_to_ExIn = 150;
    P.C_ExIn_to_Pyr = 120;
    P.C_Pyr_to_InIn = 37.5;
    P.C_InIn_to_Pyr = 37.5;
    P.C_Pyr_to_Pyr = 40;

    % Glutamate
    P.W_Ast = 53.6;
    P.W_Pyr = 53.6;
    P.w_r = 90;
    P.w_d = 33;

    P.V_Glu_e_to_Ast = 4.5;
    P.r_Glu_e_to_Ast = 0.5;
    P.theta_Glu_e_to_Ast = 9;
    P.delta_Glu_e_to_Ast = 0;

    P.V_Glu_e_to_Pyr = 0.5;
    P.r_Glu_e_to_Pyr = 0.5;
    P.theta_Glu_e_to_Pyr = 9;
    P.delta_Glu_e_to_Pyr = 0;

    P.tau_Glu_Ast = 1 / 9;

    P.m_Glu_Pyr = 2.5;
    P.r_Glu_Pyr = 0.15;
    P.theta_Glu_Pyr = 30;
    P.delta_Glu_Pyr = 0;

    P.m_Glu_InIn = 1;
    P.r_Glu_InIn = 0.15;
    P.theta_Glu_InIn = 30;
    P.delta_Glu_InIn = 0;

    P.m_Glu_Ast = 4.5;
    P.r_Glu_Ast = 0.5;
    P.theta_Glu_Ast = 9;
    P.delta_Glu_Ast = 0;

    % GABA
    P.Z_Ast = 53.6;
    P.Z_InIn = 53.6;
    P.z_r = 90;
    P.z_d = 33;

    P.V_GABA_e_to_Ast = 2;
    P.K_GABA_e_to_Ast = 8;

    P.V_GABA_e_to_InIn = 5;
    P.K_GABA_e_to_InIn = 24;

    P.tau_GABA_Ast = 1 / 9;

    P.m_GABA_Pyr = 1;
    P.r_GABA_Pyr = 0.12;
    P.theta_GABA_Pyr = 25;
    P.delta_GABA_Pyr = 0;

end
end


function S = add_to_sys_param(S, A, loop_ind)
n = fieldnames(A);
if isempty(loop_ind)
    loop_ind = ones(1, numel(n));
end
for k_n = 1 : numel(n)
    for k_ind = 1 : size(loop_ind, 1)
        switch n{k_n}
            case {
                    'omega_Pyr'; 'omega_Glu'; 'omega_GABA';
                    }
                S.(n{k_n}) = A.(n{k_n})(loop_ind(k_ind, k_n));
            case {
                    'q'; 'm_Glu_Ast'; 'r_Glu_Ast'; 'theta_Glu_Ast'; 'delta_Glu_Ast';
                    'A'; 'B'; 'a'; 'b'; 'nu_max'; 'r'; 'v_0'; 'C_Pyr_to_ExIn'; 'C_ExIn_to_Pyr';
                    'C_Pyr_to_InIn'; 'C_InIn_to_Pyr'; 'C_Pyr_to_Pyr';
                    'W_Ast'; 'W_Pyr'; 'w_r'; 'w_d'; 'V_Glu_e_to_Ast'; 'r_Glu_e_to_Ast';
                    'theta_Glu_e_to_Ast'; 'delta_Glu_e_to_Ast'; 'V_Glu_e_to_Pyr'; 'r_Glu_e_to_Pyr';
                    'theta_Glu_e_to_Pyr'; 'delta_Glu_e_to_Pyr'; 'tau_Glu_Ast'; 'm_Glu_Pyr';
                    'r_Glu_Pyr'; 'theta_Glu_Pyr'; 'delta_Glu_Pyr'; 'm_Glu_InIn'; 'r_Glu_InIn';
                    'theta_Glu_InIn'; 'delta_Glu_InIn';
                    'Z_Ast'; 'Z_InIn'; 'z_r'; 'z_d'; 'V_GABA_e_to_Ast'; 'K_GABA_e_to_Ast';
                    'V_GABA_e_to_InIn'; 'K_GABA_e_to_InIn'; 'tau_GABA_Ast'; 'm_GABA_Pyr';
                    'r_GABA_Pyr'; 'theta_GABA_Pyr'; 'delta_GABA_Pyr';
                    }
                [S.nodes.(n{k_n})] = deal(A.(n{k_n})(loop_ind(k_ind, k_n)));
        end
    end
end
end
