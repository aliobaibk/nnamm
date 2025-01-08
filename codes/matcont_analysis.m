function matcont_analysis()
% - This function uses MatCont in order to continue saddle-node and Hopf
%   branches with (q, v_Glu) as parameters.
% - The input parameters and initial states are obtained from a pre-computed
%   bifurcation diagram (e.g., after running 'bd_draw_alpha.m').
% - This is not a fully automatic continuation. Knowledge from the pre-computed
%   diagram is needed.

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


% - Sets up MatCont
matcont_bin = './bin/lib-misc/MatCont7p3';

current_dir = pwd;
cd(matcont_bin)
init
cd(current_dir)


% - Pseudo codimension-2 diagrams (q, v_Glu, C_Pyr_to_Pyr)
% * By default, outputs are saved in the same directory as the input files,
%   using names consistent with input file names.
in_dir = fullfile(proj_dir, 'data', 'bifurcation-diagrams');
in_files = [
    struct('dir', fullfile(in_dir, 'alpha'), 'name', 'C_Pyr_to_Pyr-0');
    struct('dir', fullfile(in_dir, 'alpha'), 'name', 'C_Pyr_to_Pyr-7.5');
    ];


% - Continuations
for k_f = 1 : numel(in_files)
    in_dir = in_files(k_f).dir;
    bd_name = in_files(k_f).name;
    in_file = fullfile(in_dir, sprintf('%s.bd.mat', bd_name));
    out_file = fullfile(in_dir, bd_name);
    bd = load_from_mat(in_file, 'P');
    saddle_node_analysis(bd, out_file);
    hopf_analysis(bd, out_file);
end
end


%% ------------ Local functions ------------------------------------------------
function Results = saddle_node_analysis(bd, out_file)
v_Glu_min = -7;
v_Glu_max = 7;


% - Loads inputs
% Parameter names
par_names = namm_matcont;
par_names = par_names{end};

% Starts in the vicinity of a saddle-node bifurcation point in a 'JR-diagram'
K = find(strcmp(cellfun(@(c) c.bif_label.name, bd(:), ...
    'UniformOutput', false), 'sshhh'));
bd = bd{K(round(end / 2))};

K = bd.bif_label.k_sing(strcmp(bd.bif_label.k_sing(:, 1), 'saddle'), 2);
K = K{1} - 3;

x0_0 = [
    bd.E_Pyr(K);
    bd.E_ExIn_and_Pyr(K);
    bd.E_InIn(K);
    0;
    0;
    0;
    ];

bd.param.q = bd.q(K);
P0 = cellfun(@(c) bd.param.(c), par_names(:), 'UniformOutput', true);


% - Equilibrium curve continuation (forward)
% The goal here is to find the closest saddle-node bifurcation point

Results = [];
R = struct('proc', '', 'X', [], 'V', [], 'S', [], 'H', [], 'F', []);

global cds %#ok

ap_1 = find(strcmp(par_names, 'q'), 1);
[x0, v0] = init_EP_EP(@namm_matcont, x0_0, P0, ap_1);

opt = contset();
opt = contset(opt, 'Eigenvalues', 1);
opt = contset(opt, 'MaxNumPoints', 50);
opt = contset(opt, 'InitStepsize', 10^(-2));
opt = contset(opt, 'MaxStepsize', 10^(-1));
opt = contset(opt, 'MinStepsize', 10^(-3));
opt = contset(opt, 'Singularities', 1);
opt = contset(opt, 'SymDerivative', 3);
opt = contset(opt, 'SymDerivativeP', 2);
opt = matcont_set_cont_dir(@equilibrium, x0, v0, opt, 1);

R.proc = 'ep';
[R.X, R.V, R.S, R.H, R.F] = cont(@equilibrium, x0, v0, opt);
Results = [Results; R];


% - Saddle-node curve continuation
% The goal here is to continue a saddle-node branch for increasing 'v_Glu' until
% finding a BT bifurcation, then a CP bifurcation, then a BT bifurcation. From
% that point onwards, the saddle-node branch can be continued backward for
% decreasing 'v_Glu' until the specified bound.

ap_2 = [ap_1, find(strcmp(par_names, 'v_Glu'), 1)];

opt = contset(opt, 'MaxNumPoints', 1500);

k_res = 1; % uses results from equilibrium curve continuation
R_ = Results(k_res);

k = find(strcmp({R_.S.label}, 'LP'), 1);
x0 = R_.X(1 : (end - 1), R_.S(k).index);
P = P0;
P(ap_2(1)) = R_.X(end, R_.S(k).index);

[x0, v0] = init_LP_LP(@namm_matcont, x0, P, ap_2);

[R.X, R.V, R.S, R.H, R.F] = cont(@limitpoint, x0, v0, opt);
while true
    if any(R.X(end, (end - opt.MaxNumPoints + 1) : end) < v_Glu_min)
        break
    end
    [R.X, R.V, R.S, R.H, R.F] = cont(R.X, R.V, R.S, R.H, R.F, cds);
end
R.proc = 'lp-f';
Results = [Results; R];

opt = contset(opt, 'Backward', 1);
[R.X, R.V, R.S, R.H, R.F] = cont(@limitpoint, x0, v0, opt);
while true
    if any(R.X(end, (end - opt.MaxNumPoints + 1) : end) < v_Glu_min)
        break
    end
    [R.X, R.V, R.S, R.H, R.F] = cont(R.X, R.V, R.S, R.H, R.F, cds);
end
R.proc = 'lp-b';
Results = [Results; R];


% - Hopf curve continuation from BT bifurcation
% The goal here is to continue Hopf branches from the first BT bifurcation point
% obtained earlier until the specified bounds.

ap_2 = [ap_1, find(strcmp(par_names, 'v_Glu'), 1)];

opt = contset(opt, 'MaxNumPoints', 10000);

BT = extract_bt(Results);
P = P0;
P(ap_2) = BT([end - 1, end], 1);

[x0, v0] = init_BT_H(@namm_matcont, BT(1 : (end - 2), 1), P, ap_2);
opt = contset(opt, 'Backward', 0);
[R.X, R.V, R.S, R.H, R.F] = cont(@hopf, x0, v0, opt);
while true
    if any(R.X(end - 1, (end - opt.MaxNumPoints + 1) : end) < v_Glu_min)
        break
    end
    [R.X, R.V, R.S, R.H, R.F] = cont(R.X, R.V, R.S, R.H, R.F, cds);
end
R.proc = 'bt-f';
Results = [Results; R];

[x0, v0] = init_BT_H(@namm_matcont, BT(1 : (end - 2), 1), P, ap_2);
opt = contset(opt, 'Backward', 1);
[R.X, R.V, R.S, R.H, R.F] = cont(@hopf, x0, v0, opt);
while true
    if any(R.X(end - 1, (end - opt.MaxNumPoints + 1) : end) > v_Glu_max)
        break
    end
    [R.X, R.V, R.S, R.H, R.F] = cont(R.X, R.V, R.S, R.H, R.F, cds);
end
R.proc = 'bt-b';
Results = [Results; R];


% - Saves curves
out_file = sprintf('%s.saddle-node.mat', out_file);
save(out_file, 'Results', '-v7.3')
end


function Results = hopf_analysis(bd, out_file)
v_Glu_min = -7;
v_Glu_max = 7;


% - Loads inputs
% Parameter names
par_names = namm_matcont;
par_names = par_names{end};

% Starts in the vicinity of a stable Hopf bifurcation point in a 'JR-diagram'
K = find(strcmp(cellfun(@(c) c.bif_label.name, bd(:), ...
    'UniformOutput', false), 'sshhh'));
bd = bd{K(round(end / 2))};

K = bd.bif_label.k_sing(strcmp(bd.bif_label.k_sing(:, 1), 'hopf'), 2);
K = K{end} - 3;

x0_0 = [
    bd.E_Pyr(K);
    bd.E_ExIn_and_Pyr(K);
    bd.E_InIn(K);
    0;
    0;
    0;
    ];

bd.param.q = bd.q(K);
P0 = cellfun(@(c) bd.param.(c), par_names(:), 'UniformOutput', true);


% - Equilibrium curve continuation (forward)
% The goal here is to find the closest Hopf bifurcation point

Results = [];
R = struct('proc', '', 'X', [], 'V', [], 'S', [], 'H', [], 'F', []);

global cds %#ok

ap_1 = find(strcmp(par_names, 'q'), 1);
[x0, v0] = init_EP_EP(@namm_matcont, x0_0, P0, ap_1);

opt = contset();
opt = contset(opt, 'Eigenvalues', 1);
opt = contset(opt, 'MaxNumPoints', 50);
opt = contset(opt, 'InitStepsize', 10^(-2));
opt = contset(opt, 'MaxStepsize', 10^(-1));
opt = contset(opt, 'MinStepsize', 10^(-3));
opt = contset(opt, 'Singularities', 1);
opt = contset(opt, 'SymDerivative', 3);
opt = contset(opt, 'SymDerivativeP', 2);
opt = matcont_set_cont_dir(@equilibrium, x0, v0, opt, 1);

R.proc = 'ep';
[R.X, R.V, R.S, R.H, R.F] = cont(@equilibrium, x0, v0, opt);
Results = [Results; R];


% - Hopf curve continuation
% The goal here is to continue the Hopf branch for increasing 'v_Glu' until the
% specified bound then continue the Hopf branch for decreasing 'v_Glu' until
% both a GH bifurcation and the specified bound.

ap_2 = [ap_1, find(strcmp(par_names, 'v_Glu'), 1)];

k_res = 1; % uses results from equilibrium curve continuation
R_ = Results(k_res);

k = find(strcmp({R_.S.label}, 'H '), 1);
x0 = R_.X(1 : (end - 1), R_.S(k).index);
P = P0;
P(ap_2(1)) = R_.X(end, R_.S(k).index);

[x0, v0] = init_H_H(@namm_matcont, x0, P, ap_2);

opt = contset(opt, 'MaxNumPoints', 1500);
[R.X, R.V, R.S, R.H, R.F] = cont(@hopf, x0, v0, opt);
while true
    if any(R.X(end - 1, (end - opt.MaxNumPoints + 1) : end) > v_Glu_max)
        break
    end
    [R.X, R.V, R.S, R.H, R.F] = cont(R.X, R.V, R.S, R.H, R.F, cds);
end
R.proc = 'hp-f';
Results = [Results; R];

opt = contset(opt, 'MaxNumPoints', 10000);
opt = contset(opt, 'Backward', 1);
[R.X, R.V, R.S, R.H, R.F] = cont(@hopf, x0, v0, opt);
while true
    if (ismember('GH', {R.S.label}) && any(R.X(end - 1, (end - ...
            opt.MaxNumPoints + 1) : end) < v_Glu_min))
        break
    end
    [R.X, R.V, R.S, R.H, R.F] = cont(R.X, R.V, R.S, R.H, R.F, cds);
end
R.proc = 'hp-b';
Results = [Results; R];


% - Saves curves
out_file = sprintf('%s.hopf.mat', out_file);
save(out_file, 'Results', '-v7.3')
end


function BT = extract_bt(SN)
BT = [];
for k_sn = 1 : numel(SN)
    K = find(strcmp({SN(k_sn).S.label}, 'BT'));
    BT = [BT, SN(k_sn).X(:, [SN(k_sn).S(K).index])]; %#ok
end
end
