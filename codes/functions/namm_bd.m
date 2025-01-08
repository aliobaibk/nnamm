function BD = namm_bd(d_param, e_param, flag_prod, h_E_Pyr, ...
    flag_disp_summary, out_file)
% Simplified codimension-1 bifurcation diagram with 'q' as the parameter


% - Default parameters
BD.default_param = default_param();
if ~isempty(d_param)
    discard_fields = setdiff(fieldnames(d_param), fieldnames(BD.default_param));
    BD.default_param = rmfield(d_param, discard_fields);
end


% - Parameters to explore
[P, BD.explored_param, BD.loop_ind] = exploration_param(BD.default_param, ...
    e_param, flag_prod);


% - Draws and classifies diagram
n = max(1, size(BD.loop_ind, 1));
sp = cell(n, 1);
for k = 1 : n
    sp{k} = singular_points(h_E_Pyr, P(k));
end
BD.P = sp;

if flag_disp_summary
    display_summary(BD.P);
end

BD.h_E_Pyr = h_E_Pyr;
BD.n_ode = 6;


% - Saves
if ~isempty(out_file)
    create_dir(fileparts(out_file), false, false)
    save(out_file, '-struct', 'BD', '-v7.3')
end
end


%% ------------ Local functions ------------------------------------------------
function P = default_param()
% Jansen–Rit parameterization
P.A = 3.25;
P.B = 22;
P.a = 100;
P.b = 50;
P.nu_max = 5;
P.r = 0.56;
P.v_0 = 6;
P.v_GABA = 0;
P.v_Glu = 0;
P.mu_Glu_InIn_by_Pyr = 0.5;
P.C_Pyr_to_ExIn = 135;
P.C_ExIn_to_Pyr = 108;
P.C_Pyr_to_InIn = 33.75;
P.C_InIn_to_Pyr = 33.75;
P.C_Pyr_to_Pyr = 0;
end


function SP = singular_points(h_E_Pyr, P)
SP.E_Pyr = (h_E_Pyr : h_E_Pyr : (P.A * P.nu_max / P.a - h_E_Pyr)).';

v = exp(P.r * [P.v_0, P.v_0 - P.mu_Glu_InIn_by_Pyr * P.v_Glu] - (SP.E_Pyr * ...
    (P.r * [P.C_Pyr_to_ExIn, P.C_Pyr_to_InIn])));
S = ([P.A / P.a * P.C_ExIn_to_Pyr, P.B / P.b * P.C_InIn_to_Pyr] * ...
    P.nu_max) ./ (1 + v);
SP.E_ExIn_and_Pyr = S(:, 1) + P.C_Pyr_to_Pyr * SP.E_Pyr;
SP.E_InIn = S(:, 2);

SP.q = (P.a / P.A) * (P.v_0 + P.v_GABA - P.v_Glu - (log((P.A * P.nu_max / ...
    P.a) ./ SP.E_Pyr - 1) / P.r + SP.E_ExIn_and_Pyr - SP.E_InIn));
SP.E_ExIn_and_Pyr = SP.E_ExIn_and_Pyr + (P.A / P.a) * SP.q;

SP.LFP = SP.E_ExIn_and_Pyr - SP.E_InIn;

J = P.r * (S .* v ./ (1 + v));
J = jacobian(P, J, SP.LFP);

SP.eigen_values = arrayfun(@(k) eig(J(:, :, k)), 1 : size(SP.E_Pyr, 1), ...
    'UniformOutput', false);
SP.eigen_values = cat(2, SP.eigen_values{:}).';
SP.n_neg_real_part = sum(real(SP.eigen_values) < 0, 2);
SP.bif_label = bd_label(SP.n_neg_real_part, 6);

SP.param = P;
end


function J = jacobian(P, d_S, LFP)
v = exp(P.r * (P.v_0 + P.v_GABA - P.v_Glu - LFP));
d = (P.A * P.a * P.nu_max * P.r) * (v ./ ((1 + v) .^ 2));

J = repmat([ ...
    0, 0, 0, 1, 0, 0;
    0, 0, 0, 0, 1, 0;
    0, 0, 0, 0, 0, 1;
    -P.a * P.a, 0, 0, -2 * P.a, 0, 0;
    0, 0, 0, 0, -2 * P.a, 0;
    0, 0, -P.b * P.b, 0, 0, -2 * P.b; ...
    ], [1, 1, size(d, 1)]);

J(4, 2, :) = d;
J(4, 3, :) = -d;
J(5, 1, :) = (P.a * P.a * P.C_Pyr_to_ExIn) * d_S(:, 1);
J(5, 2, :) = P.C_Pyr_to_Pyr * d - P.a * P.a;
J(5, 3, :) = (-P.C_Pyr_to_Pyr) * d;
J(6, 1, :) = (P.b * P.b * P.C_Pyr_to_InIn) * d_S(:, 2);
end


function [P, e_param, loop_ind] = exploration_param(d_param, a_param, flag_prod)
if isempty(a_param)
    P = d_param;
    e_param = [];
    loop_ind = table();
else
    e_param = a_param(:, 1);
    N = cellfun(@numel, a_param(:, 2), 'UniformOutput', false);
    if flag_prod
        loop_ind = cartesian_product(N{:});
    else
        loop_ind = cell2mat(cellfun(@(c) 1 : c, N(:), 'UniformOutput', false)).';
    end
    P = duplicate_object(d_param, size(loop_ind, 1));
    for k = 1 : size(a_param, 1)
        [P.(a_param{k, 1})] = cell_expand(num2cell(a_param{k, 2}(loop_ind(:, k))));
    end
    loop_ind = array2table(loop_ind, 'VariableNames', e_param);
end
end


function b = bd_label(a, n)
d = diff(a);
[k, ~, d] = find(d);
d = abs(d);

if isempty(d)
    if (a(1) == n)
        b.name = 'all-stable';
        b.n_hopfs = 0;
        b.n_saddles = 0;
    elseif (a(1) ~= n)
        b.name = 'all-unstable';
        b.n_hopfs = 0;
        b.n_saddles = 0;
    end
    b.k_sing = {'', []};
else
    hopfs = find(d == 2);
    saddles = find(d == 1);
    b.name = char(zeros(1, numel(d)));
    b.name(hopfs) = 'h';
    b.name(saddles) = 's';
    b.n_hopfs = numel(hopfs);
    b.n_saddles = numel(saddles);
    b.k_sing = {};
    for kk = 1 : numel(d)
        if (d(kk) == 2)
            b.k_sing{end + 1, 1} = 'hopf';
            b.k_sing{end, 2} = k(kk) + 1;
        elseif (d(kk) == 1)
            b.k_sing{end + 1, 1} = 'saddle';
            b.k_sing{end, 2} = k(kk) + 1;
        end
    end
end
end


function display_summary(P)
names = cellfun(@(c) c.bif_label.name, P, 'UniformOutput', false);
cellfun(@(c) fprintf('%s: %d%%\n', c, round(sum(cellfun(@(cc) strcmp(cc, c), ...
    names, 'UniformOutput', true)) / numel(names) * 100)), unique(names));
end
