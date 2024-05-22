function G = namm_setup_exploration_grid(d_param, v_GABA, n_omega_GABA, ...
    v_Glu, n_omega_Glu)
% * Using steady-state calculations, it is possible to derive a grid for
%   (omega_GABA, omega_Glu) based on a grid with (v_GABA, v_Glu).
% * Note that it is not always trivial to account for non-null network
%   contributions.
% * Using (v_GABA, v_Glu) at the onset is convenient because we are then able to
%   sample more or less uniformly the complex dynamical landscape of interest.
%   Although it is desirable to further enrich the sampling by using (GABA_e,
%   Glu_e) directly.

% - Input parameters for drawing bifurcation diagram
% * Neuronal network feedback can be accounted for by using the parameter
%   C_Pyr_to_Pyr and assuming a homogeneous parameterization of the mass models
%   making up the network.
if isfield(d_param, 'omega_Pyr')
    d_param.C_Pyr_to_Pyr = d_param.C_Pyr_to_Pyr + d_param.omega_Pyr;
end


% - Desired exploration grid from bifurcation parameters
v_GABA = linspace(v_GABA(1), v_GABA(2), n_omega_GABA);
v_Glu = linspace(v_Glu(1), v_Glu(2), n_omega_Glu);


% - Draws the bifurcation diagram for the pairs (v_GABA, v_Glu)
e_param = {
    'v_GABA', v_GABA(:);
    'v_Glu', v_Glu(:);
    };

G.BD = namm_bd(d_param, e_param, true, 10^(-4), false, '');

G.loop_ind = G.BD.loop_ind;
G.v_GABA = cellfun(@(d) d.param.v_GABA, G.BD.P(:), 'UniformOutput', true);
G.v_Glu = cellfun(@(d) d.param.v_Glu, G.BD.P(:), 'UniformOutput', true);


% - Extracts information from the bifurcation diagram
% NOTES:
% * It is assumed that the chosen bifurcation landscape (E_Pyr = f(v_GABA,
%   v_Glu)) is injective within the targeted domain. This is important as
%   singular points will be read from the pre-computed diagram. The code would
%   need revisions when a pair (v_GABA, v_Glu) can map to multiple values of
%   E_Pyr. For instance, in the case of non-injection, criteria based on a
%   reference state could help targetting automatically the singular points of
%   interest.
% * It is also assumed that a codimension-1 diagram with 'q' as parameter was
%   drawn for each pair (v_GABA, v_Glu).
q = cell(n_omega_GABA * n_omega_Glu, 1);
E = cell(n_omega_GABA * n_omega_Glu, 1);
for k_d = 1 : (n_omega_GABA * n_omega_Glu)
    bd = G.BD.P{k_d};
    [~, K] = min(abs(d_param.q.mean - bd.q));
    q{k_d} = bd.q(K);
    E{k_d} = [bd.E_Pyr(K), bd.E_ExIn_and_Pyr(K), bd.E_InIn(K)];
end
G.q = cat(1, q{:});
G.E = cat(1, E{:});


% - Inferred exploration grid for (omega_GABA, omega_Glu)
G.var = {'omega_GABA'; 'omega_Glu'}; % consistent with variable 'loop_ind'

F_Pyr = sigmoid(G.E(:, 2) - G.E(:, 3), d_param.nu_max, d_param.r, ...
    d_param.v_0 + G.v_GABA - G.v_Glu, 0);
F_InIn = sigmoid(d_param.C_Pyr_to_InIn * G.E(:, 1), d_param.nu_max, ...
    d_param.r, d_param.v_0 - d_param.mu_Glu_InIn_by_Pyr * G.v_Glu, 0);

G.Glu_e = sigmoid_inv(G.v_Glu, d_param.m_Glu_Pyr, d_param.r_Glu_Pyr, ...
    d_param.theta_Glu_Pyr, d_param.delta_Glu_Pyr);
G.GABA_e = sigmoid_inv(G.v_GABA, d_param.m_GABA_Pyr, d_param.r_GABA_Pyr, ...
    d_param.theta_GABA_Pyr, d_param.delta_GABA_Pyr);

Glu_upt_r_Ast = sigmoid(G.Glu_e, d_param.V_Glu_e_to_Ast, ...
    d_param.r_Glu_e_to_Ast, d_param.theta_Glu_e_to_Ast, ...
    d_param.delta_Glu_e_to_Ast);
GABA_upt_r_Ast = mich_menten_kin(G.GABA_e, d_param.V_GABA_e_to_Ast, ...
    d_param.K_GABA_e_to_Ast);

G.Glu_upt_r = Glu_upt_r_Ast + sigmoid(G.Glu_e, d_param.V_Glu_e_to_Pyr, ...
    d_param.r_Glu_e_to_Pyr, d_param.theta_Glu_e_to_Pyr, ...
    d_param.delta_Glu_e_to_Pyr);
G.GABA_upt_r = GABA_upt_r_Ast + mich_menten_kin(G.GABA_e, ...
    d_param.V_GABA_e_to_InIn, d_param.K_GABA_e_to_InIn);

G.Glu_Ast = Glu_upt_r_Ast * d_param.tau_Glu_Ast;
G.GABA_Ast = GABA_upt_r_Ast * d_param.tau_GABA_Ast;

Q_Ast = sigmoid(G.Glu_e, d_param.m_Glu_Ast, d_param.r_Glu_Ast, ...
    d_param.theta_Glu_Ast, d_param.delta_Glu_Ast);

G.omega_Glu = (d_param.w_d * G.Glu_upt_r - d_param.W_Pyr * F_Pyr) ./ ( ...
    d_param.W_Ast * Q_Ast);
G.omega_GABA = (d_param.z_d * G.GABA_upt_r - d_param.Z_InIn * F_InIn) ./ ...
    (d_param.Z_Ast * Q_Ast);
end
