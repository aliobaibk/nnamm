function B = namm_setup_baseline(d_param, v_GABA, omega_GABA, v_Glu, omega_Glu)
% * Using steady-state calculations, it is possible to define a baseline state.
% * The main goal of this function is to define values for W_* and Z_* given a
%   target neuronal state defined by (v_GABA, v_Glu).
% * Note that it is not always trivial to account for non-null network
%   contributions.
% * Using (v_GABA, v_Glu) at the onset is convenient to directly map a state on
%   the bifurcation diagram.
% * It is assumed that neurons and astrocytes contribute equally likely to
%   releases (i.e., W_Ast = W_Pyr = W and Z_Ast = Z_InIn = Z). This assumption
%   could be relaxed by assuming a certain ratio between W_Ast and W_Pyr or
%   Z_Ast and Z_InIn.

% - Input parameters for drawing bifurcation diagram
% * Neuronal network feedback can be accounted for by using the parameter
%   C_Pyr_to_Pyr and assuming a homogeneous parameterization of the mass models
%   making up the network.
if isfield(d_param, 'omega_Pyr')
    d_param.C_Pyr_to_Pyr = d_param.C_Pyr_to_Pyr + d_param.omega_Pyr;
end


% - Draws the bifurcation diagram for the pair (v_GABA, v_Glu)
e_param = {
    'v_GABA', v_GABA;
    'v_Glu', v_Glu;
    };

B.BD = namm_bd(d_param, e_param, true, 10^(-4), false, '');

B.v_GABA = v_GABA;
B.omega_GABA = omega_GABA;
B.v_Glu = v_Glu;
B.omega_Glu = omega_Glu;


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
[~, K] = min(abs(d_param.q.mean - B.BD.P{1}.q));
B.q = B.BD.P{1}.q(K);
B.E = [B.BD.P{1}.E_Pyr(K), B.BD.P{1}.E_ExIn_and_Pyr(K), B.BD.P{1}.E_InIn(K)];


% - Inferred parameters for targetting the specified baseline state
F_Pyr = sigmoid(B.E(2) - B.E(3), d_param.nu_max, d_param.r, d_param.v_0 + ...
    B.v_GABA - B.v_Glu, 0);
F_InIn = sigmoid(d_param.C_Pyr_to_InIn * B.E(1), d_param.nu_max, ...
    d_param.r, d_param.v_0 - d_param.mu_Glu_InIn_by_Pyr * B.v_Glu, 0);

B.Glu_e = sigmoid_inv(B.v_Glu, d_param.m_Glu_Pyr, d_param.r_Glu_Pyr, ...
    d_param.theta_Glu_Pyr, d_param.delta_Glu_Pyr);
B.GABA_e = sigmoid_inv(B.v_GABA, d_param.m_GABA_Pyr, d_param.r_GABA_Pyr, ...
    d_param.theta_GABA_Pyr, d_param.delta_GABA_Pyr);

Glu_upt_r_Ast = sigmoid(B.Glu_e, d_param.V_Glu_e_to_Ast, ...
    d_param.r_Glu_e_to_Ast, d_param.theta_Glu_e_to_Ast, ...
    d_param.delta_Glu_e_to_Ast);
GABA_upt_r_Ast = mich_menten_kin(B.GABA_e, d_param.V_GABA_e_to_Ast, ...
    d_param.K_GABA_e_to_Ast);

B.Glu_upt_r = Glu_upt_r_Ast + sigmoid(B.Glu_e, d_param.V_Glu_e_to_Pyr, ...
    d_param.r_Glu_e_to_Pyr, d_param.theta_Glu_e_to_Pyr, ...
    d_param.delta_Glu_e_to_Pyr);
B.GABA_upt_r = GABA_upt_r_Ast + mich_menten_kin(B.GABA_e, ...
    d_param.V_GABA_e_to_InIn, d_param.K_GABA_e_to_InIn);

B.Glu_Ast = Glu_upt_r_Ast * d_param.tau_Glu_Ast;
B.GABA_Ast = GABA_upt_r_Ast * d_param.tau_GABA_Ast;

Q_Ast = sigmoid(B.Glu_e, d_param.m_Glu_Ast, d_param.r_Glu_Ast, ...
    d_param.theta_Glu_Ast, d_param.delta_Glu_Ast);

B.W = B.Glu_upt_r * d_param.w_d / (F_Pyr + B.omega_Glu * Q_Ast);
B.Z = B.GABA_upt_r * d_param.z_d / (F_InIn + B.omega_GABA * Q_Ast);
end
