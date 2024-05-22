function out = namm_matcont()
out{1} = @init;
out{2} = @fun_eval;
out{3} = @jacobian;
out{4} = @jacobianp;
out{5} = @hessians;
out{6} = @hessiansp;
out{7} = @der3;
out{8} = [];
out{9} = [];
out{10} = {
    'A';
    'a';
    'B';
    'b';
    'C_Pyr_to_Pyr';
    'C_Pyr_to_ExIn';
    'C_ExIn_to_Pyr';
    'C_Pyr_to_InIn';
    'C_InIn_to_Pyr';
    'nu_max';
    'r';
    'v_0';
    'v_Glu';
    'v_GABA';
    'mu_Glu_InIn_by_Pyr';
    'q';
    };
end

% ------------------------------------------------------------------------------
function dYdt = fun_eval(t, Y, A, a, B, b, C_Pyr_to_Pyr, C_Pyr_to_ExIn, ...
    C_ExIn_to_Pyr, C_Pyr_to_InIn, C_InIn_to_Pyr, nu_max, r, v_0, v_Glu, v_GABA, ...
    mu_Glu_InIn_by_Pyr, q) %#ok
S = (nu_max * [1, C_ExIn_to_Pyr, C_InIn_to_Pyr]) ./ (1 + exp(r * [v_0 - ...
    v_Glu + v_GABA - Y(2) + Y(3), v_0 - C_Pyr_to_ExIn * Y(1), v_0 - ...
    mu_Glu_InIn_by_Pyr * v_Glu - C_Pyr_to_InIn * Y(1)]));
dYdt = [ ...
    ... % E_Pyr
    Y(4); ...
    ... % E_ExIn_and_Pyr
    Y(5); ...
    ... % E_InIn
    Y(6); ...
    ... % d[E_Pyr]
    a * (A * S(1) - 2 * Y(4) - a * Y(1)); ...
    ... % d[E_ExIn_and_Pyr]
    a * (A * (q + C_Pyr_to_Pyr * S(1) + S(2)) - 2 * Y(5) - a * Y(2)); ...
    ... % d[E_InIn]
    b * (B * S(3) - 2 * Y(6) - b * Y(3)); ...
    ];
end

% ------------------------------------------------------------------------------
function [tspan, y0, options] = init
handles = feval(namm_matcont);
y0 = [0, 0, 0, 0, 0, 0];
options = odeset('Jacobian', handles(3), 'JacobianP', handles(4), ...
    'Hessians', handles(5), 'HessiansP', handles(6));
tspan = [0, 10];
end

% ------------------------------------------------------------------------------
function jac = jacobian(t, Y, A, a, B, b, C_Pyr_to_Pyr, C_Pyr_to_ExIn, ...
    C_ExIn_to_Pyr, C_Pyr_to_InIn, C_InIn_to_Pyr, nu_max, r, v_0, v_Glu, v_GABA, ...
    mu_Glu_InIn_by_Pyr, q) %#ok
S = exp(r * ([v_0 - v_Glu + v_GABA - Y(2) + Y(3), v_0 - C_Pyr_to_ExIn * Y( ...
    1), v_0 - mu_Glu_InIn_by_Pyr * v_Glu - C_Pyr_to_InIn * Y(1)]));
S = (nu_max * r) * ([A * a, A * a * C_Pyr_to_ExIn * C_ExIn_to_Pyr, B * b * ...
    C_Pyr_to_InIn * C_InIn_to_Pyr] .* S ./ ((S + 1).^2));
jac = [ ...
    0, 0, 0, 1, 0, 0; ...
    0, 0, 0, 0, 1, 0; ...
    0, 0, 0, 0, 0, 1; ...
    -a * a, S(1), -S(1), -2 * a, 0, 0; ...
    S(2), C_Pyr_to_Pyr * S(1) - a * a, -C_Pyr_to_Pyr * S(1), 0, -2 * a, 0; ...
    S(3), 0, -b * b, 0, 0, -2 * b; ...
    ];
end

% ------------------------------------------------------------------------------
function jacp = jacobianp(t, Y, A, a, B, b, C_Pyr_to_Pyr, C_Pyr_to_ExIn, ...
    C_ExIn_to_Pyr, C_Pyr_to_InIn, C_InIn_to_Pyr, nu_max, r, v_0, v_Glu, v_GABA, ...
    mu_Glu_InIn_by_Pyr, q) %#ok
v = exp(r * [v_0 - v_Glu + v_GABA - Y(2) + Y(3), v_0 - C_Pyr_to_ExIn * ...
    Y(1), v_0 - mu_Glu_InIn_by_Pyr * v_Glu - C_Pyr_to_InIn * Y(1)]);
s = 1 ./ (v + 1);
S = v .* s .* s;
jacp = [ ...
    ... % d[E_Pyr]
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; ...
    ... % d[E_ExIn_and_Pyr]
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; ...
    ... % d[E_InIn]
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; ...
    ... % d[d[E_Pyr]]
    nu_max * a * s(1), ... % A
    nu_max * A * s(1) - 2 * Y(4) - 2 * a * Y(1), ... % a
    0, 0, 0, 0, 0, 0, 0, ...
    A * a * s(1), ... % nu_max
    -nu_max * A * a * (v_0 - v_Glu + v_GABA - Y(2) + Y(3)) * S(1), ... % r
    -nu_max * r * A * a * S(1), ... % v_0
    nu_max * r * A * a * S(1), ... % v_Glu
    -nu_max * r * A * a * S(1), ... % v_GABA
    0, 0; ...
    ... % d[d[E_ExIn_and_Pyr]]
    nu_max * a * (C_ExIn_to_Pyr * s(2) + C_Pyr_to_Pyr * s(1)) + a * q, ... % A
    nu_max * A * (C_ExIn_to_Pyr * s(2) + C_Pyr_to_Pyr * s(1)) - 2 * Y( ...
    5) - 2 * a * Y(2) + A * q, ... % a
    0, 0, ...
    nu_max * A * a * s(1), ... % C_Pyr_to_Pyr
    nu_max * r * A * a * C_ExIn_to_Pyr * Y(1) * S(2), ... % C_Pyr_to_ExIn
    nu_max * A * a * s(2), ... % C_ExIn_to_Pyr
    0, 0, ...
    A * a * (C_ExIn_to_Pyr * s(2) + C_Pyr_to_Pyr * s(1)), ... % nu_max
    -nu_max * A * a * (C_ExIn_to_Pyr * (v_0 - C_Pyr_to_ExIn * Y(1)) * S(2) + ...
    C_Pyr_to_Pyr * (v_0 - v_Glu + v_GABA - Y(2) + Y(3)) * S(1)), ... % r
    -nu_max * r * A * a * (C_ExIn_to_Pyr * S(2) + C_Pyr_to_Pyr * S(1)), ... % v_0
    nu_max * r * A * a * C_Pyr_to_Pyr * S(1), ... % v_Glu
    -nu_max * r * A * a * C_Pyr_to_Pyr * S(1), ... % v_GABA
    0, ...
    A * a; ... % q
    ... % d[d[E_InIn]]
    0, 0, ...
    nu_max * b * C_InIn_to_Pyr * s(3), ... % B
    nu_max * B * C_InIn_to_Pyr * s(3) - 2 * Y(6) - 2 * b * Y(3), ... % b
    0, 0, 0, ...
    nu_max * r * B * b * C_InIn_to_Pyr * Y(1) * S(3), ... % C_Pyr_to_InIn
    nu_max * B * b * s(3), ... % C_InIn_to_Pyr
    B * b * C_InIn_to_Pyr * s(3), ... % nu_max
    -nu_max * B * b * C_InIn_to_Pyr * (v_0 - mu_Glu_InIn_by_Pyr * v_Glu - ...
    C_Pyr_to_InIn * Y(1)) * S(3), ... % r
    -nu_max * r * B * b * C_InIn_to_Pyr * S(3), ... % v_0
    nu_max * r * B * b * C_InIn_to_Pyr * mu_Glu_InIn_by_Pyr * S(3), ... % v_Glu
    0, ...
    nu_max * r * B * b * C_InIn_to_Pyr * v_Glu * S(3), ... % mu_Glu_InIn_by_Pyr
    0; ...
    ];
end

% ------------------------------------------------------------------------------
function hess = hessians(t, Y, A, a, B, b, C_Pyr_to_Pyr, C_Pyr_to_ExIn, ...
    C_ExIn_to_Pyr, C_Pyr_to_InIn, C_InIn_to_Pyr, nu_max, r, v_0, v_Glu, v_GABA, ...
    mu_Glu_InIn_by_Pyr, q) %#ok
S = exp(r * ([v_0 - v_Glu + v_GABA - Y(2) + Y(3), v_0 - C_Pyr_to_ExIn * Y( ...
    1), v_0 - mu_Glu_InIn_by_Pyr * v_Glu - C_Pyr_to_InIn * Y(1)]));
S = (nu_max * r * r) * ([A * a, A * a * C_Pyr_to_ExIn * C_Pyr_to_ExIn * ...
    C_ExIn_to_Pyr, B * b * C_Pyr_to_InIn * C_Pyr_to_InIn * C_InIn_to_Pyr] .* ...
    S .* (S - 1) ./ ((S + 1).^3));
hess = zeros(6, 6, 6);
hess(4, 2, 2) = S(1);
hess(4, 3, 3) = S(1);
hess(4, 2, 3) = -S(1);
hess(4, 3, 2) = -S(1);
hess(5, 1, 1) = S(2);
hess(5, 2, 2) = C_Pyr_to_Pyr * S(1);
hess(5, 3, 3) = C_Pyr_to_Pyr * S(1);
hess(5, 2, 3) = -C_Pyr_to_Pyr * S(1);
hess(5, 3, 2) = -C_Pyr_to_Pyr * S(1);
hess(6, 1, 1) = S(3);
end

% ------------------------------------------------------------------------------
function hessp = hessiansp(t, Y, A, a, B, b, C_Pyr_to_Pyr, C_Pyr_to_ExIn, ...
    C_ExIn_to_Pyr, C_Pyr_to_InIn, C_InIn_to_Pyr, nu_max, r, v_0, v_Glu, v_GABA, ...
    mu_Glu_InIn_by_Pyr, q) %#ok
v = exp(r * [v_0 - v_Glu + v_GABA - Y(2) + Y(3), v_0 - C_Pyr_to_ExIn * Y(1), ...
    v_0 - mu_Glu_InIn_by_Pyr * v_Glu - C_Pyr_to_InIn * Y(1)]);
s = v ./ ((v + 1) .^ 2);
S = s ./ (v + 1);
hessp = zeros(6, 6, 16);

% A
hessp(4, 2, 1) = nu_max * r * a * s(1);
hessp(4, 3, 1) = -hessp(4, 2, 1);
hessp(5, 1, 1) = nu_max * r * a * C_Pyr_to_ExIn * C_ExIn_to_Pyr * s(2);
hessp(5, 2, 1) = C_Pyr_to_Pyr * hessp(4, 2, 1);
hessp(5, 3, 1) = -hessp(5, 2, 1);

% a
hessp(4, 1, 2) = -2 * a;
hessp(4, 2, 2) = nu_max * r * A * s(1);
hessp(4, 3, 2) = -hessp(4, 2, 2);
hessp(4, 4, 2) = -2;
hessp(5, 1, 2) = nu_max * r * A * C_Pyr_to_ExIn * C_ExIn_to_Pyr * s(2);
hessp(5, 2, 2) = C_Pyr_to_Pyr * hessp(4, 2, 2) - 2 * a;
hessp(5, 3, 2) = C_Pyr_to_Pyr * hessp(4, 3, 2);
hessp(5, 5, 2) = -2;

% B
hessp(6, 1, 3) = nu_max * r * b * C_Pyr_to_InIn * C_InIn_to_Pyr * s(3);

% b
hessp(6, 1, 4) = nu_max * r * B * C_Pyr_to_InIn * C_InIn_to_Pyr * s(3);
hessp(6, 3, 4) = -2 * b;
hessp(6, 6, 4) = -2;

% C_Pyr_to_Pyr
hessp(5, 2, 5) = nu_max * r * A * a * s(1);
hessp(5, 3, 5) = -hessp(5, 2, 5);

% C_Pyr_to_ExIn
hessp(5, 1, 6) = nu_max * r * A * a * C_ExIn_to_Pyr * (1 - r * ...
    C_Pyr_to_ExIn * Y(1) + (1 + r * C_Pyr_to_ExIn * Y(1)) * v(2)) * S(2);

% C_ExIn_to_Pyr
hessp(5, 1, 7) = nu_max * r * A * a * C_Pyr_to_ExIn * s(2);

% C_Pyr_to_InIn
hessp(6, 1, 8) = nu_max * r * B * b * C_InIn_to_Pyr * (1 - r * ...
    C_Pyr_to_InIn * Y(1) + (1 + r * C_Pyr_to_InIn * Y(1)) * v(3)) * S(3);

% C_InIn_to_Pyr
hessp(6, 1, 9) = nu_max * r * B * b * C_Pyr_to_InIn * s(3);

% nu_max
hessp(4, 2, 10) = r * A * a * s(1);
hessp(4, 3, 10) = -hessp(4, 2, 10);
hessp(5, 1, 10) = r * A * a * C_Pyr_to_ExIn * C_ExIn_to_Pyr * s(2);
hessp(5, 2, 10) = C_Pyr_to_Pyr * hessp(4, 2, 10);
hessp(5, 3, 10) = -hessp(5, 2, 10);
hessp(6, 1, 10) = r * B * b * C_Pyr_to_InIn * C_InIn_to_Pyr * s(3);

% r
hessp(4, 2, 11) = nu_max * A * a * (1 + r * (v_0 - v_Glu + v_GABA - Y(2) + ...
    Y(3)) + (1 - r * (v_0 - v_Glu + v_GABA - Y(2) + Y(3))) * v(1)) * S(1);
hessp(4, 3, 11) = -hessp(4, 2, 11);
hessp(5, 1, 11) = nu_max * A * a * C_Pyr_to_ExIn * C_ExIn_to_Pyr * (1 + r * ...
    (v_0 - C_Pyr_to_ExIn * Y(1)) + (1 - r * (v_0 - C_Pyr_to_ExIn * Y(1))) * v( ...
    2)) * S(2);
hessp(5, 2, 11) = C_Pyr_to_Pyr * hessp(4, 2, 11);
hessp(5, 3, 11) = -hessp(5, 2, 11);
hessp(6, 1, 11) = nu_max * B * b * C_Pyr_to_InIn * C_InIn_to_Pyr * (1 + r * ...
    (v_0 - mu_Glu_InIn_by_Pyr * v_Glu - C_Pyr_to_InIn * Y(1)) + (1 - r * (v_0 - ...
    mu_Glu_InIn_by_Pyr * v_Glu - C_Pyr_to_InIn * Y(1))) * v(3)) * S(3);

% v_0
hessp(4, 2, 12) = nu_max * r * A * a * r * (1 - v(1)) * S(1);
hessp(4, 3, 12) = -hessp(4, 2, 12);
hessp(5, 1, 12) = nu_max * r * A * a * C_Pyr_to_ExIn * C_ExIn_to_Pyr * r * ( ...
    1 - v(2)) * S(2);
hessp(5, 2, 12) = C_Pyr_to_Pyr * hessp(4, 2, 12);
hessp(5, 3, 12) = -hessp(5, 2, 12);
hessp(6, 1, 12) = nu_max * r * B * b * C_Pyr_to_InIn * C_InIn_to_Pyr * r * ( ...
    1 - v(3)) * S(3);

% v_Glu
hessp(4, 2, 13) = -nu_max * r * A * a * r * (1 - v(1)) * S(1);
hessp(4, 3, 13) = -hessp(4, 2, 13);
hessp(5, 2, 13) = C_Pyr_to_Pyr * hessp(4, 2, 13);
hessp(5, 3, 13) = -hessp(5, 2, 13);
hessp(6, 1, 13) = -nu_max * r * B * b * C_Pyr_to_InIn * C_InIn_to_Pyr * ...
    mu_Glu_InIn_by_Pyr * r * (1 - v(3)) * S(3);

% v_GABA
hessp(4, 2, 14) = nu_max * r * A * a * r * (1 - v(1)) * S(1);
hessp(4, 3, 14) = -hessp(4, 2, 14);
hessp(5, 2, 14) = C_Pyr_to_Pyr * hessp(4, 2, 14);
hessp(5, 3, 14) = -hessp(5, 2, 14);

% mu_Glu_InIn_by_Pyr
hessp(6, 1, 15) = -nu_max * r * B * b * C_Pyr_to_InIn * C_InIn_to_Pyr * ...
    v_Glu * r * (1 - v(3)) * S(3);

% q

end

% ------------------------------------------------------------------------------
function tens3 = der3(t, Y, A, a, B, b, C_Pyr_to_Pyr, C_Pyr_to_ExIn, ...
    C_ExIn_to_Pyr, C_Pyr_to_InIn, C_InIn_to_Pyr, nu_max, r, v_0, v_Glu, v_GABA, ...
    mu_Glu_InIn_by_Pyr, q) %#ok
S = exp(r * ([v_0 - v_Glu + v_GABA - Y(2) + Y(3), v_0 - C_Pyr_to_ExIn * Y( ...
    1), v_0 - mu_Glu_InIn_by_Pyr * v_Glu - C_Pyr_to_InIn * Y(1)]));
S = (nu_max * (r^3)) * ([A * a, A * a * (C_Pyr_to_ExIn^3) * C_ExIn_to_Pyr, ...
    B * b * (C_Pyr_to_InIn^3) * C_InIn_to_Pyr] .* S .* ((S - 2).^2 - 3) ./ ((S + ...
    1).^4));
tens3 = zeros(6, 6, 6, 6);
tens3(5, 1, 1, 1) = S(2);
tens3(6, 1, 1, 1) = S(3);
tens3(4, 2, 2, 2) = S(1);
tens3(4, 3, 3, 2) = S(1);
tens3(4, 2, 3, 2) = -S(1);
tens3(4, 3, 2, 2) = -S(1);
tens3(5, 2, 2, 2) = C_Pyr_to_Pyr * S(1);
tens3(5, 3, 3, 2) = C_Pyr_to_Pyr * S(1);
tens3(5, 2, 3, 2) = -C_Pyr_to_Pyr * S(1);
tens3(5, 3, 2, 2) = -C_Pyr_to_Pyr * S(1);
tens3(:, :, :, 3) = -tens3(:, :, :, 2);
end

% ------------------------------------------------------------------------------
function tens4 = der4(t, Y, A, a, B, b, C_Pyr_to_Pyr, C_Pyr_to_ExIn, ...
    C_ExIn_to_Pyr, C_Pyr_to_InIn, C_InIn_to_Pyr, nu_max, r, v_0, v_Glu, v_GABA, ...
    mu_Glu_InIn_by_Pyr, p) %#ok
end

% ------------------------------------------------------------------------------
function tens5 = der5(t, Y, A, a, B, b, C_Pyr_to_Pyr, C_Pyr_to_ExIn, ...
    C_ExIn_to_Pyr, C_Pyr_to_InIn, C_InIn_to_Pyr, nu_max, r, v_0, v_Glu, v_GABA, ...
    mu_Glu_InIn_by_Pyr, p) %#ok
end
