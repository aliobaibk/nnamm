function [R, var_z] = matcont_hopf(bd_dir, bd_name, flag_no_neutral, var_z, ...
    BD, H, SN)
% Full Hopf curve

if (nargin < 3)
    flag_no_neutral = true;
end
if ((nargin < 4) || isempty(var_z))
    var_z = 1;
end
if (nargin < 5)
    BD = [];
end
if (nargin < 6)
    H = [];
end
if (nargin < 7)
    SN = [];
end

if isempty(BD)
    in_file = fullfile(bd_dir, sprintf('%s.bd.mat', bd_name));
    BD = load_from_mat(in_file, {'n_ode'});
end

var = [BD.n_ode + [1; 2]; var_z(:)];
n_v = numel(var);

if isempty(H)
    in_file = fullfile(bd_dir, sprintf('%s.hopf.mat', bd_name));
    R = load_from_mat(in_file, 'Results');
else
    R = H;
end
R = [R(strcmp({R.proc}, 'hp-b')); R(strcmp({R.proc}, 'hp-f'))];
s = [fliplr(R(1).H(end, :)), R(2).H(end, :); fliplr(R(1).F), R(2).F];
s = matcont_hopf_stability(s(1, :), s(2 : end, :)).';
if flag_no_neutral
    mask = ~strcmp(s, 'N');
else
    mask = true(size(s));
end
S = s(mask);
XYZ = reshape(mask_where([fliplr(R(1).X(var, :)), R(2).X(var, :)], repmat( ...
    mask, n_v, 1)), n_v, sum(mask));

if isempty(SN)
    in_file = fullfile(bd_dir, sprintf('%s.saddle-node.mat', bd_name));
    R = load_from_mat(in_file, 'Results');
else
    R = SN;
end
R = [R(strcmp({R.proc}, 'bt-f')); R(strcmp({R.proc}, 'bt-b'))];
s = [fliplr(R(1).H(end, :)), R(2).H(end, :); fliplr(R(1).F), R(2).F];
s = matcont_hopf_stability(s(1, :), s(2 : end, :)).';
if flag_no_neutral
    mask = ~strcmp(s, 'N');
else
    mask = true(size(s));
end
S = [S, {''}, s(mask)].';
XYZ = [XYZ, NaN(n_v, 1), reshape(mask_where([fliplr(R(1).X(var, :)), R(2).X( ...
    var, :)], repmat(mask, n_v, 1)), n_v, sum(mask))].';

R = struct('X', XYZ(:, 1), 'Y', XYZ(:, 2), 'Z', XYZ(:, 3 : n_v), 'S', {S});
end
