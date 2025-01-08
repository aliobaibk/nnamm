function [R, var_z] = matcont_saddle_node(bd_dir, bd_name, var_z, BD, SN)
% Full Hopf curve

if ((nargin < 3) || isempty(var_z))
    var_z = 1;
end
if (nargin < 4)
    BD = [];
end
if (nargin < 5)
    SN = [];
end

if isempty(BD)
    in_file = fullfile(bd_dir, sprintf('%s.bd.mat', bd_name));
    BD = load_from_mat(in_file, {'n_ode'});
end

var = [BD.n_ode + [1; 2]; var_z(:)];
n_v = numel(var);

if isempty(SN)
    in_file = fullfile(bd_dir, sprintf('%s.saddle-node.mat', bd_name));
    R = load_from_mat(in_file, 'Results');
else
    R = SN;
end
R = [R(strcmp({R.proc}, 'lp-f')); R(strcmp({R.proc}, 'lp-b'))];
XYZ = [fliplr(R(1).X(var, :)), R(2).X(var, :)].';

R = struct('X', XYZ(:, 1), 'Y', XYZ(:, 2), 'Z', XYZ(:, 3 : n_v));
end
