function [R, var_z] = matcont_special_points(H, SN, BD, var_z)
% Full Hopf curve

if ((nargin < 4) || isempty(var_z))
    var_z = 1;
end

if isempty(BD)
    in_file = fullfile(bd_dir, sprintf('%s.bd.mat', bd_name));
    BD = load_from_mat(in_file, {'n_ode'});
end

var = [BD.n_ode + [1; 2]; var_z(:)];
n_v = numel(var);

R = cell(0, 2);
for k = 2 : numel(H)
    s = H(k).S(2 : (end - 1));
    R = [R; arrayfun(@(a) {a.label, H(k).X(var, a.index).'}, s(:), ...
        'UniformOutput', false)]; %#ok
end
for k = 2 : numel(SN)
    s = SN(k).S(2 : (end - 1));
    R = [R; arrayfun(@(a) {a.label, SN(k).X(var, a.index).'}, s(:), ...
        'UniformOutput', false)]; %#ok
end
R = cat(1, R{:});
if isempty(R)
    return
end
R = {R(:, 1), cat(1, R{:, 2})};
R = table(R{1}, R{2}(:, 1), R{2}(:, 2), R{2}(:, 3 : n_v), ...
    'VariableNames', {'label', 'X', 'Y', 'Z'});
end
