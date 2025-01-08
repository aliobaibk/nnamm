function D = namm_setup_bd(D, par, var)
% * Extracts (X: 'q', Y: par.name, Z: var).
% * X is adjusted with 'v_GABA' if needed.
D = D(:).';
X = cell2mat(cellfun(@(d) d.q, D, 'UniformOutput', false));
Y = cell2mat(cellfun(@(d) d.param.(par.name) * ones(size(d.E_Pyr)), D, ...
    'UniformOutput', false));
Z = struct( ...
    'E_Pyr', cell2mat(cellfun(@(d) d.E_Pyr, D, 'UniformOutput', false)), ...
    'E_ExIn_and_Pyr', cell2mat(cellfun(@(d) d.E_ExIn_and_Pyr, D, ...
    'UniformOutput', false)), ...
    'E_InIn', cell2mat(cellfun(@(d) d.E_InIn, D, 'UniformOutput', false)), ...
    'LFP', cell2mat(cellfun(@(d) d.LFP, D, 'UniformOutput', false)));
if isfield(par, 'lim')
    mask = is_in(par.lim, '<=', Y(1, :), '<=');
else
    mask = true(1, size(Y, 2));
end
D = struct('X', X(:, mask), 'Y', Y(:, mask), 'Z', Z.(var)(:, mask));
end
