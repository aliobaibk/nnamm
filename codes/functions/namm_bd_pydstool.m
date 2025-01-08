function BD = namm_bd_pydstool(BD, par_name, v)
% From a pseudo codimension-2 bifurcation diagram, filters number of diagrams to
% speed up PyDSTool
F = cellfun(@(bd) bd.param.(par_name), BD.P, 'UniformOutput', true);
F = find(is_in(v, '<=', F, '<='));
BD = struct('P', {arrayfun(@(k) setfield(rmfield(BD.P{k}, {'eigen_values', ...
    'n_neg_real_part'}), 'original', k), F, 'UniformOutput', false)});
end
