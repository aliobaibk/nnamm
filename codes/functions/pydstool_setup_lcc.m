function lcc = pydstool_setup_lcc(bd_file, BD, v_par, opt)
opt.var = [strcat(opt.var, {'_max'; '_min'}); {'period'}];
mask = cellfun(@(c) c.param.(opt.par), BD.P, 'UniformOutput', true);
[~, mask] = min(abs(mask(:) - v_par(:).'));
L = cell2mat(pydstool_load_lcc(bd_file, BD, mask, opt));
[L.Y] = cell_expand(arrayfun(@(a) a.Y * ones(size(a.X)), L, ...
    'UniformOutput', false));
n_c = numel(L);
lcc = cell(n_c, 1);
for k_c = 1 : n_c
    lcc{k_c} = struct('X', L(k_c).X, 'Y', L(k_c).Y, 'Z', L(k_c).Z(:, [1, 2]), ...
        'T', L(k_c).Z(:, 3));
end
lcc = cat(1, lcc{:});
end
