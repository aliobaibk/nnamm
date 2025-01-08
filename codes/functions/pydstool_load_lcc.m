function [L, discard] = pydstool_load_lcc(bd_file, BD, mask_bd, opt)
% Loads LCC saved by PyDSTool
% opt.var: variables to extract from LCC (L.X: {'q'}; L.Y: opt.var])
% opt.par: parameters to extract from BD (L.Y: opt.par)
% opt.ver: LCC directory name suffix

% Usage example:
% > bd_file = 'bd.mat';
% > BD = load_from_mat(bd_file, {'explored_param', 'P'});
% > mask_bd = [1, 4, 5];
% > opt = struct('var', {{'E_Pyr_min', 'E_Pyr_max'}}, 'par', {'v_Glu'}, ...
%   'ver', 'lcc');
% > L = pydstool_load_lcc(bd_file, BD, mask_bd, opt);

[in_dir, name] = fileparts(bd_file);
if isempty(BD)
    BD = load_from_mat(bd_file, {'explored_param', 'P'});
end

if ~isfield(opt, 'par')
    opt.par = '';
end
if (~isempty(opt.par) && ~ismember(opt.par, BD.explored_param))
    error('Parameter ''%s'' was not explored...', opt.par)
end

BD = BD.P;
if iscell(BD)
    BD = cell2mat(BD);
end
if ~isempty(mask_bd)
    BD = BD(mask_bd);
else
    mask_bd = (1 : numel(BD));
end
BD = struct('label', {BD.bif_label}, 'par', {BD.param});

n = numel(BD);
L = cell(n, 1);
lcc = struct('X', [], 'Y', [], 'Z', [], 'stab', {{}}, 'type', {{}}, ...
    'diagram', [], 'name', '', 'hopf', []);
for k_d = 1 : n
    label = BD(k_d).label.name;
    mask_h = find(label == 'h');
    n_h = numel(mask_h);
    if (n_h == 0)
        continue
    end

    L{k_d} = duplicate_object(lcc, n_h);
    discard = false(n_h, 1);
    for k_h = 1 : n_h
        lcc_file = fullfile(in_dir, sprintf('%s_%s', name, opt.ver), sprintf( ...
            '%d-%d.lcc', mask_bd(k_d), mask_h(k_h)));
        if ~isfile(lcc_file)
            discard(k_h) = true;
            continue
        end
        D = readtable(lcc_file, 'FileType', 'text', 'VariableNamingRule', 'preserve');

        L{k_d}(k_h).X = D.q;
        if ~isempty(opt.par)
            L{k_d}(k_h).Y = BD(k_d).par.(opt.par);
        end
        L{k_d}(k_h).Z = D(:, opt.var).Variables;
        L{k_d}(k_h).stab = D.stability;
        L{k_d}(k_h).type = D.type;
        L{k_d}(k_h).diagram = mask_bd(k_d);
        L{k_d}(k_h).name = label;
        L{k_d}(k_h).hopf = [L{k_d}(k_h).X(1), L{k_d}(k_h).Y, L{k_d}(k_h).Z(1, :)];
    end
    L{k_d}(discard) = [];
end
discard = cellfun(@isempty, L, 'UniformOutput', true);
L(discard) = [];
end
