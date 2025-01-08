function namm_bd_display_3d(h_tile, bd, c, opt)
% Partial pseudo codimension-2 bifurcation diagram ('bd.q', 'opt.param')
% * Only the following points are displayed: equilibrium, Hopf, saddle-node
n_ode = 6;

bd = bd(:).';

colors = 3;
if (isfield(opt, 'show_bp') && opt.show_bp)
    colors = colors + 1;
end
colors = lines(colors);
if (isfield(opt, 'colors') && ~isempty(opt.colors))
    colors = opt.colors;
end

% 'x' = 'q', 'y' = 'opt.param', 'z' = 'opt.var'
x = cell2mat(cellfun(@(d) d.q, bd, 'UniformOutput', false));
y = cell2mat(cellfun(@(d) d.param.(opt.param) * ones(size(d.q)), bd, ...
    'UniformOutput', false));
switch opt.var
    case {'E_Pyr'; 'E_ExIn_and_Pyr'; 'E_InIn'; 'LFP'}
        z = cell2mat(cellfun(@(d) d.(opt.var), bd, 'UniformOutput', false));
    otherwise
        error('Unknown state variable ''%s''...', opt.var)
end


% - Equilibrium and bifurcation points
hold(h_tile, 'on')

if isempty(c)
    K = n_ode - cell2mat(cellfun(@(d) d.n_neg_real_part, bd, ...
        'UniformOutput', false));
    c = reshape(colors(K + 1, :), [size(z), 3]);

    if (isfield(opt, 'show_bp') && opt.show_bp)
        K = arrayfun(@(k) [cell2mat(bd{k}.bif_label.k_sing(:, 2)), repmat(k, size( ...
            bd{k}.bif_label.k_sing, 1), 1)], 1 : numel(bd), 'UniformOutput', false);
        K = cat(1, K{:});
        K = sub2ind(size(z), K(:, 1), K(:, 2)) + [0, 1, 2] * numel(z);
        c(K) = repmat(colors(end, :), [size(K, 1), 1]);
    end
end
surf(h_tile, x, y, z, c, 'EdgeColor', 'interp', 'FaceColor', 'interp', ...
    'LineStyle', 'none')

hold(h_tile, 'off')


% - Axis formatting
xlabel(h_tile, 'q (Hz)', 'FontSize', 9)
ylabel(h_tile, opt.lgd_param, 'FontSize', 9)
ytickformat(h_tile, '%.2e')
switch opt.var
    case 'E_Pyr'
        zlabel(h_tile, 'E_{Pyr} (mV)', 'FontSize', 9)
    case 'E_ExIn_and_Pyr'
        zlabel(h_tile, 'E_{ExIn∪Pyr} (mV)', 'FontSize', 9)
    case 'E_InIn'
        zlabel(h_tile, 'E_{InIn} (mV)', 'FontSize', 9)
    case 'LFP'
        zlabel(h_tile, 'LFP (mV)', 'FontSize', 9)
end
ztickformat(h_tile, '%.2e')

axis(h_tile, 'tight')
box(h_tile, 'on')

h_tile.TickDir = 'both';

if (isfield(opt, 'title') && ~isempty(opt.title))
    title(h_tile, opt.title, 'FontSize', 12)
end
end
