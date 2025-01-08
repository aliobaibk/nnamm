function bd_draw_alpha()
% Draws multiple pseudo codimension-2 bifurcation diagrams, independently with
% respect to the parameter 'C_Pyr_to_Pyr', i.e., draws diagrams with (q, v_Glu,
% C_Pyr_to_Pyr)

proj_dir = './data+codes_pcbi2025_aliobaibk';
if ~isfolder(proj_dir)
    error(['The "proj_dir" variable is not set to an existing directory.\n', ...
        'Ensure "proj_dir" points to the directory where the repository was ', ...
        'cloned.\nFor example, if you cloned the repository into ', ...
        '"./data+codes_pcbi2025_aliobaibk",\n', ...
        'use: proj_dir = ''./data+codes_pcbi2025_aliobaibk''.'], '')
end

o_dir = fullfile(proj_dir, 'data', 'bifurcation-diagrams');

addpath(genpath(fullfile(proj_dir, 'codes', 'util')))
addpath(genpath(fullfile(proj_dir, 'codes', 'functions')))

flag_overwrite = false;


% - Main
create_dir(o_dir, false, false)

[~, out_name] = fileparts(mfilename('fullpath'));
out_name = strrep(out_name, 'bd_draw_', '');

W = init(out_name, o_dir);
main(W.d_param, W.e_param, W.par, W.h_E_Pyr, W.opt_bd_fig, W.out_dir, ...
    flag_overwrite)
end


%% ------------ Local functions ------------------------------------------------
function W = init(out_name, out_dir)
W.out_dir = fullfile(out_dir, out_name);
create_dir(out_dir, false, false)

W.d_param = { % Jansen–Rit parameterization
    'A', 3.25;
    'B', 22;
    'a', 100;
    'b', 50;
    'nu_max', 5;
    'r', 0.56;
    'v_0', 6;
    'v_GABA', 0;
    'mu_Glu_InIn_by_Pyr', 0.5;
    'C_Pyr_to_ExIn', 135;
    'C_ExIn_to_Pyr', 108;
    'C_Pyr_to_InIn', 33.75;
    'C_InIn_to_Pyr', 33.75;
    };

W.e_param = { % 'name', full, reduced
    'v_Glu', -2 : 0.005 : 2, [-1.5, 0.505];
    };

W.par = {
    'C_Pyr_to_Pyr', [0, 7.5];
    };

W.h_E_Pyr = 10^(-4);

W.opt_bd_fig = struct('var', 'E_Pyr', 'param', 'v_Glu', ...
    'lgd_param', 'v_{Glu} (mV)', 'lgd_title', 'C^{Pyr\rightarrowPyr}');
end


function main(d_param, e_param, par, h_E_Pyr, fig_opt, out_dir, flag_overwrite)
h_fig = figure('Color', [1, 1, 1], 'InvertHardcopy', 'off', ...
    'Position', [0, 0, 500, 400], 'Visible', 'off');
h_tile = nexttile(tiledlayout(h_fig, 1, 1, 'Padding', 'compact'));

d_param = cell2struct(d_param(:, 2), d_param(:, 1), 1);
for k = 1 : numel(par{2})
    cla(h_tile, 'reset')

    out_name = sprintf('%s-%g.bd', par{1}, par{2}(k));

    bif_file = fullfile(out_dir, sprintf('%s.mat', out_name));
    d_param.(par{1}) = par{2}(k);
    if (flag_overwrite || ~isfile(bif_file))
        BD = namm_bd(d_param, e_param([1, 2]), false, h_E_Pyr, false, bif_file);
    else
        BD = load(bif_file);
    end

    fig_file = fullfile(out_dir, sprintf('%s.png', out_name));
    fig_opt.title = sprintf('%s = %g', fig_opt.lgd_title, par{2}(k));
    namm_bd_display_3d(h_tile, BD.P, [], fig_opt)
    view(h_tile, -45, 30)
    exportgraphics(h_fig, fig_file, 'Resolution', 150)

    bif_file_py = fullfile(out_dir, sprintf('%s.pydstool.mat', out_name));
    if (flag_overwrite || ~isfile(bif_file_py))
        BD = namm_bd_pydstool(BD, e_param{1}, e_param{3});
        save(bif_file_py, '-struct', 'BD', '-v7.3')
    end
end

close(h_fig)
end
