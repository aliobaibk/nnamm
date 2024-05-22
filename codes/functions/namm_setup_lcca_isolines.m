function C = namm_setup_lcca_isolines(bd_dir, bd_name, opt, LCC)
% Contour lines of limit cycle peak–peak amplitudes

if (nargin < 4)
    LCC = [];
end

if isempty(LCC)
    in_file = fullfile(bd_dir, sprintf('%s.lcc.mat', bd_name));
    LCC = load(in_file);
end

if ~ismember(opt.var, LCC.var)
    error('Unknown state variable: %s...', opt.var)
elseif ~ismember(opt.x_par, LCC.par)
    error('Unknown x-parameter: %s...', opt.x_par)
elseif ~ismember(opt.y_par, LCC.par)
    error('Unknown y-parameter: %s...', opt.y_par)
end

X_q = LCC.(opt.var).(opt.x_par);
Y_q = LCC.(opt.var).(opt.y_par);

P = [opt.f_x_inv(opt.levels(:, 1)), opt.f_y_inv(opt.levels(:, 2))];

z = LCC.(opt.var).max - LCC.(opt.var).min;
[~, m] = min((X_q(:) - P(:, 1).').^2 + (Y_q(:) - P(:, 2).').^2, [], 1);
levels = unique(z(m), 'sorted');
C = contourc(X_q(1, :), Y_q(:, 1), z, levels);
[~, C_ind] = contours_to_curves(C, numel(levels));
C(:, C_ind(:, 1) - 1) = NaN;

z = [LCC.(opt.var).max(:), LCC.(opt.var).min(:)];
[~, m] = min((X_q(:) - C(1, :)).^2 + (Y_q(:) - C(2, :)).^2, [], 1);
C = [C.', z(m, :)];

C = struct('curves', C, 'levels', levels, 'ind', C_ind);
end
