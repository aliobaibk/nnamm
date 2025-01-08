function opt = matcont_set_cont_dir(curve_name, x0, v0, opt, p_sign)
n = opt.MaxNumPoints;
opt = contset(opt, 'Backward', 0);
opt = contset(opt, 'MaxNumPoints', 2);
x = cont(curve_name, x0, v0, opt);
opt = contset(opt, 'MaxNumPoints', n);
if (sign(x(end, 2) - x(end, 1)) ~= p_sign)
    opt = contset(opt, 'Backward', 1);
end
end
