function S = matcont_hopf_stability(H, F)
mask_supc = H(end, :) < 0;

[r, c] = find(triu(ones(size(F, 1)), 1) == 1);
mask_n_inv = F;
mask_n_inv(abs(imag(mask_n_inv)) > 10^(-9)) = NaN;
mask_n_inv = round(real(mask_n_inv), 6);
mask_n_inv = (sum(mask_n_inv(r, :) + mask_n_inv(c, :) == 0) ~= 1);

S = repmat({'N'}, numel(mask_n_inv), 1);
S(mask_n_inv & mask_supc) = {'S'};
S(mask_n_inv & ~mask_supc) = {'U'};
end
