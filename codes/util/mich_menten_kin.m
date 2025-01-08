function MM = mich_menten_kin(X, v_max, K)
MM = v_max .* X ./ (K + X);
end
