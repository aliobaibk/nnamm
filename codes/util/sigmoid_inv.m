function X = sigmoid_inv(Y, s_max, s_rate, s_mode, s_off)
X = s_mode - log(s_max ./ (Y - s_off) - 1) / s_rate;
X(~is_in([s_off, s_off + s_max], '<', Y, '<')) = NaN;
end
