function Y = sigmoid(X, s_max, s_rate, s_mode, s_off)
Y = s_max ./ (1 + exp(s_rate .* (s_mode - X))) + s_off;
end
