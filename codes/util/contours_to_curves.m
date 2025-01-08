function [c, K] = contours_to_curves(c, n)
% c: contour matrix, e.g., obtained using countour()
% n: number of contour levels
c = c.';
K = [2, c(1, 2) + 1; zeros(n - 1, 2)];
for k = 1 : (n - 2)
    K(k + 1, :) = K(k, 2) + [2, c(K(k, 2) + 1, 2) + 1];
end
K(n, :) = [K(n - 1, 2) + 2, size(c, 1)];
c = arrayfun(@(a) c(K(a, 1) : K(a, 2), :), 1 : n, 'UniformOutput', false);
end
