function P = cartesian_product(varargin)
% * Usage example:
% cartesian_product(1, 3, 2)
N = cellfun(@(n) 1 : n, varargin, 'UniformOutput', false);
P = cell(1, nargin);
[P{1 : nargin}] = ndgrid(N{:});
P = cell2mat(cellfun(@(e) e(:), P, 'UniformOutput', false));
end
