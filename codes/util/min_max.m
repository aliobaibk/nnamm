function I = min_max(varargin)
I = zeros(1, 2);
[I(1), I(2)] = bounds(varargin{:});
end
