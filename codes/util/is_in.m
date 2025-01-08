function tf = is_in(lub, f_lb, x, f_ub)
% * Usage example:
% is_in([1, 3], @lt, 1 : 5, '<=')
tf = @(arg) (isequal(arg, '<') || isequal(arg, '<=') || isequal(arg, @lt) || ...
    isequal(arg, @le));
if ((numel(lub) ~= 2) || any(isnan(lub)) || (diff(lub) <= 0))
    error('Expected two increasing scalars for the lower and upper bounds...')
elseif (~tf(f_lb) || ~tf(f_ub))
    error('Accepted lower and upper bound operators are {<, <=, @lt, @le}...')
end
f = op_char_to_function_handle({f_lb, f_ub});
tf = f{1}(lub(1), x) & f{2}(x, lub(2));
end


%% ------------ Local functions ------------------------------------------------
function f = op_char_to_function_handle(f)
for k = 1 : numel(f)
    if ~ischar(f{k})
        continue
    end
    switch f{k}
        case '<'
            f{k} = @lt;
        case '<='
            f{k} = @le;
        otherwise
            error('Accepted operators are < and <=...')
    end
end
end
