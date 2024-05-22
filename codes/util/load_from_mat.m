function S = load_from_mat(filepath, fields)
if ischar(fields)
    S = getfield(load(filepath, fields), fields);
elseif iscell(fields)
    S = load(filepath, fields{:});
else
    error('''fields'' is not a character array or cell array of character vectors')
end
end
