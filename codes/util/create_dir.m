function create_dir(dir_path, verbose, warning_on)
if ((nargin < 2) || isempty(verbose))
    verbose = true;
end
if ((nargin < 3) || isempty(warning_on))
    warning_on = true;
end
if ~exist(dir_path, 'dir')
    if verbose
        fprintf('- Creating the non existing output directory:\n%s\n', dir_path)
    end
    [status, msg] = mkdir(dir_path);
    if ~status
        error('- Failed to create the output directory, \n%s\n', msg)
    elseif verbose
        fprintf('- Output directory successfully created\n\n');
    end
elseif warning_on
    warning(['\n- Old files might get replaced in the already existing output ', ...
        'directory:\n%s\n\n'], dir_path);
end
end
