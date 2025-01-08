function A = mask_where(A, mask, val)
if (nargin < 3)
    A = A(mask);
else
    if ((ischar(mask) || isstring(mask)) && strcmp(mask, 'all'))
        A(:) = val;
    else
        A(mask) = val;
    end
end
end
