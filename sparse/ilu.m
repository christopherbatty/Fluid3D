function [l, u] = ilu(a, level, omega, storage)
%ILU   Incomplete LU and Cholesky factorization by level-of-fill.
%   [L, U] = ILU(A, LEVEL, OMEGA, STORAGE) computes an
%   incomplete factorization of A based on level-of-fill LEVEL.
%   L is lower triangular, and U is unit upper triangular.
%
%   If only one output argument is used, then the matrix is assumed to be 
%   SPD and the lower triangular incomplete Cholesky factor is returned.
%   The Cholesky version computes the nonsymmetric factorization internally.
%
%   The parameters OMEGA and STORAGE are optional, but both must be
%   specified if STORAGE is specified.
%
%   OMEGA is the relaxation factor for relaxed factorizations.
%   A value of 0 is the standard factorization (default) and a value
%   of 1 produces a modified factorization.  0.95 is often a good value
%   for some elliptic PDE problems.
%
%   STORAGE is the number of nonzeros allocated for each factor L and U.
%   This parameter should be specified if the function fails and indicates
%   that the default value is too small.
%
%   Examples:
%   [l u] = ilu(a, 0);
%   l = ilu(a, 0);
%   [l u] = ilu(a, n, 0, n*(n+1)/2);

%   Version of 2-12-02, Edmond Chow

if (nargin < 3)
    omega = 0;
end

if (nargin < 4)
    % set default value
    storage = max(nnz(tril(a)),nnz(triu(a))) + length(a);
    if (level > 0)
        storage = storage*(level+2);
    end
end

if (issparse(a) == 0)
    a = sparse(a); 
end

% call the mex file to perform the factorization
[l u] = ilu_mex(a, level, omega, storage);

if (nargout == 1) 

    % incomplete Cholesky case

    % l is lower triangular and u is unit upper triangular
    d = diag(l);
    if (sum(d <= 0) > 0) error('IC factorization: matrix is not SPD'); end
    d = diag(sqrt(d));
    l = l*inv(d);

    %use to check that l=u'
    %u = u + speye(size(a));
    %u = d*u;

else 

    % incomplete LU case

    u = u + speye(size(a));

end
