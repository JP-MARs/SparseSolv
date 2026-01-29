function [x, flag, relres, iter, resvec] = sparsesolv(A, b, varargin)
% SPARSESOLV Solve sparse symmetric linear system using SparseSolv library
%
% Solves the linear system Ax = b where A is a symmetric (or complex symmetric)
% sparse matrix using iterative methods from the SparseSolv library.
%
% SYNTAX:
%   x = sparsesolv(A, b)
%   x = sparsesolv(A, b, tol)
%   x = sparsesolv(A, b, tol, maxit)
%   x = sparsesolv(A, b, tol, maxit, method)
%   x = sparsesolv(A, b, tol, maxit, method, accel)
%   x = sparsesolv(A, b, tol, maxit, method, accel, x0)
%   [x, flag] = sparsesolv(...)
%   [x, flag, relres] = sparsesolv(...)
%   [x, flag, relres, iter] = sparsesolv(...)
%   [x, flag, relres, iter, resvec] = sparsesolv(...)
%
% INPUTS:
%   A      - Symmetric sparse matrix (n x n)
%   b      - Right-hand side vector (n x 1)
%   tol    - Convergence tolerance (default: 1e-6)
%   maxit  - Maximum iterations (default: 1000)
%   method - Solver method: 'iccg' (default), 'iccocg', 'icmrtr', 'sgsmrtr'
%   accel  - Acceleration/shift parameter (default: 1.0, negative for auto)
%   x0     - Initial guess (default: zeros)
%
% OUTPUTS:
%   x      - Solution vector
%   flag   - Convergence flag: 0 = converged, 1 = not converged
%   relres - Final relative residual ||b - Ax|| / ||b||
%   iter   - Number of iterations performed
%   resvec - Residual history at each iteration
%
% SOLVER METHODS:
%   'iccg'    - Incomplete Cholesky Conjugate Gradient (real symmetric)
%   'iccocg'  - IC + Conjugate Orthogonal CG (complex symmetric)
%   'icmrtr'  - IC + MRTR method (real symmetric)
%   'sgsmrtr' - Symmetric Gauss-Seidel + MRTR (real symmetric)
%
% OPTIONS (via name-value pairs after positional arguments):
%   'DiagScaling'  - Enable diagonal scaling (default: false)
%   'SaveBest'     - Save best result for divergence recovery (default: true)
%   'DivergeType'  - Divergence detection type (default: 1)
%   'DivergeVal'   - Divergence threshold multiplier (default: 10.0)
%   'DivergeCount' - Divergence count threshold (default: 10)
%
% EXAMPLES:
%   % Simple usage
%   A = gallery('poisson', 50);
%   b = ones(size(A, 1), 1);
%   x = sparsesolv(A, b);
%
%   % With options
%   [x, flag, relres, iter] = sparsesolv(A, b, 1e-8, 500, 'iccg', 1.05);
%
%   % Complex symmetric matrix
%   A = A + 1i * speye(size(A));
%   x = sparsesolv(A, b, 1e-6, 500, 'iccocg');
%
% COMSOL COMPATIBILITY:
%   This function is designed to be compatible with COMSOL LiveLink for MATLAB.
%   The input matrix A can come directly from COMSOL's mphmatrix function.
%
% See also: PCG, BICGSTAB, SYMMETRIC_SPARSE_TO_CSR

    % Input validation
    narginchk(2, inf);

    if ~issparse(A)
        A = sparse(A);
    end

    [m, n] = size(A);
    if m ~= n
        error('sparsesolv:notSquare', 'Matrix A must be square');
    end

    if length(b) ~= n
        error('sparsesolv:sizeMismatch', 'Vector b must have length equal to matrix dimension');
    end
    b = b(:);  % Ensure column vector

    % Parse optional arguments
    p = inputParser;
    addOptional(p, 'tol', 1e-6, @isnumeric);
    addOptional(p, 'maxit', 1000, @isnumeric);
    addOptional(p, 'method', 'iccg', @ischar);
    addOptional(p, 'accel', 1.0, @isnumeric);
    addOptional(p, 'x0', [], @isnumeric);
    addParameter(p, 'DiagScaling', false, @(x) islogical(x) || isnumeric(x));
    addParameter(p, 'SaveBest', true, @(x) islogical(x) || isnumeric(x));
    addParameter(p, 'DivergeType', 1, @isnumeric);
    addParameter(p, 'DivergeVal', 10.0, @isnumeric);
    addParameter(p, 'DivergeCount', 10, @isnumeric);

    parse(p, varargin{:});

    tol = p.Results.tol;
    maxit = p.Results.maxit;
    method = lower(p.Results.method);
    accel = p.Results.accel;
    x0 = p.Results.x0;

    % Build options structure for MEX
    options = struct();
    options.diag_scaling = logical(p.Results.DiagScaling);
    options.save_best = logical(p.Results.SaveBest);
    options.diverge_type = p.Results.DivergeType;
    options.diverge_val = p.Results.DivergeVal;
    options.diverge_count = p.Results.DivergeCount;

    % Validate method
    valid_methods = {'iccg', 'iccocg', 'icmrtr', 'sgsmrtr'};
    if ~ismember(method, valid_methods)
        error('sparsesolv:invalidMethod', ...
            'Invalid method ''%s''. Valid methods: %s', ...
            method, strjoin(valid_methods, ', '));
    end

    % Auto-detect complex and switch to ICCOCG if needed
    is_complex = ~isreal(A) || ~isreal(b);
    if is_complex && strcmp(method, 'iccg')
        method = 'iccocg';
        warning('sparsesolv:autoComplex', ...
            'Complex input detected, switching to ICCOCG method');
    end

    % Convert to CSR format (lower triangular only)
    [vals, col_ind, row_ptr] = symmetric_sparse_to_csr(A);

    % Prepare initial guess
    if isempty(x0)
        x0 = [];
    else
        x0 = x0(:);
        if length(x0) ~= n
            error('sparsesolv:x0Size', 'Initial guess x0 must have length %d', n);
        end
    end

    % Call MEX function
    try
        [x, flag, relres, iter, resvec] = sparsesolv_mex(...
            method, vals, col_ind, row_ptr, b, tol, maxit, accel, x0, options);
    catch ME
        % If MEX not found, try to build it
        if contains(ME.identifier, 'UndefinedFunction') || ...
           contains(ME.message, 'sparsesolv_mex')
            error('sparsesolv:mexNotFound', ...
                ['MEX file not found. Run build_mex.m first to compile.\n' ...
                 'Error: %s'], ME.message);
        else
            rethrow(ME);
        end
    end

    % Ensure output is column vector
    x = x(:);
end
