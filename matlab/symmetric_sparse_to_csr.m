function [vals, col_ind, row_ptr] = symmetric_sparse_to_csr(A)
% SYMMETRIC_SPARSE_TO_CSR Convert symmetric sparse matrix to CSR format
%
% Extracts only the lower triangular part of a symmetric sparse matrix
% and converts it to Compressed Sparse Row (CSR) format.
%
% SYNTAX:
%   [vals, col_ind, row_ptr] = symmetric_sparse_to_csr(A)
%
% INPUTS:
%   A - Symmetric sparse matrix (full or sparse MATLAB matrix)
%
% OUTPUTS:
%   vals    - Non-zero values (lower triangular part only)
%   col_ind - Column indices (0-based for C++ compatibility)
%   row_ptr - Row pointers (0-based)
%
% NOTES:
%   - Input matrix must be square and symmetric
%   - Only lower triangular part (including diagonal) is extracted
%   - Output indices are 0-based for direct use with C++ MEX functions
%
% EXAMPLE:
%   A = sparse([4 -1 0; -1 4 -1; 0 -1 4]);
%   [vals, col_ind, row_ptr] = symmetric_sparse_to_csr(A);
%
% See also: SPARSESOLV, SPARSESOLV_MEX, TRIL

    [m, n] = size(A);

    % Check square matrix
    if m ~= n
        error('symmetric_sparse_to_csr:notSquare', ...
            'Input matrix must be square (got %dx%d)', m, n);
    end

    % Check symmetry (optional, can be disabled for performance)
    if ~issparse(A)
        A = sparse(A);
    end

    % Symmetry check (with tolerance for numerical errors)
    tol = 1e-12 * max(abs(nonzeros(A)));
    if norm(A - A.', 'fro') > tol
        warning('symmetric_sparse_to_csr:notSymmetric', ...
            'Input matrix may not be symmetric');
    end

    % Extract lower triangular part
    L = tril(A);

    % Get COO format (row, col, val)
    [rows, cols, values] = find(L);

    % Sort by row (MATLAB's find already returns sorted by column, then row)
    % We need to sort by row first
    [rows_sorted, sort_idx] = sort(rows);
    cols_sorted = cols(sort_idx);
    vals = values(sort_idx);

    % Convert to 0-based indexing
    col_ind = cols_sorted - 1;

    % Build row_ptr
    row_ptr = zeros(m + 1, 1);
    for k = 1:length(rows_sorted)
        row_ptr(rows_sorted(k) + 1) = row_ptr(rows_sorted(k) + 1) + 1;
    end
    row_ptr = cumsum(row_ptr);

    % Ensure column vector output
    vals = vals(:);
    col_ind = col_ind(:);
    row_ptr = row_ptr(:);
end
