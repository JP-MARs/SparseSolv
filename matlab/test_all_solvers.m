% Test all solvers
fprintf('\n=== Real Symmetric Matrix ===\n');

n = 30;
A = gallery('poisson', n);
b = ones(size(A, 1), 1);

% ICCG
[~, flag, relres, iter] = sparsesolv(A, b, 1e-8, 500, 'iccg');
fprintf('ICCG: flag=%d, relres=%.2e, iter=%d\n', flag, relres, iter);

% ICMRTR
[~, flag, relres, iter] = sparsesolv(A, b, 1e-8, 500, 'icmrtr');
fprintf('ICMRTR: flag=%d, relres=%.2e, iter=%d\n', flag, relres, iter);

% SGSMRTR
[~, flag, relres, iter] = sparsesolv(A, b, 1e-8, 500, 'sgsmrtr');
fprintf('SGSMRTR: flag=%d, relres=%.2e, iter=%d\n', flag, relres, iter);

fprintf('\n=== Complex Symmetric Matrix ===\n');

A_c = A + 0.1i * speye(size(A));
b_c = b + 0.5i * b;

% ICCOCG
[~, flag, relres, iter] = sparsesolv(A_c, b_c, 1e-6, 500, 'iccocg');
fprintf('ICCOCG: flag=%d, relres=%.2e, iter=%d\n', flag, relres, iter);

% ICMRTR (complex)
[~, flag, relres, iter] = sparsesolv(A_c, b_c, 1e-6, 500, 'icmrtr');
fprintf('ICMRTR (complex): flag=%d, relres=%.2e, iter=%d\n', flag, relres, iter);

% SGSMRTR (complex)
[~, flag, relres, iter] = sparsesolv(A_c, b_c, 1e-6, 500, 'sgsmrtr');
fprintf('SGSMRTR (complex): flag=%d, relres=%.2e, iter=%d\n', flag, relres, iter);

fprintf('\nAll solvers tested!\n');
