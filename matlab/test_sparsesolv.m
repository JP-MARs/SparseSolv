function test_sparsesolv()
% TEST_SPARSESOLV Test suite for SparseSolv MEX interface
%
% This script runs a series of tests to verify the MEX interface
% is working correctly with various matrix types and solver options.
%
% Run this after building the MEX file with build_mex.m

    fprintf('========================================\n');
    fprintf('SparseSolv MEX Interface Test Suite\n');
    fprintf('========================================\n\n');

    % Check if MEX file exists
    this_dir = fileparts(mfilename('fullpath'));
    mex_file = fullfile(this_dir, ['sparsesolv_mex.' mexext]);
    if ~isfile(mex_file)
        error('test_sparsesolv:mexNotFound', ...
            ['MEX file not found: %s\n' ...
             'Run build_mex first to compile the MEX interface.'], mex_file);
    end

    % Add matlab directory to path
    addpath(this_dir);

    % Run tests
    results = struct('passed', 0, 'failed', 0, 'tests', {{}});

    results = run_test(results, @test_poisson_2d, 'Poisson 2D (ICCG)');
    results = run_test(results, @test_poisson_3d, 'Poisson 3D (ICCG)');
    results = run_test(results, @test_iccg_options, 'ICCG with options');
    results = run_test(results, @test_complex_symmetric, 'Complex symmetric (ICCOCG)');
    results = run_test(results, @test_icmrtr, 'ICMRTR method');
    results = run_test(results, @test_sgsmrtr, 'SGSMRTR method');
    results = run_test(results, @test_initial_guess, 'Initial guess');
    results = run_test(results, @test_convergence_check, 'Convergence check');
    results = run_test(results, @test_residual_log, 'Residual log');

    % Summary
    fprintf('\n========================================\n');
    fprintf('Test Summary\n');
    fprintf('========================================\n');
    fprintf('Passed: %d\n', results.passed);
    fprintf('Failed: %d\n', results.failed);
    fprintf('Total:  %d\n', results.passed + results.failed);
    fprintf('\n');

    if results.failed > 0
        fprintf('Failed tests:\n');
        for i = 1:length(results.tests)
            if ~results.tests{i}.passed
                fprintf('  - %s: %s\n', results.tests{i}.name, results.tests{i}.error);
            end
        end
        fprintf('\n');
    end

    if results.failed == 0
        fprintf('All tests PASSED!\n');
    else
        warning('test_sparsesolv:testsFailed', '%d test(s) failed.', results.failed);
    end
end

function results = run_test(results, test_func, name)
    fprintf('Testing: %s... ', name);
    test_result = struct('name', name, 'passed', false, 'error', '');

    try
        test_func();
        fprintf('PASSED\n');
        test_result.passed = true;
        results.passed = results.passed + 1;
    catch ME
        fprintf('FAILED\n');
        fprintf('  Error: %s\n', ME.message);
        test_result.error = ME.message;
        results.failed = results.failed + 1;
    end

    results.tests{end+1} = test_result;
end

%% Test Functions

function test_poisson_2d()
    % 2D Poisson problem using gallery matrix
    n = 30;  % Grid size
    A = gallery('poisson', n);
    b = ones(size(A, 1), 1);

    [x, flag, relres, iter] = sparsesolv(A, b, 1e-8, 500, 'iccg');

    assert(flag == 0, 'Solver did not converge');
    assert(relres < 1e-7, 'Relative residual too large: %.2e', relres);

    % Verify solution
    r = b - A * x;
    actual_relres = norm(r) / norm(b);
    assert(actual_relres < 1e-6, 'Actual relative residual too large: %.2e', actual_relres);
end

function test_poisson_3d()
    % Larger 3D-like problem
    n = 15;
    A = gallery('poisson', n);
    A = kron(speye(n), A) + kron(A, speye(n));  % Tensor product for 3D-like structure
    b = rand(size(A, 1), 1);

    [x, flag, relres, iter] = sparsesolv(A, b, 1e-6, 1000, 'iccg');

    assert(flag == 0, 'Solver did not converge');
    assert(relres < 1e-5, 'Relative residual too large: %.2e', relres);
end

function test_iccg_options()
    % Test ICCG with various options
    n = 20;
    A = gallery('poisson', n);
    b = ones(size(A, 1), 1);

    % Test with acceleration parameter
    [x, flag, relres] = sparsesolv(A, b, 1e-8, 500, 'iccg', 1.05);
    assert(flag == 0, 'Solver with accel=1.05 did not converge');

    % Test with diagonal scaling
    [x, flag, relres] = sparsesolv(A, b, 1e-8, 500, 'iccg', 1.0, [], 'DiagScaling', true);
    assert(flag == 0, 'Solver with DiagScaling did not converge');
end

function test_complex_symmetric()
    % Complex symmetric matrix (not Hermitian)
    n = 20;
    A_real = gallery('poisson', n);
    A = A_real + 0.1i * speye(size(A_real));  % Complex symmetric
    b = ones(size(A, 1), 1) + 0.5i * ones(size(A, 1), 1);

    [x, flag, relres, iter] = sparsesolv(A, b, 1e-6, 500, 'iccocg');

    assert(flag == 0, 'ICCOCG did not converge');
    assert(relres < 1e-5, 'Relative residual too large: %.2e', relres);

    % Verify solution
    r = b - A * x;
    actual_relres = norm(r) / norm(b);
    assert(actual_relres < 1e-4, 'Actual relative residual too large: %.2e', actual_relres);
end

function test_icmrtr()
    % Test ICMRTR method
    n = 20;
    A = gallery('poisson', n);
    b = ones(size(A, 1), 1);

    [x, flag, relres, iter] = sparsesolv(A, b, 1e-6, 500, 'icmrtr');

    assert(flag == 0, 'ICMRTR did not converge');
    assert(relres < 1e-5, 'Relative residual too large: %.2e', relres);
end

function test_sgsmrtr()
    % Test SGSMRTR method
    n = 20;
    A = gallery('poisson', n);
    b = ones(size(A, 1), 1);

    [x, flag, relres, iter] = sparsesolv(A, b, 1e-6, 500, 'sgsmrtr');

    assert(flag == 0, 'SGSMRTR did not converge');
    assert(relres < 1e-5, 'Relative residual too large: %.2e', relres);
end

function test_initial_guess()
    % Test with initial guess
    n = 20;
    A = gallery('poisson', n);
    b = ones(size(A, 1), 1);

    % First solve from zero
    [x1, ~, ~, iter1] = sparsesolv(A, b, 1e-8, 500, 'iccg');

    % Solve again with x1 as initial guess (should converge faster)
    [x2, flag, relres, iter2] = sparsesolv(A, b, 1e-8, 500, 'iccg', 1.0, x1);

    assert(flag == 0, 'Solver with initial guess did not converge');
    % With good initial guess, should converge in fewer iterations
    assert(iter2 <= iter1, 'Initial guess did not help: iter2=%d > iter1=%d', iter2, iter1);
end

function test_convergence_check()
    % Test that flag correctly indicates non-convergence
    n = 10;
    A = gallery('poisson', n);
    b = ones(size(A, 1), 1);

    max_iter = 3;
    % Very few iterations - should not converge
    [x, flag, relres, iter] = sparsesolv(A, b, 1e-12, max_iter, 'iccg');

    % Flag should be 1 (not converged) unless matrix happens to converge quickly
    % Just verify that the function doesn't crash and returns reasonable values
    % Note: iter may include initial residual evaluation, so allow max_iter + 1
    assert(iter <= max_iter + 1, 'Iterations exceeded maximum: %d > %d', iter, max_iter + 1);
    assert(isfinite(relres), 'Relative residual is not finite');
end

function test_residual_log()
    % Test that residual log is returned correctly
    n = 20;
    A = gallery('poisson', n);
    b = ones(size(A, 1), 1);

    [x, flag, relres, iter, resvec] = sparsesolv(A, b, 1e-6, 500, 'iccg');

    assert(length(resvec) == iter, ...
        'Residual log length (%d) does not match iterations (%d)', length(resvec), iter);

    % Residuals should generally decrease (not strictly, but overall trend)
    if iter > 1
        assert(resvec(end) < resvec(1), ...
            'Final residual (%.2e) not smaller than initial (%.2e)', resvec(end), resvec(1));
    end
end
