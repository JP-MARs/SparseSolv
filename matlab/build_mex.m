function build_mex(varargin)
% BUILD_MEX Build the SparseSolv MEX interface
%
% SYNTAX:
%   build_mex           - Build with default options
%   build_mex('debug')  - Build with debug flags
%   build_mex('clean')  - Clean build files
%   build_mex('verbose') - Build with verbose output
%
% REQUIREMENTS:
%   - MATLAB with supported C++ compiler
%   - Windows: Visual Studio 2019 or later recommended
%   - Linux/Mac: GCC 7+ or Clang 7+
%
% The MEX file will be created in the current directory.

    % Parse options
    debug_mode = any(strcmpi(varargin, 'debug'));
    clean_mode = any(strcmpi(varargin, 'clean'));
    verbose_mode = any(strcmpi(varargin, 'verbose'));

    % Get paths
    this_dir = fileparts(mfilename('fullpath'));
    parent_dir = fileparts(this_dir);
    sparse_solv_dir = fullfile(parent_dir, 'SparseSolv');

    % Clean build
    if clean_mode
        fprintf('Cleaning build files...\n');
        delete(fullfile(this_dir, '*.obj'));
        delete(fullfile(this_dir, '*.o'));
        delete(fullfile(this_dir, ['sparsesolv_mex.' mexext]));
        fprintf('Clean complete.\n');
        return;
    end

    % Check if SparseSolv directory exists
    if ~isfolder(sparse_solv_dir)
        error('build_mex:notFound', ...
            'SparseSolv directory not found at: %s', sparse_solv_dir);
    end

    fprintf('========================================\n');
    fprintf('Building SparseSolv MEX Interface\n');
    fprintf('========================================\n');
    fprintf('SparseSolv directory: %s\n', sparse_solv_dir);
    fprintf('Output directory: %s\n', this_dir);
    fprintf('\n');

    % Source files to compile
    mex_source = fullfile(this_dir, 'sparsesolv_mex.cpp');

    % SparseSolv source files needed for the MEX interface
    sparse_solv_sources = {
        'SparseMat.cpp'
        'SparseMatC.cpp'
        'SparseMatOperators.cpp'
        'MatSolvers_Base.cpp'
        'MatSolvers_ICCG.cpp'
        'MatSolvers_ICCOCG.cpp'
        'MatSolvers_ICMRTR.cpp'
        'MatSolvers_SGSMRTR.cpp'
    };

    % Build full paths
    source_files = cell(1, length(sparse_solv_sources) + 1);
    source_files{1} = mex_source;
    for i = 1:length(sparse_solv_sources)
        source_files{i+1} = fullfile(sparse_solv_dir, sparse_solv_sources{i});

        % Verify file exists
        if ~isfile(source_files{i+1})
            error('build_mex:sourceNotFound', ...
                'Source file not found: %s', source_files{i+1});
        end
    end

    fprintf('Source files:\n');
    for i = 1:length(source_files)
        fprintf('  [%d] %s\n', i, source_files{i});
    end
    fprintf('\n');

    % Compiler flags
    include_flag = ['-I' sparse_solv_dir];

    % Common flags
    cxx_flags = {};

    % Platform-specific flags
    if ispc
        % Windows (MSVC)
        % /openmp for OpenMP support, /EHsc for exception handling
        cxx_flags = [cxx_flags, {'COMPFLAGS="$COMPFLAGS /std:c++17 /EHsc /W3 /openmp"'}];
        if debug_mode
            cxx_flags = [cxx_flags, {'COMPFLAGS="$COMPFLAGS /Od /Zi"'}];
            cxx_flags = [cxx_flags, {'-g'}];
        else
            cxx_flags = [cxx_flags, {'COMPFLAGS="$COMPFLAGS /O2"'}];
        end
    else
        % Linux/Mac (GCC/Clang)
        % -fopenmp for OpenMP support
        cxx_flags = [cxx_flags, {'CXXFLAGS="$CXXFLAGS -std=c++17 -fopenmp"'}];
        cxx_flags = [cxx_flags, {'LDFLAGS="$LDFLAGS -fopenmp"'}];
        if debug_mode
            cxx_flags = [cxx_flags, {'CXXFLAGS="$CXXFLAGS -O0 -g"'}];
        else
            cxx_flags = [cxx_flags, {'CXXFLAGS="$CXXFLAGS -O3"'}];
        end
    end

    % Build MEX arguments
    mex_args = {};

    if verbose_mode
        mex_args = [mex_args, {'-v'}];
    end

    % Use R2017b API for complex number compatibility
    % This allows mxGetPr/mxGetPi to work in newer MATLAB versions
    mex_args = [mex_args, {'-R2017b'}];

    % Output file
    output_file = fullfile(this_dir, 'sparsesolv_mex');
    mex_args = [mex_args, {'-output', output_file}];

    % Include paths
    mex_args = [mex_args, {include_flag}];
    % Add parent directory for Eigen (000_thirdparty/Eigen)
    eigen_include_flag = ['-I' parent_dir];
    mex_args = [mex_args, {eigen_include_flag}];

    % Compiler flags
    mex_args = [mex_args, cxx_flags];

    % Source files
    mex_args = [mex_args, source_files];

    % Display build command (for debugging)
    fprintf('MEX build arguments:\n');
    for i = 1:length(mex_args)
        fprintf('  %s\n', mex_args{i});
    end
    fprintf('\n');

    % Build
    fprintf('Compiling... (this may take a moment)\n');
    try
        tic;
        mex(mex_args{:});
        elapsed = toc;
        fprintf('\n');
        fprintf('========================================\n');
        fprintf('Build successful!\n');
        fprintf('========================================\n');
        fprintf('Output: %s.%s\n', output_file, mexext);
        fprintf('Build time: %.1f seconds\n', elapsed);
        fprintf('\n');
        fprintf('To test, run:\n');
        fprintf('  test_sparsesolv\n');
    catch ME
        fprintf(2, '\n');
        fprintf(2, '========================================\n');
        fprintf(2, 'Build FAILED!\n');
        fprintf(2, '========================================\n');
        fprintf(2, 'Error: %s\n', ME.message);
        fprintf(2, '\n');
        fprintf(2, 'Troubleshooting:\n');
        fprintf(2, '  1. Ensure a supported C++ compiler is installed\n');
        fprintf(2, '     Run: mex -setup C++\n');
        fprintf(2, '  2. Check that all source files exist\n');
        fprintf(2, '  3. Try verbose mode: build_mex(''verbose'')\n');
        rethrow(ME);
    end
end
