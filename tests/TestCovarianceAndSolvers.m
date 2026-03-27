classdef TestCovarianceAndSolvers < matlab.unittest.TestCase
    methods (TestClassSetup)
        function addRoot(testCase)
            root = fileparts(fileparts(mfilename("fullpath")));
            addpath(root);
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(root));
        end
    end

    methods (Test)
        function testCovarianceModels(testCase)
            x = linspace(0, 10, 8);
            y = linspace(0, 4, 6);
            [X, Y] = meshgrid(x, y);
            p = [X(:), Y(:)];
            models = {'gaussian', 'exponential', 'spherical', 'matern'};
            for i = 1:numel(models)
                params = struct("variance", 2.0, "corrLength", [3 2], "smoothness", 1.1, "nugget", 1e-10);
                C = rf2d.CovarianceModel.pairwise(models{i}, params, p, p);
                testCase.verifySize(C, [numel(p(:, 1)), numel(p(:, 1))]);
                testCase.verifyGreaterThanOrEqual(min(diag(C)), 1.9);
            end
        end

        function testSolverOutputs(testCase)
            x = linspace(0, 20, 16);
            y = linspace(0, 10, 8);
            params = struct("variance", 1.5, "corrLength", [5 2], "smoothness", 1.2, "nugget", 1e-8);
            solvers = {'cholesky', 'kl', 'circulant', 'approximate'};
            for i = 1:numel(solvers)
                g = rf2d.createGenerator(x, y, "matern", params, "Solver", solvers{i}, "Seed", 17, "KLModes", 32, "MaxRank", 32);
                f = g.realize();
                testCase.verifySize(f, [numel(y), numel(x)]);
            end
        end
    end
end
