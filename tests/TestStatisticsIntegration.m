classdef TestStatisticsIntegration < matlab.unittest.TestCase
    methods (TestClassSetup)
        function addRoot(testCase)
            root = fileparts(fileparts(mfilename("fullpath")));
            addpath(root);
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(root));
        end
    end

    methods (Test)
        function testEnsembleStatistics(testCase)
            x = linspace(0, 15, 10);
            y = linspace(0, 12, 10);
            params = struct("variance", 1.0, "corrLength", [4 3], "smoothness", 1.5, "nugget", 1e-10);
            g = rf2d.createGenerator(x, y, "gaussian", params, "Solver", "cholesky", "Seed", 123);
            fields = g.generateBatch(20000, false, false);
            empiricalMean = mean(fields(:));
            empiricalVar = var(fields(:), 1);
            testCase.verifyLessThanOrEqual(abs(empiricalMean - 0), 0.01);
            testCase.verifyLessThanOrEqual(abs(empiricalVar - params.variance) / params.variance, 0.01);

            [GX, GY] = meshgrid(x, y);
            pts = [GX(:), GY(:)];
            idxA = 11;
            idxB = 11;
            theor = rf2d.CovarianceModel.pairwise("gaussian", params, pts(idxA, :), pts(idxB, :));
            covEmp = mean(fields(idxA, :) .* fields(idxB, :)) - mean(fields(idxA, :)) * mean(fields(idxB, :));
            testCase.verifyLessThanOrEqual(abs(covEmp - theor) / max(abs(theor), eps), 0.011);
        end

        function testExportUtilities(testCase)
            x = 0:1:15;
            y = 0:1:9;
            params = struct("variance", 1.2, "corrLength", [3 2], "smoothness", 1.0, "nugget", 1e-8);
            g = rf2d.createGenerator(x, y, "exponential", params, "Solver", "circulant", "Seed", 33);
            f = g.realize();
            tmp = tempname;
            csvFile = [tmp '.csv'];
            vtkFile = [tmp '.vtk'];
            matFile = [tmp '.mat'];
            g.exportCSV(f, csvFile);
            g.exportVTK(f, vtkFile);
            g.exportMAT(f, matFile, "fieldA");
            testCase.verifyTrue(isfile(csvFile));
            testCase.verifyTrue(isfile(vtkFile));
            testCase.verifyTrue(isfile(matFile));
        end
    end
end
