classdef TestApiCoverage < matlab.unittest.TestCase
    methods (TestClassSetup)
        function addRoot(testCase)
            root = fileparts(fileparts(mfilename("fullpath")));
            addpath(root);
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(root));
        end
    end

    methods (Test)
        function testGeneratorMethods(testCase)
            x = linspace(0, 20, 16);
            y = linspace(0, 10, 8);
            params = struct("variance", 1.0, "corrLength", [5 3], "smoothness", 1.2, "nugget", 1e-8);
            g = rf2d.createGenerator(x, y, "matern", params, "Solver", "auto", "Seed", 11);
            g.reseed(99);
            f = g.realize();
            batch = g.generateBatch(4, false, false);
            testCase.verifySize(f, [numel(y), numel(x)]);
            testCase.verifySize(batch, [numel(x) * numel(y), 4]);

            fig1 = figure("Visible", "off"); g.contourPlot(f, 12); close(fig1);
            fig2 = figure("Visible", "off"); g.surfacePlot(f); close(fig2);
            fig3 = figure("Visible", "off"); g.histogramPlot(f, 20); close(fig3);
            fig4 = figure("Visible", "off"); g.variogramCloud(f, 200); close(fig4);
        end

        function testValidationBranches(testCase)
            testCase.verifyError(@() rf2d.Validation.mustBeSupportedModel("abc"), "rf2d:InvalidModel");
            testCase.verifyError(@() rf2d.Validation.mustBeSupportedSolver("abc"), "rf2d:InvalidSolver");
            testCase.verifyError(@() rf2d.Validation.mustBePositive(0, "v"), "rf2d:InvalidParameter");
            testCase.verifyError(@() rf2d.Validation.mustBeNonnegative(-1, "v"), "rf2d:InvalidParameter");
            testCase.verifyError(@() rf2d.Validation.mustBeRegularGrid([0 1 3], [0 1 2]), "rf2d:Grid");
        end
    end
end
