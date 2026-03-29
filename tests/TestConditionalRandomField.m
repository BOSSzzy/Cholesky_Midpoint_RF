classdef TestConditionalRandomField < matlab.unittest.TestCase
    methods (TestClassSetup)
        function addRoot(testCase)
            root = fileparts(fileparts(mfilename("fullpath")));
            addpath(root);
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(root));
        end
    end

    methods (Test)
        function testFromTables(testCase)
            [gridTable, obsTable] = TestConditionalRandomField.buildInputTables();
            config = struct();
            config.prior.mean = 5.0;
            config.prior.std = 1.4;
            config.covModel.type = "exponential";
            config.covModel.rangeX = 4.0;
            config.covModel.rangeY = 3.0;
            config.covModel.nugget = 0;
            config.sim.nRealizations = 16;
            config.sim.randomSeed = 77;
            config.io.writeOutput = false;

            result = rf2d.ConditionalRandomField.fromTables(gridTable, obsTable, config);
            condTable = result.conditional.conditionalField;
            rfVars = startsWith(string(condTable.Properties.VariableNames), "RF_");
            rfData = table2array(condTable(:, rfVars));

            testCase.verifySize(result.unconditional.unconditionalField, [height(gridTable), 3 + config.sim.nRealizations]);
            testCase.verifySize(condTable, [height(gridTable), 3 + config.sim.nRealizations]);
            testCase.verifySize(result.conditional.krigingGain, [height(gridTable), height(obsTable)]);
            testCase.verifySize(rfData, [height(gridTable), config.sim.nRealizations]);

            obsIdx = result.conditional.observationIndex;
            obsDelta = abs(rfData(obsIdx, :) - obsTable.obsValue);
            testCase.verifyLessThanOrEqual(max(obsDelta(:)), 1e-5);
        end

        function testRunWithFileIO(testCase)
            [gridTable, obsTable] = TestConditionalRandomField.buildInputTables();
            tmpDir = tempname;
            mkdir(tmpDir);
            gridPath = fullfile(tmpDir, "grid.csv");
            obsPath = fullfile(tmpDir, "observation.csv");
            uncondPath = fullfile(tmpDir, "unconditional_field.csv");
            condPath = fullfile(tmpDir, "conditional_field.csv");
            writetable(gridTable, gridPath);
            writetable(obsTable, obsPath);

            config = struct();
            config.grid.dx = 1.0;
            config.grid.dy = 1.0;
            config.covModel.type = "gaussian";
            config.covModel.rangeX = 5.0;
            config.covModel.rangeY = 5.0;
            config.prior.mean = 0.0;
            config.prior.std = 1.0;
            config.sim.nRealizations = 8;
            config.sim.randomSeed = 19;
            config.io.dataDir = tmpDir;
            config.io.gridTablePath = gridPath;
            config.io.observationPath = obsPath;
            config.io.unconditionalPath = uncondPath;
            config.io.conditionalPath = condPath;
            config.logging.enabled = true;
            config.logging.level = "error";

            result = rf2d.runConditionalRandomField(config);
            testCase.verifyTrue(isfile(uncondPath));
            testCase.verifyTrue(isfile(condPath));
            testCase.verifySize(result.conditional.conditionalField, [height(gridTable), 3 + config.sim.nRealizations]);

            wrapper = conrandex(config);
            testCase.verifySize(wrapper.conditionalField, [height(gridTable), 3 + config.sim.nRealizations]);
        end
    end

    methods (Static, Access = private)
        function [gridTable, obsTable] = buildInputTables()
            [jGrid, kGrid] = meshgrid(1:5, 1:4);
            id = (1:numel(jGrid))';
            gridTable = table(id, jGrid(:), kGrid(:), 'VariableNames', {'id', 'j', 'k'});

            obsTable = table([1; 3; 5], [1; 2; 4], [4.5; 5.2; 6.0], [1e-8; 1e-8; 1e-8], ...
                'VariableNames', {'j', 'k', 'obsValue', 'obsVar'});
        end
    end
end
