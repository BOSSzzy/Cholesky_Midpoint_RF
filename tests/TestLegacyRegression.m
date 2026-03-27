classdef TestLegacyRegression < matlab.unittest.TestCase
    methods (TestClassSetup)
        function addRoot(testCase)
            root = fileparts(fileparts(mfilename("fullpath")));
            addpath(root);
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(root));
        end
    end

    methods (Test)
        function testLegacyWrapperCallable(testCase)
            x = repmat((1:20)', 60, 1);
            y = kron((1:60)', ones(20, 1));
            tmp = [tempname '.csv'];
            writematrix([x y], tmp);
            [RFC, RFPHI, c, phi] = RandomField2DCholMethod(tmp, 1, 0);
            testCase.verifySize(RFC, [60, 20]);
            testCase.verifySize(RFPHI, [60, 20]);
            testCase.verifySize(c, [1200, 1]);
            testCase.verifySize(phi, [1200, 1]);
        end
    end
end
