classdef ConditionalRandomField
    methods (Static)
        function result = run(config)
            if nargin < 1
                config = struct();
            end
            config = rf2d.ConditionalRandomField.normalizeConfig(config);
            rf2d.ConditionalRandomField.log(config, "info", "Starting CRF generation.");
            try
                if config.io.writeOutput && ~exist(config.io.dataDir, "dir")
                    mkdir(config.io.dataDir);
                end
                rng(config.sim.randomSeed, "twister");
                gridTable = readtable(config.io.gridTablePath);
                obsTable = readtable(config.io.observationPath);
                result = rf2d.ConditionalRandomField.fromTables(gridTable, obsTable, config);
                if config.io.writeOutput
                    writetable(result.unconditional.unconditionalField, config.io.unconditionalPath);
                    writetable(result.conditional.conditionalField, config.io.conditionalPath);
                end
            catch ME
                rf2d.ConditionalRandomField.log(config, "error", ME.message);
                rethrow(ME)
            end
            rf2d.ConditionalRandomField.log(config, "info", "Finished CRF generation.");
        end

        function result = fromTables(gridTable, obsTable, config)
            if nargin < 3
                config = struct();
            end
            config = rf2d.ConditionalRandomField.normalizeConfig(config);
            rf2d.ConditionalRandomField.validateGridTable(gridTable);
            rf2d.ConditionalRandomField.validateObservationTable(obsTable);
            rf2d.Validation.mustBePositive(config.prior.std, "prior.std");
            rf2d.Validation.mustBePositive(config.sim.nRealizations, "sim.nRealizations");
            rf2d.Validation.mustBeNonnegative(obsTable.obsVar, "obsVar");

            uncond = rf2d.ConditionalRandomField.runUnconditional(config, gridTable);
            cond = rf2d.ConditionalRandomField.runConditional(config, uncond.unconditionalField, obsTable);
            result = struct("config", config, "unconditional", uncond, "conditional", cond);
        end
    end

    methods (Static, Access = private)
        function uncond = runUnconditional(config, gridTable)
            coordMat = rf2d.ConditionalRandomField.gridToCoordinates(gridTable, config.grid.dx, config.grid.dy);
            covParams = rf2d.ConditionalRandomField.toCovarianceParameters(config);
            covMat = rf2d.CovarianceModel.pairwise(config.covModel.type, covParams, coordMat, coordMat);
            jitter = config.sim.jitter * eye(size(covMat, 1));
            try
                cholLower = chol(covMat + jitter, "lower");
            catch ME
                error("rf2d:CRFCholeskyFailed", "CRF covariance factorization failed: %s", ME.message);
            end

            nRealizations = round(config.sim.nRealizations);
            simStdNormal = randn(size(coordMat, 1), nRealizations);
            simGaussian = config.prior.mean + cholLower * simStdNormal;

            varNames = compose("RF_%04d", 1:nRealizations);
            simTable = array2table(simGaussian, "VariableNames", cellstr(varNames));
            uncondTable = [gridTable(:, ["id", "j", "k"]), simTable];
            uncondTable = sortrows(uncondTable, "id", "ascend");

            uncond = struct();
            uncond.coordinate = coordMat;
            uncond.covariance = covMat;
            uncond.choleskyLower = cholLower;
            uncond.unconditionalField = uncondTable;
        end

        function cond = runConditional(config, uncondTable, obsTable)
            simVarMask = startsWith(string(uncondTable.Properties.VariableNames), "RF_");
            simVarNames = uncondTable.Properties.VariableNames(simVarMask);
            priorField = table2array(uncondTable(:, simVarNames));

            gridCoord = rf2d.ConditionalRandomField.gridToCoordinates(uncondTable, config.grid.dx, config.grid.dy);
            idxObs = rf2d.ConditionalRandomField.locateObservationIndices(uncondTable, obsTable);
            obsCoord = gridCoord(idxObs, :);
            obsValue = obsTable.obsValue;
            obsVar = obsTable.obsVar;

            covParams = rf2d.ConditionalRandomField.toCovarianceParameters(config);
            covGO = rf2d.CovarianceModel.pairwise(config.covModel.type, covParams, gridCoord, obsCoord);
            covOO = rf2d.CovarianceModel.pairwise(config.covModel.type, covParams, obsCoord, obsCoord);
            gainMat = covGO / (covOO + diag(obsVar));

            priorObs = priorField(idxObs, :);
            innovation = obsValue - priorObs;
            correction = gainMat * innovation;
            condField = priorField + correction;

            condTable = [uncondTable(:, ["id", "j", "k"]), array2table(condField, "VariableNames", simVarNames)];
            condTable = sortrows(condTable, "id", "ascend");

            cond = struct();
            cond.observationIndex = idxObs;
            cond.krigingGain = gainMat;
            cond.conditionalField = condTable;
            cond.conditionalMean = mean(condField, 2);
            cond.conditionalStd = std(condField, 0, 2);
        end

        function coordMat = gridToCoordinates(gridTable, dx, dy)
            coordX = (gridTable.j - 0.5) * dx;
            coordY = (gridTable.k - 0.5) * dy;
            coordMat = [coordX, coordY];
        end

        function idxObs = locateObservationIndices(gridTable, obsTable)
            obsKey = strcat(string(obsTable.j), "_", string(obsTable.k));
            gridKey = strcat(string(gridTable.j), "_", string(gridTable.k));
            [isFound, idxObs] = ismember(obsKey, gridKey);
            if any(~isFound)
                error("rf2d:CRFObservationNotFound", "Observation (j,k) is not found in grid table.");
            end
        end

        function params = toCovarianceParameters(config)
            params = struct();
            params.variance = config.prior.std^2;
            params.corrLength = [config.covModel.rangeX, config.covModel.rangeY];
            params.smoothness = config.covModel.smoothness;
            params.nugget = config.covModel.nugget;
        end

        function validateGridTable(gridTable)
            required = ["id", "j", "k"];
            names = string(gridTable.Properties.VariableNames);
            if ~all(ismember(required, names))
                error("rf2d:CRFInvalidGridTable", "Grid table must include columns: id, j, k.");
            end
        end

        function validateObservationTable(obsTable)
            required = ["j", "k", "obsValue", "obsVar"];
            names = string(obsTable.Properties.VariableNames);
            if ~all(ismember(required, names))
                error("rf2d:CRFInvalidObservationTable", "Observation table must include columns: j, k, obsValue, obsVar.");
            end
        end

        function config = normalizeConfig(config)
            if ~isfield(config, "grid"), config.grid = struct(); end
            if ~isfield(config.grid, "dx"), config.grid.dx = 1.0; end
            if ~isfield(config.grid, "dy"), config.grid.dy = 1.0; end

            if ~isfield(config, "covModel"), config.covModel = struct(); end
            if ~isfield(config.covModel, "type"), config.covModel.type = "exponential"; end
            if ~isfield(config.covModel, "nugget"), config.covModel.nugget = 0.0; end
            if ~isfield(config.covModel, "rangeX"), config.covModel.rangeX = 20.0; end
            if ~isfield(config.covModel, "rangeY"), config.covModel.rangeY = 10.0; end
            if ~isfield(config.covModel, "smoothness"), config.covModel.smoothness = 1.5; end
            rf2d.Validation.mustBeSupportedModel(config.covModel.type);
            rf2d.Validation.mustBeNonnegative(config.covModel.nugget, "covModel.nugget");
            rf2d.Validation.mustBePositive(config.covModel.rangeX, "covModel.rangeX");
            rf2d.Validation.mustBePositive(config.covModel.rangeY, "covModel.rangeY");

            if ~isfield(config, "prior"), config.prior = struct(); end
            if ~isfield(config.prior, "mean"), config.prior.mean = 10.0; end
            if ~isfield(config.prior, "std"), config.prior.std = 2.0; end

            if ~isfield(config, "sim"), config.sim = struct(); end
            if ~isfield(config.sim, "nRealizations"), config.sim.nRealizations = 300; end
            if ~isfield(config.sim, "randomSeed"), config.sim.randomSeed = 20260329; end
            if ~isfield(config.sim, "jitter"), config.sim.jitter = 1e-10; end
            rf2d.Validation.mustBePositive(config.sim.jitter, "sim.jitter");

            if ~isfield(config, "io"), config.io = struct(); end
            if ~isfield(config.io, "dataDir"), config.io.dataDir = ".\data"; end
            if ~isfield(config.io, "gridTablePath"), config.io.gridTablePath = ".\data\grid.csv"; end
            if ~isfield(config.io, "observationPath"), config.io.observationPath = ".\data\observation.csv"; end
            if ~isfield(config.io, "unconditionalPath"), config.io.unconditionalPath = ".\data\unconditional_field.csv"; end
            if ~isfield(config.io, "conditionalPath"), config.io.conditionalPath = ".\data\conditional_field.csv"; end
            if ~isfield(config.io, "writeOutput"), config.io.writeOutput = true; end

            if ~isfield(config, "logging"), config.logging = struct(); end
            if ~isfield(config.logging, "enabled"), config.logging.enabled = false; end
            if ~isfield(config.logging, "level"), config.logging.level = "info"; end
        end

        function log(config, level, message)
            if ~isfield(config, "logging") || ~isfield(config.logging, "enabled") || ~config.logging.enabled
                return
            end
            levels = ["debug", "info", "warn", "error"];
            req = find(levels == lower(string(config.logging.level)), 1);
            cur = find(levels == lower(string(level)), 1);
            if isempty(req), req = 2; end
            if isempty(cur), cur = 2; end
            if cur >= req
                fprintf("[rf2d.crf][%s] %s\n", upper(char(string(level))), char(string(message)));
            end
        end
    end
end
