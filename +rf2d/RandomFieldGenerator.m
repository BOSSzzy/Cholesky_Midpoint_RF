classdef RandomFieldGenerator < handle
    properties (SetAccess = private)
        X
        Y
        GridX
        GridY
        Points
        CovarianceModelName (1,1) string
        CovarianceParameters (1,1) struct
        Options (1,1) struct
        Solver
        SolverName (1,1) string
        Stream
    end

    methods
        function obj = RandomFieldGenerator(x, y, covarianceModel, covarianceParameters, options)
            obj.X = x(:)';
            obj.Y = y(:)';
            [obj.GridX, obj.GridY] = meshgrid(obj.X, obj.Y);
            obj.Points = [obj.GridX(:), obj.GridY(:)];
            obj.CovarianceModelName = lower(string(covarianceModel));
            obj.CovarianceParameters = rf2d.CovarianceModel.normalizeParameters(obj.CovarianceModelName, covarianceParameters);
            obj.Options = options;
            obj.Stream = RandStream("mt19937ar", "Seed", options.Seed);
            [obj.Solver, obj.SolverName] = rf2d.solvers.SolverFactory.create(obj.X, obj.Y, obj.Points, obj.CovarianceModelName, obj.CovarianceParameters, options);
            obj.prepareSolverWithFallback();
        end

        function reseed(obj, seed)
            obj.Options.Seed = seed;
            obj.Stream = RandStream("mt19937ar", "Seed", seed);
        end

        function field = realize(obj, options)
            arguments
                obj
                options.UseGPU (1,1) logical = obj.Options.UseGPU
                options.UseParallel (1,1) logical = false
            end
            fields = obj.generateBatch(1, options.UseGPU, options.UseParallel);
            field = reshape(fields(:, 1), numel(obj.Y), numel(obj.X));
        end

        function fields = generateBatch(obj, nFields, useGPU, useParallel)
            arguments
                obj
                nFields (1,1) double {mustBeInteger, mustBePositive}
                useGPU (1,1) logical = obj.Options.UseGPU
                useParallel (1,1) logical = false
            end
            oldStream = RandStream.getGlobalStream;
            cleanupObj = onCleanup(@() RandStream.setGlobalStream(oldStream));
            RandStream.setGlobalStream(obj.Stream);
            if useParallel
                fields = rf2d.RandomFieldGenerator.sampleWithSPMD(obj.Solver, nFields, useGPU);
                if isempty(fields)
                    fields = obj.Solver.sample(nFields, useGPU, true);
                end
            else
                fields = obj.Solver.sample(nFields, useGPU, false);
            end
            obj.Stream = RandStream.getGlobalStream;
            clear cleanupObj
        end

        function contourPlot(obj, field, levels)
            if nargin < 3
                levels = 20;
            end
            rf2d.viz.FieldPlotter.contour(obj.X, obj.Y, field, levels);
        end

        function surfacePlot(obj, field)
            rf2d.viz.FieldPlotter.surface(obj.X, obj.Y, field);
        end

        function histogramPlot(~, field, nBins)
            if nargin < 3
                nBins = 40;
            end
            rf2d.viz.FieldPlotter.histogram(field, nBins);
        end

        function variogramCloud(obj, field, sampleCount)
            if nargin < 3
                sampleCount = 5000;
            end
            rf2d.viz.FieldPlotter.variogramCloud(obj.Points, field(:), sampleCount);
        end

        function exportCSV(obj, field, filePath)
            rf2d.io.FieldExporter.toCSV(obj.GridX, obj.GridY, field, filePath);
        end

        function exportVTK(obj, field, filePath)
            rf2d.io.FieldExporter.toVTK(obj.GridX, obj.GridY, field, filePath);
        end

        function exportMAT(obj, field, filePath, variableName)
            if nargin < 5
                variableName = "field";
            end
            rf2d.io.FieldExporter.toMAT(obj.GridX, obj.GridY, field, filePath, variableName);
        end
    end

    methods (Access = private)
        function prepareSolverWithFallback(obj)
            try
                obj.Solver.prepare();
            catch ME
                if ~strcmp(char(obj.SolverName), "approximate")
                    obj.Solver = rf2d.solvers.ApproximateNystromSolver(obj.Points, obj.CovarianceModelName, obj.CovarianceParameters, obj.Options.Mean, obj.Options.MaxRank, obj.Options.Jitter);
                    obj.SolverName = "approximate";
                    obj.Solver.prepare();
                else
                    rethrow(ME)
                end
            end
        end
    end

    methods (Static, Access = private)
        function fields = sampleWithSPMD(solver, nFields, useGPU)
            fields = [];
            pool = gcp("nocreate");
            if isempty(pool) || pool.NumWorkers < 2
                return
            end
            spmd
                localCount = floor(nFields / numlabs) + (labindex <= mod(nFields, numlabs));
                localFields = solver.sample(localCount, useGPU, false);
            end
            chunks = cell(1, numel(localFields));
            for k = 1:numel(localFields)
                chunks{k} = localFields{k};
            end
            fields = cat(2, chunks{:});
        end
    end
end
