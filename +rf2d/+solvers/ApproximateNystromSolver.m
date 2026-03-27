classdef ApproximateNystromSolver < rf2d.solvers.BaseSolver
    properties
        Points
        Model (1,1) string
        Parameters (1,1) struct
        Rank (1,1) double = 512
        Jitter (1,1) double = 1e-10
        Basis
    end

    methods
        function obj = ApproximateNystromSolver(points, model, parameters, meanValue, rankValue, jitter)
            obj.Points = points;
            obj.Model = lower(string(model));
            obj.Parameters = parameters;
            obj.NumPoints = size(points, 1);
            obj.Mean = meanValue;
            obj.Rank = min(rankValue, obj.NumPoints);
            obj.Jitter = jitter;
        end

        function prepare(obj)
            n = obj.NumPoints;
            idx = randperm(n, obj.Rank);
            landmarks = obj.Points(idx, :);
            Kmm = rf2d.CovarianceModel.pairwise(obj.Model, obj.Parameters, landmarks, landmarks);
            Knm = rf2d.CovarianceModel.pairwise(obj.Model, obj.Parameters, obj.Points, landmarks);
            Kmm = (Kmm + Kmm') / 2 + obj.Jitter * eye(obj.Rank);
            [L, p] = chol(Kmm, "lower");
            if p ~= 0
                error("rf2d:NumericalStability", "Nyström landmark matrix is singular.");
            end
            obj.Basis = Knm / L';
        end

        function fields = sample(obj, nFields, useGPU, useParallel)
            r = size(obj.Basis, 2);
            B = obj.Basis;
            if useGPU
                B = gpuArray(B);
                Z = randn(r, nFields, "gpuArray");
            else
                Z = randn(r, nFields);
            end
            if useParallel && ~useGPU && nFields > 1
                fields = zeros(obj.NumPoints, nFields);
                parfor i = 1:nFields
                    fields(:, i) = obj.Mean + B * Z(:, i);
                end
            else
                fields = obj.Mean + B * Z;
            end
            if useGPU
                fields = gather(fields);
            end
        end
    end
end
