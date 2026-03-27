classdef KLSolver < rf2d.solvers.BaseSolver
    properties
        CovarianceMatrix
        NumModes (1,1) double = 256
        EigenVectors
        EigenValues
    end

    methods
        function obj = KLSolver(C, meanValue, numModes)
            obj.CovarianceMatrix = C;
            obj.NumPoints = size(C, 1);
            obj.Mean = meanValue;
            obj.NumModes = min(numModes, obj.NumPoints);
        end

        function prepare(obj)
            C = (obj.CovarianceMatrix + obj.CovarianceMatrix') / 2;
            [V, D] = eigs(C, obj.NumModes, "largestreal", "Tolerance", 1e-7);
            d = max(real(diag(D)), 0);
            keep = d > 0;
            obj.EigenVectors = V(:, keep);
            obj.EigenValues = d(keep);
        end

        function fields = sample(obj, nFields, useGPU, useParallel)
            m = numel(obj.EigenValues);
            Z = randn(m, nFields);
            S = sqrt(obj.EigenValues);
            if useGPU
                Z = gpuArray(Z);
                V = gpuArray(obj.EigenVectors);
                S = gpuArray(S);
            else
                V = obj.EigenVectors;
            end
            if useParallel && ~useGPU && nFields > 1
                fields = zeros(obj.NumPoints, nFields);
                parfor i = 1:nFields
                    fields(:, i) = obj.Mean + V * (S .* Z(:, i));
                end
            else
                fields = obj.Mean + V * (S .* Z);
            end
            if useGPU
                fields = gather(fields);
            end
        end
    end
end
