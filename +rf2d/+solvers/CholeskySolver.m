classdef CholeskySolver < rf2d.solvers.BaseSolver
    properties
        CovarianceMatrix
        Jitter (1,1) double = 1e-10
        Factor
    end

    methods
        function obj = CholeskySolver(C, meanValue, jitter)
            obj.CovarianceMatrix = C;
            obj.NumPoints = size(C, 1);
            obj.Mean = meanValue;
            obj.Jitter = jitter;
        end

        function prepare(obj)
            C = (obj.CovarianceMatrix + obj.CovarianceMatrix') / 2;
            n = size(C, 1);
            jit = obj.Jitter;
            success = false;
            for k = 1:8
                [L, p] = chol(C + jit * speye(n), "lower");
                if p == 0
                    obj.Factor = L;
                    success = true;
                    break
                end
                jit = jit * 10;
            end
            if ~success
                error("rf2d:NumericalStability", "Cholesky factorization failed.");
            end
        end

        function fields = sample(obj, nFields, useGPU, useParallel)
            Z = randn(obj.NumPoints, nFields);
            if useGPU
                Z = gpuArray(Z);
                F = gpuArray(obj.Factor);
            else
                F = obj.Factor;
            end
            if useParallel && ~useGPU && nFields > 1
                fields = zeros(obj.NumPoints, nFields);
                parfor i = 1:nFields
                    fields(:, i) = obj.Mean + F * Z(:, i);
                end
            else
                fields = obj.Mean + F * Z;
            end
            if useGPU
                fields = gather(fields);
            end
        end
    end
end
