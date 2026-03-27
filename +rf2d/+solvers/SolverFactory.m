classdef SolverFactory
    methods (Static)
        function [solver, solverName] = create(x, y, points, model, params, options)
            n = size(points, 1);
            requested = lower(char(string(options.Solver)));
            rf2d.Validation.mustBeSupportedSolver(requested);
            solverName = requested;
            if strcmp(requested, "auto")
                if n <= options.MaxDirectPoints
                    solverName = 'cholesky';
                elseif rf2d.solvers.SolverFactory.isRegularGrid(x, y)
                    solverName = 'circulant';
                else
                    solverName = 'approximate';
                end
            end

            if strcmp(solverName, "cholesky")
                if n > options.MaxDirectPoints
                    solver = rf2d.solvers.ApproximateNystromSolver(points, model, params, options.Mean, options.MaxRank, options.Jitter);
                    solverName = 'approximate';
                else
                    C = rf2d.CovarianceModel.pairwise(model, params, points, points);
                    solver = rf2d.solvers.CholeskySolver(C, options.Mean, options.Jitter);
                end
            elseif strcmp(solverName, "kl")
                if n > options.MaxDirectPoints
                    solver = rf2d.solvers.ApproximateNystromSolver(points, model, params, options.Mean, options.MaxRank, options.Jitter);
                    solverName = 'approximate';
                else
                    C = rf2d.CovarianceModel.pairwise(model, params, points, points);
                    solver = rf2d.solvers.KLSolver(C, options.Mean, options.KLModes);
                end
            elseif strcmp(solverName, "circulant")
                solver = rf2d.solvers.CirculantEmbeddingSolver(x, y, model, params, options.Mean);
            elseif strcmp(solverName, "approximate")
                solver = rf2d.solvers.ApproximateNystromSolver(points, model, params, options.Mean, options.MaxRank, options.Jitter);
            else
                error("rf2d:InvalidSolver", "Unsupported solver.");
            end
        end
    end

    methods (Static, Access = private)
        function tf = isRegularGrid(x, y)
            tf = true;
            if numel(x) > 2
                dx = diff(x);
                tf = tf && max(abs(dx - mean(dx))) < 1e-10 * max(1, abs(mean(dx)));
            end
            if numel(y) > 2
                dy = diff(y);
                tf = tf && max(abs(dy - mean(dy))) < 1e-10 * max(1, abs(mean(dy)));
            end
        end
    end
end
