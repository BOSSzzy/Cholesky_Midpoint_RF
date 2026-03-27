classdef Validation
    methods (Static)
        function mustBeSupportedModel(model)
            valid = {'gaussian', 'exponential', 'spherical', 'matern'};
            if ~any(strcmpi(char(string(model)), valid))
                error("rf2d:InvalidModel", "Model must be Gaussian, Exponential, Spherical, or Matern.");
            end
        end

        function mustBeSupportedSolver(solver)
            valid = {'auto', 'cholesky', 'kl', 'circulant', 'approximate'};
            if ~any(strcmpi(char(string(solver)), valid))
                error("rf2d:InvalidSolver", "Unsupported solver.");
            end
        end

        function mustBePositive(value, name)
            if any(value <= 0)
                error("rf2d:InvalidParameter", "%s must be positive.", name);
            end
        end

        function mustBeNonnegative(value, name)
            if any(value < 0)
                error("rf2d:InvalidParameter", "%s must be nonnegative.", name);
            end
        end

        function mustBeRegularGrid(x, y)
            dx = diff(x);
            dy = diff(y);
            if any(abs(dx - mean(dx)) > 1e-10 * max(1, abs(mean(dx))))
                error("rf2d:Grid", "x must be regular for circulant embedding.");
            end
            if any(abs(dy - mean(dy)) > 1e-10 * max(1, abs(mean(dy))))
                error("rf2d:Grid", "y must be regular for circulant embedding.");
            end
        end
    end
end
