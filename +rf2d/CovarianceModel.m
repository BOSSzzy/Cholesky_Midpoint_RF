classdef CovarianceModel
    methods (Static)
        function params = normalizeParameters(model, params)
            model = lower(char(string(model)));
            if ~isfield(params, "variance")
                params.variance = 1;
            end
            if ~isfield(params, "corrLength")
                error("rf2d:InvalidParameter", "corrLength is required.");
            end
            if numel(params.corrLength) == 1
                params.corrLength = [params.corrLength params.corrLength];
            end
            if ~isfield(params, "smoothness")
                params.smoothness = 0.5;
            end
            if ~isfield(params, "nugget")
                params.nugget = 0;
            end
            rf2d.Validation.mustBeSupportedModel(model);
            rf2d.Validation.mustBePositive(params.variance, "variance");
            rf2d.Validation.mustBePositive(params.corrLength(1), "corrLength(1)");
            rf2d.Validation.mustBePositive(params.corrLength(2), "corrLength(2)");
            if strcmp(model, "matern")
                rf2d.Validation.mustBePositive(params.smoothness, "smoothness");
            end
            rf2d.Validation.mustBeNonnegative(params.nugget, "nugget");
        end

        function C = pairwise(model, params, pointsA, pointsB)
            model = lower(char(string(model)));
            params = rf2d.CovarianceModel.normalizeParameters(model, params);
            xa = pointsA(:, 1);
            ya = pointsA(:, 2);
            xb = pointsB(:, 1)';
            yb = pointsB(:, 2)';
            dx = abs(xa - xb) ./ params.corrLength(1);
            dy = abs(ya - yb) ./ params.corrLength(2);
            h = sqrt(dx.^2 + dy.^2);
            C = rf2d.CovarianceModel.core(model, h, params);
            if isequal(size(pointsA), size(pointsB)) && all(pointsA(:) == pointsB(:))
                C = C + params.nugget * eye(size(C), "like", C);
            end
        end

        function c = lagged(model, params, dx, dy)
            model = lower(char(string(model)));
            params = rf2d.CovarianceModel.normalizeParameters(model, params);
            hx = abs(dx) ./ params.corrLength(1);
            hy = abs(dy) ./ params.corrLength(2);
            h = sqrt(hx.^2 + hy.^2);
            c = rf2d.CovarianceModel.core(model, h, params);
        end
    end

    methods (Static, Access = private)
        function C = core(model, h, params)
            if strcmp(model, "gaussian")
                kernel = exp(-(h.^2));
            elseif strcmp(model, "exponential")
                kernel = exp(-h);
            elseif strcmp(model, "spherical")
                kernel = zeros(size(h), "like", h);
                m = h <= 1;
                hm = h(m);
                kernel(m) = 1 - 1.5 * hm + 0.5 * hm.^3;
            elseif strcmp(model, "matern")
                nu = params.smoothness;
                z = sqrt(2 * nu) * max(h, eps);
                pref = (2^(1 - nu)) / gamma(nu);
                kernel = pref * (z.^nu) .* besselk(nu, z);
                kernel(h == 0) = 1;
            else
                error("rf2d:InvalidModel", "Unsupported covariance model.");
            end
            C = params.variance * kernel;
        end
    end
end
