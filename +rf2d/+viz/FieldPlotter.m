classdef FieldPlotter
    methods (Static)
        function contour(x, y, field, levels)
            f = rf2d.viz.FieldPlotter.reshapeField(field, x, y);
            contourf(x, y, f, levels, "LineColor", "none");
            axis tight
            colorbar
        end

        function surface(x, y, field)
            f = rf2d.viz.FieldPlotter.reshapeField(field, x, y);
            surf(x, y, f, "EdgeColor", "none");
            view(35, 35)
            axis tight
            colorbar
        end

        function histogram(field, nBins)
            histogram(field(:), nBins, "Normalization", "pdf");
        end

        function variogramCloud(points, values, sampleCount)
            n = size(points, 1);
            sampleCount = min(sampleCount, n * (n - 1) / 2);
            idx1 = randi(n, sampleCount, 1);
            idx2 = randi(n, sampleCount, 1);
            d = sqrt(sum((points(idx1, :) - points(idx2, :)).^2, 2));
            g = 0.5 * (values(idx1) - values(idx2)).^2;
            scatter(d, g, 8, "filled", "MarkerFaceAlpha", 0.3);
            xlabel("Distance")
            ylabel("Semivariance")
        end
    end

    methods (Static, Access = private)
        function f = reshapeField(field, x, y)
            if isvector(field)
                f = reshape(field, numel(y), numel(x));
            else
                f = field;
            end
        end
    end
end
