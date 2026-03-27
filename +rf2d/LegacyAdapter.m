classdef LegacyAdapter
    methods (Static)
        function [RFC, RFPHI, c, phi] = runLegacy(coordFile, acfType, seed)
            if nargin < 1 || strlength(string(coordFile)) == 0
                coordFile = "Coord1.xlsx";
            end
            if nargin < 2
                acfType = 1;
            end
            if nargin < 3
                seed = 0;
            end
            coord = readmatrix(coordFile);
            xUnique = unique(coord(:, 1))';
            yUnique = unique(coord(:, 2))';
            modelMap = {'exponential', 'gaussian', 'exponential', 'matern', 'spherical'};
            model = modelMap{max(1, min(5, acfType))};
            mu = [10; 30];
            covv = [0.3; 0.2];
            sigma = mu .* covv;
            corrLength = [40 4];
            params = struct("variance", 1, "corrLength", corrLength, "smoothness", 1.5, "nugget", 1e-10);
            g = rf2d.createGenerator(xUnique, yUnique, model, params, "Solver", "cholesky", "Seed", seed, "Mean", 0);
            f1 = g.realize();
            f2 = g.realize();
            sLn = sqrt(log(1 + (sigma ./ mu).^2));
            mLn = log(mu) - sLn.^2 / 2;
            c = exp(mLn(1) + sLn(1) * f1(:));
            phi = exp(mLn(2) + sLn(2) * f2(:));
            RFC = reshape(c, numel(yUnique), numel(xUnique));
            RFPHI = reshape(phi, numel(yUnique), numel(xUnique));
        end
    end
end
