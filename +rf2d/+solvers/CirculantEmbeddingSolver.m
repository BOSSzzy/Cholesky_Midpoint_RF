classdef CirculantEmbeddingSolver < rf2d.solvers.BaseSolver
    properties
        X
        Y
        Model (1,1) string
        Parameters (1,1) struct
        Spectrum
        EmbedSize
    end

    methods
        function obj = CirculantEmbeddingSolver(x, y, model, parameters, meanValue)
            obj.X = x(:)';
            obj.Y = y(:)';
            obj.Model = lower(string(model));
            obj.Parameters = parameters;
            obj.NumPoints = numel(x) * numel(y);
            obj.Mean = meanValue;
        end

        function prepare(obj)
            rf2d.Validation.mustBeRegularGrid(obj.X, obj.Y);
            nx = numel(obj.X);
            ny = numel(obj.Y);
            hx = mean(diff(obj.X));
            hy = mean(diff(obj.Y));
            ix = [0:nx-1 nx:-1:1];
            iy = [0:ny-1 ny:-1:1];
            [DX, DY] = meshgrid(ix * hx, iy * hy);
            covLag = rf2d.CovarianceModel.lagged(obj.Model, obj.Parameters, DX', DY');
            S = real(fft2(covLag));
            if min(S(:)) < -1e-8
                error("rf2d:NumericalStability", "Circulant embedding is not positive semidefinite.");
            end
            obj.Spectrum = max(S, 0);
            obj.EmbedSize = size(obj.Spectrum);
        end

        function fields = sample(obj, nFields, useGPU, useParallel)
            nx = numel(obj.X);
            ny = numel(obj.Y);
            m = obj.EmbedSize(1);
            n = obj.EmbedSize(2);
            S = obj.Spectrum;
            if useGPU
                S = gpuArray(S);
            end
            if useParallel && ~useGPU && nFields > 1
                fields = zeros(nx * ny, nFields);
                parfor k = 1:nFields
                    w = randn(m, n) + 1i * randn(m, n);
                    if useGPU
                        w = gpuArray(w);
                    end
                    z = ifft2(sqrt(S / (m * n)) .* w);
                    f = real(z(1:nx, 1:ny));
                    fields(:, k) = obj.Mean + f(:);
                end
            else
                fields = zeros(nx * ny, nFields);
                for k = 1:nFields
                    w = randn(m, n) + 1i * randn(m, n);
                    if useGPU
                        w = gpuArray(w);
                    end
                    z = ifft2(sqrt(S / (m * n)) .* w);
                    f = real(z(1:nx, 1:ny));
                    if useGPU
                        f = gather(f);
                    end
                    fields(:, k) = obj.Mean + f(:);
                end
            end
        end
    end
end
