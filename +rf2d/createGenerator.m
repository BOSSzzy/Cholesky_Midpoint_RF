function generator = createGenerator(x, y, covarianceModel, covarianceParameters, options)
arguments
    x {mustBeNumeric, mustBeVector}
    y {mustBeNumeric, mustBeVector}
    covarianceModel (1,1) string
    covarianceParameters (1,1) struct
    options.Solver (1,1) string = "auto"
    options.Seed (1,1) double = 1
    options.UseGPU (1,1) logical = false
    options.Mean (1,1) double = 0
    options.MaxDirectPoints (1,1) double = 2.5e4
    options.MaxRank (1,1) double = 512
    options.KLModes (1,1) double = 256
    options.Jitter (1,1) double = 1e-10
end
generator = rf2d.RandomFieldGenerator(x, y, covarianceModel, covarianceParameters, options);
end
