clc; clear;

rootDir = fileparts(fileparts(mfilename("fullpath")));
demoDir = fullfile(rootDir, "demo");

config = struct();
config.grid.dx = 1.0;
config.grid.dy = 1.0;
config.covModel.type = "exponential";
config.covModel.nugget = 0.0;
config.covModel.rangeX = 20.0;
config.covModel.rangeY = 10.0;
config.prior.mean = 10.0;
config.prior.std = 2.0;
config.sim.nRealizations = 100;
config.sim.randomSeed = 20260329;
config.io.dataDir = demoDir;
config.io.gridTablePath = fullfile(demoDir, "crf_grid.csv");
config.io.observationPath = fullfile(demoDir, "crf_observation.csv");
config.io.unconditionalPath = fullfile(demoDir, "crf_unconditional_field.csv");
config.io.conditionalPath = fullfile(demoDir, "crf_conditional_field.csv");
config.logging.enabled = true;
config.logging.level = "info";

result = rf2d.runConditionalRandomField(config);
disp(result.conditional.conditionalField(1:5, 1:8));
