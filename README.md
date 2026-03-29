# Cholesky_Midpoint_RF

Cholesky_Midpoint_RF is a MATLAB project for generating 2D Gaussian random fields on structured grids.  
It provides multiple covariance kernels, solver backends, visualization/export utilities, benchmarking scripts, and a backward-compatible legacy entry point.

## Core capabilities

- Covariance models: `gaussian`, `exponential`, `spherical`, `matern`
- Solver backends: `cholesky`, `kl`, `circulant`, `approximate`, and `auto` selection
- Generation modes: single realization, batch generation, optional parallel and GPU workflows
- Conditional random fields (CRF): unconditional + conditioned realizations from grid and observation tables
- Visualization: contour, surface, histogram, variogram cloud
- Export: CSV, VTK, MAT
- Compatibility: legacy wrapper `RandomField2DCholMethod.m`

## Repository structure

- `+rf2d/`: main package
  - `createGenerator.m`: public factory API
  - `RandomFieldGenerator.m`: high-level generator object
  - `CovarianceModel.m`: covariance kernels and pairwise/lagged covariance
  - `ConditionalRandomField.m`: CRF workflow for unconditional + conditional fields
  - `runConditionalRandomField.m`: public CRF entry point
  - `Validation.m`: model/solver/grid/parameter validation
  - `LegacyAdapter.m`: bridge for the legacy function interface
  - `+solvers/`: solver implementations and auto-selection factory
  - `+viz/`: plotting helpers
  - `+io/`: export helpers
- `RandomField2DCholMethod.m`: historical API wrapper
- `randex.m` / `conrandex.m`: CRF compatibility wrappers
- `demo/`: benchmark and Live Script conversion utility
- `tests/`: MATLAB unit tests
- `+docs/generateDocs.m`: HTML documentation publishing
- `+toolbox/buildToolbox.m`: `.mltbx` packaging helper
- `toolbox/rf2d_toolbox.prj`: toolbox project file

## Requirements

- MATLAB R2020b or newer
- Statistics and Machine Learning Toolbox
- Parallel Computing Toolbox (optional, for `parfor`, `spmd`, and GPU workflows)

## Setup

Clone the repository and add it to your MATLAB path:

```matlab
addpath(genpath(pwd))
```

## Quick start

```matlab
x = linspace(0, 100, 128);
y = linspace(0, 40, 64);
params = struct("variance", 2.0, "corrLength", [20 8], "smoothness", 1.2, "nugget", 1e-8);
g = rf2d.createGenerator(x, y, "matern", params, "Solver", "auto", "Seed", 1234);
f = g.realize();
figure; g.contourPlot(f, 20);
```

## API overview

`rf2d.createGenerator(x, y, covarianceModel, covarianceParameters, Name=Value)` supports:

- `Solver`: `"auto"` | `"cholesky"` | `"kl"` | `"circulant"` | `"approximate"`
- `Seed`: random seed
- `UseGPU`: logical flag for solver sampling
- `Mean`: scalar mean value of the field
- `MaxDirectPoints`: threshold used by automatic solver routing
- `MaxRank`: rank for Nyström approximation
- `KLModes`: number of KL modes
- `Jitter`: numerical stabilization term

Common generator methods:

- `realize()` or `realize("UseGPU", ..., "UseParallel", ...)`
- `generateBatch(nFields, useGPU, useParallel)`
- `reseed(seed)`
- `contourPlot`, `surfacePlot`, `histogramPlot`, `variogramCloud`
- `exportCSV`, `exportVTK`, `exportMAT`

CRF API:

- `rf2d.runConditionalRandomField(config)`: file-based CRF pipeline
- `rf2d.ConditionalRandomField.fromTables(gridTable, observationTable, config)`: in-memory CRF pipeline
- `randex(config)` and `conrandex(config)`: compatibility entry points

Minimal CRF config fields:

- `grid.dx`, `grid.dy`
- `covModel.type`, `covModel.rangeX`, `covModel.rangeY`, `covModel.nugget`
- `prior.mean`, `prior.std`
- `sim.nRealizations`, `sim.randomSeed`, `sim.jitter`
- `io.gridTablePath`, `io.observationPath`, `io.unconditionalPath`, `io.conditionalPath`
- `logging.enabled`, `logging.level`

CRF input tables:

- Grid table columns: `id`, `j`, `k`
- Observation table columns: `j`, `k`, `obsValue`, `obsVar`

## Solver selection behavior

When `Solver="auto"`:

- uses Cholesky for smaller problems (`n <= MaxDirectPoints`)
- uses circulant embedding on regular grids for larger problems
- falls back to approximate Nyström for large non-regular cases

If solver preparation fails for numerical reasons, the generator falls back to the approximate solver.

## Legacy entry point

The historical function is preserved:

```matlab
[RFC, RFPHI, c, phi] = RandomField2DCholMethod("Coord1.xlsx", 1, 0);
```

It delegates to `rf2d.LegacyAdapter` and keeps output compatibility.

## Demo and benchmark

Run:

```matlab
demo.demo_benchmark_rf2d
```

Run the CRF demo with project-level sample tables (outside any `data` subfolder):

```matlab
demo.demo_crf_rf2d
```

To convert the demo script into a Live Script:

```matlab
demo.createLiveScript
```

## Documentation generation

Run:

```matlab
docs.generateDocs
```

Generated HTML is written to `docs/html/`.

## Toolbox packaging

Run:

```matlab
toolbox.buildToolbox
```

This creates an installable `.mltbx` under `dist/`.

## Testing

Run all tests:

```matlab
results = runtests("tests", "IncludeSubfolders", true);
table(results)
```

Run tests with coverage output:

```matlab
import matlab.unittest.TestRunner
import matlab.unittest.plugins.CodeCoveragePlugin
import matlab.unittest.plugins.codecoverage.CoverageReport

suite = testsuite("tests", "IncludeSubfolders", true);
runner = TestRunner.withTextOutput;
runner.addPlugin(CodeCoveragePlugin.forFolder(pwd, "Producing", CoverageReport("tests/coverage")));
results = runner.run(suite);
```
