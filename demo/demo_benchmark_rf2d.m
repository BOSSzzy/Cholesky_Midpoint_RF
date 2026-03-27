x = linspace(0, 1000, 256);
y = linspace(0, 1000, 256);
params = struct("variance", 1.0, "corrLength", [80 80], "smoothness", 1.5, "nugget", 1e-8);
solvers = ["cholesky", "kl", "circulant", "approximate"];
sizes = [256 512 1024 2048];
results = struct("solver", [], "n", [], "time", [], "memoryGB", []);
row = 1;
for n = sizes
    xi = linspace(0, 1000, n);
    yi = linspace(0, 1000, n);
    for s = solvers
        tStart = tic;
        g = rf2d.createGenerator(xi, yi, "matern", params, "Solver", s, "Seed", 1, "KLModes", 256, "MaxRank", 1024, "MaxDirectPoints", 2e5);
        f = g.realize();
        elapsed = toc(tStart);
        mem = whos("f");
        results(row).solver = char(s);
        results(row).n = n;
        results(row).time = elapsed;
        results(row).memoryGB = mem.bytes / 1024^3;
        row = row + 1;
    end
end
T = struct2table(results);
disp(T);
g = rf2d.createGenerator(x, y, "matern", params, "Solver", "circulant", "Seed", 7);
f = g.realize();
figure; g.contourPlot(f, 30);
figure; g.surfacePlot(f);
figure; g.histogramPlot(f, 50);
figure; g.variogramCloud(f, 5000);
