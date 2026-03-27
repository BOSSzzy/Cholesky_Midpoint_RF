function createLiveScript()
sourceFile = fullfile(pwd, "demo", "demo_benchmark_rf2d.m");
targetFile = fullfile(pwd, "demo", "demo_benchmark_rf2d.mlx");
if ~isfile(sourceFile)
    error("rf2d:Demo", "Source script not found.");
end
if exist("matlab.internal.liveeditor.openAndConvert", "file")
    matlab.internal.liveeditor.openAndConvert(sourceFile, targetFile);
else
    error("rf2d:Demo", "Live Editor conversion API is unavailable in this MATLAB version.");
end
end
