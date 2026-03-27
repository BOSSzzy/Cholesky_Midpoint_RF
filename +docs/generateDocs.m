function generateDocs(outputDir)
if nargin < 1
    outputDir = fullfile(pwd, "docs", "html");
end
if ~exist(outputDir, "dir")
    mkdir(outputDir);
end
targets = {
    fullfile(pwd, "+rf2d", "createGenerator.m")
    fullfile(pwd, "+rf2d", "RandomFieldGenerator.m")
    fullfile(pwd, "+rf2d", "CovarianceModel.m")
    fullfile(pwd, "RandomField2DCholMethod.m")
    };
for i = 1:numel(targets)
    publish(targets{i}, struct("format", "html", "outputDir", outputDir, "showCode", true));
end
end
