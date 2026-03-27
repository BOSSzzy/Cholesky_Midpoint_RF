function outputFile = buildToolbox(projectFile, outputFolder)
if nargin < 1
    projectFile = fullfile(pwd, "toolbox", "rf2d_toolbox.prj");
end
if nargin < 2
    outputFolder = fullfile(pwd, "dist");
end
if ~exist(outputFolder, "dir")
    mkdir(outputFolder);
end
if exist(projectFile, "file")
    outputFile = matlab.addons.toolbox.packageToolbox(projectFile, outputFolder);
else
    error("rf2d:Toolbox", "Toolbox project file not found at %s", projectFile);
end
end
