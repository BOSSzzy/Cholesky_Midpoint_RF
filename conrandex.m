function result = conrandex(config)
if nargin < 1
    config = struct();
end
result = rf2d.runConditionalRandomField(config);
result = result.conditional;
end
