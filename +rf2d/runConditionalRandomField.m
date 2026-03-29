function result = runConditionalRandomField(config)
if nargin < 1
    config = struct();
end
result = rf2d.ConditionalRandomField.run(config);
end
