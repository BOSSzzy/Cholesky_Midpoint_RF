function result = randex(config)
if nargin < 1
    config = struct();
end
result = rf2d.runConditionalRandomField(config);
end
