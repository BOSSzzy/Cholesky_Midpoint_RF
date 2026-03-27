classdef (Abstract) BaseSolver < handle
    properties
        NumPoints (1,1) double
        Mean (1,1) double = 0
    end

    methods (Abstract)
        prepare(obj)
        fields = sample(obj, nFields, useGPU, useParallel)
    end
end
