classdef FieldExporter
    methods (Static)
        function toCSV(X, Y, field, filePath)
            filePath = char(filePath);
            if isvector(field)
                F = reshape(field, size(X));
            else
                F = field;
            end
            T = table(X(:), Y(:), F(:), 'VariableNames', {'x', 'y', 'value'});
            writetable(T, filePath);
        end

        function toMAT(X, Y, field, filePath, variableName)
            filePath = char(filePath);
            if isvector(field)
                F = reshape(field, size(X));
            else
                F = field;
            end
            S.X = X;
            S.Y = Y;
            S.(char(variableName)) = F;
            save(filePath, "-struct", "S");
        end

        function toVTK(X, Y, field, filePath)
            filePath = char(filePath);
            if isvector(field)
                F = reshape(field, size(X));
            else
                F = field;
            end
            nx = size(X, 2);
            ny = size(X, 1);
            dx = mean(diff(X(1, :)));
            dy = mean(diff(Y(:, 1)));
            fid = fopen(filePath, "w");
            if fid < 0
                error("rf2d:IO", "Cannot open file.");
            end
            c = onCleanup(@() fclose(fid));
            fprintf(fid, "# vtk DataFile Version 3.0\n");
            fprintf(fid, "rf2d\n");
            fprintf(fid, "ASCII\n");
            fprintf(fid, "DATASET STRUCTURED_POINTS\n");
            fprintf(fid, "DIMENSIONS %d %d 1\n", nx, ny);
            fprintf(fid, "ORIGIN %g %g 0\n", X(1, 1), Y(1, 1));
            fprintf(fid, "SPACING %g %g 1\n", dx, dy);
            fprintf(fid, "POINT_DATA %d\n", nx * ny);
            fprintf(fid, "SCALARS field double 1\n");
            fprintf(fid, "LOOKUP_TABLE default\n");
            fprintf(fid, "%g\n", F');
            clear c
        end
    end
end
