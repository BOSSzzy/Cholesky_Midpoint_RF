function [RFC, RFPHI, c, phi] = RandomField2DCholMethod(coordFile, acfType, seed)
if nargin < 1
    coordFile = "Coord1.xlsx";
end
if nargin < 2
    acfType = 1;
end
if nargin < 3
    seed = 0;
end
[RFC, RFPHI, c, phi] = rf2d.LegacyAdapter.runLegacy(coordFile, acfType, seed);
if nargout == 0
    assignin("base", "RFC", RFC);
    assignin("base", "RFPHI", RFPHI);
    assignin("base", "c", c);
    assignin("base", "phi", phi);
end
end
