function [ cellCentres ] = getCellCentres( path )

% files
fhx = fopen([path 'ccx'],'r');
fhy = fopen([path 'ccy'],'r');
fhz = fopen([path 'ccz'],'r');

% jump some lines
for i=1:20
    fgets(fhx);
    fgets(fhy);
    fgets(fhz);
end

% number of cells
numCells = str2double(fgets(fhx)); fgets(fhy); fgets(fhz);

% jump one line
fgets(fhx); fgets(fhy); fgets(fhz);

% initialize cellCentres
cellCentres = zeros(numCells, 3);

% read in data
for i=1:numCells;
    cellCentres(i,1) = str2double(fgets(fhx));
    cellCentres(i,2) = str2double(fgets(fhy));
    cellCentres(i,3) = str2double(fgets(fhz));
end

end