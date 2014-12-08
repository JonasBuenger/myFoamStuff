function [ cellVolumes ] = getCellVolumes( path )

% files
fh = fopen(path,'r');

% jump some lines
for i=1:20
    fgets(fh);
end

% number of cells
numCells = str2double(fgets(fh));

% jump one line
fgets(fh);

% initialize cellCentres
cellVolumes = zeros(numCells, 1);

% read in data
for i=1:numCells;
    cellVolumes(i,1) = str2double(fgets(fh));
end

end

