function [ Theta ] = getTheta( path )

% files
fh = fopen( path );

% jump some lines
for i=1:20
    fgetl(fh);
end

% number of cells
numCells = str2double(fgetl(fh));

% jump one line
fgetl(fh);

% initialize cellCentres
Theta = zeros(numCells, 1);

% read in data
for i=1:numCells;
    Theta(i,1) = str2double(fgetl(fh));
end

end

