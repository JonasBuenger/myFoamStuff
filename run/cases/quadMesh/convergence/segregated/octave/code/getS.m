function [ s ] = getS( path )

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
s = zeros(numCells, 3);

% read in data
for i=1:numCells;
    [s(i,1), s(i,2), s(i,3)] = strread(fgetl(fh),'(%f %f %f)');
end

end

