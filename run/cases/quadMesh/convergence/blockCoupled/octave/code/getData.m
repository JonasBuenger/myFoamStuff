function [ cc cv s T ] = getData( path2data )

cc = getCellCentres([path2data '0/']);
cv = getCellVolumes([path2data 'constant/polyMesh/Volumes']);

s  = getS([path2data '120/s']);
T  = getTheta([path2data '120/Theta']);

end

