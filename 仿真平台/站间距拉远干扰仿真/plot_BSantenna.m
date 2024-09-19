% Define pattern parameters
clear all;
fq = 4.9e9;
azvec = -180:180;
elvec = 0:90;
antenna = antennas.BSAntenna(6, 8, 1, 0);

% Define antenna pattern
[az,el] = meshgrid(azvec,elvec);
combinedMagPattern = zeros(size(az));

for i = 1:size(az, 1)
    for j = 1:size(az, 2)
      combinedMagPattern(i,j) = antenna.elementPattern(el(i,j), az(i,j), 0);
    end
end

phasepattern = zeros(size(combinedMagPattern));

% Create antenna element
antennaElement = phased.CustomAntennaElement(...
    'AzimuthAngles',azvec, ...
    'ElevationAngles',elvec, ...
    'MagnitudePattern',combinedMagPattern, ...
    'PhasePattern',phasepattern);
% Display radiation pattern
f = figure;
pattern(antennaElement,fq);