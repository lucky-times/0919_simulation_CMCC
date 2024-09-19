% Define pattern parameters
clear all;
fq = 4.9e9;
azvec = -180:180;
elvec = -90:90;
Am = 30; % Maximum attenuation (dB)

tilt = 0; % gamma角度

phi_3dB = 120; % 3 dB bandwidth in azimuth
theta_3dB = 6; % 3 dB bandwidth in elevation

[Azvec, Elvec]= meshgrid(azvec, elvec);


BS_anten = BSAntenna(16, 8, 1, phi_3dB, theta_3dB, tilt);
combinedMagPattern = BS_anten.elementPattern(Elvec, Azvec);

phasepattern = zeros(size(combinedMagPattern));
antennaElement = phased.CustomAntennaElement(...
    'AzimuthAngles',azvec, ...
    'ElevationAngles',elvec, ...
    'MagnitudePattern',combinedMagPattern, ...
    'PhasePattern',phasepattern);

f = figure;
pattern(antennaElement,fq);
title("天线类型II，无上倾");
