fclose('all'); clear all; close all force; profile off; clc; format long g; rng default;
addpath tests data functions

%% Load MCD

MCD = load('MCDMeanAtmosphere');
altitude = MCD.hs_interp;
longitude = MCD.lon_interp;
latitude = MCD.lat_interp;
tabularTimeAvgFull = MCD.tabular_interp;

%% Interpolate to Conditions

clc;

%...Get conditions
mode = 3;
switch mode
    case 0
        altitudeInput = altitude(1);
        longitudeInput = longitude(1);
        latitudeInput = latitude(1);
    case 1
        altitudeInput = 236.9862;
        longitudeInput = 72.98632;
        latitudeInput = -65.9762;
    case 2
        altitudeInput = 6.6202093e2;
        longitudeInput = -1.685714e+02;
        latitudeInput = -7.851064e+01;
    case 3
        altitudeInput = altitude(1);
        longitudeInput = 0;
        latitudeInput = latitude(end);
end

%...Interpolate at given conditions
for i = 1:5
    result( i ) = interp3(longitude,latitude,altitude,...
        tabularTimeAvgFull(:,:,:,i),longitudeInput,latitudeInput,altitudeInput,'linear');
end
result
sqrt( prod( result( 3:end ) ) )