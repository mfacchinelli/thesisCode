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
mode = 0;
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

%% Test Mean Mean Mean Atmosphere

%...Open actual atmosphere
filename = '/Users/Michele/GitHub/tudat/tudatBundle/tudat/Tudat/External/AtmosphereTables/MCDMeanAtmosphere.dat';
fileID = fopen(filename,'r');
atmosphere = textscan(fileID,'%f %f %f %f %f %f','CollectOutput',true); atmosphere = atmosphere{1};
fclose(fileID);

%...Load MCD mat-file
MCD = load('/Users/Michele/Library/Mobile Documents/com~apple~CloudDocs/University/Master Thesis/Code/SPARTA/MCDEnvir.mat');

%...Open Tudat interpolated atmosphere
filename = '/Users/Michele/GitHub/tudat/tudatBundle/tudatApplications/Test/SimulationOutput/atmosphere.dat';
fileID = fopen(filename,'r');
atmosphereTudat = textscan(fileID,'%f,%f,%f,%f','CollectOutput',true); atmosphereTudat = atmosphereTudat{1};
fclose(fileID);

%...Interpolate to match Tudat conditions
altitude = atmosphereTudat(:,1);
atmosphereInterpolated = interp1(atmosphere(:,1),atmosphere,altitude,'spline');
MCDInterpolated(:,1) = interp1(MCD.altitude*1e3,MCD.density,altitude,'spline');
MCDInterpolated(:,2) = interp1(MCD.altitude*1e3,MCD.pressure,altitude,'spline');
MCDInterpolated(:,3) = interp1(MCD.altitude*1e3,MCD.temperature,altitude,'spline');
MCDInterpolated(:,4) = interp1(MCD.altitude*1e3,MCD.numberDensity,altitude,'spline');
MCDInterpolated(:,5) = interp1(MCD.altitude*1e3,MCD.speedOfSound,altitude,'spline');

%...Compare plots
figure;
for i = 1:3
    subplot(1,3,i)
    hold on
    plot(altitude,abs(atmosphereInterpolated(:,i)-atmosphereTudat(:,i+1)),'LineWidth',1.25)
    plot(altitude,abs(MCDInterpolated(:,i)-atmosphereTudat(:,i+1)),'LineWidth',1.25)
    plot(altitude,abs(atmosphereInterpolated(:,i)-MCDInterpolated(:,i)),'LineWidth',1.25)
    hold off
    grid on
    set(gca,'FontSize',15,'YScale','log')
    legend('File-Tudat','Mat-Tudat','File-Mat')
end

%...Tables
array2table(vertcat(horzcat(atmosphereInterpolated(find(altitude == 100e3,1),2:4),...
    6.022140857e23/8.3144598 * atmosphereInterpolated(find(altitude == 100e3,1),2) * atmosphereInterpolated(find(altitude == 100e3,1),5),...
    sqrt(atmosphereInterpolated(find(altitude == 100e3,1),6)*atmosphereInterpolated(find(altitude == 100e3,1),5)*atmosphereInterpolated(find(altitude == 100e3,1),4))),...
    MCDInterpolated(find(altitude == 100e3,1),:)),'VariableNames',...
    {'Density','Pressure','Temperature','NumberDensity','SpeedOfSound'},'RowNames',{'File','Mat'})

%% Number Density

%...Find number density with equation
numberDensity = 6.022140857e23/8.3144598 * MCD.density .* MCD.gasConstant;

%...Plot comparison
figure;
subplot(1,2,1)
hold on
plot(MCD.altitude,MCD.numberDensity,'LineWidth',1.25)
plot(MCD.altitude,numberDensity,'LineWidth',1.25)
hold off
grid on
legend('Actual','Computed')
set(gca,'FontSize',15,'XScale','log','YScale','log')

subplot(1,2,2)
plot(MCD.altitude,abs(MCD.numberDensity-numberDensity),'LineWidth',1.25)
grid on
set(gca,'FontSize',15,'XScale','log','YScale','log')
