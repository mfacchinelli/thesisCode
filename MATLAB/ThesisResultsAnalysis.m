fclose('all'); clear all; close all force; profile off; clc; format long g; rng default;
addpath functions

%% Settings

%...Figure settings
showFigure = true;
saveFigure = false;
[figSizeLarge,figSizeMedium,figSizeSmall] = saveFigureSettings(saveFigure);

%...Labels
timeConversion = 3600 * 24;
timeLabel = 'Time [d]';
CartesianLabels = {'x [km]','y [km]','z [km]','v_x [km s^{-1}]','v_y [km s^{-1}]','v_z [km s^{-1}]'};
KeplerianLabels = {'a [km]','e [-]','i [deg]','\Omega [deg]','\omega [deg]','\vartheta [deg]'};
rotationLabels = {'\eta [-]','\epsilon_1 [-]','\epsilon_2 [-]','\epsilon_3 [-]','Norm Offset [-]',...
    '\omega_1 [deg s^{-1}]','\omega_2 [deg s^{-1}]','\omega_3 [deg s^{-1}]'};

%% Load C++ Results

%...Load translational motion
filename = '/Users/Michele/GitHub/tudat/tudatBundle/tudatApplications/Thesis/SimulationOutput/CartesianTranslational.dat';
fileID = fopen(filename,'r');
CartesianTranslationalMotion = textscan(fileID,repmat('%f',[1,7]),'Delimiter',',','CollectOutput',true);
simulationTime = ( CartesianTranslationalMotion{1}(:,1) - CartesianTranslationalMotion{1}(1) ) / timeConversion; 
CartesianTranslationalMotion = CartesianTranslationalMotion{1}(:,2:end) / 1e3;
fclose(fileID);

%...Load Keplerian translational motion
filename = '/Users/Michele/GitHub/tudat/tudatBundle/tudatApplications/Thesis/SimulationOutput/KeplerianTranslational.dat';
fileID = fopen(filename,'r');
KeplerianTranslationalMotion = textscan(fileID,repmat('%f',[1,7]),'Delimiter',',','CollectOutput',true);
KeplerianTranslationalMotion = KeplerianTranslationalMotion{1}(:,2:end);
KeplerianTranslationalMotion(:,1) = KeplerianTranslationalMotion(:,1) / 1e3;
KeplerianTranslationalMotion(:,3:end) = rad2deg(KeplerianTranslationalMotion(:,3:end));
fclose(fileID);

%...Load rotational motion
filename = '/Users/Michele/GitHub/tudat/tudatBundle/tudatApplications/Thesis/SimulationOutput/rotational.dat';
fileID = fopen(filename,'r');
rotationalMotion = textscan(fileID,repmat('%f',[1,8]),'Delimiter',',','CollectOutput',true);
rotationalMotion = rotationalMotion{1}(:,2:end);
rotationalMotion(:,6:8) = rad2deg(rotationalMotion(:,5:7));
rotationalMotion(:,5) = 1.0 - sqrt( sum( rotationalMotion(:,1:4).^2, 2 ) );
fclose(fileID);

%...Clean up
clear filename fileID

%% Plot States Over Time

%...Plot Cartesian translational motion
F = figure('rend','painters','pos',figSizeLarge);
for i = 1:size(CartesianTranslationalMotion,2)
    subplot(2,3,i)
    plot(simulationTime,CartesianTranslationalMotion(:,i),'LineWidth',1.25)
    xlabel(timeLabel)
    ylabel(CartesianLabels{i})
    set(gca,'FontSize',15)
    grid on
end

%...Plot Keplerian translational motion
F = figure('rend','painters','pos',figSizeLarge);
for i = 1:size(KeplerianTranslationalMotion,2)
    subplot(2,3,i)
    plot(simulationTime,KeplerianTranslationalMotion(:,i),'LineWidth',1.25)
    xlabel(timeLabel)
    ylabel(KeplerianLabels{i})
    set(gca,'FontSize',15)
    grid on
end

%...Plot rotational motion
F = figure('rend','painters','pos',figSizeLarge);
for i = 1:size(rotationalMotion,2)
    subplot(2,4,i)
    plot(simulationTime,rotationalMotion(:,i),'LineWidth',1.25)
    xlabel(timeLabel)
    ylabel(rotationLabels{i})
    set(gca,'FontSize',15)
    grid on
end