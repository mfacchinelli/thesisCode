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

%% Load C++ Results For Propagation

%...Load translational motion
filename = '/Users/Michele/GitHub/tudat/tudatBundle/tudatApplications/Thesis/SimulationOutput/cartesianPropagated.dat';
fileID = fopen(filename,'r');
CartesianPropagatedResults = textscan(fileID,repmat('%f',[1,7]),'Delimiter',',','CollectOutput',true);
simulationTime = ( CartesianPropagatedResults{1}(:,1) - CartesianPropagatedResults{1}(1) ) / timeConversion; 
CartesianPropagatedResults = CartesianPropagatedResults{1}(:,2:end) / 1e3;
fclose(fileID);

%...Load Keplerian translational motion
filename = '/Users/Michele/GitHub/tudat/tudatBundle/tudatApplications/Thesis/SimulationOutput/keplerianPropagated.dat';
fileID = fopen(filename,'r');
KeplerianPropagatedResults = textscan(fileID,repmat('%f',[1,7]),'Delimiter',',','CollectOutput',true);
KeplerianPropagatedResults = KeplerianPropagatedResults{1}(:,2:end);
KeplerianPropagatedResults(:,1) = KeplerianPropagatedResults(:,1) / 1e3;
KeplerianPropagatedResults(:,3:end) = rad2deg(KeplerianPropagatedResults(:,3:end));
fclose(fileID);

%...Load rotational motion
filename = '/Users/Michele/GitHub/tudat/tudatBundle/tudatApplications/Thesis/SimulationOutput/rotationalPropagated.dat';
fileID = fopen(filename,'r');
rotationalPropagatedResults = textscan(fileID,repmat('%f',[1,8]),'Delimiter',',','CollectOutput',true);
rotationalPropagatedResults = rotationalPropagatedResults{1}(:,2:end);
rotationalPropagatedResults(:,6:8) = rad2deg(rotationalPropagatedResults(:,5:7));
rotationalPropagatedResults(:,5) = 1.0 - sqrt( sum( rotationalPropagatedResults(:,1:4).^2, 2 ) );
fclose(fileID);

%...Clean up
clear filename fileID

%% Load C++ Results For Estimation

%...Load translational motion
filename = '/Users/Michele/GitHub/tudat/tudatBundle/tudatApplications/Thesis/SimulationOutput/cartesianEstimated.dat';
fileID = fopen(filename,'r');
CartesianEstimatedResults = textscan(fileID,repmat('%f',[1,7]),'Delimiter',',','CollectOutput',true);
onboardTime = ( CartesianEstimatedResults{1}(:,1) - CartesianEstimatedResults{1}(1) ) / timeConversion; 
CartesianEstimatedResults = CartesianEstimatedResults{1}(:,2:end) / 1e3;
fclose(fileID);

%...Load Keplerian translational motion
filename = '/Users/Michele/GitHub/tudat/tudatBundle/tudatApplications/Thesis/SimulationOutput/keplerianEstimated.dat';
fileID = fopen(filename,'r');
KeplerianEstimatedResults = textscan(fileID,repmat('%f',[1,7]),'Delimiter',',','CollectOutput',true);
KeplerianEstimatedResults = KeplerianEstimatedResults{1}(:,2:end);
KeplerianEstimatedResults(:,1) = KeplerianEstimatedResults(:,1) / 1e3;
KeplerianEstimatedResults(:,3:end) = rad2deg(KeplerianEstimatedResults(:,3:end));
fclose(fileID);

%...Load rotational motion
filename = '/Users/Michele/GitHub/tudat/tudatBundle/tudatApplications/Thesis/SimulationOutput/rotationalEstimated.dat';
fileID = fopen(filename,'r');
rotationalEstimatedResults = textscan(fileID,repmat('%f',[1,8]),'Delimiter',',','CollectOutput',true);
rotationalEstimatedResults = rotationalEstimatedResults{1}(:,2:end);
rotationalEstimatedResults(:,6:8) = rad2deg(rotationalEstimatedResults(:,5:7));
rotationalEstimatedResults(:,5) = 1.0 - sqrt( sum( rotationalEstimatedResults(:,1:4).^2, 2 ) );
fclose(fileID);

%...Clean up
clear filename fileID

%% Plot States Over Time

%...Plot Cartesian translational motion
F = figure('rend','painters','pos',figSizeLarge);
for i = 1:size(CartesianPropagatedResults,2)
    subplot(2,3,i)
    hold on
    plot(simulationTime,CartesianPropagatedResults(:,i),'LineWidth',1.25)
    plot(onboardTime,CartesianEstimatedResults(:,i),'LineWidth',1.25)
    hold off
    xlabel(timeLabel)
    ylabel(CartesianLabels{i})
    set(gca,'FontSize',15)
    grid on
end
subplotLegend({'Actual','Estimated'})

%...Plot Keplerian translational motion
F = figure('rend','painters','pos',figSizeLarge);
for i = 1:size(KeplerianPropagatedResults,2)
    subplot(2,3,i)
    hold on
    plot(simulationTime,KeplerianPropagatedResults(:,i),'LineWidth',1.25)
    plot(onboardTime,KeplerianEstimatedResults(:,i),'LineWidth',1.25)
    hold off
    xlabel(timeLabel)
    ylabel(KeplerianLabels{i})
    set(gca,'FontSize',15)
    grid on
end
subplotLegend({'Actual','Estimated'})

%...Plot rotational motion
F = figure('rend','painters','pos',figSizeLarge);
for i = 1:size(rotationalPropagatedResults,2)
    subplot(2,4,i)
    hold on
    plot(simulationTime,rotationalPropagatedResults(:,i),'LineWidth',1.25)
    plot(onboardTime,rotationalEstimatedResults(:,i),'LineWidth',1.25)
    hold off
    xlabel(timeLabel)
    ylabel(rotationLabels{i})
    set(gca,'FontSize',15)
    grid on
end
subplotLegend({'Actual','Estimated'})