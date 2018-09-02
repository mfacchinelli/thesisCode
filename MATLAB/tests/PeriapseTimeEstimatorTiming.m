fclose('all'); clear all; close all force; profile off; clc; format long g; rng default;
addpath functions tests

%% Settings

%...Figure settings
showFigure = true;
saveFigure = false;
[figSizeLarge,figSizeMedium,figSizeSmall] = saveFigureSettings(saveFigure);

%...Constants
marsRadius = 3389526.666666667;
marsGravitationalParameter = 4.282e13;
marsAtmosphericInterface = 5000;

%...Labels
timeConversion = 3600 * 24;
timeLabel = 'Time [d]';
CartesianLabels = {'x [km]','y [km]','z [km]','v_x [m s^{-1}]','v_y [m s^{-1}]','v_z [m s^{-1}]'};
CartesianLabelsDifference = {'\Delta x [km]','\Delta y [km]','\Delta z [km]',...
    '\Delta v_x [m s^{-1}]','\Delta v_y [m s^{-1}]','\Delta v_z [m s^{-1}]'};
KeplerianLabels = {'a [km]','e [-]','i [deg]','\Omega [deg]','\omega [deg]','\vartheta [deg]'};
rotationLabels = {'\eta [-]','\epsilon_1 [-]','\epsilon_2 [-]','\epsilon_3 [-]','Norm Offset [-]',...
    '\omega_1 [deg s^{-1}]','\omega_2 [deg s^{-1}]','\omega_3 [deg s^{-1}]'};

mainPath = '/Users/Michele/Library/Mobile Documents/com~apple~CloudDocs/University/Master Thesis/Code/MATLAB/data/';

%% Load C++ Results For Propagation

%...Load translational motion
filename = fullfile(mainPath,'trajectory.dat');
fileID = fopen(filename,'r');
CartesianResults = textscan(fileID,repmat('%f',[1,7]),'Delimiter',',','CollectOutput',true);
simulationTime = ( CartesianResults{1}(:,1) - CartesianResults{1}(1) ) / timeConversion; 
CartesianResults = CartesianResults{1}(:,2:end);
CartesianResults(:,1:3) = CartesianResults(:,1:3) / 1e3;
fclose(fileID);

%...Load Keplerian translational motion
filename = fullfile(mainPath,'orbit.dat');
fileID = fopen(filename,'r');
KeplerianResults = textscan(fileID,repmat('%f',[1,7]),'Delimiter',',','CollectOutput',true);
KeplerianResults = KeplerianResults{1}(:,2:end);
KeplerianResults(:,1) = KeplerianResults(:,1) / 1e3;
KeplerianResults(:,3:end) = rad2deg(KeplerianResults(:,3:end));
fclose(fileID);

%...Load accerations
filename = fullfile(mainPath,'dependent.dat');
fileID = fopen(filename,'r');
accelerationResults = textscan(fileID,repmat('%f',[1,46]),'Delimiter',',','CollectOutput',true);
accelerationResults = accelerationResults{1}(:,2:end);
fclose(fileID);

%...Compute supporting variables
altitude = sqrt(sum(CartesianResults(:,1:3).^2,2)) - marsRadius / 1e3;
centralGravity = accelerationResults(:,1:3);
sphericalHarmonics = accelerationResults(:,4:end);

%...Clean up
clear filename fileID

%% Plot States Over Time

%...Plot Cartesian translational motion
F = figure('rend','painters','pos',figSizeLarge);
for i = 1:size(CartesianResults,2)
    subplot(2,3,i)
    plot(simulationTime,CartesianResults(:,i),'LineWidth',1.25)
    xlabel(timeLabel)
    ylabel(CartesianLabels{i})
    set(gca,'FontSize',15)
    grid on
end

%...Plot Keplerian translational motion
F = figure('rend','painters','pos',figSizeLarge);
for i = 1:size(KeplerianResults,2)
    subplot(2,3,i)
    plot(simulationTime,KeplerianResults(:,i),'LineWidth',1.25)
    xlabel(timeLabel)
    ylabel(KeplerianLabels{i})
    set(gca,'FontSize',15)
    grid on
end

%% Plot PTE Timing Criteria

plotTrueAnomaly = false;

%...Get time and location of variables below atmospheric interface
reducedTime = simulationTime(1:end-1);
locAbove = altitude(1:end-1) >= 0.25 * ( KeplerianResults(1:end-1,1) .* ( 1 + KeplerianResults(1:end-1,2) ) );
locBelow = altitude(1:end-1) <= 0.25 * ( KeplerianResults(1:end-1,1) .* ( 1 + KeplerianResults(1:end-1,2) ) );

semiMajorAxisDerivative = diff(KeplerianResults(:,1))./diff(simulationTime);
semiMajorAxisDerivativeBelow = semiMajorAxisDerivative;
semiMajorAxisDerivativeBelow(locAbove) = NaN;
semiMajorAxisDerivativeAbove = semiMajorAxisDerivative;
semiMajorAxisDerivativeAbove(locBelow) = NaN;

trueAnomaly = KeplerianResults(1:end-1,6);
trueAnomalyBelow = trueAnomaly;
trueAnomalyBelow(locAbove) = NaN;
trueAnomalyAbove = trueAnomaly;
trueAnomalyAbove(locBelow) = NaN;

%...Plot Cartesian translational motion
F = figure('rend','painters','pos',figSizeSmall);

yyaxis left
hold on
plot(simulationTime,sqrt(sum(centralGravity.^2,2)),'LineWidth',1.25)
plot(simulationTime,sqrt(sum(sphericalHarmonics.^2,2)),'LineWidth',1.25)
hold off
ylabel('Acceleration [m s^{-2}]')
set(gca,'YScale','log')

if plotTrueAnomaly
    yyaxis right
    hold on
    plot(reducedTime,trueAnomalyAbove,'LineWidth',1.25)
    plot(reducedTime,trueAnomaly,'LineWidth',1.25)
    plot([simulationTime(1),simulationTime(end)],[180,180])
    hold off
    ylabel('True Anomaly [deg]')
else
    yyaxis right
    hold on
    plot(reducedTime,semiMajorAxisDerivativeAbove,'LineWidth',1.25)
    plot(reducedTime,semiMajorAxisDerivative,'LineWidth',1.25)
    hold off
    ylabel('Semi-major Axis Derivative [m s^{-1}]')
end

xlabel(timeLabel)
set(gca,'FontSize',15)
if plotTrueAnomaly
    legend('Central','Spherical Harmonics','Above','Below','Location','Best')
else
    L = legend('Central','Spherical Harmonics','Above','Below','Location','NE');
end
grid on

if ~plotTrueAnomaly && saveFigure
    xlim([13,13.8])
    saveas(F,'../../Report/figures/pte_timing_high_ecc','epsc')
    
    L.Position = [0.15 0.464285714285714 0.298214285714286 0.165476190476191];
    xlim([144.84,144.94])
    saveas(F,'../../Report/figures/pte_timing_low_ecc','epsc')
end