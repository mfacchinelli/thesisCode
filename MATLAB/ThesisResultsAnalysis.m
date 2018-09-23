fclose('all'); clear all; close all force; profile off; clc; format long g; rng default;
addpath functions

%% Settings

%...Figure settings
showFigure = true;
saveFigure = false;
[figSizeLarge,figSizeMedium,figSizeSmall] = saveFigureSettings(saveFigure);

%...Constants
marsRadius = 3389526.666666667;
marsGravitationalParameter = 42828375815756.1;
marsAtmosphericInterface = 175;

%...Plot settings
loadMeasurements = false;
loadFilter = false;
applyInterpolation = false;

%...Labels
timeConversion = 3600 * 24;
timeLabel = 'Time [d]';
CartesianLabels = {'x [km]','y [km]','z [km]','v_x [m s^{-1}]','v_y [m s^{-1}]','v_z [m s^{-1}]'};
CartesianLabelsDifference = {'\Delta x [km]','\Delta y [km]','\Delta z [km]',...
    '\Delta v_x [m s^{-1}]','\Delta v_y [m s^{-1}]','\Delta v_z [m s^{-1}]'};
KeplerianLabels = {'a [km]','e [-]','i [deg]','\Omega [deg]','\omega [deg]','\vartheta [deg]'};
rotationLabels = {'\eta [-]','\epsilon_1 [-]','\epsilon_2 [-]','\epsilon_3 [-]','Norm Offset [-]',...
    '\omega_1 [deg s^{-1}]','\omega_2 [deg s^{-1}]','\omega_3 [deg s^{-1}]'};

%% Load C++ Results For Propagation

outputFolder = 'SimulationOutputTransOnlyIMAN/';
% outputFolder = 'SimulationOutputTransOnlyAltimeter/';

%...Load translational motion
filename = ['/Users/Michele/GitHub/tudat/tudatBundle/tudatApplications/Thesis/',outputFolder,'/cartesianPropagated.dat'];
fileID = fopen(filename,'r');
CartesianPropagatedResults = textscan(fileID,repmat('%f',[1,7]),'Delimiter',',','CollectOutput',true);
initialTime = CartesianPropagatedResults{1}(1);
simulationTime = ( CartesianPropagatedResults{1}(:,1) - initialTime ) / timeConversion;
CartesianPropagatedResults = CartesianPropagatedResults{1}(:,2:end);
CartesianPropagatedResults(:,1:3) = CartesianPropagatedResults(:,1:3) / 1e3;
fclose(fileID);

%...Load Keplerian translational motion
filename = ['/Users/Michele/GitHub/tudat/tudatBundle/tudatApplications/Thesis/',outputFolder,'/keplerianPropagated.dat'];
fileID = fopen(filename,'r');
KeplerianPropagatedResults = textscan(fileID,repmat('%f',[1,7]),'Delimiter',',','CollectOutput',true);
KeplerianPropagatedResults = KeplerianPropagatedResults{1}(:,2:end);
KeplerianPropagatedResults(:,1) = KeplerianPropagatedResults(:,1) / 1e3;
KeplerianPropagatedResults(:,3:end) = rad2deg(KeplerianPropagatedResults(:,3:end));
fclose(fileID);

% %...Load rotational motion
% filename = ['/Users/Michele/GitHub/tudat/tudatBundle/tudatApplications/Thesis/',outputFolder,'/rotationalPropagated.dat'];
% fileID = fopen(filename,'r');
% rotationalPropagatedResults = textscan(fileID,repmat('%f',[1,8]),'Delimiter',',','CollectOutput',true);
% rotationalPropagatedResults = rotationalPropagatedResults{1}(:,2:end);
% rotationalPropagatedResults(:,6:8) = rad2deg(rotationalPropagatedResults(:,5:7));
% rotationalPropagatedResults(:,5) = 1.0 - quatnorm(rotationalPropagatedResults(:,1:4));
% fclose(fileID);

%...Clean up
clear filename fileID

%% Load C++ Results For Estimation

%...Load translational motion
filename = ['/Users/Michele/GitHub/tudat/tudatBundle/tudatApplications/Thesis/',outputFolder,'/cartesianEstimated.dat'];
fileID = fopen(filename,'r');
CartesianEstimatedResults = textscan(fileID,repmat('%f',[1,7]),'Delimiter',',','CollectOutput',true);
onboardTime = ( CartesianEstimatedResults{1}(:,1) - CartesianEstimatedResults{1}(1) ) / timeConversion; 
CartesianEstimatedResults = CartesianEstimatedResults{1}(:,2:end);
CartesianEstimatedResults(:,1:3) = CartesianEstimatedResults(:,1:3) / 1e3;
fclose(fileID);

%...Load Keplerian translational motion
filename = ['/Users/Michele/GitHub/tudat/tudatBundle/tudatApplications/Thesis/',outputFolder,'/keplerianEstimated.dat'];
fileID = fopen(filename,'r');
KeplerianEstimatedResults = textscan(fileID,repmat('%f',[1,7]),'Delimiter',',','CollectOutput',true);
KeplerianEstimatedResults = KeplerianEstimatedResults{1}(:,2:end);
KeplerianEstimatedResults(:,1) = KeplerianEstimatedResults(:,1) / 1e3;
KeplerianEstimatedResults(:,3:end) = rad2deg(KeplerianEstimatedResults(:,3:end));
fclose(fileID);

% %...Load rotational motion
% filename = ['/Users/Michele/GitHub/tudat/tudatBundle/tudatApplications/Thesis/',outputFolder,'/rotationalEstimated.dat'];
% fileID = fopen(filename,'r');
% rotationalEstimatedResults = textscan(fileID,repmat('%f',[1,8]),'Delimiter',',','CollectOutput',true);
% rotationalEstimatedResults = rotationalEstimatedResults{1}(:,2:end);
% rotationalEstimatedResults(:,6:8) = rad2deg(rotationalEstimatedResults(:,5:7));
% rotationalEstimatedResults(:,5) = 1.0 - quatnorm(rotationalEstimatedResults(:,1:4));
% fclose(fileID);

%...Clean up
clear filename fileID

%% Load C++ Results For Navigation Filter

%...Only if filtering is toggled
if loadFilter
    %...Load translational motion
    filename = ['/Users/Michele/GitHub/tudat/tudatBundle/tudatApplications/Thesis/',outputFolder,'/filterStateEstimates.dat'];
    fileID = fopen(filename,'r');
    filterStateEstimatedResults = textscan(fileID,repmat('%f',[1,13]),'Delimiter',',','CollectOutput',true);
    filterTime = ( filterStateEstimatedResults{1}(:,1) - initialTime ) / timeConversion;
    filterStateEstimatedResults = filterStateEstimatedResults{1}(:,2:end);
    filterStateEstimatedResults(:,1:3) = filterStateEstimatedResults(:,1:3)/1e3;
    fclose(fileID);
    
    %...Load Keplerian translational motion
    filename = ['/Users/Michele/GitHub/tudat/tudatBundle/tudatApplications/Thesis/',outputFolder,...
        '/filterCovarianceEstimates.dat'];
    fileID = fopen(filename,'r');
    filterCovarianceEstimatedResults = textscan(fileID,repmat('%f',[1,13]),'Delimiter',',','CollectOutput',true);
    filterCovarianceEstimatedResults = filterCovarianceEstimatedResults{1}(:,2:end);
    filterCovarianceEstimatedResults(:,1:3) = filterCovarianceEstimatedResults(:,1:3)/1e3;
    fclose(fileID);
       
    %...Clean up
    clear filename fileID
end

%% Load C++ Results For Dependent Variables

%...Load dependent variables
filename = ['/Users/Michele/GitHub/tudat/tudatBundle/tudatApplications/Thesis/',outputFolder,'/dependentVariables.dat'];
fileID = fopen(filename,'r');
dependentVariables = textscan(fileID,repmat('%f',[1,13]),'Delimiter',',','CollectOutput',true);
dependentVariables = dependentVariables{1}(:,2:end); dependentVariables(:,7:9) = rad2deg(dependentVariables(:,7:9));
fclose(fileID);

% %...Load control torques
% filename = ['/Users/Michele/GitHub/tudat/tudatBundle/tudatApplications/Thesis/',outputFolder,'/controlTorques.dat'];
% fileID = fopen(filename,'r');
% controlTorques = textscan(fileID,repmat('%f',[1,4]),'Delimiter',',','CollectOutput',true);
% controlTorques = controlTorques{1}(:,2:end);
% fclose(fileID);

%...Clean up
clear filename fileID

%% Load C++ Results For Measurements

%...Only if measurements are toggled
if loadMeasurements
    %...Load IMU measurements
    filename = ['/Users/Michele/GitHub/tudat/tudatBundle/tudatApplications/Thesis/',outputFolder,...
        '/accelerometerMeasurements.dat'];
    fileID = fopen(filename,'r');
    accelerometerMeasurements = textscan(fileID,repmat('%f',[1,7]),'Delimiter',',','CollectOutput',true);
    measurementTime = ( accelerometerMeasurements{1}(:,1) - initialTime ) / timeConversion;
    accelerometerMeasurements = accelerometerMeasurements{1}(:,2:end);
    fclose(fileID);
    
    %...Load expected measurements
    filename = ['/Users/Michele/GitHub/tudat/tudatBundle/tudatApplications/Thesis/',outputFolder,'/expectedlMeasurements.dat'];
    fileID = fopen(filename,'r');
    expectedMeasurements = textscan(fileID,repmat('%f',[1,4]),'Delimiter',',','CollectOutput',true);
    expectedMeasurements = expectedMeasurements{1}(:,2:end);
    fclose(fileID);
    
    %...Clean up
    clear filename fileID
end

%% Interpolate Results to Match Times

%...Set interpolation time
interpolatedTime = simulationTime;

%...Interpolate
if applyInterpolation
    %...Interpolate propagation results
    CartesianPropagatedResults = interp1( simulationTime, CartesianPropagatedResults, interpolatedTime, 'spline' );
    KeplerianPropagatedResults = interp1( simulationTime, KeplerianPropagatedResults, interpolatedTime, 'spline' );
    
    %...Interpolate estimation results
    CartesianEstimatedResults = interp1( onboardTime, CartesianEstimatedResults, interpolatedTime, 'spline' );
    KeplerianEstimatedResults = interp1( onboardTime, KeplerianEstimatedResults, interpolatedTime, 'spline' );
    
    %..Interpolate filter results
    if loadFilter
        filterStateEstimatedResults = interp1( filterTime, filterStateEstimatedResults, interpolatedTime, 'linear', NaN );
        filterCovarianceEstimatedResults = interp1( filterTime, filterCovarianceEstimatedResults, ...
            interpolatedTime, 'linear', NaN );
    end
    
    %...Interpolate other results
%     dependentVariables = interp1( simulationTime, dependentVariables, interpolatedTime, 'spline' );
    if loadMeasurements
        accelerometerMeasurements = interp1( measurementTime, accelerometerMeasurements, interpolatedTime, 'spline' );
        expectedMeasurements = interp1( measurementTime, expectedMeasurements, interpolatedTime, 'spline' );
    end
end

%% Plot 3D Orbit

%...Plot trajectory
F = figure('rend','painters','pos',figSizeLarge);
hold on
plot3(CartesianPropagatedResults(:,1),CartesianPropagatedResults(:,2),CartesianPropagatedResults(:,3),'LineWidth',1.5)
plot3(CartesianEstimatedResults(:,1),CartesianEstimatedResults(:,2),CartesianEstimatedResults(:,3),'LineWidth',1.5)
[x,y,z] = sphere; surf(marsRadius/1e3*x,marsRadius/1e3*y,marsRadius/1e3*z)
hold off
xlabel('x [km]'), ylabel('y [km]'), zlabel('z [km]')
grid on
axis equal tight
set(gca,'FontSize',15)

%% Plot Apses Altitudes

%...Retrieve apo- and periapses altitudes
actualAltitude = sqrt( sum( CartesianPropagatedResults(:,1:3).^2, 2 ) ) - marsRadius/1e3;
actualApoapses = findpeaks(actualAltitude);
actualPeriapses = -findpeaks(-actualAltitude);

estimatedAltitude = sqrt( sum( CartesianEstimatedResults(:,1:3).^2, 2 ) ) - marsRadius/1e3;
estimatedApoapses = findpeaks(estimatedAltitude);
[estimatedPeriapses,periapsisLocations] = findpeaks(-estimatedAltitude);
estimatedPeriapses = -estimatedPeriapses;

maximumDensityPerOrbit = findpeaks(dependentVariables(:,10),'MinPeakHeight',1.0e-13);

%...Plot altitude
F = figure('rend','painters','pos',figSizeSmall);
yyaxis left
hold on
plot(interpolatedTime,actualAltitude,'LineWidth',1.25)
plot(interpolatedTime,estimatedAltitude,'LineWidth',1.25)
hold off
ylabel('Altitude [km]')
yyaxis right
semilogy(interpolatedTime,abs(estimatedAltitude-actualAltitude)*1e3,'LineWidth',1.25)
xlabel(timeLabel)
ylabel('Altitude Difference [m]')
legend('Actual','Estimated','Location','Best')
set(gca,'FontSize',15);
grid on

%...Plot periapses
F = figure('rend','painters','pos',figSizeSmall);
yyaxis left
hold on
scatter(1:length(actualPeriapses),actualPeriapses,'filled','s')
scatter(1:length(estimatedPeriapses),estimatedPeriapses,'filled','d')
% scatter(1:length(targetedPeriapses),targetedPeriapses,'filled')
ylabel('Altitude [km]')
hold off
yyaxis right
scatter(1:length(maximumDensityPerOrbit),maximumDensityPerOrbit,'filled','v')
ylabel('Density [kg m^{-3}]')
set(gca,'YScale','log')
xlabel('Orbit Number [-]')
legend('Actual Periapses','Estimated Periapses','Density')%,'Commanded Periapses')
set(gca,'FontSize',15);
grid on

%% RMS Error

%...Compute RMS error in position and velocity
rmsPositionError = rms( sqrt(sum(CartesianEstimatedResults(:,1:3).^2,2)) - ...
    sqrt(sum(CartesianPropagatedResults(:,1:3).^2,2)) ) * 1e3;
rmsVelocityError = rms( sqrt(sum(CartesianEstimatedResults(:,4:6).^2,2)) - ...
    sqrt(sum(CartesianPropagatedResults(:,4:6).^2,2)) );

%...Show errors
table(rmsPositionError,rmsVelocityError)

%% Plot States Over Time

%...Plot Cartesian translational motion
F = figure('rend','painters','pos',figSizeLarge);
for i = 1:size(CartesianPropagatedResults,2)
    subplot(2,3,i)
    hold on
    plot(interpolatedTime,CartesianPropagatedResults(:,i),'LineWidth',1.25)
    plot(interpolatedTime,CartesianEstimatedResults(:,i),'LineWidth',1.25)
    hold off
    xlabel(timeLabel)
    ylabel(CartesianLabels{i})
    set(gca,'FontSize',15)
    grid on
end
subplotLegend({'Actual','Estimated'})

%...Plot error in Cartesian translational motion
F = figure('rend','painters','pos',figSizeLarge);
for i = 1:size(CartesianPropagatedResults,2)
    subplot(2,3,i)
    plot(interpolatedTime,CartesianEstimatedResults(:,i)-CartesianPropagatedResults(:,i),'LineWidth',1.25)
    xlabel(timeLabel)
    ylabel(CartesianLabelsDifference{i},'LineWidth',1.25)
    set(gca,'FontSize',15)
    grid on
end

%...Plot Keplerian translational motion
F = figure('rend','painters','pos',figSizeLarge);
for i = 1:size(KeplerianPropagatedResults,2)
    subplot(2,3,i)
    hold on
    plot(interpolatedTime,KeplerianPropagatedResults(:,i),'LineWidth',1.25)
    plot(interpolatedTime,KeplerianEstimatedResults(:,i),'LineWidth',1.25)
    hold off
    xlabel(timeLabel)
    ylabel(KeplerianLabels{i})
    set(gca,'FontSize',15)
    grid on
end
subplotLegend({'Actual','Estimated'})

%...Plot error in Keplerian translational motion
F = figure('rend','painters','pos',figSizeLarge);
for i = 1:size(KeplerianPropagatedResults,2)
    subplot(2,3,i)
    plot(interpolatedTime,KeplerianEstimatedResults(:,i)-KeplerianPropagatedResults(:,i),'LineWidth',1.25)
    xlabel(timeLabel)
    ylabel(KeplerianLabels{i})
    set(gca,'FontSize',15)
    grid on
end

% %...Plot rotational motion
% F = figure('rend','painters','pos',figSizeLarge);
% for i = size(rotationalPropagatedResults,2):-1:1
%     subplot(2,4,i)
%     hold on
%     plot(interpolatedTime,rotationalPropagatedResults(:,i),'LineWidth',1.25)
%     plot(interpolatedTime,rotationalEstimatedResults(:,i),'LineWidth',1.25)
%     if i <= 4
%         plot(onboardTime,commandedRotationalState(:,i),'LineWidth',1.25)
%     end
%     hold off
%     xlabel(timeLabel)
%     ylabel(rotationLabels{i})
%     set(gca,'FontSize',15)
%     grid on
% end
% subplotLegend({'Actual','Estimated','Commanded'})

%% Plot Filter States

%...Only if filtering is toggled
if loadFilter
    %...Plot Cartesian translational motion
    F = figure('rend','painters','pos',figSizeLarge);
    for i = 1:6
        subplot(2,3,i)
        hold on
        plot(interpolatedTime,filterStateEstimatedResults(:,i)-CartesianPropagatedResults(:,i),'LineWidth',1.25)
        plot(interpolatedTime(2:end),sqrt(filterCovarianceEstimatedResults(2:end,i)),'LineWidth',1.25,'LineStyle','--')
        plot(interpolatedTime(2:end),-sqrt(filterCovarianceEstimatedResults(2:end,i)),'LineWidth',1.25,'LineStyle','--')
        hold off
        xlabel(timeLabel)
        set(gca,'FontSize',15)
        grid on
    end
    
    %...Plot instrument errors
    actual = [-0.000123974  5.23187e-06   9.7221e-06 -4.07472e-05 -0.000166043  -8.9237e-05];
    F = figure('rend','painters','pos',figSizeLarge);
    for i = 7:12
        subplot(2,3,i-6)
        hold on
        plot(interpolatedTime,filterStateEstimatedResults(:,i)-actual(i-6),'LineWidth',1.25)
        plot(interpolatedTime(2:end),sqrt(filterCovarianceEstimatedResults(2:end,i)),'LineWidth',1.25,'LineStyle','--')
        plot(interpolatedTime(2:end),-sqrt(filterCovarianceEstimatedResults(2:end,i)),'LineWidth',1.25,'LineStyle','--')
        hold off
        xlabel(timeLabel)
        set(gca,'FontSize',15)
        grid on
    end
end

%% Plot Other Results

%...Only if measurements are toggled
if loadMeasurements
    %...Plot measurements
    F = figure('rend','painters','pos',figSizeLarge);
    
    subplot(1,2,1)
    hold on
    for i = 1:3
        plot(interpolatedTime,accelerometerMeasurements(:,i),'LineWidth',1.0)
        plot(interpolatedTime,smooth(accelerometerMeasurements(:,i),25),'LineWidth',1.25)
        plot(interpolatedTime,expectedMeasurements(:,i),'LineWidth',1.25,'LineStyle','--')
        plot(interpolatedTime,dependentVariables(:,i+3),'LineWidth',1.25,'LineStyle',':')
    end
    hold off
    xlabel(timeLabel)
    ylabel('Translational Acceleration [m s^{-1}]')
    set(gca,'FontSize',15)
    grid on
    legend('x_{IMU}','x_{smooth}','x_{exp}','x_{real}',...
        'y_{IMU}','y_{smooth}','y_{exp}','y_{real}',...
        'z_{IMU}','z_{smooth}','z_{exp}','z_{real}')
    % legend('x_{IMU}','x_{real}','y_{IMU}','y_{real}','z_{IMU}','z_{real}')
    
    subplot(1,2,2)
    plot(interpolatedTime,accelerometerMeasurements(:,4:6),'LineWidth',1.25)
    hold off
    xlabel(timeLabel)
    ylabel('Rotational Velocity [rad s^{-1}]')
    set(gca,'FontSize',15)
    grid on
    legend('x','y','z')
end

% %...Plot control torques
% figure;
% plot(simulationTime,controlTorques(:,1:3),'LineWidth',1.25)
% xlabel(timeLabel)
% ylabel('Control Torque [N m]')
% set(gca,'FontSize',15)
% grid on
% legend('x','y','z','Location','Best')

%...Plot aerodynamic angles
figure;
plot(interpolatedTime,dependentVariables(:,7:9),'LineWidth',1.25)
xlabel(timeLabel)
ylabel('Angle [deg]')
set(gca,'FontSize',15)
grid on
legend('Attack','Side-slip','Bank','Location','Best')

% %...Plot aerodynamic acceleration
% figure;
% ax = polaraxes;
% polarplot(deg2rad(KeplerianPropagatedResults(:,6)),dependentVariables(:,10),'LineWidth',1.25)
% grid on
% ax.RAxis.Label.String = 'Aerodynamic Acceleration [m s^{-2}]';
% ax.ThetaAxis.Label.String = 'True Anomaly [deg]';
% set(gca,'FontSize',15)
% 
% figure;
% plot(KeplerianPropagatedResults(:,6),dependentVariables(:,10),'LineWidth',1.25)
% grid on
% xlabel('True Anomaly [deg]')
% ylabel('Aerodynamic Acceleration [m s^{-2}]')
% set(gca,'FontSize',15)

%...Plot heating conditions
figure;
yyaxis left
hold on
plot(interpolatedTime,dependentVariables(:,11),'LineWidth',1.25)
plot([interpolatedTime(1),interpolatedTime(end)],[0.19,0.19],'LineWidth',1.25,'LineStyle','--')
set(gca,'YScale','log')
hold off
ylabel('Dynamic Pressure [kg m^{-2}]')
ylims = ylim; ylim([1e-6,ylims(2)]);
yyaxis right
hold on
plot(interpolatedTime,dependentVariables(:,12),'LineWidth',1.25)
plot([interpolatedTime(1),interpolatedTime(end)],[2800,2800],'LineWidth',1.25,'LineStyle','--')
set(gca,'YScale','log')
hold off
ylabel('Heat Rate [W m^{-2}]')
ylims = ylim; ylim([1e-2,ylims(2)]);
xlabel(timeLabel)
set(gca,'FontSize',15)
grid on

%% Plot Density

%...Only if density is toggled
if loadMeasurements
    %...Find pericenter height
    estimatedAltitude = sqrt( sum( CartesianEstimatedResults(:,1:3).^2, 2 ) ) - marsRadius / 1e3;
    pericenter = min( estimatedAltitude );
    
    %...Retireve density
    aerodynamicAcceleration = sqrt( sum( accelerometerMeasurements(:,1:3).^2, 2 ) );
    atmosphericDensity = 2 * 1000 / 37.5 / 1.87 ./ ...
        sum( ( CartesianEstimatedResults(:,4:6) ).^2, 2 ) .* aerodynamicAcceleration;
    
    %...Data below atmospheric interface
%     loc = ( estimatedAltitude <= marsAtmosphericInterface );
%     loc(find(estimatedAltitude==pericenter)+1:end) = false;
    atmospherePhaseTime = simulationTime;%simulationTime( loc );
    % estimatedAltitude = estimatedAltitude( loc );
    % atmosphericDensity = atmosphericDensity( loc( 1:end-1 ) );
    
    %...Plot density
    figure;
    hold on
    plot(atmospherePhaseTime,atmosphericDensity,'LineWidth',1.25)
%     plot(atmospherePhaseTime,dependentVariables(:,10),'LineWidth',1.25)
    hold off
    xlabel(timeLabel)
    ylabel('Density [kg m^{-3}]')
    set(gca,'FontSize',15,'YScale','log')
    grid on
    % legend('Measured','Actual','Location','Best')
end

%%

x = [133653 1.75013e-09     5652.26   -0.626817   0.0248434   0.0492727];
x = [138087 7.85814e-10     3945.28   -0.585279  0.00602754  0.00572154];
x = [137946 9.55613e-10     1910.18   -0.207692    0.007515  0.00841763];
x = [137350 7.59367e-10      1860.1    -0.27187  0.00221564   0.0042782];
dens_func = @(h) x(2) * exp( x(4) * ( h*1e3 - x(1) ) / x(3) + x(5) * cos( 2*pi*( ( h*1e3 - x(1) ) / x(3) ) ) + ...
    x(6) * sin( 2*pi*( ( h*1e3 - x(1) ) / x(3) ) ) );

figure
hold on
% plot(dens_func(125:1000),125:1000,'LineWidth',1.25)
plot(dependentVariables(:,10),actualAltitude,'LineWidth',1.25)
plot(dens_func(actualAltitude),actualAltitude,'LineWidth',1.25)
hold off
xlabel('Density [kg m^{-3}]')
ylabel('Altitude [km]')
ylim([125,175])
grid on
set(gca,'FontSize',15,'XScale','log')

%% Clean Up

if ~showFigure, close all, end