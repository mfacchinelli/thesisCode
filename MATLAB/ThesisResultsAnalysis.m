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
onboardAtmosphereModel = 0; % 0: exp, 1: wave, 2: 5-param

%...Output folder
outputFolder = 'SimulationOutput';
% outputFolder = [outputFolder,'TransOnly'];
% outputFolder = [outputFolder,'TransOnlyReduced'];
% outputFolder = [outputFolder,'TransOnly/server/without_force'];
% outputFolder = [outputFolder,'TransGuidOnly'];
% outputFolder = [outputFolder,'TransOnlyIMAN/high_ecc'];
outputFolder = [outputFolder,'TransOnlyIMANRMS/high_ecc/7/10'];

%...Plot settings
loadRotational = false;
loadDependent = false;
loadFilter = true;
loadMeasurements = false;
applyInterpolation = true;

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

%...Load rotational motion
if loadRotational
    filename = ['/Users/Michele/GitHub/tudat/tudatBundle/tudatApplications/Thesis/',outputFolder,'/rotationalPropagated.dat'];
    fileID = fopen(filename,'r');
    rotationalPropagatedResults = textscan(fileID,repmat('%f',[1,8]),'Delimiter',',','CollectOutput',true);
    rotationalPropagatedResults = rotationalPropagatedResults{1}(:,2:end);
    rotationalPropagatedResults(:,6:8) = rad2deg(rotationalPropagatedResults(:,5:7));
    rotationalPropagatedResults(:,5) = 1.0 - quatnorm(rotationalPropagatedResults(:,1:4));
    fclose(fileID);
end

%...Clean up
clear filename fileID

%% Load C++ Results For Estimation

%...Load translational motion
filename = ['/Users/Michele/GitHub/tudat/tudatBundle/tudatApplications/Thesis/',outputFolder,'/cartesianEstimated.dat'];
fileID = fopen(filename,'r');
CartesianEstimatedResults = textscan(fileID,repmat('%f',[1,7]),'Delimiter',',','CollectOutput',true);
onboardTime = ( CartesianEstimatedResults{1}(:,1) - initialTime ) / timeConversion; 
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

%...Load rotational motion
if loadRotational
    filename = ['/Users/Michele/GitHub/tudat/tudatBundle/tudatApplications/Thesis/',outputFolder,'/rotationalEstimated.dat'];
    fileID = fopen(filename,'r');
    rotationalEstimatedResults = textscan(fileID,repmat('%f',[1,8]),'Delimiter',',','CollectOutput',true);
    rotationalEstimatedResults = rotationalEstimatedResults{1}(:,2:end);
    rotationalEstimatedResults(:,6:8) = rad2deg(rotationalEstimatedResults(:,5:7));
    rotationalEstimatedResults(:,5) = 1.0 - quatnorm(rotationalEstimatedResults(:,1:4));
    fclose(fileID);
end

%...Clean up
clear filename fileID

%% Load C++ Results For Navigation Filter

%...Only if filtering is toggled
if loadFilter
    %...Load translational motion
    filename = ['/Users/Michele/GitHub/tudat/tudatBundle/tudatApplications/Thesis/',outputFolder,'/filterStateEstimates.dat'];
    fileID = fopen(filename,'r');
    filterStateEstimatedResults = textscan(fileID,repmat('%f',[1,11]),'Delimiter',',','CollectOutput',true);
    filterTime = ( filterStateEstimatedResults{1}(:,1) - initialTime ) / timeConversion;
    filterStateEstimatedResults = filterStateEstimatedResults{1}(:,2:end);
    filterStateEstimatedResults(:,1:3) = filterStateEstimatedResults(:,1:3)/1e3;
    fclose(fileID);
    
    %...Load Keplerian translational motion
    filename = ['/Users/Michele/GitHub/tudat/tudatBundle/tudatApplications/Thesis/',outputFolder,...
        '/filterCovarianceEstimates.dat'];
    fileID = fopen(filename,'r');
    filterCovarianceEstimatedResults = textscan(fileID,repmat('%f',[1,11]),'Delimiter',',','CollectOutput',true);
    filterCovarianceEstimatedResults = filterCovarianceEstimatedResults{1}(:,2:end);
    filterCovarianceEstimatedResults(:,1:3) = filterCovarianceEstimatedResults(:,1:3)/1e3;
    fclose(fileID);
    
    %...Remove first entry
    filterTime = filterTime(2:end);
    filterStateEstimatedResults = filterStateEstimatedResults(2:end,:);
    filterCovarianceEstimatedResults = filterCovarianceEstimatedResults(2:end,:);
    
    %...Clean up
    clear filename fileID
end

%% Load C++ Results For Dependent Variables

if loadDependent
    %...Load dependent variables
    filename = ['/Users/Michele/GitHub/tudat/tudatBundle/tudatApplications/Thesis/',outputFolder,'/dependentVariables.dat'];
    fileID = fopen(filename,'r');
    dependentVariables = textscan(fileID,repmat('%f',[1,13]),'Delimiter',',','CollectOutput',true);
    dependentVariables = dependentVariables{1}(:,2:end); dependentVariables(:,7:9) = rad2deg(dependentVariables(:,7:9));
    fclose(fileID);
    
    %...Load control torques
    if loadRotational && loadMeasurements
        filename = ['/Users/Michele/GitHub/tudat/tudatBundle/tudatApplications/Thesis/',outputFolder,'/controlTorques.dat'];
        fileID = fopen(filename,'r');
        controlTorques = textscan(fileID,repmat('%f',[1,4]),'Delimiter',',','CollectOutput',true);
        controlTorques = controlTorques{1}(:,2:end);
        fclose(fileID);
    end

    %...Load periapsis corridors
    filename = ['/Users/Michele/GitHub/tudat/tudatBundle/tudatApplications/Thesis/',outputFolder,'/periapsisCorridors.dat'];
    fileID = fopen(filename,'r');
    periapsisCorridors = textscan(fileID,repmat('%f',[1,3]),'Delimiter',',','CollectOutput',true);
    periapsisCorridors = periapsisCorridors{1};
    fclose(fileID);

    %...Load apopasis maneuvers
    filename = ['/Users/Michele/GitHub/tudat/tudatBundle/tudatApplications/Thesis/',outputFolder,'/apoapsisManeuver.dat'];
    fileID = fopen(filename,'r');
    apoapsisManeuver = textscan(fileID,repmat('%f',[1,2]),'Delimiter',',','CollectOutput',true);
    apoapsisManeuver = apoapsisManeuver{1};
    fclose(fileID);
    
    %...Load atmosphere parameters
    if onboardAtmosphereModel == 0, limit = 4; else, limit = 7; end
    filename = ['/Users/Michele/GitHub/tudat/tudatBundle/tudatApplications/Thesis/',outputFolder,'/atmosphericParameters.dat'];
    fileID = fopen(filename,'r');
    atmosphericParameters = textscan(fileID,repmat('%f',[1,limit]),'Delimiter',',','CollectOutput',true);
    atmosphericParameters = atmosphericParameters{1};
    fclose(fileID);
    
    %...Clean up
    clear filename fileID limit
else
    dependentVariables = NaN(length(simulationTime),12);
end

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
interpolatedTime = onboardTime(1:end-1);

%...Interpolate
if applyInterpolation
    %...Interpolate propagation results
    CartesianPropagatedResults = interp1( simulationTime, CartesianPropagatedResults, interpolatedTime, 'linear' );
    KeplerianPropagatedResults = interp1( simulationTime, KeplerianPropagatedResults, interpolatedTime, 'linear' );
    if loadRotational
        rotationalPropagatedResults = interp1( simulationTime, rotationalPropagatedResults, interpolatedTime, 'linear' );
    end
    
    %...Interpolate estimation results
    CartesianEstimatedResults = interp1( onboardTime, CartesianEstimatedResults, interpolatedTime, 'linear' );
    KeplerianEstimatedResults = interp1( onboardTime, KeplerianEstimatedResults, interpolatedTime, 'linear' );
    if loadRotational
        rotationalEstimatedResults = interp1( onboardTime, rotationalEstimatedResults, interpolatedTime, 'linear' );
    end
    
    %..Interpolate filter results
    if loadFilter
        filterStateEstimatedResults = interp1( filterTime, filterStateEstimatedResults, interpolatedTime, 'linear', NaN );
        filterCovarianceEstimatedResults = interp1( filterTime, filterCovarianceEstimatedResults, ...
            interpolatedTime, 'linear', NaN );
    end
    
    %...Interpolate other results
    if loadDependent
        dependentVariables = interp1( simulationTime, dependentVariables, interpolatedTime, 'linear' );
    end
    if loadMeasurements
        accelerometerMeasurements = interp1( measurementTime, accelerometerMeasurements, interpolatedTime, 'linear' );
        expectedMeasurements = interp1( measurementTime, expectedMeasurements, interpolatedTime, 'linear' );
    end
end

%% Plot 3D Orbit

%...Plot trajectory
F = figure('rend','painters','pos',figSizeLarge);
hold on
plot3(CartesianPropagatedResults(1:100:end,1),CartesianPropagatedResults(1:100:end,2),...
    CartesianPropagatedResults(1:100:end,3),'LineWidth',1.5)
% plot3(CartesianEstimatedResults(1:100:end,1),CartesianEstimatedResults(1:100:end,2),...
%     CartesianEstimatedResults(1:100:end,3),'LineWidth',1.5)
[x,y,z] = sphere; surf(marsRadius/1e3*x,marsRadius/1e3*y,marsRadius/1e3*z)
hold off
xlabel('x [km]'), ylabel('y [km]'), zlabel('z [km]')
view([0,0])
grid on
axis equal tight
set(gca,'FontSize',15)

%% Plot Apses Altitudes

%...Retrieve apo- and periapses altitudes
actualAltitude = sqrt( sum( CartesianPropagatedResults(:,1:3).^2, 2 ) ) - marsRadius/1e3;
actualApoapses = findpeaks(actualAltitude);
actualPeriapses = -findpeaks(-actualAltitude); %actualPeriapses(8) = [];

estimatedAltitude = sqrt( sum( CartesianEstimatedResults(:,1:3).^2, 2 ) ) - marsRadius/1e3;
estimatedApoapses = findpeaks(estimatedAltitude);
[estimatedPeriapses,periapsisLocations] = findpeaks(-estimatedAltitude);
estimatedPeriapses = -estimatedPeriapses;

% %...Plot altitude
% F = figure('rend','painters','pos',figSizeSmall);
% yyaxis left
% hold on
% plot(interpolatedTime,actualAltitude,'LineWidth',1.25)
% plot(interpolatedTime,estimatedAltitude,'LineWidth',1.25)
% hold off
% ylabel('Altitude [km]')
% yyaxis right
% semilogy(interpolatedTime,abs(estimatedAltitude-actualAltitude)*1e3,'LineWidth',1.25)
% xlabel(timeLabel)
% ylabel('Altitude Difference [m]')
% legend('Actual','Estimated','Location','Best')
% set(gca,'FontSize',15);
% grid on

%...Plot periapses
if loadDependent
    F = figure('rend','painters','pos',figSizeSmall);
    hold on
    scatter(135+(1:length(actualPeriapses)),actualPeriapses,50)
    plot(periapsisCorridors(2:end,1),periapsisCorridors(1:end-1,2)/1e3,'LineWidth',1.25,'LineStyle','--')
    plot(periapsisCorridors(2:end,1),periapsisCorridors(1:end-1,3)/1e3,'LineWidth',1.25,'LineStyle','--')
    hold off
    xlabel('Orbit Number [-]')
    ylabel('Altitude [km]')
    legend('Actual Periapses','Lower Bound','Upper Bound','Location','Best')
    set(gca,'FontSize',15)
    grid on
end

%% Plot Densities

if loadDependent
    %...Find maximum density places
    [maximumDensityPerOrbit,locMaxDens] = findpeaks(dependentVariables(:,10),'MinPeakHeight',1.0e-13);
    estimatedMaxDensity = zeros(size(maximumDensityPerOrbit));
    
    %...Plot real and estimated densities
    figure
    hold on
    loc = actualAltitude < 175;
    reducedAltitude = actualAltitude(loc);
    plot(dependentVariables(loc,10),reducedAltitude,'LineWidth',1.25)
    
    i = 0;
    for x = atmosphericParameters(:,2:end)'
        i = i+1;
        if onboardAtmosphereModel == 0
            dens_func = @(h) x(2) * exp( - ( h*1e3 - x(1) ) / x(3) );
        else
            dens_func = @(h) x(2) * exp( x(4) * ( h*1e3 - x(1) ) / x(3) + x(5) * cos( 2*pi*( ( h*1e3 - x(1) ) / x(3) ) ) + ...
                x(6) * sin( 2*pi*( ( h*1e3 - x(1) ) / x(3) ) ) );
        end
%         estimatedMaxDensity(i) = dens_func(actualAltitude(locMaxDens(i)));
        plot(dens_func(reducedAltitude),reducedAltitude,'LineWidth',1.25)
    end
    hold off
    xlabel('Density [kg m^{-3}]')
    ylabel('Altitude [km]')
    ylim([125,175])
    grid on
    set(gca,'FontSize',15,'XScale','log')
    
    %...Plot real and estimated peak densities
    figure
    hold on
    scatter(135+(1:length(maximumDensityPerOrbit)),maximumDensityPerOrbit,50)
%     scatter(atmosphericParameters(:,1),estimatedMaxDensity,50,'s')    
    hold off
    ylabel('Density [kg m^{-3}]')
    xlabel('Orbit Number [-]')
    set(gca,'YScale','log')
    legend('Actual','Estimated')
    set(gca,'FontSize',15)
    grid on
end

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

%...Plot rotational motion
if loadRotational
    F = figure('rend','painters','pos',figSizeLarge);
    for i = size(rotationalPropagatedResults,2):-1:1
        subplot(2,4,i)
        hold on
        plot(interpolatedTime,rotationalPropagatedResults(:,i),'LineWidth',1.25)
        plot(interpolatedTime,rotationalEstimatedResults(:,i),'LineWidth',1.25)
        hold off
        xlabel(timeLabel)
        ylabel(rotationLabels{i})
        set(gca,'FontSize',15)
        grid on
    end
    subplotLegend({'Actual','Estimated'})
end

%% Plot Peri- and Apoapsis Convergence

%...Compute apsis altitudes
apoapsisAltitude = KeplerianPropagatedResults(:,1) .* ...
    (1.0 + KeplerianPropagatedResults(:,2)) - marsRadius/1e3;
periapsisAltitude = KeplerianPropagatedResults(:,1) .* ...
    (1.0 - KeplerianPropagatedResults(:,2)) - marsRadius/1e3;

%...Compute DAIA
percentage = 0.5 * ( 1 + ( 1 - KeplerianPropagatedResults(:,1) / 26021000.0 ) );
dynamicAtmosphericInterfaceAltitude = percentage .* ( KeplerianPropagatedResults(:,1) - marsRadius/1e3 );

%...Plot convergence
figure
hold on
plot(interpolatedTime,apoapsisAltitude,'LineWidth',1.25)
plot(interpolatedTime,periapsisAltitude,'LineWidth',1.25)
plot(interpolatedTime,sqrt(sum(CartesianPropagatedResults(:,1:3).^2,2)) - marsRadius/1e3,'LineWidth',1.25,'LineStyle','--')
plot(interpolatedTime,dynamicAtmosphericInterfaceAltitude,'LineWidth',1.25,'LineStyle',':')
hold off
xlabel(timeLabel)
ylabel('Altitude [km]')
legend('Apoapsis','Periapsis','Altitude','DAIA','Location','Best')
grid on
set(gca,'FontSize',15,'YScale','log')

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
    actual = [-0.000123974  5.23187e-06   9.7221e-06];
    F = figure('rend','painters','pos',figSizeLarge);
    for i = 7:9
        subplot(1,3,i-6)
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
if loadMeasurements && loadDependent
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

%...Plot control torques
if loadRotational && loadMeasurements
    figure;
    plot(simulationTime,controlTorques(:,1:3),'LineWidth',1.25)
    xlabel(timeLabel)
    ylabel('Control Torque [N m]')
    set(gca,'FontSize',15)
    grid on
    legend('x','y','z','Location','Best')
end

%...Plot aerodynamic angles
if loadDependent
    figure;
    plot(interpolatedTime,dependentVariables(:,7:9),'LineWidth',1.25)
    xlabel(timeLabel)
    ylabel('Angle [deg]')
    set(gca,'FontSize',15)
    grid on
    legend('Attack','Side-slip','Bank','Location','Best')
end

% %...Plot aerodynamic acceleration
% if loadDependent
%     figure;
%     ax = polaraxes;
%     polarplot(deg2rad(KeplerianPropagatedResults(:,6)),dependentVariables(:,10),'LineWidth',1.25)
%     grid on
%     ax.RAxis.Label.String = 'Aerodynamic Acceleration [m s^{-2}]';
%     ax.ThetaAxis.Label.String = 'True Anomaly [deg]';
%     set(gca,'FontSize',15)
%     
%     figure;
%     plot(KeplerianPropagatedResults(:,6),dependentVariables(:,10),'LineWidth',1.25)
%     grid on
%     xlabel('True Anomaly [deg]')
%     ylabel('Aerodynamic Acceleration [m s^{-2}]')
%     set(gca,'FontSize',15)
% end

%...Plot heating conditions
if loadDependent
    %...Find peak dynamic pressure
    [maximumDynamicPressure,locDynPressPeak] = findpeaks(dependentVariables(:,11));
    loc = maximumDynamicPressure > 1e-2;
    maximumDynamicPressure = maximumDynamicPressure(loc);
    locDynPressPeak = locDynPressPeak(loc);
    
    %...Find peak heat rate
    [maximumHeatRate,locHeatRatePeak] = findpeaks(dependentVariables(:,12));
    loc = maximumHeatRate > 1e-2;
    maximumHeatRate = maximumHeatRate(loc);
    locHeatRatePeak = locHeatRatePeak(loc);
    
    %...Compute peak heat load
    locHeatLoadPeak = locHeatRatePeak;
    maximumHeatLoad = zeros(size(maximumHeatRate));
    for i = 1:length(locHeatLoadPeak)
        if i == 1
            l = locHeatLoadPeak(i)/2;
        else
            l = (locHeatLoadPeak(i-1) + locHeatLoadPeak(i))/2;
        end
        if i == length(locHeatLoadPeak)
            u = (locHeatLoadPeak(i) + length(interpolatedTime))/2;
        else
            u = (locHeatLoadPeak(i) + locHeatLoadPeak(i+1))/2;
        end
        l = floor(l); u = floor(u);
        maximumHeatLoad(i) = trapz(interpolatedTime(l:u)*timeConversion,dependentVariables(l:u,12)) / 1e3;
    end
    
    %...Plot data
    F = figure;
    
    subplot(1,2,1)
    hold on
    scatter(interpolatedTime(locDynPressPeak),maximumDynamicPressure,50)
    plot([interpolatedTime(1),interpolatedTime(end)],[0.19,0.19],'LineWidth',1.25,'LineStyle','--')
    hold off
    xlabel(timeLabel)
    ylabel('Dynamic Pressure [kg m^{-2}]')
    legend('Dyn. Press.','Min. Dyn. Press.','Location','NW')
    set(gca,'FontSize',15)
    grid on
    
%     hold on
%     for maneuver = apoapsisManeuver'
%         if maneuver(1) == 0
%             time = interpolatedTime(1);
%         elseif maneuver(1) == length(locHeatRatePeak)
%             loc = floor((locHeatRatePeak(maneuver(1)) + length(interpolatedTime))/2);
%             time = interpolatedTime(loc);
%         else
%             loc = floor((locHeatRatePeak(maneuver(1)) + locHeatRatePeak(maneuver(1)+1))/2);
%             time = interpolatedTime(loc);
%         end
%         plot([time,time],ylim,'k','LineWidth',0.5,'LineStyle','-','Marker','none')
%     end
%     hold off
    
    subplot(1,2,2)
    yyaxis left
    hold on
    scatter(interpolatedTime(locHeatRatePeak),maximumHeatRate,50)
    plot([interpolatedTime(1),interpolatedTime(end)],[2800,2800],'LineWidth',1.25,'LineStyle','--')
    hold off
    ylabel('Heat Rate [W m^{-2}]')
    
    yyaxis right
    hold on
    scatter(interpolatedTime(locHeatLoadPeak),maximumHeatLoad,50,'s')
    plot([interpolatedTime(1),interpolatedTime(end)],[500,500],'LineWidth',1.25,'LineStyle','-.')
    hold off
    ylabel('Heat Load [kJ m^{-2}]')
    xlabel(timeLabel)
    legend('Heat Rate','Max. Heat Rate','Heat Load','Max. Heat Load','Location','NW')
    set(gca,'FontSize',15)
    grid on
end

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
    atmospherePhaseTime = interpolatedTime;
    
    %...Plot density
    figure;
    hold on
    plot(atmospherePhaseTime,atmosphericDensity,'LineWidth',1.25)
    plot(atmospherePhaseTime,dependentVariables(:,10),'LineWidth',1.25)
    hold off
    xlabel(timeLabel)
    ylabel('Density [kg m^{-3}]')
    ylim([1e-11,1e-8])
    grid on
    legend('Measured','Actual','Location','Best')
    set(gca,'FontSize',15,'YScale','log')
    
    %...Plot altitude vs. acceleration
    F = figure('rend','painters','pos',figSizeSmall);
    hold on
    scatter(estimatedAltitude,aerodynamicAcceleration,'LineWidth',1.25)
    plot(estimatedAltitude,sqrt( sum( dependentVariables(:,4:6).^2, 2 ) ),'LineWidth',1.25)
    plot([150,150],[1e-5,1e-2],'LineWidth',1.25,'LineStyle','--')
    hold off
    xlabel('Altitude [km]')
    ylabel('Acceleration [kg m^{-3}]')
    xlim([100,250])
    ylim([1e-5,1e-2])
    grid on
    legend('Noisy','Ideal','Reduced Interface','Location','NE')
    set(gca,'FontSize',15,'YScale','log')
    if saveFigure, saveas(F,'../../Report/figures/imu_acc_raia','epsc'), end
end

%% Plot Aerobraking Evolution

if saveFigure && strcmp(outputFolder,'SimulationOutputTransGuidOnly/')
    %...Determine change in apoapsis altitude drop rate
    altitudes = sqrt(sum(CartesianPropagatedResults(:,1:3).^2,2)) - marsRadius/1e3;
    [apoapsesAltitudes,apoapsesLocs] = findpeaks(altitudes);
    apoapsesAltitudes = [altitudes(1);apoapsesAltitudes];
    apoapsesTimes = [0;interpolatedTime(apoapsesLocs)];
    
    %...Reduce to exclude last day
    loc = apoapsesTimes < (interpolatedTime(end)-1); loc(find(loc==0,1)-1) = false;
    apoapsesTimes = apoapsesTimes(loc);
    apoapsesAltitudes = apoapsesAltitudes(loc);
    
    apoapsesPercentDiff = diff(apoapsesAltitudes)./apoapsesAltitudes(1:end-1) * 100;
    
    loc = interpolatedTime <= apoapsesTimes(end);
    plotTime = interpolatedTime(loc);
    plotEccentricity = KeplerianPropagatedResults(loc,2);

    %...Plot results
    F = figure('rend','painters','pos',figSizeSmall);
    yyaxis left
    plot(plotTime,plotEccentricity,'LineWidth',1.25)
    ylabel(KeplerianLabels{2})
    yyaxis right
    scatter(apoapsesTimes(2:end),apoapsesPercentDiff,50)
    ylabel('Apoapsis Altitude Change [%]')
    xlabel(timeLabel)
    grid on
    legend('Eccentricity','Apoapsis Altitude Change','Location','SW')
    set(gca,'FontSize',15)
    saveas(F,'../../Report/figures/aerobrake_cont_evol','epsc')
end

%% Clean Up

if ~showFigure, close all, end