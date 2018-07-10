fclose('all'); clear all; close all force; profile off; clc; format long g; rng default;
addpath functions

%% Settings

%...Figure settings
showFigure = true;
saveFigure = false;
[figSizeLarge,figSizeMedium,figSizeSmall] = saveFigureSettings(saveFigure);

%...Constants
marsRadius = 3.396e6;
marsGravitationalParameter = 4.282e13;
marsAtmosphericInterface = 175;

%...Labels
timeConversion = 3600 * 24;
timeLabel = 'Time [d]';
CartesianLabels = {'x [km]','y [km]','z [km]','v_x [km s^{-1}]','v_y [km s^{-1}]','v_z [km s^{-1}]'};
CartesianLabelsDifference = {'\Delta x [km]','\Delta y [km]','\Delta z [km]',...
    '\Delta v_x [km s^{-1}]','\Delta v_y [km s^{-1}]','\Delta v_z [km s^{-1}]'};
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
rotationalPropagatedResults(:,5) = 1.0 - quatnorm(rotationalPropagatedResults(:,1:4));
% [sigma,beta,alpha] = quat2angle(rotationalPropagatedResults(:,1:4),'XZY');
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
rotationalEstimatedResults(:,5) = 1.0 - quatnorm(rotationalEstimatedResults(:,1:4));
fclose(fileID);

%...Clean up
clear filename fileID

%...Compute commanded directon cosine matrix
commandedDirectionCosineMatrix = zeros(3,3,length(onboardTime));
commandedDirectionCosineMatrix(:,1,:) = ( CartesianEstimatedResults(:,4:6) ./ ...
    sqrt( sum( CartesianEstimatedResults(:,4:6).^2, 2 ) ) )';
commandedDirectionCosineMatrix(:,3,:) = ( CartesianEstimatedResults(:,1:3) ./ ...
    sqrt( sum( CartesianEstimatedResults(:,1:3).^2, 2 ) ) )';
for i = 1:length(onboardTime)
    commandedDirectionCosineMatrix(:,3,i) = commandedDirectionCosineMatrix(:,3,i) - ...
        dot( commandedDirectionCosineMatrix(:,3,i), commandedDirectionCosineMatrix(:,1,i) ) * ...
        commandedDirectionCosineMatrix(:,1,i);
    commandedDirectionCosineMatrix(:,2,i) = cross( commandedDirectionCosineMatrix(:,3,i), ...
        commandedDirectionCosineMatrix(:,1,i) );
    commandedDirectionCosineMatrix(:,:,i) = commandedDirectionCosineMatrix(:,:,i)';
end

%...Compute commanded quaternions
commandedRotationalState = dcm2quat(commandedDirectionCosineMatrix);

%% Load C++ Results For Dependent Variables

%...Load dependent variables
filename = '/Users/Michele/GitHub/tudat/tudatBundle/tudatApplications/Thesis/SimulationOutput/dependentVariables.dat';
fileID = fopen(filename,'r');
dependentVariables = textscan(fileID,repmat('%f',[1,11]),'Delimiter',',','CollectOutput',true);
dependentVariables = dependentVariables{1}(:,2:end); dependentVariables(:,7:9) = rad2deg(dependentVariables(:,7:9));
fclose(fileID);

%...Load control torques
filename = '/Users/Michele/GitHub/tudat/tudatBundle/tudatApplications/Thesis/SimulationOutput/controlTorques.dat';
fileID = fopen(filename,'r');
controlTorques = textscan(fileID,repmat('%f',[1,4]),'Delimiter',',','CollectOutput',true);
controlTorques = controlTorques{1}(:,2:end);
fclose(fileID);

%...Clean up
clear filename fileID

%% Load C++ Results For Measurements

%...Load IMU measurements
filename = '/Users/Michele/GitHub/tudat/tudatBundle/tudatApplications/Thesis/SimulationOutput/inertialMeasurements.dat';
fileID = fopen(filename,'r');
inertialMeasurements = textscan(fileID,repmat('%f',[1,7]),'Delimiter',',','CollectOutput',true);
measurementTime = ( inertialMeasurements{1}(:,1) - inertialMeasurements{1}(1) ) / timeConversion; 
inertialMeasurements = inertialMeasurements{1}(:,2:end);
fclose(fileID);

%...Load expected measurements
filename = '/Users/Michele/GitHub/tudat/tudatBundle/tudatApplications/Thesis/SimulationOutput/expectedlMeasurements.dat';
fileID = fopen(filename,'r');
expectedMeasurements = textscan(fileID,repmat('%f',[1,4]),'Delimiter',',','CollectOutput',true);
expectedMeasurements = expectedMeasurements{1}(:,2:end);
fclose(fileID);

%...Clean up
clear filename fileID

%% Plot 3D Orbit

%...Plot trajectory
F = figure('rend','painters','pos',figSizeLarge);
hold on
quiver3(CartesianPropagatedResults(1,1),CartesianPropagatedResults(1,2),CartesianPropagatedResults(1,3),...
    commandedDirectionCosineMatrix(1,1,1),commandedDirectionCosineMatrix(1,2,1),commandedDirectionCosineMatrix(1,3,1),...
    'AutoScaleFactor',1e4,'LineWidth',1)
quiver3(CartesianPropagatedResults(1,1),CartesianPropagatedResults(1,2),CartesianPropagatedResults(1,3),...
    commandedDirectionCosineMatrix(2,1,1),commandedDirectionCosineMatrix(2,2,1),commandedDirectionCosineMatrix(2,3,1),...
    'AutoScaleFactor',1e4,'LineWidth',1)
quiver3(CartesianPropagatedResults(1,1),CartesianPropagatedResults(1,2),CartesianPropagatedResults(1,3),...
    commandedDirectionCosineMatrix(3,1,1),commandedDirectionCosineMatrix(3,2,1),commandedDirectionCosineMatrix(3,3,1),...
    'AutoScaleFactor',1e4,'LineWidth',1)
plot3(CartesianPropagatedResults(:,1),CartesianPropagatedResults(:,2),CartesianPropagatedResults(:,3),'LineWidth',1.5)
plot3(CartesianEstimatedResults(:,1),CartesianEstimatedResults(:,2),CartesianEstimatedResults(:,3),'LineWidth',1.5)
[x,y,z] = sphere; surf(marsRadius/1e3*x,marsRadius/1e3*y,marsRadius/1e3*z)
hold off
xlabel('x [km]'), ylabel('y [km]'), zlabel('z [km]')
grid on
axis equal tight
set(gca,'FontSize',15)

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

%...Plot error in Cartesian translational motion
F = figure('rend','painters','pos',figSizeLarge);
for i = 1:size(CartesianPropagatedResults,2)
    subplot(2,3,i)
    hold on
    plot(simulationTime,CartesianEstimatedResults(:,i)-CartesianPropagatedResults(:,i),'LineWidth',1.25)
    hold off
    xlabel(timeLabel)
    ylabel(CartesianLabelsDifference{i})
    set(gca,'FontSize',15)
    grid on
end

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
for i = size(rotationalPropagatedResults,2):-1:1
    subplot(2,4,i)
    hold on
    plot(simulationTime,rotationalPropagatedResults(:,i),'LineWidth',1.25)
    plot(onboardTime,rotationalEstimatedResults(:,i),'LineWidth',1.25)
    if i <= 4
        plot(onboardTime,commandedRotationalState(:,i),'LineWidth',1.25)
    end
    hold off
    xlabel(timeLabel)
    ylabel(rotationLabels{i})
    set(gca,'FontSize',15)
    grid on
end
subplotLegend({'Actual','Estimated','Commanded'})

%% Plot Other Results

%...Plot measurements
F = figure('rend','painters','pos',figSizeLarge);

subplot(1,2,1)
hold on
for i = 1:3
    plot(measurementTime,inertialMeasurements(:,i),'LineWidth',1.25)
    plot(onboardTime(1:end-1),expectedMeasurements(:,i),'LineWidth',1.25,'LineStyle','--')
    plot(simulationTime,dependentVariables(:,i+3),'LineWidth',1.25,'LineStyle',':')
end
hold off
xlabel(timeLabel)
ylabel('Translational Acceleration [m s^{-1}]')
set(gca,'FontSize',15)
grid on
legend('x_{IMU}','x_{exp}','x_{real}',...
    'y_{IMU}','y_{exp}','y_{real}',...
    'z_{IMU}','z_{exp}','z_{real}')

subplot(1,2,2)
plot(measurementTime,inertialMeasurements(:,4:6),'LineWidth',1.25)
hold off
xlabel(timeLabel)
ylabel('Rotational Velocity [rad s^{-1}]')
set(gca,'FontSize',15)
grid on
legend('x','y','z')

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
plot(simulationTime,dependentVariables(:,7:9),'LineWidth',1.25)
xlabel(timeLabel)
ylabel('Angle [deg]')
set(gca,'FontSize',15)
grid on
legend('Attack','Side-slip','Bank','Location','Best')

%% Plot Density

%...Find pericenter height
estimatedAltitude = sqrt( sum( CartesianEstimatedResults(:,1:3).^2, 2 ) ) - marsRadius / 1e3;
pericenter = min( estimatedAltitude );

%...Retireve density
aerodynamicAcceleration = sqrt( sum( inertialMeasurements(:,1:3).^2, 2 ) );
atmosphericDensity = 2 * 1000 / 37.5 ./ sum( ( CartesianEstimatedResults(1:end-1,4:6) * 1e3 ).^2, 2 ) / ...
    sqrt( 1.87^2 + 0.013^2 ) .* aerodynamicAcceleration;

%...Data below atmospheric interface
loc = ( estimatedAltitude <= marsAtmosphericInterface );
loc(find(estimatedAltitude==pericenter)+1:end) = false;
atmospherePhaseTime = simulationTime;%simulationTime( loc );
% estimatedAltitude = estimatedAltitude( loc );
% atmosphericDensity = atmosphericDensity( loc( 1:end-1 ) );

%...Plot density
figure;
hold on
plot(atmospherePhaseTime(1:end-1),atmosphericDensity,'LineWidth',1.25)
plot(atmospherePhaseTime,dependentVariables(:,10),'LineWidth',1.25)
hold off
xlabel(timeLabel)
ylabel('Density [kg m^{-3}]')
set(gca,'FontSize',15,'YScale','log')
grid on
legend('Measured','Actual','Location','Best')