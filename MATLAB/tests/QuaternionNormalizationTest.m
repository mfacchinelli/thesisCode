fclose('all'); clear all; close all force; profile off; clc; format long g; rng default;
addpath tests data functions

%% Settings

%...Constants
radius = 6378.1363;
numberOfNormalizationMethods = 6;

%...Result repository
TudatApplicationOutput = '/Users/Michele/GitHub/tudat/tudatBundle/tudatApplications/Test/SimulationOutput/Quaternions';

%...Figure settings
showFigure = true;
saveFigure = false;
[figSizeLarge,figSizeMedium,figSizeSmall] = saveFigureSettings(saveFigure);

%...Labels
timeConversion = 3600 * 24;
timeLabel = 'Time [d]';
CartesianLabels = {'x [km]','y [km]','z [km]','v_x [km s^{-1}]','v_y [km s^{-1}]','v_z [km s^{-1}]'};
CartDiffLabels = {'\Delta x [km]','\Delta y [km]','\Delta z [km]','\Delta r [km]',...
    '\Delta v_x [km/s]','\Delta v_y [km/s]','\Delta v_z [km/s]','\Delta v [km/s]'};
USM7Labels = {'C Hodograph [km/s]','R_1 Hodograph [km/s]','R_2 Hodograph [km/s]',...
    'q_1 [-]','q_2 [-]','q_3 [-]','q_4 [-]','Quaternion Norm Offset [-]'};
legendLabels = {'1) None','2) Phillips, Variable','3) Phillips, Constant',...
    '4) Fukushima, Quaternion','5) Fukushima, Derivative','6) Fukushima, Both'};

%...Computation times
computationTimes = [ 2.70521, 3.11181, 2.89613, 2.81819, 2.8117, 2.80242 ]; % average of 100 runs

%% Load C++ Results

%...Allocate cell sizes
simulationTime = cell(1,numberOfNormalizationMethods);
CartesianTranslationalMotion = cell(1,numberOfNormalizationMethods);
USM7TranslationalMotion = cell(1,numberOfNormalizationMethods);
functionEvaluations = cell(1,numberOfNormalizationMethods);

%...Loop over quaternion normalization methods
for i = 1:numberOfNormalizationMethods
    %...Load Cartesian results
    filename = fullfile( TudatApplicationOutput, ['cartesianResult_',num2str(i-1),'.dat'] );
    fileID = fopen(filename,'r');
    fileOutput = textscan(fileID,repmat('%f',[1,7]),'Delimiter',',','CollectOutput',true);
    simulationTime{i} = ( fileOutput{1}(:,1) - fileOutput{1}(1) ) / timeConversion;
    CartesianTranslationalMotion{i} = fileOutput{1}(:,2:end) / 1e3; % convert to position and velocity
    fclose(fileID);
    
    %...Load USM7 results
    filename = fullfile( TudatApplicationOutput, ['usm7Result_',num2str(i-1),'.dat'] );
    fileID = fopen(filename,'r');
    fileOutput = textscan(fileID,repmat('%f',[1,8]),'Delimiter',',','CollectOutput',true);
    USM7TranslationalMotion{i} = fileOutput{1}(:,2:end);
    USM7TranslationalMotion{i}(:,1:3) = USM7TranslationalMotion{i}(:,1:3) / 1e3; % convert to velocity
    USM7TranslationalMotion{i}(:,end+1) = abs( 1.0 - ...
        sqrt( sum( USM7TranslationalMotion{i}(:,4:end) .^ 2, 2 ) ) ); % add quaternion norm offset
    fclose(fileID);
    
    %...Load function evaluation statistics
    filename = fullfile( TudatApplicationOutput, ['functionEvaluations_',num2str(i-1),'.dat'] );
    fileID = fopen(filename,'r');
    fileOutput = textscan(fileID,repmat('%f',[1,2]),'Delimiter',',','CollectOutput',true);
    functionEvaluations{i} = fileOutput{1}(:,2);
    fclose(fileID);
end

%...Load reference Cartesian results
filename = fullfile( TudatApplicationOutput, ['cartesianResult_',num2str(numberOfNormalizationMethods),'.dat'] );
fileID = fopen(filename,'r');
fileOutput = textscan(fileID,repmat('%f',[1,7]),'Delimiter',',','CollectOutput',true);
referenceTime = ( fileOutput{1}(:,1) - fileOutput{1}(1) ) / timeConversion;
CartesianTranslationalMotionReference = fileOutput{1}(:,2:end) / 1e3; % convert to position and velocity
fclose(fileID);

%% Interpolate Results

%...Loop over quaternion normalization methods
CartesianTranslationalMotionInterpolated = cell(1,numberOfNormalizationMethods);
for i = 1:numberOfNormalizationMethods
    CartesianTranslationalMotionInterpolated{i} = interp1( simulationTime{i}, CartesianTranslationalMotion{i}, ...
        referenceTime, 'spline' );
end

%% Plot Results

%...Plot Cartesian position
F = figure('rend','painters','pos',figSizeLarge);
[x,y,z] = sphere;
hold on
plot3(CartesianTranslationalMotionReference(:,1),CartesianTranslationalMotionReference(:,2),...
    CartesianTranslationalMotionReference(:,3),'LineWidth',1.5)
surf(radius*x,radius*y,radius*z)
hold off
xlabel('x [km]'), ylabel('y [km]'), zlabel('z [km]')
grid on
axis equal tight
set(gca,'FontSize',15)
clear x y z

%...Plot Cartesian elements difference
F = figure('rend','painters','pos',figSizeLarge);
for i = 1:8
    subplot(2,4,i)
    hold on
    for j = 1:numberOfNormalizationMethods
        if i <= 3
            plot(referenceTime,abs( CartesianTranslationalMotionInterpolated{j}(:,i) - ...
                CartesianTranslationalMotionReference(:,i) ),'LineWidth',1.1)
        elseif i == 4
            plot(referenceTime,abs( sqrt( sum( CartesianTranslationalMotionInterpolated{j}(:,1:3) .^ 2, 2) ) - ...
                sqrt( sum( CartesianTranslationalMotionReference(:,1:3) .^ 2, 2 ) ) ),'LineWidth',1.1)
        elseif i == 8
            plot(referenceTime,abs( sqrt( sum( CartesianTranslationalMotionInterpolated{j}(:,4:6) .^2, 2 ) ) - ...
                sqrt( sum( CartesianTranslationalMotionReference(:,4:6) .^ 2, 2) ) ),'LineWidth',1.1)
        else
            plot(referenceTime,abs( CartesianTranslationalMotionInterpolated{j}(:,i-1) - ...
                CartesianTranslationalMotionReference(:,i-1) ),'LineWidth',1.1)
        end
    end
    hold off
    xlabel(timeLabel)
    ylabel(CartDiffLabels{i})
    set(gca,'FontSize',15,'YScale','log')
    grid on
    xlim([referenceTime(end)-1,referenceTime(end)])
end
subplotLegend(legendLabels)

%% Compute RMS Error w.r.t. Reference

%...Difference in Cartesian elements
rmsError = zeros(8,numberOfNormalizationMethods);
for j = 1:numberOfNormalizationMethods
    for i = 1:8
        if i <= 3
            rmsError(i,j) = rms( CartesianTranslationalMotionInterpolated{j}(:,i) - ...
                CartesianTranslationalMotionReference(:,i) );
        elseif i == 4
            rmsError(i,j) = rms( sqrt( sum( CartesianTranslationalMotionInterpolated{j}(:,1:3) .^ 2, 2) ) - ...
                sqrt( sum( CartesianTranslationalMotionReference(:,1:3) .^ 2, 2 ) ) );
        elseif i == 8
            rmsError(i,j) = rms( sqrt( sum( CartesianTranslationalMotionInterpolated{j}(:,4:6) .^2, 2 ) ) - ...
                sqrt( sum( CartesianTranslationalMotionReference(:,4:6) .^ 2, 2) ) );
        else
            rmsError(i,j) = rms(CartesianTranslationalMotionInterpolated{j}(:,i-1) - ...
                CartesianTranslationalMotionReference(:,i-1));
        end
    end
end

%% Plot Results

%...Plot quaternion norm offset for each type
styles = repmat({'-',':','--','-.'},[1,3]);
F = figure('rend','painters','pos',figSizeSmall);
hold on
for i = 1:numberOfNormalizationMethods
    plot(simulationTime{i},USM7TranslationalMotion{i}(:,end),'LineWidth',1.25,'LineStyle',styles{i})
end
hold off
xlabel(timeLabel)
ylabel(USM7Labels{end})
legend(legendLabels{:},'Location',[0.6,0.25,0.2,0.2])
set(gca,'FontSize',15,'XScale','log','YScale','log')
grid on
if saveFigure, saveas(F,'../../Report/figures/quat_offset','epsc'), end

%...Plot RMS error
styles = {'o','d','s','v','h','p','^'};
F = figure('rend','painters','pos',figSizeSmall);
hold on
for i = 1:numberOfNormalizationMethods
    offsetInFunctionEvaluations = ( functionEvaluations{i}(end) - functionEvaluations{1}(end) ) / ...
        functionEvaluations{1}(end) * 100;
    scatter(offsetInFunctionEvaluations,rmsError(4,i)*1e3,250,'filled',styles{i})
end
for i = 1:numberOfNormalizationMethods
    offsetInFunctionEvaluations = ( functionEvaluations{i}(end) - functionEvaluations{1}(end) ) / ...
        functionEvaluations{1}(end) * 100;
    offsetInComputationTime = ( computationTimes(i) - computationTimes(1) ) / ...
        computationTimes(1) * 100;
    text(offsetInFunctionEvaluations,rmsError(4,i)*1e3,[num2str(offsetInComputationTime','%.1f'),' %  '],...
        'FontSize',15,'HorizontalAlignment','right','VerticalAlignment','bottom');
end
hold off
xlabel('Offset in Function Evaluations [%]')
ylabel('RMS Position Error [m]')
[L,icons] = legend(legendLabels{:},'Location','NE');
for i = 1:numberOfNormalizationMethods
    icons(i).FontSize = 12.5;
    icons(i+numberOfNormalizationMethods).Children.MarkerSize = 10;
end
set(gca,'FontSize',15)
grid on
if saveFigure, saveas(F,'../../Report/figures/quat_error','epsc'), end