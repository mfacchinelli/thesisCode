fclose('all'); clear all; close all force; profile off; clc; format long g; rng default;
addpath functions tests

%% Settings

%...Figure settings
showFigure = false;
showAllFigures = false;
saveFigure = false;
[figSizeLarge,figSizeMedium,figSizeSmall] = saveFigureSettings(saveFigure);

%...Constants
marsRadius = 3389526.666666667;
marsGravitationalParameter = 42828375815756.1;
marsAtmosphericInterface = 175;

%...Plot settings
loadMeasurements = false;
loadFilter = true;
applyInterpolation = true;

%...Labels
timeConversion = 3600 * 24;
timeLabel = 'Time [d]';
CartesianLabels = {'x [km]','y [km]','z [km]','v_x [m s^{-1}]','v_y [m s^{-1}]','v_z [m s^{-1}]'};
CartesianLabelsDifference = {'\Delta x [km]','\Delta y [km]','\Delta z [km]',...
    '\Delta v_x [m s^{-1}]','\Delta v_y [m s^{-1}]','\Delta v_z [m s^{-1}]'};
KeplerianLabels = {'a [km]','e [-]','i [deg]','\omega [deg]','\Omega [deg]','\vartheta [deg]'};
rotationLabels = {'\eta [-]','\epsilon_1 [-]','\epsilon_2 [-]','\epsilon_3 [-]','Norm Offset [-]',...
    '\omega_1 [deg s^{-1}]','\omega_2 [deg s^{-1}]','\omega_3 [deg s^{-1}]'};

%...Main output folder
extension = 'high_ecc/';
mainOutputFolder = ['SimulationOutputTransOnlyIMANRMS/',extension];

%% Loop Over Each Settings

%...Predefine variables
simulations = 0:7;
if strcmp(extension,'high_ecc/')
    ratios = [10,100,1000];
elseif strcmp(extension,'low_ecc/')
    ratios = [10,100];
else
    error('Input folder not recognized.')
end
rmsErrors = cell(length(simulations),length(ratios));
finalErrors = cell(length(simulations),length(ratios));
standardDeviations = cell(length(simulations),length(ratios));

%...Loop
for simulation = simulations
    for ratio = ratios
        %...Output folder
        outputFolder = fullfile(mainOutputFolder,num2str(simulation),num2str(ratio));
        
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
        
        %...Clean up
        clear filename fileID
        
        %% Load C++ Results For Estimation
        
        %...Load translational motion
        filename = ['/Users/Michele/GitHub/tudat/tudatBundle/tudatApplications/Thesis/',outputFolder,'/cartesianEstimated.dat'];
        fileID = fopen(filename,'r');
        CartesianEstimatedResults = textscan(fileID,repmat('%f',[1,7]),'Delimiter',',','CollectOutput',true);
        onboardTime = ( CartesianEstimatedResults{1}(1:end-1,1) - initialTime ) / timeConversion;
        CartesianEstimatedResults = CartesianEstimatedResults{1}(1:end-1,2:end);
        CartesianEstimatedResults(:,1:3) = CartesianEstimatedResults(:,1:3) / 1e3;
        fclose(fileID);
        
        %...Load Keplerian translational motion
        filename = ['/Users/Michele/GitHub/tudat/tudatBundle/tudatApplications/Thesis/',outputFolder,'/keplerianEstimated.dat'];
        fileID = fopen(filename,'r');
        KeplerianEstimatedResults = textscan(fileID,repmat('%f',[1,7]),'Delimiter',',','CollectOutput',true);
        KeplerianEstimatedResults = KeplerianEstimatedResults{1}(1:end-1,2:end);
        KeplerianEstimatedResults(:,1) = KeplerianEstimatedResults(:,1) / 1e3;
        KeplerianEstimatedResults(:,3:end) = rad2deg(KeplerianEstimatedResults(:,3:end));
        fclose(fileID);
        
        %...Clean up
        clear filename fileID
        
        %% Load C++ Results For Navigation Filter
        
        %...Only if filtering is toggled
        if loadFilter
            %...Load translational motion
            filename = ['/Users/Michele/GitHub/tudat/tudatBundle/tudatApplications/Thesis/',...
                outputFolder,'/filterStateEstimates.dat'];
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
        
        %% Interpolate Results to Match Times

        %...Set interpolation time
        interpolatedTime = onboardTime;
        
        %...Interpolate
        if applyInterpolation
            %...Interpolate propagation results
            CartesianPropagatedResults = interp1( simulationTime, CartesianPropagatedResults, interpolatedTime, 'linear' );
            KeplerianPropagatedResults = interp1( simulationTime, KeplerianPropagatedResults, interpolatedTime, 'linear' );
            
            %...Interpolate estimation results
            CartesianEstimatedResults = interp1( onboardTime, CartesianEstimatedResults, interpolatedTime, 'linear' );
            KeplerianEstimatedResults = interp1( onboardTime, KeplerianEstimatedResults, interpolatedTime, 'linear' );
            
            %..Interpolate filter results
            if loadFilter
                filterStateEstimatedResults = interp1( filterTime, filterStateEstimatedResults, interpolatedTime, 'linear', NaN );
                filterCovarianceEstimatedResults = interp1( filterTime, filterCovarianceEstimatedResults, ...
                    interpolatedTime, 'linear', NaN );
            end
        end
        
        %% Plot 3D Orbit
        
        %...Plot trajectory
        if showAllFigures
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
        
        if showAllFigures
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
        end
        
        %% Plot Filter States
        
        %...Only if filtering is toggled
        if showFigure && loadFilter
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
            actual = [-1.21966e-05  7.44836e-05  0.000102096];
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
            
            %...Plot drag coefficient
            actual = 1.9;
            F = figure('rend','painters','pos',figSizeLarge);
            hold on
            plot(interpolatedTime,filterStateEstimatedResults(:,10)-actual,'LineWidth',1.25)
            plot(interpolatedTime(2:end),sqrt(filterCovarianceEstimatedResults(2:end,10)),'LineWidth',1.25,'LineStyle','--')
            plot(interpolatedTime(2:end),-sqrt(filterCovarianceEstimatedResults(2:end,10)),'LineWidth',1.25,'LineStyle','--')
            hold off
            xlabel(timeLabel)
            set(gca,'FontSize',15)
            grid on
        end
        
        %% Store Results
        
        %...RMS Error
        rmsErrors{simulation+1,log10(ratio)} = [rmsPositionError,rmsVelocityError];
        
        %...Final errors
        finalErrors{simulation+1,log10(ratio)} = CartesianEstimatedResults(end,:) - CartesianPropagatedResults(end,:);
        
        %...STD values
        if loadFilter
            standardDeviations{simulation+1,log10(ratio)} = {interpolatedTime,sqrt(filterCovarianceEstimatedResults)};
        end
        
        %% Clean Up
        
        if ~showFigure, close all, end
    end
end

%...Clean up
clc
clearvars -except simulations ratios rmsErrors finalErrors standardDeviations ...
    showFigure saveFigure figSizeLarge figSizeMedium figSizeSmall extension ...
    marsRadius marsGravitationalParameter marsAtmosphericInterface ...
    loadMeasurements loadFilter applyInterpolation ...
    timeConversion timeLabel CartesianLabels CartesianLabelsDifference KeplerianLabels rotationLabels

%% Analyze Results

%...Plot STD
if loadFilter
    F = figure('rend','painters','pos',figSizeLarge);
    for i = 1:3
        for j = 1:length(simulations)
            subplot(8,3,i+3*(j-1))
            hold on
            for k = 1:length(ratios)
                plot(standardDeviations{j,k}{1},standardDeviations{j,k}{2}(:,i),'LineWidth',1.25)
            end
            hold off
            xlabel(timeLabel)
            ylabel(CartesianLabels{i})
            grid on
            set(gca,'FontSize',15)
        end
    end
end

%...LaTeX tables
rates = 2000 ./ ratios;
fprintf(['& ',repmat('{%d} & ',[1,length(simulations)-1]),'{%d} \\\\\n'],1:length(simulations))
fprintf('\\midrule\n')
if strcmp(extension,'high_ecc/')
    fprintf('\\num{%d} & %.3f & %.2f & %.1f & %.1f & %.1f & %.1f & %.2f & %.2f \\\\\n',[rates',cellfun(@(x)x(1),rmsErrors)']')
elseif strcmp(extension,'low_ecc/')
    fprintf('\\num{%d} & %.3f & %.2f & %.3f & %.2f & %.2f & %.2f & %.2f & %.2f \\\\\n',[rates',cellfun(@(x)x(1),rmsErrors)']')
end