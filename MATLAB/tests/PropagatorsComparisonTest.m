fclose('all'); clear all; close all force; profile off; clc; format long g; rng default;
addpath tests data functions

%% Settings

%...Figure settings
plotAllFigures = false;
saveFigure = false;
[figSizeLarge,~,figSizeSmall] = saveFigureSettings(saveFigure);

%...Result repository
TudatApplicationOutput = '/Users/Michele/GitHub/tudat/tudatBundle/tudatApplications/Test/SimulationOutput';
repository = fullfile(TudatApplicationOutput,'Propagators');

%...Select test case
%	0: Aerocapture
%	1: Full aerobraking
% 	2: Interplanetary trajectory
%	3: Circular orbit at LEO (Low Earth Orbit)
% 	4: Molniya orbit
%	5: Low-thrust trajectory
testCase = 1;
switch testCase
    case 0 % aerocapture
        R = 3.396e3;
        simulationDuration = 0.625;
        scaling = 3600;
        timeLabel = 'Time [h]';
        timeStepSeconds = 10;
        constantStepSizes = [1.0, 5.0, 10.0, 25.0, 50.0, 75.0, 100.0, 150.0, 200.0];
        repository = fullfile(repository,'aero');
    case 1 % full aerobraking
        R = 3.396e3;
        simulationDuration = 140.0;
        scaling = 24*3600;
        timeLabel = 'Time [day]';
        timeStepSeconds = 250;
        constantStepSizes = [20.0, 30.0, 40.0, 50.0, 75.0, 100.0, 150.0, 200.0, 250.0, 300.0];
        repository = fullfile(repository,'aero_full');
    case 2 % interplanetary trajectory
        R = 695.508e3;
        simulationDuration = 260.0;
        scaling = 24*3600;
        timeLabel = 'Time [day]';
        timeStepSeconds = 350;
        constantStepSizes = [50.0, 100.0, 200.0, 400.0, 800.0, 1200.0, 1500.0];
        repository = fullfile(repository,'inter');
    case 3 % circular orbit
        R = 6378.1363;
        simulationDuration = 10.0;
        scaling = 24*3600;
        timeLabel = 'Time [day]';
        timeStepSeconds = 50;
        constantStepSizes = [20.0, 30.0, 40.0, 50.0, 75.0, 100.0, 150.0, 200.0];
        repository = fullfile(repository,'circ');
    case 4 % Molniya orbit
        R = 6378.1363;
        simulationDuration = 25.0;
        scaling = 24*3600;
        timeLabel = 'Time [day]';
        timeStepSeconds = 50;
        constantStepSizes = [20.0, 30.0, 40.0, 50.0, 75.0, 100.0, 150.0, 200.0, 250.0, 300.0];
        repository = fullfile(repository,'moln');
    case 5 % low-thurst orbit
        R = 6378.1363;
        simulationDuration = 10.0;
        scaling = 24*3600;
        timeLabel = 'Time [day]';
        timeStepSeconds = 50;
        constantStepSizes = [20.0, 30.0, 40.0, 50.0, 75.0, 100.0, 150.0, 200.0, 250.0, 300.0];
        repository = fullfile(repository,'low_thrust');
end

values = [7,length(constantStepSizes)];
standardFunctionEvals = [8,13;4,0]; % [RK5(6),RK7(8);RK4,-]

%...Labels
keplerLabels = {'Semi-major Axis [km]','Eccentricity [-]','Inclination [deg]',...
    'Right Ascension of Ascending Node [deg]','Argument of Perigee [deg]','True Anomaly [deg]'};
keplerDiffLabels = {'\Delta a [km]','\Delta e [-]','\Delta i [deg]',...
    '\Delta\Omega [deg]','\Delta\omega [deg]','\Delta\vartheta [deg]'};
cartDiffLabels = {'\Delta x [km]','\Delta y [km]','\Delta z [km]','\Delta r [km]',...
    '\Delta v_x [km/s]','\Delta v_y [km/s]','\Delta v_z [km/s]','\Delta v [km/s]'};
usm7Labels = {'C Hodograph [km/s]','R_1 Hodograph [km/s]','R_2 Hodograph [km/s]',...
    'q_1 [-]','q_2 [-]','q_3 [-]','q_4 [-]','Norm Offset [-]'};
usm6Labels = {'C Hodograph [km/s]','R_1 Hodograph [km/s]','R_2 Hodograph [km/s]',...
    '\sigma_1 [-]','\sigma_2 [-]','\sigma_3 [-]'};
usmemLabels = {'C Hodograph [km/s]','R_1 Hodograph [km/s]','R_2 Hodograph [km/s]',...
    'e_1 [-]','e_2 [-]','e_3 [-]'};

%% File Data

propagators = {'cowell','usm7','usm6','usmem','ref'};
propagatorNames = {'Cowell','USM7','USM6','USMEM','Reference'};
integrators = {'var','const'};
integratorNames = {'Variable Step Size','Constant Step Size'};
results = repmat(struct('prop',[],'int',[],'val',[],'kepler',[],'cartesian',[],...
    'timeOut',[],'keplerOut',[],'cartesianOut',[],'evaluations',[]),[length(propagators),1]);

%% Read Propagated Orbits

%...Loop over propagators
for j = 1:length(integrators)
    for k = 1:values(j)
        if plotAllFigures
            F = figure('rend','painters','pos',figSizeLarge);
            hold on
        end
        for i = 1:length(propagators)
            if ~( i == length(propagators) && j == length(integrators) ) && ...
                    ~( i == length(propagators) && k ~= 1 )
%                 if i == length(propagators) % <<<<<<<------ USMEM as reference
%                     repository_old = repository;
%                     [repository,folder] = fileparts(repository);
%                     repository = fullfile(repository,'USMEM',folder);
%                 end
                
                %...Get orbital data
                fileName = fullfile(repository,['orbit_',...
                    propagators{i},'_',integrators{j},'_1_',num2str(k),'.dat']);
                fileID = fopen(fileName);
                kepler = textscan(fileID,repmat('%f ',[1,7]),'CollectOutput',true,'Delimiter',',');
                time = (kepler{1}(:,1)-kepler{1}(1))/scaling; kepler = kepler{1}(:,2:end);
                kepler(:,1) = kepler(:,1)/1e3; kepler(:,3:end) = rad2deg(kepler(:,3:end));
                fclose(fileID);
                
                %...Get trajecotry data
                fileName = fullfile(repository,['trajectory_',...
                    propagators{i},'_',integrators{j},'_1_',num2str(k),'.dat']);
                fileID = fopen(fileName);
                cartesian = textscan(fileID,repmat('%f ',[1,8]),'CollectOutput',true,'Delimiter',',');
                cartesian = cartesian{1}(:,2:7); cartesian = cartesian/1e3;
                fclose(fileID);
                
                %...Get time step data
                fileName = fullfile(repository,['eval_',...
                    propagators{i},'_',integrators{j},'_1_',num2str(k),'.dat']);
                fileID = fopen(fileName);
                l = 1; if i == length(propagators), l = 2; end
                evaluations = textscan(fileID,repmat('%f ',[1,2]),'CollectOutput',true,'Delimiter',',');
                evaluations = evaluations{1}; evaluations(:,1) = (evaluations(:,1)-evaluations(1,1))/scaling;
                evaluations = evaluations(ismembertol(evaluations(:,1),time),2);
                evaluations = vertcat(evaluations(1),diff(evaluations));
                evaluations = round( evaluations(1:end-1) / standardFunctionEvals(j,l) );
                fclose(fileID);
                
                %...Interpolate data
                timeConstantStepSize = 0:timeStepSeconds/scaling:simulationDuration;
                results(i,j,k).prop = propagatorNames{i};
                results(i,j,k).int = integratorNames{j};
                results(i,j,k).val = k;
                results(i,j,k).timeOut = time;
                results(i,j,k).keplerOut = kepler;
                results(i,j,k).cartesianOut = cartesian;
                results(i,j,k).kepler = interp1(time,kepler,timeConstantStepSize,'spline',NaN);
                results(i,j,k).cartesian = interp1(time,cartesian,timeConstantStepSize,'spline',NaN);
                results(i,j,k).evaluations = evaluations;
                
                %...Plot time steps
                if i ~= length(propagators) && j == 1 && plotAllFigures
                    subplot(2,2,i)
                    yyaxis right
                    scatter(time(1:end-1),evaluations,25,'filled')
                    ylabel('Evaluations [-]')
                    yyaxis left
                    plot(time(1:end-1),diff(time)*scaling,'LineWidth',1.1)
                    ylabel('Time Step [s]')
                    xlabel(timeLabel)
                    grid on
                    set(gca,'FontSize',15,'YScale','log')
                    title(propagatorNames{i})
                end
                
%                 if i == length(propagators), repository = repository_old; end % <<<<<<<------ USMEM as reference
            end
        end
    end
    if plotAllFigures
        hold off
    end
end
clear fileName fileID k kepler cartesian evaluations time

if plotAllFigures
    F = figure('rend','painters','pos',figSizeLarge);
    yyaxis right
    scatter(results(end,1).timeOut(1:end-1),results(end,1).evaluations,35,'filled')
    ylabel('Evaluations [-]')
    yyaxis left
    plot(results(end,1,1).timeOut(1:end-1),diff(results(end,1,1).timeOut)*scaling,'LineWidth',1.1)
    ylabel('Time Step [s]')
    xlabel(timeLabel)
    grid on
    set(gca,'FontSize',15,'YScale','log')
    if ~saveFigure, title('Reference'), end
end

%% Computation Times

%...Settings
numberOfConstantSteps = [9,10,7,8,10,10];
numberOfVariableSteps = 7;
fullConstantStepSizes = {[1.0, 5.0, 10.0, 25.0, 50.0, 75.0, 100.0, 150.0, 200.0],...
    [20.0, 30.0, 40.0, 50.0, 75.0, 100.0, 150.0, 200.0, 250.0, 300.0],...
    [50.0, 100.0, 200.0, 400.0, 800.0, 1200.0, 1500.0],...
    [20.0, 30.0, 40.0, 50.0, 75.0, 100.0, 150.0, 200.0],...
    [20.0, 30.0, 40.0, 50.0, 75.0, 100.0, 150.0, 200.0, 250.0, 300.0],...
    [20.0, 30.0, 40.0, 50.0, 75.0, 100.0, 150.0, 200.0, 250.0, 300.0]};

%...Open file
filename = fullfile(fullfile(TudatApplicationOutput,'Propagators'),'computationTimes.dat');
fileID = fopen(filename,'r');

%...Read lines seprately
i = 0;
computationTimes = cell(length(propagators)-1,length(integrators),6);
while ~feof(fileID)
    i = i+1;
    currentLine = cell2mat(cellfun(@str2num,split(fgetl(fileID)),'UniformOutput',false))';
    init = 1; fin = numberOfVariableSteps;
    for j = 1:length(propagators)-1
        if j ~= 1, init = fin + 1; fin = fin + numberOfVariableSteps; end
        computationTimes{j,1,i} = currentLine(init:fin);
        init = fin + 1; fin = fin + numberOfConstantSteps(i);
        computationTimes{j,2,i} = currentLine(init:fin);
    end
end

%...Close file
fclose(fileID);

% %...Plot computation times for each propagator
% styles = {'-o','-d','-s','-v','-p','-h','-*','-x','-^','-o','-d'};
% for i = 1:6*length(integrators)
%     F = figure('rend','painters','pos',figSizeSmall);
%     hold on
%     if mod(i,2) == 0
%         j = 2;
%         for k = 1:length(propagators)-1
%             plot(fullConstantStepSizes{ceil(i/2)},vertcat(computationTimes{k,j,ceil(i/2)}),...
%                 styles{k},'LineWidth',1.5,'MarkerSize',10)
%         end
%         xlabel('Constant Step Size [-]')
%     else
%         j = 1;
%         for k = 1:length(propagators)-1
%             plot(10.^(-13:-7),vertcat(computationTimes{k,j,ceil(i/2)}),styles{k},'LineWidth',1.5,'MarkerSize',10)
%         end
%         xlabel('Integration Tolerance [-]')
%         set(gca,'XScale','log')
%     end
%     hold off
%     ylabel('Computation Time [s]')
%     legend(propagatorNames(1:end-1),'Location','Best')
%     grid on
%     set(gca,'FontSize',15,'YScale','log')
%     if saveFigure, saveas(F,['../../Report/figures/comp_time_',num2str(ceil(i/2)),'_',num2str(j)],'epsc'), end
% end

%% Plot Reference Trajectory

ref = length(propagators);
reference = results(ref,1,1);
time = timeConstantStepSize;

%...Plot Kepler elements
F = figure('rend','painters','pos',figSizeLarge);
for i = 1:6
    subplot(2,3,i)
    hold on
    plot(reference.timeOut,reference.keplerOut(:,i),'LineWidth',1.1)
    if i == 1
        plot([reference.timeOut(1),reference.timeOut(end)],R*ones(1,2),'LineStyle','--')
    end
    hold off
    if testCase == 0 && i == 1
        ylim([-2e4,1.5e4])
    end
    xlabel(timeLabel)
    ylabel(keplerLabels{i})
    grid on
    set(gca,'FontSize',15)
end

%...Plot Cartesian position
[x,y,z] = sphere;
F = figure('rend','painters','pos',figSizeLarge);
hold on
plot3(reference.cartesianOut(:,1),reference.cartesianOut(:,2),...
    reference.cartesianOut(:,3),'LineWidth',1.5)
surf(R*x,R*y,R*z)
hold off
xlabel('x [km]'), ylabel('y [km]'), zlabel('z [km]')
grid on
axis equal tight
set(gca,'FontSize',15)

%% Plot Function Evaluations

if plotAllFigures
    %...Plot true anomaly
    for k = 1:values
        F = figure('rend','painters','pos',figSizeLarge);
        for i = 1:length(propagators)-1
            subplot(2,2,i)
            scatter(results(i,1,k).timeOut(1:end-1),results(i,1,k).keplerOut(1:end-1,6),...
                25,results(i,1,k).evaluations,'filled')
            c = colorbar; c.Label.String = 'Evaluations [-]'; c.Location = 'southoutside';
            colormap jet
            timeArray = 0:15:simulationDuration; timeArray(end) = simulationDuration;
            xticks(timeArray), xlim([0,simulationDuration])
            yticks(0:90:360), ylim([0,360])
            xlabel(timeLabel), ylabel(keplerLabels{end})
            grid on
            set(gca,'FontSize',15)
            title(propagatorNames{i})
        end
    end
    
    %...Plot area of interest
    F = figure('rend','painters','pos',figSizeLarge);
    offsetTrueAnomaly = 30;
    loc = ( reference.keplerOut(:,6) > ( 360 - offsetTrueAnomaly ) ) | ...
        ( reference.keplerOut(:,6) < offsetTrueAnomaly ); loc(end) = false;
    scatter3(reference.cartesianOut(loc,1),reference.cartesianOut(loc,2),...
        reference.cartesianOut(loc,3),50,reference.evaluations(loc(1:end-1)),'filled')
    c = colorbar; c.Label.String = 'Evaluations [-]';
    colormap jet
    xlabel('x [km]'), ylabel('y [km]'), zlabel('z [km]')
    grid on
    axis equal tight
    view([-30,25])
    set(gca,'FontSize',15)
    title('Reference')
end

%% Compare Results of Cowell Propagator and USM

%...Difference in Cartesian elements
rmsError = zeros(length(propagators),8,length(integrators),max(values));
for k = 1:length(integrators)
    for l = 1:values(k)
        if plotAllFigures, F = figure('rend','painters','pos',figSizeLarge); end
        for i = 1:8
            if plotAllFigures
                subplot(2,4,i)
                hold on
            end
            for j = 1:length(propagators)-1
                if i < 4
                    rmsError(j,i,k,l) = rms(results(j,k,l).cartesian(:,i) - reference.cartesian(:,i));
                    if plotAllFigures
                        plot(time,abs(results(j,k,l).cartesian(:,i)-reference.cartesian(:,i)),'LineWidth',1.1)
                    end
                elseif i > 4 && i < 8
                    h = i - 1;
                    rmsError(j,i,k,l) = rms(results(j,k,l).cartesian(:,h) - reference.cartesian(:,h));
                    if plotAllFigures
                        plot(time,abs(results(j,k,l).cartesian(:,k)-reference.cartesian(:,k)),'LineWidth',1.1)
                    end
                elseif i == 4
                    rmsError(j,i,k,l) = rms(sqrt(sum(results(j,k,l).cartesian(:,1:3).^2,2)) - ...
                        sqrt(sum(reference.cartesian(:,1:3).^2,2)));
                    if plotAllFigures
                        plot(time,abs(sqrt(sum(results(j,k,l).cartesian(:,1:3).^2,2)) - ...
                        sqrt(sum(reference.cartesian(:,1:3).^2,2))),'LineWidth',1.1)
                    end
                else
                    rmsError(j,i,k,l) = rms(sqrt(sum(results(j,k,l).cartesian(:,4:6).^2,2)) - ...
                        sqrt(sum(reference.cartesian(:,4:6).^2,2)));
                    if plotAllFigures
                        plot(time,abs(sqrt(sum(results(j,k,l).cartesian(:,4:6).^2,2)) - ...
                        sqrt(sum(reference.cartesian(:,4:6).^2,2))),'LineWidth',1.1)
                    end
                end
            end
            if plotAllFigures
                hold off
                xlabel(timeLabel)
                ylabel(cartDiffLabels{i})
                grid on
                set(gca,'FontSize',15,'YScale','log')
            end
        end
        if plotAllFigures
            subplotLegend(propagatorNames(1:end-1))
            if ~saveFigure, subplotTitle(integratorNames{k}), end
        end
    end
end

%% Root-mean-squared Error

%...Table of RMS error
if plotAllFigures
    format short g
    for l = 1:values
        rmsErrorTableVar = array2table([rmsError(1:end-1,:,1,l),arrayfun(@(i)sum(results(i,1,l).evaluations),...
            1:size(results,1)-1)'],'VariableNames',{'x','y','z','r','v_x','v_y','v_z','v','eval'},...
            'RowNames',propagatorNames(1:end-1))
        rmsErrorTableCont = array2table([rmsError(1:end-1,:,2,l),arrayfun(@(i)sum(results(i,2,l).evaluations),...
            1:size(results,1)-1)'],'VariableNames',{'x','y','z','r','v_x','v_y','v_z','v','eval'},...
            'RowNames',propagatorNames(1:end-1))
    end
    format long g
end

styles = {'-o','-d','-s','-v','-p','-h','-*','-x','-^','-o','-d'};

%...Plot RMS error for variable step size
F = figure('rend','painters','pos',figSizeSmall);
hold on
for i = 1:length(propagators)-1
    plot(arrayfun(@(j)sum(results(i,1,j).evaluations),1:values(1)),...
        squeeze(rmsError(i,4,1,1:values(1)))*1e3,styles{i},'LineWidth',1.5,'MarkerSize',10)
end
hold off
xlabel('Function Evaluations [-]')
ylabel('RMS Position Error [m]')
legend(propagatorNames(1:end-1),'Location','Best')
if testCase == 2
    set(gca,'FontSize',15,'XScale','log','YScale','log')
else
    set(gca,'FontSize',15,'YScale','log')
end
grid on
if saveFigure, saveas(F,['../../Report/figures/rms_var_',num2str(testCase)],'epsc'), end

F = figure('rend','painters','pos',figSizeSmall);
hold on
for i = 1:length(propagators)-1
    plot(computationTimes{i,1,testCase+1},squeeze(rmsError(i,4,1,1:values(1)))*1e3,styles{i},'LineWidth',1.5,'MarkerSize',10)
end
hold off
xlabel('Computation Time [s]')
ylabel('RMS Position Error [m]')
legend(propagatorNames(1:end-1),'Location','Best')
if testCase == 2
    set(gca,'FontSize',15,'XScale','log','YScale','log')
else
    set(gca,'FontSize',15,'YScale','log')
end
grid on
if saveFigure, saveas(F,['../../Report/figures/rms_var_time_',num2str(testCase)],'epsc'), end
    
%...Plot RMS error for constant step size
F = figure('rend','painters','pos',figSizeSmall);
hold on
for i = 1:length(propagators)-1
    plot(constantStepSizes,squeeze(rmsError(i,4,2,1:values(2)))*1e3,styles{i},'LineWidth',1.5,'MarkerSize',10)
end
hold off
xlabel('Constant Time Step [s]')
ylabel('RMS Position Error [m]')
legend(propagatorNames(1:end-1),'Location','Best')
set(gca,'FontSize',15,'YScale','log')
grid on
if saveFigure, saveas(F,['../../Report/figures/rms_const_',num2str(testCase)],'epsc'), end

F = figure('rend','painters','pos',figSizeSmall);
hold on
for i = 1:length(propagators)-1
    plot(computationTimes{i,2,testCase+1},squeeze(rmsError(i,4,2,1:values(2)))*1e3,styles{i},'LineWidth',1.5,'MarkerSize',10)
end
hold off
xlabel('Computation Time [s]')
ylabel('RMS Position Error [m]')
legend(propagatorNames(1:end-1),'Location','Best')
set(gca,'FontSize',15,'YScale','log')
grid on
if saveFigure, saveas(F,['../../Report/figures/rms_const_time_',num2str(testCase)],'epsc'), end

%% Close All Figures

if saveFigure, close all, end