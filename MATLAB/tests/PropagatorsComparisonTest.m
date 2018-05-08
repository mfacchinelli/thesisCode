fclose('all'); clear all; close all force; profile off; clc; format long g; rng default;
addpath tests data functions

%% Settings

%...Figure settings
saveFigure = false;
figSize = saveFigureSettings(saveFigure);

%...Result repository
TudatApplicationOutput = '/Users/Michele/GitHub/tudat/tudatBundle/tudatApplications/Test/SimulationOutput/';
repository = [TudatApplicationOutput,'Propagators'];

%...Select test case
testCase = 1;
switch testCase
    case 0 % single aerobraking sweep
        simulationDuration = 24;
        scaling = 3600;
        timeLabel = 'Time [h]';
        timeStepSeconds = 10;
        repository = fullfile(repository,'aero');
    case 1 % full aerobraking
        simulationDuration = 150;
        scaling = 24*3600;
        timeLabel = 'Time [day]';
        timeStepSeconds = 100;
        repository = fullfile(repository,'aero_full');
    case 2 % circular orbit
        simulationDuration = 27.5;
        scaling = 24*3600;
        timeLabel = 'Time [day]';
        timeStepSeconds = 100;
        repository = fullfile(repository,'circ');
    case 3 % interplanetary trajectory
        simulationDuration = 30;
        scaling = 24*3600;
        timeLabel = 'Time [day]';
        timeStepSeconds = 100;
        repository = fullfile(repository,'inter');
end

standardFunctionEvals = [8,13;4,0]; % RK5(6),RK7(8),RK4,NaN

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

propagators = {'cowell','encke','gauss','usm7','usm6','usmem','ref'};
propagatorNames = {'Cowell','Encke','Gauss','USM7','USM6','USMEM','Reference'};
integrators = {'var','const'};
integratorNames = {'Variable Step Size','Constant Step Size'};
results = repmat(struct('prop',[],'int',[],'kepler',[],'cartesian',[],...
    'timeOut',[],'keplerOut',[],'cartesianOut',[],'evaluations',[]),[length(propagators),1]);

%% Read Propagated Orbits

%...Mars data
R = 3.390e3;

%...Loop over propagators
F = figure('rend','painters','pos',figSize);
hold on
for i = 1:length(propagators)
    for j = 1:length(integrators)
        if ~( i == length(propagators) && j == 2 )
            %...Get orbital data
            fileName = fullfile(repository,['orbit_',propagators{i},'_',integrators{j},'.dat']);
            fileID = fopen(fileName);
            kepler = textscan(fileID,repmat('%f ',[1,7]),'CollectOutput',true,'Delimiter',',');
            time = (kepler{1}(:,1)-kepler{1}(1))/scaling; kepler = kepler{1}(:,2:end);
            kepler(:,1) = kepler(:,1)/1e3; kepler(:,3:end) = rad2deg(kepler(:,3:end));
            fclose(fileID);
            
            %...Get trajecotry data
            fileName = fullfile(repository,['trajectory_',propagators{i},'_',integrators{j},'.dat']);
            fileID = fopen(fileName);
            cartesian = textscan(fileID,repmat('%f ',[1,8]),'CollectOutput',true,'Delimiter',',');
            cartesian = cartesian{1}(:,2:7); cartesian = cartesian/1e3;
            fclose(fileID);
            
            %...Get time step data
            fileName = fullfile(repository,['eval_',propagators{i},'_',integrators{j},'.dat']);
            fileID = fopen(fileName);
            k = 1; if i == length(propagators), k = 2; end
            evaluations = textscan(fileID,repmat('%f ',[1,2]),'CollectOutput',true,'Delimiter',',');
            evaluations = evaluations{1}; evaluations(:,1) = (evaluations(:,1)-evaluations(1,1))/scaling;
            evaluations = evaluations(ismembertol(evaluations(:,1),time),2);
            evaluations = vertcat(evaluations(1),diff(evaluations)); 
            evaluations = round( evaluations(1:end-1) / standardFunctionEvals(j,k) );
            fclose(fileID);
            
            %...Interpolate data
            timeConstantStepSize = 0:timeStepSeconds/scaling:simulationDuration;
            results(i,j).prop = propagatorNames{i};
            results(i,j).int = integratorNames{j};
            results(i,j).timeOut = time;
            results(i,j).keplerOut = kepler;
            results(i,j).cartesianOut = cartesian;
            results(i,j).kepler = interp1(time,kepler,timeConstantStepSize,'spline',NaN);
            results(i,j).cartesian = interp1(time,cartesian,timeConstantStepSize,'spline',NaN);
            results(i,j).evaluations = evaluations;
            
            %...Plot time steps
            if i ~= length(propagators) && j == 1
                subplot(2,3,i)
                yyaxis right
                scatter(time(1:end-1),evaluations,'LineWidth',1.1)
                ylabel('Evaluations [-]')
                yyaxis left
                plot(time(1:end-1),diff(time)*scaling,'LineWidth',1.1)
                ylabel('Time Step [s]')
                xlabel(timeLabel)
                grid on
                set(gca,'FontSize',15,'YScale','log')
                title(propagatorNames{i})
            end
        end
    end
end
hold off
if saveFigure, saveas(F,['../../Report/figures/prop_comp_dt_',testCase],'epsc'), end
clear fileName fileID k kepler cartesian evaluations time

F = figure('rend','painters','pos',figSize);
yyaxis right
% scatter(results(end,1).timeOut(1:end-1),results(end,1).evaluations,'LineWidth',1.1)
ylabel('Evaluations [-]')
yyaxis left
plot(results(end,1).timeOut(1:end-1),diff(results(end,1).timeOut)*scaling,'LineWidth',1.1)
ylabel('Time Step [s]')
xlabel(timeLabel)
grid on
set(gca,'FontSize',15,'YScale','log')
if saveFigure, saveas(F,['../../Report/figures/prop_ref_dt_',testCase],'epsc'), 
else, title('Reference'), end

%% Plot Reference Trajectory

ref = length(propagators);
reference = results(ref,1);
time = timeConstantStepSize;

%...Plot Kepler elements
F = figure('rend','painters','pos',figSize);
for i = 1:6
    subplot(2,3,i)
    hold on
    plot(reference.timeOut,reference.keplerOut(:,i),'LineWidth',1.1)
    if i == 1
        plot([reference.timeOut(1),reference.timeOut(end)],R*ones(1,2),'LineStyle','--')
    end
    hold off
    xlabel(timeLabel)
    ylabel(keplerLabels{i})
    grid on
    set(gca,'FontSize',15)
end
if saveFigure, saveas(F,['../../Report/figures/prop_kepl_ref_',testCase],'epsc'), end

%...Plot Cartesian position
[x,y,z] = sphere;
F = figure('rend','painters','pos',figSize);
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

%...Plot true anomaly
F = figure('rend','painters','pos',figSize);
for i = 1:length(propagators)-1
    subplot(2,3,i)
    scatter(results(i,1).timeOut(1:end-1),results(i,1).keplerOut(1:end-1,6),...
        25,results(i,1).evaluations,'filled')
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
if saveFigure, saveas(F,['../../Report/figures/prop_peri_',testCase],'epsc'), end

%...Plot area of interest
F = figure('rend','painters','pos',figSize);
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

%% Compare Results of Existing Propagators and USM

% %...Difference in Keplerian elements
% for k = 1:length(integrators)
%     F = figure('rend','painters','pos',figSize);
%     for i = 1:6
%         subplot(2,3,i)
%         hold on
%         for j = 1:length(propagators)-1
%             plot(time,abs(results(j,k).kepler(:,i)-reference.kepler(:,i)),'LineWidth',1.1)
%         end
%         hold off
%         xlabel(timeLabel)
%         ylabel(keplerDiffLabels{i})
%         grid on
%         set(gca,'FontSize',15,'YScale','log')
%     end
%     subplotLegend({'Cowell','Encke','Gauss','USM7','USM6','USMEM'})
%     subplotTitle(integratorNames{k})
% end

%...Difference in Cartesian elements
rmsError = zeros(length(propagators),8,length(integrators));
for k = 1:length(integrators)
    F = figure('rend','painters','pos',figSize);
    for i = 1:8
        subplot(2,4,i)
        hold on
        for j = 1:length(propagators)-1
            if i < 4
                rmsError(j,i,k) = rms(results(j,k).cartesian(:,i) - reference.cartesian(:,i));
                plot(time,abs(results(j,k).cartesian(:,i)-reference.cartesian(:,i)),'LineWidth',1.1)
            elseif i > 4 && i < 8
                h = i - 1;
                rmsError(j,i,k) = rms(results(j,k).cartesian(:,h) - reference.cartesian(:,h));
                plot(time,abs(results(j,k).cartesian(:,k)-reference.cartesian(:,k)),'LineWidth',1.1)
            elseif i == 4
                rmsError(j,i,k) = rms(sqrt(sum(results(j,k).cartesian(:,1:3).^2,2)) - ...
                    sqrt(sum(reference.cartesian(:,1:3).^2,2)));
                plot(time,abs(sqrt(sum(results(j,k).cartesian(:,1:3).^2,2)) - ...
                    sqrt(sum(reference.cartesian(:,1:3).^2,2))),'LineWidth',1.1)
            else
                rmsError(j,i,k) = rms(sqrt(sum(results(j,k).cartesian(:,4:6).^2,2)) - ...
                    sqrt(sum(reference.cartesian(:,4:6).^2,2)));
                plot(time,abs(sqrt(sum(results(j,k).cartesian(:,4:6).^2,2)) - ...
                    sqrt(sum(reference.cartesian(:,4:6).^2,2))),'LineWidth',1.1)
            end
        end
        hold off
        xlabel(timeLabel)
        ylabel(cartDiffLabels{i})
        grid on
        set(gca,'FontSize',15,'YScale','log')
    end
    subplotLegend({'Cowell','Encke','Gauss','USM7','USM6','USMEM'})
    if saveFigure, saveas(F,['../../Report/figures/prop_diff_cart_',testCase],'epsc')
    else, subplotTitle(integratorNames{k}), end
end

%...Table of RMS
format short g
rmsErrorTableVar = array2table([rmsError(1:end-1,:,1),arrayfun(@(i)sum(results(i,1).evaluations),1:size(results,1)-1)'],...
    'VariableNames',{'x','y','z','r','v_x','v_y','v_z','v','eval'},...
    'RowNames',{'Cowell','Encke','Gauss','USM7','USM6','USMEM'})
rmsErrorTableCont = array2table([rmsError(1:end-1,:,2),arrayfun(@(i)sum(results(i,2).evaluations),1:size(results,1)-1)'],...
    'VariableNames',{'x','y','z','r','v_x','v_y','v_z','v','eval'},...
    'RowNames',{'Cowell','Encke','Gauss','USM7','USM6','USMEM'})
format long g

%% Plot USM Elements

%...Loop over propagators
for i = 4:6
    switch propagators{i}
        case 'usm7'
            limit = 8; spyValue = 4;
            correctState = @(x) horzcat(x,abs(1-sqrt(sum(x(:,4:7).^2,2))));
            labels = usm7Labels;
        case 'usm6'
            limit = 8; spyValue = 3;
            correctState = @(x)x;
            labels = usm6Labels;
        case 'usmem'
            limit = 7; spyValue = 3;
            correctState = @(x)x;
            labels = usmemLabels;
        otherwise
            error('USM not recognized')
    end
    
    %...Get USM data
    fileName = fullfile(repository,['usm_',propagators{i},'_',integrators{1},'.dat']);
    fileID = fopen(fileName);
    usm = textscan(fileID,repmat('%f ',[1,limit]),'CollectOutput',true,'Delimiter',',');
    time = (usm{1}(:,1)-usm{1}(1))/scaling; usm = usm{1}(:,2:end); usm(:,1:3) = usm(:,1:3)/1e3;
    usm = correctState(usm);
    fclose(fileID);
    
    %...Plot USM elements
    F = figure('rend','painters','pos',figSize);
    for j = 1:length(labels)
        subplot(2,spyValue,j)
        hold on
        plot(time,usm(:,j),'LineWidth',1.1)
        if strcmpi(propagators{i},'usm6')
            if j >= 4
                plot(time,usm(:,7),'LineStyle','--')
            end
        end
        hold off
        xlabel(timeLabel)
        ylabel(labels{j})
        grid on
        set(gca,'FontSize',15)
    end
end