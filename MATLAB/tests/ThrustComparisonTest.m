fclose('all'); clear all; close all force; profile off; clc; format long g; rng default;
addpath tests data functions

%% Settings

%...Figure settings
saveFigure = false;
[figSizeLarge,~,figSizeSmall] = saveFigureSettings(saveFigure);

%...Result repository
TudatApplicationOutput = '/Users/Michele/GitHub/tudat/tudatBundle/tudatApplications/Test/SimulationOutput';
repository = fullfile(TudatApplicationOutput,'Propagators');

%...Simulation values
R = 6378.1363;
simulationDuration = 10.0;
scaling = 24*3600;
timeLabel = 'Time [day]';
timeStepSeconds = 50;
constantStepSizes = [20.0, 30.0, 40.0, 50.0, 75.0, 100.0, 150.0, 200.0, 250.0, 300.0];
thrustValues = 5 * 10.^(-3:1);
repository = fullfile(repository,'low_thrust');

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
    for l = 1:length(thrustValues)
        for k = 1:values(j)
            for i = 1:length(propagators)
                if ~( i == length(propagators) && j == length(integrators) ) && ...
                        ~( i == length(propagators) && k ~= 1 )
                    %...Get orbital data
                    fileName = fullfile(repository,['orbit_',...
                        propagators{i},'_',integrators{j},'_',num2str(l),'_',num2str(k),'.dat']);
                    fileID = fopen(fileName);
                    kepler = textscan(fileID,repmat('%f ',[1,7]),'CollectOutput',true,'Delimiter',',');
                    time = (kepler{1}(:,1)-kepler{1}(1))/scaling; kepler = kepler{1}(:,2:end);
                    kepler(:,1) = kepler(:,1)/1e3; kepler(:,3:end) = rad2deg(kepler(:,3:end));
                    fclose(fileID);
                    
                    %...Get trajecotry data
                    fileName = fullfile(repository,['trajectory_',...
                        propagators{i},'_',integrators{j},'_',num2str(l),'_',num2str(k),'.dat']);
                    fileID = fopen(fileName);
                    cartesian = textscan(fileID,repmat('%f ',[1,8]),'CollectOutput',true,'Delimiter',',');
                    cartesian = cartesian{1}(:,2:7); cartesian = cartesian/1e3;
                    fclose(fileID);
                    
                    %...Get time step data
                    fileName = fullfile(repository,['eval_',...
                        propagators{i},'_',integrators{j},'_',num2str(l),'_',num2str(k),'.dat']);
                    fileID = fopen(fileName);
                    m = 1; if i == length(propagators), m = 2; end
                    evaluations = textscan(fileID,repmat('%f ',[1,2]),'CollectOutput',true,'Delimiter',',');
                    evaluations = evaluations{1}; evaluations(:,1) = (evaluations(:,1)-evaluations(1,1))/scaling;
                    evaluations = evaluations(ismembertol(evaluations(:,1),time),2);
                    evaluations = vertcat(evaluations(1),diff(evaluations));
                    evaluations = round( evaluations(1:end-1) / standardFunctionEvals(j,m) );
                    fclose(fileID);
                    
                    %...Interpolate data
                    timeConstantStepSize = 0:timeStepSeconds/scaling:simulationDuration;
                    results(i,j,k,l).prop = propagatorNames{i};
                    results(i,j,k,l).int = integratorNames{j};
                    results(i,j,k,l).val = k;
                    results(i,j,k,l).timeOut = time;
                    results(i,j,k,l).keplerOut = kepler;
                    results(i,j,k,l).cartesianOut = cartesian;
                    results(i,j,k,l).kepler = interp1(time,kepler,timeConstantStepSize,'spline',NaN);
                    results(i,j,k,l).cartesian = interp1(time,cartesian,timeConstantStepSize,'spline',NaN);
                    results(i,j,k,l).evaluations = evaluations;
                end
            end
        end
    end
end
clear fileName fileID k kepler cartesian evaluations time

%% Plot Reference Trajectory

ref = length(propagators);
reference = results(ref,1,1,end);
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

%% Compare Results of Cowell Propagator and USM

%...Difference in Cartesian elements
rmsError = zeros(length(propagators),8,length(integrators),max(values),length(thrustValues));
for k = 1:length(integrators)
    for m = 1:length(thrustValues)
        for l = 1:values(k)
            for i = 1:8
                for j = 1:length(propagators)-1
                    if i < 4
                        rmsError(j,i,k,l,m) = rms(results(j,k,l,m).cartesian(:,i) - results(end,1,1,m).cartesian(:,i));
                    elseif i > 4 && i < 8
                        h = i - 1;
                        rmsError(j,i,k,l,m) = rms(results(j,k,l,m).cartesian(:,h) - results(end,1,1,m).cartesian(:,h));
                    elseif i == 4
                        rmsError(j,i,k,l,m) = rms(sqrt(sum(results(j,k,l,m).cartesian(:,1:3).^2,2)) - ...
                            sqrt(sum(results(end,1,1,m).cartesian(:,1:3).^2,2)));
                    else
                        rmsError(j,i,k,l,m) = rms(sqrt(sum(results(j,k,l,m).cartesian(:,4:6).^2,2)) - ...
                            sqrt(sum(results(end,1,1,m).cartesian(:,4:6).^2,2)));
                    end
                end
            end
        end
    end
end

%% Root-mean-squared Error

styles = {'-o','-d','-s','-v','-p','-h','-*','-x','-^','-o','-d'};

%...Plot RMS error for variable step size
for m = 1:length(thrustValues)
    F = figure('rend','painters','pos',figSizeSmall);
    hold on
    for i = 1:length(propagators)-1
        plot(arrayfun(@(j)sum(results(i,1,j,m).evaluations),1:values(1)),...
            squeeze(rmsError(i,4,1,1:values(1),m))*1e3,styles{i},'LineWidth',1.5,'MarkerSize',10)
    end
    hold off
    xlabel('Function Evaluations [-]')
    ylabel('RMS Position Error [m]')
    legend(propagatorNames(1:end-1),'Location','Best')
    set(gca,'FontSize',15,'YScale','log')
    grid on
    if saveFigure, saveas(F,['../../Report/figures/rms_var_thrust_',num2str(m)],'epsc'),
    else, title(['Thrust ',num2str(thrustValues(m))]), end
    
    %...Plot RMS error for constant step size
    F = figure('rend','painters','pos',figSizeSmall);
    hold on
    for i = 1:length(propagators)-1
        plot(constantStepSizes,squeeze(rmsError(i,4,2,1:values(2),m))*1e3,styles{i},'LineWidth',1.5,'MarkerSize',10)
    end
    hold off
    xlabel('Constant Time Step [s]')
    ylabel('RMS Position Error [m]')
    legend(propagatorNames(1:end-1),'Location','Best')
    set(gca,'FontSize',15,'YScale','log')
    grid on
    if saveFigure, saveas(F,['../../Report/figures/rms_const_thrust_',num2str(m)],'epsc'),
    else, title(['Thrust ',num2str(thrustValues(m))]), end
end