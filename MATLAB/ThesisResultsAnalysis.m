fclose('all'); clear all; close all force; profile off; clc; format long g; rng default;
addpath functions

%% Settings

%...Figure settings
showFigures = false;
saveFigures = false;
[figSizeLarge,figSizeMedium,figSizeSmall] = saveFigureSettings(saveFigures);
figSizeWideLAR = [150,150,1250,475]; % low aspect ratio
figSizeWideHAR = [150,150,1250,400]; % high aspect ratio

%...Constants
marsRadius = 3389526.666666667;
marsGravitationalParameter = 42828375815756.1;
marsAtmosphericInterface = 200;
onboardAtmosphereModel = 0; % 0: exp, 1: wave, 2: 5-param
spacecraftInertia = zeros(3,3);
spacecraftInertia(1,1) = 5750.0;
spacecraftInertia(2,2) = 1215.0;
spacecraftInertia(3,3) = 5210.0;

%...Output folder
outputFolder = 'SimulationOutput';
% outputFolder = [outputFolder,'/low'];
% outputFolder = [outputFolder,'TransOnly'];
outputFolder = [outputFolder,'TransOnly/full'];
% outputFolder = [outputFolder,'TransOnly/server'];
% outputFolder = [outputFolder,'TransOnlyReduced'];
% outputFolder = [outputFolder,'TransGuidOnly/no_sun'];
% outputFolder = [outputFolder,'TransOnlyIMAN/high_ecc_double'];
% outputFolder = [outputFolder,'TransOnlyIMANLoop/0-0--1/'];
% outputFolder = [outputFolder,'TransOnlyIMANRMS/high_ecc/6/100'];

%...Plot settings
loadRotational = false;
loadDependent = true;
loadEstimated = true;
loadFilter = false;
loadMeasurements = false;
applyInterpolation = true;

%...Labels
timeConversion = 3600 * 24;
timeLabel = 'Time [d]';
CartesianLabels = {'x [km]','y [km]','z [km]','v_x [m s^{-1}]','v_y [m s^{-1}]','v_z [m s^{-1}]'};
CartesianLabelsDifference = {'\Delta x [km]','\Delta y [km]','\Delta z [km]',...
    '\Delta v_x [m s^{-1}]','\Delta v_y [m s^{-1}]','\Delta v_z [m s^{-1}]'};
CartesianLabelsAbsoluteDifference = {'|\Delta x| [km]','|\Delta y| [km]','|\Delta z| [km]',...
    '|\Delta v_x| [m s^{-1}]','|\Delta v_y| [m s^{-1}]','|\Delta v_z| [m s^{-1}]'};
KeplerianLabels = {'a [km]','e [-]','i [deg]','\omega [deg]','\Omega [deg]','\vartheta [deg]'};
KeplerianLabelsDifference = {'\Delta a [km]','\Delta e [-]','\Delta i [deg]',...
    '\Delta\omega [deg]','\Delta\Omega [deg]','\Delta\vartheta [deg]'};
rotationLabels = {'\eta [-]','\epsilon_1 [-]','\epsilon_2 [-]','\epsilon_3 [-]','Norm Offset [-]',...
    '\omega_1 [deg s^{-1}]','\omega_2 [deg s^{-1}]','\omega_3 [deg s^{-1}]'};
rotationLabelsDifference = {'\Delta\eta [-]','\Delta\epsilon_1 [-]','\Delta\epsilon_2 [-]','\Delta\epsilon_3 [-]'};

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

locOmega = KeplerianPropagatedResults(:,4) > 180;
KeplerianPropagatedResults(locOmega,4) = KeplerianPropagatedResults(locOmega,4) - 360;

originalKeplerianPropagatedResults = KeplerianPropagatedResults; % create copy 

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

loc = diff(onboardTime) == 0;
onboardTime(loc) = [];
CartesianEstimatedResults(loc,:) = [];
KeplerianEstimatedResults(loc,:) = [];

%...Load rotational data
if loadRotational
    %...Load rotational motion
    filename = ['/Users/Michele/GitHub/tudat/tudatBundle/tudatApplications/Thesis/',outputFolder,'/rotationalEstimated.dat'];
    fileID = fopen(filename,'r');
    rotationalEstimatedResults = textscan(fileID,repmat('%f',[1,8]),'Delimiter',',','CollectOutput',true);
    rotationalEstimatedResults = rotationalEstimatedResults{1}(:,2:end);
    rotationalEstimatedResults(:,6:8) = rad2deg(rotationalEstimatedResults(:,5:7));
    rotationalEstimatedResults(:,5) = 1.0 - quatnorm(rotationalEstimatedResults(:,1:4));
    fclose(fileID);
    
    %...Load commanded quaternions
    filename = ['/Users/Michele/GitHub/tudat/tudatBundle/tudatApplications/Thesis/',outputFolder,'/commandedQuaternions.dat'];
    fileID = fopen(filename,'r');
    commandedQuaternions = textscan(fileID,repmat('%f',[1,5]),'Delimiter',',','CollectOutput',true);
    commandedQuaternions = commandedQuaternions{1};
    fclose(fileID);
end

%...Clean up
clear filename fileID

%% Load C++ Results For Navigation Filter

%...Only if filtering is toggled
if loadFilter
    if contains(outputFolder,'IMAN')
        fliterLimit = 11;
    elseif loadRotational
        fliterLimit = 17;
    else
        fliterLimit = 10;
    end
    
    %...Load translational motion
    filename = ['/Users/Michele/GitHub/tudat/tudatBundle/tudatApplications/Thesis/',outputFolder,'/filterStateEstimates.dat'];
    fileID = fopen(filename,'r');
    filterStateEstimatedResults = textscan(fileID,repmat('%f',[1,fliterLimit]),'Delimiter',',','CollectOutput',true);
    filterTime = ( filterStateEstimatedResults{1}(:,1) - initialTime ) / timeConversion;
    filterStateEstimatedResults = filterStateEstimatedResults{1}(:,2:end);
    filterStateEstimatedResults(:,1:3) = filterStateEstimatedResults(:,1:3)/1e3;
    fclose(fileID);
    
    %...Load Keplerian translational motion
    filename = ['/Users/Michele/GitHub/tudat/tudatBundle/tudatApplications/Thesis/',outputFolder,...
        '/filterCovarianceEstimates.dat'];
    fileID = fopen(filename,'r');
    filterCovarianceEstimatedResults = textscan(fileID,repmat('%f',[1,fliterLimit]),'Delimiter',',','CollectOutput',true);
    filterCovarianceEstimatedResults = filterCovarianceEstimatedResults{1}(:,2:end);
    filterCovarianceEstimatedResults(:,1:3) = filterCovarianceEstimatedResults(:,1:3)/1e3^2;
    fclose(fileID);
    
    %...Clean up
    clear filename fileID
end

%% Load C++ Results For Dependent Variables
    
%...Load control torques
if loadRotational
    filename = ['/Users/Michele/GitHub/tudat/tudatBundle/tudatApplications/Thesis/',outputFolder,'/controlTorques.dat'];
    fileID = fopen(filename,'r');
    controlTorques = textscan(fileID,repmat('%f',[1,4]),'Delimiter',',','CollectOutput',true);
    controlTorques = controlTorques{1}(:,2:end);
    fclose(fileID);
end

%...Load dependent variables
if loadDependent
    filename = ['/Users/Michele/GitHub/tudat/tudatBundle/tudatApplications/Thesis/',outputFolder,'/dependentVariables.dat'];
    fileID = fopen(filename,'r');
    dependentVariables = textscan(fileID,repmat('%f',[1,16]),'Delimiter',',','CollectOutput',true);
    dependentVariables = dependentVariables{1}(:,2:end); dependentVariables(:,7:9) = rad2deg(dependentVariables(:,7:9));
    fclose(fileID);
    
%     locAlpha = dependentVariables(:,7) > 0 & simulationTime > 0.8;
%     dependentVariables(locAlpha,7) = dependentVariables(locAlpha,7) - 360;
    locSigma = dependentVariables(:,9) > 0;
    dependentVariables(locSigma,9) = dependentVariables(locSigma,9) - 360;
else
    dependentVariables = NaN(length(simulationTime),12);
end

if loadEstimated
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
    
    %...Load PTE estimates
    filename = ['/Users/Michele/GitHub/tudat/tudatBundle/tudatApplications/Thesis/',outputFolder,'/pteEstimates.dat'];
    fileID = fopen(filename,'r');
    pteEstimates = textscan(fileID,repmat('%f',[1,7]),'Delimiter',',','CollectOutput',true);
    pteEstimates = pteEstimates{1};
    fclose(fileID);
end
    
%...Clean up
clear filename fileID limit

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
    filename = ['/Users/Michele/GitHub/tudat/tudatBundle/tudatApplications/Thesis/',outputFolder,'/expectedMeasurements.dat'];
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
        %...Interpolate
        filterStateEstimatedResults = interp1( filterTime, filterStateEstimatedResults, interpolatedTime, 'linear', NaN );
        filterCovarianceEstimatedResults = interp1( filterTime, filterCovarianceEstimatedResults, ...
            interpolatedTime, 'linear', NaN );
        
        %...Add NaNs between non-filter moments
        locAddNan = ~ismembertol(interpolatedTime,filterTime,2.5e-5);
        filterStateEstimatedResults(locAddNan,:) = NaN;
        filterCovarianceEstimatedResults(locAddNan,:) = NaN;
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
[x,y,z] = sphere; surf(marsRadius/1e3*x,marsRadius/1e3*y,marsRadius/1e3*z)
hold off
xlabel('x [km]'), ylabel('y [km]'), zlabel('z [km]')
view([0,0])
grid on
axis equal tight
set(gca,'FontSize',17.5)
if saveFigures, saveas(F,'../../Report/figures/aero_traj','epsc'), end

%% Plot Apses Altitudes

if loadEstimated
    %...Retrieve apo- and periapses altitudes
    actualAltitude = sqrt( sum( CartesianPropagatedResults(:,1:3).^2, 2 ) ) - marsRadius/1e3;
    [actualApoapses,locActualApoapses] = findpeaks(actualAltitude,'MinPeakProminence',100);
    actualPeriapses = -findpeaks(-actualAltitude,'MinPeakProminence',100);
    
    estimatedAltitude = sqrt( sum( CartesianEstimatedResults(:,1:3).^2, 2 ) ) - marsRadius/1e3;
    estimatedApoapses = findpeaks(estimatedAltitude,'MinPeakProminence',100);
    estimatedPeriapses = -findpeaks(-estimatedAltitude,'MinPeakProminence',100);
    
    if ~isempty(actualPeriapses)
        %...Analyze maneuvers
        maneuverLoc = apoapsisManeuver(:,1) + 0.5;
        maneuverMagnitude = apoapsisManeuver(:,2);
        
        locExtra = maneuverLoc > length(actualPeriapses);
        maneuverMangitude = maneuverMagnitude(~locExtra);
        maneuverLoc = maneuverLoc(~locExtra);
        
        maneuverValue = actualPeriapses(maneuverLoc(2:end)-0.5);
        maneuverValue = [actualPeriapses(1);maneuverValue];
        
        locUp = apoapsisManeuver(1:end-2,2) > 0;
        upManeuverLoc = maneuverLoc(locUp);
        upManeuverValue = maneuverValue(locUp);
        
        locDown = apoapsisManeuver(1:end-2,2) < 0;
        downManeuverLoc = maneuverLoc(locDown);
        downManeuverValue = maneuverValue(locDown);
        
        %...Analyze periapses
        periapsesLocs = 1:length(actualPeriapses);
        locFail = ( actualPeriapses < periapsisCorridors(1:end-2,2)/1e3 ) | ...
            ( actualPeriapses > periapsisCorridors(1:end-2,3)/1e3 );
%         locFail = ( actualPeriapses < periapsisCorridors(:,2)/1e3 ) | ( actualPeriapses > periapsisCorridors(:,3)/1e3 );
        periapsesLocsFail = periapsesLocs(locFail); actualPeriapsesFail = actualPeriapses(locFail);
        periapsesLocs(locFail) = []; actualPeriapses(locFail) = [];
        
        %...Plot periapses
        F = figure('rend','painters','pos',figSizeWideLAR);
        hold on
        scatter(periapsesLocs,actualPeriapses,50)
        scatter((1:length(estimatedPeriapses)),estimatedPeriapses,50,'s')
        scatter(upManeuverLoc,upManeuverValue,50,'^','filled')
        scatter(downManeuverLoc,downManeuverValue,50,'v','filled')
        plot(periapsisCorridors(2:end,1),periapsisCorridors(1:end-1,2)/1e3,'LineWidth',1.25,'LineStyle','--')
        plot(periapsisCorridors(2:end,1),periapsisCorridors(1:end-1,3)/1e3,'LineWidth',1.25,'LineStyle','-.')
        scatter(periapsesLocsFail,actualPeriapsesFail,50,[0,0.447,0.741],'o','filled')
        hold off
        xlabel('Orbit Number [-]')
        ylabel('Altitude [km]')
        xlim([0,210])
        set(gca,'FontSize',15)
        [~,icons] = legend('Actual','Estimated','Up Maneuver','Down Maneuver','Lower Bound','Upper Bound','Location','Best');
        for i = 1:14
            if i >= 7 && i <= 10
                icons(i).Children.MarkerSize = 10;
            end
        end
        grid on
        if saveFigures
            saveas(F,'../../Report/figures/aero_peri_corr_full','epsc')
            xlim([110,160])
            saveas(F,'../../Report/figures/aero_peri_corr_term','epsc')
        end
    end
end

%% Plot PTE Estimates

if loadEstimated
    modLocActualApoapses = [1;locActualApoapses];
    
    %...Find difference between semi-major axis and eccentricity between apoapses
    semiMajorAxesPTE = KeplerianPropagatedResults(modLocActualApoapses+1,1);
    eccentricitiesPTE = KeplerianPropagatedResults(modLocActualApoapses+1,2);
    
    diffSMAPTE = diff(semiMajorAxesPTE);
    diffEccPTE = diff(eccentricitiesPTE);
    
    %...Verify
    figure
    subplot(1,2,1)
    hold on
    scatter(interpolatedTime(modLocActualApoapses+1),semiMajorAxesPTE)
    plot(interpolatedTime,KeplerianPropagatedResults(:,1))
    hold off
    subplot(1,2,2)
    hold on
    scatter(interpolatedTime(modLocActualApoapses+1),eccentricitiesPTE)
    plot(interpolatedTime,KeplerianPropagatedResults(:,2))
    hold off
    
    %...Plot compared to PTE estimates
    F = figure('rend','painters','pos',figSizeSmall);
    hold on
    scatter(1:length(diffSMAPTE),-diffSMAPTE,50)
    scatter(pteEstimates(1:end-1,1),-pteEstimates(1:end-1,2)/1e3,50,'s')
    hold off
    xlabel('Orbit Number [-]')
    ylabel('- \Delta a [km]')
    legend('Actual','Estimated','Location','Best')
    set(gca,'FontSize',15,'YScale','log')
    grid on
    if saveFigures, saveas(F,'../../Report/figures/vv_pte_sma','epsc'), end
    
    F = figure('rend','painters','pos',figSizeSmall);
    hold on
    scatter(1:length(diffEccPTE),-diffEccPTE,50)
    scatter(pteEstimates(1:end-1,1),-pteEstimates(1:end-1,3),50,'s')
    hold off
    xlabel('Orbit Number [-]')
    ylabel('- \Delta e [-]')
    legend('Actual','Estimated','Location','Best')
    set(gca,'FontSize',15,'YScale','log')
    grid on
    if saveFigures, saveas(F,'../../Report/figures/vv_pte_ecc','epsc'), end    
end

%% Analyze Maneuvers

if loadEstimated
    reducedManeuver = apoapsisManeuver(1:end-2,:); % exclude termination and aerobraking complete phases
    
    locWalkIn = reducedManeuver(:,1) < 7;
    locWalkOut = KeplerianEstimatedResults(modLocActualApoapses(reducedManeuver(:,1)+1),2) < 0.3;
    locMain = ~locWalkIn & ~locWalkOut;
    
    walkInDV = sum(abs(reducedManeuver(locWalkIn,2)));
    mainDV = sum(abs(reducedManeuver(locMain,2)));
    walkOutDV = sum(abs(reducedManeuver(locWalkOut,2)));
    terminationDV = sum(abs(apoapsisManeuver(end-1:end,2)));
    totDV = sum(abs(apoapsisManeuver(:,2)));
    DeltaVs = [walkInDV,mainDV,walkOutDV,terminationDV];
    
    table(walkInDV,mainDV,walkOutDV,terminationDV,totDV)
    fprintf('\\num{%.2f} & \\num{%.2f} & \\num{%.1f} & \\num{%.1f} \\\\\n',DeltaVs)
end

%% Plot Target vs. Actual Periapsis

if loadEstimated
    %...Find periapses
    [actualPeriapses,locActualPeriapses] = findpeaks(-actualAltitude,'MinPeakProminence',100);
    actualPeriapses = -actualPeriapses;
    
    locTargerPeriapses = reducedManeuver(:,1)+1;
    actualTargetedPeriapses = actualPeriapses(locTargerPeriapses);
    targetPeriapses = mean(periapsisCorridors(locTargerPeriapses,2:3),2)/1e3;
    
    %...Plot target vs actual
    F = figure('rend','painters','pos',figSizeWideHAR);
    hold on
    scatter(reducedManeuver(:,1)+1,actualTargetedPeriapses,50)
    scatter(reducedManeuver(:,1)+1,targetPeriapses,50,'s')
    plot(periapsisCorridors(2:end,1),periapsisCorridors(1:end-1,2)/1e3,'LineWidth',1.25,'LineStyle','--')
    plot(periapsisCorridors(2:end,1),periapsisCorridors(1:end-1,3)/1e3,'LineWidth',1.25,'LineStyle','-.')
    hold off
    xlabel('Orbit Number [-]')
    ylabel('Altitude [km]')
    legend('Actual','Target','Lower Bound','Upper Bound','Location','Best')
    set(gca,'FontSize',15)
    grid on
    if saveFigures, saveas(F,'../../Report/figures/vv_me','epsc'), end 
end

%% Plot Densities

if loadEstimated
    %...Find maximum density places
    [maximumDensityPerOrbit,locMaxDens] = findpeaks(dependentVariables(:,10),'MinPeakHeight',1.0e-10);
    estimatedMaxDensity = zeros(size(atmosphericParameters,1),1);
    meanEstMaxDensity = zeros(size(atmosphericParameters,1),1);
    
    %...Moving average
    meanAP = [smooth(atmosphericParameters(:,2),7),smooth(atmosphericParameters(:,3),7),...
        smooth(atmosphericParameters(:,4),7)];
    mean_dens_func = @(i,h) meanAP(i,2) * exp( - ( h*1e3 - meanAP(i,1) ) / meanAP(i,3) );
    
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
        estimatedMaxDensity(i) = dens_func(actualAltitude(locMaxDens(i)));
        meanEstMaxDensity(i) = mean_dens_func(i,actualAltitude(locMaxDens(i)));
        plot(dens_func(reducedAltitude),reducedAltitude,'LineWidth',1.25)
    end
    hold off
    xlabel('Density [kg m^{-3}]')
    ylabel('Altitude [km]')
    ylim([125,175])
    grid on
    set(gca,'FontSize',15,'XScale','log')
    
    %...Plot real and estimated peak densities
    F = figure('rend','painters','pos',figSizeSmall);
    hold on
    scatter((1:length(maximumDensityPerOrbit)),maximumDensityPerOrbit,50)
    scatter(atmosphericParameters(:,1),estimatedMaxDensity,50,'s') 
    scatter(atmosphericParameters(:,1),meanEstMaxDensity,50,'^')
    hold off
    xlabel('Orbit Number [-]')
    ylabel('Density [kg m^{-3}]')
    legend('Actual','Estimated','Moving Average','Location','Best')
    set(gca,'FontSize',15,'YScale','log')
    grid on
    if saveFigures, saveas(F,'../../Report/figures/aero_peak_dens','epsc'), end
    
    densityOffset = (meanEstMaxDensity-maximumDensityPerOrbit)./maximumDensityPerOrbit*100;
    
    F = figure('rend','painters','pos',figSizeSmall);
    hold on
    scatter((1:length(maximumDensityPerOrbit)),densityOffset,50)
    hold off
    xlabel('Orbit Number [-]')
    ylabel('Offset [%]')
    set(gca,'FontSize',15)
    grid on
%     if saveFigure, saveas(F,'../../Report/figures/aero_peak_dens','epsc'), end
    
    minDensityError = min((estimatedMaxDensity-maximumDensityPerOrbit)./maximumDensityPerOrbit*100);
    maxDensityError = max((estimatedMaxDensity-maximumDensityPerOrbit)./maximumDensityPerOrbit*100);
    table(minDensityError,maxDensityError)
    
    meanMinDensityError = min((meanEstMaxDensity-maximumDensityPerOrbit)./maximumDensityPerOrbit*100);
    meanMaxDensityError = max((meanEstMaxDensity-maximumDensityPerOrbit)./maximumDensityPerOrbit*100);
end

%% RMS Error

%...Compute RMS error in position and velocity
rmsPositionError = rms( rssq(CartesianEstimatedResults(:,1:3),2) - rssq(CartesianPropagatedResults(:,1:3),2) ) * 1e3;
rmsVelocityError = rms( rssq(CartesianEstimatedResults(:,4:6),2) - rssq(CartesianPropagatedResults(:,4:6),2) );

%...Show errors
table(rmsPositionError,rmsVelocityError)

%...Compute STD of error in position and velocity for aided navigation
locTime = interpolatedTime < 0.5;
stdErrorPositionAided = std( rssq(CartesianEstimatedResults(locTime,1:3),2) - ...
    rssq(CartesianPropagatedResults(locTime,1:3),2) ) * 1e3;
stdErrorVelocityAided = std( rssq(CartesianEstimatedResults(locTime,4:6),2) - ...
    rssq(CartesianPropagatedResults(locTime,4:6),2) );
% stdErrorPositionAided = std( CartesianEstimatedResults(locTime,1:3) - CartesianPropagatedResults(locTime,1:3) ) * 1e3;
% stdErrorVelocityAided = std( CartesianEstimatedResults(locTime,4:6) - CartesianPropagatedResults(locTime,4:6) );

%...Show errors
table(stdErrorPositionAided,stdErrorVelocityAided)

%...Compute STD of error in position and velocity for unaided navigation
locTime = ( interpolatedTime > 0.5448 & interpolatedTime < 0.6872 ) | ...
    ( interpolatedTime > 0.7822 & interpolatedTime < 0.9236 );
stdErrorPositionUnaided = std( rssq(CartesianEstimatedResults(locTime,1:3),2) - ...
    rssq(CartesianPropagatedResults(locTime,1:3),2) ) * 1e3;
stdErrorVelocityUnaided = std( rssq(CartesianEstimatedResults(locTime,4:6),2) - ...
    rssq(CartesianPropagatedResults(locTime,4:6),2) );

%...Show errors
table(stdErrorPositionUnaided,stdErrorVelocityUnaided)

%...Rotational velocity and angular momentum
if loadRotational
    angularMomentum = spacecraftInertia * deg2rad(rotationalPropagatedResults(:,6:8))';
    maxAngularMomentum = max(max(angularMomentum,[],2));
    table(maxAngularMomentum)
end

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

if saveFigures
    F = figure('rend','painters','pos',figSizeWideHAR);
    for i = 1:2
        subplot(1,2,i)
        plot(interpolatedTime,KeplerianPropagatedResults(:,i),'LineWidth',1.5)
        xlabel(timeLabel)
        ylabel(KeplerianLabels{i})
        set(gca,'FontSize',16.5)
        grid on
    end
    if saveFigures, saveas(F,'../../Report/figures/aero_kepl_ae','epsc'), end
    
    F = figure('rend','painters','pos',figSizeWideHAR);
    for i = 3:5
        subplot(1,3,i-2)
        plot(interpolatedTime,KeplerianPropagatedResults(:,i),'LineWidth',1.5)
        xlabel(timeLabel)
        ylabel(KeplerianLabels{i})
        set(gca,'FontSize',16.5)
        grid on
    end
    if saveFigures, saveas(F,'../../Report/figures/aero_kepl_iOo','epsc'), end
end

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
        if i < 5
            plot([interpolatedTime(1),interpolatedTime(end)],commandedQuaternions(1,i+1)*ones(1,2),...
                'LineWidth',1.25,'LineStyle','--')
        end
        hold off
        xlabel(timeLabel)
        ylabel(rotationLabels{i})
        set(gca,'FontSize',15)
        grid on
    end
    subplotLegend({'Actual','Estimated','Commanded'})
end

if loadRotational && saveFigures
    F = figure('rend','painters','pos',figSizeSmall);
    for i = 1:4
        subplot(2,2,i)
        plot(interpolatedTime,rotationalPropagatedResults(:,i)-commandedQuaternions(1,i+1),'LineWidth',1.25)
%         plot(interpolatedTime,rotationalEstimatedResults(:,i)-rotationalPropagatedResults(:,i),'LineWidth',1.25)
        xlabel(timeLabel)
        ylabel(rotationLabelsDifference{i})
        set(gca,'FontSize',15)
        grid on
    end
    if saveFigures, saveas(F,'../../Report/figures/aero_rot_diff','epsc'), end
    
    F = figure('rend','painters','pos',figSizeSmall);
    for i = 1:3
        subplot(3,1,i)
        plot(interpolatedTime,rotationalPropagatedResults(:,i+5),'LineWidth',1.25)
        xlabel(timeLabel)
        ylabel(rotationLabels{i+5})
        set(gca,'FontSize',15)
        grid on
    end
    if saveFigures, saveas(F,'../../Report/figures/aero_rot_vel','epsc'), end
end

%% Plot Peri- and Apoapsis During Aerobraking

%...Compute apsis altitudes
apoapsisAltitude = KeplerianPropagatedResults(:,1) .* ...
    (1.0 + KeplerianPropagatedResults(:,2)) - marsRadius/1e3;
periapsisAltitude = KeplerianPropagatedResults(:,1) .* ...
    (1.0 - KeplerianPropagatedResults(:,2)) - marsRadius/1e3;

%...Compute DAIA
percentage = 0.5 * ( 1 + ( 1 - KeplerianPropagatedResults(:,1) / originalKeplerianPropagatedResults(1) ) );
dynamicAtmosphericInterfaceAltitude = percentage .* ( KeplerianPropagatedResults(:,1) - marsRadius/1e3 );

%...Plot convergence during aerobraking
F = figure('rend','painters','pos',figSizeWideLAR);
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
if saveFigures, saveas(F,'../../Report/figures/aero_apo_peri','epsc'), end

%% Plot Peri- and Apoapsis After Aerobraking

locStable = ( simulationTime > ( interpolatedTime(end) - 0.0385 ) ) & ( simulationTime < interpolatedTime(end) );
if any(locStable)
    stableTime = simulationTime(locStable); 
    initialStableTime = stableTime(1);
    stableTime = stableTime - initialStableTime;
    stableKeplerianElements = originalKeplerianPropagatedResults(locStable,:);
    
    %...Compute apsis altitudes
    stableApoapsisAltitude = stableKeplerianElements(:,1) .* ...
        (1.0 + stableKeplerianElements(:,2)) - marsRadius/1e3;
    stablePeriapsisAltitude = stableKeplerianElements(:,1) .* ...
        (1.0 - stableKeplerianElements(:,2)) - marsRadius/1e3;
    
    %...Compute DAIA
    percentage = 0.5 * ( 1 + ( 1 - stableKeplerianElements(:,1) / originalKeplerianPropagatedResults(1) ) );
    stableDynamicAtmosphericInterfaceAltitude = percentage .* ( stableKeplerianElements(:,1) - marsRadius/1e3 );
    
    %...Compute altitude
    stableAltitude = stableKeplerianElements(:,1) .* ( 1.0 - stableKeplerianElements(:,2).^2 ) ./ ...
        ( 1.0 + stableKeplerianElements(:,2) .* cosd(stableKeplerianElements(:,6)) ) - marsRadius/1e3;
    
    %...Plot stable conditions after aerobraking
    F = figure('rend','painters','pos',figSizeSmall);
    hold on
    plot(stableTime,stableApoapsisAltitude,'LineWidth',1.25)
    plot(stableTime,stablePeriapsisAltitude,'LineWidth',1.25)
    plot(stableTime,stableAltitude,'LineWidth',1.25,'LineStyle','--')
    plot(stableTime,stableDynamicAtmosphericInterfaceAltitude,'LineWidth',1.25,'LineStyle',':')
    hold off
    xlabel(timeLabel)
    ylabel('Altitude [km]')
    legend('Apoapsis','Periapsis','Altitude','DAIA','Location','Best')
    grid on
    set(gca,'FontSize',15,'YScale','log')
    if saveFigures, saveas(F,'../../Report/figures/aero_apo_peri_stable','epsc'), end
    
    F = figure('rend','painters','pos',figSizeLarge);
    for i = 1:size(KeplerianPropagatedResults,2)
        subplot(2,3,i)
        plot(stableTime,stableKeplerianElements(:,i),'LineWidth',1.25)
        xlabel(timeLabel)
        ylabel(KeplerianLabels{i})
        set(gca,'FontSize',15)
        grid on
    end
    
    %...Mean elements
    stableKeplerianElements(stableKeplerianElements(:,4)>180,4) = ...
        stableKeplerianElements(stableKeplerianElements(:,4)>180,4) - 360;
    mean([stableKeplerianElements(:,1:5),stableApoapsisAltitude,stablePeriapsisAltitude])
    
    %...Offset
    offsetApoapsis = ( mean(stableApoapsisAltitude) - 320 ) / 320 * 100;
    offsetPeriapsis = ( mean(stablePeriapsisAltitude) - 255 ) / 255 * 100;
end

%% Plot Filter States

%...Only if filtering is toggled
saveIMAN = false;
saveRot = false;
if loadFilter
    plotFreq = 1;
    if saveRot, plotFreq = 10; end
    
    %...Plot Cartesian translational motion
    F = figure('rend','painters','pos',figSizeLarge);
    for i = 1:6
        subplot(2,3,i)
        hold on
        if saveIMAN
            plot(interpolatedTime(1:plotFreq:end),abs(filterStateEstimatedResults(1:plotFreq:end,i)-...
                CartesianPropagatedResults(1:plotFreq:end,i)),'LineWidth',1.25)
            plot(interpolatedTime(3:plotFreq:end),sqrt(filterCovarianceEstimatedResults(3:plotFreq:end,i)),...
                'LineWidth',1.25,'LineStyle','--')
        else
            plot(interpolatedTime(1:plotFreq:end),filterStateEstimatedResults(1:plotFreq:end,i)-...
                CartesianPropagatedResults(1:plotFreq:end,i),'LineWidth',1.25)
            plot(interpolatedTime(3:plotFreq:end),sqrt(filterCovarianceEstimatedResults(3:plotFreq:end,i)),...
                'LineWidth',1.25,'LineStyle','--')
            plot(interpolatedTime(3:plotFreq:end),-sqrt(filterCovarianceEstimatedResults(3:plotFreq:end,i)),...
                'LineWidth',1.25,'LineStyle','--','Color',[0.85,0.325,0.098])
        end
        hold off
        xlabel(timeLabel)
        if saveIMAN
            ylabel(CartesianLabelsAbsoluteDifference{i})
        else
            ylabel(CartesianLabelsDifference{i})
        end
        if saveIMAN
            set(gca,'FontSize',15,'YScale','log')
        else
            set(gca,'FontSize',15)
        end
        grid on
    end
    subplotLegend({'Difference','STD'})
    if ( saveFigures && saveIMAN ), saveas(F,'../../Report/figures/aero_filt_cart','epsc'),
    elseif ( saveFigures && saveRot ), saveas(F,'../../Report/figures/aero_withrot_cart','epsc'), end
    
    %...Plot instrument errors
    actual = [0.0001239735903658918  5.231874102765194e-06   9.72210316125009e-06];
    F = figure('rend','painters','pos',figSizeMedium);
    for i = 7:9
        subplot(1,3,i-6)
        hold on
        if saveIMAN
            plot(interpolatedTime(1:plotFreq:end),abs(filterStateEstimatedResults(1:plotFreq:end,i)-actual(i-6)),'LineWidth',1.25)
            plot(interpolatedTime(3:plotFreq:end),sqrt(filterCovarianceEstimatedResults(3:plotFreq:end,i)),...
                'LineWidth',1.25,'LineStyle','--')
        else
            plot(interpolatedTime(1:plotFreq:end),filterStateEstimatedResults(1:plotFreq:end,i)-actual(i-6),'LineWidth',1.25)
            plot(interpolatedTime(3:plotFreq:end),sqrt(filterCovarianceEstimatedResults(3:plotFreq:end,i)),...
                'LineWidth',1.25,'LineStyle','--')
            plot(interpolatedTime(3:plotFreq:end),-sqrt(filterCovarianceEstimatedResults(3:plotFreq:end,i)),...
                'LineWidth',1.25,'LineStyle','--','Color',[0.85,0.325,0.098])
        end
        hold off
        xlabel(timeLabel)
        if saveIMAN
            ylabel(['|\Delta b_{t,',num2str(i-6),'}| [-]'])
        else
            ylabel(['\Delta b_{t,',num2str(i-6),'} [-]'])
        end
        if saveIMAN
            set(gca,'FontSize',15,'YScale','log')
        else
            set(gca,'FontSize',15)
        end
        grid on
    end
    
    %...Plot drag coefficient
    if contains(outputFolder,'IMAN')
        actual = 1.9;
        F = figure('rend','painters','pos',figSizeSmall);
        hold on
        plot(interpolatedTime(1:plotFreq:end),abs(filterStateEstimatedResults(1:plotFreq:end,10)-actual),'LineWidth',1.25)
        plot(interpolatedTime(3:plotFreq:end),sqrt(filterCovarianceEstimatedResults(3:plotFreq:end,10)),...
            'LineWidth',1.25,'LineStyle','--')
        %         plot(interpolatedTime(2:end),-sqrt(filterCovarianceEstimatedResults(2:end,i)),'LineWidth',1.25,'LineStyle','--')
        hold off
        xlabel(timeLabel)
        ylabel('|\Delta C_D| [-]')
        ylim([1e-5,1e10])
        legend('Difference','STD','Location','Best')
        set(gca,'FontSize',15,'YScale','log')
        grid on
        if ( saveFigures && saveIMAN ), saveas(F,'../../Report/figures/aero_filt_cd','epsc'), end
    end
    
    %...Plot rotational motion
    if loadRotational
        %...Plot quaternions
        F = figure('rend','painters','pos',figSizeLarge);
        for i = 10:13
            subplot(2,2,i-9)
            hold on
            plot(interpolatedTime(1:plotFreq:end),filterStateEstimatedResults(1:plotFreq:end,i)-...
                rotationalPropagatedResults(1:plotFreq:end,i-9),'LineWidth',1.25)
            plot(interpolatedTime(3:plotFreq:end),sqrt(filterCovarianceEstimatedResults(3:plotFreq:end,i)),...
                'LineWidth',1.25,'LineStyle','--')
            plot(interpolatedTime(3:plotFreq:end),-sqrt(filterCovarianceEstimatedResults(3:plotFreq:end,i)),...
                'LineWidth',1.25,'LineStyle','--','Color',[0.85,0.325,0.098])
            hold off
            xlabel(timeLabel)
            ylabel(rotationLabelsDifference{i-9})
            set(gca,'FontSize',15)
            grid on
        end
        subplotLegend({'Difference','STD'})
        if ( saveFigures && saveRot ), saveas(F,'../../Report/figures/aero_withrot_rot','epsc'), end
        
        %...Plot instrument errors
        actual = [8.060047168656996e-09  1.998853308544197e-09 -1.143908123010419e-09];
        F = figure('rend','painters','pos',figSizeMedium);
        for i = 14:16
            subplot(1,3,i-13)
            hold on
            plot(interpolatedTime(1:plotFreq:end),filterStateEstimatedResults(1:plotFreq:end,i)-actual(i-13),'LineWidth',1.25)
            plot(interpolatedTime(3:plotFreq:end),sqrt(filterCovarianceEstimatedResults(3:plotFreq:end,i)),...
                'LineWidth',1.25,'LineStyle','--')
            plot(interpolatedTime(3:plotFreq:end),-sqrt(filterCovarianceEstimatedResults(3:plotFreq:end,i)),...
                'LineWidth',1.25,'LineStyle','--','Color',[0.85,0.325,0.098])
            hold off
            xlabel(timeLabel)
            ylabel(['\Delta b_{r,',num2str(i-13),'} [-]'])
            if saveIMAN
                set(gca,'FontSize',15,'YScale','log')
            else
                set(gca,'FontSize',15)
            end
            grid on
        end
    end
end

%% Plot Other Results

%...Only if measurements are toggled
if loadMeasurements && loadDependent
    %...Plot measurements
    F = figure('rend','painters','pos',figSizeLarge);
    
    if loadRotational, subplot(1,2,1), end
    hold on
    for i = 1:3
        scatter(interpolatedTime,accelerometerMeasurements(:,i))
        plot(interpolatedTime,smooth(accelerometerMeasurements(:,i),50),'LineWidth',1.25)
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
    
    if loadRotational
        subplot(1,2,2)
        scatter(interpolatedTime,accelerometerMeasurements(:,4:6),'LineWidth',1.25)
        hold off
        xlabel(timeLabel)
        ylabel('Rotational Velocity [rad s^{-1}]')
        set(gca,'FontSize',15)
        grid on
        legend('x','y','z')
    end
end

%...Plot control torques
if loadRotational
    names = {'x_B','y_B','z_B'};
    styles = {'--','-','-.'};
    
    F = figure('rend','painters','pos',figSizeSmall);
    hold on
    for i = 1:3
        plot(onboardTime,controlTorques(:,i),'LineWidth',1.25,'LineStyle',styles{i})
    end
    hold off
    xlabel(timeLabel)
    ylabel('Torque [N m]')
    set(gca,'FontSize',15)
    grid on
    legend(names{:},'Location','Best')
end

%...Plot aerodynamic angles
if loadDependent && loadRotational
    styles = {'-','--','-.'};
    aeroAngles = dependentVariables(:,7:9);
    
    F = figure('rend','painters','pos',figSizeSmall);
    hold on
    for i = 1:3
        plot(interpolatedTime,aeroAngles(:,i),'LineWidth',1.25,'LineStyle',styles{i})
    end
    hold off
    xlabel(timeLabel)
    ylabel('Aerodynamic Angle [deg]')
    set(gca,'FontSize',15)
    grid on
    legend('Attack','Side-slip','Bank','Location','Best')
    if saveFigures, saveas(F,'../../Report/figures/aero_angles','epsc'), end
end

%...Plot aerodynamic coefficients
if loadDependent && loadRotational
    styles = {'-','--','-.'};
    aeroCoeff = dependentVariables(:,13:15);
    
    F = figure('rend','painters','pos',figSizeSmall);
    hold on
    for i = 1:3
        plot(interpolatedTime,aeroCoeff(:,i),'LineWidth',1.25,'LineStyle',styles{i})
    end
    hold off
    xlabel(timeLabel)
    ylabel('Aerodynamic Coefficient [-]')
    set(gca,'FontSize',15)
    grid on
    legend('Drag','Side','Lift','Location','Best')
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
    F = figure('rend','painters','pos',figSizeSmall);
    hold on
    scatter(1:length(maximumDynamicPressure),maximumDynamicPressure,50)
    plot(xlim,[0.19,0.19],'LineWidth',1.25,'LineStyle','--')
    hold off
    xlabel('Orbit Number [-]')
    ylabel('Dynamic Pressure [N m^{-2}]')
    legend('Dyn. Press.','Min. Dyn. Press.','Location','NE')
    set(gca,'FontSize',15)
    grid on
    if saveFigures, saveas(F,'../../Report/figures/aero_press_cond','epsc'), end
    
    F = figure('rend','painters','pos',figSizeSmall);
    hold on
    scatter(1:length(maximumHeatRate),maximumHeatRate,50)
    plot(xlim,[2800,2800],'LineWidth',1.25,'LineStyle','--')
    hold off
    xlabel('Orbit Number [-]')
    ylabel('Heat Rate [W m^{-2}]')
    legend('Heat Rate','Max. Heat Rate','Location','Best')
    set(gca,'FontSize',15)
    grid on
    if saveFigures, saveas(F,'../../Report/figures/aero_heat_rate','epsc'), end
    
    F = figure('rend','painters','pos',figSizeSmall);
    hold on
    scatter(1:length(maximumHeatLoad),maximumHeatLoad,50)
    plot(xlim,[500,500],'LineWidth',1.25,'LineStyle','--')
    hold off
    xlabel('Orbit Number [-]')
    ylabel('Heat Load [kJ m^{-2}]')
    legend('Heat Load','Max. Heat Load','Location','Best')
    set(gca,'FontSize',15)
    grid on
    if saveFigures, saveas(F,'../../Report/figures/aero_heat_load','epsc'), end
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
    ylabel('Acceleration [m s^{-2}]')
    xlim([100,250])
    ylim([1e-5,1e-2])
    grid on
    legend('Noisy','Ideal','Reduced Interface','Location','NE')
    set(gca,'FontSize',15,'YScale','log')
    if saveFigures, saveas(F,'../../Report/figures/imu_acc_raia','epsc'), end
end

%% Plot Aerobraking Evolution

if saveFigures && strcmp(outputFolder,'SimulationOutputTransGuidOnly/')
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

if ~showFigures, close all, end
if saveFigures, close all, end