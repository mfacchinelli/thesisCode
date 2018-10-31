fclose('all'); clear all; close all force; profile off; clc; format long g; rng default;
addpath functions

%...Figure settings
showFigure = true;
saveFigure = false;
[figSizeLarge,figSizeMedium,figSizeSmall] = saveFigureSettings(saveFigure);
figSizeWide = [150,150,1250,475];

%...Constants
marsRadius = 3389526.666666667;
marsGravitationalParameter = 42828375815756.1;
marsReducedAtmosphericInterface = 150;

minDynPress = 0.19;
maxHeatRate = 2800 / 2;
maxHeatLoad = 500 / 2;

%...Settings
loadPTE = true;
loadAE = true;
loadME = true;
loadPS = true;
loadMC = true;

timeLabel = 'Time [d]';
KeplerianLabels = {'a [km]','e [-]','i [deg]','\omega [deg]','\Omega [deg]','\vartheta [deg]'};

output = '/Users/Michele/GitHub/tudat/tudatBundle/tudat/bin/unit_tests/TestingResults/';

%% Read Data

if loadAE
    %...Retrieve altitude
    filename = [output,'nsAltitude.dat'];
    fileID = fopen(filename,'r');
    estimatedAltitude = textscan(fileID,'%f','Delimiter',',','CollectOutput',true);
    estimatedAltitude = estimatedAltitude{1}/1e3;
    fclose(fileID);
    
    %...Retrieve acceleration
    filename = [output,'nsAcceleration.dat'];
    fileID = fopen(filename,'r');
    measuredAcceleration = textscan(fileID,'%f','Delimiter',',','CollectOutput',true);
    measuredAcceleration = measuredAcceleration{1};
    fclose(fileID);
    
    %...Retrieve density
    filename = [output,'nsDensity.dat'];
    fileID = fopen(filename,'r');
    estimatedDensity = textscan(fileID,'%f','Delimiter',',','CollectOutput',true);
    estimatedDensity = estimatedDensity{1};
    fclose(fileID);
    
    %...Retrieve actual Keplerian elements
    filename = [output,'nsKepler_act.dat'];
    fileID = fopen(filename,'r');
    orbitData = textscan(fileID,repmat('%f ',[1,7]),'Delimiter',',','CollectOutput',true);
    actTime = orbitData{1}(:,1); actTime = actTime-actTime(1);
    KeplerianActualResults = orbitData{1}(:,2:end);
    fclose(fileID);
    
    KeplerianActualResultsPlot = KeplerianActualResults;
    KeplerianActualResultsPlot(:,1) = KeplerianActualResults(:,1) / 1e3;
    KeplerianActualResultsPlot(:,3:end) = rad2deg(KeplerianActualResults(:,3:end));
end

if loadPTE    
    %...Retrieve estimated Keplerian elements
    filename = [output,'nsKepler_est.dat'];
    fileID = fopen(filename,'r');
    orbitData = textscan(fileID,repmat('%f ',[1,7]),'Delimiter',',','CollectOutput',true);
    pteTime = orbitData{1}(:,1); pteTime = pteTime-pteTime(1);
    pteKeplerianEstimatedResults = orbitData{1}(:,2:end);
    fclose(fileID);
    
    pteKeplerianEstimatedResultsPlot = pteKeplerianEstimatedResults;
    pteKeplerianEstimatedResultsPlot(:,1) = pteKeplerianEstimatedResults(:,1) / 1e3;
    pteKeplerianEstimatedResultsPlot(:,3:end) = rad2deg(pteKeplerianEstimatedResults(:,3:end));
    pteKeplerianEstimatedResultsPlot(:,end) = wrapTo360(pteKeplerianEstimatedResultsPlot(:,end));
    
    %...Retrieve acceleration
    filename = [output,'nsAero.dat'];
    fileID = fopen(filename,'r');
    pteAcceleration = textscan(fileID,'%f','Delimiter',',','CollectOutput',true);
    pteAcceleration = pteAcceleration{1};
    fclose(fileID);
end

if loadME
    %...Retrieve heating conditions
    [meIndependentVariables,mePeakDynPress] = readNdData([output,'gsMePeakDynPress.dat']);
    [~,mePeakHeatRate] = readNdData([output,'gsMePeakHeatRate.dat']);
    [~,meHeatLoad] = readNdData([output,'gsMeHeatLoad.dat']);
    
    %...Retrieve orbit conditions
    [~,meNominalManeuver] = readNdData([output,'gsMeNominalManeuver.dat']);
    [~,meTargetPeriapsis] = readNdData([output,'gsMeTargetPeriapsis.dat']);
    [~,meActualPeriapsis] = readNdData([output,'gsMeActualPeriapsis.dat']);
end

if loadPS
    %...Retrieve heating conditions
    [psIndependentVariables,psPeakDynPress] = readNdData([output,'gsPsPeakDynPress.dat']);
    [~,psPeakHeatRate] = readNdData([output,'gsPsPeakHeatRate.dat']);
    [~,psHeatLoad] = readNdData([output,'gsPsHeatLoad.dat']);
    
    %...Retrieve orbit conditions
    [~,psNominalTarget] = readNdData([output,'gsPsNominalTarget.dat']);
    [~,psAdjustedTarget] = readNdData([output,'gsPsAdjustedTarget.dat']);
    [~,psActualPeriapsis] = readNdData([output,'gsPsActualPeriapsis.dat']);
end

if loadMC
    %...Retrieve heating conditions
    [mcIndependentVariables,mcPeakDynPress] = readNdData([output,'gsMcPeakDynPress.dat']);
    [~,mcPeakHeatRate] = readNdData([output,'gsMcPeakHeatRate.dat']);
    [~,mcHeatLoad] = readNdData([output,'gsMcHeatLoad.dat']);
    
    %...Retrieve orbit conditions
    [~,mcNominalTarget] = readNdData([output,'gsMcNominalTarget.dat']);
    [~,mcAdjustedTarget] = readNdData([output,'gsMcAdjustedTarget.dat']);
    [~,mcXError] = readNdData([output,'gsMcXError.dat']);
    [~,mcYError] = readNdData([output,'gsMcYError.dat']);
    [~,mcZError] = readNdData([output,'gsMcZError.dat']);
    [~,mcVxError] = readNdData([output,'gsMcVxError.dat']);
    [~,mcVyError] = readNdData([output,'gsMcVyError.dat']);
    [~,mcVzError] = readNdData([output,'gsMcVzError.dat']);
end
    
if loadPTE && loadAE
    %...Plot Keplerian translational motion
    figure
    for i = 1:size(KeplerianActualResults,2)
        subplot(2,3,i)
        hold on
        plot(actTime,KeplerianActualResultsPlot(:,i),'LineWidth',1.25)
        plot(pteTime,pteKeplerianEstimatedResultsPlot(:,i),'LineWidth',1.25)
        hold off
        xlabel(timeLabel)
        ylabel(KeplerianLabels{i})
        set(gca,'FontSize',15)
        grid on
    end
    subplotLegend({'Actual','Estimated'})
end

%% Run PTE

if loadPTE
    clc
    
    %...Run Periapse Time Estimator function
    [tp,dtheta,DV,Da,De] = PTE(pteTime,pteKeplerianEstimatedResults,pteAcceleration)
    
    keplerAltitude = pteKeplerianEstimatedResults(:,1) .* ( 1 - pteKeplerianEstimatedResults(:,2).^2 ) ./ ...
        ( 1 + pteKeplerianEstimatedResults(:,2) .* cos( pteKeplerianEstimatedResults(:,6) ) );
    keplerAltitude = ( keplerAltitude - marsRadius ) / 1e3;
    keplerAltitude = keplerAltitude( keplerAltitude <= marsReducedAtmosphericInterface );
    
    correctedAltitude = pteKeplerianEstimatedResults(:,1) .* ( 1-pteKeplerianEstimatedResults(:,2).^2 ) ./ ...
        ( 1 + pteKeplerianEstimatedResults(:,2) .* cos( pteKeplerianEstimatedResults(:,6)+deg2rad(dtheta) ) ) - marsRadius;
    correctedAltitude = correctedAltitude(correctedAltitude<150e3)/1e3;
    
    figure
    hold on
    plot(keplerAltitude)
    plot(correctedAltitude)
    if loadAE, plot(estimatedAltitude), end
    hold off
    grid on
end

%% Fit Atmospheric Data

clc

if loadAE
    %...Find pericenter height
    pericenter = min( estimatedAltitude );
    
    %...Data below atmospheric interface
    marsAtmosphericInterface = 150;
    loc = ( estimatedAltitude <= marsAtmosphericInterface );
    loc(1:find(estimatedAltitude==pericenter)+1:end) = false;
    altitude = estimatedAltitude( loc );
    atmosphericDensity = estimatedDensity( loc );
    
    %...Least squares of data
    designMatrix = [ ones( size( atmosphericDensity ) ), ( pericenter - altitude ) * 1e3 ];
    weightMatrix = diag(1./altitude.^2);
    updatedEstimate = ( designMatrix' * weightMatrix * designMatrix ) \ designMatrix' * weightMatrix * log( atmosphericDensity );
    densityFunction = @(h) exp( updatedEstimate( 1 ) ) * exp( updatedEstimate( 2 ) * ( pericenter * 1e3 - h ) );
    
    referenceDensity = @(h) 2.42386294453635e-08 * exp( ( pericenter * 1e3 - h ) / 6532.91178699257 );
    
    %...MATLAB least squares
    ft = fittype('rho0 + altitude / H',...
        'independent',{'altitude'},...
        'dependent',{'rho'},...
        'coefficients',{'rho0','H'});
    [lsq,G,O] = fit(( pericenter - altitude ) * 1e3,log(atmosphericDensity),ft,...
        'StartPoint',[log(2.42386294453635e-08),6532.91178699257],'Robust','LAR','Algorithm','Trust-Region');
    matlabDensityFunction = @(h) exp(lsq.rho0) * exp( h / lsq.H );
    
    %...Show comparison table
    table( [ 2.42386294453635e-08; exp( updatedEstimate( 1 ) ); exp( lsq.rho0 ) ], ...
        [ 6532.91178699257; 1 / updatedEstimate( 2 ); lsq.H ], ...
        'VariableNames', {'DensityAtReferenceAltitude', 'ScaleHeight'} )
    
    %...Plot result
    % F = figure('rend','painters','pos',figSizeLarge);
    figure
    hold on
    scatter( atmosphericDensity , altitude(1:end) )
    plot( densityFunction( altitude * 1e3 ), altitude, 'LineWidth',1.25,'LineStyle','--')
    plot( matlabDensityFunction( ( pericenter - altitude ) * 1e3 ), altitude, 'LineWidth',1.25,'LineStyle','-.')
    % plot( referenceDensity( altitude * 1e3 ), altitude, 'LineWidth',1.25)
    hold off
    set(gca,'FontSize',15,'XScale','log')
    grid on
    legend('Measured','Own','MATLAB')%,'Reference')
end

%% Fit Other Model

if loadAE
    clc
    
    %...Initial values
    x = [log( 2.424e-08 ), 1.0 / 6533.0, -1.0, 0.0, 0.0]';
    refAlt = ( altitude - min(altitude) ) * 1e3;
    twoPiRefAlt = 2*pi*refAlt;
    
    tau = 1e-6;
    lambda = 0;
    nu = 2;
    lsqFunc = @(x) x(1) + ( x(3) * x(2) * refAlt + x(4) * cos( twoPiRefAlt * x(2) ) + x(5) * sin( twoPiRefAlt * x(2) ) );
    
    %...Iterate over values
    for i = 1:10
        estDens = lsqFunc(x);
        A = [ ones(size( refAlt )), ...
            x(3) * refAlt - x(4) * sin( twoPiRefAlt * x(2) ) .* twoPiRefAlt + x(5) * cos( twoPiRefAlt * x(2) ) .* twoPiRefAlt, ...
            x(2) * refAlt, cos( twoPiRefAlt * x(2) ), sin( twoPiRefAlt * x(2) ) ];
        
        if i == 1
            lambda = tau * max(diag(A'*A));
        end
        
        offset = ( log( atmosphericDensity ) - estDens );
        dx = (A'*A + lambda * eye(5,5)) \ A' * offset;
        
        rho = ( lsqFunc(x)'*lsqFunc(x) - lsqFunc(x+dx)'*lsqFunc(x+dx) ) / ( dx' * ( lambda * dx - A' * offset ) );
        if rho > 0
            lambda = lambda * max( 1/3, 1-(2*rho-1)^3 );
            nu = 2;
            x = x + dx;
        else
            lambda = lambda * nu;
            nu = nu * 2;
        end
        %     [lambda,rho,dx']
    end
    newDensityFunction = @(h) exp(x(1)) * exp( x(3) * h * x(2) + x(4) * cos( 2*pi*h*x(2) ) + ...
        x(5) * sin( 2*pi*h*x(2) ) );
    
    %...MATLAB least squares
    ft = fittype('rho0 + K1 * altitude / H + K2 * cos( 2 * pi * altitude / H ) + K3 * sin( 2 * pi * altitude / H )',...
        'independent',{'altitude'},...
        'dependent',{'rho'},...
        'coefficients',{'rho0','H','K1','K2','K3'});
    [lsq,G,O] = fit(refAlt,log(atmosphericDensity),ft,...
        'StartPoint',[log(2.42386294453635e-08),6532.91178699257,-1,0,0],'Robust','LAR','Algorithm','Trust-Region');
    matlabDensityFunction = @(h) exp(lsq.rho0) * exp( lsq.K1 * h / lsq.H + lsq.K2 * cos( 2*pi*h/lsq.H ) + ...
        lsq.K3 * sin( 2*pi*h/lsq.H ) );
    
    %...Show comparison table
    table( [ 2.42386294453635e-08; exp( updatedEstimate(1) ); exp( x(1) ); exp( lsq.rho0 ) ], ...
        [ 6532.91178699257; 1 / updatedEstimate(2); 1/x(2); lsq.H ], ...
        [ - 1; NaN; x(3); lsq.K1 ],...
        [ 0; NaN; x(4); lsq.K2 ],...
        [ 0; NaN; x(5); lsq.K3 ],...
        'VariableNames',{'DensityAtReferenceAltitude','ScaleHeight','Kappa1','Kappa2','Kappa3'} )
    
    %...C++ results
    % cppFit = [132831 2.37122e-09     3991.26   -0.478632   0.0050322  0.00490655];
    cppFit = [110122 5.62421e-08     6254.71    -1.13478    0.158959    0.175122];
    cppDensityFunction = @(h) cppFit(2) * exp( cppFit(4) * ( h*1e3 - cppFit(1) ) / cppFit(3) + ...
        cppFit(5) * cos( 2*pi*( ( h*1e3 - cppFit(1) ) / cppFit(3) ) ) + ...
        cppFit(6) * sin( 2*pi*( ( h*1e3 - cppFit(1) ) / cppFit(3) ) ) );
    
    %...Plot result
    figure
    hold on
    scatter( atmosphericDensity, altitude )
    plot( newDensityFunction( refAlt ), altitude, 'LineWidth',1.25,'LineStyle','--')
    plot( matlabDensityFunction( refAlt ), altitude, 'LineWidth',1.25,'LineStyle',':')
    % plot( cppDensityFunction( altitude ), altitude, 'LineWidth',1.25,'LineStyle','-.')
    hold off
    set(gca,'FontSize',15,'XScale','log')
    grid on
    legend('Measured','Own','MATLAB')%,'C++')
end

%% Plot Periapsis Sensitivity

if loadPS
    styles = {'-o','-d','-s','-v','-p','-h','-*','-x','-^','-o','-d'};
    altitudeOffset = 1.5 * psIndependentVariables{3};
    
    %...Plot heat rate
    F = figure('rend','painters','pos',figSizeSmall);
    hold on
    for i = 1:length(psIndependentVariables{1})
        plot(altitudeOffset,squeeze(psPeakHeatRate(i,2,:))-maxHeatRate,styles{i},'LineWidth',1.25,'MarkerSize',10)
    end
    plot(xlim,[0,0],'LineWidth',1.25,'LineStyle','--')
    hold off
    xlabel('Altitude Offset [km]')
    ylabel('Offset in Heat Rate [W m^{-2}]')
    set(gca,'FontSize',15)
    grid on
    legend('Walk-in','Main','Walk-out','Threshold','Location','NE')
    if saveFigure, saveas(F,'../../../Report/figures/sens_heat_rate','epsc'), end
    
    %...Plot heat load
    F = figure('rend','painters','pos',figSizeSmall);
    hold on
    for i = 1:length(psIndependentVariables{1})
        plot(altitudeOffset,squeeze(psHeatLoad(i,2,:))/1e3-maxHeatLoad,styles{i},'LineWidth',1.25,'MarkerSize',10)
    end
    plot(xlim,[0,0],'LineWidth',1.25,'LineStyle','--')
    hold off
    xlabel('Altitude Offset [km]')
    ylabel('Offset in Heat Load [kJ m^{-2}]')
    set(gca,'FontSize',15)
    grid on
    legend('Walk-in','Main','Walk-out','Threshold','Location','NE')
    if saveFigure, saveas(F,'../../../Report/figures/sens_heat_load','epsc'), end
end

%% Plot Maneuver Sensitivity

if loadME
    styles = {'-o','-d','-s','-v','-p','-h','-*','-x','-^','-o','-d'};
    magnitudeOffset = 3 * meIndependentVariables{3} * 1e-1;
    
    %...Plot heat rate
    F = figure('rend','painters','pos',figSizeSmall);
    hold on
    for i = 1:length(meIndependentVariables{1})
        plot(magnitudeOffset,squeeze(mePeakHeatRate(i,2,:))-maxHeatRate,styles{i},'LineWidth',1.25,'MarkerSize',10)
    end
    plot(xlim,[0,0],'LineWidth',1.25,'LineStyle','--')
    hold off
    xlabel('Maneuver Magnitude Offset [m s^{-1}]')
    ylabel('Offset in Heat Rate [W m^{-2}]')
    set(gca,'FontSize',15)
    grid on
    legend('Walk-in','Main','Walk-out','Threshold','Location','NE')
    if saveFigure, saveas(F,'../../../Report/figures/sens_heat_rate','epsc'), end
    
    %...Plot heat load
    F = figure('rend','painters','pos',figSizeSmall);
    hold on
    for i = 1:length(meIndependentVariables{1})
        plot(magnitudeOffset,squeeze(meHeatLoad(i,2,:))/1e3-maxHeatLoad,styles{i},'LineWidth',1.25,'MarkerSize',10)
    end
    plot(xlim,[0,0],'LineWidth',1.25,'LineStyle','--')
    hold off
    xlabel('Maneuver Magnitude Offset [m s^{-1}]')
    ylabel('Offset in Heat Load [kJ m^{-2}]')
    set(gca,'FontSize',15)
    grid on
    legend('Walk-in','Main','Walk-out','Threshold','Location','NE')
    if saveFigure, saveas(F,'../../../Report/figures/sens_heat_load','epsc'), end
end

if loadPS
    %...Plot maneuver targeting
    F = figure('rend','painters','pos',figSizeSmall);
    hold on
    for i = 1:length(psIndependentVariables{1})
        plot(psIndependentVariables{3},(squeeze(psActualPeriapsis(i,2,:))-squeeze(psAdjustedTarget(i,2,:)))/1e3,...
            styles{i},'LineWidth',1.25,'MarkerSize',10)
    end
    hold off
    xlabel('Altitude Offset [km]')
    ylabel('Offset in Target [km]')
    set(gca,'FontSize',15)
    grid on
    legend('Walk-in','Main','Walk-out','Location','Best')
    if saveFigure, saveas(F,'../../../Report/figures/sens_maneuver','epsc'), end
end

%% Plot Monte Carlo

if loadMC
    styles = {'o','s','^'};
    colors = {[0,0.447,0.741],[0.85,0.325,0.098],[0.929,0.694,0.125]};
    
    %...Plot heat rate
    F = figure('rend','painters','pos',figSizeSmall);
    hold on
    for i = 1:length(mcIndependentVariables{1})
        for j = 1:length(mcIndependentVariables{2})
            scatter(mcPeakHeatRate(i,j)-maxHeatRate,mcHeatLoad(i,j)/1e3-maxHeatLoad,50,colors{j},styles{j})
        end
    end
    hold off
    xlabel('Offset in Heat Rate [W m^{-2}]')
    ylabel('Offset in Heat Load [kJ m^{-2}]')
    set(gca,'FontSize',15)
    grid on
    legend('Init. Cond. 1','Init. Cond. 2','Init. Cond. 3','Location','SE')
    if saveFigure, saveas(F,'../../../Report/figures/mc_heat','epsc'), end
end

%%

mu = 4.282e13;
R = 3.39e3;

ra = 40000e3 + R;

lowRp = 100e3 + R;
highRp = 175e3 + R;
dtheta = deg2rad(0.1);

maxHeatRate = 2800e3;
maxHeatLoad = 500e3;
minDynPress = 0.19;

rho0 = 7.00137e-10;
h0 = 148162;
H = 9160.63;
dens = @(r) rho0 * exp( - ( r - R - h0 ) / H );
vel = @(r,a) sqrt( mu * ( 2 ./ r - 1 / a ) );

%%

periapsisAltitudes = lowRp:1e3:highRp;
dynPress = [];
heatRate = [];
psHeatLoad = [];

for rp = periapsisAltitudes
    a = 0.5 * ( ra + rp );
    e = (ra - rp)/(ra + rp);
    p = a * ( 1 - e^2 );
    
    atmTheta = acos( ( p / (250e3+R) - 1 ) / e );
    
    theta = -atmTheta:dtheta:atmTheta;
    
    radius = p ./ ( 1 + e * cos( theta ) );
    velocity = vel( radius, a );
    
    dyn = 0.5 * dens( radius ) .* velocity.^2;
    heat = dyn .* velocity;
    dtime = abs( dtheta * radius.^2 / sqrt( mu * p ) .* theta );
    time = cumsum(dtime);
    
    dynPress = horzcat( dynPress, max(dyn) );
    heatRate = horzcat( heatRate, max(heat) );
    psHeatLoad = horzcat( psHeatLoad, trapz(time,heat) );
end

%%

figure
hold on
scatter3(heatRate/1e3,psHeatLoad/1e3,dynPress/1e3,75,(periapsisAltitudes-R)/1e3,'Filled')
xlimits = xlim;
ylimits = ylim;
zlimits = zlim;
plot3(xlimits(2)*ones(size(heatRate)),psHeatLoad/1e3,dynPress/1e3,'LineWidth',1.1,'LineStyle','--')
plot3(heatRate/1e3,ylimits(2)*ones(size(psHeatLoad)),dynPress/1e3,'LineWidth',1.1,'LineStyle','--')
plot3(heatRate/1e3,psHeatLoad/1e3,zlimits(2)*ones(size(dynPress)),'LineWidth',1.1,'LineStyle','--')
hold off
c = colorbar;
c.Label.String = 'Periapsis Altitudes [km]';
xlabel('Heat Rate [kW m^{-2}]')
ylabel('Heat Load [kJ m^{-2}]')
zlabel('Dynamic Pressure [kN m^{-2}]')
xlim(xlimits)
ylim(ylimits)
zlim(zlimits)
grid on
axis on
set(gca,'FontSize',15,'XScale','log','YScale','log','ZScale','log')

%%

heatRateOffset = ( heatRate - maxHeatRate ) / maxHeatRate;
heatLoadOffset = ( psHeatLoad - maxHeatLoad ) / maxHeatLoad;
dynPressOffset = minDynPress - dynPress;

sign = ones(size(heatRateOffset));
sign(heatRateOffset < 0 & heatLoadOffset < 0) = -1;

figure
hold on
plot( periapsisAltitudes/1e3, max(heatRateOffset,heatLoadOffset),'LineWidth',1.1)
% plot( periapsisAltitudes/1e3, heatRateOffset,'LineWidth',1.1)
% plot( periapsisAltitudes/1e3, heatLoadOffset,'LineWidth',1.1)
plot( periapsisAltitudes/1e3, dynPressOffset,'LineWidth',1.1)
plot( [lowRp,highRp]/1e3, [0,0],'LineWidth',1.1,'LineStyle','--')
hold off
xlabel('Periapsis Altitudes [km]')
ylabel('Corridor Condition Number [-]')
% legend('Lower Altitude Bound','Heat Rate','Heat Load','Higher Altitude Bound')
grid on
set(gca,'FontSize',15)

%%

figure;
yyaxis left
plot(rad2deg(theta),theta)
ylim([min(theta),max(theta)])
yyaxis right
plot(rad2deg(theta),time)
ylim([min(time),max(time)])
grid on
% xlim([0,360])

time(end)
2 * pi * sqrt( a^3 / mu ) / 3600

%% Altimeter

hRange = 100e3:5e3:1e8;
theta = 2*atan(marsRadius./(hRange+marsRadius));
figure;
semilogx(hRange/1e3,rad2deg(theta))
grid on

fullFOV = 30e3 / 400e3;
pixelFOV = 6 / 400e3; % from MRO CTX

figure
hold on
plot(hRange/1e3,marsRadius./tan((theta+2*pixelFOV)/2)-marsRadius-hRange)
plot(hRange/1e3,marsRadius./tan((theta-2*pixelFOV)/2)-marsRadius-hRange)
hold off
xlabel('Altitude [km]')
ylabel('Error in Altitude [m]')
set(gca,'FontSize',15,'XScale','log')
grid on

%%

V = sqrt(marsGravitationalParameter * (2./(marsRadius+hRange) - 2./(2*marsRadius+hRange(1)+hRange(end))));
tBeamTravel = 2*hRange / 299792458;

figure
plot(hRange/1e3,tBeamTravel.*V,'LineWidth',1.5)
xlabel('Altitude [km]')
ylabel('Distance Travelled [m]')
set(gca,'FontSize',15,'XScale','log')
grid on

%%

pseudoAltitudeFunction = @(distance,theta) distance .* cos(theta) - ...
    sqrt( marsRadius^2 - distance.^2 * sin( theta ).^2 );

figure
alt = 4.50065e+07;%hRange(end)
distance = alt+marsRadius;
theta = 0.001:0.001:asin(marsRadius/distance);
pseudoAltitude = pseudoAltitudeFunction(distance,theta);
%marsRadius ./ sin( theta ) .* sin( - theta + asin( distance/marsRadius * sin(theta) ) );
hold on
plot(rad2deg(theta),pseudoAltitude,'LineWidth',1.5)
plot(rad2deg([theta(1),theta(end)]),[alt,alt],'LineWidth',1.5,'LineStyle','--')
plot([3.5,3.5],ylim,'LineWidth',1.5,'LineStyle',':')
hold off
xlabel('Pointing Angle [deg]')
ylabel('Altitude [km]')
set(gca,'FontSize',15)
grid on

% figure
% altitude = marsRadius * ( sin( pi - theta - asin( pseudoAltitude / marsRadius .* sin( theta ) ) ) ./ sin( theta ) - 1.0 ); 
% hold on
% plot(rad2deg(theta),altitude)
% plot(rad2deg([theta(1),theta(end)]),[hRange(end),hRange(end)],'LineStyle','--')
% hold off
% grid on

%%

figure
hold on
plot(hRange/1e3,asind(marsRadius./(hRange+marsRadius)),'LineWidth',1.5)
plot([5e3,5e3],ylim,'LineWidth',1.5,'LineStyle','--')
hold off
xlabel('Altitude [km]')
ylabel('Maximum Pointing Angle [deg]')
set(gca,'FontSize',15,'XScale','log')
grid on

%%

altitudeInformation = [90 92.7385 0.970471
    130 132.758 0.979228
    110 112.78 0.975352
    100 102.776 0.972989
    105 107.768 0.974316
    107.5 110.328 0.974365
    108.75 111.519 0.975171
    108.125 110.917 0.974828
    107.812 110.605 0.974757
    100 102.776 0.972989
    130 132.758 0.979228
    115 117.89 0.975489
    107.5 110.328 0.974365
    111.25 114.119 0.974858
    113.125 115.891 0.976131
    114.062 116.826 0.976346
    114.531 117.4 0.975566];
% altitudeInformation = [90 92.4593 0.973402
%     130 132.491 0.981201
%     110 112.526 0.977551
%     100 102.506 0.975554
%     105 107.501 0.976733
%     107.5 110.059 0.976747
%     108.75 111.251 0.977522
%     108.125 110.649 0.977191
%     107.812 110.352 0.976986
%     100 102.506 0.975554
%     130 132.491 0.981201
%     115 117.664 0.97736
%     107.5 110.059 0.976747
%     111.25 113.853 0.977135
%     113.125 115.625 0.97838
%     114.062 116.57 0.978492
%     114.531 117.139 0.977737];

altitudePlot = linspace(min(altitudeInformation(:,1)),max(altitudeInformation(:,1)),15);

designMatrix = ones(size(altitudeInformation,1),2);
designMatrix(:,2) = altitudeInformation(:,1);

x = (designMatrix'*designMatrix)\designMatrix'*altitudeInformation(:,3);

F = figure('rend','painters','pos',figSizeSmall);
hold on
scatter(altitudeInformation(:,1),altitudeInformation(:,3),50)
plot(altitudePlot,x(1) + x(2) * altitudePlot,'LineWidth',1.5,'LineStyle','--')
hold off
xlabel('Altitude Guess [km]')
ylabel('Ratio [-]')
set(gca,'FontSize',15)
legend('Ratio','Linear Fit','Location','Best')
grid on
if saveFigure, saveas(F,'../../Report/figures/guid_altitude_ratio','epsc'), end