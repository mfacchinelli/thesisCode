%% Fit Atmospheric Data

%...Find pericenter height
% estimatedAltitude = sqrt( sum( CartesianEstimatedResults(:,1:3).^2, 2 ) ) - marsRadius / 1e3;
estimatedAltitude = altitude;
pericenter = min( estimatedAltitude );

%...Retireve density
% aerodynamicAcceleration = sqrt( sum( inertialMeasurements(:,1:3).^2, 2 ) );
% atmosphericDensity = 2 * 1000 / 37.5 ./ sum( ( CartesianEstimatedResults(1:end-1,4:6) * 1e3 ).^2, 2 ) / ...
%     sqrt( 1.87^2 + 0.013^2 ) .* aerodynamicAcceleration;
atmosphericDensity = accelerationResults;

%...Data below atmospheric interface
marsAtmosphericInterface = 175;
loc = ( estimatedAltitude <= marsAtmosphericInterface );
loc(1:find(estimatedAltitude==pericenter)+1:end) = false;
estimatedAltitude = estimatedAltitude( loc );
atmosphericDensity = atmosphericDensity( loc );

%...Least squares of data
designMatrix = [ ones( size( atmosphericDensity ) ), ( pericenter - estimatedAltitude ) * 1e3 ]; 
weightMatrix = diag(1./estimatedAltitude.^2);
updatedEstimate = ( designMatrix' * weightMatrix * designMatrix ) \ designMatrix' * weightMatrix * log( atmosphericDensity );
densityFunction = @(h) exp( updatedEstimate( 1 ) ) * exp( updatedEstimate( 2 ) * ( pericenter * 1e3 - h ) );

referenceDensity = @(h) 2.42386294453635e-08 * exp( ( pericenter * 1e3 - h ) / 6532.91178699257 );

%...MATLAB least squares
ft = fittype('rho0 + altitude / H',...
    'independent',{'altitude'},...
    'dependent',{'rho'},...
    'coefficients',{'rho0','H'});
lsq = fit(( pericenter - estimatedAltitude ) * 1e3,log(atmosphericDensity),ft,...
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
scatter( atmosphericDensity , estimatedAltitude(1:end) )
plot( referenceDensity( estimatedAltitude * 1e3 ), estimatedAltitude, 'LineWidth',1.25)
plot( densityFunction( estimatedAltitude * 1e3 ), estimatedAltitude, 'LineWidth',1.25,'LineStyle','--')
plot( matlabDensityFunction( ( pericenter - estimatedAltitude ) * 1e3 ), estimatedAltitude, 'LineWidth',1.25,'LineStyle','-.')
hold off
set(gca,'FontSize',15,'XScale','log')
grid on
legend('Measured','Reference','Own','MATLAB')

%% Fit Other Model

clc

%...Initial values
x = [log( 2.424e-08 ), 1.0 / 6533.0, -1.0, 0.0, 0.0]';
refAlt = ( estimatedAltitude - min(estimatedAltitude) ) * 1e3;
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
    [lambda,rho,dx']
end
newDensityFunction = @(h) exp(x(1)) * exp( x(3) * h * x(2) + x(4) * cos( 2*pi*h*x(2) ) + ...
    x(5) * sin( 2*pi*h*x(2) ) );

%...MATLAB least squares
ft = fittype('rho0 + K1 * altitude / H + K2 * cos( 2 * pi * altitude / H ) + K3 * sin( 2 * pi * altitude / H )',...
    'independent',{'altitude'},...
    'dependent',{'rho'},...
    'coefficients',{'rho0','H','K1','K2','K3'});
lsq = fit(refAlt,log(atmosphericDensity),ft,...
    'StartPoint',[log(2.42386294453635e-08),6532.91178699257,-1,0,0],'Robust','LAR','Algorithm','Trust-Region');
matlabDensityFunction = @(h) exp(lsq.rho0) * exp( lsq.K1 * h / lsq.H + lsq.K2 * cos( 2*pi*h/lsq.H ) + ...
    lsq.K3 * sin( 2*pi*h/lsq.H ) );

%...Show comparison table
table( [ 2.42386294453635e-08; exp( updatedEstimate(1) ); exp( x(1) ); exp( lsq.rho0 ) ], ...
    [ 6532.91178699257; 1 / updatedEstimate(2); 1/x(2); lsq.H ], ...
    [ - 1; NaN; x(3); lsq.K1 ],...
    [ 0; NaN; x(4); lsq.K2 ],...
    [ 0; NaN; x(5); lsq.K3 ],...
    'VariableNames', {'DensityAtReferenceAltitude', 'ScaleHeight','Kappa1','Kappa2','Kappa3'} )

%...Plot result
figure
hold on
scatter( atmosphericDensity, estimatedAltitude(1:end) )
plot( newDensityFunction( refAlt ), estimatedAltitude, 'LineWidth',1.25,'LineStyle','--')
plot( matlabDensityFunction( refAlt ), estimatedAltitude, 'LineWidth',1.25,'LineStyle','-.')
hold off
set(gca,'FontSize',15,'XScale','log')
grid on
legend('Measured','Own','MATLAB')

%% Run PTE

pteKeplerianEstimatedResults = KeplerianEstimatedResults(1:end-1,:);
pteKeplerianEstimatedResults(:,1) = pteKeplerianEstimatedResults(:,1) * 1e3;
pteKeplerianEstimatedResults(:,3:end) = deg2rad(pteKeplerianEstimatedResults(:,3:end));
[tp,DV,Da,DP] = PTE(onboardTime(1:end-1)*timeConversion,sqrt(sum(inertialMeasurements(:,1:3).^2,2)),...
    CartesianEstimatedResults(1:end-1,:)*1e3,pteKeplerianEstimatedResults,sqrt(sum(dependentVariables(1:end-1,4:6).^2,2)))
[tp,DV,Da,DP] = pte2(onboardTime(1:end-1)*timeConversion,sqrt(sum(inertialMeasurements(:,1:3).^2,2)),...
    CartesianEstimatedResults(1:end-1,:)*1e3,pteKeplerianEstimatedResults,sqrt(sum(dependentVariables(1:end-1,4:6).^2,2)))

%%

mu = 4.282e13;
R = 3.39e3;

ra = 40000e3 + R;

lowRp = 75e3 + R;
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
heatLoad = [];

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
    heatLoad = horzcat( heatLoad, trapz(time,heat) );
end

%%

figure
hold on
scatter3(heatRate/1e3,heatLoad/1e3,dynPress/1e3,75,(periapsisAltitudes-R)/1e3,'Filled')
xlimits = xlim;
ylimits = ylim;
zlimits = zlim;
plot3(xlimits(2)*ones(size(heatRate)),heatLoad/1e3,dynPress/1e3,'LineWidth',1.1,'LineStyle','--')
plot3(heatRate/1e3,ylimits(2)*ones(size(heatLoad)),dynPress/1e3,'LineWidth',1.1,'LineStyle','--')
plot3(heatRate/1e3,heatLoad/1e3,zlimits(2)*ones(size(dynPress)),'LineWidth',1.1,'LineStyle','--')
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
heatLoadOffset = ( heatLoad - maxHeatLoad ) / maxHeatLoad;
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

hRange = 100e3:5e3:40e6;
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
plot(hRange/1e3, tBeamTravel.*V)
xlabel('Altitude [km]')
ylabel('Distance Travelled [m]')
set(gca,'FontSize',15,'XScale','log')
grid on

%%

figure
distance = hRange(end)+marsRadius;
theta = 0.001:0.001:asin(marsRadius/distance);
pseudoAltitude = marsRadius ./ sin( theta ) .* sin( - theta + asin( distance/marsRadius * sin(theta) ) );
hold on
plot(rad2deg(theta),pseudoAltitude)
plot(rad2deg([theta(1),theta(end)]),[hRange(end),hRange(end)],'LineStyle','--')
hold off
grid on

figure
altitude = marsRadius * ( sin( pi - theta - asin( pseudoAltitude / marsRadius .* sin( theta ) ) ) ./ sin( theta ) - 1.0 ); 
hold on
plot(rad2deg(theta),altitude)
plot(rad2deg([theta(1),theta(end)]),[hRange(end),hRange(end)],'LineStyle','--')
hold off
grid on