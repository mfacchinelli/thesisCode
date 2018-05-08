fclose('all'); clear all; close all force; profile off; clc; format long g; rng default;
addpath tests data functions

%% State Derivative USM7

clc;

%...Constants
mu = 4.28284e+13;

%...USM state
USMState = [2549.57,79.2391,-2160.97,-0.670883,-0.275833,-0.0126141,-0.688239]';
C = USMState(1);
R = USMState(2:3);
eps = USMState(4:6);
eta = USMState(7);

%...Acceleration
acc = zeros(3,1);
acc = [5.01045e-07;2.26855e-06;-3.01176e-06];

%...Compute right ascension of latitude
lambda = atan2(2*eps(3)*eta/(eps(3)^2+eta^2),(eta^2-eps(3)^2)/(eps(3)^2+eta^2));
sinLambda = sin(lambda)
cosLambda = cos(lambda)

%...Compute gamma
gamma = (eps(1)*eps(3) - eps(2)*eta)/(eps(3)^2+eta^2)

%...Compute rotational velocity
omega = zeros(3,1);
omega(1) = acc(3)/(C - R(1)*sinLambda + R(2)*cosLambda);
omega(3) = C/mu * (C - R(1)*sinLambda + R(2)*cosLambda)^2;
omega

%...P parameter
p = 1/(C - R(1)*sinLambda + R(2)*cosLambda) * [C;R(2);R(1)]

%...Hodograph matrix
hodograph = [0, -p(1), 0;
    cosLambda, -(1+p(1)) * sinLambda, - gamma * p(2);
    sinLambda, (1+p(1)) * cosLambda, gamma * p(3)]

%...Augmented matrix
augmented = [0,omega(3),0,omega(1);
    -omega(3),0,omega(1),0;
    0,-omega(1),0,omega(3);
    -omega(1),0,-omega(3),0]

%...Derivative
[hodograph*acc;0.5*augmented*[eps;eta]]

%% State Derivative USM6

clc;
skew = @(x) [0,-x(3),x(2);x(3),0,-x(1);-x(2),x(1),0];

%...Constants
mu = 4.28284e+13;

%...USM state
USMState = [2494.82    79.0582   -2156.03   0.397386   0.163385 0.00747172,true]';
C = USMState(1);
R = USMState(2:3);
mrp = USMState(4:6);
flag = USMState(7);
sigma = norm(mrp)

%...Acceleration
acc = zeros(3,1);
acc = [7.75075e-07 2.99033e-06  -3.964e-06];

%...Compute right ascension of latitude
den = (1-norm(mrp)^2)^2+4*mrp(3)^2
lambda = atan2(4*mrp(3)*(1-norm(mrp)^2)/den,((1-norm(mrp)^2)^2-4*mrp(3)^2)/den);
sinLambda = sin(lambda)
cosLambda = cos(lambda)

%...Compute gamma
gamma = 2*(mrp(2)^3+mrp(1)^2*mrp(2)+mrp(2)*mrp(3)^2+2*mrp(1)*mrp(3)-mrp(2))/...
    (mrp(1)^4+mrp(2)^4+mrp(3)^4+2*(mrp(1)^2*mrp(2)^2 + mrp(1)^2*mrp(3)^2 + mrp(2)^2*mrp(3)^2 - mrp(1)^2-mrp(2)^2+mrp(3)^2)+1)
2*(mrp(2)*(sigma^2-1)+2*mrp(1)*mrp(3))/...
    (mrp(1)^4+mrp(2)^4+mrp(3)^4+2*(mrp(1)^2* ( mrp(2)^2 + mrp(3)^2 - 1 ) + ( mrp(2)^2 + 1 ) * ( mrp(3)^2 - 1 ) + 1 ) + 1 )

%...Compute rotational velocity
omega = zeros(3,1);
omega(1) = acc(3)/(C - R(1)*sinLambda + R(2)*cosLambda);
omega(3) = C/mu * (C - R(1)*sinLambda + R(2)*cosLambda)^2;
omega

%...P parameter
p = 1/(C - R(1)*sinLambda + R(2)*cosLambda) * [C;R(2);R(1)]

%...Hodograph matrix
hodograph = [0, -p(1), 0;
    cosLambda, -(1+p(1)) * sinLambda, - gamma * p(2);
    sinLambda, (1+p(1)) * cosLambda, gamma * p(3)]

%...Augmented matrix
mrpDot = 1/4 * ( (1-sigma^2) * eye(3) + 2 * skew(mrp) + 2 * mrp * mrp' ) * omega

%% State Derivative USMEM

clc;
skew = @(x) [0,-x(3),x(2);x(3),0,-x(1);-x(2),x(1),0];

%...Constants
mu = 4.28284e+13;

%...USM state
USMState = [2541.35,78.8923,-2151.51,-0.750645,-0.308627,-0.0141137]';
C = USMState(1);
R = USMState(2:3);
expVec = USMState(4:6)

%...Acceleration
acc = zeros(3,1);
acc = [1.10672e-06,6.52374e-07,2.06722e-06]';

%...Quaternions
expMag = norm(expVec)
eta = cos(expMag/2);
eps = expVec/expMag * sin(expMag/2);

%...Compute right ascension of latitude
lambda = atan2(2*eps(3)*eta/(eps(3)^2+eta^2),(eta^2-eps(3)^2)/(eps(3)^2+eta^2));
sinLambda = sin(lambda)
cosLambda = cos(lambda)

%...Compute gamma
gamma = (eps(1)*eps(3) - eps(2)*eta)/(eps(3)^2+eta^2)

%...Compute rotational velocity
omega = zeros(3,1);
omega(1) = acc(3)/(C - R(1)*sinLambda + R(2)*cosLambda);
omega(3) = C/mu * (C - R(1)*sinLambda + R(2)*cosLambda)^2;
omega

%...P parameter
p = 1/(C - R(1)*sinLambda + R(2)*cosLambda) * [C;R(2);R(1)]

%...Hodograph matrix
hodograph = [0, -p(1), 0;
    cosLambda, -(1+p(1)) * sinLambda, - gamma * p(2);
    sinLambda, (1+p(1)) * cosLambda, gamma * p(3)]

%...Augmented matrix
expDot = (eye(3) + skew(expVec)/2 + 1/expMag^2 * ( 1-expMag/2 * cot(expMag/2) ) * (skew(expVec))^2 ) * omega

%...Derivative
[hodograph*acc;expDot]

%% Read Kepler Orbit

test = false;
kepler_orbit = false;

%...Mars data
R = 3.390e3;

%...Get orbital data
fileID = fopen('/Users/Michele/GitHub/tudat/tudatBundle/tudatApplications/Test/SimulationOutput/test_orbit.dat');
kepler = textscan(fileID,repmat('%f ',[1,7]),'CollectOutput',true,'Delimiter',',');
time = (kepler{1}(:,1)-kepler{1}(1))/3600; kepler = kepler{1}(:,2:end);
kepler(:,1) = kepler(:,1)/1e3; kepler(:,3:end) = rad2deg(kepler(:,3:end));
if kepler_orbit, kepler(:,1:5) = abs( kepler(:,1:5) - kepler(1,1:5) ); end % get difference
fclose(fileID);

%...Get trajecotry data
fileID = fopen('/Users/Michele/GitHub/tudat/tudatBundle/tudatApplications/Test/SimulationOutput/test_trajectory.dat');
cartesian = textscan(fileID,repmat('%f ',[1,8]),'CollectOutput',true,'Delimiter',',');
cartesian = cartesian{1}(:,2:7); cartesian = cartesian/1e3;
fclose(fileID);

%...Interpolate data
timeConstantStepSize = 0:1/3600:48;
kepler = interp1(time,kepler,timeConstantStepSize,'spline');
cartesian = interp1(time,cartesian,timeConstantStepSize,'spline');

%...Plot Kepler elements
figure;
labels = {'Semi-major Axis [km]','Eccentricity [-]','Inclination [deg]',...
    'Right Ascension of Ascending Node [deg]','Argument of Perigee [deg]','True Anomaly [deg]'};
for i = 1:6
    subplot(2,3,i)
    hold on
    plot(timeConstantStepSize,kepler(:,i),'LineWidth',1.1)
    if i == 1 && ( ~test || ~kepler_orbit )
        plot([timeConstantStepSize(1),timeConstantStepSize(end)],R*ones(1,2),'LineStyle','--')
    end
    hold off
    xlabel('Time [h]')
    ylabel(labels{i})
    grid on
    if i ~= 6 && test && kepler_orbit, set(gca,'FontSize',15,'YScale','log'), else, set(gca,'FontSize',15), end
end

%...Plot Cartesian elements
[x,y,z] = sphere;
figure;
hold on
plot3(cartesian(:,1),cartesian(:,2),cartesian(:,3),'LineWidth',1.1)
if test, scatter3(cartesian(:,1),cartesian(:,2),cartesian(:,3),30), end
surf(R*x,R*y,R*z)
hold off
xlabel('x [km]')
ylabel('y [km]')
zlabel('z [km]')
grid on
axis equal
set(gca,'FontSize',15)

%% Compare Results of Cowell and USM7

% save('data/XXX_solution','cartesian','kepler')

%...Load results
time = 0:86400*2;
cowell = load('cowell_solution');
encke = load('encke_solution');
usm7 = load('usm7_solution');
usm6 = load('usm6_solution');
usmem = load('usmem_solution');
gauss_kepler = load('gauss_kepler_solution');

reference = usmem;

%...Difference in Keplerian elements
figure;
labels = {'\Delta a [km]','\Delta e [-]','\Delta i [deg]',...
    '\Delta\Omega [deg]','\Delta\omega [deg]','\Delta\vartheta [deg]'};
for i = 1:6
    subplot(2,3,i)
    hold on
    plot(time,abs(cowell.kepler(:,i)-reference.kepler(:,i)),'LineWidth',1.1)
    plot(time,abs(encke.kepler(:,i)-reference.kepler(:,i)),'LineWidth',1.1)
    plot(time,abs(usm7.kepler(:,i)-reference.kepler(:,i)),'LineWidth',1.1)
    plot(time,abs(usm6.kepler(:,i)-reference.kepler(:,i)),'LineWidth',1.1)
    plot(time,abs(usmem.kepler(:,i)-reference.kepler(:,i)),'LineWidth',1.1)
    plot(time,abs(gauss_kepler.kepler(:,i)-reference.kepler(:,i)),'LineWidth',1.1)
%     plot(time,abs(kepler(:,i)-reference.kepler(:,i)),'LineWidth',1.1)
    hold off
    xlabel('Time [h]')
    ylabel(labels{i})
    grid on
    set(gca,'FontSize',15,'YScale','log')
end
SubplotLegend({'Cowell','Encke','USM7','USM6','USMEM','Gauss'}) %,'Current'

%...Difference in Cartesian elements
figure;
labels = {'\Delta x [km]','\Delta y [km]','\Delta z [km]',...
    '\Delta v_x [km/s]','\Delta v_y [km/s]','\Delta v_z [km/s]'};
for i = 1:6
    subplot(2,3,i)
    hold on
    plot(time,abs(cowell.cartesian(:,i)-reference.cartesian(:,i)),'LineWidth',1.1)
    plot(time,abs(encke.cartesian(:,i)-reference.cartesian(:,i)),'LineWidth',1.1)
    plot(time,abs(usm7.cartesian(:,i)-reference.cartesian(:,i)),'LineWidth',1.1)
    plot(time,abs(usm6.cartesian(:,i)-reference.cartesian(:,i)),'LineWidth',1.1)
    plot(time,abs(usmem.cartesian(:,i)-reference.cartesian(:,i)),'LineWidth',1.1)
    plot(time,abs(gauss_kepler.cartesian(:,i)-reference.cartesian(:,i)),'LineWidth',1.1)
%     plot(time,abs(cartesian(:,i)-reference.cartesian(:,i)),'LineWidth',1.1)
    hold off
    xlabel('Time [h]')
    ylabel(labels{i})
    grid on
    set(gca,'FontSize',15,'YScale','log')
end
SubplotLegend({'Cowell','Encke','USM7','USM6','USMEM','Gauss'}) %,'Current'

%% Plot Together

figure;
labels = {'a [km]','e [-]','i [deg]',...
    '\Omega [deg]','\omega [deg]','\vartheta [deg]'};
for i = 1:6
    subplot(2,3,i)
    hold on
    plot(time,kepler(:,i),'LineWidth',1.1)
    plot(time,usm7.kepler(:,i),'LineWidth',1.1)
    hold off
    xlabel('Time [h]')
    ylabel(labels{i})
    grid on
    set(gca,'FontSize',15)
end
SubplotLegend({'Cowell','USMEM'})

%% Plot USM Elements

usm_type = 'usmem';
load_results = false;
fileName = ''; % 

switch usm_type
    case {'usm7','usm6'}
        limit = 8;
    case 'usmem'
        limit = 7;
    otherwise
        error('USM not recognized')
end

%...Get USM data
% save('data/quat_norm_prop_1_omega','time','usm')
if load_results
    data = load(fileName);
    time = data.time;
    usm = data.usm;
else
    fileID = fopen('/Users/Michele/GitHub/tudat/tudatBundle/tudatApplications/Test/SimulationOutput/test_propagation.dat');
    usm = textscan(fileID,repmat('%f ',[1,limit]),'CollectOutput',true,'Delimiter',',');
    time = (usm{1}(:,1)-usm{1}(1))/3600; usm = usm{1}(:,2:end); usm(:,1:3) = usm(:,1:3)/1e3;
    fclose(fileID);
end

switch usm_type
    case 'usm7'
        usm(:,8) = abs( 1 - sqrt( sum( usm( :, 4:7 ).^2, 2 ) ) );
        labels = {'C Hodograph [km/s]','R_1 Hodograph [km/s]','R_2 Hodograph [km/s]',...
            'q_1 [-]','q_2 [-]','q_3 [-]','q_4 [-]','Norm Offset [-]'};
    case 'usm6'
        usm(usm(:,7) == 0,7) = -1;
        labels = {'C Hodograph [km/s]','R_1 Hodograph [km/s]','R_2 Hodograph [km/s]',...
            '\sigma_1 [-]','\sigma_2 [-]','\sigma_3 [-]'};
    case 'usmem'
        labels = {'C Hodograph [km/s]','R_1 Hodograph [km/s]','R_2 Hodograph [km/s]',...
            'e_1 [-]','e_2 [-]','e_3 [-]'};
end

%...Plot USM elements
figure;
for i = 1:length(labels)
    subplot(3,3,i)
    hold on
    plot(time,usm(:,i),'LineWidth',1.1)
    if strcmp(usm_type,'usm6')
        if i >= 4
            plot(time,usm(:,7),'LineStyle','--')
        end
    end
    hold off
    xlabel('Time [h]')
    ylabel(labels{i})
    grid on
    set(gca,'FontSize',15)
end