fclose('all'); clear all; close all force; profile off; clc; format long g; rng default;
addpath data functions

%% Initial Conditions

%...Planet
radius = 3.390e6;
mu = 4.282e13;
omega_mars = 2*pi/(24*3600 + 39*60 + 35.244);
atm_interface = 200e3;

%...Satellite parameters
m = 1000;
A = 37.5;
CD = 2.2;
CL = 0;
K = m/A/CD;

%% Access Data

start = 6050;
stop = 6350;

%...Get orbital data
fileID = fopen('/Users/Michele/GitHub/tudat/tudatBundle/tudatApplications/Test/SimulationOutput/Aerobraking/orbit.dat');
orbit = textscan(fileID,repmat('%f ',[1,7]),'CollectOutput',true,'Delimiter',',');
time = (orbit{1}(start:stop,1)-orbit{1}(1))/3600; orbit = orbit{1}(start:stop,2:end); 
fclose(fileID);

%...Get acceleration data
fileID = fopen('/Users/Michele/GitHub/tudat/tudatBundle/tudatApplications/Test/SimulationOutput/Aerobraking/dependent.dat');
dependent = textscan(fileID,repmat('%f ',[1,9]),'CollectOutput',true,'Delimiter',',');
dependent = dependent{1};
dependent(abs(dependent)<1e-15) = 0;
sphe = dependent(start:stop,2);
aero = dependent(start:stop,3);
rad = dependent(start:stop,4);
sun = dependent(start:stop,5);
ven = dependent(start:stop,6);
ear = dependent(start:stop,7);
jup = dependent(start:stop,8);
sat = dependent(start:stop,9);
fclose(fileID);

%...Get altitude
radial = orbit(:,1).*(1-orbit(:,2).^2)./(1+orbit(:,2).*cos(orbit(:,6)));
altitude = radial - radius;
velocity = sqrt(mu*(2./radial-1./orbit(:,1)));

%...Adapt
loc = true(size(altitude));%altitude <= atm_interface;
time = time(loc);
orbit = orbit(loc,:);
sphe = sphe(loc); aero = aero(loc); rad = rad(loc); sun = sun(loc);
ven = ven(loc); ear = ear(loc); jup = jup(loc); sat = sat(loc);
radial = radial(loc);
velocity = velocity(loc,:);
altitude = altitude(loc,:);

figure;
accs = {sphe,aero,rad,sun,ven,ear,jup,sat};
hold on
for i = 1:length(accs)
    plot(time,accs{i})
end
hold off
xlim([time(1),time(end)])
set(gca,'yscale','log')
legend('sphe','aero','rad','sun','ven','ear','jup','sat')
grid on

%...IMU model
imu_model = 'KVH';
% mis = @(x) [0,x(1),x(2);x(3),0,x(4);x(5),x(6),0];
switch imu_model
    case 'custom'
        std_imu = 1e-3;
        bias = [0.1;-0.2;0.05];
        scale = [0.025;0.1;-0.075];
        misalign = [0.01;-0.02;0.01;-0.05;0.03;0.01];
    case 'KVH'
        g0 = 9.81;
        std_imu = 0.024e-3*g0;
        bias = 0.4e-3*g0*randn(3,1)*3;
        scale = 0.8e-2*randn(3,1)*3;
        misalign = 1e-3*randn(6,1)*3;
end        
% imu_error = [bias;scale;misalign];
% acc_imu = bias + (eye(3) + diag(scale) + mis(misalign)) * (aero+rad) + std_imu*randn(size(aero));
% acc_imu = (eye(3) + diag(scale) + mis(misalign))\(acc_imu-bias);
acc_imu = abs(sum([aero,rad,sun,ven,ear,jup,sat],2) + std_imu*randn(size(aero)));

%% Find Area Division

%...Start plotting
figure;
hold on
plot(time,acc_imu)
% plot(time,aero,'LineWidth',2)

%...First iteration
a0 = find(max(acc_imu)==acc_imu)-100;
b0 = find(max(acc_imu)==acc_imu)+100;

%...Iterative method
c = areaBisection(time,acc_imu,a0,b0);

%...Improvement
fprintf('Before: %.3f.\n',time(max(acc_imu)==acc_imu) * 60)
fprintf('After: %.3f.\n',time(c) * 60)

%...Finish plotting
plot([time(c),time(c)],ylim,'LineWidth',2,'LineStyle','--')
hold off
xlabel('Time [h]')
ylabel('Acceleration Norm [m/s^2]')
legend('IMU','Actual','Barycenter')
grid on

%% Delta V and period change

figure;
yyaxis left
plot(time,orbit(:,1),'LineWidth',1.1)
ylabel('Semi-major axis [m]')
yyaxis right
plot(time,aero,'LineWidth',1.1)
ylabel('Acceleration Norm [m/s^2]')
xlabel('Time [h]')
set(gca,'YScale','log')
grid on

%...Find Delta V
DV_imu = -trapz(time*3600,acc_imu); % time is in hours
DV = -trapz(time*3600,aero);

%...Find Delta a
a0 = orbit(1,1); a1 = orbit(end,1); n0 = sqrt(mu/a0^3);
e0 = orbit(1,2); e1 = orbit(end,2);
Da_imu = 2/n0*sqrt((1+e0)/(1-e0))*DV_imu;
De = 2*sqrt(a0*(1-e0^2)/mu)*DV_imu;

Da = 2/n0*sqrt((1+e0)/(1-e0))*DV;
Da_actual = a1 - a0;
De_actual = e1 - e0;

DV_actual = Da_actual/(2/n0*sqrt((1+e0)/(1-e0)));

table(DV_imu,DV,DV_actual)
table(Da_imu,Da,Da_actual)

%...Find Delta Period
kappa = 0.955;
DP_imu = 2*pi*kappa*((a0+Da_imu)^(3/2)-a0^(3/2))/sqrt(mu)/3600;

P = 2*pi*sqrt(a0^3/mu)/3600;
F_MGS = -7.097 + 5.4596 * P + 0.8577 * P^2 - 0.00361 * P^3;
DP_MGS = 0.96 * F_MGS * DV_imu / 3600;

DP_actual = 2*pi*(sqrt(a1^3/mu)-sqrt(a0^3/mu))/3600;
table(DP_imu,DP_MGS,DP_actual)

%% Plot

%...Plot orbit
figure;
for i = 1:6
    subplot(2,3,i)
    plot(time,orbit(:,i),'LineWidth',1.1)
    xlabel('Time [h]')
    grid on
end

%...Plot altitude
% figure;
% plot(time,altitude,'LineWidth',1.1)
% ylabel('Altitude [m]')
% xlabel('Time [h]')
% grid on

%...Smoothed data
% figure;
% plot(time,aero_imu,'LineWidth',1.1)
% ylabel('Acceleration Norm [m/s^2]')
% xlabel('Time [h]')
% grid on

%% Acceleration norm and true anomaly

figure;
hold on
plot(time,aero,'LineWidth',1.1)
a = 1:50:500;
for i = a
    plot(time,smooth(acc_imu,i),'LineWidth',1.1)
end
hold off
legend('True',cellstr(num2str(a')))
grid on

figure;
yyaxis left
hold on
plot(time,aero,'LineWidth',1.1)
plot(time,smooth(acc_imu,200),'LineWidth',1.1)
hold off
ylabel('Acceleration Norm [m/s^2]')
% set(gca,'yscale','log')
yyaxis right
plot(time,orbit(:,6),'LineWidth',1.1,'LineStyle','-.')
ylabel('True Anomaly [deg]')
xlabel('Time [h]')
legend('Actual','Smoothed','Anomaly')
grid on

%% Fit Density

%...Get before peak
loc = ( altitude <= atm_interface ) & ( true(size(find(altitude==min(altitude)),1)) );
z = altitude(loc);

%...Get density
V = velocity(loc);
rho = acc_imu(loc)*2*K./V.^2/sqrt(1+(CL/CD)^2);

x0 = [log(0.02);1/11.1e3];
A = zeros(length(z),2);
for i = 1:length(z)
    A(i,:) = - [x0(2),x0(1)] .* z(i);
end
x_hat = (A'*A)\A'*log(rho);

table([exp(x_hat(1));1/x_hat(2)],[0.02;11.1e3])

rho_lsq = @(h) exp(x_hat(1)) * exp( - h * x_hat(2) );

%...Plot
figure;
hold on
plot(z,rho_lsq(z))
scatter(z,rho)
% plot(z,rho_true(z))
hold off
set(gca,'yscale','log')
grid on

% %...Fit
% ft = fittype('rho0 - altitude / H',...
%     'independent',{'altitude'},...
%     'dependent',{'rho'},...
%     'coefficients',{'rho0','H'});
% lsq = fit(z,log(rho),ft,'StartPoint',[log(0.01),10e3],'Lower',[log(0),0],'Robust',...
%     'LAR','Algorithm','Trust-Region');
% rho_lsq = @(h) exp(lsq.rho0) * exp( - h / lsq.H );
% 
% table([exp(lsq.rho0);lsq.H],[0.02;11.1e3])

% %...Plot
% figure;
% hold on
% plot(z,rho_lsq(z))
% scatter(z,rho)
% % plot(z,rho_true(z))
% hold off
% set(gca,'yscale','log')
% grid on