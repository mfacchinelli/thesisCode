fclose('all'); clear all; close all force; profile off; clc; format long g;

%...General
radius = 3.396e6;
mu = 4.282e13;
% radius = 695.508e6;
% mu = 1.32712440019e20;

%% Initial Conditions

clc;

%...Elliptical orbit
periapsis = 165e3+radius;
apoapsis = 47500e3+radius;
% periapsis = 200e3+radius;
% apoapsis = 250e3+radius;
% periapsis = 100e6+radius;
% apoapsis = 200e6+radius;

a = (apoapsis+periapsis)/2;
e = (apoapsis-periapsis)/(apoapsis+periapsis);

%...Print
fprintf('Semi-major axis: %f\n',a)
fprintf('Eccentricity: %f\n',e)

%% Hyperbolic Initial Conditions

clc;

%...Hyperbolic orbit
periapsis = 65e3+radius;
e = 1.2;

%...Print
a = periapsis/(1-e);
t = acos(-1/e);
fprintf('Semi-major axis: %f\n',a)
fprintf('Eccentricity: %f\n',e)
fprintf('Asymptotyc True Anomaly: %f\n',rad2deg(t))
sqrt(mu*(2/(a*(1-1.2))-1/a))

%% Other Conditions

clc;

%...Elliptical orbit
semiMajorAxis = 31241000;
eccentricity = 0.89;

rp = (semiMajorAxis*(1-eccentricity) - radius)/1e3;
ra = (semiMajorAxis*(1+eccentricity) - radius)/1e3;

%...Print
fprintf('Periapsis: %f\n',rp)
fprintf('Apoapsis: %f\n',ra)

%% Access Data

start = 5950;
stop = 6450;
start = 1;
stop = 18802;
% start = 1;
% stop = 518401;

%...Get orbital data
fileID = fopen('/Users/Michele/GitHub/tudat/tudatBundle/tudatApplications/Test/SimulationOutput/Aerobraking/orbit.dat');
orbit = textscan(fileID,repmat('%f ',[1,7]),'CollectOutput',true,'Delimiter',',');
time = (orbit{1}(start:stop,1)-orbit{1}(1))/3600/24; orbit = orbit{1}(start:stop,2:end);
orbit(:,1) = orbit(:,1)/1e3; orbit(:,3:end) = rad2deg(orbit(:,3:end));
fclose(fileID);

%...Get acceleration data
fileID = fopen('/Users/Michele/GitHub/tudat/tudatBundle/tudatApplications/Test/SimulationOutput/Aerobraking/dependent.dat');
dependent = textscan(fileID,repmat('%f ',[1,9]),'CollectOutput',true,'Delimiter',',');
dependent = dependent{1}(start:stop,:);
dependent(abs(dependent)<1e-15) = 0;
fclose(fileID);

altitude = 1e3*orbit(:,1).*(1-orbit(:,2).^2)./(1+orbit(:,2).*cosd(orbit(:,6))) - radius;

%% Plot orbit

figure;
labels = {'Semi-major Axis [km]','Eccentricity [-]','Inclination [deg]',...
    'Right Ascension of Ascending Node [deg]','Argument of Perigee [deg]','True Anomaly [deg]'};
for i = 1:6
    subplot(2,3,i)
    plot(time,orbit(:,i),'LineWidth',1.1)
    xlabel('Time [h]')
    ylabel(labels{i})
    grid on
    set(gca,'FontSize',15)
end

%% Plot accelerations

figure;
styles = {'-.','-','--',':','-','-.','--',':'};
hold on
for i = 1:size(dependent,2)-1
    plot(time,dependent(:,i+1),'LineWidth',1.5,'LineStyle',styles{i})
end
% annotation('textarrow',[0.6,0.5],[0.895,0.865],'String','Periareion','FontSize',25)
hold off
xlabel('Time [h]')
ylabel('Acceleration [m s^{-2}]')
legend('Mars (gravity)','Mars (aerodynamic)','Sun (radiation)','Sun (gravity)',...
    'Venus (gravity)','Earth (gravity)','Jupiter (gravity)','Saturn (gravity)','Location','Best')
grid on
set(gca,'FontSize',25,'YScale','log')

%% Aerodynamic Acceleration

figure;
yyaxis left
plot(time,dependent(:,3),'LineWidth',1.5)
ylabel('Acceleration [m/s^2]')
set(gca,'FontSize',25,'YScale','log')
yyaxis right
plot(time(1:end-1),abs(diff(dependent(:,3))./diff(time)),'LineWidth',1.5)
ylabel('Time Derivative of Acceleration [m/s^3]')
xlabel('Time [h]')
set(gca,'FontSize',25,'YScale','log')
grid on

%% Plot Accelerations Until Periapsis

reduce = find(altitude == min(altitude));

figure;
styles = {'-.','-','--',':','-','-.','--',':'};
hold on
for i = 1:8
    plot(altitude(1:reduce)/1e3,dependent(1:reduce,i+1),'LineWidth',1.5,'LineStyle',styles{i})
end
hold off
xlabel('Altitude [km]')
ylabel('Acceleration [m/s^2]')
legend('Mars (gravity)','Mars (aerodynamic)','Sun (radiation)','Sun (gravity)',...
    'Venus (gravity)','Earth (gravity)','Jupiter (gravity)','Saturn (gravity)','Location','Best')
grid on
set(gca,'FontSize',25,'XScale','log','YScale','log')