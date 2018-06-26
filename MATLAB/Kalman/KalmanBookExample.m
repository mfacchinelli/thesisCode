fclose('all'); clear all; close all force; profile off; clc; format long g; rng default;
addpath Kalman

N_red = 301;

%% Load C++ EKF Output

filename = '/Users/Michele/GitHub/tudat/tudatBundle/tudatApplications/Test/SimulationOutput/KFBook/EKFActualStateHistory.dat';
fileID = fopen(filename,'r');
cpp_result = textscan(fileID,'%f %f %f %f','CollectOutput',true,'Delimiter',',');
state_actual_EKF = cpp_result{1}(:,2:end)';
fclose(fileID);

filename = '/Users/Michele/GitHub/tudat/tudatBundle/tudatApplications/Test/SimulationOutput/KFBook/EKFEstimatedStateHistory.dat';
fileID = fopen(filename,'r');
cpp_result = textscan(fileID,'%f %f %f %f','CollectOutput',true,'Delimiter',',');
state_extimated_ekf = cpp_result{1}(:,2:end)';
fclose(fileID);

filename = ['/Users/Michele/GitHub/tudat/tudatBundle/tudatApplications/Test/'...
    'SimulationOutput/KFBook/EKFEstimatedCovarianceHistory.dat'];
fileID = fopen(filename,'r');
cpp_result = textscan(fileID,'%f %f %f %f %f %f %f %f %f %f','CollectOutput',true,'Delimiter',',');
covariance_extimated_ekf = cpp_result{1}(:,2:end)';
fclose(fileID);

filename = '/Users/Michele/GitHub/tudat/tudatBundle/tudatApplications/Test/SimulationOutput/KFBook/systemNoise.dat';
fileID = fopen(filename,'r');
cpp_result = textscan(fileID,'%f %f %f','CollectOutput',true);
system_noise = cpp_result{1}';
fclose(fileID);

filename = '/Users/Michele/GitHub/tudat/tudatBundle/tudatApplications/Test/SimulationOutput/KFBook/measurementNoise.dat';
fileID = fopen(filename,'r');
cpp_result = textscan(fileID,'%f','CollectOutput',true);
measurement_noise = cpp_result{1}';
fclose(fileID);

%% Load C++ UKF Output

filename = '/Users/Michele/GitHub/tudat/tudatBundle/tudatApplications/Test/SimulationOutput/KFBook/UKFActualStateHistory.dat';
fileID = fopen(filename,'r');
cpp_result = textscan(fileID,'%f %f %f %f','CollectOutput',true,'Delimiter',',');
state_actual_UKF = cpp_result{1}(:,2:end)';
fclose(fileID);

filename = '/Users/Michele/GitHub/tudat/tudatBundle/tudatApplications/Test/SimulationOutput/KFBook/UKFEstimatedStateHistory.dat';
fileID = fopen(filename,'r');
cpp_result = textscan(fileID,'%f %f %f %f','CollectOutput',true,'Delimiter',',');
state_estimated_ukf = cpp_result{1}(:,2:end)';
fclose(fileID);

filename = ['/Users/Michele/GitHub/tudat/tudatBundle/tudatApplications/Test/'...
    'SimulationOutput/KFBook/UKFEstimatedCovarianceHistory.dat'];
fileID = fopen(filename,'r');
cpp_result = textscan(fileID,'%f %f %f %f %f %f %f %f %f %f','CollectOutput',true,'Delimiter',',');
covariance_extimated_ukf = cpp_result{1}(:,2:end)';
fclose(fileID);

%% Plot States

%...Plot both
figure;
subplot(1,3,1)
hold on
plot(1:N_red,state_extimated_ekf(1,1:N_red),'LineWidth',1.25)
plot(1:N_red,state_estimated_ukf(1,1:N_red),'LineWidth',1.25)
plot(1:N_red,state_actual_EKF(1,1:N_red),'LineWidth',1.25)
hold off
grid on
legend('C++ EKF','C++ UKF','Actual','Location','Best')
% legend('C++ EKF','Actual','Location','Best')
set(gca,'FontSize',15)
title('x_1')

subplot(1,3,2)
hold on
plot(1:N_red,state_extimated_ekf(2,1:N_red),'LineWidth',1.25)
plot(1:N_red,state_estimated_ukf(2,1:N_red),'LineWidth',1.25)
plot(1:N_red,state_actual_EKF(2,1:N_red),'LineWidth',1.25)
hold off
grid on
legend('C++ EKF','C++ UKF','Actual','Location','Best')
% legend('C++ EKF','Actual','Location','Best')
set(gca,'FontSize',15)
title('x_2')

subplot(1,3,3)
hold on
plot(1:N_red,state_extimated_ekf(3,1:N_red),'LineWidth',1.25)
plot(1:N_red,state_estimated_ukf(3,1:N_red),'LineWidth',1.25)
plot(1:N_red,state_actual_EKF(3,1:N_red),'LineWidth',1.25)
hold off
grid on
legend('C++ EKF','C++ UKF','Actual','Location','Best')
% legend('C++ EKF','Actual','Location','Best')
set(gca,'FontSize',15)
title('x_3')

%% Plot State Differences

%...Plot both
figure;
subplot(1,3,1)
hold on
plot(1:N_red,state_extimated_ekf(1,1:N_red)-state_actual_EKF(1,1:N_red),'LineWidth',1.25)
plot(1:N_red,state_estimated_ukf(1,1:N_red)-state_actual_EKF(1,1:N_red),'LineWidth',1.25)
plot(1:N_red,sqrt(covariance_extimated_ekf(1,1:N_red)),'LineWidth',1.25,'Color',[0.929,0.694,0.125])
plot(1:N_red,-sqrt(covariance_extimated_ekf(1,1:N_red)),'LineWidth',1.25,'Color',[0.929,0.694,0.125])
hold off
grid on
legend('C++ EKF','C++ UKF','STD','Location','Best')
% legend('C++ EKF','STD','Location','Best')
set(gca,'FontSize',15)
title('x_1')

subplot(1,3,2)
hold on
plot(1:N_red,state_extimated_ekf(2,1:N_red)-state_actual_EKF(2,1:N_red),'LineWidth',1.25)
plot(1:N_red,state_estimated_ukf(2,1:N_red)-state_actual_EKF(2,1:N_red),'LineWidth',1.25)
plot(1:N_red,sqrt(covariance_extimated_ekf(5,1:N_red)),'LineWidth',1.25,'Color',[0.929,0.694,0.125])
plot(1:N_red,-sqrt(covariance_extimated_ekf(5,1:N_red)),'LineWidth',1.25,'Color',[0.929,0.694,0.125])
hold off
grid on
legend('C++ EKF','C++ UKF','STD','Location','Best')
% legend('C++ EKF','STD','Location','Best')
set(gca,'FontSize',15)
title('x_2')

subplot(1,3,3)
hold on
plot(1:N_red,state_extimated_ekf(3,1:N_red)-state_actual_EKF(3,1:N_red),'LineWidth',1.25)
plot(1:N_red,state_estimated_ukf(3,1:N_red)-state_actual_EKF(3,1:N_red),'LineWidth',1.25)
plot(1:N_red,sqrt(covariance_extimated_ekf(9,1:N_red)),'LineWidth',1.25,'Color',[0.929,0.694,0.125])
plot(1:N_red,-sqrt(covariance_extimated_ekf(9,1:N_red)),'LineWidth',1.25,'Color',[0.929,0.694,0.125])
hold off
grid on
legend('C++ EKF','C++ UKF','STD','Location','Best')
% legend('C++ EKF','STD','Location','Best')
set(gca,'FontSize',15)
title('x_3')