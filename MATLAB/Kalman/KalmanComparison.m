fclose('all'); clear all; close all force; profile off; clc; format long g; rng default;
addpath Kalman

%% Load MATLAB EKF and UKF Output

mat_ekf = load('/Users/Michele/Desktop/KF/ekf.mat');
mat_ukf = load('/Users/Michele/Desktop/KF/ukf.mat');

%% Load C++ EKF Output

fileID = fopen('/Users/Michele/Desktop/KF/EKFActualStateHistory.dat','r');
cpp_result = textscan(fileID,'%f %f %f','CollectOutput',true,'Delimiter',',');
x_cpp = cpp_result{1}(:,2:end)';
fclose(fileID);

fileID = fopen('/Users/Michele/Desktop/KF/EKFEstimatedStateHistory.dat','r');
cpp_result = textscan(fileID,'%f %f %f','CollectOutput',true,'Delimiter',',');
x_cor_cpp = cpp_result{1}(:,2:end)';
fclose(fileID);

fileID = fopen('/Users/Michele/Desktop/KF/EKFMeasurementHistory.dat','r');
cpp_result = textscan(fileID,'%f %f','CollectOutput',true,'Delimiter',',');
z_cpp = cpp_result{1}(:,2:end)';
fclose(fileID);

%% Load C++ UKF Output

% fileID = fopen('/Users/Michele/Desktop/KF/UKFActualStateHistory.dat','r');
% cpp_result = textscan(fileID,'%f %f %f','CollectOutput',true,'Delimiter',',');
% x_cpp_ukf = cpp_result{1}(:,2:end)';
% fclose(fileID);

fileID = fopen('/Users/Michele/Desktop/KF/UKFEstimatedStateHistory.dat','r');
cpp_result = textscan(fileID,'%f %f %f','CollectOutput',true,'Delimiter',',');
x_cor_cpp_ukf = cpp_result{1}(:,2:end)';
fclose(fileID);

% fileID = fopen('/Users/Michele/Desktop/KF/UKFMeasurementHistory.dat','r');
% cpp_result = textscan(fileID,'%f %f','CollectOutput',true,'Delimiter',',');
% z_cpp_ukf = cpp_result{1}(:,2:end)';
% fclose(fileID);

%% Plots

N_red = 999;

%...Plot both
figure;
subplot(1,2,1)
hold on
plot(1:N_red,mat_ekf.x_cor(1,1:N_red),'LineWidth',1.25)
plot(1:N_red,mat_ukf.xx(1,1:N_red),'LineWidth',1.25)
plot(1:N_red,x_cor_cpp(1,1:N_red),'LineWidth',1.25)
plot(1:N_red,x_cor_cpp_ukf(1,1:N_red),'LineWidth',1.25)
plot(1:N_red,x_cpp(1,1:N_red),'LineWidth',1.25)
hold off
grid on
legend('MATLAB EKF','MATLAB UKF','C++ EKF','C++ UKF','Actual','Location','Best')
set(gca,'FontSize',15)
title('x_1')

subplot(1,2,2)
hold on
plot(1:N_red,mat_ekf.x_cor(2,1:N_red),'LineWidth',1.25)
plot(1:N_red,mat_ukf.xx(2,1:N_red),'LineWidth',1.25)
plot(1:N_red,x_cor_cpp(2,1:N_red),'LineWidth',1.25)
plot(1:N_red,x_cor_cpp_ukf(2,1:N_red),'LineWidth',1.25)
plot(1:N_red,x_cpp(2,1:N_red),'LineWidth',1.25)
hold off
grid on
legend('MATLAB EKF','MATLAB UKF','C++ EKF','C++ UKF','Actual','Location','Best')
set(gca,'FontSize',15)
title('x_2')