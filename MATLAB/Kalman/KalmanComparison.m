fclose('all'); clear all; close all force; profile off; clc; format long g; rng default;
addpath Kalman

%% Load MATLAB EKF and UKF Output

mat_ekf = load('/Users/Michele/GitHub/tudat/tudatBundle/tudatApplications/Test/SimulationOutput/KF/ekf.mat');
mat_ukf = load('/Users/Michele/GitHub/tudat/tudatBundle/tudatApplications/Test/SimulationOutput/KF/ukf.mat');

%% General Output

filename = '/Users/Michele/GitHub/tudat/tudatBundle/tudatApplications/Test/SimulationOutput/KF/actualStateHistory.dat';
fileID = fopen(filename,'r');
cpp_result = textscan(fileID,'%f %f %f','CollectOutput',true,'Delimiter',',');
x_cpp = cpp_result{1}(:,2:end)';
fclose(fileID);

filename = '/Users/Michele/GitHub/tudat/tudatBundle/tudatApplications/Test/SimulationOutput/KF/measurementHistory.dat';
fileID = fopen(filename,'r');
cpp_result = textscan(fileID,'%f %f','CollectOutput',true,'Delimiter',',');
z_cpp = cpp_result{1}(:,2:end)';
fclose(fileID);

%% Load C++ EKF Output

filename = '/Users/Michele/GitHub/tudat/tudatBundle/tudatApplications/Test/SimulationOutput/KF/EKFEstimatedStateHistory.dat';
% filename = '/Users/Michele/Desktop/Filter/EKFEstimatedStateHistory.dat';
fileID = fopen(filename,'r');
cpp_result = textscan(fileID,'%f %f %f','CollectOutput',true,'Delimiter',',');
x_cor_cpp = cpp_result{1}(:,2:end)';
fclose(fileID);

%% Load C++ UKF Output

filename = '/Users/Michele/GitHub/tudat/tudatBundle/tudatApplications/Test/SimulationOutput/KF/UKFEstimatedStateHistory.dat';
% filename = '/Users/Michele/Desktop/Filter/UKFEstimatedStateHistory.dat';
fileID = fopen(filename,'r');
cpp_result = textscan(fileID,'%f %f %f','CollectOutput',true,'Delimiter',',');
x_cor_cpp_ukf = cpp_result{1}(:,2:end)';
fclose(fileID);

%% Plots

N_red = 999;

%...Plot both
figure;
for i = 1:2
    subplot(1,2,i)
    hold on
    plot(1:N_red,mat_ekf.x_cor(i,1:N_red)-x_cpp(i,1:N_red),'LineWidth',1.25)
    plot(1:N_red,mat_ukf.xx(i,1:N_red)-x_cpp(i,1:N_red),'LineWidth',1.25)
    plot(1:N_red,x_cor_cpp(i,1:N_red)-x_cpp(i,1:N_red),'LineWidth',1.25)
    plot(1:N_red,x_cor_cpp_ukf(i,1:N_red)-x_cpp(i,1:N_red),'LineWidth',1.25)
    hold off
    grid on
    legend('MATLAB EKF','MATLAB UKF','C++ EKF','C++ UKF','Location','Best')
    set(gca,'FontSize',15)
    title(['x_',num2str(i),])
end