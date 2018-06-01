% EKF.m
% Extended Kalman Filter (UKF)

%...Clean up
fclose('all'); clear all; close all force; profile off; clc; format long g; rng default;
addpath data other Kalman

%% Access Data

fileID = fopen('/Users/Michele/Desktop/KF/EKFMeasurementHistory.dat','r');
cpp_result = textscan(fileID,'%f %f','CollectOutput',true,'Delimiter',',');
z_k = cpp_result{1}(:,2:end)';
fclose(fileID);

%% Extended Kalman Filter (UKF)

clc;

%...Initial conditions
Ex_0 = 10; stdx_0 = 10; P_0 = stdx_0^2; x_0 = [3; -0.3];

%...Simulation data
dt = 0.01; N = 1000;
t_end = N*dt;

%...System noise
Ew = 0; stdw = 10; 
Q = stdw^2;

%...Measurement noise
Ev = 0; stdv = 10; 
R = stdv^2;

%...Input of constants
Ts = dt;
x_k1_k1 = [10;-3];          % x(0|0) = E(x_0)
P_k1_k1 = diag([100,100]);	% P(0|0) = P(0)

%...Extended Kalman Filter (EKF)
ti = 0; tf = Ts;  
n = length(x_k1_k1); % n: state dimension
for k = 1:N
    %...First
    x_cor(:,k) = x_k1_k1;
    stdx_cor(:,k) = sqrt(diag(P_k1_k1));	% standard deviation
    
    %...Prediction
%     [t,x] = ode45(@state,[ti tf],x_k1_k1);  % predicted states
%     x_k_k1 = x(length(t),:)'               % x(k|k-1)
    x_k_k1 = x_k1_k1 + state(ti,x_k1_k1) * Ts;
    x_pred(:,k) = x_k_k1; 
    z_k_k1 = measurement(x_k_k1,0);       	% z(k|k-1)
    z_pred(:,k) = z_k_k1;
    
    [Phi,Gamma] = stateJacobian(x_k1_k1,Ts);  	% Phi(k|k-1) & Gamma(k|k-1)
    P_k_k1 = Phi*P_k1_k1*Phi'+Gamma*Q*Gamma';	% P(k|k-1) (prediction) 
    P_pred(:,k) = diag(P_k_k1);
    stdx_pred(:,k) = sqrt(diag(P_k_k1));
    
    %...Correction
    H = measurementJacobian(x_k_k1,n);
    Ve = (H*P_k_k1*H'+R);               % covariance matrix of innovation
    Vek(:,k) = sqrt(diag(Ve)); 
    K = P_k_k1*H'/Ve;                   % K(k) (Kalman gain)
    
    x_k_k = x_k_k1+K*(z_k(:,k)-z_k_k1);                 % x(k|k) (correction)
    P_k_k = (eye(n)-K*H)*P_k_k1*(eye(n)-K*H)'+K*R*K';   % P(k|k) (correction)
    
    %...Next step
    x_k1_k1 = x_k_k; 
    P_k1_k1 = P_k_k; 
    ti = tf; tf = tf+Ts;
end

%% Final State

save('/Users/Michele/Desktop/KF/ekf.mat','x_cor')

%...Show final state
x_k1_k1

%% Functions

function x_dot = state(~,x)
    x_dot(1,1) = x(2)*cos(x(1))^3; 
    x_dot(2,1) = sin(x(1));
end

function z = measurement(x,~)
    z(1,1) = x(1)^3;
end

function [Phi,Gamma,F,G] = stateJacobian(x,Ts)
    F(1,1) = x(2)*(-3*cos(x(1))^2*sin(x(1))); 
    F(1,2) = cos(x(1))^3;
    F(2,1) = cos(x(1));
    F(2,2) = 0;
    G = [1;1];
    [Phi,Gamma] = c2dLocal(F,G,Ts);
end

function H = measurementJacobian(x,~) 
    H = [3*x(1)^2,0];
end

function [Phi,Gamma] = c2dLocal(F,G,Ts)
    Mat = [F,G;zeros(1,3)];
    MatExp = eye(3) + Mat * Ts + Mat^2 * Ts^2 / 2 + Mat^3 * Ts^3 / 6;
    Phi = MatExp(1:2,1:2);
    Gamma = MatExp(1:2,3);
end