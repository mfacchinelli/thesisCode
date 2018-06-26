% UKF.m
% Unscented Kalman Filter (UKF)

%...Clean up
fclose('all'); clear all; close all force; profile off; clc; format long g; rng default;
addpath data other Kalman

%% Access Data

filename = '/Users/Michele/GitHub/tudat/tudatBundle/tudatApplications/Test/SimulationOutput/KF/EKFMeasurementHistory.dat';
fileID = fopen(filename,'r');
cpp_result = textscan(fileID,'%f %f','CollectOutput',true,'Delimiter',',');
z_k = cpp_result{1}(:,2:end)';
fclose(fileID);

%% UKF Data

%...Simulation data
N = 1000; 	% N: number of epochs
n = 2;      % n: state dimension
m = 1;      % m: observation dimension
L = n*2+m;	% L: augmented state dimension
S = 2*L+1;	% S: number of sigma points
dt = 0.01;

%...Scaling
model = 'mst';
switch model
    case 'mst'
        alpha = 0.003;
        kappa = 0;
    case 'art'
        alpha = 0.05;
        kappa = 3-L;
    case 'ieee'
        alpha = 0.001;
        kappa = 1;
    case 'test'
        alpha = 1;
        kappa = 1;
end
lambda = alpha^2*(L+kappa) - L;
gamma = sqrt(L+lambda);

%...Weights
beta = 2;
W_x = horzcat(lambda/(L+lambda),1/2/(L+lambda)*ones(1,S-1));
W_P = horzcat(lambda/(L+lambda)+(1-alpha^2+beta),1/2/(L+lambda)*ones(1,S-1));

%...System noise
Ew = zeros(n,1);
stdw = 1e1 * ones(n,1);
Q = diag(stdw.^2); 

%...Measurement noise
Ev = zeros(m,1);
stdv = 1e1 * ones(m,1);
R = diag(stdv.^2);

%% Unscented Kalman Filter (UKF)

clc;

%...Initial conditions
x_k1_k1 = [10;-3];	% x(0|0): augmented state vector
Pxx_k1_k1 = Q;      % P(0|0): covariance matrix

%...Pre-allocate arrays
x = 1:n; w = n+1:2*n; v = 2*n+1:L;
x_est = zeros(n,N-1);
stdx_est = zeros(n,N-1);
xs = zeros(n,S,N-1);
xx = zeros(n,N-1);
zz = zeros(m,N-1);

%...Loop over time
tic
for k = 1:N-1
    %...Augmented parameters
    xa_k1_k1 = vertcat(x_k1_k1,Ew,Ev); % augmented state vector
    P_k1_k1 = blkdiag(Pxx_k1_k1,Q,R); % augmented covariance matrix
    
    %...Sigma points
    xa_k1_k1 = sigmaPoint(gamma,xa_k1_k1,P_k1_k1);
        
    %...State prediction
    x_k_k1 = xa_k1_k1(x,:) + ( state(xa_k1_k1(x,:)) + xa_k1_k1(w,:) ) * dt;  % x(k|k-1): a-priori
    xx_k_k1 = sum(W_x.*x_k_k1,2);
    
    e_x = x_k_k1 - xx_k_k1;
    Pxx_k_k1 = (W_P.*e_x)*e_x';     % P(k|k-1): a-priori
    
    %...Sigma points (again)
    xa_k_k1 = sigmaPoint(gamma,vertcat(xx_k_k1,Ew,Ev),blkdiag(Pxx_k_k1,Q,R));
    
    %...Measurement prediction
    z_k_k1 = measurement(xa_k_k1(x,:)) + xa_k_k1(v,:);
    zz_k_k1 = sum(W_x.*z_k_k1,2);
    
    e_z = z_k_k1 - zz_k_k1;
    Pzz = (W_P.*e_z)*e_z'; 	% innovation
    
    %...State-measurement covariance
    Pxz = (W_P.*e_x)*e_z';  % cross-correlation
    
    %...Kalman gain
    K = Pxz/Pzz;	% K(k): Kalman gain
        
    %...Correction
    x_k_k = xx_k_k1 + K*(z_k(k)-zz_k_k1);                 % x(k|k): a-posteriori
    Pxx_k_k = Pxx_k_k1 - K*Pzz*K';                      % P(k|k): a-posteriori
    
    %...Store values
    x_est(:,k) = x_k_k;                     % estimated state
    stdx_est(:,k) = sqrt(diag(Pxx_k_k));	% standard deviation
    zz(:,k) = zz_k_k1;
    xx(:,k) = xx_k_k1;
    xs(:,:,k) = xa_k1_k1(x,:);%x_k_k1;
    
    %...Next step
    x_k1_k1 = x_k_k;
    Pxx_k1_k1 = Pxx_k_k;
end
toc

%% Final State

save('/Users/Michele/GitHub/tudat/tudatBundle/tudatApplications/Test/SimulationOutput/KF/ukf.mat','xx')

%...Show final state
x_k1_k1

%% Plot Sigma Points

offset = 10;
y1 = round(length(xs)/2)-offset;
y2 = round(length(xs)/2)+offset;

figure; 
hold on
plot(xx(1,y1:y2),xx(2,y1:y2),'LineWidth',1.25)
plot(x_est(1,y1:y2),x_est(2,y1:y2),'LineWidth',1.25,'LineStyle','--')
for i = y1:y2
    scatter(squeeze(xs(1,:,i)),squeeze(xs(2,:,i)),'filled')
end
hold off
xlabel('x [m]'); ylabel('y [m]');
legend('A-priori Mean','A-posteriori','Location','Best')
grid on
set(gca,'FontSize',15)

%% Functions

function sigma = sigmaPoint(g,x,P)
    sqrtP = g * real(sqrtm(P));
    sigma = x + horzcat(zeros(size(x)),sqrtP,-sqrtP);
end

function x_dot = state(x)
    x_dot( 1, : ) = x(2,:) .* cos( x(1,:) ).^3;
    x_dot( 2, : ) = sin( x(1,:) );
end

function z = measurement(x)
    z = x( 1, : ).^3;
end