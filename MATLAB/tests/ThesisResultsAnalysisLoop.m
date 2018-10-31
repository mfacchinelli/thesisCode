fclose('all'); clear all; close all force; profile off; clc; format long g; rng default;
addpath functions tests

%% Settings

%...Constants
marsRadius = 3389526.666666667;
marsGravitationalParameter = 42828375815756.1;
marsAtmosphericInterface = 175;

%...Plot settings
applyInterpolation = true;

%...Labels
timeConversion = 3600 * 24;
timeLabel = 'Time [d]';
CartesianLabels = {'x [km]','y [km]','z [km]','v_x [m s^{-1}]','v_y [m s^{-1}]','v_z [m s^{-1}]'};
CartesianLabelsDifference = {'\Delta x [km]','\Delta y [km]','\Delta z [km]',...
    '\Delta v_x [m s^{-1}]','\Delta v_y [m s^{-1}]','\Delta v_z [m s^{-1}]'};
KeplerianLabels = {'a [km]','e [-]','i [deg]','\omega [deg]','\Omega [deg]','\vartheta [deg]'};
rotationLabels = {'\eta [-]','\epsilon_1 [-]','\epsilon_2 [-]','\epsilon_3 [-]','Norm Offset [-]',...
    '\omega_1 [deg s^{-1}]','\omega_2 [deg s^{-1}]','\omega_3 [deg s^{-1}]'};

ps = -1:3;
vs = -3:0;
ms = -1:2;

results = zeros(length(ps),length(vs),length(ms));
for p = ps
    for v = vs
        for m = ms
            %% Load C++ Results For Propagation
            
            outputFolder = ['SimulationOutputTransOnlyIMANLoop/',num2str(p),'-',num2str(v),'-',num2str(m),'/'];
            
            %...Load translational motion
            filename = ['/Users/Michele/GitHub/tudat/tudatBundle/tudatApplications/Thesis/',outputFolder,'/cartesianPropagated.dat'];
            fileID = fopen(filename,'r');
            CartesianPropagatedResults = textscan(fileID,repmat('%f',[1,7]),'Delimiter',',','CollectOutput',true);
            initialTime = CartesianPropagatedResults{1}(1);
            simulationTime = ( CartesianPropagatedResults{1}(:,1) - initialTime ) / timeConversion;
            CartesianPropagatedResults = CartesianPropagatedResults{1}(:,2:end);
            CartesianPropagatedResults(:,1:3) = CartesianPropagatedResults(:,1:3) / 1e3;
            fclose(fileID);
            
            %...Load Keplerian translational motion
            filename = ['/Users/Michele/GitHub/tudat/tudatBundle/tudatApplications/Thesis/',outputFolder,'/keplerianPropagated.dat'];
            fileID = fopen(filename,'r');
            KeplerianPropagatedResults = textscan(fileID,repmat('%f',[1,7]),'Delimiter',',','CollectOutput',true);
            KeplerianPropagatedResults = KeplerianPropagatedResults{1}(:,2:end);
            KeplerianPropagatedResults(:,1) = KeplerianPropagatedResults(:,1) / 1e3;
            KeplerianPropagatedResults(:,3:end) = rad2deg(KeplerianPropagatedResults(:,3:end));
            fclose(fileID);
            
            % %...Load rotational motion
            % filename = ['/Users/Michele/GitHub/tudat/tudatBundle/tudatApplications/Thesis/',outputFolder,'/rotationalPropagated.dat'];
            % fileID = fopen(filename,'r');
            % rotationalPropagatedResults = textscan(fileID,repmat('%f',[1,8]),'Delimiter',',','CollectOutput',true);
            % rotationalPropagatedResults = rotationalPropagatedResults{1}(:,2:end);
            % rotationalPropagatedResults(:,6:8) = rad2deg(rotationalPropagatedResults(:,5:7));
            % rotationalPropagatedResults(:,5) = 1.0 - quatnorm(rotationalPropagatedResults(:,1:4));
            % fclose(fileID);
            
            %...Clean up
            clear filename fileID
            
            %% Load C++ Results For Estimation
            
            %...Load translational motion
            filename = ['/Users/Michele/GitHub/tudat/tudatBundle/tudatApplications/Thesis/',outputFolder,'/cartesianEstimated.dat'];
            fileID = fopen(filename,'r');
            CartesianEstimatedResults = textscan(fileID,repmat('%f',[1,7]),'Delimiter',',','CollectOutput',true);
            onboardTime = ( CartesianEstimatedResults{1}(:,1) - CartesianEstimatedResults{1}(1) ) / timeConversion;
            CartesianEstimatedResults = CartesianEstimatedResults{1}(:,2:end);
            CartesianEstimatedResults(:,1:3) = CartesianEstimatedResults(:,1:3) / 1e3;
            fclose(fileID);
            
            %...Load Keplerian translational motion
            filename = ['/Users/Michele/GitHub/tudat/tudatBundle/tudatApplications/Thesis/',outputFolder,'/keplerianEstimated.dat'];
            fileID = fopen(filename,'r');
            KeplerianEstimatedResults = textscan(fileID,repmat('%f',[1,7]),'Delimiter',',','CollectOutput',true);
            KeplerianEstimatedResults = KeplerianEstimatedResults{1}(:,2:end);
            KeplerianEstimatedResults(:,1) = KeplerianEstimatedResults(:,1) / 1e3;
            KeplerianEstimatedResults(:,3:end) = rad2deg(KeplerianEstimatedResults(:,3:end));
            fclose(fileID);
            
            % %...Load rotational motion
            % filename = ['/Users/Michele/GitHub/tudat/tudatBundle/tudatApplications/Thesis/',outputFolder,'/rotationalEstimated.dat'];
            % fileID = fopen(filename,'r');
            % rotationalEstimatedResults = textscan(fileID,repmat('%f',[1,8]),'Delimiter',',','CollectOutput',true);
            % rotationalEstimatedResults = rotationalEstimatedResults{1}(:,2:end);
            % rotationalEstimatedResults(:,6:8) = rad2deg(rotationalEstimatedResults(:,5:7));
            % rotationalEstimatedResults(:,5) = 1.0 - quatnorm(rotationalEstimatedResults(:,1:4));
            % fclose(fileID);
            
            %...Clean up
            clear filename fileID
            
            %% Interpolate Results to Match Times
            
            %...Set interpolation time
            interpolatedTime = onboardTime(1:end-1);
            
            %...Interpolate
            if applyInterpolation
                %...Interpolate propagation results
                CartesianPropagatedResults = interp1( simulationTime, CartesianPropagatedResults, interpolatedTime, 'linear' );
                KeplerianPropagatedResults = interp1( simulationTime, KeplerianPropagatedResults, interpolatedTime, 'linear' );
                
                %...Interpolate estimation results
                CartesianEstimatedResults = interp1( onboardTime, CartesianEstimatedResults, interpolatedTime, 'linear' );
                KeplerianEstimatedResults = interp1( onboardTime, KeplerianEstimatedResults, interpolatedTime, 'linear' );
            end
            
            %% RMS Error
            
            %...Compute RMS error in position and velocity
            rmsPositionError = rms( sqrt(sum(CartesianEstimatedResults(:,1:3).^2,2)) - ...
                sqrt(sum(CartesianPropagatedResults(:,1:3).^2,2)) ) * 1e3;
            rmsVelocityError = rms( sqrt(sum(CartesianEstimatedResults(:,4:6).^2,2)) - ...
                sqrt(sum(CartesianPropagatedResults(:,4:6).^2,2)) );
            
            %...Show errors
            currentValue = {[num2str(p),'-',num2str(v),'-',num2str(m)]};
            table(currentValue,rmsPositionError,rmsVelocityError)
            
            %...Store errors
            results(p-ps(1)+1,v-vs(1)+1,m-ms(1)+1) = rmsPositionError;
        end
    end
end

results