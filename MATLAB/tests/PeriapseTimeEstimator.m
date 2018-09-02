fclose('all'); clear all; close all force; profile off; clc; format long g; rng default;
addpath functions tests

%% Settings

%...Constants
marsRadius = 3389526.666666667;
marsGravitationalParameter = 4.282e13;
marsAtmosphericInterface = 175;

%% Data

orbitNumber = 3;

%...Load translational motion
filename = ['/Users/Michele/Desktop/aero_',num2str(orbitNumber),'.dat'];
fileID = fopen(filename,'r');
accelerationResults = textscan(fileID,'%f','Delimiter',',','CollectOutput',true);
accelerationResults = accelerationResults{1};
fclose(fileID);

%...Load Keplerian translational motion
filename = ['/Users/Michele/Desktop/kepler_',num2str(orbitNumber),'.dat'];
fileID = fopen(filename,'r');
KeplerianResults = textscan(fileID,repmat('%f',[1,7]),'Delimiter',',','CollectOutput',true);
atmosphericTime = ( KeplerianResults{1}(:,1) - KeplerianResults{1}(1) );
KeplerianResults = KeplerianResults{1}(:,2:end);
KeplerianResults(:,1) = KeplerianResults(:,1) / 1e3;
KeplerianResults(:,3:end) = rad2deg(KeplerianResults(:,3:end));
fclose(fileID);

KeplerianResults(:,6) = KeplerianResults(:,6) - ( -12.5 );

%% Plot

altitude = KeplerianResults(:,1) .* ( 1 - KeplerianResults(:,2).^2 ) ./ ...
    ( 1 + KeplerianResults(:,2) .* cosd( KeplerianResults(:,6) ) ) - marsRadius / 1e3;

try
    [tp,DV,Da,DP] = PTE(atmosphericTime,KeplerianResults,accelerationResults);
    atPeak = find( atmosphericTime == tp, 1 );
catch
    atPeak = find(accelerationResults==max(accelerationResults),1);
end

beforeZero = find(KeplerianResults(:,6)<=0,1,'last');
afterZero = find(KeplerianResults(:,6)>=0,1,'first');

figure
yyaxis left
hold on
plot(atmosphericTime,accelerationResults)
plot([atmosphericTime(atPeak),atmosphericTime(atPeak)],[min(accelerationResults),max(accelerationResults)])
hold off
yyaxis right
hold on
plot(atmosphericTime,KeplerianResults(:,6))
plot([atmosphericTime(beforeZero),atmosphericTime(beforeZero)],[min(KeplerianResults(:,6)),max(KeplerianResults(:,6))])
plot([atmosphericTime(afterZero),atmosphericTime(afterZero)],[min(KeplerianResults(:,6)),max(KeplerianResults(:,6))])
hold off
grid on

%% Plot

figure
plot(accelerationResults(altitude<=marsAtmosphericInterface),altitude(altitude<=marsAtmosphericInterface))
grid on
xlabel('Acceleration [m s^{-2}]')
ylabel('Altitude [km]')
set(gca,'FontSize',15,'XScale','log')