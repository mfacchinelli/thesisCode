fclose('all'); clear all; close all force; profile off; clc; format long g; rng default;

%%

% fileName = '/Users/Michele/GitHub/tudat/tudatBundle/tudatApplications/Test/SimulationOutput/translational.dat';
% fileID = fopen(fileName,'r');
% translational = textscan(fileID,'%f %f %f %f %f %f %f %f','CollectOutput',true,'Delimiter',',');
% translational = translational{1}(:,2:end);
% fclose(fileID);

fileName = '/Users/Michele/GitHub/tudat/tudatBundle/tudatApplications/Test/SimulationOutput/rotational.dat';
% fileName = '/Users/Michele/Desktop/state.dat';
fileID = fopen(fileName,'r');
rotational = textscan(fileID,'%f %f %f %f %f %f %f %f','CollectOutput',true,'Delimiter',',');
time = rotational{1}(:,1) - rotational{1}(1,1); 
rotational = rotational{1}(:,2:end);
fclose(fileID);

fileName = '/Users/Michele/GitHub/tudat/tudatBundle/tudatApplications/Test/SimulationOutput/rotationalCopy.dat';
% fileName = '/Users/Michele/Desktop/processed_state.dat';
fileID = fopen(fileName,'r');
rotationalCopy = textscan(fileID,'%f %f %f %f %f %f %f %f','CollectOutput',true,'Delimiter',',');
rotationalCopy = rotationalCopy{1}(:,2:end);
fclose(fileID);

% time = time(1:30);
% rotational = rotational(1:30,:);

figure;
for i = 1:4
    subplot(2,4,i)
    hold on
    plot(time,rotational(:,i))
    plot(time,rotationalCopy(:,i))
    hold off
    grid on
end
subplot(2,4,5)
plot(time,sqrt(sum(rotational(:,1:4).^2,2)) - 1)
grid on
for i = 6:8
    subplot(2,4,i)
    plot(time,rotational(:,i-1))
    grid on
end

% figure;
% for i = 1:4
%     subplot(2,4,i)
%     plot(time,rotationalCopy(:,i))
%     grid on
% end
% subplot(2,4,5)
% plot(time,sqrt(sum(rotationalCopy(:,1:4).^2,2)) - 1)
% grid on
% for i = 6:8
%     subplot(2,4,i)
%     plot(time,rotationalCopy(:,i-1))
%     grid on
% end

%%

% locs = zeros(length(time)-2,1);
% for i = 1:4
%     deriv = diff(rotational(:,i))./diff(time);
%     diffDeriv = diff(deriv);
%     threshold = std(diffDeriv);
%     locs = locs | abs(diffDeriv) > threshold;
% end
% locs = find(locs)
% 
% offset = 1;
% for i = 1:length(locs)
%     if i ~= length(locs)
%         if ( locs( i + 1 ) - locs( i ) ) > 3
%             for j = locs(i):length(deriv)
%                 rotational(j+offset,1:4) = -rotational(j+offset,1:4);
%             end
%         end
%     else
%         for j = locs(i):length(deriv)
%             rotational(j+offset,1:4) = -rotational(j+offset,1:4);
%         end
%     end
% end

figure;
for i = 1:4
    subplot(2,2,i)
    plot(time(1:end-1),diff(rotational(:,i))./diff(time))
    grid on
end

% figure;
% for i = 1:4
%     subplot(2,4,i)
%     plot(time,rotational(:,i))
%     grid on
% end
% subplot(2,4,5)
% plot(time,sqrt(sum(rotational(:,1:4).^2,2)) - 1)
% grid on
% for i = 6:8
%     subplot(2,4,i)
%     plot(time,rotational(:,i-1))
%     grid on
% end

%%

fileName = '/Users/Michele/Desktop/state_0.dat';
fileID = fopen(fileName,'r');
case1 = textscan(fileID,'%f %f %f %f %f %f %f %f','CollectOutput',true,'Delimiter',',');
case1 = case1{1}; time1 = case1(:,1) - case1(1,1); case1 = case1(:,2:end);
fclose(fileID);

fileName = '/Users/Michele/Desktop/state_1.dat';
fileID = fopen(fileName,'r');
case2 = textscan(fileID,'%f %f %f %f %f %f %f %f','CollectOutput',true,'Delimiter',',');
case2 = case2{1}; time2 = case2(:,1) - case2(1,1); case2 = case2(:,2:end);
fclose(fileID);

fileName = '/Users/Michele/Desktop/state_2.dat';
fileID = fopen(fileName,'r');
case3 = textscan(fileID,'%f %f %f %f %f %f %f %f','CollectOutput',true,'Delimiter',',');
case3 = case3{1}; time3 = case3(:,1) - case3(1,1); case3 = case3(:,2:end);
fclose(fileID);

% threshold = std(diff(case3(:,1))./diff(time3));
% locs = abs(diff(diff(case3(:,1))./diff(time3))) > threshold;
% locs = find(locs);
% locs(find(diff(locs)<5)+1) = [];
% for i = locs'
% %     case3(i+2:end,1:4) = -case3(i+2:end,1:4);
% end
% locs

figure;
for i = 1:4
    subplot(2,2,i)
    plot(time3(1:end-1),diff(case3(:,i))./diff(time3))
    grid on
end

%%

% figure;
% for i = 1:4
%     subplot(2,4,i)
%     plot(time1,case1(:,i))
%     grid on
% end
% subplot(2,4,5)
% plot(time1,sqrt(sum(case1(:,1:4).^2,2)) - 1)
% grid on
% for i = 6:8
%     subplot(2,4,i)
%     plot(time1,case1(:,i-1))
%     grid on
% end
% 
% figure;
% for i = 1:4
%     subplot(2,4,i)
%     plot(time2,case2(:,i))
%     grid on
% end
% subplot(2,4,5)
% plot(time2,sqrt(sum(case2(:,1:4).^2,2)) - 1)
% grid on
% for i = 6:8
%     subplot(2,4,i)
%     plot(time2,case2(:,i-1))
%     grid on
% end

figure;
for i = 1:4
    subplot(2,4,i)
    plot(time3,case3(:,i))
    grid on
end
subplot(2,4,5)
plot(time3,sqrt(sum(case3(:,1:4).^2,2)) - 1)
grid on
for i = 6:8
    subplot(2,4,i)
    plot(time3,case3(:,i-1))
    grid on
end
% figure;
% for i = 1:3
%     subplot(2,4,i)
%     plot(time3,case3(:,i))
%     grid on
% end
% subplot(2,4,4)
% plot(time3,sqrt(sum(case3(:,1:3).^2,2)))
% grid on
% for i = 5:7
%     subplot(2,4,i)
%     plot(time3,case3(:,i-1))
%     grid on
% end

%% 

time = 0:10:864000.0;
int1 = interp1(time1,case1,time,'spline');
int2 = interp1(time2,case2,time,'spline');
int3 = interp1(time3,case3,time,'spline');

% int1 = diff(int1)./diff(time)';
% int2 = diff(int2)./diff(time)';
% int3 = diff(int3)./diff(time)';
% time = time(1:end-1);

figure;
for i = 1:4
    subplot(1,4,i)
    hold on
    plot( time, int1(:,i) )
    plot( time, int2(:,i) )
    plot( time, int3(:,i) )
    plot( [260640,260640], ylim )
    hold off
    legend('Quat','MRP','EXP')
    grid on
%     set(gca,'YScale','log')
end

meanMotion = sqrt(4.282e13/9376.0e3^3);
eulerFrequency = ( 0.5024 - 0.4265 ) / 0.4265 * meanMotion;
expected = [ (0.1 * meanMotion * cos( eulerFrequency * time ) )' , ...
             (0.1 * meanMotion * sin( eulerFrequency * time ) )' , ...
              meanMotion * ones(length(time),1) ];

figure;
for i = 1:3
    subplot(1,3,i)
    hold on
    plot( time, abs( int1(:,4+i) - expected(:,i) ) )
    plot( time, abs( int2(:,4+i) - expected(:,i) ) )
    plot( time, abs( int3(:,4+i) - expected(:,i) ) )
    hold off
    legend('Quat','MRP','EXP')
    grid on
    set(gca,'YScale','log')
end