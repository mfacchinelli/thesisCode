function [tp,DV,Da,DP] = pte2(time,aero,cart,kepl,actual_aero)

%...Constants
mu = 4.282e13;      % Mars gravitational parameter
Rm = 3.390e6;       % Mars radius
atm_int = 500e3;	% atmospheric interface altitude

%...Reduce data
h = sqrt(sum(cart.^2,2)) - Rm;
loc = h < atm_int;
time = time(loc,:);
theta = kepl(loc,end); theta(theta>pi) = theta(theta>pi) - 2*pi; % use true anomaly as independent variable
aero = sqrt(sum(aero(loc,:).^2,2));
kepl = kepl(loc,:);
actual_aero = actual_aero(loc,:);

%% Periapse Time

%...Start plotting
figure;
hold on
plot(theta,aero,'LineWidth',1.5)

%...First iteration
a0 = find(max(aero)==aero)-100;
b0 = find(max(aero)==aero)+100;

%...Iterative method
c = areaBisection(theta,aero,a0,b0);
theta(c)
tp = interp1(theta,time,theta(c),'spline');

%...Finish plotting
plot(theta,actual_aero,'LineWidth',1.5,'LineStyle','--')
plot([theta(c),theta(c)],ylim,'LineWidth',2,'LineStyle',':')
hold off
xlabel('\theta [rad]')
ylabel('Acceleration Norm [m/s^2]')
legend('IMU','Actual','Barycenter')
grid on

%% 

a0 = kepl(1,1); n0 = sqrt(mu/a0^3);
e0 = kepl(1,2);
dM = tp * n0;
[dE,dT] = ma2eta(dM,e0);

%% Period Change

%...Find Delta V
DV = -trapz(time,aero);

%...Find Delta a
a0 = kepl(1,1); n0 = sqrt(mu/a0^3);
e0 = kepl(1,2);
Da = 2/n0*sqrt((1+e0)/(1-e0))*DV;

%...Find Delta Period
kappa = 0.955;
DP = 2*pi*kappa*((a0+Da)^(3/2)-a0^(3/2))/sqrt(mu);