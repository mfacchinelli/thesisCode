function [tp,DV,Da,DP] = PTE(time,aero,cart,kepl)

%...Constants
mu = 4.282e13;      % Mars gravitational parameter
Rm = 3.390e6;       % Mars radius
atm_int = 250e3;	% atmospheric interface altitude

%...Reduce data
h = sqrt(sum(cart.^2,2)) - Rm;
loc = h < atm_int;
time = time(loc);
aero = aero(loc);
kepl = kepl(loc,:);

%% Periapse Time

%...Start plotting
figure;
hold on
plot(time,aero,'LineWidth',1.5)

%...First iteration
a0 = find(max(aero)==aero)-100;
b0 = find(max(aero)==aero)+100;

%...Iterative method
c = areaBisection(time,aero,a0,b0);
tp = time(c);

%...Finish plotting
plot([time(c),time(c)],ylim,'LineWidth',2,'LineStyle','--')
hold off
xlabel('Time [min]')
ylabel('Acceleration Norm [m/s^2]')
legend('IMU','Barycenter')
xlim(fix([time(1),time(end)]))
set(gca,'FontSize',25)
grid on

%% Period Change

%...Find Delta V
DV = -trapz(time,aero); % time is in minutes

%...Find Delta a
a0 = kepl(1,1); n0 = sqrt(mu/a0^3);
e0 = kepl(1,2);
Da = 2/n0*sqrt((1+e0)/(1-e0))*DV;

%...Find Delta Period
kappa = 0.955;
DP = 2*pi*kappa*((a0+Da)^(3/2)-a0^(3/2))/sqrt(mu);