%  MATLAB Function < ma2eta >
%
%  Author:      Michele Facchinelli
%  Date:        10th April 2017
%  Purpose:     compute eccentric and true anomalies from mean anomaly
%  Input:
%   - MA:       mean anomaly
%   - e:        eccentricity
% Output:
%   - EA:       eccentric anomaly
%   - TA:       true anomaly

function [EA,TA] = ma2eta(MA,e)

%...Initial guess
EA = MA;
EA0 = zeros(size(MA));

%...Iterate until convergence to find eccentric anomaly
iter = 0;
while any(abs(EA-EA0)>1e-10) && iter <= 10 % iterative process (with safety break)
    iter = iter+1;
    EA0 = EA;
    EA = EA0 + (MA-EA0+e.*sin(EA0))./(1-e.*cos(EA0));       % [rad]   eccentric anomaly
end

%...Find true anomaly
if nargout == 2
    TA = mod(2*atan(sqrt((1+e)./(1-e)).*tan(EA/2)),2*pi);   % [rad]   true anomaly
end