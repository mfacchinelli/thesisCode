%  Author:     Michele Facchinelli
%  Date:       3rd March 2017
%  Purpose:    The function 'acos2' computes the value of the arccosine (in radians)  
%              of the input 'cosine', based on the hemisphere function 'H', which
%              specifies the quadrant. If H is 1 or 0, then 'angle' is 
%              in the first or second quadrant. Otherwise, 'angle' is in 
%              the third or forth quadrant.
%  Inputs:
%   - cosine:  cosine of the angle to be computed (within +/- 1)
%   - H:       hemisphere function of the angle (+/- 1 or 0)
%  Outputs:
%   - angle:   angle with correct quadrant information (from 0 to 2*pi radian)

function angle = acos2(cosine,H)

%...Check number of input arguments
if nargin ~= 2
    error('Wrong number of input arguments.')
end

%...Check number of output arguments
if nargout > 1
    error('Wrong number of output arguments.')
end

%...Check for value of hemisphere function
check = (abs(H)~=ones(size(H))) + (H~=zeros(size(H)));
if ~isempty(find(check~=1,1))
    error('Hemisphere function H must be +/- 1 or 0.')
end

%...Check for value of cosine (with small adjustment)
margin = 1e-15;                                  % introduce margin for ...
cosine(cosine>1) = cosine(cosine>1)-margin;      % ... possible floating point ... 
cosine(cosine<-1) = cosine(cosine<-1)+margin;    % ... approximation errors
check = abs(cosine)<=(ones(size(cosine)));
if ~isempty(find(check~=1,1))
    error('Cosine value must be within +/- 1.')
end

%...Compute angle value
angle = wrapTo2Pi(2*pi + H.*acos(cosine));