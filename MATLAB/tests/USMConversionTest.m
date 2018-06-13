fclose('all'); clear all; close all force; profile off; clc; format long g; rng default;
addpath tests data functions

%% Read Output Of C++ Conversion

clc;

%...Tolerance in difference
tolerance = 1e-7;

for i = 1:4
    %...Print mode
    fprintf([newline,'Mode: %d',newline],i)
    
    %...Keplerian input values
    fileName = ['/Users/Michele/GitHub/tudat/tudatBundle/tudatApplications/Test/SimulationOutput/'...
        'Conversions/kepler_input_',num2str(i-1),'.dat'];
    fileID = fopen(fileName,'r');
    keplerInput = textscan(fileID,repmat('%f ',[1,7]),'CollectOutput',true,'Delimiter',',');
    keplerInput = keplerInput{1};
    fclose(fileID);
    
    %...Keplerian output values
    fileName = ['/Users/Michele/GitHub/tudat/tudatBundle/tudatApplications/Test/SimulationOutput/'...
        'Conversions/kepler_usm7_output_',num2str(i-1),'.dat'];
    fileID = fopen(fileName,'r');
    keplerUSM7Output = textscan(fileID,repmat('%f ',[1,7]),'CollectOutput',true,'Delimiter',',');
    keplerUSM7Output = keplerUSM7Output{1};
    fclose(fileID);
    fileName = ['/Users/Michele/GitHub/tudat/tudatBundle/tudatApplications/Test/SimulationOutput/'...
        'Conversions/kepler_usm6_output_',num2str(i-1),'.dat'];
    fileID = fopen(fileName,'r');
    keplerUSM6Output = textscan(fileID,repmat('%f ',[1,7]),'CollectOutput',true,'Delimiter',',');
    keplerUSM6Output = keplerUSM6Output{1};
    fclose(fileID);
    fileName = ['/Users/Michele/GitHub/tudat/tudatBundle/tudatApplications/Test/SimulationOutput/'...
        'Conversions/kepler_usmem_output_',num2str(i-1),'.dat'];
    fileID = fopen(fileName,'r');
    keplerUSMEMOutput = textscan(fileID,repmat('%f ',[1,7]),'CollectOutput',true,'Delimiter',',');
    keplerUSMEMOutput = keplerUSMEMOutput{1};
    fclose(fileID);
    
    %...USM output values
    fileName = ['/Users/Michele/GitHub/tudat/tudatBundle/tudatApplications/Test/SimulationOutput/'...
        'Conversions/usm7_output_',num2str(i-1),'.dat'];
    fileID = fopen(fileName,'r');
    usm7Output = textscan(fileID,repmat('%f ',[1,8]),'CollectOutput',true,'Delimiter',',');
    usm7Output = usm7Output{1};
    fclose(fileID);
    fileName = ['/Users/Michele/GitHub/tudat/tudatBundle/tudatApplications/Test/SimulationOutput/'...
        'Conversions/usm6_output_',num2str(i-1),'.dat'];
    fileID = fopen(fileName,'r');
    usm6Output = textscan(fileID,repmat('%f ',[1,8]),'CollectOutput',true,'Delimiter',',');
    usm6Output = usm6Output{1};
%     usm6Output = [ usm6Output, ...
%         ( cos( 0.5 * keplerInput( :, 4 ) ) .* cos( 0.5 * sum( keplerInput( :, 5 : 7 ), 2 ) ) ) < 0 ];
    fclose(fileID);
    fileName = ['/Users/Michele/GitHub/tudat/tudatBundle/tudatApplications/Test/SimulationOutput/'...
        'Conversions/usmem_output_',num2str(i-1),'.dat'];
    fileID = fopen(fileName,'r');
    usmemOutput = textscan(fileID,repmat('%f ',[1,8]),'CollectOutput',true,'Delimiter',',');
    usmemOutput = usmemOutput{1};
    fclose(fileID);
    
    %...Compare Keplerian input and output
    differenceUSM7 = abs( keplerUSM7Output - keplerInput );
    locLargeDifferenceUSM7 = any( ( ( differenceUSM7 > tolerance ) & ( abs( differenceUSM7 - 2 * pi ) > tolerance ) ) ...
        | isnan( keplerUSM7Output ), 2 );
    differenceUSM6 = abs( keplerUSM6Output - keplerInput );
    locLargeDifferenceUSM6 = any( ( ( differenceUSM6 > tolerance ) & ( abs( differenceUSM6 - 2 * pi ) > tolerance ) ) ...
        | isnan( keplerUSM6Output ), 2 );
    differenceUSMEM = abs( keplerUSMEMOutput - keplerInput );
    locLargeDifferenceUSMEM = any( ( ( differenceUSMEM > tolerance ) & ( abs( differenceUSMEM - 2 * pi ) > tolerance ) ) ...
        | isnan( keplerUSMEMOutput ), 2 );
    
    %...Percentage of wrong conversions
    failureRateUSM7 = length( find( locLargeDifferenceUSM7 ) ) / size( keplerInput, 1 ) * 100
    failureRateUSM6 = length( find( locLargeDifferenceUSM6 ) ) / size( keplerInput, 1 ) * 100
    failureRateUSMEM = length( find( locLargeDifferenceUSMEM ) ) / size( keplerInput, 1 ) * 100
    
    %...NaN check
    nansUSM7 = any( any( isnan( usm7Output ), 2 ) )
    nansUSM6 = any( any( isnan( usm6Output ), 2 ) )
    nansUSMEM = any( any( isnan( usmemOutput ), 2 ) )
    
    %...Maximum difference in SMA
    maxDiffSMAUSM7 = max( differenceUSM7( :, 2 ) )
    maxDiffSMAUSM6 = max( differenceUSM6( :, 2 ) )
    maxDiffSMAUSMEM = max( differenceUSMEM( :, 2 ) )
    
    %...Only keep elements where large difference
    keplerInput( :, end + 1 ) = wrapTo2Pi( sum( keplerInput( :, 5 : 7 ), 2 ) ); % add lambda
    
    keplerUSM7Input = keplerInput( locLargeDifferenceUSM7, : );
    keplerUSM7Output = keplerUSM7Output( locLargeDifferenceUSM7, : );
    usm7Output = usm7Output( locLargeDifferenceUSM7, : );
    differenceUSM7 = differenceUSM7( locLargeDifferenceUSM7, : );
    differenceUSM7( differenceUSM7 < 1e-10 ) = 0;
    
    keplerUSM6Input = keplerInput( locLargeDifferenceUSM6, : );
    keplerUSM6Output = keplerUSM6Output( locLargeDifferenceUSM6, : );
    usm6Output = usm6Output( locLargeDifferenceUSM6, : );
    differenceUSM6 = differenceUSM6( locLargeDifferenceUSM6, : );
    differenceUSM6( differenceUSM6 < 1e-10 ) = 0;
    
    keplerUSMEMInput = keplerInput( locLargeDifferenceUSMEM, : );
    keplerUSMEMOutput = keplerUSMEMOutput( locLargeDifferenceUSMEM, : );
    usmemOutput = usmemOutput( locLargeDifferenceUSMEM, : );
    differenceUSMEM = differenceUSMEM( locLargeDifferenceUSMEM, : );
    differenceUSMEM( differenceUSMEM < 1e-10 ) = 0;
    
    %...Convert to degrees
    keplerInput( :, 4 : end ) = rad2deg( keplerInput( :, 4 : end ) );
    
    keplerUSM7Input( :, 4 : end ) = rad2deg( keplerUSM7Input( :, 4 : end ) );
    keplerUSM7Output( :, 4 : end ) = rad2deg( keplerUSM7Output( :, 4 : end ) );
    differenceUSM7( :, 4 : end ) = rad2deg( differenceUSM7( :, 4 : end ) );
    
    keplerUSM6Input( :, 4 : end ) = rad2deg( keplerUSM6Input( :, 4 : end ) );
    keplerUSM6Output( :, 4 : end ) = rad2deg( keplerUSM6Output( :, 4 : end ) );
    differenceUSM6( :, 4 : end ) = rad2deg( differenceUSM6( :, 4 : end ) );
    
    keplerUSMEMInput( :, 4 : end ) = rad2deg( keplerUSMEMInput( :, 4 : end ) );
    keplerUSMEMOutput( :, 4 : end ) = rad2deg( keplerUSMEMOutput( :, 4 : end ) );
    differenceUSMEM( :, 4 : end ) = rad2deg( differenceUSMEM( :, 4 : end ) );
end

%%

clc;

centralBodyGravitationalParameter = 1.32712440018e20;

mode = 4;
switch mode
    case 1
        i = deg2rad(50); O = deg2rad(15); o = deg2rad(350); TA = deg2rad(10);
    case 2
        i = deg2rad(170); O = deg2rad(15); o = deg2rad(350); TA = deg2rad(10);
    case 3
        i = deg2rad(170); O = deg2rad(15); o = deg2rad(350); TA = deg2rad(170);
    case 4
        i = pi; O = deg2rad(15); o = deg2rad(350); TA = deg2rad(170);
    case 5
        i = deg2rad(0); O = deg2rad(0); o = deg2rad(0); TA = deg2rad(170);
    case 6
        i = deg2rad(0); O = deg2rad(0); o = deg2rad(15); TA = deg2rad(345);
end

u = o + TA;

quat = [ sin( i / 2 ) * cos( ( O - u ) / 2 ); sin( i / 2 ) * sin( ( O - u ) / 2 ); 
    cos( i / 2 ) * sin( ( O + u ) / 2 ); cos( i / 2 ) * cos( ( O + u ) / 2 ) ]'

% denominator = 1.0 + cos( 0.5 * i ) * cos( 0.5 * ( O + u ) )
% shadowFlag = false;
% if abs( denominator ) < 1e-5
%     shadowFlag = true
% else
%     normMrp = ( 1.0 + power( cos( 0.5 * i ), 2 ) * ( power( sin( 0.5 * ( O + u ) ), 2 ) - 1.0 ) ) / denominator / denominator
%     if normMrp > 1.0
%         shadowFlag = true
%     end
% end

% if shadowFlag
%     mrp = quat( 1 : 3 ) / ( quat( 4 ) - 1 )
%     
%     den = (norm(mrp)^2-1)^2 + 4*mrp(3)^2
%     sinl = 4 * mrp(3) * (1-norm(mrp)^2)
%     cosl = (norm(mrp)^2-1)^2 - 4*mrp(3)^2
% else
%     mrp = quat( 1 : 3 ) / ( 1 + quat( 4 ) )
%     
%     den = (1-norm(mrp)^2)^2 + 4*mrp(3)^2
%     sinl = 4 * mrp(3) * (1-norm(mrp)^2)
%     cosl = (1-norm(mrp)^2)^2 - 4*mrp(3)^2
% end
% lambda = atan2d( sinl/den, cosl/den )
% 
% arccosineArgument = ( 4 * mrp(3)^2 - mrp(1)^2 - mrp(2)^2 + ( 1-norm(mrp)^2 )^2 ) / ( 1+norm(mrp)^2 )^2;
% if ( ( abs( arccosineArgument ) - 1.0 ) > 0.0 )
%     if arccosineArgument > 0.0
%         arccosineArgument = 1.0
%     else
%         arccosineArgument = - 1.0
%     end
% end
% incl = acos( arccosineArgument )

%% Convert Quaternions to Exponential Map

clc;

%...Quaternions
% quat = [-0.419002703925548,-0.0551627524676706,-0.118296904421275,-0.898554198280556;
%     -0.987672114350896,-0.130029500651719,-0.0113761072309622,-0.0864101132863834;
%     -0.299561523151596,0.95008776981561,-0.0870727897926938,-0.00380168010402369;
%     -0.300705799504273,0.953716950748227,-6.11740603377039e-17,-2.67091715588637e-18;
%     0,0,0.996194698091746,0.0871557427476581]';
% quat = [1 5e-05 -1.08793e-11  3.39399e-13]';
quat = [0.688239  0.670883  0.275833 0.0126141]';
quat = vertcat(quat(2:4,:),quat(1,:));

%...Check norm
sum( quat.^2, 1 )

%...Convert
exp = quaternion2exponentialMap( quat )'
sqrt( sum( exp.^2, 2 ) )

mrp = quaternion2modifiedRodriguesParameter( quat )'
sqrt( sum( mrp.^2, 2 ) )

%% Convert Exponential Map to Quaternions

clc;

%...Exponential map
exp = [0.266321,-0.0350618,0.780813];
exp = [1.70094  -3.59025 -0.453745];
exp = [-1.51134      1.15969      1.75999];

%...Convert
quat = exponentialMap2quaternion( exp )'
quat.^2
dcm_quat = quat2dcm([quat(4),quat(1:3)])
trace(dcm_quat)

pos = [-3.10872e+10 -1.79482e+10 -1.33967e+11]';
vel = [-26301.8 -15185.3  10326.2]';
h = cross(pos,vel)

den = quat(3)^2 + quat(4)^2
cosl = (- quat(3)^2 + quat(4)^2)/den
sinl = 2*quat(3)*quat(4)/den
lambda = atan2d( sinl, cosl )

ve1 = norm( dot(pos,vel) / norm(pos)^2 * pos )
ve2 = norm( vel - ( dot(pos,vel) / norm(pos)^2 * pos ) )

C = 1.32712440018e20 / norm( h )

if ( ( lambda > 0 ) || ( lambda < 180 ) )
    Rf1 = ve1 * sinl + ( ve2 - C ) * cosl
    Rf2 = ve1 * cosl - ( ve2 - C ) * sinl
else
    Rf1 = ve1 * cosl - ( ve2 - C ) * sinl
    Rf2 = ve1 * sinl + ( ve2 - C ) * cosl
end

[ - 0.1 * C * sind( 210 + 330 ), 0.1 * C * cosd( 210 + 330 )]

%%

%...Function handle
skew = @(x) [0,-x(3),x(2);x(3),0,-x(1);-x(2),x(1),0];

%...Symbols
O = sym('O','real');
u = sym('u','real');
i = sym('i','real');
l = sym('l','real');

% l = O+u;

%...First rotation
n = [cos(O);sin(O);0];
dcm1 = eye(3)*cos(i)+(1-cos(i))*(n*n')-skew(n)*sin(i);

%...Second rotation
dcm2 = [cos(l),sin(l),0;-sin(l),cos(l),0;0,0,1];

% %...Quaternion
% quat = [sin(i/2)*cos((O-u)/2);sin(i/2)*sin((O-u)/2);cos(i/2)*sin(l/2);cos(i/2)*cos(l/2)];
%
% %...Modified Rodrigues Parameters
% mrp = quat(1:3)/(1+quat(4));

%...Full rotation
dcm = simplify(dcm2*dcm1);
% dcm_quat = simplify(eye(3)*(quat(4)^2-quat(1:3)'*quat(1:3)) + 2*(quat(1:3)*quat(1:3)') - 2*quat(4)*skew(quat(1:3)));
%
% %...Eigenaxis rotation (quaternion)
% eta = simplify(1/2*sqrt(trace(dcm)+1));
% epsilon = simplify(1/4/eta*simplify([dcm(2,3)-dcm(3,2);dcm(3,1)-dcm(1,3);dcm(1,2)-dcm(2,1)]));
% quat_dcm = [epsilon;eta];

% xi_quat = simplify(2*acos(quat_dcm(4)));
% e_quat = simplify(1/(sin(xi_quat/2))*quat_dcm(1:3));
% em_quat = simplify(xi_quat*e_quat);

%...Eigenaxis rotation (DCM)
Xi = simplify(trace(dcm)-1);
xi = simplify(acos(1/2*Xi));
e = simplify(1/2/sin(xi)*simplify([dcm(2,3)-dcm(3,2);dcm(3,1)-dcm(1,3);dcm(1,2)-dcm(2,1)]));

%...Exponential map
em = simplify(xi*e);

%% Supporting Functions

function e = quaternion2exponentialMap(q)
    %...Find where eta is close to 1
    locSing = abs(abs(q(4,:)) - 1) <= 1e-6;

    %...Compute exponential map, ref: https://en.wikipedia.org/wiki/Axis-angle_representation
    e = zeros(3,size(q,2));
    xi(~locSing) = 2 * atan2( sqrt( sum( q(1:3,~locSing).^2, 1 ) ), q(4,~locSing) );
    e(:,~locSing) = xi .* q(1:3,~locSing) ./ sin( xi / 2 );
    if any(locSing)
        e(:,locSing) = zeros(3,1);
    end
end

function q = exponentialMap2quaternion(e)
    %...Compute eigenaxis
    xi = norm( e );
    eig = e / xi;
    
    %...Compute quaternion
    q = zeros( 4, 1 );
    q( 1 : 3 ) = eig * sin( xi / 2 );
    q( 4 ) = cos( xi / 2 );
end

function s = quaternion2modifiedRodriguesParameter( q )
    %...Find exponential map
    e = quaternion2exponentialMap(q);
    xi = sqrt( sum( e.^2, 1 ) );
    
    %...Find MRP
    shadow = xi > 3/2*pi;
    s = zeros( 3, length( shadow ) );
    s( : , ~shadow ) = q( 1 : 3, ~shadow ) ./ ( 1 + q( 4, ~shadow ) );
    s( : , shadow ) = - q( 1 : 3, shadow ) ./ ( 1 - q( 4, shadow ) );
end