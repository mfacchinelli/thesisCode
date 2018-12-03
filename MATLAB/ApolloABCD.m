fclose('all'); clear all; close all force; profile off; clc; format long g; rng default;
%
% General
%
degrad = pi/180;
raddeg = 180/pi;
Rcm    = [0.0;0.0;0.1375];
Rmrc   = Rcm;
dR     = Rmrc - Rcm;
%
% Load variables
%
Tsin = 0.2;
aero = load('data/aero');
alphas = aero.simAnglesOfAttack;
betas = aero.simAnglesOfSideSlip;
Machs = aero.MachNumber;
CDtab = aero.dragCoefficients;
CLtab = aero.liftCoefficients;
CStab = aero.sideCoefficients;
Cltab = aero.xMomentCoefficients;
Cmtab = aero.yMomentCoefficients;
Cntab = aero.zMomentCoefficients;
[atab,btab,Machtab] = meshgrid(alphas,betas,Machs);
P = [2,1,3];
CDtab = permute(CDtab,P);
CLtab = permute(CLtab,P);
CStab = permute(CStab,P);
Cltab = permute(Cltab,P);
Cmtab = permute(Cmtab,P);
Cntab = permute(Cntab,P);
atmos = load('data/MCDMeanAtmosphere');
atmTab = squeeze(mean(mean(atmos.tabular_interp,1),2));
dens = [atmos.hs_interp',atmTab(:,1)];
%
% Define constant vehicle parameters.
%
m    = 1000.0;
Ixx  = 5750.0;
Iyy  = 1215.0;
Izz  = 5210.0;
Ixy  = 0.0;
Iyz  = 0.0;
Ixz  = 0.0;
dref = 2.5;
cref = dref;
bref = dref;
Sref = 37.5;
ch1  = 0.005;
ch2  = 0.05;
Txmax = 22; % 22 N is max thrust
Tymax = 22;
Tzmax = 22;
%
h       = 40e6;%115e3;%
r       = 3389526.666666667 + h;
V       = sqrt( 42828375815756.1 / r );
L       = 0.0;
rho     = interp1(dens(:,1),dens(:,2),h/1e3,'spline');
qdyn    = 0.5 * rho * V^2;
D       = qdyn * 1.9 * Sref;
Mach    = V / 180;
g       = sqrt( L^2 + D^2 ) / m;
mV      = m * V;
mV2     = mV * V;
qS      = qdyn * Sref;
qsc     = qS * cref;
qsb     = qS * bref;
gamma   = degrad * 0.0;
cgamma  = cos(gamma);
sgamma  = sin(gamma);
tgamma  = tan(gamma);
chi     = pi/2;
alpha   = degrad * -180.0;
calpha  = cos(alpha);
salpha  = sin(alpha);
beta    = 0;
sbeta   = sin(beta);
cbeta   = cos(beta);
sigma   = degrad * 0.0;
csigma  = cos(sigma);
ssigma  = sin(sigma);
%
% Commanded angular rate.
%
c1 = g/V*cgamma*ssigma;
c2 = L/mV*tgamma*ssigma;
pc = c1*salpha + c2*calpha;
qc = L/mV - g/V*cgamma*csigma;
rc = -c1*calpha + c2*salpha;
%
% Calculate aerodynamic derivatives
%
alpha1  = alpha - 1*degrad;
salpha1 = sin(alpha1);
calpha1 = cos(alpha1);
alpha2  = alpha + 1*degrad;
salpha2 = sin(alpha2);
calpha2 = cos(alpha2);
Mach1   = min(0.95*Mach,max(max(max(Machtab))));
Mach2   = min(1.05*Mach,max(max(max(Machtab))));
if ((Mach2-Mach1) < 1e-10), Mach1 = 0.95*Mach2; end
%
% Mach
%
CD1       = interp3(atab,btab,Machtab,CDtab,alpha,beta,Mach1,'linear');
CD2       = interp3(atab,btab,Machtab,CDtab,alpha,beta,Mach2,'linear');
CL1       = interp3(atab,btab,Machtab,CLtab,alpha,beta,Mach1,'linear');
CL2       = interp3(atab,btab,Machtab,CLtab,alpha,beta,Mach2,'linear');
dCDdM     = (CD2-CD1)/(Mach2-Mach1);
dCLdM     = (CL2-CL1)/(Mach2-Mach1);
%
% Alpha
%
CD1    = interp3(atab,btab,Machtab,CDtab,alpha1,beta,min(Mach,Machtab(end)),'linear');
CL1    = interp3(atab,btab,Machtab,CLtab,alpha1,beta,min(Mach,Machtab(end)),'linear');
Cm1    = interp3(atab,btab,Machtab,Cmtab,alpha1,beta,min(Mach,Machtab(end)),'linear');
%
CD2    = interp3(atab,btab,Machtab,CDtab,alpha2,beta,min(Mach,Machtab(end)),'linear');
CL2    = interp3(atab,btab,Machtab,CLtab,alpha2,beta,min(Mach,Machtab(end)),'linear');
Cm2    = interp3(atab,btab,Machtab,Cmtab,alpha2,beta,min(Mach,Machtab(end)),'linear');
%
dCDdalpha = (CD2-CD1)/(alpha2-alpha1);
dCLdalpha = (CL2-CL1)/(alpha2-alpha1);
dCmdalpha = (Cm2-Cm1)/(alpha2-alpha1);
%
% Sideslip derivatives, dCSdbeta, dCldbeta and dCndbeta.
%
beta1  = beta - 1*degrad;
sbeta1 = sin(beta1);
cbeta1 = cos(beta1);
beta2  = beta + 1*degrad;
sbeta2 = sin(beta2);
cbeta2 = cos(beta2);
%
% Beta
%
CS1    = interp3(atab,btab,Machtab,CStab,alpha,beta1,min(Mach,Machtab(end)),'linear');
Cl1    = interp3(atab,btab,Machtab,Cltab,alpha,beta1,min(Mach,Machtab(end)),'linear');
Cn1    = interp3(atab,btab,Machtab,Cntab,alpha,beta1,min(Mach,Machtab(end)),'linear');
%
CS2    = interp3(atab,btab,Machtab,CStab,alpha,beta2,min(Mach,Machtab(end)),'linear');
Cl2    = interp3(atab,btab,Machtab,Cltab,alpha,beta2,min(Mach,Machtab(end)),'linear');
Cn2    = interp3(atab,btab,Machtab,Cntab,alpha,beta2,min(Mach,Machtab(end)),'linear');
%
dCSdbeta = (CS2-CS1)/(beta2-beta1);
dCldbeta = (Cl2-Cl1)/(beta2-beta1);
dCndbeta = (Cn2-Cn1)/(beta2-beta1);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         ROTATIONAL MOTION        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Compute Amp matrix of the longitudinal state-space equations.
%
a11 =  0;
a12 =  dCmdalpha*qS*cref/Iyy;
a21 =  1.0;
a22 =  -dCLdalpha*qS/mV;
Amp = [[a11 a12]
       [a21 a22]];
%
% Prepare input lqr controller
%
dq     = 1.5*degrad;
dalpha = 1*degrad;
dp     = 1.5*degrad;
dsigma = 4*degrad;
dr     = 1.5*degrad;
dbeta  = 1*degrad;

dq     = 2*degrad;
dalpha = 1*degrad;
dp     = 2*degrad;
dsigma = 3*degrad;
dr     = 2*degrad;
dbeta  = 2*degrad;
%
% Compute gains for reaction control only.
%
Qm = [[1/dq^2 0]
      [0      1/dalpha^2]];
Bmp  = [1/Iyy;0];
Rm   = 1/(Tymax^2);
Dmp  = 0;
[K,~,E] = lqr(Amp,Bmp,Qm,Rm);
Kqy  = K(1);
Kay  = K(2);
%
% Prepare lqr input lateral controller.
%
% Compute Aml and Cml matrices of the lateral state-space equations.
%
a11 = 0.0;
a12 = 0.0;
a13 = dCldbeta*qS*bref/Ixx;
a14 = 0.0;
a21 = 0.0;
a22 = 0.0;
a23 = dCndbeta*qS*bref/Izz;
a24 = 0.0;
a31 = salpha;
a32 = -calpha;
a33 = -1/mV*dCSdbeta*qS;
a34 = -g/V*cgamma*csigma;
a41 = -calpha;
a42 = -salpha;
a43 = tgamma*csigma/mV*dCSdbeta*qS-L/mV+g/V*cgamma*csigma;
a44 = tgamma/mV*csigma*L;
Aml = [[a11 a12 a13 a14]
       [a21 a22 a23 a24]
       [a31 a32 a33 a34]
       [a41 a42 a43 a44]];
%
% Compute gains; reaction roll control only.
%
Qm = [[1/dp^2 0      0       0]
      [0      1/dr^2 0       0]
      [0      0      1/dbeta^2 0]
      [0      0      0         1/dsigma^2]];
%
% Compute gains.
%
Bml     = [1/Ixx 0; 0 1/Izz;0 0;0 0];
Rm      = [[1/(Txmax^2) 0]
           [0           1/(Tzmax^2)]];
[K,S,E] = lqr(Aml,Bml,Qm,Rm);
Kpx     = K(1,1); Krx = K(1,2); Kbx = K(1,3); Ksx = K(1,4);
Kpz     = K(2,1); Krz = K(2,2); Kbz = K(2,3); Ksz = K(2,4);
Km      = [Kpx 0   Krx 0   Kbx Ksx; ...
           0   Kqy 0   Kay 0   0; ...
           Kpz 0   Krz 0   Kbz Ksz];
%
% Complete A matrix
%
Istar     = Ixx*Izz-Ixz*Ixz;
Ip1       = (Ixx-Iyy+Izz)/Istar;
Ip2       = ((Iyy-Izz)*Izz-Ixz*Ixz)/Istar;
Amc2      =  zeros(6,6);
Amc2(1,1) =  Ip1*qc;
Amc2(1,2) =  Ip1*pc + Ip2*rc;
Amc2(1,3) =  Ip2*qc;
Amc2(1,5) =  dCldbeta*qS*bref/Ixx;
Amc2(2,1) =  -2*Ixz/Iyy*pc + (Izz-Ixx)/Iyy*qc;
Amc2(2,2) =  0.0;
Amc2(2,3) =  (Izz-Ixx)/Iyy*pc + 2*Ixz/Iyy*rc;
Amc2(2,4) =  dCmdalpha*qS*cref/Iyy;
Ir1       = ((Ixx-Iyy)*Ixx+Ixz*Ixz)/Istar;
Ir2       = ((-Ixx+Iyy-Izz)*Ixz)/Istar;
Amc2(3,1) =  Ir1*qc;
Amc2(3,2) =  Ir1*pc + Ir2*rc;
Amc2(3,3) =  Ir2*qc;
Amc2(3,5) =  dCndbeta*qS*bref/Izz;
Amc2(4,2) =  1.0;
Amc2(4,4) =  -dCLdalpha*qS/mV;
Amc2(5,1) =  salpha;
Amc2(5,3) = -calpha;
Amc2(5,5) = -dCSdbeta*qS/mV;
Amc2(5,6) = -g*cgamma*csigma/V;
Amc2(6,1) = -calpha;
Amc2(6,3) = -salpha;
Amc2(6,4) =  tgamma*ssigma*dCLdalpha*qS/mV;
Amc2(6,5) =  tgamma*csigma*dCSdbeta*qS/mV-L/mV+g*cgamma*csigma/V;
Amc2(6,6) =  tgamma*csigma*L/mV;
%
Bmc2      = zeros(6,3);
Bmc2(1,1) = Izz/Istar;
Bmc2(1,3) = Ixz/Istar;
Bmc2(2,2) = 1/Iyy;
Bmc2(3,1) = Ixz/Istar;
Bmc2(3,3) = Ixx/Istar;
%
% Compute gains.
%
Qmc2 = [[1/dp^2 0      0       0          0         0]
        [0      1/dq^2 0       0          0         0]
        [0      0      1/dr^2  0          0         0]
        [0      0      0       1/dalpha^2 0         0]
        [0      0      0       0          1/dbeta^2 0]
        [0      0      0       0          0         1/dsigma^2]];
Rmc2 = [[1/(Txmax^2) 0           0]
        [0           1/(Tymax^2) 0]
        [0           0           1/(Tzmax^2)]];
[Kmc,S,E] = lqr(Amc2,Bmc2,Qmc2,Rmc2);
fprintf([repmat('%f, ',[1,5]),'%f,\n'],Kmc)
%
% Concatenate and digitise reference model rotational motion
%
nmstate  = 6;
nrcont   = 3;
nrout    = 6;
Amc      = zeros(nmstate,nmstate);
Bmc      = zeros(nmstate,nrcont);
Amc(1,1) = Aml(1,1);
Amc(1,3) = Aml(1,2);
Amc(1,5) = Aml(1,3);
Amc(1,6) = Aml(1,4);
Amc(2,2) = Amp(1,1);
Amc(2,4) = Amp(1,2);
Amc(3,1) = Aml(2,1);
Amc(3,3) = Aml(2,2);
Amc(3,5) = Aml(2,3);
Amc(3,6) = Aml(2,4);
Amc(4,2) = Amp(2,1);
Amc(4,4) = Amp(2,2);
Amc(5,1) = Aml(3,1);
Amc(5,3) = Aml(3,2);
Amc(5,5) = Aml(3,3);
Amc(5,6) = Aml(3,4);
Amc(6,1) = Aml(4,1);
Amc(6,3) = Aml(4,2);
Amc(6,5) = Aml(4,3);
Amc(6,6) = Aml(4,4);

Bmc(1,1) = 1/Ixx;
Bmc(2,2) = 1/Iyy;
Bmc(3,3) = 1/Izz;
Cmc    = eye(nrout);
Dmc    = zeros(nrout,nrcont);
sysm   = ss(Amc,Bmc,Cmc,Dmc);
sysd   = c2d(sysm,Tsin,'zoh');
Amc    = sysd.a;
Bmc    = sysd.b;
Cmc    = sysd.c;
