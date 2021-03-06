function jac_imu = imu_jac(AP,CC1,CC2,CC3,CC4,CC5,CC6,H,R,h0,rho0,s1,s2,s3)
%IMU_JAC
%    JAC_IMU = IMU_JAC(AP,CC1,CC2,CC3,CC4,CC5,CC6,H,R,H0,RHO0,S1,S2,S3)

%    This function was generated by the Symbolic Math Toolbox version 8.2.
%    02-Oct-2018 14:42:00

t2 = 1.0./H;
t3 = CC1.^2;
t4 = CC2.^2;
t5 = CC3.^2;
t6 = t3+t4+t5;
t7 = sqrt(t6);
t8 = R+h0-t7;
t9 = t2.*t8;
t10 = exp(t9);
t11 = abs(CC4);
t12 = s1+1.0;
t13 = 1.0./sqrt(t6);
t14 = sign(CC4);
t15 = abs(CC5);
t16 = s2+1.0;
t17 = sign(CC5);
t18 = abs(CC6);
t19 = s3+1.0;
t20 = sign(CC6);
jac_imu = reshape([AP.*CC1.*CC4.*rho0.*t2.*t10.*t11.*t12.*t13,AP.*CC1.*CC5.*rho0.*t2.*t10.*t13.*t15.*t16,AP.*CC1.*CC6.*rho0.*t2.*t10.*t13.*t18.*t19,AP.*CC2.*CC4.*rho0.*t2.*t10.*t11.*t12.*t13,AP.*CC2.*CC5.*rho0.*t2.*t10.*t13.*t15.*t16,AP.*CC2.*CC6.*rho0.*t2.*t10.*t13.*t18.*t19,AP.*CC3.*CC4.*rho0.*t2.*t10.*t11.*t12.*t13,AP.*CC3.*CC5.*rho0.*t2.*t10.*t13.*t15.*t16,AP.*CC3.*CC6.*rho0.*t2.*t10.*t13.*t18.*t19,-(AP.*CC4.*rho0.*t10.*t12.*(t14.^2+1.0))./t14,0.0,0.0,0.0,-(AP.*CC5.*rho0.*t10.*t16.*(t17.^2+1.0))./t17,0.0,0.0,0.0,-(AP.*CC6.*rho0.*t10.*t19.*(t20.^2+1.0))./t20,1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0,-AP.*CC4.*rho0.*t10.*t11,0.0,0.0,0.0,-AP.*CC5.*rho0.*t10.*t15,0.0,0.0,0.0,-AP.*CC6.*rho0.*t10.*t18],[3,12]);
