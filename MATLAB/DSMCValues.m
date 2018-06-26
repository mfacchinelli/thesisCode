fclose('all'); clear all; close all force; profile off; clc; format long g; rng default;

%% Translational Temperature

clc;

k = 1.38064852e-23; % boltzmann constant
c = 3000; % velocity

m = [0.03995,0.02801,0.04401,0.00101,0.00202,0.01802,0.004,0.02801,0.016,0.032,0.048];

T_tr = m * c^2 / ( 3 * k )

%% Rotational Degrees of Freedom

clc;

T_rot = [3103.0,6159.0,3371.0,2256.0,3103.0,6159.0,3371.0,2256.0];
T_ref = 500;

Z_rot = 2 * T_rot ./ T_ref ./ ( exp( T_rot ./ T_ref ) - 1 )

%% Rotational Relaxation Number

clc;

%...Parker
Z_inf_1 = 15.7;
T_star_1 = 91.5;

%...Lordi-Mates
Z_inf_2 = 23.0;
T_star_2 = 80.0;

Z = @(Z_inf, T_star, T_tr) Z_inf ./ ( 1 + pi^(3/2)/2 * (T_star./T_tr).^(1/2) + ...
    ( pi + pi^2/4 ) * (T_star./T_tr) );
T_plot = 0:0.1:1200;

figure
hold on
plot(T_plot,Z(Z_inf_1,T_star_1,T_plot),'LineWidth',1.25)
plot(T_plot,Z(Z_inf_2,T_star_2,T_plot),'LineWidth',1.25,'LineStyle','--')
hold off
legend('Parker','Lordi-Mates','Location','Best')
set(gca,'FontSize',15)
grid on

T_ref = [3103.0,3339.0,6159.0,3371.0,2256.0];
1./Z(Z_inf_1,T_star_1,T_ref)
1./Z(Z_inf_2,T_star_2,T_ref)

%% Vibrational Collision Number

clc;

gas = 'CO2';
switch gas
    case 'N2'
        C1 = 9.1;
        C2 = 220.0;
        T_ref = 3371.0;
        omega = 0.74;
    case 'O2'
        C1 = 56.5;
        C2 = 153.5;
        T_ref = 2256.0;
        omega = 0.77;
    case 'CO'
        C1 = 37.7;
        C2 = 175.0;
        T_ref = 3103.0;
        omega = 0.73;
    case 'CO2'
        C1 = 37.7;
        C2 = 175.0;
        T_ref = 3339.0;
        omega = 0.93;
end
Z = @(T) C1 ./ T.^omega .* exp( C2 ./ power( T, 1/3 ) );

1./Z(T_ref)