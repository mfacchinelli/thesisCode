fclose('all'); clear all; close all force; profile off; clc; format long g; rng default;

sparta = false;
showFigure = true;
saveFigure = true;

%% Satellite Dimensions

D_ant = 3; % antenna diameter [m]
A_sa = 12; % solar panel area [m^2]
A_tot = 37.5; % total frontal area [m^2]

A_ant = pi*D_ant^2/4; % antenna frontal area [m^2]
A_sa_ant = 2*A_sa + A_ant; % frontal area of solar panels and antenna

L_bus = sqrt(A_tot - A_sa_ant); % bus side length [m]
L_sa = A_sa/L_bus; % solar panel extension [m]

L_bus^2+2*(L_sa*L_bus)+pi*D_ant^2/4;

%% Satellite Mass Properties

%...Dimensions
L_bus = 2.5; 

L_sa = 4.75;
cx_sa = 0.5; % x-location coefficient
y_sa = 0.5*L_bus + 0.5*L_sa; %x_sa1 = 0; z_sa1 = 0;

L_ant = 3;
cx_ant = 0.5; cy_ant = 0.5; cz_ant = 0.5; % x-, y- and z-location coefficient
th1 = 2/6*pi+pi; th2 = pi/2+pi; th3 = th2-th1;
n = 20; r = 0.5*L_ant/sin(th3); h = r*(1-cos(th3));
z_ant = 0.5*L_bus + 0.5*L_ant; % modeled as disk above CM

%...Mass distribution
m_tot = 1000;
f_bus = 0.65; f_sa = 0.15; f_ant = 1-f_bus-2*f_sa;
m_bus = f_bus*m_tot;
m_sa = f_sa*m_tot;
m_ant = f_ant*m_tot;

%...Center of mass
%   d = sum(d_i*m_i)/sum(m_i)
x_cm = 0; % symmetry
y_cm = 0; % symmetry
z_cm = (m_ant*z_ant)/m_tot;

%...Moment of inertia
%   I = sum(I_cm_i + m_i*d_i^2)
I = @(I,m,d) I + m*d^2;

I_bus = 1/6*m_bus*L_bus^2; % equal for all axes
I_sa_xx = 1/12*m_sa*(L_sa^2 + L_bus^2); I_sa_yy = 1/12*m_sa*L_bus^2; I_sa_zz = 1/12*m_sa*L_sa^2;
I_ant_xx = 1/8*m_ant*L_ant^2; I_ant_yy = 1/16*m_ant*L_ant^2; I_ant_zz = 1/16*m_ant*L_ant^2;

I_xx = I(I_bus,m_bus,z_cm) + 2*I(I_sa_xx,m_sa,y_sa) + I(I_ant_xx,m_ant,z_ant-z_cm);

I_yy = I(I_bus,m_bus,z_cm) + 2*I(I_sa_yy,m_sa,0) + I(I_ant_yy,m_ant,z_ant-z_cm);

I_zz = I(I_bus,m_bus,0) + 2*I(I_sa_zz,m_sa,y_sa) + I(I_ant_zz,m_ant,0);

diag([I_xx,I_yy,I_zz])
round(diag([I_xx,I_yy,I_zz]),-1)

%% Satellite Parametrization

% L_bus = 2.5;
% L_sa = 4.75;
% L_ant = 3;

%...Bus
x_bus = [0;0;1;1;0;0;1;1]*L_bus;
y_bus = -[0;1;1;0;0;1;1;0]*L_bus;
z_bus = [0;0;0;0;1;1;1;1]*L_bus;
faces_bus = [1,2,6,5;4,3,2,1;2,3,7,6;3,4,8,7;4,1,5,8;5,6,7,8];
T_bus = vertcat(faces_bus(:,[1,2,3]),faces_bus(:,[1,3,4]));

%...Solar arrays
x_sa = [cx_sa;cx_sa;cx_sa;cx_sa]*L_bus;
y_sa = -[0;0;-1;-1]*L_sa;
z_sa = [0;1;0;1]*L_bus;
faces_sa = [3,1,2,4];
T_sa = vertcat(faces_sa(:,[1,2,3]),faces_sa(:,[1,3,4]));

%...Antenna
theta = linspace(-pi,pi,n);
phi = linspace(th1,th2,n)';

x_ant = r*sin(phi)*ones(1,n);
y_ant = -r*cos(phi)*cos(theta);
z_ant = r*cos(phi)*sin(theta);
if sparta
    [T_ant,x_ant,y_ant,z_ant] = generateAntenna(r,x_ant,y_ant,z_ant);
end
x_ant = x_ant+(r+cy_ant*L_bus)-cy_ant*r*(1-cos(th3));
y_ant = y_ant-cx_ant*L_bus;
z_ant = z_ant+(cz_ant*L_ant+L_bus);

%..Triangulation
triangulation = @(t,p) struct('vertices',p,'faces',t,'facevertexcdata',t(:,2));
Tri_bus = triangulation(T_bus,[x_bus,y_bus,z_bus]);
Tri_sa1 = triangulation(T_sa,[x_sa,y_sa,z_sa]);
Tri_sa2 = triangulation(T_sa,[x_sa,y_sa-L_sa-L_bus,z_sa]);
if sparta
    Tri_ant = triangulation(T_ant,[x_ant,y_ant,z_ant]);
else
    Tri_ant = surf2patch(x_ant,y_ant,z_ant,'triangles');
end
Tri_ant.facevertexcdata = 4*(sqrt(sum(Tri_ant.vertices.^2,2))-4);

%...Plot
if showFigure
    F = figure;
    hold on
    quiver3( L_bus/2, L_sa - 1, L_bus + L_ant/2, 1, 0, 0, 'LineWidth', 2 )
    quiver3( L_bus/2, L_sa - 1, L_bus + L_ant/2, 0, 1, 0, 'LineWidth', 2 )
    quiver3( L_bus/2, L_sa - 1, L_bus + L_ant/2, 0, 0, 1, 'LineWidth', 2 )
    patch(Tri_bus)
    patch(Tri_sa1)
    patch(Tri_sa2)
    patch(Tri_ant,'LineStyle','none')
    shading faceted
    hold off
    axis equal
    xlabel('x'); ylabel('y'); zlabel('z')
    xlim([-1,4]), ylim([-8.5,6]), zlim([-1,7])
    view([127.5,30])
    if saveFigure
        axis off
        figView = {[127.5,30],[90,0],[180,0],[90,90]};
        figName = {'iso','x','y','z'};
        for i = 1:length(figView)
            if i == 1
                legend('x','y','z','Location','South','Orientation','Horizontal')
            else
                legend('hide')
            end
            view(figView{i})
%             pause(2.5)
            saveas(F,['../../Report/figures/mro_design_',figName{i}],'epsc')
        end
    end
end

%% Other

function [T_ant,x_ant,y_ant,z_ant] = generateAntenna(r,x_ant,y_ant,z_ant)
    [SPARTApoints,SPARTAtriangles] = ReadSPARTASphere();
    SPARTApoints = r*SPARTApoints{1};
    SPARTAtriangles = SPARTAtriangles{1};

    locX = SPARTApoints(:,1) <= max(max(x_ant)) & SPARTApoints(:,1) >= min(min(x_ant));
    locY = SPARTApoints(:,2) <= max(max(y_ant)) & SPARTApoints(:,2) >= min(min(y_ant));
    locZ = SPARTApoints(:,3) <= max(max(z_ant)) & SPARTApoints(:,3) >= min(min(z_ant));
    loc_p = locX & locY & locZ;
    SPARTApoints = SPARTApoints(loc_p,:);

    loc_t = any(ismember(SPARTAtriangles,find(~loc_p)),2);
    T_ant = SPARTAtriangles(~loc_t,:);
    [~,T_ant] = ismember(T_ant,find(loc_p));

    x_ant = SPARTApoints(:,1); 
    y_ant = SPARTApoints(:,2);
    z_ant = SPARTApoints(:,3);
end