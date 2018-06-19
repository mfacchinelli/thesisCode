function [MROExtent,points,connectivity,solarPanelElements,antennaElement] = ...
    SPARTAModel(MRODataFile,rotateMRO,showFigure,saveFigure)
%% Verification Settings

%...Toggle verification
verification = false;

%% Satellite Properties

%...Bus dimensions
L_bus = 2.5;

%...Solar array dimensions
L_sa = 4.75;
cx_sa = 0.5; % x-location coefficient
ct_sa = 0.02; % thichness coefficient

%...Antenna dimensions
L_ant = 3;
cx_ant = 0.5; % x-location coefficient
ct_ant = 0.02; % thichness coefficient

%% Satellite Parametrization

%...Function handles to turn faces in triangles
getTriangles = @(faces) vertcat(faces(:,[1,2,5]),faces(:,[2,3,5]),faces(:,[3,4,5]),faces(:,[4,1,5]));

%...Bus
x_bus = [0;0;cx_sa-ct_sa;cx_sa+ct_sa;1;1;cx_sa+ct_sa;cx_sa-ct_sa;...
    0;(cx_sa-ct_sa)/2;1-(cx_sa-ct_sa)/2;1;1-(cx_sa-ct_sa)/2;(cx_sa-ct_sa)/2;...
    (cx_sa-ct_sa)/2;1-(cx_sa-ct_sa)/2;1-(cx_sa-ct_sa)/2;cx_sa;(cx_sa-ct_sa)/2;...
    0;0;cx_sa-ct_sa;cx_sa+ct_sa;1;1;cx_sa+ct_sa;cx_sa-ct_sa]*L_bus;
y_bus = -[0;1;1;1;1;0;0;0;...
    0.5;1;1;0.5;0;0;...
    0.5;0.5;0.5;0.5;0.5;...
    0;1;1;1;1;0;0;0]*L_bus;
z_bus = [0;0;0;0;0;0;0;0;...
    0.5;0.5;0.5;0.5;0.5;0.5;...
    1;1;0;0;0;...
    1;1;1;1;1;1;1;1]*L_bus;
faces_bus = [1,2,21,20,9;2,3,22,21,10;4,5,24,23,11;5,6,25,24,12;6,7,26,25,13;8,1,20,27,14;...
    20,21,22,27,15;26,23,24,25,16;6,5,4,7,17;7,4,3,8,18;8,3,2,1,19];
T_bus = getTriangles(faces_bus);

%...Solar arrays
x_sa = [cx_sa-ct_sa;cx_sa-ct_sa;cx_sa+ct_sa;cx_sa+ct_sa;...
    cx_sa-ct_sa;cx_sa;cx_sa+ct_sa;cx_sa;cx_sa;cx_sa;...
    cx_sa-ct_sa;cx_sa-ct_sa;cx_sa+ct_sa;cx_sa+ct_sa]*L_bus;
y_sa = -[0;-1;-1;0;...
    -0.5;-1;-0.5;0;-0.5;-0.5;...
    0;-1;-1;0]*L_sa;
z_sa = [0;0;0;0;...
    0.5;0.5;0.5;0.5;1;0;...
    1;1;1;1]*L_bus;
faces_sa1 = [11,12,2,1,5;1,2,3,4,10;12,13,3,2,6;13,14,4,3,7;14,13,12,11,9];
T_sa1 = getTriangles(faces_sa1);
faces_sa2 = [11,12,2,1,5;1,2,3,4,10;14,11,1,4,8;13,14,4,3,7;14,13,12,11,9];
T_sa2 = getTriangles(faces_sa2);

%...Antenna
x_ant = [cx_ant-ct_ant;cx_ant+ct_ant;cx_ant+ct_ant;cx_ant-ct_ant;...
    cx_ant-ct_ant;cx_ant;cx_ant+ct_ant;cx_ant;cx_ant;...
    cx_ant-ct_ant;cx_ant+ct_ant;cx_ant+ct_ant;cx_ant-ct_ant]*L_bus;
y_ant = -[1;1;0;0;...
    0.5;1;0.5;0;0.5;...
    1;1;0;0]*L_bus;
z_ant = [0;0;0;0;...
    0.5;0.5;0.5;0.5;1;...
    1;1;1;1]*L_ant+L_bus;
faces_ant = [1,2,11,10,6;2,3,12,11,7;4,1,10,13,5;3,4,13,12,8;10,11,12,13,9];
T_ant = getTriangles(faces_ant);

%..Triangulation
triangulation = @(t,p) struct('vertices',p,'faces',t,'facevertexcdata',t(:,3));
Tri_bus = triangulation(T_bus,[x_bus,y_bus,z_bus]);
Tri_sa1 = triangulation(T_sa1,[x_sa,y_sa,z_sa]);
Tri_sa2 = triangulation(T_sa2,[x_sa,y_sa-L_sa-L_bus,z_sa]);
Tri_ant = triangulation(T_ant,[x_ant,y_ant,z_ant]);

%% Full Model

%...Combine triangulations
Tri = triangulation(vertcat(Tri_bus.faces,...
                            Tri_sa1.faces + size(Tri_bus.vertices,1),...
                            Tri_sa2.faces + size(Tri_bus.vertices,1) + size(Tri_sa1.vertices,1),...
                            Tri_ant.faces + size(Tri_bus.vertices,1) + 2*size(Tri_sa1.vertices,1)),...
                    vertcat(Tri_bus.vertices,Tri_sa1.vertices,Tri_sa2.vertices,Tri_ant.vertices));

%...Remove duplicate points
Tri = reducepatch(Tri,size(Tri.faces,1));
Tri.facevertexcdata = Tri.faces(:,3);

%...Get solar panel elements
solarPanelElements.vertices = ismember(Tri.vertices,Tri_sa1.vertices,'rows');
solarPanelElements.faces = all(ismember(Tri.faces,find(solarPanelElements.vertices)),2);

%...Get antenna elements
antennaElement.vertices = ismember(Tri.vertices,Tri_ant.vertices,'rows');
antennaElement.faces = all(ismember(Tri.faces,find(antennaElement.vertices)),2);

%...Shift position to center of bus
center = L_bus/2*[1,-1,1];
Tri.vertices = Tri.vertices - center;

%...Accomodate rotation
if rotateMRO
    %...Transform points
    Tri.vertices = cell2mat(arrayfun(@(i)roty(-90)*Tri.vertices(i,:)',...
        1:size(Tri.vertices,1),'UniformOutput',false))';
end

%% Figures

if showFigure
    if verification
        figure
        patch(Tri_bus)
        shading faceted
        xlabel('x'), ylabel('y'), zlabel('z')
        title('Bus')

        figure
        hold on
        patch(Tri_sa1)
        patch(Tri_sa2)
        hold off
        shading faceted
        xlabel('x'), ylabel('y'), zlabel('z')
        title('Solar Array')

        figure
        patch(Tri_ant)
        shading faceted
        xlabel('x'), ylabel('y'), zlabel('z')
        title('Antenna')
    end
    
    F = figure;
    hold on
    quiver3( 0, L_bus/2 + L_sa - 1, L_bus + L_ant/2, 1, 0, 0, 'LineWidth',2,'MaxHeadSize',0.5 )
    quiver3( 0, L_bus/2 + L_sa - 1, L_bus + L_ant/2, 0, 1, 0, 'LineWidth',2,'MaxHeadSize',0.5 )
    quiver3( 0, L_bus/2 + L_sa - 1, L_bus + L_ant/2, 0, 0, 1, 'LineWidth',2,'MaxHeadSize',0.5 )
    text( 1, L_bus/2 + L_sa - 1, L_bus + L_ant/2, 'x','FontWeight','bold',...
        'VerticalAlignment','mid','HorizontalAlignment','right')
    text( 0, L_bus/2 + L_sa, L_bus + L_ant/2, 'y','FontWeight','bold',...
        'VerticalAlignment','mid','HorizontalAlignment','left')
    text( 0, L_bus/2 + L_sa - 1, L_bus + L_ant/2 + 1, 'z','FontWeight','bold',...
        'VerticalAlignment','bottom','HorizontalAlignment','center')
    patch(Tri)
    hold off
    shading faceted
    view([127.5,30])
    if saveFigure
        axis equal off
        saveas(F,'../../Report/figures/mro_design_sparta','epsc')
    else
        axis equal
        xlabel('x'), ylabel('y'), zlabel('z')
        title('MRO')
    end
end

%% Save File For SPARTA

%...Combine triangulation points
points = Tri.vertices;
pointsList = horzcat(transpose(1:size(points,1)),points);

%...Combine triangulation triangles
connectivity = Tri.faces;
connectivityList = horzcat(transpose(1:size(connectivity,1)),connectivity);

%...Open file
if exist(MRODataFile,'file') == 2, delete(MRODataFile), end % delete if already exists
fileID = fopen(MRODataFile,'a');

%...Get and add header
header = sprintf('# surf file from MATLAB\n\n%d points\n%d triangles\n\n',...
                 size(pointsList,1),size(connectivityList,1));
fprintf(fileID,header);

%...Add points
fprintf(fileID,'Points\n\n');
fprintf(fileID,'%d %.3f %.3f %.3f\n',pointsList');

%...Add triangles
fprintf(fileID,'\nTriangles\n\n');
fprintf(fileID,'%d %d %d %d\n',connectivityList');

%...Parameters for SPARTA input file
MROExtent = zeros(3,2);
MROExtent(:,1) = floor(min(points) - 3);
MROExtent(:,2) = ceil(max(points) + 3);