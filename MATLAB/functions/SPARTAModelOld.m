function [MROExtent,points,connectivity,solarPanelElements] = SPARTAModelOld(MRODataFile,rotateMRO,showFigure,saveFigure)
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

%...Bus
x_bus = [0;0;cx_sa-ct_sa;cx_sa+ct_sa;1;1;cx_sa+ct_sa;cx_sa-ct_sa;...
    0;0;cx_sa-ct_sa;cx_sa+ct_sa;1;1;cx_sa+ct_sa;cx_sa-ct_sa]*L_bus;
y_bus = -[0;1;1;1;1;0;0;0;0;1;1;1;1;0;0;0]*L_bus;
z_bus = [0;0;0;0;0;0;0;0;1;1;1;1;1;1;1;1]*L_bus;
faces_bus = [1,2,10,9;2,3,11,10;4,5,13,12;5,6,14,13;6,7,15,14;8,1,9,16;...
    9,10,11,16;15,12,13,14;6,5,4,7;7,4,3,8;8,3,2,1];
T_bus = vertcat(faces_bus(:,[1,2,3]),faces_bus(:,[1,3,4]));

%...Solar arrays
x_sa = [cx_sa-ct_sa;cx_sa-ct_sa;cx_sa+ct_sa;cx_sa+ct_sa;...
    cx_sa-ct_sa;cx_sa-ct_sa;cx_sa+ct_sa;cx_sa+ct_sa]*L_bus;
y_sa = -[0;-1;-1;0;0;-1;-1;0]*L_sa;
z_sa = [0;0;0;0;1;1;1;1]*L_bus;
faces_sa1 = [5,6,2,1;1,2,3,4;6,7,3,2;7,8,4,3;8,7,6,5];
T_sa1 = vertcat(faces_sa1(:,[1,2,3]),faces_sa1(:,[1,3,4]));
faces_sa2 = [5,6,2,1;1,2,3,4;8,5,1,4;7,8,4,3;8,7,6,5];
T_sa2 = vertcat(faces_sa2(:,[1,2,3]),faces_sa2(:,[1,3,4]));

%...Antenna
x_ant = [cx_ant-ct_ant;cx_ant+ct_ant;cx_ant+ct_ant;cx_ant-ct_ant;...
    cx_ant-ct_ant;cx_ant+ct_ant;cx_ant+ct_ant;cx_ant-ct_ant]*L_bus;
y_ant = -[1;1;0;0;1;1;0;0]*L_bus;
z_ant = [0;0;0;0;L_ant;L_ant;L_ant;L_ant]+L_bus;
faces_ant = [1,2,6,5;2,3,7,6;4,1,5,8;3,4,8,7;5,6,7,8];
T_ant = vertcat(faces_ant(:,[1,2,3]),faces_ant(:,[1,3,4]));

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
    patch(Tri)
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