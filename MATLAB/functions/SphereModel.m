function [sphereExtent,points,connectivity] = SphereModel(sphereDataFile,showFigure)
%% Generate Sphere

%...Radius and discretization
L = 3;
n = 20;

%...Spherical coordinates
theta = linspace(-pi,pi,n);
phi = linspace(pi/2,-pi/2,n)';

%...Cartesian coordinates
x = L*sin(phi)*ones(1,n);
y = L*cos(phi)*cos(theta);
z = L*cos(phi)*sin(theta);

%% Triangulation

%...Generate triangulation
T = surf2patch(x,y,z,'triangles');

%...Remove duplicates
T = reducepatch(T,size(T.faces,1));
T.facevertexcdata = T.vertices(:,3);

%...Patch to surface
if showFigure
    figure;
    patch(T)
    shading faceted
    axis equal
end

%% Save To File

%...Combine triangulation points
points = T.vertices;
pointsList = horzcat(transpose(1:size(points,1)),points);

%...Combine triangulation triangles
connectivity = T.faces;
connectivityList = horzcat(transpose(1:size(connectivity,1)),connectivity);

%...Open file
if exist(sphereDataFile,'file') == 2, delete(sphereDataFile), end % delete if already exists
fileID = fopen(sphereDataFile,'a');

%...Get and add header
header = sprintf('surf file from MATLAB\n\n%d points\n%d triangles\n\n',...
                 size(pointsList,1),size(connectivityList,1));
fprintf(fileID,header);

%...Add points
fprintf(fileID,'Points\n\n');
fprintf(fileID,'%d %.3f %.3f %.3f\n',pointsList');

%...Add triangles
fprintf(fileID,'\nTriangles\n\n');
fprintf(fileID,'%d %d %d %d\n',connectivityList');

%...Parameters for SPARTA input file
sphereExtent = zeros(3,2);
sphereExtent(:,1) = floor(min(points)-1.5);
sphereExtent(:,2) = ceil(max(points)+1.5);