function [sphereExtent,points,triangles] = ReadSPARTASphere(sphereDataFile,showFigure)
%% Get Data

fileID = fopen(sphereDataFile,'r');
points = textscan(fileID,'%*d %f %f %f','HeaderLines',7,'CollectOutput',true); points = points{1};
triangles = textscan(fileID,'%*d %d %d %d','HeaderLines',2,'CollectOutput',true); triangles = triangles{1};
fclose(fileID);

%...Parameters for SPARTA input file
sphereExtent = zeros(3,2);
sphereExtent(:,1) = floor(min(points)-1.5);
sphereExtent(:,2) = ceil(max(points)+1.5);

%% Triangulation

%...Get surface object
T.vertices = points;
T.faces = triangles;
T.facevertexcdata = T.faces(:,3);

%...Patch to surface
if showFigure
    figure;
    patch(T)
    shading faceted
    axis equal
end

%% Create New File

% %...Extract new surface information
% points = T_red.vertices;
% points = horzcat(transpose(1:size(points,1)),points);
% connectivityList = T_red.faces;
% connectivityList = horzcat(transpose(1:size(connectivityList,1)),connectivityList);
% 
% %...Open file
% sphereInputFile = '/Users/Michele/AE Software/SPARTA/examples/sphere/data.sphere';
% if exist(sphereInputFile,'file') == 2, delete(sphereInputFile), end % delete if already exists
% fileID = fopen(sphereInputFile,'a');
% 
% %...Get header
% header = sprintf('surf file from MATLAB\n\n%d points\n%d triangles\n\n',size(points,1),size(connectivityList,1));
% fprintf(fileID,header);
% 
% %...Add points
% fprintf(fileID,'Points\n\n');
% fprintf(fileID,'%d %.12f %.12f %.12f\n',points');
% 
% %...Add triangles
% fprintf(fileID,'\nTriangles\n\n');
% fprintf(fileID,'%d %d %d %d\n',connectivityList');