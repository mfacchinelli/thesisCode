fclose('all'); clear all; close all force; profile off; clc; format long g; rng default;

%% SPARTA and Modeling Variables

%...Rotate by 180 degrees
rotateMRO = false;
if rotateMRO
    repository = 'rotated';
else
    repository = 'nominal';
end

%...Host file settings
useHostFile = false; % run SPARTA simulation with 4 cores (slows down computer a lot)

%...SPARTA
SPARTAExec = '/home/michele/sparta/spa_mpi'; % path to SPARTA executable

%...Repositories
AeroRepository = '/home/michele/aero'; % path to SPARTA folder
MRORepository = fullfile(AeroRepository,'mro'); % path to MRO folder
OutputRepository = fullfile(MRORepository,repository,'data'); % path to SPARTA results

%...Files
MROInputFile = fullfile(MRORepository,'in.mro'); % path to MRO input file
MRODataFile = fullfile(MRORepository,'data.mro'); % path to MRO model file

%...Templates
MROInputTemplate = fullfile(AeroRepository,'SPARTAInputTemplateServer.txt');

%...Mars Climate Database
MCDData = fullfile(AeroRepository,'MCDEnvir');

%% Spacecraft Model

%...Compute model based on MRO discretization
[MROExtent,points,triangles,solarPanelElements] = SPARTAModelServer(MRODataFile,rotateMRO);

%...Add more space in x-direction
MROExtent(1,1) = MROExtent(1,1) - 0.5;
MROExtent(1,2) = MROExtent(1,2) + 0.5;

%...Get variables for simulation
simBox = reshape(MROExtent',1,numel(MROExtent));

%% Constants

%...Martian environment
MarsGravityParameter = 4.282e13;
MarsRadius = 3.390e6;

%...Grid specifications
gridSpacing = 0.25; % grid size
simGrid = diff(MROExtent,[],2)/gridSpacing;
r2sOffset = 50; % number of simulated particles per cell

%...Simulation conditions
simAnglesOfAttack = linspace(-30,30,13); % angles of attack for simulation
%[linspace(-75,-30,4),linspace(-25,25,11),linspace(30,75,4)];
simAltRarefied = 125;%[100,125,150,200,300,500]; % altitudes for rarefied regime simulations

%...Time settings
simTimeStep = round(0.1*diff(MROExtent(1,:))/3500,5); % time step (1/10-th box traverse time)
simSteps = 1000; % simulation time (specified as number of steps)

%% Mars Climate Database Model

%...Load MCD parameters
%   Should be struct containing 3 fields:
%       - altitude:         array of altitudes
%       - density:          atmospheric density (per altitude)
%       - specificHeatRatio:            specific heat ratio (per altitude)
%       - gasStr:           cell array of gasses in mixture
%       - gasRatio:         matrix of gas ratios (per altitude)
%       - knudsenNumber:    array of Knudsen numbers (per altitude)
%       - numberDensity:    array of number density (per altitude)
%       - pressure:         atmospheric pressure (per altitude)
%       - speedOfSound:     speed of sound (per altitude)
%   Note that values are averaged over latitude, longitude and time.
MCD = load(MCDData);

%...Extract MCD parameters
altitude = MCD.altitude;
density = MCD.density;
pressure = MCD.pressure;
temperature = MCD.temperature;
gamma = MCD.specificHeatRatio;
speedOfSound = MCD.speedOfSound;
gasRatios = MCD.gasRatio';
gasNames = MCD.gasStr;
numberDensity = MCD.numberDensity;

%...Interpolate to find data for rarefied regime
[density,pressure,temperature,gamma,speedOfSound,gasRatios,numberDensity] = interpolate(altitude,simAltRarefied,...
    density,pressure,temperature,gamma,speedOfSound,gasRatios,numberDensity);

%...Create cell arrays of gas information for each altitude
frac = cell(size(simAltRarefied));
gas = cell(size(simAltRarefied));
gasList = cell(size(simAltRarefied));
for h = 1:length(simAltRarefied)
    loc = gasRatios(h,:) >= 0.05; % only consider gases with more than 5 % traces
    frac{h} = round(gasRatios(h,loc)/sum(gasRatios(h,loc)),2);
    frac{h} = round(frac{h}/sum(frac{h}),2);
    gas{h} = gasNames(loc);
    gasList{h} = join(gas{h},' '); gasList{h} = gasList{h}{1};
end
gasRatios = frac; 
% gasRatios{4} = [0.12,0.28,0.18,0.42]; % make sure it all sums up to one
gasNames = cellfun(@(x)erase(x,'_'),gas,'UniformOutput',false);
gasNamesString = cellfun(@(x)erase(x,'_'),gasList,'UniformOutput',false);
clear frac gas

%...Check if fractions sum up to one
if any( abs( cellfun(@sum,gasRatios) - 1 ) > 10 * eps )
    error('Gas fractions do not sum up to 100 %.')
end

%...Set ratio of real to simulated particles
real2sim = numberDensity / r2sOffset * gridSpacing^3;

%...Find circular velocity at altitude
MarsCircularVelocity = @(altitude) sqrt( MarsGravityParameter ./ ( MarsRadius + altitude * 1e3 ) );
streamVelocity = zeros(length(simAltRarefied),3);
streamVelocity(:,1) = - MarsCircularVelocity(simAltRarefied)';

%...Find molecular speed ratio and Mach number
molecularSpeedRatio = sqrt(sum(streamVelocity.^2,2))'./speedOfSound.*sqrt(gamma/2);
MachNumber = sqrt(sum(streamVelocity.^2,2))'./speedOfSound;

%% Generate SPARTA Input File

%...Get input template
fileID = fopen(MROInputTemplate,'r');
template = textscan(fileID,'%s','Delimiter','\n','CommentStyle','#');
template = join(template{1},'\n'); template = template{1};
fclose(fileID);

%% Run SPARTA

%...Function handle to adapt paths to UNIX
adapt2UNIX = @(x) replace(x,' ','\ ');

%...Folder names for SPARTA output
dataFolder = cell(size(simAltRarefied));

%...Start timer
tic;

%...Clean up folders
system(['rm ',fullfile(adapt2UNIX(OutputRepository),'*/*.coeff.*')]);

%...Generate command
if useHostFile
    commandString = ['cd ',adapt2UNIX(MRORepository),';',...
        'mpirun -hostfile hostfile -np 14 ',adapt2UNIX(SPARTAExec),' -in ',...
        adapt2UNIX(MROInputFile)];
else
    commandString = ['cd ',adapt2UNIX(MRORepository),';',...
        'mpirun -np 14 ',adapt2UNIX(SPARTAExec),' -in ',...
        adapt2UNIX(MROInputFile)];
end

%...Loop over altitudes
status = cell(size(simAltRarefied));
outcome = cell(size(simAltRarefied));
for h = 1:length(simAltRarefied)
    %...Add folders to list
    dataFolder{h} = fullfile(OutputRepository,num2str(simAltRarefied(h)));
    
    %...Create folder if inexistent
    if ~exist(dataFolder{h},'dir'), mkdir(dataFolder{h}), end
    
    %...Gas fraction specifications
    gasFractions = '';
    for g = 1:length(gasNames{h})
        gasFractions = horzcat(gasFractions,...
            sprintf('mixture %s frac %s\n',gasNames{h}{g},num2str(gasRatios{h}(g),'%.2f ')'));
    end
    gasFractions = gasFractions(1:end-1);
    
    %...Generate input file based on template
    fileID = fopen(MROInputFile,'w');
    fprintf(fileID,template,simBox,simGrid,numberDensity(h),real2sim(h),gasNamesString{h},...
        gasFractions,gasNamesString{h},streamVelocity(h,:),gasNamesString{h},...
        temperature(h),num2str(simAnglesOfAttack,'%.0f '),simTimeStep,...
        erase(dataFolder{h},[MRORepository,'/']),...
        simSteps);
    fclose(fileID);
    
    %...Run command via Terminal
    [status{h},outcome{h}] = system(commandString,'-echo');
    if status{h}
        system(['open ',adapt2UNIX(fullfile(MRORepository,'log.sparta')),' -a textedit']);
        error('SPARTA simulation failed. Opening log file.')
    end
end
clear commandString gasFractions fileID

%...Stop timer
toc;

%% Area Analysis

%...Compute area and normal to surface of triangles and surface areas
%   around each axis for whole spacecraft
[trianglesArea,trianglesNormal,crossSectionalArea] = computeTriangleAreaNormal(points,triangles);
if rotateMRO, referenceArea = 3; else, referenceArea = 1; end
referenceLength = 2.5;
momentArm = computeMomentArm(points,triangles);

%...Compute same data only for solar panels
[solarPanelTrianglesArea,solarPanelTrianglesNormal,solarPanelCrossSectionalArea] = ...
    computeTriangleAreaNormal(points,triangles(solarPanelElements.faces,:));
momentArmSolarPanel = computeMomentArm(points,triangles(solarPanelElements.faces,:));

%% Analyze SPARTA Rarefied Results

%...Get results for rarefied flow
pressureCoeffRarefied = cell(length(simAnglesOfAttack),length(simAltRarefied));
frictionCoeffRarefied = cell(length(simAnglesOfAttack),length(simAltRarefied));
aeroCoeffRarefied = cell(length(simAnglesOfAttack),length(simAltRarefied));
bendingMomentRarefied = cell(length(simAnglesOfAttack),length(simAltRarefied));
for h = 1:length(simAltRarefied)
    for a = 1:length(simAnglesOfAttack)
        %...Get (average) distribution
        distribution = cell(size(simAltRarefied));
        files = dir(fullfile(dataFolder{h},[num2str(simAnglesOfAttack(a)),'.coeff.*']));
        for i = 1:length(files)
            fileID = fopen(fullfile(files(i).folder,files(i).name),'r');
            result = textscan(fileID,repmat('%f ',[1,7]),'HeaderLines',9,'CollectOutput',true);
            result = sortrows(result{1},1);
            distribution{i} = result(:,2:end);
            fclose(fileID);
        end
        
        %...Take median to minimize effect of outliers
        %   Note that the first 2 files are NOT used, since the simulation
        %   is still starting up
        averageDistribution = median(cat(3,distribution{3:end}),3);
        
        %...Compute bending moment
        bendingMomentRarefied{a,h} = sum(cross(momentArmSolarPanel,(averageDistribution(solarPanelElements.faces,1:3) + ...
            averageDistribution(solarPanelElements.faces,4:6))) .* solarPanelTrianglesArea);
        
        %...Compute pressure and friction coefficients
        dynamicPressure = 1/2 * density(h) * norm(streamVelocity(h,:))^2;
        pressureCoeffRarefied{a,h} = (averageDistribution(:,1:3) - pressure(h) * trianglesNormal * ...
            roty(simAnglesOfAttack(a))) / dynamicPressure; % rotate normal to account for angle of attack
        frictionCoeffRarefied{a,h} = averageDistribution(:,4:6) / dynamicPressure;
        
        %...Integrate to find force coefficients
        forceCoefficients = sum((pressureCoeffRarefied{a,h} + frictionCoeffRarefied{a,h}) .* ...
            trianglesArea) / crossSectionalArea(referenceArea);
        
        %...Integrate to find moment coefficients
        momentCoefficients = sum(cross(momentArm,pressureCoeffRarefied{a,h} + frictionCoeffRarefied{a,h}) .* ...
            trianglesArea) / crossSectionalArea(referenceArea) / referenceLength;
        
        %...Find side, drag and lift coefficients
        %   Note that a negative sign is added since aerodynamic forces are
        %   positive in the negative direction
        aeroCoeffRarefied{a,h} = - [forceCoefficients,momentCoefficients]';
    end
end
clear distribution files fileID result averageDistribution dynamicPressure forceCoefficients momentCoefficients

%% Save Text File With Coefficients

%...File names
coefficient = {'MRODragCoefficients.txt','','MROLiftCoefficients.txt','','MROMomentCoefficients.txt',''};

%...Loop over settings
for i = 1:2:6
    %...Set file name and open
    fileName = fullfile(AeroRepository,coefficient{i});
    fileID = fopen(fileName,'w');
    
    %...Write number of independent variables
    fprintf(fileID,'%d\n',2);
    fprintf(fileID,'\n'); % separator
    
    %...Add independent variables
    fprintf(fileID,[repmat('%.10f\t ',[1,length(simAnglesOfAttack)]),'\n'],deg2rad(simAnglesOfAttack));
    fprintf(fileID,[repmat('%.3f\t ',[1,length(simAltRarefied)]),'\n'],simAltRarefied*1e3);
    fprintf(fileID,'\n'); % separator
    
    %...Add coefficients
    for j = 1:length(simAnglesOfAttack)
        fprintf(fileID,[repmat('%.6f\t ',[1,length(simAltRarefied)]),'\n'],cellfun(@(x)x(i),aeroCoeffRarefied(j,:)));
    end
    
    %...Close file
    fclose(fileID);
end

%% Close All Figures

close all;

%% Supporting Functions

function varargout = interpolate(h,hq,varargin)
    %...Loop over various input
    varargout = cell(size(varargin));
    for i = 1:length(varargin)
        %...Interpolate
        varargout{i} = interp1(h,varargin{i}',hq,'spline');
    end
end

function [trianglesArea,trianglesNormal,crossSectionalArea] = computeTriangleAreaNormal(points,triangles)
    %...Function handle for surface normal
    surfaceNormal = @(p) cross(p(2,:)-p(1,:),p(3,:)-p(1,:));

    %...Get triangle vertices
    vertices = arrayfun(@(i)points(triangles(i,:),:),1:size(triangles,1),'UniformOutput',false);

    %...Compute areas
    trianglesNormal = cellfun(@(x)surfaceNormal(x),vertices,'UniformOutput',false)'; % compute normal vector
    trianglesArea = cellfun(@(x)norm(x)/2,trianglesNormal);

    %...Compute surface normal
    trianglesNormal = cellfun(@(x)x/norm(x),trianglesNormal,'UniformOutput',false); % normalize vector
    trianglesNormal = cell2mat(trianglesNormal);
    
    %...Cross sectional area
    crossSectionalArea = 0.5 * sum(abs(trianglesNormal) .* trianglesArea);
end

function momentArm = computeMomentArm(points,triangles)
    %...Function handle for surface normal
    centroid = @(p) [sum(p(:,1)),sum(p(:,2)),sum(p(:,3))]/3;

    %...Get triangle vertices
    vertices = arrayfun(@(i)points(triangles(i,:),:),1:size(triangles,1),'UniformOutput',false);

    %...Compute surface normal
    trianglesCentroid = cellfun(@(x)centroid(x),vertices,'UniformOutput',false)'; % compute centroid
    trianglesCentroid = cell2mat(trianglesCentroid);
    
    %...Find distance to center
    momentArm = trianglesCentroid - [0,0,1.25+0.1375];
end