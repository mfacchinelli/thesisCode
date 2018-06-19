fclose('all'); clear all; close all force; profile off; clc; format long g; rng default;

%% SPARTA and Modeling Variables

%...Select test case
%       0: nominal case (all altitudes)
%       1: both real2sim and gridSpacing
%       2: only real2sim
%       3: only gridSpacing
% testCase = 2;

%...Rotate by 180 degrees
rotateMRO = false;
if rotateMRO
    repository = 'rotated';
else
    repository = 'nominal';
end

%...Host file settings
useHostFile = false; % run SPARTA simulation with 30 cores

%...SPARTA
SPARTAExec = '/home/michele/sparta/src/spa_mpi'; % path to SPARTA executable

%...Repositories
AeroRepository = '/home/michele/aero'; % path to SPARTA folder
MRORepository = fullfile(AeroRepository,'mro'); % path to MRO folder
switch testCase
    case 0
        OutputRepository = fullfile(MRORepository,repository,'data'); % path to SPARTA results
    case 1
        OutputRepository = fullfile(MRORepository,repository,'data_both'); % path to SPARTA results
    case 2
        OutputRepository = fullfile(MRORepository,repository,'data_ratio'); % path to SPARTA results
    case 3
        OutputRepository = fullfile(MRORepository,repository,'data_grid'); % path to SPARTA results
end

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
switch testCase
    case 0
        r2sOffset = 20; % number of simulated particles per cell
        gridSpacing = 0.25; % grid size
    case 1
        r2sOffset = 35; % number of simulated particles per cell
        gridSpacing = 0.15; % grid size
    case 2
        r2sOffset = 50; % number of simulated particles per cell
        gridSpacing = 0.25; % grid size
    case 3
        r2sOffset = 17.5; % number of simulated particles per cell
        gridSpacing = 0.1; % grid size
end
simGrid = diff(MROExtent,[],2)/gridSpacing;

%...Simulation conditions
switch testCase
    case 0
        simAnglesOfAttack = [linspace(-75,-30,4),-25:5:25,linspace(30,75,4)]; % angles of attack for simulation
        simAltRarefied = [100,125,150,200,300,500]; % altitudes for rarefied regime simulations
    otherwise
        simAnglesOfAttack = -30:5:30; % angles of attack for simulation
        simAltRarefied = 125; % altitudes for rarefied regime simulations
end

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
if testCase == 0
    gasRatios{4} = [0.12,0.28,0.18,0.42]; % make sure it all sums up to one
end
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
        'nice mpirun -np 30 ',adapt2UNIX(SPARTAExec),' -in ',...
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
        erase(dataFolder{h},[MRORepository,'/']),simSteps);
    fclose(fileID);
    
    %...Run command via Terminal
    status = system(commandString,'-echo');
end
clear commandString gasFractions fileID

%...Stop timer
toc;

%% Supporting Functions

function varargout = interpolate(h,hq,varargin)
    %...Loop over various input
    varargout = cell(size(varargin));
    for i = 1:length(varargin)
        %...Interpolate
        varargout{i} = interp1(h,varargin{i}',hq,'spline');
    end
end