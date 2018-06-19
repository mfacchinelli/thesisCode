fclose('all'); clear all; close all force; profile off; clc; format long g; rng default;

%% SPARTA and Modeling Variables

%...SPARTA
SPARTAExec = '/home/michele/sparta/src/spa_mpi'; % path to SPARTA executable

%...Repositories
AeroRepository = '/home/michele/aero'; % path to SPARTA folder
MRORepository = fullfile(AeroRepository,'gallis'); % path to MRO folder
OutputRepository = fullfile(MRORepository,'data'); % path to SPARTA results

%...Files
MROInputFile = fullfile(MRORepository,'in.mro'); % path to MRO input file
MRODataFile = fullfile(MRORepository,'data.mro'); % path to MRO model file

%...Templates
MROInputTemplate = fullfile(AeroRepository,'SPARTAInputTemplateServer.txt');

%...Mars Climate Database
MCDData = fullfile(AeroRepository,'MCDEnvir');

%% Spacecraft Model

%...Compute model based on MRO discretization
[MROExtent,points,triangles,solarPanelElements] = SPARTAModelServer(MRODataFile,false);

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
r2sOffset = 15; % number of simulated particles per cell

%...Simulation conditions
simAnglesOfAttack = 0; % angles of attack for simulation
simAltRarefied = linspace(100,150,11); % altitudes for rarefied regime simulations
simGases = {'CO2'};
gasRatios = {1.0};

%...Time settings
simTimeStep = round(0.1*diff(MROExtent(1,:))/3500,5); % time step (1/10-th box traverse time)
simSteps = 1000; % simulation time (specified as number of steps)

%% Mars Climate Database Model

%...Load MCD parameters
%   Should be struct containing 3 fields:
%       - altitude:         array of altitudes
%       - density:          atmospheric density (per altitude)
%       - gamma:            specific heat ratio (per altitude)
%       - gasStr:           cell array of gasses in mixture
%       - gasRatio:         matrix of gas ratios (per altitude)
%       - knudsenNumber:    array of Knudsen numbers (per altitude)
%       - numberDensity:    array of number density (per altitude)
%       - pressure:         atmospheric pressure (per altitude)
%       - sos:              speed of sound (per altitude)
%   Note that values are averaged over latitude, longitude and time.
MCD = load(MCDData);

%...Extract MCD parameters
altitude = MCD.altitude;
density = MCD.density;
pressure = MCD.pressure;
temperature = 150;
gamma = MCD.specificHeatRatio;
speedOfSound = MCD.speedOfSound;

%...Interpolate to find data for rarefied regime
[density,pressure,gamma,speedOfSound] = interpolate(altitude,simAltRarefied,density,...
    pressure,gamma,speedOfSound);

%...Number density
avogadrosNumber = 6.022140857e23;
molarMass = 44.01e-3;
numberDensity = avogadrosNumber / molarMass * density;
real2sim = numberDensity / r2sOffset * gridSpacing^3; % overwrite real-to-simulated-particles ratio

%...Find circular velocity at altitude
MarsCircularVelocity = @(altitude) sqrt( MarsGravityParameter ./ ( MarsRadius + altitude * 1e3 ) );
streamVelocity = zeros(length(simAltRarefied),3);
streamVelocity(:,1) = - MarsCircularVelocity(simAltRarefied)';

%...Create cell arrays of gas
gasNamesString = cell(1,length(simGases));
gasFractions = cell(1,length(simGases));
molecularSpeedRatio = zeros(length(simGases),length(simAltRarefied));
MachNumber = zeros(length(simGases),length(simAltRarefied));
for g = 1:length(simGases)
    gasNamesString{g} = simGases{g};
    
    splitGases = split(simGases{g});
    gasFraction = [];
    for i = 1:length(splitGases)
        gasFraction = horzcat(gasFraction,...
            sprintf('mixture %s frac %s\n',splitGases{i},num2str(gasRatios{g}(i),'%.2f ')'));
    end
    gasFractions{g} = gasFraction(1:end-1);
    
    %...Find molecular speed ratio and Mach number
    molecularSpeedRatio(g,:) = sqrt(sum(streamVelocity.^2,2))'./speedOfSound.*sqrt(gamma/2);
    MachNumber(g,:) = sqrt(sum(streamVelocity.^2,2))'./speedOfSound;
end

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
dataFolder = cell(length(simGases),length(simAltRarefied));
imageFolder = cell(length(simGases),length(simAltRarefied));

%...Start timer
tic;

%...Clean up folders
system(['rm ',fullfile(adapt2UNIX(OutputRepository),'*/*/*.coeff.*')]);

%...Generate command
commandString = ['cd ',adapt2UNIX(MRORepository),';',...
    'nice mpirun -np 14 ',adapt2UNIX(SPARTAExec),' -in ',...
    adapt2UNIX(MROInputFile)];

%...Loop over altitudes
status = cell(length(simGases),length(simAltRarefied));
outcome = cell(length(simGases),length(simAltRarefied));
for g = 1:length(simGases)
    for h = 1:length(simAltRarefied)
        %...Add folders to list
        dataFolder{g,h} = fullfile(OutputRepository,erase(simGases{g},' '),num2str(simAltRarefied(h)));
        
        %...Create folder if inexistent
        if ~exist(dataFolder{g,h},'dir'), mkdir(dataFolder{g,h}), end
        if ~exist(imageFolder{g,h},'dir'), mkdir(imageFolder{g,h}), end
        
        %...Generate input file based on template
        fileID = fopen(MROInputFile,'w');
        fprintf(fileID,template,simBox,simGrid,numberDensity(g,h),real2sim(g,h),gasNamesString{g},...
            gasFractions{g},gasNamesString{g},streamVelocity(h,:),gasNamesString{g},temperature,...
            num2str(simAnglesOfAttack,'%.0f '),simTimeStep,...
            erase(dataFolder{g,h},[MRORepository,'/']),simSteps);
        fclose(fileID);
        
        %...Run command via Terminal
        status{g,h} = system(commandString,'-echo');
    end
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