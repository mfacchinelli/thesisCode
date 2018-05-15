fclose('all'); clear all; close all force; profile off; clc; format long g; rng default;
addpath functions tests

%...Set system path
setenv('PATH',[getenv('PATH'),':/Users/Michele/Software/ImageMagick-7.0.7/bin:/opt/local/bin:',...
    '/opt/local/sbin:/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin:/Library/TeX/texbin:/opt/X11/bin']);

%...Set ImageMagick path
setenv('MAGICK_HOME','/Users/Michele/Software/ImageMagick-7.0.7')
setenv('DYLD_LIBRARY_PATH','$MAGICK_HOME/ImageMagick-7.0.7/lib/')

%% SPARTA and Modeling Variables

%...Figures and tables setting
showFigure = true;
saveFigure = false;
[figSizeLarge,~,figSizeSmall] = saveFigureSettings(saveFigure);

%...SPARTA
SPARTAExec = '/Users/Michele/AE Software/SPARTA/src/spa_mac_mpi'; % path to SPARTA executable

%...Repositories
SPARTARepository = ['/Users/Michele/Library/Mobile Documents/com~apple~CloudDocs/University/Master ',...
    'Thesis/Code/SPARTA/']; % path to SPARTA folder
MRORepository = fullfile(SPARTARepository,'gallis'); % path to MRO folder
OutputRepository = fullfile(MRORepository,'data'); % path to SPARTA results
ImageRepository = fullfile(MRORepository,'figures'); % path to SPARTA figures
VideoRepository = fullfile(MRORepository,'videos'); % path to SPARTA videos

%...Files
MROInputFile = fullfile(MRORepository,'in.mro'); % path to MRO input file
MRODataFile = fullfile(MRORepository,'data.mro'); % path to MRO model file

%...Templates
MROInputTemplate = fullfile(SPARTARepository,'SPARTAInputTemplate.txt');

%...Mars Climate Database
MCDData = fullfile(SPARTARepository,'MCDEnvir');

%% Spacecraft Model

%...Compute model based on MRO discretization
[MROExtent,points,triangles,solarPanelElements] = ...
    SPARTAModel(MRODataFile,false,false,false); % output min and max dimensions of MRO

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
simGases = {'CO2','H','O'};
gasRatios = {1.0,1.0,1.0};

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
gamma = MCD.gamma;
speedOfSound = MCD.sos;

%...Interpolate to find data for rarefied regime
[density,pressure,gamma,speedOfSound] = interpolate(altitude,simAltRarefied,density,...
    pressure,gamma,speedOfSound);
numberDensity = repmat(interpolate(altitude,simAltRarefied,MCD.numberDensity),[length(simGases),1]); % overwrite number density
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

%...Ask user to proceed
selection = questdlg('Proceeding with the analysis will overwrite all previous data.',...
    'Warning','Proceed','Stop','Stop');

%...Folder names for SPARTA output
dataFolder = cell(length(simGases),length(simAltRarefied));
imageFolder = cell(length(simGases),length(simAltRarefied));

%...Switch based on user input
switch lower(selection)
    case 'proceed'
        %...Start timer
        tic;
        
        %...Clean up folders
        system(['rm ',fullfile(adapt2UNIX(OutputRepository),'*/*/*.coeff.*')]);
        system(['rm ',fullfile(adapt2UNIX(ImageRepository),'*/*/*.image.*')]);
        
        %...Generate command
        commandString = ['cd ',adapt2UNIX(MRORepository),';',...
            'mpirun -np 2 ',adapt2UNIX(SPARTAExec),' -in ',...
            adapt2UNIX(MROInputFile)];
        
        %...Loop over altitudes
        status = cell(length(simGases),length(simAltRarefied));
        outcome = cell(length(simGases),length(simAltRarefied));
        for g = 1:length(simGases)
            for h = 1:length(simAltRarefied)
                %...Add folders to list
                dataFolder{g,h} = fullfile(OutputRepository,erase(simGases{g},' '),num2str(simAltRarefied(h)));
                imageFolder{g,h} = fullfile(ImageRepository,erase(simGases{g},' '),num2str(simAltRarefied(h)));
                
                %...Create folder if inexistent
                if ~exist(dataFolder{g,h},'dir'), mkdir(dataFolder{g,h}), end
                if ~exist(imageFolder{g,h},'dir'), mkdir(imageFolder{g,h}), end
                
                %...Generate input file based on template
                fileID = fopen(MROInputFile,'w');
                fprintf(fileID,template,simBox,simGrid,numberDensity(g,h),real2sim(g,h),gasNamesString{g},...
                    gasFractions{g},gasNamesString{g},streamVelocity(h,:),gasNamesString{g},temperature,...
                    num2str(simAnglesOfAttack,'%.0f '),simTimeStep,...
                    erase(dataFolder{g,h},[MRORepository,'/']),erase(imageFolder{g,h},[MRORepository,'/']),...
                    simSteps);
                fclose(fileID);
                
                %...Run command via Terminal
                [status{g,h},outcome{g,h}] = system(commandString,'-echo');
                if status{g,h}
                    system(['open ',adapt2UNIX(fullfile(MRORepository,'log.sparta')),' -a textedit']);
                    error('SPARTA simulation failed. Opening log file.')
                end
            end
        end
        
        %...Stop timer
        toc;
    otherwise
        %...Confirm to user
        warning('SPARTA analysis skipped.')
        
        %...Add folders to list
        for g = 1:length(simGases)
            for h = 1:length(simAltRarefied)
                dataFolder{g,h} = fullfile(OutputRepository,erase(simGases{g},' '),num2str(simAltRarefied(h)));
                imageFolder{g,h} = fullfile(ImageRepository,erase(simGases{g},' '),num2str(simAltRarefied(h)));
            end
        end
end
clear commandString gasFractions fileID

%% Area Analysis

%...Compute area and normal to surface of triangles and surface areas
%   around each axis for whole spacecraft
[trianglesArea,trianglesNormal,crossSectionalArea] = computeTriangleAreaNormal(points,triangles);

%...Compute same data only for solar panels
[solarPanelTrianglesArea,solarPanelTrianglesNormal,solarPanelCrossSectionalArea] = ...
    computeTriangleAreaNormal(points,triangles(solarPanelElements.faces,:));

%% Analyze SPARTA Results

%...Get results for rarefied flow
aeroCoeffRarefied = cell(length(simGases),length(simAltRarefied),length(simAnglesOfAttack));
for g = 1:length(simGases)
    for h = 1:length(simAltRarefied)
        for a = 1:length(simAnglesOfAttack)
            %...Get (average) distribution
            distribution = cell(length(simGases),length(simAltRarefied));
            files = dir(fullfile(dataFolder{g,h},[num2str(simAnglesOfAttack(a)),'.coeff.*']));
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
            
            %...Compute pressure and friction coefficients
            dynamicPressure = 1/2 * density(h) * norm(streamVelocity(h,:))^2;
            pressureCoeffRarefied = (averageDistribution(:,1:3) - pressure(h) * trianglesNormal * ...
                roty(simAnglesOfAttack(a))) / dynamicPressure; % rotate normal to account for angle of attack
            frictionCoeffRarefied = averageDistribution(:,4:6) / dynamicPressure;
            
            %...Integrate to find force coefficients
            forceCoefficients = sum((pressureCoeffRarefied + frictionCoeffRarefied) .* ...
                trianglesArea) ./ crossSectionalArea(1);
            
            %...Find side, drag and lift coefficients
            %   Note that a negative sign is added since aerodynamic forces are
            %   positive in the negative direction
            aeroCoeffRarefied{g,h,a} = - forceCoefficients';
        end
    end
end
clear distribution files fileID result averageDistribution dynamicPressure pressureCoeffRarefied 
clear frictionCoeffRarefied forceCoefficients transformation

%% Validation

%...Get density
densityValues = interpolate(altitude,simAltRarefied,MCD.density);

%...Equation
dragDensity = @( d ) 1.47952 - 0.032578 * log( d );

%...Plot
if showFigure
    F = figure('rend','painters','pos',figSizeSmall);
    hold on
    for g = 1:length(simGases)
        plot(densityValues,cellfun(@(x)x(1),aeroCoeffRarefied(g,:)),'-o','LineWidth',1.25,'MarkerSize',10)
    end
    plot(densityValues,dragDensity(densityValues),'--','LineWidth',1.25)
    hold off
    xlabel('Density [kg m^{-3}]')
    ylabel('Drag Coefficient [-]')
    set(gca,'FontSize',15,'XScale','log')
    grid on
    legend(simGases{:},'Reference')
    if saveFigure, saveas(F,['../../Report/figures/valid_gallis'],'epsc'), end
end

%...Print LaTeX table
percentageOffset = ( cellfun(@(x)x(1),aeroCoeffRarefied(g,:)) - dragDensity(densityValues) ) ./ ...
    dragDensity(densityValues) * 100;
fprintf('\\num{%.0f} & \\num{%.1e} & %.2f & %.2f & %.1f \\\\ \n',...
    [simAltRarefied',densityValues',cellfun(@(x)x(1),aeroCoeffRarefied(g,:))',...
    dragDensity(densityValues)',percentageOffset']')

%% GIF Commands

if strcmpi(selection,'proceed')
    %...Open terminal
    ! open /Applications/Utilities/Terminal.app/
    
    %...Create command string
    GIFCommand = ['export MAGICK_HOME="/Users/Michele/Software/ImageMagick-7.0.7"; '...
        'export PATH="$MAGICK_HOME/bin:$PATH"; export DYLD_LIBRARY_PATH="$MAGICK_HOME/lib/"; '...
        'cd ',adapt2UNIX(MRORepository)];
    for h = simAltRarefied
        hs = num2str(h);
        GIFCommand = [GIFCommand,'; convert '];
        for a = simAnglesOfAttack
            as = num2str(a);
            GIFCommand = [GIFCommand,sprintf('%s ',fullfile('figures',hs,[as,'.image.*']))];
        end
        GIFCommand = [GIFCommand,fullfile('videos',['movie_',hs,'.gif'])];
    end
    
    %...Add command string to clipboard
    clipboard('copy',GIFCommand)
end

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