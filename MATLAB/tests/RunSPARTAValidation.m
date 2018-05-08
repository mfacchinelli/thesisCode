fclose('all'); clear all; close all force; profile off; clc; format long g; rng default;
addpath functions other

%...Set system path
setenv('PATH',[getenv('PATH'),':/Users/Michele/Software/ImageMagick-7.0.7/bin:/opt/local/bin:',...
    '/opt/local/sbin:/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin:/Library/TeX/texbin:/opt/X11/bin']);

%...Set ImageMagick path
setenv('MAGICK_HOME','/Users/Michele/Software/ImageMagick-7.0.7')
setenv('DYLD_LIBRARY_PATH','$MAGICK_HOME/ImageMagick-7.0.7/lib/')

%% SPARTA and Modeling Variables

%...Figures setting
showFigure = true;

%...SPARTA
SPARTAExec = '/Users/Michele/AE Software/SPARTA/src/spa_mac_mpi'; % path to SPARTA executable

%...Repositories
SPARTARepository = ['/Users/Michele/Library/Mobile Documents/com~apple~CloudDocs/University/Master ',...
    'Thesis/Code/SPARTA/']; % path to SPARTA folder
MRORepository = fullfile(SPARTARepository,'validation'); % path to MRO folder
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
[MROExtent,points,triangles,solarPanelElements] = SPARTAModel(MRODataFile,false,false,false); % output min and max dimensions of MRO

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
r2sOffset = 17.5; % number of simulated particles per cell

%...Simulation conditions
simAnglesOfAttack = 0; % angles of attack for simulation
simDensity = logspace(-6,-15,10); % altitudes for rarefied regime simulations

%...Time settings
simTimeStep = round(0.1*diff(MROExtent(1,:))/4500,5); % time step (1/10-th box traverse time)
simSteps = 1000; % simulation time (specified as number of steps)

%% Mars Climate Database Model

%...Create cell arrays of gas
gasRatios = [0.95,0.05];
gasNames = {'CO2','N2'};
gasNamesString = 'CO2 N2';

%...Set ratio of real to simulated particles
molarMasses = ( 44.01 * gasRatios(1) + 28.01 * gasRatios(2) ) * 1e-3;
numberDensity = 6.0221e23 / molarMasses * simDensity;
real2sim = numberDensity / (r2sOffset / gridSpacing^3);

%...Find circular velocity at altitude
streamVelocity = zeros(length(simDensity),3);
streamVelocity(simDensity<100e-9,1) = - 4811;
streamVelocity(simDensity>=100e-9,1) = - 3611;

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
dataFolder = cell(size(simDensity));
imageFolder = cell(size(simDensity));

%...Switch based on user input
switch lower(selection)
    case 'proceed'
        %...Start timer
        tic;
        
        %...Clean up folders
        system(['rm ',fullfile(adapt2UNIX(OutputRepository),'*/*.coeff.*')]);
        system(['rm ',fullfile(adapt2UNIX(ImageRepository),'*/*.image.*')]);
        system(['rm ',fullfile(adapt2UNIX(VideoRepository),'*')]);
        
        %...Generate command
        commandString = ['cd ',adapt2UNIX(MRORepository),';',...
            'mpirun -np 2 ',adapt2UNIX(SPARTAExec),' -in ',...
            adapt2UNIX(MROInputFile)];
        
        %...Loop over altitudes
        status = cell(size(simDensity));
        outcome = cell(size(simDensity));
        for d = 1:length(simDensity)
            %...Add folders to list
            dataFolder{d} = fullfile(OutputRepository,num2str(d));
            imageFolder{d} = fullfile(ImageRepository,num2str(d));
            
            %...Create folder if inexistent
            if ~exist(dataFolder{d},'dir'), mkdir(dataFolder{d}), end
            if ~exist(imageFolder{d},'dir'), mkdir(imageFolder{d}), end
            
            %...Gas fraction specifications
            gasFractions = '';
            for g = 1:length(gasNames)
                gasFractions = horzcat(gasFractions,...
                    sprintf('mixture %s frac %s\n',gasNames{g},num2str(gasRatios(g),'%.2f ')'));
            end
            gasFractions = gasFractions(1:end-1);
            
            %...Generate input file based on template
            fileID = fopen(MROInputFile,'w');
            fprintf(fileID,template,simBox,simGrid,numberDensity(d),real2sim(d),gasNamesString,...
                gasFractions,gasNamesString,streamVelocity(d,:),gasNamesString,...
                num2str(simAnglesOfAttack,'%.0f '),simTimeStep,...
                erase(dataFolder{d},[MRORepository,'/']),erase(imageFolder{d},[MRORepository,'/']),...
                simSteps);
            fclose(fileID);
            
            %...Run command via Terminal
            [status{d},outcome{d}] = system(commandString,'-echo');
            if status{d}
                system(['open ',adapt2UNIX(fullfile(MRORepository,'log.sparta')),' -a textedit']);
                error('SPARTA simulation failed. Opening log file.')
            end
        end
        
        %...Stop timer
        toc;
    otherwise
        %...Confirm to user
        warning('SPARTA analysis skipped.')
        
        %...Add folders to list
        for d = 1:length(simDensity)
            dataFolder{d} = fullfile(OutputRepository,num2str(d));
            imageFolder{d} = fullfile(ImageRepository,num2str(d));
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
momentArm = computeSolarPanelMomentArm(points,triangles(solarPanelElements.faces,:));

%% Analyze SPARTA Results

%...Get results for rarefied flow
pressureCoeffRarefied = cell(length(simAnglesOfAttack),length(simDensity));
frictionCoeffRarefied = cell(length(simAnglesOfAttack),length(simDensity));
aeroCoeffRarefied = cell(length(simAnglesOfAttack),length(simDensity));
bendingMomentRarefied = cell(length(simAnglesOfAttack),length(simDensity));
for d = 1:length(simDensity)
    for a = 1:length(simAnglesOfAttack)
        %...Get (average) distribution
        distribution = cell(size(simDensity));
        files = dir(fullfile(dataFolder{d},[num2str(simAnglesOfAttack(a)),'.coeff.*']));
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
        bendingMomentRarefied{a,d} = sum(cross(momentArm,(averageDistribution(solarPanelElements.faces,1:3) + ...
            averageDistribution(solarPanelElements.faces,4:6)) .* solarPanelTrianglesArea));
        
        %...Compute pressure and friction coefficients
        dynamicPressure = 1/2 * simDensity(d) * norm(streamVelocity(d,:))^2;
        pressureCoeffRarefied{a,d} = averageDistribution(:,1:3) / dynamicPressure;
        frictionCoeffRarefied{a,d} = averageDistribution(:,4:6) / dynamicPressure;
        
        %...Integrate to find force coefficients
        forceCoefficients = sum((pressureCoeffRarefied{a,d} + frictionCoeffRarefied{a,d}) .* ...
            trianglesArea) ./ crossSectionalArea(1);
        
        %...Find side, drag and lift coefficients
        %   Note that a negative sign is added since aerodynamic forces are
        %   positive in the negative direction
        aeroCoeffRarefied{a,d} = - forceCoefficients';
    end
end
clear distribution files fileID result averageDistribution dynamicPressure forceCoefficients transformation

%% Read Coefficients With Hess' Code

%...Get results for rarefied flow
aeroCoeffRarefiedHess = cell(length(simAnglesOfAttack),length(simDensity));
for d = 1:length(simDensity)
    for a = 1:length(simAnglesOfAttack)        
        aeroCoeffRarefiedHess{a,d} = obtainLiftDragCurveSPARTA(dataFolder{d},trianglesArea,simAnglesOfAttack,...
            simDensity(d),crossSectionalArea(1),norm(streamVelocity(d,:)));
    end
end
clear distribution files fileID result averageDistribution dynamicPressure forceCoefficients transformation

%% Bending Moment Analysis

%...Find maximum bending moment conditions for rarefied flow
bendingMomentRarefiedMagnitude = cellfun(@norm,bendingMomentRarefied);
maxBendingMomentRarefiedMagnitude = max(max(bendingMomentRarefiedMagnitude));
[i,j] = ind2sub(size(bendingMomentRarefiedMagnitude),...
    find(maxBendingMomentRarefiedMagnitude==bendingMomentRarefiedMagnitude));

%...Print results
fprintf(['Maximum bending moment in RAREFIED flow: %.3f Nm.\n ',...
    '\tConditions: AOA %.1f deg, Alt. %.2e km.\n'],maxBendingMomentRarefiedMagnitude,...
    simAnglesOfAttack(i),simDensity(j))

%% Plot Aerodynamic Coefficients

%...Plot
if showFigure
    %...Plot aerodynamic coefficients in 2D for rarefied flow
    styles = {'-o','-d','-s','-v','-p','-h','-*','-x','-^','-o','-d'}; labels = {'Drag','Side','Lift'}; locs = {'NE','NE','NW'};
    for i = 1:3
        figure;
        plot(simDensity,cellfun(@(x)x(i),aeroCoeffRarefied),'LineWidth',1.25,'MarkerSize',10)
        xlabel('Density [kg m^{-3}]')
        ylabel([labels{i},' Coefficient [-]'])
        if i == 2
            ylim([-1e-2,1e-2])
        end
        set(gca,'FontSize',15,'XScale','log')
        grid on
    end
end

%% Validation

%...Equation
dragDensity = @(d) 1.47952 - 0.032578 * log( d );

%...Plot
if showFigure
    figure;
    hold on
    plot(simDensity,cellfun(@(x)x(1),aeroCoeffRarefied(simAnglesOfAttack==0,:)),styles{1},'LineWidth',1.25,'MarkerSize',10)
    plot(simDensity,dragDensity(simDensity),'--','LineWidth',1.25)
%     plot(simDensity,1.79032 + 8.41026 ./ (log(simDensity)).^2,':','LineWidth',1.25)
    hold off
    xlabel('Density [km m^{-3}]')
    ylabel('Drag Coefficient [-]')
    set(gca,'FontSize',15,'XScale','log')
    grid on
    legend('Computed','Reference')
    ylim([1,3])
end

%...Print percentage offset
fprintf( ['Offset to validation data:\n',repmat('%.1f\t',[1,length(simDensity)]),'\n'], ...
    ( cellfun(@(x)x(1),aeroCoeffRarefied(simAnglesOfAttack==0,:)) - dragDensity(simDensity) ) ./ ...
    dragDensity(simDensity) * 100)

%% Plot Pressure Distributions

%...Create triangulation
Tri = [];
Tri.vertices = points;
Tri.faces = triangles;

%...Plot
if showFigure
    %...Plot triangulation for rarefied flow
    titles = {'Pressure','Shear'};
    figure;
    j = 0;
    color = {sqrt(sum(pressureCoeffRarefied{1,1}.^2,2)),sqrt(sum(frictionCoeffRarefied{1,1}.^2,2))};
    for i = 1:2
        j = j+1;
        subplot(1,2,j)
        Tri.facevertexcdata = color{i};
        patch(Tri), shading faceted, colormap jet
        c = colorbar; c.Label.String = 'Pa'; c.Location = 'southoutside';
        set(gca,'FontSize',15), view([127.5,30])
        axis off tight equal
        title([num2str(simAnglesOfAttack(1)),' deg'])
    end
end
clear F j color c

%% GIF Commands

if strcmpi(selection,'proceed')
    %...Open terminal
    ! open /Applications/Utilities/Terminal.app/
    
    %...Create command string
    GIFCommand = ['export MAGICK_HOME="/Users/Michele/Software/ImageMagick-7.0.7"; '...
        'export PATH="$MAGICK_HOME/bin:$PATH"; export DYLD_LIBRARY_PATH="$MAGICK_HOME/lib/"; '...
        'cd ',adapt2UNIX(MRORepository)];
    for d = 1:length(simDensity)
        ds = num2str(d);
        GIFCommand = [GIFCommand,'; convert '];
        for a = simAnglesOfAttack
            as = num2str(a);
            GIFCommand = [GIFCommand,sprintf('%s ',fullfile('figures',ds,[as,'.image.*']))];
        end
        GIFCommand = [GIFCommand,fullfile('videos',['movie_',ds,'.gif'])];
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
    axes = eye(3);
    crossSectionalArea = zeros(1,3);
    for i = 1:length(axes)
        crossSectionalArea(i) = sum(abs(trianglesNormal * axes(:,i)) .* trianglesArea)/2;
    end
end

function momentArm = computeSolarPanelMomentArm(points,triangles)
    %...Function handle for surface normal
    centroid = @(p) [sum(p(:,1)),sum(p(:,2)),sum(p(:,3))]/3;

    %...Get triangle vertices
    vertices = arrayfun(@(i)points(triangles(i,:),:),1:size(triangles,1),'UniformOutput',false);

    %...Compute surface normal
    trianglesCentroid = cellfun(@(x)centroid(x),vertices,'UniformOutput',false)'; % compute centroid
    trianglesCentroid = cell2mat(trianglesCentroid);
    
    %...Find distance to center
    momentArm = trianglesCentroid - [0,1.25,0];
end