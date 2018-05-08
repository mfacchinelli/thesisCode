fclose('all'); clear all; close all force; profile off; clc; format long g; rng default;
addpath functions tests other

%...Set system path
setenv('PATH',[getenv('PATH'),':/Users/Michele/Software/ImageMagick-7.0.7/bin:/opt/local/bin:',...
    '/opt/local/sbin:/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin:/Library/TeX/texbin:/opt/X11/bin']);

%...Set ImageMagick path
setenv('MAGICK_HOME','/Users/Michele/Software/ImageMagick-7.0.7')
setenv('DYLD_LIBRARY_PATH','$MAGICK_HOME/ImageMagick-7.0.7/lib/')

%% SPARTA and Modeling Variables

%...Select mode
mode = 'speedRatio'; % 'sphere' or 'speedRatio'

%...Figures setting
showFigure = true;
saveFigure = false;
[figSizeLarge,figSizeSmall] = saveFigureSettings(saveFigure);

%...SPARTA
SPARTAExec = '/Users/Michele/AE Software/SPARTA/src/spa_mac_mpi'; % path to SPARTA executable

%...Repositories
SPARTARepository = ['/Users/Michele/Library/Mobile Documents/com~apple~CloudDocs/University/Master ',...
    'Thesis/Code/SPARTA/']; % path to SPARTA folder
SphereRepository = fullfile(SPARTARepository,'verification'); % path to sphere folder
OutputRepository = fullfile(SphereRepository,mode,'data'); % path to SPARTA results
ImageRepository = fullfile(SphereRepository,mode,'figures'); % path to SPARTA figures

%...Files
SphereInputFile = fullfile(SphereRepository,'in.mro'); % path to MRO input file
switch mode
    case 'sphere'
        SphereDataFile = fullfile(SphereRepository,'data.sphere'); % path to MRO model file
    case 'speedRatio'
        SphereDataFile = fullfile(SphereRepository,'data.mro'); % path to MRO model file
end

%...Templates
SphereInputTemplate = fullfile(SPARTARepository,'SPARTAInputTemplateVerification.txt');

%...Mars Climate Database
MCDData = fullfile(SPARTARepository,'MCDEnvir');

%% Spacecraft Model

%...Compute model based on MRO discretization
switch mode
    case 'sphere'
        [SphereExtent,points,triangles] = ReadSPARTASphere(SphereDataFile,showFigure);
    case 'speedRatio'
        [SphereExtent,points,triangles] = SPARTAModel(SphereDataFile,false,showFigure,false);
end

%...Get variables for simulation
simBox = reshape(SphereExtent',1,numel(SphereExtent));

%% Constants

%...Grid specifications
switch mode
    case 'sphere'
        gridSpacing = 0.1; % grid size
        r2sOffset = 10; % number of simulated particles per cell
    case 'speedRatio'
        gridSpacing = 0.25; % grid size
        r2sOffset = 15; % number of simulated particles per cell
end
simGrid = diff(SphereExtent,[],2)/gridSpacing;

%...Time settings
simSteps = 1000; % simulation time (specified as number of steps)

%% Climate Model

%...MRO conditions
MRONumberDensity = [7.74999177706116e+16,3.61518614036900e+15,101127687472605,3801215185302.87];
MROMolecularSpeedRatios = [14.0299656853382,11.8991913642848,10.1037171207576,8.87895229868620];
MRODragCoefficient = [2.04746001998197,1.83088934893315,2.61653426222919,1.69965588176817];
MROLiftCoefficient = [0.0116989119634681,0.00585381331181692,0.00786808650764345,0.0101542468626223];

%...Flow constants
simNumberDensities = logspace(17,12,4);
temperature = 150; % gas mixture temperature

%...Create cell arrays of gas information for each altitude
switch mode
    case 'sphere'
        gasNames = {''};
    case 'speedRatio'
        gasNames = {'CO2','O','H'};
end
molarMasses = [44.01,1.01,16.00]*1e-3;
gamma = [1.28,1.67,1.67];
gasConstant = 8.3145 ./ molarMasses;

%...Get dependent values
density = simNumberDensities / 6.0221e23 .* molarMasses';
pressure = simNumberDensities / 6.0221e23 * 273.15 * 8.3145;
speedOfSound = sqrt(gamma .* gasConstant * temperature);

%...Angles of attack
switch mode
    case 'sphere'
        simAnglesOfAttack = [-15,0,15];
    case 'speedRatio'
        simAnglesOfAttack = -30:15:30;
end

%...Molecular speed ratio
simMolecularSpeedRatios = [1,2.5,5,10,25,50,100];
Mach = simMolecularSpeedRatios .* sqrt(2./gamma)';

%...Find circular velocity at altitude
streamVelocity = zeros(length(simMolecularSpeedRatios),3,length(gasNames));
for g = 1:length(gasNames)
    streamVelocity(:,1,g) = - Mach(g,:) * speedOfSound(g);
end
simTimeStep = round(0.09*diff(SphereExtent(1,:))./sqrt(sum(streamVelocity.^2,2)),6); % time step (1/10-th box traverse time)

%...Set ratio of real to simulated particles
real2sim = simNumberDensities / (r2sOffset / gridSpacing^3);

%% Generate SPARTA Input File

%...Get input template
fileID = fopen(SphereInputTemplate,'r');
template = textscan(fileID,'%s','Delimiter','\n','CommentStyle','#');
template = join(template{1},'\n'); template = template{1};
fclose(fileID);

%% Run SPARTA

%...Function handle to adapt paths to UNIX
adapt2UNIX = @(x) replace(x,' ','\ ');

%...Ask user to proceed
selection = questdlg('Proceeding with the analysis will overwrite all previous data.',...
    'Warning','Proceed','Stop','Stop');

%...Switch based on user input
switch lower(selection)
    case 'proceed'
        %...Start timer
        tic;
        
        %...Clean up folders
        system(['rm -rf ',fullfile(adapt2UNIX(OutputRepository),'*')]);
        system(['rm -rf ',fullfile(adapt2UNIX(ImageRepository),'*')]);
        system(['rm -rf ',fullfile(adapt2UNIX(SphereRepository),'videos/*')]);
        
        %...Generate command
        commandString = ['cd ',adapt2UNIX(SphereRepository),';',...
            'mpirun -np 2 ',adapt2UNIX(SPARTAExec),' -in ',...
            adapt2UNIX(SphereInputFile)];
        
        %...Loop over altitudes
        status = cell(length(simMolecularSpeedRatios),length(simNumberDensities),length(gasNames));
        outcome = cell(length(simMolecularSpeedRatios),length(simNumberDensities),length(gasNames));
        for g = 1:length(gasNames)
            for d = 1:length(simNumberDensities)
                for s = 1:length(simMolecularSpeedRatios)
                    %...Add folders to list
                    dataFolder = fullfile(OutputRepository,gasNames{g},...
                        num2str(log10(simNumberDensities(d)),'%.0f'),num2str(simMolecularSpeedRatios(s),'%.0f'));
                    imageFolder = fullfile(ImageRepository,gasNames{g},...
                        num2str(log10(simNumberDensities(d)),'%.0f'),num2str(simMolecularSpeedRatios(s),'%.0f'));
                    
                    %...Create folder if inexistent
                    if ~exist(dataFolder,'dir'), mkdir(dataFolder), end
                    if ~exist(imageFolder,'dir'), mkdir(imageFolder), end
                    
                    %...Generate input file based on template
                    fileID = fopen(SphereInputFile,'w');
                    fprintf(fileID,template,simBox,simGrid,simNumberDensities(d),real2sim(d),gasNames{g},...
                        gasNames{g},streamVelocity(s,:,g),gasNames{g},num2str(simAnglesOfAttack),simTimeStep(s,1,g),...
                        erase(dataFolder,[SphereRepository,'/']),erase(imageFolder,[SphereRepository,'/']),...
                        simSteps);
                    fclose(fileID);
                    
                    %...Run command via Terminal
                    [status{s,d,g},outcome{s,d,g}] = system(commandString,'-echo');
                    if status{s,d,g}
                        system(['open ',adapt2UNIX(fullfile(SphereRepository,'log.sparta')),' -a textedit']);
                        error('SPARTA simulation failed. Opening log file.')
                    end
                end
            end
        end
        
        %...Stop timer
        toc;
    otherwise
        %...Confirm to user
        warning('SPARTA analysis skipped.')
end
clear commandString dataFolder imageFolder fileID

%% Analyze SPARTA Results

%...Compute area and normal to surface of triangles and surface areas around each axis
[trianglesArea,trianglesNormal,crossSectionalArea] = computeTriangleAreaNormal(points,triangles);

%...Get results for rarefied flow
pressureCoeffRarefied = cell(length(simAnglesOfAttack),length(simMolecularSpeedRatios),...
    length(simNumberDensities),length(gasNames));
frictionCoeffRarefied = cell(length(simAnglesOfAttack),length(simMolecularSpeedRatios),...
    length(simNumberDensities),length(gasNames));
aeroCoeffRarefied = cell(length(simAnglesOfAttack),length(simMolecularSpeedRatios),...
    length(simNumberDensities),length(gasNames));
aeroCoeffRarefiedHess = cell(length(simAnglesOfAttack),length(simMolecularSpeedRatios),...
    length(simNumberDensities),length(gasNames));
for g = 1:length(gasNames)
    for d = 1:length(simNumberDensities)
        for s = 1:length(simMolecularSpeedRatios)
            for a = 1:length(simAnglesOfAttack)
                %...Get (average) distribution
                distribution = cell(size(simMolecularSpeedRatios));
                switch mode
                    case 'sphere'
                        files = dir(fullfile(OutputRepository,num2str(simMolecularSpeedRatios(s),'%.0f'),...
                            [num2str(simAnglesOfAttack(a)),'.coeff.*']));
                    case 'speedRatio'
                        files = dir(fullfile(OutputRepository,gasNames{g},num2str(log10(simNumberDensities(d)),'%.0f'),...
                            num2str(simMolecularSpeedRatios(s),'%.0f'),[num2str(simAnglesOfAttack(a)),'.coeff.*']));
                end
                for i = 1:length(files)
                    fileID = fopen(fullfile(files(i).folder,files(i).name),'r');
                    result = textscan(fileID,repmat('%f ',[1,10]),'HeaderLines',9,'CollectOutput',true);
                    result = sortrows(result{1},1);
                    distribution{i} = result(:,2:end);
                    fclose(fileID);
                end
                    
                %...Take median to minimize effect of outliers
                %   Note that the first 2 files are NOT used, since the simulation
                %   is still starting up
                averageDistribution = median(cat(3,distribution{3:end}),3);
                
                %...Compute pressure and friction coefficients
                dynamicPressure = 1/2 * density(g,d) * norm(streamVelocity(s,:,g))^2;
                pressureCoeffRarefied{a,s,d,g} = (averageDistribution(:,1:3) - pressure(d) * trianglesNormal * ...
                    roty(simAnglesOfAttack(a))) / dynamicPressure; % rotate normal to account for angle of attack
                frictionCoeffRarefied{a,s,d,g} = averageDistribution(:,4:6) / dynamicPressure;
                
                %...Integrate to find force coefficients
                %   Note that a negative sign is added since aerodynamic forces are
                %   positive in the negative direction
                forceCoefficients = sum((pressureCoeffRarefied{a,s,d,g} + frictionCoeffRarefied{a,s,d,g}) .* ...
                    trianglesArea) ./ crossSectionalArea(1);
                %     forceCoefficients = sum(averageDistribution(:,7:9))./crossSectionalArea(1)/dynamicPressure;
                
                %...Find side, drag and lift coefficients
                %   Note that MATLAB has a different definition of positive rotations
                aeroCoeffRarefied{a,s,d,g} = - forceCoefficients';
                
                %...Read coefficients with Hess' code
                aeroCoeffRarefiedHess{a,s,d,g} = obtainLiftDragCurveSPARTA(averageDistribution,...
                    trianglesArea,simAnglesOfAttack(a),density(g,d),crossSectionalArea(1),...
                    norm(streamVelocity(s,:,g)));
            end
        end
    end
end

%% Plot Aerodynamic Coefficients

%...Plot
if showFigure
    %...Plot aerodynamic coefficients in 2D for rarefied flow
    styles = {'-o','-d','-s','-v','-p','-h','-*','-x','-o','-d','-s','-v','-p','-h','-*','-x'};
    labels = {'Drag','Side','Lift'};
    if strcmpi(mode,'sphere') || ~saveFigure, F = figure('rend','painters','pos',figSizeSmall); end
    for g = 1:length(gasNames)
        for i = 1:2:3
            if strcmpi(mode,'sphere') || ~saveFigure, subplot(length(gasNames),2,i - 0.5*(i-1) + 2*(g-1)),
            else, F = figure('rend','painters','pos',figSizeSmall); end
            hold on
            surfs = cell(1,length(simAnglesOfAttack));
            colors = {'blue','red','green','yellow','magenta'};
            for a = 1:length(simAnglesOfAttack)
                surfs{a} = surf(simMolecularSpeedRatios,simNumberDensities,...
                    squeeze(cellfun(@(x)x(i),aeroCoeffRarefied(a,:,:,g)))');%,...
%                    'FaceColor',colors{a},'FaceAlpha',0.5);
            end
            hold off
            xlabel('Molecular Speed Ratio [-]')
            ylabel('Number Density [-]')
            zlabel([labels{i},' Coefficient [-]'])
            set(gca,'FontSize',15,'XScale','log','YScale','log')
%             legend([surfs{:}],split(num2str(simAnglesOfAttack)),'Location','Best');
            view([15,15])
            if ~saveFigure, title(gasNames{g}), end
            grid on
%             if saveFigure
%                 switch mode
%                     case 'sphere'
%                         saveas(F,'../../Report/figures/ver_rare_3d','epsc')
%                     case 'speedRatio'
%                         saveas(F,['../../Report/figures/speed_rare_3d_',gasNames{g},'_',lower(labels{i})],'epsc')
%                 end
%             end
        end
    end
    if ~saveFigure, subplotTitle('Rarefied'), end
    
    %...Plot aerodynamic coefficients in 2D for rarefied flow with MRO data overlayed
    if strcmpi(mode,'sphere') || ~saveFigure, F = figure('rend','painters','pos',figSizeSmall); end
    for g = 1:length(gasNames)
        for i = 1:2:3
            if strcmpi(mode,'sphere') || ~saveFigure, subplot(length(gasNames),2,i - 0.5*(i-1) + 2*(g-1)),
            else, F = figure('rend','painters','pos',figSizeSmall); end
            hold on
            surf(simMolecularSpeedRatios,simNumberDensities,...
                squeeze(cellfun(@(x)x(i),aeroCoeffRarefied(3,:,:,g)))',... % plot only 0 deg AOA
                'FaceAlpha',0.5);
            expectedValues = interp2(simMolecularSpeedRatios,simNumberDensities,squeeze(cellfun(@(x)x(i),...
                aeroCoeffRarefied(3,:,:,g)))',MROMolecularSpeedRatios,MRONumberDensity);
            if i == 1
                scatter3(MROMolecularSpeedRatios,MRONumberDensity,MRODragCoefficient,'filled','red')
                scatter3(MROMolecularSpeedRatios,MRONumberDensity,expectedValues,'filled','green')
            elseif i == 3
                scatter3(MROMolecularSpeedRatios,MRONumberDensity,MROLiftCoefficient,'filled','red')
                scatter3(MROMolecularSpeedRatios,MRONumberDensity,expectedValues,'filled','green')
            end
            hold off
            xlabel('Molecular Speed Ratio [-]')
            ylabel('Number Density [-]')
            zlabel([labels{i},' Coefficient [-]'])
            set(gca,'FontSize',15,'XScale','log','YScale','log')
            view([15,15])
            if ~saveFigure, title(gasNames{g}), end
            grid on
%             if saveFigure
%                 switch mode
%                     case 'sphere'
%                         saveas(F,'../../Report/figures/ver_rare_2d','epsc')
%                     case 'speedRatio'
%                         saveas(F,['../../Report/figures/speed_rare_2d_',gasNames{g},'_',lower(labels{i})],'epsc')
%                 end
%             end
        end
    end
    if ~saveFigure, subplotTitle('Rarefied'), end
end

%% Validation

%...Load MCD
MCD = load(MCDData);

%...Get drag coefficients at zero angle of attack
dragCoefficients = squeeze(cellfun(@(x)x(1),aeroCoeffRarefied(simAnglesOfAttack==0,:,:,:)));

%...Get interpolated values
dragValues = zeros(4,4,3);
for g = 1:length(gasNames)
    [N,M] = meshgrid(MRONumberDensity,MROMolecularSpeedRatios);
    dragValues(:,:,g) = interp2(simNumberDensities,simMolecularSpeedRatios,dragCoefficients(:,:,g),...
        N,M);
end

%...Equation
dragDensity = @( d ) 1.47952 - 0.032578 * log( d );

%...Plot
if showFigure
    for g = 1:length(gasNames)
        F = figure('rend','painters','pos',figSizeSmall);
        hold on
        plot(densityValues,dragValues(:,:,g),...
            '-o','LineWidth',1.25,'MarkerSize',10)
        plot(densityValues,dragDensity(densityValues),'--','LineWidth',1.25)
        hold off
        xlabel('Density [kg m^{-3}]')
        ylabel('Drag Coefficient [-]')
        set(gca,'FontSize',15,'XScale','log')
        grid on
        legend('Computed','Reference')
        if ~saveFigure, title(gasNames{g}), end
%         if saveFigure, saveas(F,['../../Report/figures/valid_gallis'],'epsc'), end
    end
end

%% GIF Commands

if strcmpi(selection,'proceed')
    %...Open terminal
    ! open /Applications/Utilities/Terminal.app/
    
    %...Create command string
    GIFCommand = ['export MAGICK_HOME="/Users/Michele/Software/ImageMagick-7.0.7"; '...
        'export PATH="$MAGICK_HOME/bin:$PATH"; export DYLD_LIBRARY_PATH="$MAGICK_HOME/lib/"; '...
        'cd ',adapt2UNIX(SphereRepository)];
    for g = gasNames
        for s = simMolecularSpeedRatios
            ss = num2str(s,'%.0f');
            GIFCommand = [GIFCommand,'; convert '];
            for a = simAnglesOfAttack
                as = num2str(a);
                GIFCommand = [GIFCommand,sprintf('%s ',fullfile(mode,'figures',g{1},ss,[as,'.image.*']))];
            end
            GIFCommand = [GIFCommand,fullfile(mode,'videos',g{1},['movie_',g{1},'_',ss,'.gif'])];
        end
    end
    
    %...Add command string to clipboard
    clipboard('copy',GIFCommand)
end

%% Supporting Functions

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