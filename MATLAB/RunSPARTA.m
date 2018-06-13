fclose('all'); clear all; close all force; profile off; clc; format long g; rng default;
addpath functions

%...Set system path
setenv('PATH',[getenv('PATH'),':/Users/Michele/Software/ImageMagick-7.0.7/bin:/opt/local/bin:',...
    '/opt/local/sbin:/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin:/Library/TeX/texbin:/opt/X11/bin']);

%...Set ImageMagick path
setenv('MAGICK_HOME','/Users/Michele/Software/ImageMagick-7.0.7')
setenv('DYLD_LIBRARY_PATH','$MAGICK_HOME/ImageMagick-7.0.7/lib/')

%% SPARTA and Modeling Variables

%...Rotate by 180 degrees
rotateMRO = false;
if rotateMRO
    repository = 'rotated';
else
    repository = 'nominal';
end

%...Figures and tables setting
saveTable = false;
showFigure = true;
saveFigure = false;
[figSizeLarge,~,figSizeSmall] = saveFigureSettings(saveFigure);

%...Host file settings
useHostFile = false; % run SPARTA simulation with 4 cores (slows down computer a lot)

%...SPARTA
SPARTAExec = '/Users/Michele/Software/SPARTA/src/spa_mac_mpi'; % path to SPARTA executable

%...Repositories
SPARTARepository = ['/Users/Michele/Library/Mobile Documents/com~apple~CloudDocs/University/Master ',...
    'Thesis/Code/SPARTA/']; % path to SPARTA folder
MRORepository = fullfile(SPARTARepository,'mro'); % path to MRO folder
OutputRepository = fullfile(MRORepository,repository,'data'); % path to SPARTA results
ImageRepository = fullfile(MRORepository,repository,'figures'); % path to SPARTA figures
VideoRepository = fullfile(MRORepository,repository,'videos'); % path to SPARTA videos

%...Files
MROInputFile = fullfile(MRORepository,'in.mro'); % path to MRO input file
MRODataFile = fullfile(MRORepository,'data.mro'); % path to MRO model file

%...Templates
MROInputTemplate = fullfile(SPARTARepository,'SPARTAInputTemplate.txt');

%...Mars Climate Database
MCDData = fullfile(SPARTARepository,'MCDEnvir');

%% Spacecraft Model

%...Compute model based on MRO discretization
[MROExtent,points,triangles,solarPanelElements,antennaElement] = ...
    SPARTAModel(MRODataFile,rotateMRO,showFigure,false); % output min and max dimensions of MRO

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
r2sOffset = 20; % number of simulated particles per cell

%...Simulation conditions
simAnglesOfAttack = linspace(-30,30,13); % angles of attack for simulation
%[linspace(-75,-30,4),linspace(-25,25,11),linspace(30,75,4)];
simAltRarefied = [100,125,150,200,300,500]; % altitudes for rarefied regime simulations
simAltContinuum = 65; % altitudes for continuum regime simulations

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
gasRatios{4} = [0.12,0.28,0.18,0.42]; % make sure it all sums up to one
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

%...Ask user to proceed
% selection = questdlg('Proceeding with the analysis will overwrite all previous data.',...
%     'Warning','Proceed','Stop','Stop');
selection = 'stop';

%...Folder names for SPARTA output
dataFolder = cell(size(simAltRarefied));
imageFolder = cell(size(simAltRarefied));

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
        if useHostFile
            commandString = ['cd ',adapt2UNIX(MRORepository),';',...
                'mpirun -hostfile hostfile -np 4 ',adapt2UNIX(SPARTAExec),' -in ',...
                adapt2UNIX(MROInputFile)];
        else
            commandString = ['cd ',adapt2UNIX(MRORepository),';',...
                'mpirun -np 2 ',adapt2UNIX(SPARTAExec),' -in ',...
                adapt2UNIX(MROInputFile)];
        end
        
        %...Loop over altitudes
        status = cell(size(simAltRarefied));
        outcome = cell(size(simAltRarefied));
        for h = 1:length(simAltRarefied)
            %...Add folders to list
            dataFolder{h} = fullfile(OutputRepository,num2str(simAltRarefied(h)));
            imageFolder{h} = fullfile(ImageRepository,num2str(simAltRarefied(h)));
            
            %...Create folder if inexistent
            if ~exist(dataFolder{h},'dir'), mkdir(dataFolder{h}), end
            if ~exist(imageFolder{h},'dir'), mkdir(imageFolder{h}), end
            
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
                erase(dataFolder{h},[MRORepository,'/']),erase(imageFolder{h},[MRORepository,'/']),...
                simSteps);
            fclose(fileID);
            
            %...Run command via Terminal
            [status{h},outcome{h}] = system(commandString,'-echo');
            if status{h}
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
        for h = 1:length(simAltRarefied)
            dataFolder{h} = fullfile(OutputRepository,num2str(simAltRarefied(h)));
            imageFolder{h} = fullfile(ImageRepository,num2str(simAltRarefied(h)));
        end
end
clear commandString gasFractions fileID

%% Area Analysis

%...Compute area and normal to surface of triangles and surface areas
%   around each axis for whole spacecraft
[trianglesArea,trianglesNormal,crossSectionalArea] = computeTriangleAreaNormal(points,triangles);
if rotateMRO, referenceArea = 3; else, referenceArea = 1; end
referenceLength = 2.5;
centerOfMass = [0,0,0.1375];
[momentArm,trianglesCentroid] = computeMomentArm(points,triangles,centerOfMass);

%...Compute same data only for solar panels
[solarPanelTrianglesArea,solarPanelTrianglesNormal,solarPanelCrossSectionalArea] = ...
    computeTriangleAreaNormal(points,triangles(solarPanelElements.faces,:));
solarPanelAttachmentPoint = [0,1.25,0];
momentArmSolarPanel = computeMomentArm(points,triangles(solarPanelElements.faces,:),solarPanelAttachmentPoint);

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
        solarPanelDistribution = averageDistribution(solarPanelElements.faces,1:3) + ...
            averageDistribution(solarPanelElements.faces,4:6);
        rotatedDistribution = arrayfun(@(i)roty(simAnglesOfAttack(a))' * ...
            solarPanelDistribution(i,:)',1:size(solarPanelDistribution,1),'UniformOutput',false);
        bendingMomentRarefied{a,h} = sum(cross(momentArmSolarPanel,...
            [rotatedDistribution{:}]') .* solarPanelTrianglesArea);
        
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
clear distribution files fileID result averageDistribution solarPanelDistribution rotatedDistribution ...
    dynamicPressure forceCoefficients momentCoefficients

%% Compute Continuum Flow

%...Function handles
stagnationPressure = @(gamma,Mach) 2 ./ gamma ./ Mach.^2 .* ( ( ( (gamma+1).^2 .* Mach.^2 ) ./ ...
    ( 4 * gamma .* Mach.^2 - 2 * (gamma-1) ) ).^(gamma./(gamma-1)) .* ...
    ( ( 2 * gamma .* Mach.^2 - (gamma-1) ) ./ (gamma+1) ) - 1 );
basePressure = @(gamma,Mach) 2 ./ gamma ./ Mach.^2 .* ( ( 1 ./ Mach.^2 .* ( 2 ./ (gamma+1) ) ).^(1.4) .* ...
    ( ( 2 * gamma .* Mach.^2 - (gamma-1) ) ./ (gamma+1) ) - 1 );

%...Extract MCD parameters
gamma = MCD.specificHeatRatio;
speedOfSound = MCD.speedOfSound;

%...Interpolate to find data for rarefied regime
[gamma,speedOfSound] = interpolate(altitude,simAltContinuum,gamma,speedOfSound);

%...Set up for analysis
Mach = MarsCircularVelocity(simAltContinuum) ./ speedOfSound; % Mach number
velocityVector = [-1;0;0]; % direction of incoming flow

%...Compute basis pressure coefficients
cps = stagnationPressure(gamma,Mach); % stagnation pressure
cpb = basePressure(gamma,Mach); % base pressure

%...Loop over altitudes and angles of attack
pressureCoeffContinuum = cell(length(simAnglesOfAttack),length(simAltContinuum));
aeroCoeffContinuum = cell(length(simAnglesOfAttack),length(simAltContinuum));
bendingMomentContinuum = cell(length(simAnglesOfAttack),length(simAltContinuum));
for h = 1:length(simAltContinuum)
    for a = 1:length(simAnglesOfAttack)
        %...Find velocity vector at angle of attack
        %   Note that MATLAB has a different definition of positive rotations
        transformation = roty(simAnglesOfAttack(a));
        V = transformation * velocityVector;
        
        %...Find incidence angle
        sineIncidenceAngle = - trianglesNormal * V; % dot product
        locPos = sineIncidenceAngle > 0;
        locNeg = sineIncidenceAngle <= 0;
        
        %...Find pressure coefficient
        pressureCoeffContinuum{a,h} = zeros(size(sineIncidenceAngle));
        pressureCoeffContinuum{a,h}(locPos) = cps(h) * sineIncidenceAngle(locPos).^2;
        pressureCoeffContinuum{a,h}(locNeg) = cpb(h);
        pressureCoeffContinuum{a,h} = pressureCoeffContinuum{a,h} .* trianglesNormal;
        
        %...Compute bending moment
        dynamicPressure = 1/2 * density(h) * norm(streamVelocity(h,:))^2;
        solarPanelDistribution = pressureCoeffContinuum{a,h}(solarPanelElements.faces,:);
        rotatedDistribution = arrayfun(@(i)roty(simAnglesOfAttack(a))' * ...
            solarPanelDistribution(i,:)',1:size(solarPanelDistribution,1),'UniformOutput',false);
        bendingMomentContinuum{a,h} = sum(cross(momentArmSolarPanel,...
            [rotatedDistribution{:}]') * ...
            dynamicPressure .* solarPanelTrianglesArea);
%         bendingMomentContinuum{a,h} = sum(cross(momentArmSolarPanel,...
%             pressureCoeffContinuum{a,h}(solarPanelElements.faces,:)) * ...
%             dynamicPressure .* solarPanelTrianglesArea);
        
        %...Integrate to find force coefficients
        forceCoefficients = sum(pressureCoeffContinuum{a,h} .* trianglesArea) ./ crossSectionalArea(referenceArea);
        
        %...Integrate to find moment coefficients
        momentCoefficients = sum(cross(momentArm,pressureCoeffContinuum{a,h}) .* ...
            trianglesArea) / crossSectionalArea(referenceArea) / referenceLength;
        
        %...Find aerodynamic coefficients
        aeroCoeffContinuum{a,h} = [ transformation' * forceCoefficients'; ...
            transformation' * momentCoefficients']; % inverse rotation
    end
end
clear transformation V sineIncidenceAngle locPos locNeg dynamicPressure solarPanelDistribution rotatedDistribution ...
    forceCoefficients momentCoefficients

%% Compute Transition Regime

%...Bridging function
bridge = @(Kn) sin( pi/8 * ( 3 + log10(Kn) ) ).^2;

%...Extract MCD parameters
knudsen = MCD.knudsenNumber;

%...Interpolate with combined altitudes
simAltTotal = linspace(simAltContinuum(end),simAltRarefied(2),6); % use 125 km instead of 100 km
simAltTransition = simAltTotal(2:end-1); % remove extremes
knudsen = interpolate(altitude,simAltTransition,knudsen);

%...Combine coefficients
%   Note that here, the values at 125 km are selected as reference
aeroCoeffTotal = horzcat(aeroCoeffContinuum(:,end),aeroCoeffRarefied(:,2));

%...Find transition coefficients
aeroCoeffTranstion = cell(length(simAnglesOfAttack),length(simAltTransition));
for h = 1:length(simAltTransition)
    for a = 1:length(simAnglesOfAttack)
        aeroCoeffTranstion{a,h} = (aeroCoeffTotal{a,2} - aeroCoeffTotal{a,1}) * ...
            bridge(knudsen(h)) + aeroCoeffTotal{a,1};
    end
end

%...Combine coefficients to add transition
aeroCoeffTotal = horzcat(aeroCoeffTotal(:,1),aeroCoeffTranstion,aeroCoeffTotal(:,2));

%...Combine coefficients without transition
aeroCoefficients = horzcat(aeroCoeffContinuum,aeroCoeffTranstion,aeroCoeffRarefied(:,2:end));
simAltitudes = horzcat(simAltContinuum,simAltTransition,simAltRarefied(2:end));

%% Bending Moment Analysis

%...Find maximum bending moment conditions for rarefied flow
bendingMomentRarefiedMagnitude = cellfun(@norm,bendingMomentRarefied);
maxBendingMomentRarefiedMagnitude = max(max(bendingMomentRarefiedMagnitude));
[i,j] = ind2sub(size(bendingMomentRarefiedMagnitude),...
    find(maxBendingMomentRarefiedMagnitude==bendingMomentRarefiedMagnitude));

%...Print results
fprintf(['Maximum bending moment in RAREFIED flow: %.3f Nm.\n ',...
    '\tConditions: AOA %.1f deg, Alt. %.0f km.\n'],maxBendingMomentRarefiedMagnitude,...
    simAnglesOfAttack(i),simAltRarefied(j))

%...Find maximum bending moment conditions for continuum flow
bendingMomentContinuumMagnitude = cellfun(@norm,bendingMomentContinuum);
maxBendingMomentContinuumMagnitude = max(max(bendingMomentContinuumMagnitude));
[i,j] = ind2sub(size(bendingMomentContinuumMagnitude),...
    find(maxBendingMomentContinuumMagnitude==bendingMomentContinuumMagnitude));

%...Print results
fprintf(['Maximum bending moment in CONTINUUM flow: %.3f Nm.\n ',...
    '\tConditions: AOA %.1f deg, Alt. %.0f km.\n'],maxBendingMomentContinuumMagnitude,...
    simAnglesOfAttack(i),simAltContinuum(j))

%% Plot Aerodynamic Coefficients

%...Plot
if showFigure
    styles = {'-o','-d','-s','-v','-p','-h','-*','-x','-^','-o','-d'};
    labels = {'Drag','Side','Lift','X-Moment','Y-Moment','Z-Moment'};
    if ~saveFigure
        %...Plot aerodynamic coefficients in 2D for rarefied and continuum flows
        for i = 1:6
            figure('rend','painters','pos',figSizeSmall);
            hold on
            for h = 1:length(simAltitudes)
                plot(simAnglesOfAttack,cellfun(@(x)x(i),aeroCoefficients(:,h)),styles{h},'LineWidth',1.25,'MarkerSize',10)
            end
            hold off
            xlabel('Angle of Attack [deg]')
            ylabel([labels{i},' Coefficient [-]'])
            set(gca,'FontSize',15)
            grid on
            legend(split(num2str(simAltitudes)),'Location','Best')
        end
        
        %...Plot drag polar
        figure('rend','painters','pos',figSizeSmall);
        hold on
        for h = 1:length(simAltitudes)
            plot(cellfun(@(x)x(3),aeroCoefficients(:,h)),cellfun(@(x)x(1),aeroCoefficients(:,h)),...
                styles{h},'LineWidth',1.25,'MarkerSize',10)
        end
        hold off
        xlabel('Lift Coefficient [-]')
        ylabel('Drag Coefficient [-]')
        set(gca,'FontSize',15)
        grid on
        legend(split(num2str(simAltitudes)),'Location','Best')
    else
        %...Plot aerodynamic coefficients in 2D for rarefied flow
        for i = 1:6
            F = figure('rend','painters','pos',figSizeSmall);
            hold on
            for h = 2:length(simAltRarefied) % skip 100 km
                plot(simAnglesOfAttack,cellfun(@(x)x(i),aeroCoeffRarefied(:,h)),styles{h},'LineWidth',1.25,'MarkerSize',10)
            end
            hold off
            xlabel('Angle of Attack [deg]')
            ylabel([labels{i},' Coefficient [-]'])
            if i == 2 || i == 4
                ylim([-1e-2,1e-2])
            elseif i == 6
                ylim([-1e-2,1e-2])%ylim([-1e-1,1e-1])
            end
            set(gca,'FontSize',15)
            grid on
            if i == 1
                legend(split(num2str(simAltRarefied(2:end))),'NumColumns',2,'Location','Best')
            elseif i >= 4
                legend(split(num2str(simAltRarefied(2:end))),'Location','SW')
            else
                legend(split(num2str(simAltRarefied(2:end))),'Location','Best')
            end
            if saveFigure, saveas(F,['../../Report/figures/aero_rare_2d_',lower(labels{i})],'epsc'), end
        end
        
        F = figure('rend','painters','pos',figSizeSmall);
        hold on
        plot(simAnglesOfAttack,cellfun(@(x)x(5),aeroCoeffRarefied(:,2)), ...
            styles{1},'LineWidth',1.25,'MarkerSize',10) % full moment
        plot(simAnglesOfAttack,cellfun(@(x)x(2),... % moment due to pressure only
            arrayfun(@(a)-sum(cross(momentArm,pressureCoeffRarefied{a,2}) .* ...
            trianglesArea) / crossSectionalArea(referenceArea) / referenceLength, ...
            1:length(simAnglesOfAttack),'UniformOutput',false) ), ...
            styles{2},'LineWidth',1.25,'MarkerSize',10)
        plot(simAnglesOfAttack,cellfun(@(x)x(2),... % moment due to friction only
            arrayfun(@(a)-sum(cross(momentArm,frictionCoeffRarefied{a,2}) .* ...
            trianglesArea) / crossSectionalArea(referenceArea) / referenceLength, ...
            1:length(simAnglesOfAttack),'UniformOutput',false) ), ...
            styles{3},'LineWidth',1.25,'MarkerSize',10)
        hold off
        xlabel('Angle of Attack [deg]')
        ylabel([labels{5},' Coefficient [-]'])
        set(gca,'FontSize',15)
        grid on
        legend('Combined','Pressure','Shear','Location','NE')
        if saveFigure, saveas(F,['../../Report/figures/aero_rare_2d_',lower(labels{5}),'_split'],'epsc'), end
        
        %...Plot aerodynamic coefficients in 2D for continuum flow
        F = figure('rend','painters','pos',figSizeSmall);
        yyaxis left
        plot(simAnglesOfAttack,cellfun(@(x)x(1),aeroCoeffContinuum(:,1)),styles{1},'LineWidth',1.25,'MarkerSize',10)
        ylabel('Drag Coefficients [-]')
        yyaxis right
        plot(simAnglesOfAttack,cellfun(@(x)x(3),aeroCoeffContinuum(:,1)),styles{2},'LineWidth',1.25,'MarkerSize',10)
        ylabel('Lift Coefficient [-]')
        xlabel('Angle of Attack [deg]')
        set(gca,'FontSize',15)
        grid on
        legend('Drag','Lift','Location','NW')
        if saveFigure, saveas(F,'../../Report/figures/aero_cont_2d_force','epsc'), end
        
        F = figure('rend','painters','pos',figSizeSmall);
        hold on
        plot(simAnglesOfAttack,cellfun(@(x)x(5),aeroCoeffContinuum(:,1)),styles{1},'LineWidth',1.25,'MarkerSize',10)
        hold off
        xlabel('Angle of Attack [deg]')
        ylabel('Moment Coefficient [-]')
        set(gca,'FontSize',15)
        grid on
        if saveFigure, saveas(F,'../../Report/figures/aero_cont_2d_moment','epsc'), end
    end
    
    %...Plot aerodynamic coefficients in 3D against altitude
    for i = 1:2:6
        F = figure('rend','painters','pos',figSizeSmall);
        surf(simAnglesOfAttack,simAltTotal,cellfun(@(x)x(i),aeroCoeffTotal)')
        xlabel('Angle of Attack [deg]')
        ylabel('Altitude [km]')
        zlabel([labels{i},' Coefficient [-]'])
        set(gca,'FontSize',15)
        grid on
        if i == 5
            view([35,30])
        else
            view([-35,30])
        end
        if saveFigure, saveas(F,['../../Report/figures/aero_3d_',lower(labels{i})],'epsc'), end
    end
end

%% Validation

%...Get density
densityValues = interpolate(altitude,simAltRarefied,MCD.density);

%...Equation
dragDensity = @(d) 1.47952 - 0.032578 * log( d );

%...Plot
if showFigure
    F = figure('rend','painters','pos',figSizeSmall);
    figAxes1 = gca;
    hold on
    rectangle('Position',[1e-10,1.5,2e-7,3],'FaceColor',[0.5,0.5,0.5,0.5],'LineStyle','none','Parent',figAxes1)
    line(densityValues,cellfun(@(x)x(1),aeroCoeffRarefied(simAnglesOfAttack==0,:)),...
        'LineStyle','-','LineWidth',1.25,'Marker','o','MarkerSize',10,'Parent',figAxes1)
    line(densityValues,dragDensity(densityValues),...
        'LineStyle','--','LineWidth',1.25,'Color',[0.85,0.325,0.098],'Parent',figAxes1)
    hold off
    xlabel('Density [kg m^{-3}]')
    ylabel('Drag Coefficient [-]')
    ylim([1.5,3])
    set(gca,'FontSize',15,'XScale','log')
    grid on
    legend('Computed','Reference','Location','Best')
%     
%     figAxes2 = axes('Position',figAxes1.Position,'XAxisLocation','top','Color','none');
%     line(simAltRarefied,cellfun(@(x)x(1),aeroCoeffRarefied(simAnglesOfAttack==0,:)),...
%         'LineStyle','none','Parent',figAxes2)
%     xlabel('Altitude [km]')
%     xlim([85,550])
%     ylim([1.5,3])
%     set(gca,'FontSize',15,'XScale','log','XDir','Reverse')
    
    if saveFigure, saveas(F,['../../Report/figures/valid_gallis_diff'],'epsc'), end
end

%...Print percentage offset
fprintf( ['Offset to validation data:\n',repmat('%.1f\t',[1,length(densityValues)]),'\n'], ...
    ( cellfun(@(x)x(1),aeroCoeffRarefied(simAnglesOfAttack==0,:)) - dragDensity(densityValues) ) ./ ...
    dragDensity(densityValues) * 100)

%% Plot Bending Moments

if showFigure
    %...Plot bending moment magnitude for rarefied flow
    F = figure('rend','painters','pos',figSizeSmall);
    plot(simAltRarefied,bendingMomentRarefiedMagnitude(find(simAnglesOfAttack==0,1),:),'-o','LineWidth',1.25,'MarkerSize',10)
    xlabel('Altitude [km]')
    ylabel('Bending Moment [N m]')
    xticks(simAltRarefied)
    set(gca,'FontSize',15,'XScale','log','YScale','log')
    grid on
    if saveFigure, saveas(F,'../../Report/figures/bend_rare','epsc'), 
    else, title('Rarefied'); end
    
    %...Plot bending moment magnitude for rarefied flow at 100 km
    F = figure('rend','painters','pos',figSizeSmall);
    hold on
    for i = 1:3
        plot(simAnglesOfAttack,cellfun(@(x)x(i),bendingMomentRarefied(:,1)),styles{i},'LineWidth',1.25,'MarkerSize',10)
    end
    plot(simAnglesOfAttack,bendingMomentRarefiedMagnitude(:,1),styles{i+1},'LineWidth',1.25,'MarkerSize',10)
    hold off
    xlabel('Angle of Attack [deg]')
    ylabel('Bending Moment [N m]')
    set(gca,'FontSize',15)
    legend({'x_B axis','y_B axis','z_B axis','Magnitude'},'NumColumns',2,'Location','Best')%,'Location',[0.5,0.5,0.1,0.15])%
    grid on
    if saveFigure, saveas(F,'../../Report/figures/bend_rare_max','epsc'), 
    else, title('Rarefied at 100 km'); end
    
    %...Plot bending moment magnitude for continuum flow
    F = figure('rend','painters','pos',figSizeSmall);
    plot(simAnglesOfAttack,bendingMomentContinuumMagnitude','-o','LineWidth',1.25,'MarkerSize',10)
    xlabel('Angle of Attack [deg]')
    ylabel('Bending Moment [N m]')
    set(gca,'FontSize',15)
    grid on
    if saveFigure, saveas(F,'../../Report/figures/bend_cont','epsc'), 
    else, title('Rarefied'); end
end

%% Plot Pressure Distributions

%...Create triangulation
Tri.vertices = points;
Tri.faces = triangles;

%...Angle to show
% a = length(simAnglesOfAttack);
a = find(simAnglesOfAttack==0,1);
h = 2;

%...Plot
if showFigure
    %...Plot triangulation for rarefied flow
    titles = {'Pressure','Friction','Moment','Pressure','','Moment'};
    color = {sqrt(sum(pressureCoeffRarefied{a,h}.^2,2)),...
        sqrt(sum(frictionCoeffRarefied{a,h}.^2,2)),...
        sqrt(sum(cross(momentArm,pressureCoeffRarefied{a,h} + frictionCoeffRarefied{a,h}).^2,2)), ...
        sqrt(sum(pressureCoeffContinuum{a,1}.^2,2)),...
        NaN, ...
        sqrt(sum(cross(momentArm,pressureCoeffContinuum{a,1}).^2,2))};
    F = figure('rend','painters','pos',figSizeLarge);
    for i = 1:length(color)
        if i ~= 5
            subplot(2,3,i)
            Tri.facevertexcdata = color{i};
            patch(Tri), shading faceted, colormap jet
            c = colorbar; c.Location = 'southoutside';
            set(gca,'FontSize',15), view([127.5,30])
            axis off tight equal
            title(titles{i})
        end
    end
    subplotTitle([num2str(simAnglesOfAttack(a)),' deg, ',num2str(simAltRarefied(h)),' km'])
end
clear F j c

%...Vector Plots
a = find(simAnglesOfAttack==0,1);
if showFigure
    %...Plot pressure and shear forces for rarefied flow
    F = figure('rend','painters','pos',figSizeSmall);
    subplot(1,2,1)
    Tri.facevertexcdata = color{1};
    hold on
    patch(Tri,'FaceAlpha',0.5), shading faceted, colormap jet
    quiver3(trianglesCentroid(:,1),trianglesCentroid(:,2),trianglesCentroid(:,3),...
        pressureCoeffRarefied{a,h}(:,1),pressureCoeffRarefied{a,h}(:,2),pressureCoeffRarefied{a,h}(:,3))
    hold off
    c = colorbar; c.Location = 'southoutside';
    set(gca,'FontSize',15), view([127.5,30])
    set(gca,'FontSize',15)
    axis off tight equal
    title(titles{1})
    
    subplot(1,2,2)
    Tri.facevertexcdata = color{2};
    hold on
    patch(Tri,'FaceAlpha',0.5), shading faceted, colormap jet
    quiver3(trianglesCentroid(:,1),trianglesCentroid(:,2),trianglesCentroid(:,3),...
        frictionCoeffRarefied{a,h}(:,1),frictionCoeffRarefied{a,h}(:,2),frictionCoeffRarefied{a,h}(:,3))
    hold off
    c = colorbar; c.Location = 'southoutside';
    set(gca,'FontSize',15), view([127.5,30])
    set(gca,'FontSize',15)
    axis off tight equal
    title(titles{2})
end
subplotTitle([num2str(simAnglesOfAttack(a)),' deg, ',num2str(simAltRarefied(h)),' km'])

%% Compute Moment Coefficient Without Antenna

%...Disclaimer:
%   The way the moment coefficient without antenna is computed, is by
%   removing the pressure and shear coefficients of the panels
%   corresponding to the antenna. Thus, effects of the presence of the
%   antenna will not be removed, and the resulting coefficients do NOT
%   equal the ones that would result from a SPARTA simulation without the
%   antenna.

%...Take 125 km as reference
h = 2;

%...Loop over angles of attack
momentCoefficientNoAntenna = cell(1,length(simAnglesOfAttack));
for a = 1:length(simAnglesOfAttack)
    %...Retrieve pressure and shear coefficients (except antenna)
    forceValues = pressureCoeffRarefied{a,h}(~antennaElement.faces,:) + ...
        frictionCoeffRarefied{a,h}(~antennaElement.faces,:);
    
    %...Compute moment coefficient
    momentCoefficientNoAntenna{a} = - sum(cross(computeMomentArm(points,triangles(~antennaElement.faces,:),zeros(1,3)), ...
        forceValues) .* trianglesArea(~antennaElement.faces,:) ) / ...
        crossSectionalArea(referenceArea) / referenceLength;
end

%...Plot results
if showFigure
    %...Plot aerodynamic coefficients in 2D for rarefied flow
    F = figure('rend','painters','pos',figSizeSmall);
    hold on
    plot(simAnglesOfAttack,cellfun(@(x)x(5),aeroCoeffRarefied(:,h)),styles{1},'LineWidth',1.25,'MarkerSize',10)
    plot(simAnglesOfAttack,cellfun(@(x)x(2),momentCoefficientNoAntenna),styles{2},'LineWidth',1.25,'MarkerSize',10)
    hold off
    xlabel('Angle of Attack [deg]')
    ylabel('Y-Moment Coefficient [-]')
    legend('Nominal','No Antenna')
    set(gca,'FontSize',15)
    grid on
end

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
            GIFCommand = [GIFCommand,sprintf('%s ',fullfile(repository,'figures',hs,[as,'.image.*']))];
        end
        GIFCommand = [GIFCommand,fullfile(repository,'videos',['movie_',hs,'.gif'])];
    end
    
    %...Add command string to clipboard
    clipboard('copy',GIFCommand)
end

%% Save Text File With Coefficients

%...Save to file
if saveTable
    %...File names
    coefficient = {'MRODragCoefficients.txt','','MROLiftCoefficients.txt','','MROMomentCoefficients.txt',''};
    
    %...Loop over settings
    for i = 1:2:6
        %...Set file name and open
        fileName = ['/Users/Michele/Library/Mobile Documents/com~apple~CloudDocs/',...
            'University/Master Thesis/Code/MATLAB/data/',coefficient{i}];
        fileID = fopen(fileName,'w');
        
        %...Write number of independent variables
        fprintf(fileID,'%d\n',2);
        fprintf(fileID,'\n'); % separator
        
        %...Add independent variables
        fprintf(fileID,[repmat('%.10f\t ',[1,length(simAnglesOfAttack)]),'\n'],deg2rad(simAnglesOfAttack));
        fprintf(fileID,[repmat('%.3f\t ',[1,length(simAltitudes)]),'\n'],simAltitudes*1e3);
        fprintf(fileID,'\n'); % separator
        
        %...Add coefficients
        for j = 1:length(simAnglesOfAttack)
            fprintf(fileID,[repmat('%.6f\t ',[1,length(simAltitudes)]),'\n'],cellfun(@(x)x(i),aeroCoefficients(j,:)));
        end
        
        %...Close file
        fclose(fileID);
    end
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

function [momentArm,trianglesCentroid] = computeMomentArm(points,triangles,referencePoint)
    %...Function handle for surface normal
    centroid = @(p) [sum(p(:,1)),sum(p(:,2)),sum(p(:,3))]/3;

    %...Get triangle vertices
    vertices = arrayfun(@(i)points(triangles(i,:),:),1:size(triangles,1),'UniformOutput',false);

    %...Compute surface normal
    trianglesCentroid = cellfun(@(x)centroid(x),vertices,'UniformOutput',false)'; % compute centroid
    trianglesCentroid = cell2mat(trianglesCentroid);
    
    %...Find distance to center
    momentArm = trianglesCentroid - referencePoint;
end