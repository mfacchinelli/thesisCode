fclose('all'); clear all; close all force; profile off; clc; format long g; rng default;
addpath functions tests

%% Settings

%...Figures and tables setting
showFigure = true;
saveFigure = true;
[figSizeLarge,figSizeMedium,figSizeSmall] = saveFigureSettings(saveFigure);

%...Preallocate variables
allAeroCoeff = cell(2,3);

%% Collect All Data

%...Loop over each test case
for testCase = 1:3
    %% SPARTA and Modeling Variables
    
    %...Select test case
    %       1: both real2sim and gridSpacing
    %       2: only real2sim
    %       3: only gridSpacing
    
    %...Repositories
    repositoryLocal = 'nominal';
    repositoryServer = 'nominal_eudx';
    SPARTARepository = ['/Users/Michele/Library/Mobile Documents/com~apple~CloudDocs/University/Master ',...
        'Thesis/Code/SPARTA/']; % path to SPARTA folder
    MRORepository = fullfile(SPARTARepository,'mro'); % path to MRO folder
    OutputRepositoryLocal = fullfile(MRORepository,repositoryLocal,'data'); % path to SPARTA local results
    switch testCase % path to SPARTA server results
        case 1
            OutputRepositoryServer = fullfile(MRORepository,repositoryServer,'data_both');
        case 2
            OutputRepositoryServer = fullfile(MRORepository,repositoryServer,'data_ratio');
        case 3
            OutputRepositoryServer = fullfile(MRORepository,repositoryServer,'data_grid');
    end
    
    %...Files
    MRODataFile = fullfile(MRORepository,'data.mro'); % path to MRO model file
    
    %...Mars Climate Database
    MCDData = fullfile(SPARTARepository,'MCDEnvir');
    
    %% Spacecraft Model
    
    %...Compute model based on MRO discretization
    [MROExtent,points,triangles,solarPanelElements] = ...
        SPARTAModel(MRODataFile,false,false,false); % output min and max dimensions of MRO
    
    %% Constants
    
    %...Simulation conditions
    simAnglesOfAttack = -30:5:30; % angles of attack for simulation
    %[linspace(-75,-30,4),linspace(-25,25,11),linspace(30,75,4)];
    simAltitudes = 125; % altitudes for rarefied regime simulations
    
    %% Run SPARTA
    
    %...Folder names for SPARTA output
    dataFolderLocal = cell(size(simAltitudes));
    dataFolderServer = cell(size(simAltitudes));
    
    %...Add folders to list
    for h = 1:length(simAltitudes)
        dataFolderLocal{h} = fullfile(OutputRepositoryLocal,num2str(simAltitudes(h)));
        dataFolderServer{h} = fullfile(OutputRepositoryServer,num2str(simAltitudes(h)));
    end
    
    %% Mars Climate Database Model
    
    %...Martian environment
    MarsGravityParameter = 4.282e13;
    MarsRadius = 3.390e6;
    
    %...Load MCD parameters
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
    [density,pressure,temperature,gamma,speedOfSound,gasRatios,numberDensity] = interpolate(altitude,simAltitudes,...
        density,pressure,temperature,gamma,speedOfSound,gasRatios,numberDensity);
    
    %...Find circular velocity at altitude
    MarsCircularVelocity = @(altitude) sqrt( MarsGravityParameter ./ ( MarsRadius + altitude * 1e3 ) );
    streamVelocity = zeros(length(simAltitudes),3);
    streamVelocity(:,1) = - MarsCircularVelocity(simAltitudes)';
    
    %% Area Analysis
    
    %...Compute area and normal to surface of triangles and surface areas
    %   around each axis for whole spacecraft
    [trianglesArea,trianglesNormal,crossSectionalArea] = computeTriangleAreaNormal(points,triangles);
    referenceArea = 1;
    referenceLength = 2.5;
    centerOfMass = [0,0,0.1375];
    [momentArm,trianglesCentroid] = computeMomentArm(points,triangles,centerOfMass);
    
    %% Analyze SPARTA Rarefied Results Local
    
    %...Get results for rarefied flow
    pressureCoeffLocal = cell(length(simAnglesOfAttack),length(simAltitudes));
    frictionCoeffLocal = cell(length(simAnglesOfAttack),length(simAltitudes));
    aeroCoeffLocal = cell(length(simAnglesOfAttack),length(simAltitudes));
    for h = 1:length(simAltitudes)
        for a = 1:length(simAnglesOfAttack)
            %...Get (average) distribution
            distribution = cell(size(simAltitudes));
            files = dir(fullfile(dataFolderLocal{h},[num2str(simAnglesOfAttack(a)),'.coeff.*']));
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
            
            %...Transformation from aerodynamic to body frame
            transformation = roty(simAnglesOfAttack(a));
            
            %...Compute pressure and friction coefficients
            dynamicPressure = 1/2 * density(h) * norm(streamVelocity(h,:))^2;
            pressureCoeffLocal{a,h} = (averageDistribution(:,1:3) - pressure(h) * ...
                ( transformation' * trianglesNormal' )' ) / dynamicPressure; % rotate normal to account for angle of attack
            frictionCoeffLocal{a,h} = averageDistribution(:,4:6) / dynamicPressure;
            
            %...Integrate to find force coefficients
            forceCoefficients = sum((pressureCoeffLocal{a,h} + frictionCoeffLocal{a,h}) .* ...
                trianglesArea) / crossSectionalArea(referenceArea);
            
            %...Integrate to find moment coefficients
            momentArmInAerodynamicFrame = arrayfun(@(i)transformation' * ...
                momentArm(i,:)',1:size(momentArm,1),'UniformOutput',false);
            momentCoefficients = sum(cross([momentArmInAerodynamicFrame{:}]', ...
                pressureCoeffLocal{a,h} + frictionCoeffLocal{a,h}) .* ...
                trianglesArea) / crossSectionalArea(referenceArea) / referenceLength;
            
            %...Find side, drag and lift coefficients
            %   Note that a negative sign is added since aerodynamic forces are
            %   positive in the negative direction
            aeroCoeffLocal{a,h} = - [forceCoefficients,momentCoefficients]';
        end
    end
    allAeroCoeff{1,testCase} = aeroCoeffLocal;
    clear distribution files fileID result averageDistribution transformation solarPanelDistribution ...
        solarPanelDistributionInAerodynamicFrame dynamicPressure forceCoefficients momentCoefficients
    
    %% Analyze SPARTA Rarefied Results Server
    
    %...Get results for rarefied flow
    pressureCoeffServer = cell(length(simAnglesOfAttack),length(simAltitudes));
    frictionCoeffServer = cell(length(simAnglesOfAttack),length(simAltitudes));
    aeroCoeffServer = cell(length(simAnglesOfAttack),length(simAltitudes));
    for h = 1:length(simAltitudes)
        for a = 1:length(simAnglesOfAttack)
            %...Get (average) distribution
            distribution = cell(size(simAltitudes));
            files = dir(fullfile(dataFolderServer{h},[num2str(simAnglesOfAttack(a)),'.coeff.*']));
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
            
            %...Transformation from aerodynamic to body frame
            transformation = roty(simAnglesOfAttack(a));
            
            %...Compute pressure and friction coefficients
            dynamicPressure = 1/2 * density(h) * norm(streamVelocity(h,:))^2;
            pressureCoeffServer{a,h} = (averageDistribution(:,1:3) - pressure(h) * ...
                ( transformation' * trianglesNormal' )' ) / dynamicPressure; % rotate normal to account for angle of attack
            frictionCoeffServer{a,h} = averageDistribution(:,4:6) / dynamicPressure;
            
            %...Integrate to find force coefficients
            forceCoefficients = sum((pressureCoeffServer{a,h} + frictionCoeffServer{a,h}) .* ...
                trianglesArea) / crossSectionalArea(referenceArea);
            
            %...Integrate to find moment coefficients
            momentArmInAerodynamicFrame = arrayfun(@(i)transformation' * ...
                momentArm(i,:)',1:size(momentArm,1),'UniformOutput',false);
            momentCoefficients = sum(cross([momentArmInAerodynamicFrame{:}]', ...
                pressureCoeffServer{a,h} + frictionCoeffServer{a,h}) .* ...
                trianglesArea) / crossSectionalArea(referenceArea) / referenceLength;
            
            %...Find side, drag and lift coefficients
            %   Note that a negative sign is added since aerodynamic forces are
            %   positive in the negative direction
            aeroCoeffServer{a,h} = - [forceCoefficients,momentCoefficients]';
        end
    end
    allAeroCoeff{2,testCase} = aeroCoeffServer;
    clear distribution files fileID result averageDistribution transformation solarPanelDistribution ...
        solarPanelDistributionInAerodynamicFrame dynamicPressure forceCoefficients momentCoefficients
    
    %% Plot Aerodynamic Coefficients
    
    %...Plot
    if showFigure
        %...Plot aerodynamic coefficients in 2D for rarefied flow
        styles = {'-o','-d','-s','-v','-p','-h','-*','-x','-^','-o','-d'};
        labels = {'Drag','Side','Lift','X-Moment','Y-Moment','Z-Moment'};
        locations = {'NE','','SE','','NE','NW'};
        for i = [1:2:5,6]
            F = figure('rend','painters','pos',figSizeSmall);
            hold on
            for h = 1:length(simAltitudes)
                plot(simAnglesOfAttack,cellfun(@(x)x(i),aeroCoeffLocal(:,h)), ...
                    styles{h},'LineWidth',1.25,'MarkerSize',10)
                plot(simAnglesOfAttack,cellfun(@(x)x(i),aeroCoeffServer(:,h)),...
                    styles{h+1},'LineWidth',1.25,'MarkerSize',10)
            end
            hold off
            xlabel('Angle of Attack [deg]')
            ylabel([labels{i},' Coefficient [-]'])
            if i == 6
                ylim([-1e-2,1e-2])
            end
            legend('Nominal','Accurate','Location',locations{i})
            set(gca,'FontSize',15)
            grid on
            if saveFigure, saveas(F,['../../Report/figures/aero_server_',num2str(testCase),...
                    '_',lower(labels{i})],'epsc'), end
        end
    end
    
    %% Plot Pressure Distributions
    
    %...Create triangulation
    Tri.vertices = points;
    Tri.faces = triangles;
    
    %...Angle to show
    angles = [1,find(simAnglesOfAttack==0,1)];
%     pressure = pressureCoeffServer;
%     friction = frictionCoeffServer;
    pressure = pressureCoeffLocal;
    friction = frictionCoeffLocal;
    if saveFigure
        pressure = pressureCoeffLocal;
        friction = frictionCoeffLocal;
    end
    h = 1;
    
    %...Plot
    if showFigure
        %...Plot triangulation for rarefied flow
        titles = {'Pressure','Friction','Moment'};
        views = {[90,0],[180,0],[90,90]};
        for a = angles
            color = {sqrt(sum(pressure{a,h}.^2,2)) .* trianglesArea / crossSectionalArea(referenceArea),...
                sqrt(sum(friction{a,h}.^2,2)) .* trianglesArea / crossSectionalArea(referenceArea),...
                sqrt(sum(cross(momentArm,pressure{a,h} + friction{a,h}).^2,2)) .* ...
                trianglesArea / crossSectionalArea(referenceArea) / referenceLength};
%             color = {sqrt(sum(pressure{a,h}.^2,2)),...
%                 sqrt(sum(friction{a,h}.^2,2)),...
%                 sqrt(sum(cross(momentArm,pressure{a,h} + friction{a,h}).^2,2))};
            F = figure('rend','painters','pos',figSizeLarge);
            for i = 1:length(color)
                for j = 1:3
                    subplot(length(color),3,(i-1)*3 + j)
                    Tri.facevertexcdata = color{i};
                    patch(Tri,'LineStyle','none'), shading faceted, colormap jet
                    c = colorbar; c.Location = 'southoutside';
                    set(gca,'FontSize',15), view(views{j})
                    axis off tight equal
                end
            end
            if saveFigure, saveas(F,['../../Report/figures/aero_color_',num2str(a)],'epsc'),
            else, subplotTitle([num2str(simAnglesOfAttack(a)),' deg']), end
        end
    end
    clear F j color c
end

%% Close All Figures

if saveFigure, close all, end

%% Plot All Results

%...Plot
if showFigure
    %...Plot aerodynamic coefficients in 2D for rarefied flow
    styles = {'-o','-d','-s','-v','-p','-h','-*','-x','-^','-o','-d'};
    labels = {'Drag','Side','Lift','X-Moment','Y-Moment','Z-Moment'};
    locations = {'NE','','SE','','NE','NW'};
    for i = [1:2:5,6]
        F = figure('rend','painters','pos',figSizeSmall);
        hold on
        for h = 1:length(simAltitudes)
            plot(simAnglesOfAttack,cellfun(@(x)x(i),allAeroCoeff{1,1}(:,h)), ...
                styles{h},'LineWidth',1.25,'MarkerSize',10)
            for t = 1:3
                plot(simAnglesOfAttack,cellfun(@(x)x(i),allAeroCoeff{2,t}(:,h)),...
                    styles{h+1},'LineWidth',1.25,'MarkerSize',10)
            end
        end
        hold off
        xlabel('Angle of Attack [deg]')
        ylabel([labels{i},' Coefficient [-]'])
        if i == 6
            ylim([-1e-2,1e-2])
        end
        legend('Nominal','Accurate 1','Accurate 2','Accurate 3','Location',locations{i})
        set(gca,'FontSize',15)
        grid on
        if saveFigure, saveas(F,['../../Report/figures/aero_server_all',...
                '_',lower(labels{i})],'epsc'), end
    end
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