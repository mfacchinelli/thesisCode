fclose('all'); clear all; close all force; profile off; clc; format long g; rng default;
addpath functions

%% Settings

%...Globals
global fullDatabase plotResults

%...Set settings
downloadDatabse = false;
fullDatabase = true;
if fullDatabase, folder = 'MCDFull';
else, folder = 'MCD'; end

genTable = true;
plotResults = false;
savePlots = false;
[figSizeLarge,figSizeSmall] = saveFigureSettings(savePlots);

%% Database Parameters

%...Constants
global Ru
Ru = 8.3145; % universal gas constant
Na = 6.0221e23; % Avogadro's number
masses = [44.01,39.95,28.01,28.01,48.00,16.00,32.00,1.01,2.02,4.00,18.02]*1e-3; % molar masses
collision = [3.941,3.542,3.798,3.760,3.467,3.467,3.467,2.827,2.827,2.551,2.650]*1e-10; % colision diameters

%...Time range
ts = 0:90:270; T = length(ts);

%...Altitude range
[H,hs,hs_plot,hs_loc_plot] = altitudeRange();
if plotResults
    figure; 
    yyaxis left
    scatter(1:length(hs),hs)
    set(gca,'YScale','log')
    yyaxis right
    scatter(1:length(hs)-1,diff(hs))
    set(gca,'YScale','log')
    grid on
end

%...Latitude and longitude
global latitude longitude altitude
Lat = 48; Lon = 64; N = T*H;
latitude = linspace(-90,90,Lat);
longitude = linspace(-180,180,Lon);
altitude = hs;

%...Pre-allocate variables
if fullDatabase
    dataStr = {'Density','Pressure','Temperature','Horiz. Wind','Vert. Wind','Specific Heat',...
        'CO_2','Ar','N_2','CO','O_3','O','O_2','H','H_2','He','H_2O'};
    dataLabel = {'dens','pres','temp','horz','vert','cp',...
        'co2','ar','n2','co','o3','o','o2','h','h2','he','h2o'};
    dataUnit = [{'log_{10}(kg m^{-3})','log_{10}(Pa)','K','m s^{-1}','m s^{-1}','J kg^{-1} K^{-1}'},...
        repmat({'-'},[1,11])];
else
    dataStr = {'Density','Pressure','Temperature','Horiz. Wind','Vert. Wind','Specific Heat'};
    dataLabel = {'dens','pres','temp','horz','vert','cp'};
    dataUnit = {'log_{10}(kg m^{-3})','log_{10}(Pa)','K','m s^{-1}','m s^{-1}','J kg^{-1} K^{-1}'};
end
D = length(dataStr);

%...Locators
global densLoc presLoc tempLoc cpLoc
densLoc = find(cellfun(@(x)isequal(x,'dens'),dataLabel));
presLoc = find(cellfun(@(x)isequal(x,'pres'),dataLabel));
tempLoc = find(cellfun(@(x)isequal(x,'temp'),dataLabel));
windLoc = find(cellfun(@(x)isequal(x,'horz'),dataLabel));
cpLoc = find(cellfun(@(x)isequal(x,'cp'),dataLabel));
gassesLoc = find(cellfun(@(x)length(x)<=4,dataStr));
if fullDatabase
    gasOrder = [2,4,1,8,9,11,10,3,6,7,5];
    gasStr = dataStr(gasOrder+min(gassesLoc)-1);
    gassesLoc = gassesLoc(gasOrder);
    masses = masses(gasOrder);
    collision = collision(gasOrder);
end
G = length(gassesLoc);

%% Download Data

if downloadDatabse
    delete(fullfile(folder,'*.txt'))
    tic
    GenerateMCDDatabase(ts,hs,D,folder)
    toc
end

%% Open Files

%...MCD format
formatSpec = ['%f ||',repmat('%f ',[1,Lat])];

%...Read files
i = 0;
files = cell(H,T);
for t = ts
    for h = hs
        i = i+1;
        fileID = fopen(fullfile(folder,[num2str(t),'_',num2str(h),'.txt']));
        tempData = [];
        for d = 1:D+1
            tempData = vertcat(tempData,cell2mat(textscan(fileID,formatSpec,...
                'CommentStyle','#','Delimiter','\t','CollectOutput',true,'HeaderLines',2)));
        end
        files{i} = tempData;
    end
end
clear tempData

%% Extract Data

%...Distribute
data = repmat({zeros(Lat,Lon,H,T)},[1,D]);
for t = 1:T
    for h = 1:H
        for d = 1:D
            data{d}(:,:,h,t) = files{h,t}((1+Lon*(d-1)):Lon*d,2:end)';
            if d == cpLoc % remove NaNs in Cp 
                locNaN = isnan(data{d}(:,:,h,t));
                if any(any(locNaN))
                    cpData = data{d}(:,:,h,t);
                    cpData(locNaN) = mean(mean(cpData(~locNaN)));
                    data{d}(:,:,h,t) = cpData;
                end
            end
        end
    end
end
clear locNaN cpData

%% Molecular Properties

%...Compute molecular mass and ratio
molarRatioLatLonTimeAvg = cell(size(gassesLoc));
molarRatio = 0;
molarMass = 0;
for g = gassesLoc
    molarRatioLatLonTimeAvg{g-min(gassesLoc)+1} = squeeze(mean(mean(mean(data{g},1),2),4));
    molarRatio = molarRatio + data{g};
    molarMass = molarMass + masses(g-6)*data{g};
end
data_extended = horzcat(data,molarMass);

%% Create Tabulated Atmosphere

if genTable
    %...Average over latitude, longitude, time, to get tabulated atmosphere
    %   over altitude only
    tabularLatLonTimeAvg = TabulatedMCDAtmosphere(data_extended,'f(h)');
    
    %...Average over time, to get tabulated atmosphere over altitude, 
    %   latitude and longitude
%     tabularTimeAvg = TabulatedMCDAtmosphere(data_extended,'f(h,d,l)');
end

%% Atmospheric Composition

%...Values
L_char = 2.5; % characteristic length
numberDensityFun = @(M,rho) Na*rho./M;
meanFreePathFun = @(M,rho,sigma) 1/sqrt(2)/pi./numberDensityFun(M,rho)./sigma.^2;
knudsenNumberFun = @(M,rho,sigma) meanFreePathFun(M,rho,sigma)/L_char;

%...Compute number density per altitude
if genTable && fullDatabase
    %...Determine gas mixture values
    gasRatio = zeros(length(hs),G);
    moleMass = zeros(size(hs));
    collisionDiameter = zeros(size(hs));
    density = tabularLatLonTimeAvg(:,1);
    pressure = tabularLatLonTimeAvg(:,2);
    temperature = tabularLatLonTimeAvg(:,3);
    numberDensity = zeros(size(hs));
    knudsenNumber = zeros(size(hs));
    for h = 1:length(hs)
        gasRatio(h,:) = arrayfun(@(g)molarRatioLatLonTimeAvg{g-min(gassesLoc)+1}(h),gassesLoc);
        moleMass(h) = sum(masses.*gasRatio(h,:));
        collisionDiameter(h) = sum(collision.*gasRatio(h,:));
        numberDensity(h) = numberDensityFun(moleMass(h),density(h));
        knudsenNumber(h) = knudsenNumberFun(moleMass(h),density(h),collisionDiameter(h));
    end
    
    %...Combine elements for SPARTA analysis
    R = Ru./squeeze(mean(mean(mean(molarMass,1),2),4)); % gas constant
    gamma = 1./(1-R./squeeze(mean(mean(mean(data{cpLoc},1),2),4))); % specific heat ratio
    sos = sqrt(R.*gamma.*tabularLatLonTimeAvg(:,3)); % speed of sound
    save('../SPARTA/MCDEnvir','altitude','density','pressure','temperature','gamma',...
        'gasRatio','gasStr','knudsenNumber','numberDensity','sos');
    
    if plotResults
        %...Plot gas percentage
        F = figure('rend','painters','pos',figSizeSmall);
        styles = repmat({'-',':','--','-.'},[1,3]);
        hold on
        for g = 1:G
            plot(hs,gasRatio(:,g),'LineWidth',1.75,'LineStyle',styles{g})
        end
        hold off
        xlabel('Altitude [km]')
        ylabel('Gas Fraction [-]')
        legend(gasStr,'Location','Best')
        set(gca,'FontSize',15,'XScale','log','YScale','log')
        grid on
        xlim([hs_plot(1),hs_plot(end)]), xticks([50,100,250,500,1500])
        ylim([1e-10,1])
        if ~savePlots, title('Atmospheric Composition'),
        else, saveas(F,'../../Report/figures/mars_atm_comp','epsc'), end
        
        %...Plot molar mass and collision diameter
        F = figure('rend','painters','pos',figSizeSmall);
        yyaxis left
        hold on
        plot(hs,moleMass*1e3,'LineWidth',1.75)
        plot([10,10^6],masses(3)*1e3*ones(1,2),'LineWidth',1.25,'LineStyle','--')
        plot([10,10^6],masses(4)*1e3*ones(1,2),'LineWidth',1.25,'LineStyle','-.')
        hold off
        ylabel('Molar Mass [g mol^{-1}]')
        set(gca,'FontSize',15,'XScale','log')
        yyaxis right
        plot(hs,abs(sum(gasRatio,2)-1),'LineWidth',1.75)
        ylabel('Offset in Gas Presence [-]')
        set(gca,'FontSize',15,'XScale','log','YScale','log')
        xlabel('Altitude [km]')
        legend('Mars','CO_2','H','Gas','Location','Best')
        grid on
        xlim([hs_plot(1),hs_plot(end)]), xticks([50,100,250,500,1500])
        if ~savePlots, title('Atmospheric Composition'),
        else, saveas(F,'../../Report/figures/mars_atm_mol','epsc'), end
        
        %...Plot number density and Knudsen
        F = figure('rend','painters','pos',figSizeSmall);
        yyaxis left
        loglog(hs,numberDensity,'LineWidth',1.75)
        ylabel('Number Density [-]')
        yyaxis right
        loglog(hs,knudsenNumber,'LineWidth',1.75)
        ylabel('Knudsen Number [-]')
        xlabel('Altitude [km]')
        set(gca,'FontSize',15)
        grid on
        xlim([hs_plot(1),hs_plot(end)]), xticks([50,100,250,500,1500])
        if ~savePlots, title('Number Density and Knudsen'),
        else, saveas(F,'../../Report/figures/mars_atm_num','epsc'), end
    end
end

%% Plot

if plotResults
    %...Make density and pressure logarithmic
    data{densLoc} = log10(data{densLoc});
    data{presLoc} = log10(data{presLoc});

    %...Mesh grid
    [lonH,latH,alt] = meshgrid(longitude,latitude,hs);
    [lonT,latT,time] = meshgrid(longitude,latitude,ts);
    
    %...MCD Data
    for i = 1:D-G
        if i == densLoc || i == presLoc % density and pressure
            F = figure('rend','painters','pos',figSizeLarge);
            for h = 1:length(hs_plot)
                subplot(3,2,h)
                S = slice(lonT,latT,time,squeeze(data{i}(:,:,hs_loc_plot(h),:)),[],[],ts);
                plotSettings(S,dataUnit{i},ts,hs_plot(h),false,true)
            end
            if ~savePlots, subplotTitle(dataStr{i}),
            else, saveas(F,['../../Report/figures/mars_',dataLabel{i}],'epsc'), end
        else % other data
            F = figure('rend','painters','pos',figSizeLarge);
            for t = 1:T
                subplot(2,2,t)
                S = slice(lonH,latH,alt,data{i}(:,:,:,t),[],[],hs_plot);
                plotSettings(S,dataUnit{i},hs_plot,ts(t),true,true)
            end
            if ~savePlots, subplotTitle(dataStr{i}),
            else, saveas(F,['../../Report/figures/mars_',dataLabel{i}],'epsc'), end
        end
    end
    
    %...Gas constant
    if fullDatabase
%         R = Ru./molarMass;
        R = 197*ones(1,1,1,T);
    else
        R = 197*ones(1,1,1,T);
    end
    gamma = 1./(1-R./data{cpLoc});

    %...Speed of sound
    F = figure('rend','painters','pos',figSizeLarge);
    for t = 1:T
        subplot(2,2,t)
        sos = sqrt(R(:,:,:,t).*gamma(:,:,:,t).*data{tempLoc}(:,:,:,t));
        S = slice(lonH,latH,alt,sos,[],[],hs_plot);
        plotSettings(S,dataUnit{windLoc},hs_plot,ts(t),true,true)
    end
    if ~savePlots, subplotTitle('Speed of Sound'),
    else, saveas(F,'../../Report/figures/mars_sos','epsc'), end

    %...Mean values
    if genTable
        mean_units = {'Pa','kg m^{-3}','K','J kg^{-1} K^{-1}',''};
        mean_titles = {'Pressure','Density','Temperature','Gas Constant','Specific Heat Ratio'};
        mean_labels = {'mars_pres_mean','mars_dens_mean','mars_temp_mean','mars_gas_const_mean','mars_heat_ratio_mean'};
        for i = 1:size(tabularTimeAvg,4)
            F = figure('rend','painters','pos',figSizeLarge);
            for h = 1:length(hs_plot)
                subplot(2,3,h)
                S = pcolor(lonH(:,:,1),latH(:,:,1),tabularTimeAvg(:,:,hs_loc_plot(h),i));
                plotSettings(S,mean_units{i},[],hs_plot(h),false,false)
            end
            if ~savePlots, subplotTitle(['Average ',mean_titles{i}]),
            else, saveas(F,['../../Report/figures/',mean_labels{i}],'epsc'), end
        end
    end
end

%% Velocities

%...General data
r = 3.390e6;
mu = 4.282e13;
o = 2*pi/(24.6597*3600);

%...Change altitude range
hs_tab = hs(hs<=500); hs_tab = [hs_tab(1:15:end),hs_tab(end)];
[~,h_loc_tab] = ismembertol(hs_tab,hs,1e-6);

%...Velocities
i = 0;
for h = 1:length(hs_tab)
    h = h_loc_tab(h); i = i+1;
    V_wind_max(i) = max(max(max(data{windLoc}(:,:,h,:))));
    V_wind_rms(i) = rms(rms(rms(data{windLoc}(:,:,h,:))));
    V_atm(i) = o*(r+hs(h)*1e3);
    V_circ(i) = sqrt(mu/(r+hs(h)*1e3));
    
    angle_max(i) = round(atan2d(V_atm(i)-V_wind_max(i),V_circ(i)),1);
    angle_rms(i) = round(atan2d(V_atm(i)-V_wind_rms(i),V_circ(i)),1);
    angle_no_wind(i) = round(atan2d(V_atm(i),V_circ(i)),1);
    mag_max(i) = round(sqrt(V_circ(i)^2+(V_atm(i)-V_wind_max(i))^2),0);
    mag_mrs(i) = round(sqrt(V_circ(i)^2+(V_atm(i)-V_wind_rms(i))^2),0);
    mag_no_wind(i) = round(sqrt(V_circ(i)^2+V_atm(i)^2),0);
end
table(hs_tab',round(V_circ',-1),round(V_atm',0),round(V_wind_max',0),round(V_wind_rms',0))
table(hs_tab',angle_max',angle_rms',angle_no_wind',mag_max',mag_mrs',mag_no_wind')

fulldata = vertcat(round(hs_tab,-1),round(V_circ,-1),round(V_atm,0),round(V_wind_max,0),round(V_wind_rms,0));
sprintf([repmat('\\num{%.0f} & ',[1,4]),'\\num{%0.f} \\\\\n'],fulldata)

%% TEST >>>>>>>>>>>

% R = Ru./molecular_mass;
% gamma = 1./(1-R./data{cp_loc});
% 
% mean(mean(mean(mean(molecular_ratio))))
% median(median(median(median(molecular_ratio))))
% max(max(max(max(molecular_ratio))))
% min(min(min(min(molecular_ratio))))
% 
% figure;
% i = 0;
% for t = 1:T
%     for h = 1:length(hs_plot)%1:H%4:H%
%         i = i+1;
%         subplot(4,6,i)%subplot(4,3,i)%
%         pcolor(gamma(:,:,hs_loc_plot(h),t)), colorbar, title(sprintf(...
%             '%i deg, %i km',ts(t),hs_plot(h))), axis tight
%     end
% end
% 
% mean(mean(mean(mean(R))))
% median(median(median(median(R))))
% max(max(max(max(R))))
% min(min(min(min(R))))

%%% TEST <<<<<<<<<<
%% Functions

function [H,hs,hs_plot,hs_loc] = altitudeRange()
    hs = []; h_bound = [50,10000]; h_steps = 150;
    for i = 1:length(h_bound)-1
        hs = horzcat(hs,logspace(log10(h_bound(i)),log10(h_bound(i+1)),h_steps(i)));
    end
    hs = unique(round(hs)); H = length(hs);
    hs_plot = [50,102,151,248,504,1519]; [~,hs_loc] = ismembertol(hs_plot,hs,1e-6);
end

function plotSettings(S,label,zs,w,log,three)
    global save_plots
    c = colorbar;
    c.Label.String = label;
    xlabel('Long. [deg]'), xticks(-180:90:180)
    ylabel('Lat. [deg]'), yticks(-90:45:90)
    if three
        view(21,21)
        set(S,'EdgeColor','none')
%         if save_plots
%             set(S,'EdgeColor','none')
%         else
%             set(S,'EdgeColor','none','FaceColor','interp','FaceAlpha','interp')
%             alpha('color')
%             alphamap('rampdown')
%             alphamap('increase',0.25)
%         end
    else
        set(S,'EdgeColor','none')
    end
    colormap jet
    if all( w > 1000 ), decimal = -2;
    else, decimal = -1; end
    if log
        zlabel('Alt. [km]'), zticks(zs)
        set(gca,'zscale','log','FontSize',15)
        title([num2str(round(w,decimal)),' deg'])
    else
        zlabel('Solar Long. [deg]'), zticks(zs)
        set(gca,'FontSize',15)
        title([num2str(round(w,decimal)),' km'])
    end
    axis tight
end