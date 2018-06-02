fclose('all'); clear all; close all force; profile off; clc; format long g; rng default;
addpath functions

%% Settings

%...Globals
global fullDatabase showFigures saveFigures

%...Set settings
downloadDatabse = false;
fullDatabase = true;
if fullDatabase, folder = 'MCDFull';
else, folder = 'MCD'; end

genTable = true;
showFigures = true;
saveFigures = true;
[figSizeLarge,figSizeMedium] = saveFigureSettings(saveFigures);

%% Database Parameters

%...Constants
global Ru
Ru = 8.3144598; % universal gas constant
Na = 6.022140857e23; % Avogadro's number
masses = [44.01,39.95,28.01,28.01,48.00,16.00,32.00,1.01,2.02,4.00,18.02]*1e-3; % molar masses
collision = [3.941,3.542,3.798,3.760,3.467,3.467,3.467,2.827,2.827,2.551,2.650]*1e-10; % colision diameters

%...Time range
ts = 0:90:270; T = length(ts);

%...Altitude range
[H,hs,hs_plot,hs_loc_plot] = altitudeRange();
if showFigures
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
gasesLoc = find(cellfun(@(x)length(x)<=4,dataStr));
G = length(gasesLoc);

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

%...Reorder gases
if fullDatabase
    gasOrder = [2,4,1,8,9,11,10,3,6,7,5];
    dataStr(min(gasesLoc):end) = dataStr(gasOrder+min(gasesLoc)-1);
    dataLabel(min(gasesLoc):end) = dataLabel(gasOrder+min(gasesLoc)-1);
    gasStr = dataStr(min(gasesLoc):end);
    gasesLoc = gasesLoc(gasOrder);
    masses = masses(gasOrder);
    collision = collision(gasOrder);
    oldData = data;
    for g = 0:G-1
        data{min(gasesLoc) + g} = oldData{min(gasesLoc) + gasOrder(g+1) - 1};
    end
    gasesLoc = find(cellfun(@(x)length(x)<=4,dataStr)); % redefine
end
HLoc = find(cellfun(@(x)isequal(x,'H'),dataStr));

%% Latitude-Longitude-Time-Average Atmosphere

%...Function handles for functions
L_char = 2.5; % characteristic length
numberDensityFun = @(M,rho) Na./M.*rho;
meanFreePathFun = @(M,rho,sigma) 1/sqrt(2)/pi./numberDensityFun(M,rho)./sigma.^2;
knudsenNumberFun = @(M,rho,sigma) meanFreePathFun(M,rho,sigma)/L_char;

%...Average data over latitude, longitude and time
dataLatLonTimeAvg = zeros(size(data{1},3),length(data));
for i = 1:length(data)
    dataLatLonTimeAvg(:,i) = mean(mean(mean(data{i},1),2),4);
end

%...Compute information on molecular mass, gas constant, specific heat
%   ratio and speed of sound
molarMassLatLonTimeAvg = dataLatLonTimeAvg(:,gasesLoc) * masses';
gasConstantLatLonTimeAvg = Ru./molarMassLatLonTimeAvg;
specificHeatRatioLatLonTimeAvg = 1./(1-gasConstantLatLonTimeAvg./dataLatLonTimeAvg(:,cpLoc));
specificHeatRatioLatLonTimeAvg(dataLatLonTimeAvg(:,HLoc) > 0.075) = 5/3; % correct values, since MCD gives weird results
speedOfSoundLatLonTimeAvg = sqrt(specificHeatRatioLatLonTimeAvg .* gasConstantLatLonTimeAvg .* ...
    dataLatLonTimeAvg(:,tempLoc));

%...Add new data to atmosphere table
dataLatLonTimeAvg = horzcat(dataLatLonTimeAvg,gasConstantLatLonTimeAvg,...
    specificHeatRatioLatLonTimeAvg,molarMassLatLonTimeAvg);

%...Compute other properties
numberDensityLatLonTimeAvg = numberDensityFun(molarMassLatLonTimeAvg,dataLatLonTimeAvg(:,densLoc));
collisionDiameterLatLonTimeAvg = dataLatLonTimeAvg(:,gasesLoc) * collision';
knudsenNumberLatLonTimeAvg = knudsenNumberFun(molarMassLatLonTimeAvg,dataLatLonTimeAvg(:,densLoc),...
    collisionDiameterLatLonTimeAvg);

%...Verify results
if showFigures && ~saveFigures
    figure
    yyaxis left
    hold on
    plot(altitude,gasConstantLatLonTimeAvg,'LineWidth',1.25)
    plot(altitude,dataLatLonTimeAvg(:,cpLoc),'LineWidth',1.25)
    hold off
    yyaxis right
    plot(altitude,specificHeatRatioLatLonTimeAvg,'LineWidth',1.25)
    xlabel('Altitude [km]')
    set(gca,'FontSize',12.5,'XScale','log')
    grid on
    xlim([altitude(1),altitude(end)])
end

%% Time-Average Atmosphere

%...Average data over latitude, longitude and time
dataTimeAvg = zeros(size(data{1},1),size(data{1},2),size(data{1},3),length(data));
for i = 1:length(data)
    dataTimeAvg(:,:,:,i) = mean(data{i},4);
end

%...Compute information on molecular mass, gas constant, specific heat
%   ratio and speed of sound
molarMassTimeAvg = arrayfun(@(i)dataTimeAvg(:,:,:,min(gasesLoc)+i-1) * masses(i),1:G,'UniformOutput',false);
molarMassTimeAvg = sum(cat(4,molarMassTimeAvg{:}),4);
gasConstantTimeAvg = Ru./molarMassTimeAvg;
specificHeatRatioTimeAvg = 1./(1-gasConstantTimeAvg./dataTimeAvg(:,:,:,cpLoc));
specificHeatRatioTimeAvg(dataTimeAvg(:,:,:,HLoc) > 0.075) = 5/3; % correct values, since MCD gives weird results
speedOfSoundTimeAvg = sqrt(specificHeatRatioTimeAvg .* gasConstantTimeAvg .* ...
    dataTimeAvg(:,:,:,tempLoc));

%...Add new data to atmosphere table
dataTimeAvg = cat(4,dataTimeAvg,gasConstantTimeAvg,...
    specificHeatRatioTimeAvg,molarMassTimeAvg);

%...Compute other properties
numberDensityTimeAvg = numberDensityFun(molarMassTimeAvg,dataTimeAvg(:,:,:,densLoc));
collisionDiameterTimeAvg = arrayfun(@(i)dataTimeAvg(:,:,:,min(gasesLoc)+i-1) * collision(i),1:G,'UniformOutput',false);
collisionDiameterTimeAvg = sum(cat(4,collisionDiameterTimeAvg{:}),4);
knudsenNumberTimeAvg = knudsenNumberFun(molarMassTimeAvg,dataTimeAvg(:,:,:,densLoc),...
    collisionDiameterTimeAvg);

%...Verify results
if showFigures && ~saveFigures
    [lonH,latH,altH] = meshgrid(longitude,latitude,hs);
    
    figure;
    S = slice(lonH,latH,altH,molarMassTimeAvg,[],[],hs(1:2:end));
    plotSettings(S,'',logspace(log10(hs(1)),log10(hs(end)),7),NaN,true,true)
    
    figure;
    S = slice(lonH,latH,altH,gasConstantTimeAvg,[],[],hs(1:2:end));
    plotSettings(S,'',logspace(log10(hs(1)),log10(hs(end)),7),NaN,true,true)
    
    figure;
    S = slice(lonH,latH,altH,specificHeatRatioTimeAvg,[],[],hs(1:2:end));
    plotSettings(S,'',logspace(log10(hs(1)),log10(hs(end)),7),NaN,true,true)
end

%% Atmospheric Composition

%...Compute number density per altitude
if fullDatabase
    %...Save elements for SPARTA analysis
    density = dataLatLonTimeAvg(:,densLoc);
    pressure = dataLatLonTimeAvg(:,presLoc);
    temperature = dataLatLonTimeAvg(:,tempLoc);
    gasRatio = dataLatLonTimeAvg(:,gasesLoc);
    gasConstant = gasConstantLatLonTimeAvg;
    specificHeatRatio = specificHeatRatioLatLonTimeAvg;
    speedOfSound = speedOfSoundLatLonTimeAvg;
    numberDensity = numberDensityLatLonTimeAvg;
    knudsenNumber = knudsenNumberLatLonTimeAvg;
    save('../SPARTA/MCDEnvir','altitude','density','pressure','temperature',...
        'gasRatio','gasStr','gasConstant','specificHeatRatio','speedOfSound','numberDensity','knudsenNumber');
    
    if showFigures
        %...Plot gas percentage
        F = figure('rend','painters','pos',figSizeMedium);
        styles = repmat({'-',':','--','-.'},[1,3]);
        hold on
        for g = 1:G
            plot(hs,gasRatio(:,g),'LineWidth',1.75,'LineStyle',styles{g})
        end
        hold off
        xlabel('Altitude [km]')
        ylabel('Gas Fraction [-]')
        legend(gasStr,'Location','Best')
        set(gca,'FontSize',17.5,'XScale','log','YScale','log')
        grid on
        xlim([hs_plot(1),hs_plot(end)]), xticks([50,100,250,500,1500])
        ylim([1e-10,1])
        if ~saveFigures, title('Atmospheric Composition'),
        else, saveas(F,'../../Report/figures/mars_atm_comp','epsc'), end
        
        %...Plot molar mass and collision diameter
        F = figure('rend','painters','pos',figSizeMedium);
        yyaxis left
        hold on
        plot(hs,molarMassLatLonTimeAvg*1e3,'LineWidth',1.75)
        plot([10,10^6],masses(3)*1e3*ones(1,2),'LineWidth',1.25,'LineStyle','--')
        plot([10,10^6],masses(4)*1e3*ones(1,2),'LineWidth',1.25,'LineStyle','-.')
        hold off
        ylabel('Molar Mass [g mol^{-1}]')
        set(gca,'FontSize',17.5,'XScale','log')
        yyaxis right
        plot(hs,abs(sum(gasRatio,2)-1),'LineWidth',1.75)
        ylabel('Offset in Gas Presence [-]')
        set(gca,'FontSize',15,'XScale','log','YScale','log')
        xlabel('Altitude [km]')
        legend('Mars','CO_2','H','Offset','Location','Best')
        grid on
        xlim([hs_plot(1),hs_plot(end)]), xticks([50,100,250,500,1500])
        if ~saveFigures, title('Atmospheric Composition'),
        else, saveas(F,'../../Report/figures/mars_atm_mol','epsc'), end
        
        %...Plot number density and Knudsen
        F = figure('rend','painters','pos',figSizeMedium);
        yyaxis left
        loglog(hs,numberDensity,'LineWidth',1.75)
        ylabel('Number Density [-]')
        yyaxis right
        loglog(hs,knudsenNumber,'LineWidth',1.75)
        ylabel('Knudsen Number [-]')
        xlabel('Altitude [km]')
        set(gca,'FontSize',17.5)
        grid on
        xlim([hs_plot(1),hs_plot(end)]), xticks([50,100,250,500,1500])
        if ~saveFigures, title('Number Density and Knudsen'),
        else, saveas(F,'../../Report/figures/mars_atm_num','epsc'), end
    end
    
    clear density pressure temperature gasRatio gasConstant specificHeatRatio speedOfSound numberDensity knudsenNumber
end

%% Create Tabulated Atmosphere

if genTable
    %...Average over latitude, longitude, time, to get tabulated atmosphere
    %   as a function of altitude only
    tabularLatLonTimeAvg = TabulatedMCDAtmosphere('f(h)',dataLatLonTimeAvg);
    
    %...Average over time, to get tabulated atmosphere as a function of
    %   altitude, latitude and longitude
    tabularTimeAvg = TabulatedMCDAtmosphere('f(h,d,l)',dataTimeAvg);
end

%% Plot

if showFigures
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
            if ~saveFigures, subplotTitle(dataStr{i}),
            else, saveas(F,['../../Report/figures/mars_',dataLabel{i}],'epsc'), end
        else % other data
            F = figure('rend','painters','pos',figSizeLarge);
            for t = 1:T
                subplot(2,2,t)
                S = slice(lonH,latH,alt,data{i}(:,:,:,t),[],[],hs_plot);
                plotSettings(S,dataUnit{i},hs_plot,ts(t),true,true)
            end
            if ~saveFigures, subplotTitle(dataStr{i}),
            else, saveas(F,['../../Report/figures/mars_',dataLabel{i}],'epsc'), end
        end
    end
    
    %...Gas constant
    if fullDatabase
%         R = Ru./molarMass;
        gasConstant = 197*ones(1,1,1,T);
    else
        gasConstant = 197*ones(1,1,1,T);
    end
    specificHeatRatio = 1./(1-gasConstant./data{cpLoc});

    %...Speed of sound
    F = figure('rend','painters','pos',figSizeLarge);
    for t = 1:T
        subplot(2,2,t)
        sos = sqrt(gasConstant(:,:,:,t).*specificHeatRatio(:,:,:,t).*data{tempLoc}(:,:,:,t));
        S = slice(lonH,latH,alt,sos,[],[],hs_plot);
        plotSettings(S,dataUnit{windLoc},hs_plot,ts(t),true,true)
    end
    if ~saveFigures, subplotTitle('Speed of Sound'),
    else, saveas(F,'../../Report/figures/mars_sos','epsc'), end
%%
    %...Mean values
    if genTable
        mean_plot_params = [densLoc,presLoc,tempLoc,size(tabularTimeAvg,4)-2,size(tabularTimeAvg,4)-1];
        mean_units = {'kg m^{-3}','Pa','K','J kg^{-1} K^{-1}',''};
        mean_titles = {'Density','Pressure','Temperature','Gas Constant','Specific Heat Ratio'};
        mean_labels = {'mars_dens_mean','mars_pres_mean','mars_temp_mean','mars_gas_const_mean','mars_heat_ratio_mean'};
        for i = 1:length(mean_plot_params)
            F = figure('rend','painters','pos',figSizeLarge);
            for h = 1:length(hs_plot)
                subplot(2,3,h)
                S = pcolor(lonH(:,:,1),latH(:,:,1),tabularTimeAvg(:,:,hs_loc_plot(h),mean_plot_params(i)));
                plotSettings(S,mean_units{i},[],hs_plot(h),false,false)
            end
            if ~saveFigures, subplotTitle(['Average ',mean_titles{i}]),
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

%% Close All Figures

close all;

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
    global saveFigures
    c = colorbar;
    c.Label.String = label;
    xlabel('Long. [deg]'), xticks(-180:90:180)
    ylabel('Lat. [deg]'), yticks(-90:45:90)
    if three
        view(21,21)
        set(S,'EdgeColor','none')
        if saveFigures
            set(S,'EdgeColor','none')
        else
            set(S,'EdgeColor','none','FaceColor','interp','FaceAlpha','interp')
            alpha('color')
            alphamap('rampdown')
            alphamap('increase',0.25)
        end
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