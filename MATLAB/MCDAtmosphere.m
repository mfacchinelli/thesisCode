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

genTable = false;
showFigures = false;
saveFigures = false;
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
[lonH,latH,altH] = meshgrid(longitude,latitude,hs);
if showFigures && ~saveFigures
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

%% Perturbed Atmospheric Density

clc

%...Get random coefficients
generateRandomValues = false;
x = [1.71415  0.0272539   0.110777  -0.146923  -0.055571 -0.0475874  0.0165448    1.50714   0.478782  -0.297889];
x = [1.71415  0.0272539   0.110777  -0.146923  -0.055571 -0.0475874  0.0165448    1.50714   0.478782  1];
% x = [1.76748 -0.163649  0.361739 -0.255198  -0.19689 -0.143325 0.0385251   1.48386  0.289886 -0.170874];
% x = [1.7885    0.072896    0.049453     0.39692    -0.20112     0.17416     0.20877    0.87814     1.1078    0.41708];
% x = [1.0679    0.18702    -0.048105      0.22215    -0.19121    -0.35057    -0.35559    1.2441    0.91131    0.90197];
if generateRandomValues
    a = 0.75 + 1.25 * rand(1);
    b = 0.25 * randn(1,3);
    c = 0.25 * randn(1,3);
    d = ( 1 + 0.5 * randn(1,3) );
else
    a = x(1);
    b = x([2,4,6]);
    c = x([3,5,7]);
    d = x(8:10);
end

%...Get independent variables
lonP = ( lonH - longitude(1) ) / ( longitude(end) - longitude(1) ) * 360 * d(1);
latP = ( latH - latitude(1) ) / ( latitude(end) - latitude(1) ) * 360 * d(2);
altP = ( altH - altitude(1) ) / ( 1500 - altitude(1) ) * 360 * d(3);

%...Compute perturbing term
perturbations = abs( a + ...
    b(1) * cosd(lonP) + c(1) * sind(lonP) + ...
    b(2) * cosd(latP) + c(2) * sind(latP) + ...
    b(3) * cosd(altP) + c(3) * sind(altP) );
perturbations( perturbations<0.5 ) = 0.5;

format short g, table(a,b,c,d), format long g

%...Plot results
F = figure('rend','painters','pos',figSizeLarge);
for h = 1:length(hs_loc_plot)
    subplot(2,3,h)
%     S = pcolor(lonH(:,:,1),latH(:,:,1),perturbations(:,:,hs_loc_plot(h)));
    S = pcolor(lonH(:,:,1),latH(:,:,1),perturbations(:,:,hs_loc_plot(h)) .* dataTimeAvg(:,:,hs_loc_plot(h),densLoc));
%     S = pcolor(lonH(:,:,1),latH(:,:,1),dataTimeAvg(:,:,hs_loc_plot(h),densLoc));
    plotSettings(S,'kg m^{-3}',[],hs_plot(h),false,false)
end
if ~saveFigures, subplotTitle(['Perturbed Average Density']),
else, saveas(F,'../../Report/figures/mars_dens_mean_pert','epsc'), end

figure
S = slice(lonH,latH,altH,perturbations,[],[],linspace(50,1500,10));
colormap jet
set(S,'EdgeColor','none','FaceColor','interp','FaceAlpha','interp')
alpha('color')
alphamap('rampdown')
alphamap('increase',0.25)
colorbar

%% Perturbed Atmosphere Evolution

%...Perturbation coefficients
X = [1.71415  0.0272539   0.110777  -0.146923  -0.055571 -0.0475874  0.0165448    1.50714   0.478782  -0.297889
1.68759  -0.014518  0.0674644   -0.14893 -0.0207594 -0.0471949   0.040052    1.45393    0.42858  -0.380483
1.61521 0.00597483  0.0607304  -0.152705 -0.0415833 -0.0742231  0.0450976    1.45956    0.45225  -0.416446
1.52041   0.0497252   0.0587099    -0.18566  -0.0411544  -0.0766623 -0.00665937      1.3907    0.393526   -0.367184
1.66091   0.0257186     0.11099   -0.186187 -0.00363248  -0.0912492  -0.0180747      1.4518    0.501847   -0.264378
1.61002 0.00285677  0.0970087  -0.172205 -0.0170566 -0.0588736 -0.0248429    1.43076   0.513475  -0.334797
1.64714 -0.0238725  0.0917116  -0.156685   -0.02832  0.0117141 -0.0115811    1.52585   0.533104  -0.292143
1.53346 -0.00145534    0.100958    -0.19212 -0.00974424  -0.0195981  -0.0167557     1.49716    0.515075   -0.408583
1.68134   0.0361143   0.0980484   -0.217666 -0.00382061  -0.0316979  0.00971486     1.53213    0.612151   -0.379812
1.59917  0.0332527   0.144773  -0.212507 -0.0676948 -0.0449416 0.00353183    1.50844   0.568709  -0.411243
1.71336  0.0188096   0.154981  -0.211699 -0.0354057 -0.0160633 -0.0122023    1.61909   0.697184  -0.315977
1.74304  0.0144644   0.141447   -0.19586 -0.0501727  0.0317361 -0.0126983    1.67206    0.77497   -0.30437
1.9766  0.0476793   0.152053  -0.172183 -0.0621753   0.046012  0.0219906    1.86091   0.903295  -0.181904
1.92003  0.0347641   0.145976  -0.171889 -0.0778345  0.0187404   0.030562    1.81997    0.77384  -0.252122
2.00841  0.0357696   0.164068  -0.154676 -0.0932909  0.0259087  0.0365928    1.79479   0.855815  -0.293374
2.08142 -0.00584229    0.188115   -0.153282  -0.0582121   0.0134735   0.0419491     1.83155    0.947064   -0.255647
2.17725   0.026738   0.205682   -0.18933 -0.0587854 0.00305996  0.0470307    1.89843   0.978075  -0.235294
2.19684  0.0242813   0.227447  -0.175553 -0.0862279  0.0151203  0.0786095    1.97455   0.985824  -0.143747
2.13507   0.0586551    0.205047   -0.201831  -0.0987977 -0.00923905   0.0785915     1.98567    0.977641   -0.100489
2.21379  0.0665699   0.204487   -0.17427 -0.0805926   0.010397   0.100427    2.03819    1.03607 -0.0401719
2.29134   0.103343   0.218426  -0.168765 -0.0905459 0.00618746   0.102828     2.1664    1.07855  0.0157786
2.47543   0.147864   0.211817  -0.114546 -0.0966025 0.00566185  0.0982684    2.27863    1.22824   0.144652
2.38193   0.135665   0.193998  -0.130253 -0.0908412  -0.024907  0.0811168    2.27414    1.12898   0.095198
2.20783   0.0938771    0.189059   -0.121936  -0.0734825  -0.0177815   0.0867224     2.18908     1.02261 -0.00101503
2.1563  0.0978768   0.164299 -0.0859596 -0.0922169 -0.0632552  0.0615566    2.10069   0.917216 -0.0112893
2.11552  0.0914607   0.173666 -0.0890087 -0.0918561 -0.0531196  0.0746797    2.06182   0.865214 -0.0958879
2.06344    0.10209   0.195076  -0.107493 -0.0838192 -0.0899616  0.0590123     2.0538   0.886891  -0.177832
2.26352  0.0907691    0.24423 -0.0913177 -0.0940806 -0.0690579  0.0451566    2.26333    1.07329  0.0178345
2.17106  0.0512555   0.230645  -0.101822 -0.0763309 -0.0405329  0.0412203     2.2001    1.01903 0.00723909
1.90508  0.0438242   0.224711  -0.112876  -0.103723 -0.0505243   0.037575    2.05623   0.948287  -0.103036
1.98712  0.0618706   0.234201  -0.145719 -0.0982657 -0.0578929  0.0418398     2.1059     1.0334  0.0645327
2.20999   0.101317   0.257053  -0.124575 -0.0802723 -0.0772104  0.0362539    2.14724     1.2005   0.253625
2.34383    0.11746   0.257479 -0.0936536 -0.0724077  -0.092876  0.0543474    2.22148    1.18773   0.261273
2.13707  0.105415  0.244459 -0.114416 -0.058077  -0.10522 0.0249768   2.02143   1.00383 0.0722512
2.0877    0.07303   0.281302 -0.0883171 -0.0519792  -0.146048 0.00667274    1.96602   0.953958   0.022823
2.04336  0.0679405   0.296755 -0.0804674 -0.0536847  -0.162097   0.012732    1.99754   0.965724  -0.116434
1.97923  0.0491896   0.268076  -0.101313 -0.0815783  -0.184594 0.00660692     1.9281    0.93528  -0.140771
2.0658   0.0452151    0.297388  -0.0976338  -0.0678827   -0.192627 -0.00443283     1.96792    0.955606  -0.0670719
2.17128   0.037362   0.301291  -0.101293 -0.0637327  -0.201586  0.0112422    2.11321   0.942752 -0.0248813
2.23306  0.0535108   0.302712  -0.120516 -0.0809033  -0.171896  0.0459813    2.16759    1.01234   0.112766
2.05681   0.0298561    0.315374   -0.159197  -0.0827201    -0.20122   0.0697673     2.04779    0.923758 -0.00212248
2.10258 0.0654126  0.331537 -0.174055 -0.089867 -0.198014 0.0699462   2.01952   1.02694 0.0986318
1.8965  0.0306899   0.318974  -0.191063  -0.136679  -0.198694  0.0694592    1.91051    0.89119 -0.0133946
1.74709 0.0144663  0.308267  -0.17372 -0.150479 -0.223494 0.0644668   1.84568  0.722312 -0.145387
1.89483 0.0311785  0.333935 -0.144264 -0.122728 -0.187848 0.0846128   2.01564  0.916797 0.0452709
1.74356 0.0217795  0.347691 -0.168462 -0.143638 -0.205423 0.0607542   1.83069  0.668548 -0.160549
1.96634  0.0445593    0.33836  -0.144885  -0.122402  -0.199514  0.0738161    1.99048   0.804368 -0.0957012
2.02747 0.0515756  0.342879 -0.129931 -0.102207 -0.216108  0.052036   2.04116  0.914419 -0.156028
2.09606  0.0576392   0.354626  -0.124346 -0.0906282  -0.219677  0.0494646    2.12319    1.03233  -0.129705
2.07712  0.0752892   0.354639  -0.141899 -0.0597999   -0.24415  0.0490001    2.09533   0.951896  -0.116022
2.2504 0.0758871  0.383434 -0.133553 -0.060752 -0.223853  0.102517   2.22214   1.12151 0.0287975];

%...Loop over each perturbation coefficient set
F = figure('rend','painters','pos',figSizeLarge);
for i = 1:6
    %...Get coefficients
    x = X(9*(i-1)+i,:);
    a = x(1);
    b = x([2,4,6]);
    c = x([3,5,7]);
    d = x(8:10);

    %...Get independent variables
    lonP = ( lonH - longitude(1) ) / ( longitude(end) - longitude(1) ) * 360 * d(1);
    latP = ( latH - latitude(1) ) / ( latitude(end) - latitude(1) ) * 360 * d(2);
    altP = ( altH - altitude(1) ) / ( 1500 - altitude(1) ) * 360 * d(3);
    
    %...Compute perturbing term
    perturbations = abs( a + ...
        b(1) * cosd(lonP) + c(1) * sind(lonP) + ...
        b(2) * cosd(latP) + c(2) * sind(latP) + ...
        b(3) * cosd(altP) + c(3) * sind(altP) );
    perturbations( perturbations<0.5 ) = 0.5;
    
    subplot(2,3,i)
    S = slice(lonH,latH,altH,perturbations,[],[],linspace(50,1500,5));
    plotSettings(S,'Atmosphere Scaling [-]',[50,500,1000,1500],2*(i-1)*100+1,false,true)
    colormap jet
    if saveFigures
        set(S,'EdgeColor','none')
    else
        set(S,'EdgeColor','none','FaceColor','interp','FaceAlpha','interp')
        alpha('color')
        alphamap('rampdown')
        alphamap('increase',0.25)
    end
    xlabel('Long. [deg]'), xticks(-180:90:180)
    ylabel('Lat. [deg]'), yticks(-90:45:90)
    zlabel('Alt. [km]'), zticks(linspace(50,1500,5))
    colorbar
    set(gca,'FontSize',15)
    title(num2str(30*(i-1)+1))
end
if saveFigures, saveas(F,'../../Report/figures/mars_dens_pert_coeff_evol','epsc'), end

%% Make GIF

%...Loop over each perturbation coefficient set
for i = 1:size(X,1)
    %...Get coefficients
    x = X(i,:);
    a = x(1);
    b = x([2,4,6]);
    c = x([3,5,7]);
    d = x(8:10);

    %...Get independent variables
    lonP = ( lonH - longitude(1) ) / ( longitude(end) - longitude(1) ) * 360 * d(1);
    latP = ( latH - latitude(1) ) / ( latitude(end) - latitude(1) ) * 360 * d(2);
    altP = ( altH - altitude(1) ) / ( 1500 - altitude(1) ) * 360 * d(3);
    
    %...Compute perturbing term
    perturbations = abs( a + ...
        b(1) * cosd(lonP) + c(1) * sind(lonP) + ...
        b(2) * cosd(latP) + c(2) * sind(latP) + ...
        b(3) * cosd(altP) + c(3) * sind(altP) );
    perturbations( perturbations<0.5 ) = 0.5;
    
    F = figure('rend','painters','pos',figSizeMedium);
    S = slice(lonH,latH,altH,perturbations,[],[],linspace(50,1500,10));
    colormap jet
    set(S,'EdgeColor','none','FaceColor','interp','FaceAlpha','interp')
    alpha('color')
    alphamap('rampdown')
    alphamap('increase',0.25)
    xlabel('Long. [deg]'), xticks(-180:90:180)
    ylabel('Lat. [deg]'), yticks(-90:45:90)
    zlabel('Alt. [km]'), zticks([50,500,1000,1500])
    colorbar
    title(num2str(3*(i-1)+1))
    set(gca,'FontSize',15)
    if saveFigures, saveas(F,['../../Report/figures/pert_coeff_evol/',num2str(i)],'png'), end
    close(gcf)
end

%..Create GIF command string
if saveFigures
    repository = ['"/Users/Michele/Library/Mobile Documents/com~apple~CloudDocs/University/'...
        'Master Thesis/Report/figures/pert_coeff_evol"'];
    GIFCommand = ['cd ',repository,'; convert '];
    for i = 1:size(X,1)
        is = num2str(i);
        GIFCommand = [GIFCommand,sprintf('%s.png ',is)];
    end
    GIFCommand = [GIFCommand,['movie.gif; rm *.png']];
    
    %...Add command string to clipboard
    clipboard('copy',GIFCommand)
    helpdlg(['Figures have been generated.',newline,...
        'GIF command string is ready to be pasted in the Terminal!'],'Ready to GIFy!');
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

if saveFigures, close all, end

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