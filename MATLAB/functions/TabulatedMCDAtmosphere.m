function tabular = TabulatedMCDAtmosphere(mode,tabular)
%% Main File

%...Globals
global showFigures latitude longitude altitude densLoc presLoc tempLoc
tab_atm_param = [densLoc,presLoc,tempLoc];

%...Switch mode based on selected averaging process
switch mode
    case 'f(h)'
        %...Average over latitude, longitude, time, to get tabulated atmosphere
        %   as a function of altitude only
        %   Only use information on interesting parameters
        tabular_mean = tabular(:,tab_atm_param);
        tabular_mean = horzcat(tabular_mean,tabular(:,end-2:end));
        
        %...Interpolate data
        hs_interp = logspace(log10(altitude(1)),log10(altitude(end)),1e3);
        tabular_interp = interp1(altitude,tabular_mean,hs_interp,'spline');
        tabular_interp = horzcat(hs_interp'*1e3,tabular_interp);

        %...Save to file
        path = '/Users/Michele/GitHub/tudat/tudatBundle/tudat/Tudat/External/AtmosphereTables/';
        fileID = fopen([path,'MCDMeanAtmosphere.dat'],'w');
        fprintf(fileID,[repmat('%.10e\t',[1,size(tabular_interp,2)]),'\n'],tabular_interp');
        fclose(fileID);
        
        %..Plot mean atmosphere
        referenceAltitude = 158;
        hs_interp(referenceAltitude)
        densityAtReferenceAltitude = tabular_interp(referenceAltitude,2)
        localScaleHeight = (tabular_interp(referenceAltitude,4)*tabular_interp(referenceAltitude,5)/3.711)
        if showFigures
            figure;
            for i = 1:3
                subplot(1,3,i)
                hold on
                plot(tabular_interp(:,i+1),tabular_interp(:,1)/1e3)
                scatter(tabular_mean(:,i),altitude)
                if i == 1 || i == 2
                    if i == tab_atm_param(tab_atm_param==densLoc)
                        scatter(0.02*exp(-hs_interp/11.1),hs_interp), xlim([1e-30,1e0])
                        scatter(densityAtReferenceAltitude*exp((hs_interp(referenceAltitude)-hs_interp)/...
                            (localScaleHeight/1e3)),hs_interp)
                        xlim([1e-30,1e0])
                        legend('Mean','Interpolated','h_0 = 0','h_0 = 110')
                    end
                    set(gca,'XScale','log')
                end
                hold off
                grid on
                set(gca,'YScale','log')
            end
        end
    case 'f(h,d,l)'
        %...Average over time, to get tabulated atmosphere as a function of
        %   altitude, latitude and longitude
        %   Only use information on interesting parameters
        tabular_mean = tabular(:,:,:,tab_atm_param);
        tabular_mean = cat(4,tabular_mean,tabular(:,:,:,end-2:end));

        %...Interpolate data
        hs_interp = logspace(log10(altitude(1)),log10(altitude(end)),500);
        lon_interp = linspace(longitude(1),longitude(end),150);
        lat_interp = linspace(latitude(1),latitude(end),75);
        for i = 1:size(tabular_mean,4)
            [lon_m,lat_m,alt_m] = meshgrid(lon_interp,lat_interp,hs_interp);
            tabular_interp(:,:,:,i) = interp3(longitude,latitude,altitude,tabular_mean(:,:,:,i),...
                lon_m,lat_m,alt_m,'spline');
        end
        
        %...Save to file
        path = ['/Users/Michele/GitHub/tudat/tudatBundle/tudat/Tudat/External/'...
            'AtmosphereTables/MCDMeanAtmosphereTimeAverage'];
        if ~exist(path,'dir')
            mkdir(path);
        end
        save('data/MCDMeanAtmosphere.mat','tabular_interp','lat_interp','lon_interp','hs_interp');
        
        %...Save files
        parameters = {'density','pressure','temperature','gasConstant','specificHeatRatio'};
        for i = 1:length(parameters)
            fileID = fopen(fullfile(path,sprintf('%s.dat',parameters{i})),'w');
            fprintf(fileID,'%d\n',3);
            fprintf(fileID,'\n\n');
            fprintf(fileID,[repmat('%.7e\t',[1,length(lon_interp)]),'\n'],deg2rad(lon_interp));
            fprintf(fileID,[repmat('%.7e\t',[1,length(lat_interp)]),'\n'],deg2rad(lat_interp));
            fprintf(fileID,[repmat('%.7e\t',[1,length(hs_interp)]),'\n'],hs_interp*1e3);
            for h = 1:length(hs_interp)
                fprintf(fileID,'\n\n');
                fprintf(fileID,[repmat('%.7e\t',[1,size(tabular_interp(:,:,h,i),1)]),'\n'],tabular_interp(:,:,h,i));
            end
            fclose(fileID);
        end
    otherwise
        error('Only mode supported are f(h) and f(h,d,l)')
end

end