function [tabular,tabular_interp] = TabulatedMCDAtmosphere(data,mode)
%% Main File

%...Globals
global plotResults Ru latitude longitude altitude densLoc presLoc tempLoc cpLoc
tab_atm_param = [densLoc,presLoc,tempLoc];

%...Switch mode based on selected averaging process
switch mode
    case 'f(h)'
        %...Average over latitude, longitude, time, to get tabulated atmosphere
        %   over altitude only
        j = 0;
        tabular = zeros(size(data{1},3),length(tab_atm_param));
        for i = tab_atm_param
            j = j+1;
            tabular(:,j) = mean(mean(mean(data{i},1),2),4);
        end
        
        %...Add gas constant and specific heat ratio
        tabular(:,4) = Ru./squeeze(mean(mean(mean(data{end},1),2),4));
        tabular(:,5) = 1./(1-tabular(:,4)./squeeze(mean(mean(mean(data{cpLoc},1),2),4)));
        
        %...Interpolate data
        hs_interp = logspace(log10(altitude(1)),log10(altitude(end)),1e3);
        tabular_interp = interp1(altitude,tabular,hs_interp,'linear','extrap');
        tabular_interp = horzcat(hs_interp'*1e3,tabular_interp);

        %...Save to file
        path = '/Users/Michele/GitHub/tudat/tudatBundle/tudat/Tudat/External/AtmosphereTables/';
        fileID = fopen([path,'MCDMeanAtmosphere.dat'],'w');
        fprintf(fileID,[repmat('%.10e\t',[1,size(tabular_interp,2)]),'\n'],tabular_interp');
        fclose(fileID);
        
        %..Plot mean atmosphere
        if plotResults
            figure;
            for i = 1:3
                subplot(1,3,i)
                hold on
                scatter(tabular_interp(:,i+1),tabular_interp(:,1)/1e3)
                scatter(mean(mean(mean(data{tab_atm_param(i)},1),2),4),altitude)
                if i == 1 || i == 2
                    if i == tab_atm_param(tab_atm_param==densLoc)
                        scatter(0.02*exp(-hs_interp/11.1),hs_interp), xlim([1e-30,1e0])
                    end
                    set(gca,'XScale','log')
                end
                hold off
                grid on
                set(gca,'YScale','log')
            end
        end
    case 'f(h,d,l)'
        %...Average over time, to get tabulated atmosphere over altitude, 
        %   latitude and longitude
        j = 0;
        tabular = zeros(size(data{1},1),size(data{1},2),size(data{1},3),length(tab_atm_param));
        for i = tab_atm_param
            j = j+1;
            tabular(:,:,:,j) = mean(data{i},4); % this step also changes the order of columns (see tab_atm_param)
        end
        
        %...Add gas constant and specific heat ratio
        tabular(:,:,:,4) = Ru./mean(data{end},4);
        tabular(:,:,:,5) = 1./(1-tabular(:,:,:,4)./mean(data{cpLoc},4));

        %...Interpolate data
        hs_interp = logspace(log10(altitude(1)),log10(altitude(end)),1e3);
        lon_interp = linspace(longitude(1),longitude(end),150);
        lat_interp = linspace(latitude(1),latitude(end),75);
        for i = 1:size(tabular,4)
            [lon_m,lat_m,alt_m] = meshgrid(lon_interp,lat_interp,hs_interp);
            tabular_interp(:,:,:,i) = interp3(longitude,latitude,altitude,tabular(:,:,:,i),...
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