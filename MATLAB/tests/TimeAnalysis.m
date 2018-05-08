fclose('all'); clear all; close all force; profile off; clc; format long g; rng default;
addpath tests functions

%% SPARTA and Modeling Variables

%...Select mode
mode = 'speedRatio'; % 'sphere' or 'speedRatio'

%...Figures setting
figSize = [440,378,840,630]; 

%...Repositories
SPARTARepository = ['/Users/Michele/Library/Mobile Documents/com~apple~CloudDocs/University/Master ',...
    'Thesis/Code/SPARTA/']; % path to SPARTA folder
SphereRepository = fullfile(SPARTARepository,'verification'); % path to sphere folder
OutputRepository = fullfile(SphereRepository,mode,'data'); % path to SPARTA output
ImageRepository = fullfile(SphereRepository,mode,'figures'); % path to SPARTA figures
repository = OutputRepository;

%% Retrieve Files

%...Scan image directory
gasses = dir(repository); gasses = gasses(4:end);

%...Continue to sub-directories
densities = cell(length(gasses),1);
speeds = cell(length(gasses),1);
angles = cell(length(gasses),1);
for g = 1:length(gasses)
    %...Scan gas directories
    den = dir(fullfile(repository,gasses(g).name));
    densities{g} = den(4:end);
    
    %...Continue to sub-directories
    for d = 1:length(densities{g})
        %...Scan density directories
        spe = dir(fullfile(repository,gasses(g).name,densities{g}(d).name));
        speeds{g,d} = spe(4:end);
        
        %...Continue to sub-directories
        for s = 1:length(speeds{g,d})
            ang = dir(fullfile(repository,gasses(g).name,densities{g}(d).name,speeds{g,d}(s).name));
            loc = {ang.name};
            loc = cellfun(@(x)x(1),cellfun(@(x)split(x,'.'),loc,'UniformOutput',false));
            angles{g,d,s} = ang(~cellfun(@isempty,loc));
        end
    end
end

%% Computation Times

%...Get all dates
times = [];
nums = cellfun(@(x)arrayfun(@(i)datenum(x(i).date),1:length(x))',angles,'UniformOutput',false);
for i = 1:numel(nums)
    times = vertcat(times,nums{i});
end
times = unique(times);
dates = datetime(times,'ConvertFrom','datenum');

%% Results

%...Plot
figure('rend','painters','pos',figSize);

subplot(1,2,1)
scatter(1:length(dates),dates,10)
title('Creation Dates')
grid on
set(gca,'FontSize',15)

subplot(1,2,2)
scatter(1:length(dates)-1,seconds(diff(dates)),10)
title('Creation Duration')
grid on
set(gca,'FontSize',15,'YScale','log')