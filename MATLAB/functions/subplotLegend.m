%  MATLAB Function < subplotTitle >
%
%  Author:          Michele Facchinelli
%  Date:            22nd May 2017
%  Purpose:         give overall legend to subplots; reference: goo.gl/Ga7AfW
%  Input:
%   - data:         cell array with legend names
%   - specifier:    specify location of legend
%  Output:
%   - N/A

function subplotLegend(data,specifier)

%...Default value to specifier
if nargin == 1
    specifier = 'mid';
end

%...Add legend to plot
legends = legend(data);
switch specifier
    case 'top'
        set(legends,'Position',[0.5,0.95,0.01,0.01],'Orientation','horizontal','Units','normalized');
    case 'mid'
        set(legends,'Position',[0.5,0.5,0.01,0.01],'Orientation','horizontal','Units','normalized');
    case 'low'
        set(legends,'Position',[0.5,0.05,0.01,0.01],'Orientation','horizontal','Units','normalized');
end