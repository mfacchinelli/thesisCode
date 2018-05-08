%  MATLAB Function < subplotColorbar >
%
%  Author:          Michele Facchinelli
%  Date:            22nd May 2017
%  Purpose:         give overall legend to subplots; reference: goo.gl/Ga7AfW
%  Input:
%   - data:         cell array with legend names
%   - specifier:    specify location of legend
%  Output:
%   - N/A

function subplotColorbar(string)

%...Add colormap to plot
c = colorbar('northoutside');
c.Label.String = string;
c.Position = [0.842857142857143+0.05,0.10952380952381-0.05,0.0357142857142857,0.815476190476191+0.05];