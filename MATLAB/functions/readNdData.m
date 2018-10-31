function [independentVariables,dependentVariables] = readNdData(input)

%...Open file
fileID = fopen(input,'r');
numberOfIndependentVariables = str2double(fgetl(fileID));
fgetl(fileID);

%...Read independent variables
independentVariables = cell(1,numberOfIndependentVariables);
for i = 1:numberOfIndependentVariables
    independentVariables{i} = str2num(fgetl(fileID));
end
fgetl(fileID);

%...Read dependent variables
if numberOfIndependentVariables == 2
    dependentVariables = zeros(length(independentVariables{1}),length(independentVariables{2}));
    for i = 1:length(independentVariables{1})
        dependentVariables(i,:) = str2num(fgetl(fileID));
    end
elseif numberOfIndependentVariables == 3
    dependentVariables = zeros(length(independentVariables{1}),length(independentVariables{2}),length(independentVariables{3}));
    for k = 1:length(independentVariables{end})
        for i = 1:length(independentVariables{1})
            dependentVariables(i,:,k) = str2num(fgetl(fileID));
        end
        fgetl(fileID);
    end
end
fclose(fileID);