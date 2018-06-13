fclose('all'); clear all; close all force; profile off; clc; format long g; rng default;
addpath functions tests

%%

data = [-3,-0.9999779095030014
    -2.8125,-0.9999303492456073
    -2.625,-0.9997946242638588
    -2.4375,-0.9994334567454198
    -2.25,-0.9985372834133188
    -2.0625,-0.9964637508747902
    -1.875,-0.9919900576701199
    -1.6875,-0.982989716601978
    -1.5,-0.9661051464753108
    -1.3125,-0.9365685747113888
    -1.125,-0.8883882317017078
    -0.9375,-0.8151024010343998
    -0.75,-0.7111556336535152
    -0.5625,-0.5736744566155919
    -0.375,-0.4041169094348223
    -0.1875,-0.2091176770593758
    0,0
    0.1875,0.2091176770593758
    0.375,0.4041169094348223
    0.5625,0.5736744566155919
    0.75,0.7111556336535152
    0.9375,0.8151024010343998
    1.125,0.8883882317017078
    1.3125,0.9365685747113888
    1.5,0.9661051464753108
    1.6875,0.982989716601978
    1.875,0.9919900576701199
    2.0625,0.9964637508747902
    2.25,0.9985372834133188
    2.4375,0.9994334567454198
    2.625,0.9997946242638588
    2.8125,0.9999303492456073
    3,0.9999779095030014];

interp1( data(:,1), data(:,2), 4, 'spline')

%%

clc;

independentValues = cell(1,2);
for i = 0:4
    independentValues{1}(i + 1) = 1950.0 + i * 10.0;
end
for i = 0:2
    independentValues{2}(i + 1) = 10.0 + i * 10.0;
end

dependentData = zeros(5,3);
dependentData(1,1) = 150.697;
dependentData(1,2) = 199.592;
dependentData(1,3) = 187.625;
dependentData(2,1) = 179.323;
dependentData(2,2) = 195.072;
dependentData(2,3) = 250.287;
dependentData(3,1) = 203.212;
dependentData(3,2) = 179.092;
dependentData(3,3) = 322.767;
dependentData(4,1) = 226.505;
dependentData(4,2) = 153.706;
dependentData(4,3) = 426.730;
dependentData(5,1) = 249.633;
dependentData(5,2) = 120.281;
dependentData(5,3) = 598.243;

interpValues = { [1940,15], [2000,15], [1975,7.5], [1975,32.5], [1940,32.5] };
interpValues = { [1950,15], [1990,15], [1975,10], [1975,30], [1950,30.0] };

nearestLowerIndices = { [1,1], [4,1], [3,1], [3,2], [1,2] };

for i = 1:length(interpValues)
%     performRecursiveInterpolationStep( independentValues, dependentData, 1, interpValues{ i }, ...
%         [ 0, 0 ], nearestLowerIndices{ i }, 2 )
    interp2(independentValues{2},independentValues{1},dependentData,interpValues{i}(2),interpValues{i}(1),'linear',NaN)
end

% interp2(independentValues{2},independentValues{1},dependentData,15,1975,'linear')

%%

performRecursiveInterpolationStep( 1, interpValues{ 1 }, [ 1, 1 ], [ 1, 1 ] )
