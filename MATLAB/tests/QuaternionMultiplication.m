fclose('all'); clear all; close all force; profile off; clc; format long g; rng default;

%% Test 1

clc;

q1 = [ 1, 0, 0, 0 ]
q2 = angle2quat( 0, 10, 35, 'ZXY' )

quatmultiply( q1, q2 )

%% Test 2

clc;

q1 = angle2quat( 135, 10, -90, 'XYZ' )
q2 = angle2quat( 0, 10, 35, 'YZX' )

quatmultiply( q1, q2 )

%% Test 3

clc;

q1 = angle2quat( 135, -45, 215, 'YZY' )
q2 = angle2quat( 0, 15, 300, 'ZXZ' )

quatmultiply( q1, q2 )