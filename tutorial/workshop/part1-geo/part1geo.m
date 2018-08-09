%% generate required paths for program execution
addpath(genpath('..\..\..\'));


%% load XML File
myAircraft=class_aircraft('simpleWing.xml',1);

%% plot aircraft for visual checking
figure; 
plot(myAircraft);
% or 
% myAircraft.plot;