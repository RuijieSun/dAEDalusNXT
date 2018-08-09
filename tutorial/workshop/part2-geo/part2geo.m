%% generate required paths for program execution
addpath(genpath('..\..\..\'));


%% load XML File
myAircraft=class_aircraft('part2-B787-800.xml',1);

%% plot aircraft for visual checking
figure; 
myAircraft.plot;

%% grid computation
myAircraft=myAircraft.compute_grid(); %<- explain OOP

%% grid plotting (2D grid)
figure; 
myAircraft.plot_grid

%% grid plotting (Volume Grid)
figure;
myAircraft.plot_grid_vol();

%% including the fuselage:
myAircraft.grid_settings.aerodynamic_fuselage=1; %<- for 2d grid
myAircraft=myAircraft.compute_fuselage_grid();
figure;
myAircraft.plot_grid_vol();

%% For better pictures, the geometry may also be written to a tecplot file
myAircraft.write_tecplot_volgrid([myAircraft.name '_volGrid'])

%% modifying the grid size
myAircraft.grid_settings.x_max_grid_size=1;
myAircraft.grid_settings.y_max_grid_size=1;
myAircraft=myAircraft.compute_grid(); 
figure; 
myAircraft.plot_grid