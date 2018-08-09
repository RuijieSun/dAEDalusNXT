%% generate required paths for program execution
addpath(genpath('..\..\..\'));


%% load XML File
myAircraft=class_aircraft('simpleWingWithStruct.xml',1);

%% grid settings

% set aerodynamic panel size
myAircraft.grid_settings.x_max_grid_size=myAircraft.reference.c_ref/8;
myAircraft.grid_settings.y_max_grid_size=myAircraft.reference.b_ref/30;

% set  structural grid size
myAircraft.grid_settings.dy_max_struct_grid=myAircraft.reference.b_ref/30;

myAircraft=myAircraft.compute_grid(); %<- aircraft grid is always required for structural model creation

%% create Structural Model
%options
structOpts=class_wingstructure_solver_settings;
%create
[myAircraft,myStructure]=create_structural_model(myAircraft, structOpts);

%% visual checks
figure; 
myAircraft.plot;
myStructure.plot_geometry;

%% visual checks2
figure; 
myAircraft.plot_grid;
myStructure.plot_structure;

%% clamp wing
myStructureClamped=myStructure;
myStructureClamped.beam(1)=myStructureClamped.beam(1).f_add_boundary_condition(class_boundary_condition(myStructureClamped.beam(1).nel/2+1,[1 1 1 1 1 1],[0 0 0 0 0 0]));
%plot deformations
myStructureClamped=myStructureClamped.f_solve;
figure; myStructureClamped.plot_deformations; 
%set force, solve, plot
myStructureClamped.Ftest(end-3)=100000;
myStructureClamped=myStructureClamped.f_solve;
myStructureClamped.plot_deformations;

%% solve modes
myStructureClamped=myStructureClamped.f_solve_modes;
%visual checks in paraview (arg1: which modes; arg2: scale factor; 
%arg3: folder; arg4: frames per animation
myStructureClamped.write_mode_animation_in_paraview(1:5,50,'modeshapes',30)

myStructureClamped.modefrequencies(1:20)
%% sovle free modes
myStructure=myStructure.f_solve_free_modes;
myStructure.modefrequencies(1:20)
myStructure.write_mode_animation_in_paraview(7:10,50,'modeshapes_free',30)

