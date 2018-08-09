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
%clamp
myStructure.beam(1)=myStructure.beam(1).f_add_boundary_condition(class_boundary_condition(myStructure.beam(1).nel/2+1,[1 1 1 1 1 1],[0 0 0 0 0 0]));
myStructure=myStructure.f_solve;
%set tip force
myStructure.Ftest(end-3)=100000;
% solve
myStructure=myStructure.f_solve;

%% create coupling matrices
myAircraft=myAircraft.compute_force_interpolation_matrix(myStructure);

%% visual checks
% where does panel x belong to
myAircraft.check_panel_to_beam_element(myStructure,10)
% which panels belong to beam element 10 of beam 1
myAircraft.check_panel_to_beam_element2(myStructure,1,10)


%% deform aerodynamic grid 
deflections=myStructure.f_get_deflections;
myAircraft=myAircraft.compute_deflected_grid(deflections);
%% plot deflected grid
figure;
myAircraft.plot_grid;
myAircraft.plot_grid_deflected;

