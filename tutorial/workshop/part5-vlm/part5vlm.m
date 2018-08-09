%% generate required paths for program execution
addpath(genpath('..\..\..\'));


%% load XML File
myAC=class_aircraft('part5-B787-800.xml',1);

%% grid settings

% set aerodynamic panel size
myAC.grid_settings.x_max_grid_size=myAC.reference.c_ref/5;
myAC.grid_settings.y_max_grid_size=myAC.reference.b_ref/20;
myAC.grid_settings.aerodynamic_fuselage=0;

% set  structural grid size
myAC.grid_settings.dy_max_struct_grid=myAC.reference.b_ref/30;

%% grid computation
% myAC=myAC.compute_shell_coords();
myAC=myAC.compute_grid();


%% Setup aerodynamic state  V   AoA° SideSlip°   Ma   dens
myState=class_aero_state(  300, 1,    0,         0,  1.225);
%% Setup VLM Solver

myVLM=class_VLM_solver(myAC.grid, myAC.te_idx,myAC.panels,myState,myAC.reference, myAC.control_surfaces_parents);

myVLM=myVLM.f_solve_std;
myVLM=myVLM.plot_cp;

%% solve derivatives wrt to aoa/beta p q r
myVLM=myVLM.f_solve_full;

%% deflect control surface
myAC=myAC.f_set_control_surface('outerFlap',10);

myAC=myAC.compute_grid();

myVLM2=class_VLM_solver(myAC.grid, myAC.te_idx,myAC.panels,myState,myAC.reference, myAC.control_surfaces_parents);
myVLM2=myVLM2.f_solve_std;
figure;
myVLM2=myVLM2.plot_cp;