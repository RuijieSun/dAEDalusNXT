%% generate required paths for program execution
addpath(genpath('..\..\..\'));


%% load XML File
myAC=class_aircraft('part6-simple787.xml',1);

%% grid settings

% set aerodynamic panel size
myAC.grid_settings.x_max_grid_size=myAC.reference.c_ref/8;
myAC.grid_settings.y_max_grid_size=myAC.reference.b_ref/20;
myAC.grid_settings.aerodynamic_fuselage=0;

% set  structural grid size
myAC.grid_settings.dy_max_struct_grid=myAC.reference.b_ref/30;

%% grid computation
myAC=myAC.compute_grid();

%% create Structural Model
%options
structOpts=class_wingstructure_solver_settings;
%create
[myAC,myStr]=create_structural_model(myAC, structOpts);
myStr=myStr.f_solve_free_modes;
%% Init Coupling
myAC=myAC.compute_force_interpolation_matrix(myStr);
%% fligth state definition

ref_state=critical_ref_state(myAC,0.4,0); %Mach Number and Altitude

%% aeroelastic trim

aeSolverSettings=class_aeroelastic_solver_settings;
[myAC,myStr,myVLM,trimmedState]=trim_elastic_aircraft(myAC,myStr,ref_state,aeSolverSettings,0);


%% check result

figure; myVLM.plot_cp;
hold on; myAC.plot_grid;



