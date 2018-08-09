%% generate required paths for program execution
addpath(genpath('..\..\..\'));


%% load XML File
myAC=class_aircraft('part7-simple787.xml',1);

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
refState=critical_ref_state(myAC,0.4,0); %Mach Number and Altitude


%% aerostructural layout
%define critical state
criticalState=critical_g_maneuver_state(refState,2.5);
aeSolverSettings=class_aeroelastic_solver_settings;
[myAC,myStr,myVLM]= structural_sizing_loop(myAC,myStr,criticalState,myAC.weights,aeSolverSettings);

%% check result
figure;
myStr.plot_structure

%% add more cases
criticalStates(1)=critical_g_maneuver_state(refState,2.5);
criticalStates(1).loadcase_index=1;
criticalStates(2)=critical_g_maneuver_state(refState,-1.0);
criticalStates(2).loadcase_index=2;
criticalStates(3)=critical_sideslip_maneuver_state(refState,10);
criticalStates(3).loadcase_index=3;
criticalStates(4)=critical_aileron_roll_load_state(refState,2.5*2/3,'aileron',-25);
criticalStates(4).loadcase_index=4;
[myAC,myStr,myVLM2]= critical_case_layout(myAC,myStr,criticalStates(1:4),myAC.weights,aeSolverSettings);

%% check result
figure;
myStr.plot_structure
