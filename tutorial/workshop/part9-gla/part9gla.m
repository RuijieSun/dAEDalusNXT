%% generate required paths for program execution
addpath(genpath('..\..\..\'));


%% load XML File
myAC=class_aircraft('part9-simple787.xml',1);

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

%% 1g trim
refState=critical_ref_state(myAC,0.4,0); %Mach Number and Altitude

[myAC,myStr,myVLM,trimmedState]=trim_elastic_aircraft(myAC,myStr,refState,aeSolverSettings,0);
nomLoads=myStr.beam(1).node_loadings_loc;

%% create aeSSM
%aircraft needs UVLM grid information (wake)
myAC.grid_settings.wake=2;
myAC=myAC.compute_grid;

% initialize aeSSM
aeSSMSettings=AeSSMSettings;
aeSSMSettings.modes=[7:16]; % specify elastic modes
aeSSMSettings.wakeLengthFactor=0.3; %performance improvement, better set to >1
aeSSMSettings.strLoadStations=[1:61];
aeSSMSettings.strPropDamping=0.02;
myAeSSM=AeSSM(myAC,myStr,refState.aerodynamic_state, aeSSMSettings);


%% determining gust loads in open and closed loop:

myAeSSM=AeSSM(myAC,myStr,refState.aerodynamic_state, aeSSMSettings);

% run open loop gusts
myAeSSM=myAeSSM.runGustAnalysis(refState,1, 5);
myAeSSM=myAeSSM.runContTurbAnalysis(refState,1);
olGustData=myAeSSM.gustData;

% specify gla surfaces in gustState
refState.glaSurf={'aileron1','aileron2','elevator'};
refState.glaUB=[25 25 25];
refState.glaLB=[-25 -25 -25];

refState.glaSurf={'aileron1'};
refState.glaUB=[25 ];
refState.glaLB=[-25 ];

refState.glaSurf={'aileron2','elevator'};
refState.glaUB=[ 25 25];
refState.glaLB=[ -25 -25];
%% gla synthesis 
myAeSSM=myAeSSM.optimFFController(refState,1,10,nomLoads);

% rerun gust/contTurb with gla
myAeSSM=myAeSSM.runGustAnalysis(refState,1, 5);
myAeSSM=myAeSSM.runContTurbAnalysis(refState,1);
clGustData=myAeSSM.gustData;
%%
loadIdx=(myStr.beam(1).nel/2)*6+4;
figure;
hold on; 
title('Discrete Gust Loads New Implementation')
plot(olGustData.timeVec,  nomLoads(loadIdx)+ squeeze(olGustData.loads(:,loadIdx,:)),'--')
plot(olGustData.timeVec,  nomLoads(loadIdx)-  squeeze(olGustData.loads(:,loadIdx,:)),'--')
plot(clGustData.timeVec,  nomLoads(loadIdx)+ squeeze(clGustData.loads(:,loadIdx,:)))
plot(clGustData.timeVec,  nomLoads(loadIdx)-  squeeze(clGustData.loads(:,loadIdx,:)))
xlabel('time [s]')
ylabel('WRBM [Nm]')

%%

figure;
hold on; 
title('Continuous Gust WRBM response New Implementation')

plot(olGustData.contTurb.timeVec, (nomLoads(loadIdx))'+olGustData.contTurb.loads(:,loadIdx),'b')
plot(olGustData.contTurb.timeVec, (nomLoads(loadIdx))'-olGustData.contTurb.loads(:,loadIdx),'b')
plot(clGustData.contTurb.timeVec, (nomLoads(loadIdx))'+clGustData.contTurb.loads(:,loadIdx),'r')
plot(clGustData.contTurb.timeVec, (nomLoads(loadIdx))'-clGustData.contTurb.loads(:,loadIdx),'r')

xlabel('time [s]')
ylabel('WRBM [Nm]')
legend('GLA off','GLA off','GLA on','GLA on')

%%
figure;
hold on; 
title('Continuous Gust Loads New Implementation')

plot(abs(nomLoads(4:6:end))'+rms(olGustData.contTurb.loads(:,4:6:end)),'b')
plot(abs(nomLoads(4:6:end))'+rms(clGustData.contTurb.loads(:,4:6:end)),'r')

xlabel('beam node')
ylabel('Bending Moment [Nm]')
legend('GLA off','GLA on')


%% plot commanded def and approximated rate
myAeSSM.ctr.getCommands(squeeze(myAeSSM.gustData.inputs(:,1,:)),myAeSSM.gustData.timeVec)

%% plot actuator def an rate
myAeSSM.plotActDeflectionAndRates(squeeze(myAeSSM.gustData.inputs(:,1,:)),myAeSSM.gustData.timeVec)

%% plot ctr bode
myAeSSM.ctr.bode