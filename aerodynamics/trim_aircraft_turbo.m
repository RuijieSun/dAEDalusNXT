function [ac, flightState,vlm]=trim_aircraft_turbo(ac, flightState, vlm, def, rollAlloc)

ac=ac.compute_deflected_grid(def);
ClTarget=flightState.get_Cl(ac.reference.S_ref);
%init alpha
a=flightState.aerodynamic_state.alpha;
%vertical trim id
vertTrimSurfId=find(strcmp(ac.control_surfaces,'elevator'));
dVert=flightState.aircraft_state.control_deflections{vertTrimSurfId};
if flightState.p==0 %vertical case
    Ctarget=[ClTarget 0];
   
    [C0,ac,vlm]=f_rigid_aero_vert(vlm,ac,def,a,dVert);
    C=C0;

    for i=1
        % aoa perturbation
        [Ca,ac,vlm]=f_rigid_aero_vert(vlm,ac,def,a+1,dVert);
        dClda=Ca(1)-C(1);
        dCmda=Ca(2)-C(2);
        % dVert perturbation
        [Cd,ac,vlm]=f_rigid_aero_vert(vlm,ac,def,a,dVert+1);
        dCldCs=Cd(1)-C(1);
        dCmdCs=Cd(2)-C(2);

        %linsolve Determination of next step
        dCdT=[dClda dCldCs; dCmda dCmdCs];
        tIncr=linsolve(dCdT,(Ctarget-C)');
        a=a+tIncr(1);
        dVert=dVert+tIncr(2);
        
        %solve with current guess
        [C,ac,vlm]=f_rigid_aero_vert(vlm,ac,def,a,dVert);
    end
elseif flightState.p~=0 %roll case
    
    Ctarget=[ClTarget 0 0 0];
    %lat trim id
    latTrimSurfId=find(strcmp(ac.control_surfaces,'rudder'));
    %lat value
    dLat=flightState.aircraft_state.control_deflections{latTrimSurfId};
    %init dRoll
    dRoll=rollAlloc*0;
    [C0,ac,vlm0]=f_rigid_aero_lat(vlm,ac,def,flightState,a,dVert,dLat,dRoll);
    C=C0;
    for i=1
        % aoa perturbation
        [Ca,ac,vlm]=f_rigid_aero_lat(vlm,ac,def,flightState,a+1,dVert,dLat,dRoll);
        dCa=Ca-C;
        
        % dVert perturbation
        [Cvert,ac,vlm]=f_rigid_aero_lat(vlm,ac,def,flightState,a,dVert+1,dLat,dRoll);
        dCdVert=Cvert-C;
        
        % dLat perturbation
        [Clat,ac,vlm]=f_rigid_aero_lat(vlm,ac,def,flightState,a,dVert,dLat+1,dRoll);
        dCdLat=Clat-C;
        
        % dRoll perturbation
        dRollPert=rollAlloc./norm(rollAlloc,1)*length(rollAlloc);
        [Croll,ac,vlm]=f_rigid_aero_lat(vlm,ac,def,flightState,a,dVert,dLat,dRoll+dRollPert);
        dCdRoll=Croll-C;
        
        %linsolve trim problem
        dCdT=[dCa' dCdVert' dCdLat' dCdRoll'];
        tIncr=linsolve(dCdT,(Ctarget-C)');
        a=a+tIncr(1);
        dVert=dVert+tIncr(2);
        dLat=dLat+tIncr(3);
        dRoll=dRoll+dRollPert*tIncr(4);
        
        %solve with current guess
        [C,ac,vlm]=f_rigid_aero_lat(vlm,ac,def,flightState,a,dVert,dLat,dRoll);
    end
    flightState.rollDef=dRoll;
end
% store required info
flightState.aerodynamic_state=vlm.state;
flightState.aircraft_state.control_deflections=ac.control_deflections;

%% aero functions
function [C,ac,vlm]=f_rigid_aero_vert(vlm,ac,def,a,d)
    %set control surfaces to given value in d
    ac=ac.f_set_control_surface('elevator',d);
    ac=ac.compute_grid();
    ac=ac.compute_deflected_grid(def);
    vlm=vlm.set_grid(ac.grid_deflected, ac.panels);
    %set alpha to given value in a
    state=vlm.state;
    state=state.set_alpha(a);
    vlm=vlm.f_set_state(state);

    vlm=vlm.f_solve_fast();

    C=[vlm.Cl, vlm.CM];
end

function [C,ac,vlm]=f_rigid_aero_lat(vlm,ac,def, flightState,a,dVert,dLat,dRoll)
    %set control surfaces to given value in dVert and dLat
    ac=ac.f_set_control_surface('elevator',dVert);
    ac=ac.f_set_control_surface('rudder',dLat);
    ac=ac.compute_grid();
    ac=ac.compute_deflected_grid(def);
    vlm=vlm.set_grid(ac.grid_deflected, ac.panels);
    %set control surfaces for roll by normal vector deflection in vlm
    for iRollCs=1:length(flightState.rollSurf)
        vlm=vlm.deflect_nvec(flightState.rollSurf(iRollCs), dRoll(iRollCs),2,1);
    end
    
    %set alpha to given value in a
    state=vlm.state;
    state=state.set_alpha(a);
    vlm=vlm.f_set_state(state);
    %
    vlm=vlm.f_solve_fast_roll(flightState.p);

    C=[vlm.Cl,  vlm.CL,vlm.CM,vlm.CN];
end
end