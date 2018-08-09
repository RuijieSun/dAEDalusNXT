%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%
function [aircraft,flight_state,wingaero ] = trim_aircraft_fast(aircraft,flight_state,wingaero,def)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
%thrust moment to be considered?
consider_thrust=1;

% set the deflections of all control surfaces
for i=1:length(flight_state.aircraft_state.control_deflections)
    %if (flight_state.aircraft_state.control_deflections{i}>=1E-4)
        aircraft=aircraft.f_set_control_surface(flight_state.aircraft_state.control_surfaces{i},flight_state.aircraft_state.control_deflections{i});
   % end
end

%init observation and convergence variables

elev=0;
d_elev=-8;
no_redeform=0;
trim_name=[];
trim_name2=[];
trim_name3=[];
trim_name4=[];

%set trim ref to current cg
trim_ref=aircraft.reference;
trim_ref.p_ref=flight_state.aircraft_state.CG_ref;

trim_idx=1;

% find elevator
if isempty(aircraft.control.trim_surfaces)
    trim_name='elevator';
    for i=1:length(aircraft.control_surfaces)
        if strcmp(trim_name,aircraft.control_surfaces{i});
            trim_idx=i;
        end
    end
else
    trim_name=aircraft.control.trim_surfaces{1};
    for i=1:length(aircraft.control_surfaces)
        if strcmp(trim_name,aircraft.control_surfaces{i});
            trim_idx=i;
        end
    end
end

elev=aircraft.control_deflections{trim_idx};
if elev~=0
    d_elev=-0.1;
end
if length(aircraft.control.trim_surfaces)==2
    trim_name2=aircraft.control.trim_surfaces{2};
end

if length(aircraft.control.trim_surfaces)==4
    trim_name2=aircraft.control.trim_surfaces{2};
    trim_name3=aircraft.control.trim_surfaces{3};
    trim_name4=aircraft.control.trim_surfaces{4};
end
%set flight state cg to current cg
flight_state.aerodynamic_state.p_ref=flight_state.aircraft_state.CG_ref;
aircraft=aircraft.f_set_control_surface(trim_name,elev);
if ~isempty(trim_name2)
   aircraft=aircraft.f_set_control_surface(trim_name2,-elev); 
end
if ~isempty(trim_name3)
   aircraft=aircraft.f_set_control_surface(trim_name3,elev); 
end

if ~isempty(trim_name4)
   aircraft=aircraft.f_set_control_surface(trim_name4,elev); 
end

%recompute grid
aircraft=aircraft.compute_grid();
if isempty(def)
    %set grid
    wingaero=wingaero.set_grid(aircraft.grid, aircraft.panels);
else
    aircraft=aircraft.compute_deflected_grid(def);	
    %set grid
    wingaero=wingaero.set_grid(aircraft.grid_deflected, aircraft.panels);
end
%solve quickly
wingaero=wingaero.f_solve_for_Cl_fast(flight_state.get_Cl(aircraft.reference.S_ref));


%compute drag to compute thrust moment
aircraft=aircraft.compute_CD_f(flight_state.aerodynamic_state,aircraft.reference.S_ref);
if consider_thrust==1
    wingaero=wingaero.f_solve_full();
    T_pe=((aircraft.CD_f+wingaero.Cdi)*1/2*flight_state.aerodynamic_state.rho_air*norm(wingaero.Uinf)^2*aircraft.reference.S_ref)/length(aircraft.engines)/cos(atan(wingaero.Uinf(3)/wingaero.Uinf(1)));   % Mistake in resolution of forces corrected - 3 occurences
    delta_CM=0;
    delta_Cl=0;
    for i=1:length(aircraft.engines)
        aircraft.engines(i).delta_t=T_pe/aircraft.engines(i).thrust;
        delta_CM=delta_CM-(aircraft.engines(i).cg_pos(3)-flight_state.aerodynamic_state.p_ref(3))*T_pe/(1/2*flight_state.aerodynamic_state.rho_air*norm(wingaero.Uinf)^2*aircraft.reference.S_ref*aircraft.reference.c_ref);
        delta_Cl=delta_Cl+T_pe*sin(atan(wingaero.Uinf(3)/wingaero.Uinf(1)))/(1/2*flight_state.aerodynamic_state.rho_air*norm(wingaero.Uinf)^2*aircraft.reference.S_ref);
    end
else
    delta_CM=0;
    delta_Cl=0;
end
CM0=wingaero.CM+delta_CM;

%deflect elevator with additional increment
elev=elev+d_elev;
aircraft=aircraft.f_set_control_surface(trim_name,elev);
if ~isempty(trim_name2)
   aircraft=aircraft.f_set_control_surface(trim_name2,-elev); 
end
if ~isempty(trim_name3)
   aircraft=aircraft.f_set_control_surface(trim_name3,elev); 
end
if ~isempty(trim_name4)
   aircraft=aircraft.f_set_control_surface(trim_name4,elev); 
end


%recompute grid (because of control deflections)
aircraft=aircraft.compute_grid();

if isempty(def)
    %set grid
    wingaero=wingaero.set_grid(aircraft.grid, aircraft.panels);
else
    aircraft=aircraft.compute_deflected_grid(def);	
    %set grid
    wingaero=wingaero.set_grid(aircraft.grid_deflected, aircraft.panels);
end


%solve quickly
wingaero=wingaero.f_solve_for_Cl_fast(flight_state.get_Cl(aircraft.reference.S_ref));

%compute drag to compute thrust moment
if consider_thrust==1
    wingaero=wingaero.f_solve_full();
    T_pe=((aircraft.CD_f+wingaero.Cdi)*1/2*flight_state.aerodynamic_state.rho_air*norm(wingaero.Uinf)^2*aircraft.reference.S_ref)/length(aircraft.engines)/cos(atan(wingaero.Uinf(3)/wingaero.Uinf(1)));
    delta_CM=0;
    delta_Cl=0;
    for i=1:length(aircraft.engines)
        aircraft.engines(i).delta_t=T_pe/aircraft.engines(i).thrust;
        delta_CM=delta_CM-(aircraft.engines(i).cg_pos(3)-flight_state.aerodynamic_state.p_ref(3))*T_pe/(1/2*flight_state.aerodynamic_state.rho_air*norm(wingaero.Uinf)^2*aircraft.reference.S_ref*aircraft.reference.c_ref);
        delta_Cl=delta_Cl+T_pe*sin(atan(wingaero.Uinf(3)/wingaero.Uinf(1)))/(1/2*flight_state.aerodynamic_state.rho_air*norm(wingaero.Uinf)^2*aircraft.reference.S_ref);
    end
else
    delta_CM=0;
    delta_Cl=0;
end
    
    
CM0p=wingaero.CM+delta_CM;

dCMelev=(CM0p-CM0)/d_elev;
k=1;
while(abs(wingaero.CM+delta_CM)>1E-5)
    CM_prev=wingaero.CM+delta_CM;
    d_elev=-CM_prev/dCMelev;
    elev=elev+d_elev;
    aircraft=aircraft.f_set_control_surface(trim_name,elev);
    if ~isempty(trim_name2)
        aircraft=aircraft.f_set_control_surface(trim_name2,-elev);
    end
    if ~isempty(trim_name3)
        aircraft=aircraft.f_set_control_surface(trim_name3,elev);
    end
    if ~isempty(trim_name4)
        aircraft=aircraft.f_set_control_surface(trim_name4,elev);
    end
   
    %recompute grid
    aircraft=aircraft.compute_grid();
 
    if isempty(def)
        %set grid
        wingaero=wingaero.set_grid(aircraft.grid, aircraft.panels);
    else
        aircraft=aircraft.compute_deflected_grid(def);	
        %set grid
        wingaero=wingaero.set_grid(aircraft.grid_deflected, aircraft.panels);
    end
    %solve quickly
    wingaero=wingaero.f_solve_for_Cl_fast(flight_state.get_Cl(aircraft.reference.S_ref));
    
    if consider_thrust==1
        wingaero=wingaero.f_solve_full();
        T_pe=((aircraft.CD_f+wingaero.Cdi)*1/2*flight_state.aerodynamic_state.rho_air*norm(wingaero.Uinf)^2*aircraft.reference.S_ref)/length(aircraft.engines)/cos(atan(wingaero.Uinf(3)/wingaero.Uinf(1)));
        delta_CM=0;
        delta_Cl=0;
        for i=1:length(aircraft.engines)
            aircraft.engines(i).delta_t=T_pe/aircraft.engines(i).thrust;
            delta_CM=delta_CM-(aircraft.engines(i).cg_pos(3)-flight_state.aerodynamic_state.p_ref(3))*T_pe/(1/2*flight_state.aerodynamic_state.rho_air*norm(wingaero.Uinf)^2*aircraft.reference.S_ref*aircraft.reference.c_ref);
            delta_Cl=delta_Cl+T_pe*sin(atan(wingaero.Uinf(3)/wingaero.Uinf(1)))/(1/2*flight_state.aerodynamic_state.rho_air*norm(wingaero.Uinf)^2*aircraft.reference.S_ref);
        end
        
    else
        delta_CM=0;
        delta_Cl=0;
    end
        
    dCMelev=(wingaero.CM+delta_CM-CM_prev)/d_elev;
    k=k+1;
    if k>10
        fprintf('Trimming not converged\n');
    end
    flight_state.aerodynamic_state=wingaero.state;
end

flight_state.aerodynamic_state=wingaero.state;
for i=1:length(flight_state.aircraft_state.control_surfaces)
    if strcmp(flight_state.aircraft_state.control_surfaces{i},trim_name)
        flight_state.aircraft_state.control_deflections{i}=elev;
    end
    if strcmp(flight_state.aircraft_state.control_surfaces{i},trim_name2) 
        flight_state.aircraft_state.control_deflections{i}=-elev;
    end
    
    if ~isempty(trim_name3)
         if strcmp(flight_state.aircraft_state.control_surfaces{i},trim_name3) 
        aircraft=aircraft.f_set_control_surface(trim_name3,elev);
         end
    end
    if ~isempty(trim_name4)
        if  strcmp(flight_state.aircraft_state.control_surfaces{i},trim_name4) 

        aircraft=aircraft.f_set_control_surface(trim_name4,elev);
                end
    end
end

flight_state.aircraft_state.control_deflections=aircraft.control_deflections;
 if consider_thrust==1
     flight_state.aircraft_state.engine_thrust=T_pe;
 end
end

