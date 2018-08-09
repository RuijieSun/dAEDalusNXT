%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%
function [aircraft,flight_state,wingaero ] = trim_wing_fast(aircraft,flight_state,wingaero)
 
%set new reference
trim_ref=aircraft.reference;
trim_ref.p_ref=flight_state.aircraft_state.CG_ref;

%set new cg
flight_state.aerodynamic_state.p_ref=flight_state.aircraft_state.CG_ref;
wingaero.state.p_ref=flight_state.aircraft_state.CG_ref;

%set the new grid
wingaero=wingaero.set_grid(aircraft.grid_deflected, aircraft.panels);

%solve quickly
wingaero=wingaero.f_solve_for_Cl_fast(flight_state.get_Cl(aircraft.reference.S_ref));

%hand over resulting state
flight_state.aerodynamic_state=wingaero.state;

end