%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%
function [aircraft,flight_state,wingaero ] = trim_wing_faster(aircraft,flight_state,wingaero)
 



% disp(['     	trimming aircraft for Cl: ' num2str(flight_state.get_Cl(aircraft.reference.S_ref))]);
wingaero=wingaero.f_solve_for_Cl_fast(flight_state.get_Cl(aircraft.reference.S_ref));

flight_state.aerodynamic_state=wingaero.state;


end