%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%
function [ref_state] =critical_ref_state(aircraft,MCruise,StructuralDesignAltitude,varargin)
%UNTITLED3 Summary of this function goes here
% varargin 1 and 2 : cg and weight
% varargin 3 : Maximum operating Mach Number (TAS)
% varargin 4 : Maximum Dive Number (TAS)
%based on: http://adg.stanford.edu/aa241/structures/placard.html
if nargin>=5
    cgIn=varargin{1};
    weightIn=varargin{2};
else
    cgIn=aircraft.reference.p_ref;
    weightIn=aircraft.weights.MTOW;
end


%aircraft_state=class_aircraft_state(aircraft,[0.083333333*aircraft.reference.c_ref 0 0]);

aircraft_state=class_aircraft_state(aircraft,cgIn,weightIn);
[~,a,press,rho_air] = atmosisa(StructuralDesignAltitude);  % To use single standard throughout the tool. 
% [rho_air,a,T,P,mu]=stdatmo(StructuralDesignAltitude);

V=[MCruise*a 0 0];
%the flight state is initialized at the speed (TAS) at which the maneuver
%has to be carried out (this is then used e.g. in VLM solvers)

if nargin>=4
    MD=varargin{3};
else
    Mc=1.06*MCruise; %maximum operating mach number, kroo
    MD=1.07*Mc; %design dive mach number, kroo
end
if(StructuralDesignAltitude<11000)
    VDTAS=MD*a;%*sqrt((288.15-0.0065*ref_state.Altitude)/288.15); % design dive TAS
    VCTAS=MCruise*a;%*sqrt((288.15-0.0065*ref_state.Altitude)/288.15); % design dive TAS
else
    VDTAS=MD*a;%*sqrt((288.15-0.0065*11000)/288.15); % design dive TAS
    VCTAS=MCruise*a;%*sqrt((288.15-0.0065*11000)/288.15); % design dive TAS
end
VCEAS=VCTAS*sqrt(rho_air/1.225);  % design dive EAS
VDEAS=VDTAS*sqrt(rho_air/1.225);  % design dive EAS 


ref_state=class_flight_state(StructuralDesignAltitude,VDTAS,0,0,aircraft_state);
ref_state.V_EAS=VDEAS;
ref_state.VD=VDTAS;
ref_state.VC=VCTAS;
ref_state.VC_EAS=VCEAS;
ref_state.VD_EAS=VDEAS;
ref_state.M_C=MCruise;
ref_state.M_D=MD;
end

