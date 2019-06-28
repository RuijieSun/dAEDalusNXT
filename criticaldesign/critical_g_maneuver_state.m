%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%
function [gmaneuverstate] =critical_g_maneuver_state(ref_state,gfactor, varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
        gmaneuverstate=ref_state;
    if isempty(varargin)
        gmaneuverstate.V=[ref_state.VD 0 0];
        gmaneuverstate.aerodynamic_state.V_A=ref_state.VD;
        gmaneuverstate.aerodynamic_state.Ma=ref_state.M_D;
        gmaneuverstate.aerodynamic_state.V_inf=[ref_state.VD 0 0];
    elseif strcmp(varargin{1},'VC')
        gmaneuverstate.V=[ref_state.VC 0 0];
        gmaneuverstate.aerodynamic_state.V_A=ref_state.VC;
        gmaneuverstate.aerodynamic_state.Ma=ref_state.M_C;
        gmaneuverstate.aerodynamic_state.V_inf=[ref_state.VC 0 0];
    end
        gmaneuverstate.load_factor=gfactor;
end

