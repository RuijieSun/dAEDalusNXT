%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%
classdef class_wingstructure_solver_settings<class_structure_solver_settings
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        fuel_mass=0;
        engine=0;
        landing_gear=0;  
        
        % information used to discretize the cross section of each beam
        % element into shell elements. One of the two should always be set
        % to zero, so that the code will automatically use the other
        % variable. The code DOES NOT check which condition is
        % dominant.
        
        nr_nodes_crosssection_element=30;  % nr of nodes per skin/spar
        ds_nodes_crosssection_element=0;  % max shell element width
    end
    
    methods
    end
    
end

