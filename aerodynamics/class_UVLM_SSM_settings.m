%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%
classdef class_UVLM_SSM_settings
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        %if 1, outputs are splitted in steady/unsteady
        splitForce =1
        %if 1, sectional pressures are output
        sectPress=1
        %reduction method 'bpod', 'balred'
        redMethod='bpod'
        %reduction order (lpvkernel); if set to zero, no reduction is
        %carried out
        redOrder=0
    end
    
    methods
    end
    
end

