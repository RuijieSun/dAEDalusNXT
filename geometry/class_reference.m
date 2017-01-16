%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%
classdef class_reference
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        S_ref;
        c_ref;
        b_ref;
        p_ref;
        
    end
    
    methods
        
        function obj=class_reference(S_ref,c_ref,b_ref,p_ref)
            obj.S_ref=S_ref;
            obj.c_ref=c_ref;
            obj.b_ref=b_ref;
            obj.p_ref=p_ref; 
        end
        
    end
    
end

