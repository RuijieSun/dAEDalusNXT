%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   This file is part of dAEDalus structures
%                   Copyright (C) 2011, Klaus Seywald
%     Author:   	Klaus Seywald
%                   klaus.seywald@mytum.de
%                   seywald@kth.se
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [obj] = f_calc_element_mass(obj)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    el_ndof=obj.el_ndof;
    %% calculate element stiffnes matrices and load vectors
    if el_ndof==6
         for i=1:obj.nel   
            %Calculate Element mass Matrices
            
            % Checks if the beam is anisotropic. If so, it runs the
            % anistropic versions of the elM_6dof and elM_lumped_6dof
            % functions.
            if obj.anisotropic == 1
                obj.beamelement(i)=obj.beamelement(i).elM_6dof_anisotropic();    
                obj.beamelement(i)=obj.beamelement(i).elM_lumped_6dof_anisotropic(); 
            else
                obj.beamelement(i)=obj.beamelement(i).elM_6dof();    
                obj.beamelement(i)=obj.beamelement(i).elM_lumped_6dof(); 
            end
         end
    else
        % number of ndof not implemented
    end
end
