%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%
%> @file class_crosssection_wingbox.m
%> @brief File containing the class for a finite wingbox crosssection
%>        This file is part of dAEDalus structures, Copyright (C) 2011,
%>        Klaus Seywald
%>   dAEDalus is published under the terms of the GNU General Public
%>   License by the Free Software Foundation;
%>   either version 2 or any later version.
%>
%>   You should have received a copy of the GNU General Public
%>   License along with dAEDalus; see the file GNU GENERAL 
%>   PUBLIC LICENSE.TXT.  If not, write to the Free Software 
%>   Foundation, 59 Temple Place -Suite 330, Boston, MA
%>   02111-1307, USA.
%>==================================================================

classdef (Abstract) class_crosssection
    
    methods (Abstract)

        obj = setGeometry(obj, varargin)
       
        obj = setMaterial(obj, varargin)
        
        dimensions = get_dimensions(obj)
        
        dm = f_calc_dm(obj)
        
        obj = f_calc_stresses(obj, Mbx, Mby, Mt, Qx, Qz, loadcase_idx, overwrite)
        
        obj = f_self_design_crosssection(obj, Mbx, Mby, Mt, Qx, Qz, loadcase_idx, overwrite)
        
        [Ix, Iy, Iz, J, A, Aenclosed] = calc_crosssection(obj)
    end   
end

