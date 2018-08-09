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
% f_initMaterialProperties :   file is part of nlFEM class_wing
%               initialize FEM model with required material properties
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function obj = f_init_material_properties(obj, structure, anisotropic_flag, layup_settings)

%%set Material properties

% anisotropic flag checks whether the wingbox is made of anisotropic
% materials. If that is the case, the code below will deal with the
% anisotropic version of class_beam_element and class_material.
switch anisotropic_flag
    case 0
        for i=1:obj.nel
            %% todo now all has Uskin.E and Uskin.G
            obj.beamelement(i).E=structure.materials.E;
            obj.beamelement(i).G=structure.materials.G;
            %Upper Skin
            E_sk_up = structure.materials.E;
            G_sk_up = structure.materials.G;
            rho_sk_up = structure.materials.rho;
            %Lower Skin
            E_sk_lo = structure.materials.E;
            G_sk_lo = structure.materials.G;
            rho_sk_lo = structure.materials.rho;
            %Spars
            rho_sp = structure.materials.rho;
            sigma_allowable_sp = structure.materials.sigma_allowable;
            sigma_allowable_sk_u=structure.materials.sigma_allowable;
            sigma_allowable_sk_l=structure.materials.sigma_allowable;
            
            
            t_min_sp=structure.t_min_sp;
            t_min_sk=structure.t_min_sk;
            
            obj.fuel_density=structure.fuel_density;
            
            obj.beamelement(i).crosssection=obj.beamelement(i).crosssection.setMaterial(E_sk_up,G_sk_up,rho_sk_up,...
                E_sk_lo,G_sk_lo,rho_sk_lo,rho_sp,sigma_allowable_sp,sigma_allowable_sk_u,sigma_allowable_sk_l,t_min_sp,t_min_sk);
            obj.beamelement(i)=obj.beamelement(i).f_calcCrossProp;
        end
        
    case 1
        % TEMPORARY CODE used to create linear distribution of
        % thickness between tip and root.
        %             thick_temp = linspace(0.0015, 0.005, obj.nel/2);
        %             thick_temp = [thick_temp, thick_temp(end:-1:1)];
        
        % Test if the wing is symmetrical to symmetrize the laminate as
        % well
        
        if obj.is_sym == 1
            for i=1:obj.nel/2
                
                obj.fuel_density=structure.fuel_density;
                
                % All the material related information is stored within
                % crosssection objects, since it will vary depending on the
                % particular skin or spar we are looking at
                obj.beamelement(i).crosssection=obj.beamelement(i).crosssection.setMaterial(structure,layup_settings);
            end

            layup_settings_sym = layup_settings;
            layup_settings_sym.skins_layup_angles = - layup_settings_sym.skins_layup_angles;
            
            for i=obj.nel/2+1:obj.nel
                
                obj.fuel_density=structure.fuel_density;
                
                % All the material related information is stored within
                % crosssection objects, since it will vary depending on the
                % particular skin or spar we are looking at
                obj.beamelement(i).crosssection=obj.beamelement(i).crosssection.setMaterial(structure,layup_settings_sym);
            end
        
        else
             for i=1:obj.nel
                
                obj.fuel_density=structure.fuel_density;
                
                % All the material related information is stored within
                % crosssection objects, since it will vary depending on the
                % particular skin or spar we are looking at
                obj.beamelement(i).crosssection=obj.beamelement(i).crosssection.setMaterial(structure,layup_settings);
            end
            
        end
        
end
end