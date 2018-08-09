%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%
function obj = f_calc_stresses( obj, overwrite )
%F_CALC_STRESSES Summary of this function goes here
%   Detailed explanation goes here
    if obj.anisotropic
        obj=obj.f_calc_stress_strain_crossmod();
    else
        Mbx = abs(obj.node_loadings_loc(4:6:end));%+abs(wing.node_loadings(6:6:end)); %bending moment about x (extend.. 
        Mbz = abs(obj.node_loadings_loc(6:6:end));
        Mt = abs(obj.node_loadings_loc(5:6:end)); %torsional moment
        Qz = abs(obj.node_loadings_loc(3:6:end));
        Qx = abs(obj.node_loadings_loc(1:6:end));

        for i = 1:1:obj.nel
            if obj.is_sym
                MBX = max(Mbx(i)*0.5+Mbx(i+1)*0.5,Mbx(end-i+1)*0.5+Mbx(end-i)*0.5);
                MBZ = max(Mbz(i)*0.5+Mbz(i+1)*0.5,Mbz(end-i+1)*0.5+Mbz(end-i)*0.5);
                MT = max(Mt(i)*0.5+Mt(i+1)*0.5,Mt(end-i+1)*0.5+Mt(end-i)*0.5);
                QX = max(Qx(i)*0.5+Qx(i+1)*0.5,Qx(end-i+1)*0.5+Qx(end-i)*0.5);
                QZ = max(Qz(i)*0.5+Qz(i+1)*0.5,Qz(end-i+1)*0.5+Qz(end-i)*0.5);
            else
                MBX = Mbx(i)*0.5+Mbx(i+1)*0.5;
                MBZ = Mbz(i)*0.5+Mbz(i+1)*0.5;
                MT = Mt(i)*0.5+Mt(i+1)*0.5;
                QX = Qx(i)*0.5+Qx(i+1)*0.5;
                QZ = Qz(i)*0.5+Qz(i+1)*0.5;
            end
            if isa(obj,'class_wing')
                obj.beamelement(i).crosssection = obj.beamelement(i).crosssection.f_calc_stresses(MBX,MBZ,MT,QX,QZ,obj.loadcase_index,overwrite);
            end
        end
    end
end

