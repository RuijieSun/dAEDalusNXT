%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%
% This function calculates the stiffness matrix of objects belonging to
% class_beamelement_anisotropic.

% It is a re-adaptation of lin_elK_6dof. This one is able to handle
% anisotropic materials (as well as isotropic ones of course, provided that
% they are properly defined).

function obj=lin_elK_6dof_anisotropic(obj)


    %% New local stiffness matrix definition
    %BEFORE THE CALCULATION, EIzx, AT, LT and FT needs to be defined.
    %AT=coupled axial-torsion stiffness
    %FT=coupled flap-torsion stiffness
    %LT=coupled lag-torsion stiffness
    % Note how all variables below belonging to the crosssection object, as
    % they were calculated by cross sectional modeller
    
    %ref: ?
    EA = obj.crosssection.Se(1,1);
%     AF = obj.crosssection.Se(1,2);
%     AL = obj.crosssection.Se(1,3);
    AT = obj.crosssection.Se(1,4);
    
    EIy = obj.crosssection.Se(2,2);
    EIyz = obj.crosssection.Se(2,3);
    FT = obj.crosssection.Se(2,4);
    
    EIz = obj.crosssection.Se(3,3);
    LT = obj.crosssection.Se(3,4);
    
    GJ = obj.crosssection.Se(4,4);
    
    L = obj.le;
    
   %%
   %anisotropic beam stiffness matrix (ref: Lopez, Ricardo de Frias; A 3D Finite Beam Element for the Modelling of Composite Wind Turbine Wings)
     K1=[   EA/L    0               0               AT/L    0               0;... 
            0       12*EIz/L^3      12*EIyz/L^3     0       -6*EIyz/L^2     6*EIz/L^2;... 
            0       12*EIyz/L^3     12*EIy/L^3      0       -6*EIy/L^2      6*EIyz/L^2;... 
            AT/L    0               0               GJ/L    FT/L            -LT/L;... 
            0       -6*EIyz/L^2     -6*EIy/L^2      FT/L    4*EIy/L         -4*EIyz/L;... 
            0       6*EIz/L^2       6*EIyz/L^2      -LT/L   -4*EIyz/L       4*EIz/L]; 
    K2=[    -EA/L   0               0               -AT/L   0               0;... 
            0       -12*EIz/L^3     -12*EIyz/L^3    0       -6*EIyz/L^2     6*EIz/L^2;... 
            0       -12*EIyz/L^3    -12*EIy/L^3     0       -6*EIy/L^2      6*EIyz/L^2;... 
            -AT/L   0               0               -GJ/L   -FT/L           LT/L;... 
            0       6*EIyz/L^2      6*EIy/L^2       -FT/L   2*EIy/L         -2*EIyz/L;... 
            0       -6*EIz/L^2      -6*EIyz/L^2     LT/L    -2*EIyz/L       2*EIz/L]; 
    K3=[    EA/L    0               0               AT/L    0               0;... 
            0       12*EIz/L^3      12*EIyz/L^3     0       6*EIyz/L^2      -6*EIz/L^2;... 
            0       12*EIyz/L^3     12*EIy/L^3      0       6*EIy/L^2       -6*EIyz/L^2;... 
            AT/L 	0               0               GJ/L    FT/L            -LT/L;... 
            0       6*EIyz/L^2      6*EIy/L^2       FT/L    4*EIy/L         -4*EIyz/L;... 
            0       -6*EIz/L^2      -6*EIyz/L^2     -LT/L   -4*EIyz/L       4*EIz/L]; 
    K=[K1 K2; K2' K3];
    %%
    %change  x -> y and y -> x axis to get the stiffness matrix in the beam element coordinate system of dAEDalus
    obj.elK=K([2 1 3 5 4 6 8 7 9 11 10 12],[2 1 3 5 4 6 8 7 9 11 10 12]);
    %flip x axis
    obj.elK=kron(ones(2,2),[-1 1 1 -1 1 1]'*[-1 1 1 -1 1 1]).*obj.elK;
    %Calculate Element Stiffness Matrices in Global Coordinates
    obj.elKglobal=obj.T'*obj.elK*obj.T;
    
end
