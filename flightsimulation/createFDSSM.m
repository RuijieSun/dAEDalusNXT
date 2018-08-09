%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%
function [sys] =  createFDSSM(m,Vb_0,Euler_0,omegab_0, Inertia, g)
    u_0=Vb_0(1);
    v_0=Vb_0(2);
    w_0=Vb_0(3);

    phi_0=Euler_0(1);
    theta_0=Euler_0(2);
    psi_0=Euler_0(3);

    p_0=omegab_0(1);
    q_0=omegab_0(2);
    r_0=omegab_0(3);




    nabThetaRDot0=[ w_0*cos(phi_0)*sin(psi_0) + v_0*sin(phi_0)*sin(psi_0) + v_0*cos(phi_0)*cos(psi_0)*sin(theta_0) - w_0*cos(psi_0)*sin(phi_0)*sin(theta_0), cos(psi_0)*(w_0*cos(phi_0)*cos(theta_0) - u_0*sin(theta_0) + v_0*cos(theta_0)*sin(phi_0)), w_0*cos(psi_0)*sin(phi_0) - v_0*cos(phi_0)*cos(psi_0) - u_0*cos(theta_0)*sin(psi_0) - w_0*cos(phi_0)*sin(psi_0)*sin(theta_0) - v_0*sin(phi_0)*sin(psi_0)*sin(theta_0);
     v_0*cos(phi_0)*sin(psi_0)*sin(theta_0) - v_0*cos(psi_0)*sin(phi_0) - w_0*cos(phi_0)*cos(psi_0) - w_0*sin(phi_0)*sin(psi_0)*sin(theta_0), sin(psi_0)*(w_0*cos(phi_0)*cos(theta_0) - u_0*sin(theta_0) + v_0*cos(theta_0)*sin(phi_0)), u_0*cos(psi_0)*cos(theta_0) - v_0*cos(phi_0)*sin(psi_0) + w_0*sin(phi_0)*sin(psi_0) + w_0*cos(phi_0)*cos(psi_0)*sin(theta_0) + v_0*cos(psi_0)*sin(phi_0)*sin(theta_0);
    cos(theta_0)*(v_0*cos(phi_0) - w_0*sin(phi_0)),            - u_0*cos(theta_0) - w_0*cos(phi_0)*sin(theta_0) - v_0*sin(phi_0)*sin(theta_0),       0];

    nabVRDot0=[ cos(psi_0)*cos(theta_0), cos(psi_0)*sin(phi_0)*sin(theta_0) - cos(phi_0)*sin(psi_0), sin(phi_0)*sin(psi_0) + cos(phi_0)*cos(psi_0)*sin(theta_0);
     cos(theta_0)*sin(psi_0), cos(phi_0)*cos(psi_0) + sin(phi_0)*sin(psi_0)*sin(theta_0), cos(phi_0)*sin(psi_0)*sin(theta_0) - cos(psi_0)*sin(phi_0);
               -sin(theta_0),                                    cos(theta_0)*sin(phi_0),                                    cos(phi_0)*cos(theta_0)];

    nabEulerEulerDot0=[ (sin(theta_0)*(q_0*cos(phi_0) - r_0*sin(phi_0)))/cos(theta_0),                (r_0*cos(phi_0) + q_0*sin(phi_0))/cos(theta_0)^2, 0;
                                 - r_0*cos(phi_0) - q_0*sin(phi_0),                                                               0, 0;
                    (q_0*cos(phi_0) - r_0*sin(phi_0))/cos(theta_0), (sin(theta_0)*(r_0*cos(phi_0) + q_0*sin(phi_0)))/cos(theta_0)^2, 0];

    nabOmegaDotEulerDot0=[ 1, (sin(phi_0)*sin(theta_0))/(cos(theta_0)*cos(phi_0)^2 + cos(theta_0)*sin(phi_0)^2), (cos(phi_0)*sin(theta_0))/(cos(theta_0)*cos(phi_0)^2 + cos(theta_0)*sin(phi_0)^2);
     0,                                          cos(phi_0)/(cos(phi_0)^2 + sin(phi_0)^2),                                         -sin(phi_0)/(cos(phi_0)^2 + sin(phi_0)^2);
     0,                sin(phi_0)/(cos(theta_0)*cos(phi_0)^2 + cos(theta_0)*sin(phi_0)^2),                cos(phi_0)/(cos(theta_0)*cos(phi_0)^2 + cos(theta_0)*sin(phi_0)^2)];

    G0=g*[0 -cos(theta_0) 0; cos(theta_0)*cos(phi_0) -sin(theta_0)*sin(phi_0) 0; -cos(theta_0)*sin(phi_0) -sin(theta_0)*cos(phi_0) 0];

    omegab0Skew=skewsym(omegab_0);
    v0Skew=skewsym(Vb_0);
    strangeterm=skewsym(Inertia*omegab_0);
    InertTens=-Inertia^-1*(omegab0Skew*Inertia-strangeterm);

    A=[ zeros(3) nabThetaRDot0 nabVRDot0 zeros(3) ;
        zeros(3)  nabEulerEulerDot0 zeros(3)  nabOmegaDotEulerDot0;
        zeros(3)  G0 -omegab0Skew v0Skew;
        zeros(3)  zeros(3)  zeros(3)  InertTens];
    B=[ zeros(3)  zeros(3) ;
        zeros(3)  zeros(3) ;
        eye(3)/m  zeros(3) ;
        zeros(3)  Inertia^-1];


    C=[eye(12);A];
    D=[zeros(12,6); B];


    sys.a=A;
    sys.b=B;
    sys.c=C;
    sys.d=D;
    sys.StateName={'rI_x';'rI_y';'rI_z';'Phi'; 'Theta'; 'Psi'; 'vB_x'; 'vB_y'; 'vB_z'; 'wB_x'; 'wB_y'; 'wB_z';};
    sys.InputName={'X';'Y';'Z';'L';'M';'N'};
    sys.InputGroup.RBMForces=1:3;
    sys.InputGroup.RBMMoments=4:6;
    sys.OutputName=[sys.StateName;{ 'vI_x';'vI_y';'vI_z';'PhiDot'; 'ThetaDot'; 'PsiDot'; 'vBDot_x'; 'vBDot_y'; 'vBDot_z'; 'wBDot_x'; 'wBDot_y'; 'wBDot_z';}];
    sys.OutputGroup.rI=1:3;
    sys.OutputGroup.Euler=4:6;
    sys.OutputGroup.vB=7:9;
    sys.OutputGroup.wB=10:12;
    sys.OutputGroup.vI=13:15;
    sys.OutputGroup.EulerDot=16:18;
    sys.OutputGroup.vBDot=19:21;
    sys.OutputGroup.wBDot=22:24;
    
end