%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%
classdef class_flight_state
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
        V % set maneuvering speed
        V_EAS% set maneuvering speed(equivalent airspeed)
        
        %> Roll rate in body COS                          double      (1x1)   [rad/s]*
        p = 0    
        %> Pitch rate in body COS                         double      (1x1)   [rad/s]*
        q = 0 
        %> Yaw rate in body COS                           double      (1x1)   [rad/s]*  
        r = 0                       
        
        %> simulation starting time
        t0=0;
        %> time of the actual state
        ts=0;               
        %> aerodynamic state
        aerodynamic_state;
        %> 
        aircraft_state;
        %> load factor
        load_factor=1;
        %> altitude
        h=0;
        
        % design speeds
        M_D;
        M_C;
        VD; %design dive speed
        VC; % design cruise speed
        VC_EAS; %design cruise speed (equivalent airspeed)
        VD_EAS; % design dive speed (equivalent airspeed)
        loadcase_index;     % loadcase index of this state (for visualization)
        % body states
        
        reference;
        % mla surfaces (array of names)
        mlaSurf;
        % mla setting - if not empty, no mla optimization is run and the
        % mlaVal value is used 
        mlaVal;
        % mla upper bounds
        mlaUB;
        % mla lower bounds
        mlaLB;
        % mla deflection symmetry (e.g. deflect aileron symmetrically)
        mlaSymDefl;
        %deformation state (used for optimization, stored within iterations
        %to speed up convergence)
        def;
        % gla surfaces (array of names)
        glaSurf;
        % gla upper bounds (deflection)
        glaUB;
        % gla lower bounds (deflection)
        glaLB;        
        % gla setting - if not empty, no gla optimization is done and the
        % glaVal is used 
        glaVal;
        % gla controller time step 
        glaTs=[]
        % roll surfaces (array of names)
        rollSurf
        % deflection used for roll 
        rollDef
        % rollAllocVal - if not empty, no rollalloc optimization is done
        % and the rollAllocVal is used
        rollAllocVal
        
    end
    
    methods
        function obj=class_flight_state(H,V,alpha,beta,aircraft_state)     
           obj.h=H;
           obj.V=V;
           [rho_air,a,T,P,mu]=stdatmo(H);
           Uinf=norm(V);
           Ma=Uinf/a;
           obj.aerodynamic_state=class_aero_state(Uinf,alpha,beta,Ma,rho_air,mu);
           obj.aircraft_state=aircraft_state;
        end
        
        function Cl=get_Cl(obj,S_ref)
            Cl=2*obj.aircraft_state.weight*9.81*obj.load_factor/(obj.aerodynamic_state.rho_air*norm(obj.V)^2*S_ref);
        end
        
        function obj=f_set_V_A(obj,V_A,varargin)
            [rho_air,a,~,P0,mu]=stdatmo(obj.h);
            if nargin==3
               V_A=correctairspeed(V_A,a,P0,varargin{1},'TAS');
            end
            Ma=V_A/a;
            obj.aerodynamic_state=class_aero_state(V_A,obj.aerodynamic_state.alpha,obj.aerodynamic_state.beta,Ma,rho_air,mu);
            obj.V_EAS=correctairspeed(V_A,a,P0,'TAS','EAS');
            obj.V=V_A;
        end
        
        
    end
    
end

