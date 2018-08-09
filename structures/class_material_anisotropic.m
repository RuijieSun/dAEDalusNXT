%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%
% This class defines anisotropic materials.

% Given the nature of anisotropic materials, all the properties below are
% ply properties.

% This class has the same purposes of class_material

classdef class_material_anisotropic
    
    properties
        %> Ply young's modulus
        E1;     % in fiber direction
        E2;     % transverse to fiber direction
        
        %> shear modulus
        G;
        
        %> material density
        rho;
        
        %> allowable yield stress
        sigma_allowable_1;
        sigma_allowable_2;
        
        %> allowable shear stress
        tau_allowable;
        
        %> ultimate tensile and compressive strengths
        ult_sigma_tens_1;   % fiber direction
        ult_sigma_tens_2;   % transverse to fiber direction
        
        ult_sigma_comp_1;   % fiber direction
        ult_sigma_comp_2;   % transverse to fiber direction
        
        %> ultimate in-plane shear stress
        ult_tau_shear_12;        
        
        %> major poisson's ratio
        poiss_ratio;          % reaction of transverse direction to change in fiber direction
        poiss_ratio_minor;    % reaction of fiber direction to change in transverse direction
        
        t_ply;  % thickness of a single ply
        
        % Material invariants used to sovlve for feasible design region
        u1;u2;u3;u4;u5;u6;
    end
    
    methods
        
        % =================================================================
        %> @brief Class constructor
        %>
        %> Initializes the coupling condition
        %>
        %> @param material_string name of the material
        %>
        %> @return instance of class_material
        % =================================================================
       
        function obj =class_material_anisotropic(material_string)

%             if strcmp('CFRP_0/90_QI',material_string)
%                 % Mech. properties taken from
%                 % http://www.performance-composites.com/carbonfibre/mechanicalproperties_2.asp
% 
%                 obj.E1 = 7.0E10;
%                 obj.E2 = 7.0E10;
%                 obj.G = 5E10;
%                 obj.rho = 1600;
% 
%                 obj.ult_sigma_tens_1 = 600.0E6;
%                 obj.ult_sigma_tens_2 = 600.0E6;
% 
%                 obj.ult_sigma_comp_1 = 570.0E6;
%                 obj.ult_sigma_comp_2 = 570.0E6;
% 
%                 obj.ult_tau_shear_12 = 90.0E6;  
%             end
% 
%             if strcmp('CFRP_UD',material_string)
%                 % Mech. properties taken from
%                 % http://www.performance-composites.com/carbonfibre/mechanicalproperties_2.asp
% 
%                 obj.E1 = 13.5E10;
%                 obj.E2 = 1.0E10;
%                 obj.G = 5E10;
%                 obj.rho = 1600;
% 
%                 obj.ult_sigma_tens_1 = 1500.0E6;
%                 obj.ult_sigma_tens_2 = 50.0E6;
% 
%                 obj.ult_sigma_comp_1 = 1200.0E6;
%                 obj.ult_sigma_comp_2 = 250.0E6;
% 
%                 obj.ult_tau_shear_12 = 70.0E6;  
%             end
            
            if strcmp('CFRP_standard',material_string)
                % Used in XML files, taken from Iksselmuiden, Abdalla and
                % Gurdal, DOI: 10.2514/1.35565, table B1

                obj.E1 = 142.0E9;
                obj.E2 = 10.3E9;
                obj.G = 7.2E9;
                obj.rho = 1320;

                obj.poiss_ratio = 0.27;
                obj.poiss_ratio_minor = 0.0195845; %poiss_ratio/(E1/E2)

                obj.t_ply = 0.000127;

                obj.ult_sigma_tens_1 = 2280.0E6;
                obj.ult_sigma_tens_2 = 57.0E6;

                obj.ult_sigma_comp_1 = 1440.0E6;
                obj.ult_sigma_comp_2 = 228.0E6;

                obj.ult_tau_shear_12 = 71.0E6;  
            end
            
            if strcmp('CFRP_Werter',material_string)
                % Used by werter for optimizations with predefined laminate
                % including knock down factors as specified in his phd
                % thesis page 210
                obj.E1 = 147.0E9;
                obj.E2 = 10.3E9;
                obj.G = 7.0E9;
                obj.rho = 1600;

                obj.poiss_ratio = 0.27;
                obj.poiss_ratio_minor = 0.018918367; %poiss_ratio/(E1/E2)

                obj.t_ply = 0.000127;

                obj.ult_sigma_tens_1 = 948.5E6;
                obj.ult_sigma_tens_2 = 23.7E6;

                obj.ult_sigma_comp_1 = 717.6E6;
                obj.ult_sigma_comp_2 = 94.8E6;

                obj.ult_tau_shear_12 = 31.6E6;  
            end
            
            if strcmp('Werter_validation',material_string)
                % Mech. properties taken from Table 3.2 of Werter's PhD thesis

                obj.E1 = 141.96E9;
                obj.E2 = 9.79E9;
                obj.G = 6E9;
                obj.rho = 1600;

                obj.poiss_ratio = 0.42;
                obj.poiss_ratio_minor = 0.02896;

                obj.t_ply = 0.000127;

                obj.ult_sigma_tens_1 = 1500.0E6;
                obj.ult_sigma_tens_2 = 50.0E6;

                obj.ult_sigma_comp_1 = 1200.0E6;
                obj.ult_sigma_comp_2 = 250.0E6;

                obj.ult_tau_shear_12 = 70.0E6;  
            end

            if strcmp('aniso_aluminum',material_string)
                % Mech. properties taken from Table 3.2 of Werter's PhD thesis

                obj.E1 = 6.89E10;
                obj.E2 = 6.89E10;
                obj.G = 2.6E10;
                obj.rho = 2800;

                obj.poiss_ratio = 0.42;
                obj.poiss_ratio_minor = 0.42;

                obj.t_ply = 0.001;

                obj.ult_sigma_tens_1 = 250E6;
                obj.ult_sigma_tens_2 = 250E6;

                obj.ult_sigma_comp_1 = 250E6;
                obj.ult_sigma_comp_2 = 250E6;

                obj.ult_tau_shear_12 = 250E6;  
            end
            
            % calculates design envelope for material specified in this
            % instance
            obj = obj.tsai_wu_fail_envelope();
        end
        
        % This function calculates the Tsai Wu failure envelope for a
        % specified ply material, according to the theory illustrated in
        % Ijsselmuiden, Abdalla  and Gurdal, 2007, DOI: 10.2514/1.35565
        function obj = tsai_wu_fail_envelope(obj)
            
            % Calculating stress based coefficients for the Tsai Wu failure criterion
            F11 = 1/(obj.ult_sigma_tens_1 * obj.ult_sigma_comp_1);
            F22 = 1/(obj.ult_sigma_tens_2 * obj.ult_sigma_comp_2);
            F1 = 1/(obj.ult_sigma_tens_1) - 1/(obj.ult_sigma_comp_1);
            F2 = 1/(obj.ult_sigma_tens_2) - 1/(obj.ult_sigma_comp_2);
            F12 = -1/sqrt(obj.ult_sigma_tens_1 * obj.ult_sigma_comp_1 * obj.ult_sigma_tens_2 * obj.ult_sigma_comp_2);
            F66 = 1/(obj.ult_tau_shear_12^2);
            
            
            % Calculating local stiffness matrix of a single ply, in the
            % ply coordinate system.
            mu_xy = obj.poiss_ratio;
            mu_yx = obj.poiss_ratio_minor;
            
            den = 1. - mu_xy*mu_yx;
            
            Q_mat = [obj.E1 / (den),    mu_xy*obj.E2/den,   0;
                     mu_yx*obj.E1/den,  obj.E2 / (den),     0;
                     0,                    0,              obj.G];
            
            % Assigning relevant elements of Q matrix to variables
            Q11 = Q_mat(1,1);
            Q22 = Q_mat(2,2);
            Q12 = Q_mat(1,2);
            Q66 = Q_mat(3,3);
            
            
            % Calculating strain based coefficients for the Tsai Wu failure criterion
            G11 = Q11^2 * F11 + Q12^2 * F22 + 2*F12 * Q11 * Q12;
            G22 = Q12^2 * F11 + Q22^2 * F22 + 2*F12 * Q12 * Q22;
            G1 = Q11 * F1 + Q12 * F2;
            G2 = Q12 * F1 + Q22 * F2;
            G12 = Q11 * Q12 * F11 + Q12 * Q22 * F22 + F12 * Q12^2 + F12 * Q11 * Q22;
            G66 = 4 * Q66^2 * F66;
            
            % Material invariants used to sovlve for feasible design region
            obj.u1 = G11 + G22 - 2*G12;
            obj.u2 = (G1 + G2) / 2;
            obj.u3 = (G11 + G22 + 2*G12) / 4;
            obj.u4 = G1 - G2;
            obj.u5 = G11 - G22;
            obj.u6 = G66;
            
        end
    end
end

