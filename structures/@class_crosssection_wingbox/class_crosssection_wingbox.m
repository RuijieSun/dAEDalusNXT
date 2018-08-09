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
%>        Klaus Seywald and Daniël de Vries
%>   dAEDalus is published under the terms of the GNU General Public
%>   License by the Free Software Foundation;
%>   either version 2 or any later version.
%>
%>   You should have received a copy of the GNU General Public
%>   License along with dAEDalus; see the file GNU GENERAL 
%>   PUBLIC LICENSE.TXT.  If not, write to the Free Software 
%>   Foundation, 59 Temple Place -Suite 330, Boston, MA
%>   02111-1307, USA.
%>

% ======================================================================
%> @brief class for a finite wingbox crosssection
%>
%> contains all GEOMETRICAL and STRUCTURAL information for a finite wingbox crosssection
% ======================================================================

classdef class_crosssection_wingbox < class_crosssection
    
    
    properties

        %> mean wingbox height                                  [m]
        h = 0;
        %> front spar height                                    [m]
        h_fs = 0;
        %> rear spar height                                     [m]
        h_rs = 0;
        %> wingbox width                                        [m]
        w = 0;
        %> chordlength of wingsection                           [m]
        c = 0;
        %> enclosed area                                        [m2]   
        A_enclosed = 0; 
        %> wetted crosssection in wing segment                  [m2]
        A_fuel = 0; 
        %> material area                                        [m2]
        A_material = 0;
        
        %% simple wingbox geometry
        %% Skin'
        %> thickness of upper skin                              [m]
        t_sk_up = 0;  
        %> loadcase index for the upper skin                    [m]
        t_sk_up_lc_idx = 0;
        %> thickness of lower skin                              [m]
        t_sk_lo = 0;    
        %> loadcase index for the lower skin                    [m]
        t_sk_lo_lc_idx = 0;

        %> thickness of upper stringers                         [m]
        t_st_up = 0;    
        %> thickness of lower stringers                         [m]
        t_st_lo = 0;     
        
        %> front spar thickness                                 [m]
        t_sp_fr = 0;
        %> rear spar thickness                                  [m]
        t_sp_re = 0;
        %> loadcase index for the front spar                    [m]
        t_sp_fr_lc_idx = 0;    
        %> loadcase index for the rear spar                     [m]
        t_sp_re_lc_idx = 0; 
        
        %> minimum allowable spar thickness
        t_min_sp = 0;
        %> minimum allowable skin thickness
        t_min_sk = 0;
        
        %% material properties at cross-section
        
        %> youngs modulus of the upper skin
        E_sk_up;
        %> shear modulus of the upper skin
        G_sk_up; 
        %> density of the upper skin
        rho_sk_up;
        %> youngs modulus of the lower skin
        E_sk_lo;
        %> shear modulus of the lower skin
        G_sk_lo;
        %> density of the lower skin
        rho_sk_lo;
        %> density of the spars
        rho_sp;
        
        %> tensile yield strength of the spars
        tensile_yield_sp;
        %> tensile yield strength of the lower skin
        tensile_yield_sk_l;
        %> tensile yield strength of the upper skin
        tensile_yield_sk_u;
        %> safety factor
        safety_factor = 1.5;
        fueling_factor = 1;
        
        %> stress upper skin
        sigma_sk_up = 0;
        %> stress lower skin
        sigma_sk_lo = 0;
        %> stress front spar
        sigma_sp_fr = 0;
        %> stress rear spar
        sigma_sp_re = 0;
        %> index of the wing_segment this crosssection belongs to
        segment_index = 0;
    end
    
    methods
        % =================================================================
        %> @brief Class constructor
        %>
        %> @return instance of the class_crosssection_wingbox
        % =================================================================
         function obj = class_crosssection_wingbox(h, w, c, t_sp_fr, t_sp_re, t_sk_up, t_sk_lo, fueling_factor, h_fs, h_rs)
                if nargin >= 3
                    obj.h = h;
                    obj.h_fs = h;
                    obj.h_rs = h;
                    obj.w = w;
                    obj.c = c;
                end
                
                if nargin >= 7
                    obj.t_sp_fr = t_sp_fr;
                    obj.t_sp_re = t_sp_re;
                    obj.t_sk_up = t_sk_up;
                    obj.t_sk_lo = t_sk_lo;
                end
                
                if nargin >= 8
                    obj.fueling_factor = fueling_factor;
                end
                
                if nargin == 10
                    obj.h_fs = h_fs;
                    obj.h_rs = h_rs;
                    obj.h = (h_fs + h_rs)/2;
                end
                                 
                obj = obj.calcArea();
         end
         
        % =================================================================
        %> @brief calculate area
        %>
        %> Calculated the enclosed and fuel area as a function of the geometry
        %>
        %> @param c wing crosssection chordlength
        %> @param h wingbox height
        %> @param w wingbox width
        %> @param h_fs front spar height
        %> @param h_rs rear spar height
        %>
        %> @return instance of the class_crosssection_wingbox
        % =================================================================
        function obj = calcArea(obj)  
            obj.A_enclosed = (obj.w - obj.t_sp_fr - obj.t_sp_re)/2 * (obj.h_fs + obj.h_rs - 2*(obj.t_sk_up + obj.t_sk_lo));
            obj.A_fuel = obj.fueling_factor * obj.A_enclosed;
            obj.A_material = obj.w/2*(obj.h_fs + obj.h_rs) - obj.A_enclosed;
        end
         
        % =================================================================
        %> @brief set geometry
        %>
        %> Sets the geometry parameters of the wingbox crosssection
        %>
        %> @param c wing crosssection chordlength
        %> @param h wingbox height
        %> @param w wingbox width
        %>
        %> @return instance of the class_crosssection_wingbox
        % =================================================================
        function obj=setGeometry(obj, varargin)
            if length(varargin) == 3
                obj.c = varargin{1};
                obj.h = varargin{2};
                obj.w = varargin{3};
            elseif length(varargin) == 4
                obj.c = varargin{1};
                obj.h_fs = varargin{2};
                obj.h_rs = varargin{3};
                obj.w = varargin{4};
                obj.h = 0.5*(obj.h_fs + obj.h_rs);
            end

            obj = obj.calcArea();
        end
  
        
        function obj=setMaterial(obj,E_sk_up,G_sk_up,rho_sk_up,E_sk_lo,G_sk_lo,rho_sk_lo,rho_sp,Ult_Tstrength_sp,Ult_Tstrength_sk_u,Ult_Tstrength_sk_l,t_min_sp,t_min_sk)
            obj.E_sk_up=E_sk_up;
            obj.G_sk_up=G_sk_up; 
            obj.rho_sk_up=rho_sk_up;
            %b) Lower Skin
            obj.E_sk_lo=E_sk_lo;
            obj.G_sk_lo=G_sk_lo;
            obj.rho_sk_lo=rho_sk_lo;
            %d) spars
            obj.rho_sp=rho_sp;
            obj.tensile_yield_sp=Ult_Tstrength_sp;
            obj.tensile_yield_sk_u=Ult_Tstrength_sk_u;
            obj.tensile_yield_sk_l=Ult_Tstrength_sk_l;
            
            obj.t_min_sk=t_min_sk;
            obj.t_min_sp=t_min_sp;
            obj.t_sk_lo=t_min_sk;
            obj.t_sk_up=t_min_sk;
            obj.t_sp_fr=t_min_sp;
            obj.t_sp_re=t_min_sp;
        end
        
        function [h,w,c, h_fs, h_rs]=get_dimensions(obj)
                h=obj.h;
                w=obj.w;
                c=obj.c;
                h_fs = obj.h_fs;
                h_rs = obj.h_rs;
        end
        
        % =================================================================
        %> @brief f_calc_skin_length
        %>
        %> Calculates the length of the skin segments assuming the wingbox
        %> is an isosceles parallelogram with the top and bottom skin
        %> panels having equal length.
        %>
        %> @return length of the skin panels, [m]
        % =================================================================
        function l_skin = f_calc_skin_length(obj)
            l_skin = sqrt((obj.h_fs - obj.h_rs)^2/4 + obj.w^2);
        end
                
        % =================================================================
        %> @brief f_calc_dm
        %>
        %> calculates the mass dm of a wingbox crosssection
        %>
        %> @return instance of the class_crosssection_wingbox
        % =================================================================
        function dm=f_calc_dm(obj)
            l_skin = obj.f_calc_skin_length();
            dm = obj.rho_sk_up*obj.t_sk_up*l_skin + obj.rho_sk_lo*obj.t_sk_lo*l_skin + obj.t_sp_fr*obj.h_fs*obj.rho_sp + obj.t_sp_re*obj.h_rs*obj.rho_sp;
        end
        
        % =================================================================
        %> @brief f_calc_shear_flows
        %>
        %> calculates the shear flows in the front/rear spars and the
        %> top/bottom skins.
        %>
        %> @return q_fs shear flow in front spar
        %> @return q_rs shear flow in rear spar
        %> @return q_ts shear flow in top skin
        %> @return q_bs shear flow in bottom skin
        % =================================================================
        function [q_fs, q_rs, q_ts, q_bs] = f_calc_shear_flows(obj, Mt, Qx, Qz)
            Sx_ts = obj.t_sk_up*obj.w*(obj.h_fs + obj.h_rs)/4; % First moment of area of top skin about the x-axis
            Sx_bs = obj.t_sk_lo*obj.w*(obj.h_fs + obj.h_rs)/4; % First moment of area of bottom skin about the x-axis
            
            Sz_fs = obj.t_sp_fr*obj.w*obj.h_fs/2; % First moment of area of front spar about the z-axis
            Sz_rs = obj.t_sp_re*obj.w*obj.h_rs/2; % First moment of area of rear spar about the z-axis
            
            [Ix, ~, Iz, ~, ~, ~, Ai] = calc_crosssection(obj); % Obtain the second moments of area about the x- and z-axis
            
            qQz_ts = abs(Qz)*Sx_ts/Ix; % Shear flow in the top skin due to the z shear force
            qQz_bs = abs(Qz)*Sx_bs/Ix; % Shear flow in the bottom skin due to the z shear force
            
            qQx_fs = abs(Qx)*Sz_fs/Iz; % Shear flow in the front spar due to the x shear force
            qQx_rs = abs(Qx)*Sz_rs/Iz; % Shear flow in the rear spar due to the x shear force
                       
            qMt = abs(Mt)/(2*Ai); % Shear flow due to the torsional moment
            
            q_fs = qQx_fs + qMt; % Total shear flow in the front spar
            q_rs = qQx_rs + qMt; % Total shear flow in the rear spar
            q_ts = qQz_ts + qMt; % Total shear flow in the top skin
            q_bs = qQz_bs + qMt; % Total shear flow in the bottom skin
        end
        
        % =================================================================
        %> @brief f_calc_stresses
        %>
        %> calculates the stresses in the front/rear spars and top/bottom
        %> skins based on the shear flows and the bending moment.
        %>
        %> @return instance of the class_crosssection_wingbox
        % =================================================================
        function obj = f_calc_stresses(obj, Mbx, Mby, Mt, Qx, Qz, loadcase_idx, overwrite)
            nt=0.8; %bending efficiency
            obj.sigma_sk_up = abs(Mbx)/(nt*obj.h*obj.w*obj.t_sk_up); % Shear stress in the top skin due to bending
            obj.sigma_sk_lo = abs(Mbx)/(nt*obj.h*obj.w*obj.t_sk_lo); % Shear stress in the bottom skin due to bending
            
            [q_fs, q_rs, q_ts, q_bs] = obj.f_calc_shear_flows(Mt, Qx, Qz); % Calculate the shear flows
                                                            
            tau_sp_fr = q_fs / obj.t_sp_fr; % Shear stress in the front spar
            tau_sp_re = q_rs / obj.t_sp_re; % Shear stress in the rear spar
            
            tau_sk_up = q_ts / obj.t_sk_up; % Shear stress in the top skin
            tau_sk_lo = q_bs / obj.t_sk_lo; % Shear stress in the bottom skin
            
            obj.sigma_sk_up = max([obj.sigma_sk_up, tau_sk_up/0.55]); 
            obj.sigma_sk_lo = max([obj.sigma_sk_lo, tau_sk_lo/0.55]);
            
            obj.sigma_sp_fr = tau_sp_fr / 0.55;
            obj.sigma_sp_re = tau_sp_re / 0.55;
        end
        
        % =================================================================
        %> @brief f_self_design_crosssection
        %>
        %> function to perform the crosssectional self-design ( this
        %> function contains the sizing rues for the crosssection)
        %>
        %> @return instance of the class_crosssection_wingbox
        % =================================================================
        function obj=f_self_design_crosssection(obj,Mbx,Mby,Mt,Qx,Qz,loadcase_idx,overwrite)
        
            nt=0.8; %bending efficiency
            
            sigma_lim=obj.tensile_yield_sk_l;
            
            t_b(1) = (abs(obj.safety_factor*Mbx)/(nt*obj.h*obj.w*sigma_lim));  % Equivalent upper skin thickness including stringer and skin thickness
            t_b(2) = (abs(obj.safety_factor*Mbx)/(nt*obj.h*obj.w*sigma_lim));  % Equivalent lower skin thickness including stringer and skin thickness

            t(1)= max(t_b(1),obj.t_min_sk);
            t(2)= max(t_b(2),obj.t_min_sk);
          
            % The absolute upper and lower wing skin thickness
            
            if ~overwrite
                if t(1)>obj.t_sk_up
                    obj.t_sk_up=t(1);
                    obj.t_sk_up_lc_idx=loadcase_idx;
                end
            else
                obj.t_sk_up=t(1);
                obj.t_sk_up_lc_idx=loadcase_idx;
            end
            
            if ~overwrite
                if t(2)>obj.t_sk_lo
                    obj.t_sk_lo=t(2);
                    obj.t_sk_lo_lc_idx=loadcase_idx;
                end
            else
                obj.t_sk_lo=t(2);
                obj.t_sk_lo_lc_idx=loadcase_idx;
            end
            
            % Sizing webs
            [q_fs, q_rs, q_ts, q_bs] = obj.f_calc_shear_flows(Mt, Qx, Qz); % calculate the shear flows                          

            shear_ult=sigma_lim*0.55;                       % limit shear stress
                                                            % load better value
                                                            % from material
                                                            % database later
                                                            
            t_fs = obj.safety_factor*q_fs/shear_ult;         % Front spar thickness at each node
            t_rs = obj.safety_factor*q_rs/shear_ult;         % Rear spar thickness at each node
            
            t_fs = max(t_fs,obj.t_min_sp);
            t_rs = max(t_rs,obj.t_min_sp);
        
            t_ts = obj.safety_factor*q_ts/shear_ult;  
            t_bs = obj.safety_factor*q_bs/shear_ult;
            
            t_ts = max(t_ts, t(1));
            t_bs = max(t_bs, t(2));
            
            if ~overwrite
                if t_ts>obj.t_sk_up
                    obj.t_sk_up=t_ts;
                    obj.t_sk_up_lc_idx=loadcase_idx;
                end
            else
                obj.t_sk_up=t_ts;
                obj.t_sk_up_lc_idx=loadcase_idx;
            end
            
            if ~overwrite
                if t_bs>obj.t_sk_lo
                    obj.t_sk_lo=t_bs;
                    obj.t_sk_lo_lc_idx=loadcase_idx;
                end
            else
                obj.t_sk_lo=t_bs;
                obj.t_sk_lo_lc_idx=loadcase_idx;
            end
               
            if ~overwrite
                if t_fs>obj.t_sp_fr
                    obj.t_sp_fr=t_fs;
                    obj.t_sp_fr_lc_idx=loadcase_idx;
                end
            else
                obj.t_sp_fr=t_fs;
                obj.t_sp_fr_lc_idx=loadcase_idx;
            end
            
            if ~overwrite
                if t_rs>obj.t_sp_re
                    obj.t_sp_re=t_rs;
                    obj.t_sp_re_lc_idx=loadcase_idx;
                end
            else
                obj.t_sp_re=t_rs;
                obj.t_sp_re_lc_idx=loadcase_idx;
            end
        end
        
        % =================================================================
        %> @brief calc_crosssection
        %>
        %> calculates the crosssectional parameterx Ix, Iy, Iy, A,
        %> Aenclosed
        %>
        %> @return Ix
        %> @return Iy
        %> @return Iz
        %> @return A
        %> @return Aenclosed
        % =================================================================
        function [Ix,Iy,Iz,J,Amaterial,Afuel,Aenclosed]=calc_crosssection(obj)          
            dh = (obj.t_sk_up + obj.t_sk_lo)/2; % The height to add/subtract to get the outer/inner heights, [m]
            dw = (obj.t_sp_fr + obj.t_sp_re)/2; % The width to add/subtract to get the outer/inner widths, [m]
            
            % Calculate the contributions of the outer section to Ix and Iz
            h_fs_o = obj.h_fs + dh; % Outer height of front spar, [m]
            h_rs_o = obj.h_rs + dh; % Outer height of the rear spar, [m]
            w_o = obj.w + dw; % Outer width, [m]
            Ix_o = w_o*(h_fs_o + h_rs_o)*(h_fs_o^2 + h_rs_o^2)/48; % Ix of the outer block shape, [m4]
            Iz_o = w_o^3*(h_fs_o^2 + 4*h_fs_o*h_rs_o + h_rs_o^2)/(36*(h_fs_o + h_rs_o)); % Iz of the outer block shape, [m4]
            
            % Calculate the contributions of the inner section to Ix and Iz
            h_fs_i = obj.h_fs - dh; % Inner height of the front spar, [m]
            h_rs_i = obj.h_rs - dh; % Inner height of the rear spar, [m]
            w_i = obj.w - dw; % Inner width, [m]
            Ix_i = w_i*(h_fs_i + h_rs_i)*(h_fs_i^2 + h_rs_i^2)/48; % Ix of the inner block shape, [m4]
            Iz_i = w_i^3*(h_fs_i^2 + 4*h_fs_i*h_rs_i + h_rs_i^2)/(36*(h_fs_i + h_rs_i)); % Iz os the inner block shape, [m4]
            
            % Calculate the total Ix and Iz using the parallel axis
            % theorem. We neglect the A*d^2 term, under the assumption that
            % the centroids of the inner and outer shape coincide.
            Ix = Ix_o - Ix_i; % Total Ix, [m4]
            Iz = Iz_o - Iz_i; % Total Iz, [m4]
            
            % The polar moment of inertia is just the sum of Ix and Iz.
            J = Ix + Iz; % Polar moment of inertia, [m4]
            
            % Calculate the areas
            obj = obj.calcArea();
            Amaterial = obj.A_material; % Area of the material, [m2]
            Afuel = obj.A_fuel; % Area of the fuel, [m2]
            Aenclosed = obj.A_enclosed; % True enclosed area of the section, [m2]
            
            % Calculate the torsion constant. Here this value is stored in
            % Iy. We use the 2nd formula of Bredt which gives the torsion
            % constant for a closed, thin walled cross section with exectly
            % one hole. This formula expresses the total torsion constant
            % as follows:
            %     It = 4*A_inner^2/(Sum d/t)
            %     where:
            %        - It is the torsion constant, [m4]
            %        - A_inner is the area truly enclosed, [m2]
            %        - d is the length of each segment, [m]
            %        - t is the thickness of each segment, [m]
            l_skin = obj.f_calc_skin_length(); % Length of the top and bottom skin pannels, [m]
            Iy = 4*obj.A_enclosed^2/(l_skin/obj.t_sk_up + l_skin/obj.t_sk_lo + obj.h_fs/obj.t_sp_fr + obj.h_rs/obj.t_sp_re); % Torsion constant, [m4]
        end
    end   
end

