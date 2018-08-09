%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%
% This class has the same purpose as the class_crosssection_wingbox class.
% It is a re-adaptation of the aforementioned class, to be used for
% anisotropic materials.

classdef class_crosssection_wingbox_anisotropic < class_crosssection
    
    
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
        
        %> youngs modulus of the ply material
        E1_ply;
        E2_ply;
        %> shear modulus of the ply material
        G12_ply; 
        %> density of the the ply material
        rho_ply;
        
        %> tensile yield strength of the ply material
        ult_sigma_tens_ply_1; %longitudinal
        ult_sigma_tens_ply_2; %transverse
        %> compresive yield strength of the ply material
        ult_sigma_comp_ply_1; %longitudinal
        ult_sigma_comp_ply_2;%transverse
        %> tensile yield strength of the ply material
        ult_tau_shear_ply_12;
        %> poisson ratios of the ply material
        mu_xy; %longitudinal
        mu_yx; %transverse
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
        
        %% variables related to cross section discretization (cross sectional modeler)
        
        % reference point used in local cross section discretization.
        % Placed in the center of the wingbox
        ref_point_yz;
        
        % number of nodes for each element in the cross section
        fs_nodes_nr;      % front spar
        rs_nodes_nr;      % rear spar
        sk_up_nodes_nr;   % upper skin
        sk_lo_nodes_nr;   % lower skin
        
        % yz coordinates of each node in discretized cross section
        % elements
        fs_nodes_yz      % front spar
        rs_nodes_yz      % rear spar
        sk_up_nodes_yz   % upper skin
        sk_lo_nodes_yz   % lower skin
        
        % Local Euler Bernoulli stiffness matrix of cross section.
        Se; %CURRENTLY ONLY Se IS USED IN THE MODEL
        % Local Euler Bernoulli compliance matrix of cross section
        Ce;
        % Local shear compliance matrix of cross section
        Cs;
        % Local coupling of euler-shear forces
        Ces;
        
        % Matrix used to recover shell strains from the nodal displacements
        % of the Euler Bernoulli beam
        Gamma_euler; %CURRENTLY ONLY Gamma_euler IS USED IN THE MODEL
        % Matrix used to recover shell strains from the shear displacements
        % of the beam
        Gamma_shear;
        % Coupling of Gamma_euler and Gamma_shear
        Gamma;
        
        % ABD stiffness matrices of each cross sectional element
%         fs_ABD;
%         rs_ABD;
%         sk_up_ABD;
%         sk_lo_ABD;
%         
%         % ABD compliance matrices of each cross sectional element
%         fs_abd;
%         rs_abd;
%         sk_up_abd;
%         sk_lo_abd;
        
        % Array containing the IDs of all shell elements of the cross
        % section (ID range from 1 to 4, representing fs,sk_up,rs,sk_lo)
        shell_id_arr;
        
        % Arrays for the strains in skins and spars. These are the strains
        % recovered from the nodal displacements through the cross
        % sectional modeler.
        strain_fs;
        strain_sk_up;
        strain_rs;
        strain_sk_lo;
        
        % Arrays for the stresses in skins and spars. These are obtained
        % from the strains mentioned above.
        stress_fs;
        stress_sk_up;
        stress_rs;
        stress_sk_lo;
        
        % Arrays with the Von Mises stresses for skins and spars, derived
        % from stresses above.
        stress_von_mis_fs;
        stress_von_mis_sk_up;
        stress_von_mis_rs;
        stress_von_mis_sk_lo;
        
        % Instances of the laminate class, representing the lamiante of
        % each respective skin or spar
        laminate_sk_up;
        laminate_sk_lo;
        laminate_fs;
        laminate_rs;
                


    end
    
    methods
        % =================================================================
        %> @brief Class constructor
        %>
        %> @return instance of the class_crosssection_wingbox
        % =================================================================
         function obj = class_crosssection_wingbox_anisotropic(h, w, c, t_sp_fr, t_sp_re, t_sk_up, t_sk_lo, fueling_factor, h_fs, h_rs)
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
         
         % cross sectional modeler which turns 2D composite cross section into
         % 1D beam (stored in separate file)
         obj= cross_sectional_modeler(obj);
         
         
         % Function which updates the cross section with the optimizer
         % output.
         
         % Input is 8 variables: 4 laminante thicknesses (in METERS), and 4 offset
         % angles (in DEGREES).
         
         % The order of the input is the following: front spar, top skin,
         % rear spar, bottom skin.
         
         function obj = crosssection_update(obj,t_fs,t_sk_up,t_rs,t_sk_lo,offset_fs,offset_sk_up,offset_rs,offset_sk_lo)
            
%              update_flag = 0;
             
             % checks if the properties stored in the laminate are
             % different from the ones inputted in this function
%              if t_fs~=obj.laminate_fs.t_lam || offset_fs~=obj.laminate_fs.offset_angle 
                 
             % Set laminate thickness to the new thickness
             obj.laminate_fs.t_lam = t_fs;
             %Set laminate offset angle to the new angle
             obj.laminate_fs.offset_angle = offset_fs;
             % Updates the laminate
             obj.laminate_fs = obj.laminate_fs.update_laminate();

             % Updates thickness in the cross section
             obj.t_sp_fr = t_fs;
                 
                 % Sets update_flag, so that cross sectional modeller can
                 % run
%                  update_flag = 1;
%              end
             
             % Check out comments above
%              if t_sk_up~=obj.laminate_sk_up.t_lam || offset_sk_up~=obj.laminate_sk_up.offset_angle 
             obj.laminate_sk_up.t_lam = t_sk_up;
             obj.laminate_sk_up.offset_angle = offset_sk_up;
             obj.laminate_sk_up = obj.laminate_sk_up.update_laminate();
             obj.t_sk_up = t_sk_up;
%                  update_flag = 1;
%              end
             
             % Check out comments above
%              if t_rs~=obj.laminate_rs.t_lam || offset_rs~=obj.laminate_rs.offset_angle 
             obj.laminate_rs.t_lam = t_rs;
             obj.laminate_rs.offset_angle = offset_rs;
             obj.laminate_rs = obj.laminate_rs.update_laminate();
             obj.t_sp_re = t_rs;
%                  update_flag = 1;
%              end
             
             % Check out comments above
%              if t_sk_lo~=obj.laminate_sk_lo.t_lam || offset_sk_lo~=obj.laminate_sk_lo.offset_angle 
             obj.laminate_sk_lo.t_lam = t_sk_lo;
             obj.laminate_sk_lo.offset_angle = offset_sk_lo;
             obj.laminate_sk_lo = obj.laminate_sk_lo.update_laminate();
             obj.t_sk_lo = t_sk_lo;
%                  update_flag = 1;
%              end



             obj = obj.calcArea();
             
             % If any of the lamiantes were updated, the cross sectional
             % modeller runs 
%              if update_flag==1
             obj = obj.cross_sectional_modeler();
%              end
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
  
        
        % The setMaterial function has been modified to use anisotropic
        % materials. It should be noted that thanks to the cross sectional
        % modeler this is not actually used.
        % It still exists in order to avoid errors when code written before
        % the anisotropic model was created is run (as old pieces of code
        % would try to set the Material).
        % The values below are not actually used, and the function should
        % be removed once the old code is revised, and calls to this
        % functions are removed when dealing with anisotropic cross
        % sections.
        function obj=setMaterial(obj, structure, layup_settings, varargin)
            
            %> youngs modulus of the ply material
            obj.E1_ply = structure.materials.E1;
            obj.E2_ply = structure.materials.E2;
            
            %> shear modulus of the ply material
            obj.G12_ply = structure.materials.G;
            
            %> density of the the ply material
            obj.rho_ply = structure.materials.rho;

            %> tensile yield strength of the ply material
            obj.ult_sigma_tens_ply_1 = structure.materials.ult_sigma_tens_1;
            obj.ult_sigma_tens_ply_2 = structure.materials.ult_sigma_tens_2; 
            
            %> compresive yield strength of the ply material
            obj.ult_sigma_comp_ply_1 = structure.materials.ult_sigma_comp_1;
            obj.ult_sigma_comp_ply_2 = structure.materials.ult_sigma_comp_1;
            
            %> tensile yield strength of the ply material
            obj.ult_tau_shear_ply_12 = structure.materials.ult_tau_shear_12;
            
            %> poisson ratios of the ply material
            obj.mu_xy = structure.materials.poiss_ratio;
            obj.mu_yx = structure.materials.poiss_ratio_minor;
            
            % Material invariants used to solve for feasible design region
            design_envelope_u_arr = [structure.materials.u1, structure.materials.u2, structure.materials.u3, structure.materials.u4, structure.materials.u5, structure.materials.u6];
            
            %> minimum wing and spar thicknesses. Used for first guess of
            %  thickness and laminate creation
            obj.t_min_sk = layup_settings.wingTminSK;
            obj.t_min_sp = layup_settings.wingTminSP;
            
            if nargin == 4
                obj.t_min_sk = varargin{1};
                obj.t_min_sp = varargin{1};
            end
            
            obj.t_sk_up = obj.t_min_sk;
            obj.t_sk_lo = obj.t_min_sk;
            
            obj.t_sp_fr = obj.t_min_sp;
            obj.t_sp_re = obj.t_min_sp;
            
            
            % runs layup generator (complex version, used only for
            % initialization) once, for the skin and spar
            
            [final_layup_sk, t_array_final_sk] = ...
                class_laminate.layup_generator(layup_settings.skins_layup_angles, layup_settings.skins_layup_fractions, obj.t_min_sk, structure.materials.t_ply);
            
            [final_layup_sp, t_array_final_sp] = ...
                class_laminate.layup_generator(layup_settings.spars_layup_angles, layup_settings.spars_layup_fractions, obj.t_min_sp, structure.materials.t_ply);
            
            % IMPORTANT: THE FINAL LAYUPS BELOW ASSUME THAT THE INPUT LAYUP
            % ORIENTATION FOLLOWS THE BEAM ELEMENT COORDINATE SYSTEM. THE
            % INPUT BELOW IS TRANSLATED IN ORDER TO MAKE SENSE IN THE CROSS
            % SECTIONAL MODELLER COORDINATE SYSTEM.
            
            % The cross sectional modeller coordinate system follows the
            % coordinate system of its shell element, and it works in the
            % same way as NASTRAN's CQUAD4.
            
            % Normal of upper skin is pointing opposite to z-axis of beam element
            % Normal of lower skin is pointing in same direction as z-axis of beam element
            % Normal of front spar is pointing in same direction as x-axis of beam element
            % Normal of lower skin is pointing opposite to z-axis of beam element
            final_layup_sk_up = final_layup_sk;
            final_layup_sk_lo = -final_layup_sk; %lower skin laminate is defined in the beam coordinate system 
            
            final_layup_sp_fr = final_layup_sp;
            final_layup_sp_re = final_layup_sp;
            
            % creating the laminate class instances belonging to each skin
            % and spar
            
            % list with ply properties needed for laminate. t_ply is not
            % actually used in the class constructor since we are feeding
            % it an array of thicknesses for each ply (t_array_final_sk and
            % t_array_final_sp.
            ply_properties = [obj.E1_ply, obj.E2_ply, obj.mu_xy, obj.mu_yx, obj.G12_ply];
            
            offset_angle_first_run = 0;
            
            obj.laminate_sk_up = class_laminate(final_layup_sk_up, ply_properties,t_array_final_sk, design_envelope_u_arr, offset_angle_first_run);
            obj.laminate_sk_lo = class_laminate(final_layup_sk_lo, ply_properties,t_array_final_sk, design_envelope_u_arr, offset_angle_first_run);
            obj.laminate_fs = class_laminate(final_layup_sp_fr, ply_properties, t_array_final_sp, design_envelope_u_arr, offset_angle_first_run);
            obj.laminate_rs = class_laminate(final_layup_sp_re, ply_properties, t_array_final_sp, design_envelope_u_arr, offset_angle_first_run);
            
            % Assigns ABD stiffness matrices to skins/spars
%             obj.fs_ABD = obj.laminate_fs.ABD_stiff;
%             obj.rs_ABD = obj.laminate_rs.ABD_stiff;
%             obj.sk_up_ABD = obj.laminate_sk_up.ABD_stiff;
%             obj.sk_lo_ABD = obj.laminate_sk_lo.ABD_stiff;
% 
%             % Assigns abd compliance matrices to skins/spars
%             obj.fs_abd = obj.laminate_fs.abd_compl;
%             obj.rs_abd = obj.laminate_rs.abd_compl;
%             obj.sk_up_abd = obj.laminate_sk_up.abd_compl;
%             obj.sk_lo_abd = obj.laminate_sk_lo.abd_compl;
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
            dm = obj.rho_ply*obj.t_sk_up*l_skin + obj.rho_ply*obj.t_sk_lo*l_skin + obj.t_sp_fr*obj.h_fs*obj.rho_ply + obj.t_sp_re*obj.h_rs*obj.rho_ply;
        end
           
        
        % ALL THE INTERIAS CALCULATED WITHIN THIS FUNCTION ARE NO LONGER
        % NEEDED, AS THE CROSS SECTIONAL MODELER GENERATES AN EULER
        % STIFFNESS MATRIX BY ITSELF.
        % =================================================================
        %> @brief calc_crosssection
        %>
        %> calculates the crosssectional parameterx Ix, Iy, Iy, A,
        %> Aenclosed
        %>
        %> @return A
        %> @return Aenclosed
        % =================================================================
        function [Amaterial,Afuel,Aenclosed]=calc_crosssection(obj)       

            
            % Calculate the areas
            obj = obj.calcArea();
            Amaterial = obj.A_material; % Area of the material, [m2]
            Afuel = obj.A_fuel; % Area of the fuel, [m2]
            Aenclosed = obj.A_enclosed; % True enclosed area of the section, [m2]

        end
        
        % Discretizes the cross section into 1D elements, in order to run
        % Abdalla's cross sectional modeler. This allows to obtain a Euler
        % Bernoulli stiffness matrix, thus representing the complex 3D
        % composite beam with a simplified 1D beam.
        
        function obj = discretize_cross_section(obj,ds_max,nr_nodes)
            
            % reference point, set in the middle of the cross ection
            obj.ref_point_yz = [obj.h/2, obj.w/2]; %Y and Z coord of ref point.
            
            
            % defining amount of nodes present in each cross section
            % element
            % if ds_max is different from zero, then the discretization
            % makes sure that the amount of nodes is enough to make all
            % shell elements shorter than ds_max
            if ds_max ~= 0
                obj.fs_nodes_nr = ceil(obj.h_fs/ds_max);
                obj.rs_nodes_nr = ceil(obj.h_rs/ds_max);
                obj.sk_up_nodes_nr = ceil(obj.w/ds_max);
                obj.sk_lo_nodes_nr = ceil(obj.w/ds_max);
            
            % if nr_nodes is different from zero, then the discretization
            % assigns nr_nodes nodes to each cross section skin/spar
            elseif nr_nodes ~= 0
                obj.fs_nodes_nr = nr_nodes;
                obj.rs_nodes_nr = nr_nodes;
                obj.sk_up_nodes_nr = nr_nodes;
                obj.sk_lo_nodes_nr = nr_nodes;
            end
            
            % defining yz coordinates of each node within each cross
            % section element. Order is counter-clockwise, starting from bottom of
            % front spar
            
            
            % Code commented out below refered to wrong version of
            % coordinates, which is also what has been used until Jan 9th
            % 2018. The fs and rs have been switched around (the front spar
            % within Daedalus has negative coordinates, meaning that they
            % should be positive within the cross sectional modeler, as the
            % flap axis is opposite between the two coord systems).
%             obj.fs_nodes_yz = [-obj.w/2 * ones(obj.fs_nodes_nr,1),linspace(-obj.h_fs/2,+obj.h_fs/2,obj.fs_nodes_nr)'];
%             obj.rs_nodes_yz = [obj.w/2 * ones(obj.rs_nodes_nr,1),linspace(+obj.h_rs/2,-obj.h_rs/2,obj.rs_nodes_nr)'];
%             obj.sk_up_nodes_yz = [linspace(-obj.w/2,obj.w/2,obj.sk_up_nodes_nr)', linspace(+obj.h_fs/2, +obj.h_rs/2, obj.sk_up_nodes_nr)'];
%             obj.sk_lo_nodes_yz = [linspace(obj.w/2,-obj.w/2,obj.sk_lo_nodes_nr)',linspace(-obj.h_rs/2, -obj.h_fs/2, obj.sk_lo_nodes_nr)'];


            % The new code below uses the right coordinates. Hopefully no
            % additional changes are required for the code to work.
            % The order in which the cross section is discretized is now
            % anti-clockwise, starting from fs, then sk_up, then rs, then sk_lo.
            % Before the same was done, but in an clockwise direction.
            
            obj.fs_nodes_yz = [obj.w/2 * ones(obj.fs_nodes_nr,1),linspace(-obj.h_fs/2,+obj.h_fs/2,obj.fs_nodes_nr)'];
            obj.rs_nodes_yz = [-obj.w/2 * ones(obj.rs_nodes_nr,1),linspace(+obj.h_rs/2,-obj.h_rs/2,obj.rs_nodes_nr)'];
            obj.sk_up_nodes_yz = [linspace(obj.w/2,-obj.w/2,obj.sk_up_nodes_nr)', linspace(+obj.h_fs/2, +obj.h_rs/2, obj.sk_up_nodes_nr)'];
            obj.sk_lo_nodes_yz = [linspace(-obj.w/2,obj.w/2,obj.sk_lo_nodes_nr)',linspace(-obj.h_rs/2, -obj.h_fs/2, obj.sk_lo_nodes_nr)'];
            
        end
        
%         % function used to assign to each spar/skin the ABD matrix of their
%         % respective laminates.
%         function obj = ABD_matrix(obj,laminate_fs,laminate_rs,laminate_sk_up,laminate_sk_lo)
% 
%             % Assigns ABD stiffness matrices to skins/spars
%             obj.fs_ABD = laminate_fs.ABD_stiff;
%             obj.rs_ABD = laminate_rs.ABD_stiff;
%             obj.sk_up_ABD = laminate_sk_up.ABD_stiff;
%             obj.sk_lo_ABD = laminate_sk_lo.ABD_stiff;
% 
%             % Assigns abd compliance matrices to skins/spars
%             obj.fs_abd = laminate_fs.abd_compl;
%             obj.rs_abd = laminate_rs.abd_compl;
%             obj.sk_up_abd = laminate_sk_up.abd_compl;
%             obj.sk_lo_abd = laminate_sk_lo.abd_compl;
%             
%         end
        
        % stresses are calculated within the anisotropic beam element using
        % function f_calc_stress_strain_crossmod.
        % The function below is defined only to avoid error due to abstract
        % methods present in class_crosssection not being present here
        function f_calc_stresses(Mbx, Mby, Mt, Qx, Qz, loadcase_idx, overwrite)
        end
        
        % the self design cross section function is not needed anymore, as
        % failure of the beam elements is calculated within the laminate
        % object, using safety factor function.
        % The function below is defined only to avoid error due to abstract
        % methods present in class_crosssection not being present here
        function f_self_design_crosssection(Mbx, Mby, Mt, Qx, Qz, loadcase_idx, overwrite)
        end
        
    end 
end

