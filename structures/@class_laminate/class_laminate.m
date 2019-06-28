%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%
% This class represents a specific laminate.

% Input: Laminate Layup and Ply properties

% Output: ABD stiffness matrix and abd compliance matrix.

classdef class_laminate
    properties
        
        %% Basic ply properties
        E_x_arr;     % Array with E_x of each ply
        E_y_arr;     % Array with E_y of each ply
        mu_xy_arr;   % Array with my_xy of each ply
        mu_yx_arr;   % Array with my_yx of each ply
        G_xy_arr;    % Array with G_xy of each ply
        thick_arr;   % Array with thickness of each ply
        phi_rad;     % Array with angle of each ply (wrt to laminate coord sys)
        z_low;       % Array with lower z coordinate of each ply (wrt to laminate mid-plane)
        z_up;        % Array with upper z coordinate of each ply (wrt to laminate mid-plane)

        %% Ply matrices
        C_loc;           % Ply stiffness matrix (in ply local coord sys)
        S_loc;           % Ply compliance matrix (in ply local coord sys)
        S_loc_inverted;  % Inverted ply compliance matrix (in ply local coord sys)
        C_phi;           % Ply stiffness matrix (in laminate coord sys)
        S_phi_inverted;  % Inverted ply compliance matrix (in laminate coord sys)
        S_phi;           % Ply compliance matrix (in local coord sys)
        
        %% Basic laminate properties
        ID;              % Array with ID of each ply (starting from the bottom, with 0 being the value of the top ply)
        t_lam;           % Thickness of the entire laminate
        ABD_stiff;       % ABD stiffness matrix of the laminate
        abd_compl;       % abd compliance matrix of the laminate
        angle_lst_deg;   % Array with laminate layup angles in degrees
        angle_lst_rad;   % Array with laminate layup angles in radians
        angle_lst_deg_original; % Array with laminate layup angles (without offset angle) in degrees
        
        
        
        %% Laminate design envelope properties
        % Material invariants used to solve for feasible design region
        u1;u2;u3;u4;u5;u6;
        % Material invariants used to solve for feasible design region as
        % defied by khani
        u1k;u2k;u3k;u4k;u5k;u6k;
        
        % Failure indices of each skin or spar, indicating how close (or
        % past) failure each skin or spar are. It includes the safety
        % factor.
        
        % The array contain the most critical failure index of
        % each shell element within the skin/spar.
        fail_idx_arr;
        
        % The properties below the most critical failure index of the whole
        % skin/spar (so the shell element which fails the worst)
        fail_idx_crit;
        
        % Angle at which the whole layup is currently rotated. THIS IS IN
        % LAMINATE COORDINATE SYSTEM
        offset_angle;
        
        % Failure envelope formulation (1 for std IJsselmuiden; 2 for
        % Khani)
        failureEnvelopeType=2;
    end
    
    methods (Static)
        
        % This function generates the layup which is used to create an
        % instance of the class_laminate.
        
        % This function is required in order to turn laminate thickness
        % into a countinuous variable.
        
        % INPUTS:
        %   * input_angles:    array with a list of angles to be included in
        %                      the layup (eg. [0,30,-30,45,-45,90]. Each angle
        %                      should only be mentioned ONCE.
        %                      WARNING: the finer the fractions (like
        %                      45.67% rather than 45%) the more plies will
        %                      be required in order to correctly represent
        %                      the fraction distribution. It is advised to
        %                      stick to integer fractions.
        %
        %   * input_fractions: these are the percentages of each angle wrt
        %                      the whole laminate. Each percentage has to
        %                      have the same index as the angle it
        %                      represents. Eg. input_angles=[0,45,-45,90]
        %                      and input_fractions=[50,20,20,10] is a
        %                      laminate with 50% of plies in the 0
        %                      direction, 20% at +45deg, 20% at -45deg and 10% at
        %                      90deg.
        %
        %  * t_laminate:       The laminate thickness. This is the continuous 
        %                      variable which will be changed by the optimization
        %                      process.
        %
        %  * t_ply_max:        The maximum ply thickness allowed. Necessary
        %                      to prevent function from making excessively
        %                      thick plies
        %
        %
        % OUTPUTS:
        %  * generated_layup: Thg layup which complies with the
        %                     angles and fractions indicated in the input.
        %
        %  * t_ply_required:  The ply thickness required for the above
        %                     layup to be exactly as thick as indicated 
        %                     in t_laminate
        
        function [final_layup, t_array_final] = layup_generator(input_angles, input_fractions, t_laminate, t_ply_max)
            
            % creates matrix where 1st column are all input angles, and 2nd
            % column their corresponding fraction
            main_list = [input_angles ; input_fractions]';
            
            % sorts matrix rows in ascending order (according to values of
            % 1st column, so angles)
            main_list_asc_angle = sortrows(main_list,1);
            
            % sorts matrix rows in ascending order (according to values of
            % 2nd column, so fractions)
            main_list_asc_frac = sortrows(main_list,2);
            
            % calculates min nr of plies required in laminate.
            nr_ply_min_thickness = floor(t_laminate/t_ply_max);
            
            % Figures out if more plies are necessary. The reason extra
            % plies may be required is if: 1) Current ply thickness is
            % too coarse to accurately represent the required orientation
            % fractions. 2) Current ply thickness does not allow to create
            % a symmetric layup
            
            % The code below finds the greatest common divisor between all
            % input fractions. Dividing 100(%) by the gcd yields the
            % minimum amount of plies required to accurately represent all 
            gcd_fraction = [];
            for i=1:length(input_fractions)
                gcd_current_row = gcd(main_list_asc_angle(:,2),main_list_asc_angle(i,2));
                gcd_fraction = [gcd_fraction; gcd_current_row];
            end
            min_gcd = min(min(gcd_fraction));
            nr_ply_min_fraction = 100/min_gcd;
            
%             nr_ply_req = max(nr_ply_min_thickness, nr_ply_min_fraction);
            if 2*nr_ply_min_fraction<nr_ply_min_thickness
                nr_ply_req = ceil(nr_ply_min_thickness/(2*nr_ply_min_fraction)  ) *(2*nr_ply_min_fraction);
            else
                nr_ply_req = 2* nr_ply_min_fraction; % added by simon as the half of the laminate needs to satisfy the fractions
            end
            % since the laminate has to be symmetric, if the required nr of
            % plies is odd, it has to be multiplied by 2, so that the
            % fraction can still be respected.
            %commented out by simon: multiplication is done previously.
%             if mod(nr_ply_req,2) ~= 0
%                 nr_ply_req = 2 * (nr_ply_req);
%             end
            
            % Since layup is always symmetric, we will only work on upper
            % half, and then mirror it.
            nr_ply_req_symm = nr_ply_req/2;
            
            
            
            % the code below checks if any extra plies of variable
            % thickness are needed in case certain orientations need a
            % non integer number of plies. This allows to comply with the
            % input fractions as well as the input maximum thickness
            % without creating a ridicolous number of ultra-thin plies.
            nr_plies_arr_asc_frac = main_list_asc_frac;
            nr_plies_arr_asc_frac(:,2) = main_list_asc_frac(:,2)*nr_ply_req_symm/100;
            excess_plies = mod(nr_plies_arr_asc_frac(:,2),1);
            
            
            % flag indicating which angles need fat plies. Follows order
            % of nr_plies_arr_asc_frac. The fat ply will be applied in
            % the form of extra thickness added to the innermost ply of
            % that specific angle in the symmetric layup.
            excess_plies_flag = excess_plies>0;
            
            % calculates nr of extra plies that need to be added for
            % perfect representation. These plies will have a thickness
            % different from the standard one.
            nr_excess_plies = 0;% (no plys have to be added as we just increase the thickness) old: sum(excess_plies>0)/2;
            
            
            % final nr of plies required in half the laminate
            nr_ply_req_symm_final = nr_ply_req_symm - nr_excess_plies;
            
            % initialize generated_layup output array
            generated_layup = ones(nr_ply_req_symm_final, 1);
            
            % calculate maximum ply thickness and initialize plt thickness array
            t_ply_max_temp = t_laminate/nr_ply_req;
            t_ply_arr = ones(nr_ply_req_symm,1) * t_ply_max_temp;
            
            
            % We want to distribute the plies as evenly as possible, while
            % complying with the following rules:
            %   
            %   1) The laminate must be symmetric
            %   2) The change in orientation between adjacent plies should
            %      be as gradual as possible
            %   3) If possible, no more than 4 adjacent plies should have
            %      the same orientation
            
            % calculates nr of plies belonging to angle with smallest
            % fraction and whether the outermost sublaminate is smaller
            % than the preceeding ones
%             nr_least_plies = floor(main_list_asc_frac(1,2)/100 * nr_ply_req_symm_final);
%             upper_sublaminate_incomplete_flag = mod((main_list_asc_frac(1,2)/100 * nr_ply_req_symm_final),1)>0;

            nr_least_plies = (main_list_asc_frac(1,2)/100 * nr_ply_req_symm);
            upper_sublaminate_incomplete_flag = mod(nr_least_plies,1)>0;
            
            
            
            % calculates how many plies of the remaining orientations
            % should be present within each sublaminate
            arr_plies_sublaminate = ceil(nr_plies_arr_asc_frac(:,2) / nr_least_plies);
            
            % calculates how many plies to fit between plies with smallest
            % fraction, so as to evenly distribute them.
            % if the direction with smallest fraction has only 1 ply, the
            % dist_least_plies is set to the nr of plies in the half
            % laminate minus 1
            if nr_least_plies == 1
                dist_least_plies = nr_ply_req_symm_final-1;
            else
%                 dist_least_plies = floor((nr_ply_req_symm_final-nr_least_plies)/(nr_least_plies-1));
                dist_least_plies = sum(arr_plies_sublaminate(2:end));
            end
            
            % distributes "rarest" plies evenly within the laminate
            if nr_least_plies == 1
                sublaminate_start_idx = [1];
            else
                sublaminate_start_idx = (1:dist_least_plies+1:length(generated_layup));
                sublaminate_start_idx = sublaminate_start_idx(1:nr_least_plies);
            end
            
            generated_layup(sublaminate_start_idx) = main_list_asc_frac(1,1);
            
            
            
            
            % loops over each other sublaminate in order to distribute all
            % other plies within the "sub-laminates" delimited by the
            % rarest plies
            if nr_least_plies == 1
                for i=1:length(sublaminate_start_idx)-upper_sublaminate_incomplete_flag

                    current_idx = sublaminate_start_idx(i);
                    prev_idx = 1;

                    for j=2:length(arr_plies_sublaminate)
                        generated_layup(current_idx+prev_idx:current_idx+arr_plies_sublaminate(j)+prev_idx-1) = nr_plies_arr_asc_frac(j,1);
                        prev_idx = arr_plies_sublaminate(j)+prev_idx;
    %                     generated_layup(current_idx+prev_idx:current_idx+arr_plies_sublaminate(j)+prev_idx) = nr_plies_arr_asc_frac(j,1);
    %                     prev_idx = arr_plies_sublaminate(j)+prev_idx;
                    end

                end
                
                
%             elseif nr_least_plies == 2
%                  
%                 i=1;
%                 current_idx = sublaminate_start_idx(i);
%                 prev_idx = 1;
% 
%                 for j=2:length(arr_plies_sublaminate)
%                     generated_layup(current_idx+prev_idx:current_idx+arr_plies_sublaminate(j)+prev_idx-1) = nr_plies_arr_asc_frac(j,1);
%                     prev_idx = arr_plies_sublaminate(j)+prev_idx;
%                 end
% 
%                 generated_layup(current_idx+prev_idx:end-1) = generated_layup(2:current_idx+prev_idx-1);
                                
                
            else
                for i=1:length(sublaminate_start_idx)-upper_sublaminate_incomplete_flag

                    current_idx = sublaminate_start_idx(i);
                    prev_idx = 1;

                    for j=2:length(arr_plies_sublaminate)
                        generated_layup(current_idx+prev_idx:current_idx+arr_plies_sublaminate(j)+prev_idx-1) = nr_plies_arr_asc_frac(j,1);
                        prev_idx = arr_plies_sublaminate(j)+prev_idx;
                    end

                end
            end
            
            
            % checks how many plies of each orientation are present, and
            % calculates the difference between how many should be present.
            % The leftover plies are used to fill the incomplete laminate
            % (if present)
            
            nr_current_plies_arr = sum(generated_layup == nr_plies_arr_asc_frac(:,1)');
            nr_leftover_plies_arr = floor(nr_plies_arr_asc_frac(:,2)') - nr_current_plies_arr;
            
            
            % if the uppermost sublaminate is incomplete, the leftover
            % plies are assigned to it.
            if upper_sublaminate_incomplete_flag == 1
                start_idx = sublaminate_start_idx(end);
                prev_idx = 1;
                
                for i = 1:length(nr_leftover_plies_arr)
                   
                    if nr_leftover_plies_arr(i)>0
                       
                        generated_layup(start_idx+prev_idx : start_idx + prev_idx + nr_leftover_plies_arr(i)-1) = nr_plies_arr_asc_frac(i,1);
                        prev_idx = nr_leftover_plies_arr(i) + prev_idx;
                       
                    % if some extra plies were erroneusly added, they will
                    % be overwritten thanks to the command below
                    elseif nr_leftover_plies_arr(i)<0
                        
                        start_idx = start_idx + nr_leftover_plies_arr(i);
                    end
                end
            end
            %change simon:
            % if the length is more than required, remove respective values
            % so that fractions match
            
            while length(generated_layup)>sum(nr_plies_arr_asc_frac(:,2))
                %get actual fractions
                actFrac=sum(generated_layup==input_angles)/length(generated_layup);
                %compare fractions to desired fractions
                devFrac=actFrac./input_fractions*100-1 ;
                %remove last entry of fraction with largest deviation
                angleToRemove=input_angles(find(devFrac==max(devFrac),1,'first'));
                
                generated_layup(find(generated_layup==angleToRemove,1,'last'))=[];
            end
            
            
            
            % creates the array with the final layup angles
            final_layup = [generated_layup(end:-1:1);generated_layup];
            
            
            % now we only have to correct the thicknesses of the plies with
            % extra thickness
            
            for i=1:length(excess_plies_flag)
                
                if excess_plies_flag(i)==1
                    
                    current_angle = nr_plies_arr_asc_frac(i,1);
                    extra_thick_ply_index = find(generated_layup == current_angle,1);
                    t_ply_arr(extra_thick_ply_index) = t_ply_arr(extra_thick_ply_index) * (1 + excess_plies(i));     
                end
            end
            
            t_array_final = [t_ply_arr(end:-1:1);t_ply_arr];
            
        end
        
        function  t_array_final = layup_generator_simple(input_t_array, new_thickness)
            
            old_thickness = sum(input_t_array);
            multiplier = new_thickness/old_thickness;
            
            t_array_final = input_t_array * multiplier;
            
        end
        

    end

    methods
        
        % Class constructor, takes layup angles (in degrees) as input, as
        % well as material properties of each ply. All methods within this
        % class are called once the constructor is called.
        
        % The class assumes that the material of each ply is always the
        % same.
        
        % material_properties_input is an array which should have the
        % following form: [E_x, E_y, mu_xy, mu_yx, G_xy, ply_thicness]
        function obj = class_laminate(angle_lst_deg_in, ply_properties_input, ply_thick_array, design_envelope_u_arr,design_envelope_u_arr_khani, offset_angle)

            % saving input layup and flipping the sign. Sign inversion done
            % for internal purposes, the final result is the correct one
            % associated with the input layup.
            obj.angle_lst_deg_original = angle_lst_deg_in;
            
%             obj.angle_lst_deg = -angle_lst_deg_in - offset_angle;
            obj.angle_lst_deg = -angle_lst_deg_in - offset_angle;
            % checks nr of plies in order to properly generate the arrays
            % below
            mat_dim = size(obj.angle_lst_deg);

            % generates arrays where the properties of each ply are stored.
            obj.E_x_arr = ones(mat_dim) * ply_properties_input(1);
            obj.E_y_arr= ones(mat_dim) * ply_properties_input(2);
            obj.mu_xy_arr= ones(mat_dim) * ply_properties_input(3);
            obj.mu_yx_arr= ones(mat_dim) * ply_properties_input(4);
            obj.G_xy_arr = ones(mat_dim) * ply_properties_input(5);
            obj.thick_arr = ply_thick_array;
            
            % assinging material invariants to laminate
            obj.u1 = design_envelope_u_arr(1);
            obj.u2 = design_envelope_u_arr(2);
            obj.u3 = design_envelope_u_arr(3);
            obj.u4 = design_envelope_u_arr(4);
            obj.u5 = design_envelope_u_arr(5);
            obj.u6 = design_envelope_u_arr(6);
            %khani invariants
            obj.u1k = design_envelope_u_arr_khani(1);
            obj.u2k = design_envelope_u_arr_khani(2);
            obj.u3k = design_envelope_u_arr_khani(3);
            obj.u4k = design_envelope_u_arr_khani(4);
            obj.u5k = design_envelope_u_arr_khani(5);
            obj.u6k = design_envelope_u_arr_khani(6);
            
            % saving offset angle
            obj.offset_angle = offset_angle;
            
            % calls next function
            obj = obj.laminate_maker();
            
        end

        function obj = ABD_maker(obj)

            sum_A = 0;
            sum_B = 0;
            sum_D = 0;

            for counter_ABD = 1:length(obj.angle_lst_rad)

                sum_A = sum_A + (obj.C_phi(:,:,counter_ABD)*(obj.z_up(counter_ABD) - obj.z_low(counter_ABD)));
                sum_B = sum_B + 1./2.*(obj.C_phi(:,:,counter_ABD) * (obj.z_up(counter_ABD)^2 - obj.z_low(counter_ABD)^2));
                sum_D = sum_D + 1./3.*(obj.C_phi(:,:,counter_ABD) * (obj.z_up(counter_ABD)^3 - obj.z_low(counter_ABD)^3));

            end

            obj.ABD_stiff = [sum_A, sum_B ; sum_B, sum_D];
            obj.abd_compl = inv(obj.ABD_stiff);

        end



        % Turns the group of plies given in the constructor into a
        % laminate. Saves their location within the laminate, calculates
        % the stiffnes and compliance matrices of each ply both in the
        % local ply coordinate system (of each respective ply) and the
        % global laminate coordinate system.
        function obj = laminate_maker(obj)

            % converts layup angles array into radians
            obj.angle_lst_rad = deg2rad(obj.angle_lst_deg);

            % calculates total thickness of layup and stores it in obj.t_lam
            thickness_tot = sum(obj.thick_arr);
            obj.t_lam = thickness_tot;

            % calculates mid-plane z-coordinates wrt the bottom of the
            % layup
            z_mid = thickness_tot/2.;

            z_prev = 0.;
            
            % loops through all plies
            for counter = 1:length(obj.angle_lst_rad)
                
                % assings an ID to each ply (bottom ply = nr of plies =1,
                % top ply = 0)
                obj.ID(counter) = length(obj.angle_lst_rad)-counter;
                
                % stores each ply angle (wrt the global laminate coord
                % system) into a separate array
                obj.phi_rad(counter) = obj.angle_lst_rad(counter);

                % Calculates bottom and top z-coordinate of each ply (wrt
                % the global lamiante coord system)
                obj.z_low(counter) = z_prev - z_mid;
                obj.z_up(counter) = obj.z_low(counter) + obj.thick_arr(counter);

                z_prev  =  z_prev + obj.thick_arr(counter);
            end
            
            % initializes stiffnes and compliance matrices of plies
            obj.C_loc = zeros(3,3,length(obj.angle_lst_rad));
            obj.C_phi = zeros(3,3,length(obj.angle_lst_rad));
            obj.S_loc = zeros(3,3,length(obj.angle_lst_rad));
            obj.S_loc_inverted = zeros(3,3,length(obj.angle_lst_rad));
            obj.S_phi_inverted = zeros(3,3,length(obj.angle_lst_rad));
            obj.S_phi = zeros(3,3,length(obj.angle_lst_rad));
            
            % loops over all plies
            for counter2 = 1:length(obj.angle_lst_rad)
                
                % creates ply stiffness matrix in local ply coordinate
                % system
                obj.C_loc(:,:,counter2) = obj.f_C_mat_local(obj.mu_xy_arr(counter2), obj.mu_yx_arr(counter2), obj.E_x_arr(counter2), obj.E_y_arr(counter2), obj.G_xy_arr(counter2));
                
                % creates ply compliance matrix in local ply coordinate
                % system, as well as its inverse
                obj.S_loc(:,:,counter2) = obj.f_S_mat_local(obj.mu_xy_arr(counter2), obj.mu_yx_arr(counter2), obj.E_x_arr(counter2), obj.E_y_arr(counter2), obj.G_xy_arr(counter2));
                obj.S_loc_inverted(:,:,counter2) = inv(obj.S_loc(:,:,counter2));
                
                % creates ply stiffness matrix in global laminate
                % coordinate system
                obj.C_phi(:,:,counter2) = obj.f_C_phi(obj.C_loc(:,:,counter2), obj.phi_rad(counter2));
                
                % creates ply compliance matrix in global laminate
                % coordinate system, as well as its inverse
                obj.S_phi_inverted(:,:,counter2) = obj.f_S_phi(obj.S_loc_inverted(:,:,counter2), obj.phi_rad(counter2));
                obj.S_phi(:,:,counter2) = obj.f_S_phi(obj.S_loc(:,:,counter2), obj.phi_rad(counter2));
            end

            % calls function which creates the ABD matrices
            obj = obj.ABD_maker();
        end

        % function used to create stress rotation matrices
        function M_mat = f_rotation_mat(obj,theta_rot)

            m = cos(theta_rot);
            n = sin(theta_rot);

            M_mat = [m*m, n*n, 2*m*n; n*n,m*m,-2*m*n; -m*n,m*n,m*m-n*n];

        end
        
        % function used to create strain rotation matrices
        function M_mat_strain = f_rotation_mat_strain(obj,theta_rot)

            m = np.cos(theta_rot);
            n = np.sin(theta_rot);

            M_mat_strain = [m*m,n*n,m*n; n*n,m*m,-m*n; -2*m*n,2*m*n,m*m-n*n];

        end

        % Function that creates the ply stiffness matrix (in local ply
        % coordinate system) with the ply properties given in the input
        function C_loc = f_C_mat_local(obj,mu_xy, mu_yx, E_x, E_y,G_xy)

            den = 1. - mu_xy*mu_yx;
            C_loc = [E_x / (den),  mu_xy*E_y/den,  0;
                     mu_yx*E_x/den,  E_y / (den),  0;
                     0,               0,          G_xy];
        end

        
        % Function that creates the ply compliance matrix (in local ply
        % coordinate system) with the ply properties given in the input
        function S_loc = f_S_mat_local(obj,mu_xy, mu_yx, E_x, E_y, G_xy)

            S_loc = [1/E_x,  -mu_yx/E_y,  0;
                     -mu_xy/E_x,  1./E_y,  0;
                     0,               0,   1./G_xy];
        end

        % Function that performs the rotation of a ply stiffness matrix
        % from the local ply coord. sys to the global laminate one
        function C_phi_loc = f_C_phi(obj,C_input, theta)

            M_loc = obj.f_rotation_mat(theta);
            M_loc_T = (M_loc)';

            C_phi_loc = mtimes(M_loc, (C_input * M_loc_T));

        end

        
        % Function that performs the rotation of a ply compliance matrix
        % from the local ply coord. sys to the global laminate one
        function S_phi_loc = f_S_phi(obj, S_input, theta)

            M_loc = obj.f_rotation_mat(theta);
            M_loc_T = M_loc';

            S_phi_loc = mtimes(M_loc, (S_input * M_loc_T));

        end
        
        % Calculates the failure index of the specific laminate instance.
        % This includes the safety factor. Since the
        % Material invariants  are the same for the whole cross section, it
        % assumes that ALL LAMINATES WITHIN THE CROSS SECTION ARE MADE OF
        % THE SAME PLY MATERIAL.
        
        function obj = calculate_safety_factor(obj, strain_in, strain_in_endA, strain_in_endB, safety_factor_input)
            
            
            critical_safety_factor_arr=zeros(2,size(strain_in,2));
                            
            C0=-(1/4)*obj.u6k^2/obj.u4k-1;
            C1=-(1/2)*obj.u3k*obj.u6k/obj.u4k+obj.u5k;
            C11=-(1/4)*obj.u3k^2/obj.u4k+obj.u1k+obj.u2k;
            C12=obj.u1k-(1/4)*obj.u3k^2/obj.u4k;
            
            
            if obj.failureEnvelopeType==1

                for i=1:1:size(strain_in,2)
                    I1 =  strain_in(1,i) + strain_in(2,i); % volumetric strain invariant
                    I2 =  sqrt(((strain_in(1,i)-strain_in(2,i))/2)^2 + strain_in(3,i)^2 ); % maximum shear strain
                    if and((I1==0),(I2==0)) % no strains
                       critical_safety_factor=Inf;
                    else
                        obj.fail_idx_crit = zeros(size(strain_in,2),1);
                        % calculating failure indices for the provided strain array
                        % according to IJsselmuiden 2008
                        % calculating coefficients to solve the 2nd and 4th order
                        % polynomial equations (to find fail index)
                        a10 = obj.u4^2 + 4*obj.u1 - 4*obj.u6;
                        a11 = -4*obj.u2 * I1 * (obj.u1 - obj.u6) + 2 * obj.u4 * obj.u5 * I1;
                        a12 = 4*obj.u6^2 * I2^2 - 4*obj.u3 * I1^2 * (obj.u1 - obj.u6) - 4 * obj.u6 * obj.u1 * I2^2 + obj.u5^2 * I1^2;
                        a20 = 1;
                        a21 = -2*obj.u2 * I1;
                        a22 = -2*obj.u3 * I1^2 + obj.u2^2 * I1^2 - I2^2 * (obj.u4^2 + 2*obj.u1);
                        a23 = 2*obj.u2 * I1^3 * obj.u3 - I2^2 * (2*obj.u4 * obj.u5 * I1 - 2*obj.u1 * obj.u2 * I1);
                        a24 = obj.u1^2 * I2^4 - I2^2 * (obj.u5^2 * I1^2 - 2*obj.u1 * obj.u3 * I1^2) + obj.u3^2 * I1^4;

                        % 1st and 2nd roots, from 2nd order polynomial
                        root_1 = (-a11 + sqrt(a11^2 - 4 * (a12 * a10))) / (2 * a12);
                        root_2 = (-a11 - sqrt(a11^2 - 4 * (a12 * a10))) / (2 * a12);

                        % 3rd to 6th roots, from 4th order polynomial
                        if any(isnan([a24,a23,a22,a21,a20]))
                            root_arr=[Inf,Inf,Inf,Inf];
                        elseif any(isinf([a24,a23,a22,a21,a20]))
                            root_arr=[Inf,Inf,Inf,Inf];
                        else
                        %                         root_arr = roots([a24,a23,a22,a21,a20]);
                            a = diag(ones(1,3),-1);
                            a(1,:) = -[a23,a22,a21,a20]./a24;
                            root_arr = eig(a);
                        end
                        root_3 = root_arr(1);
                        root_4 = root_arr(2);
                        root_5 = root_arr(3);
                        root_6 = root_arr(4);

                        % storing all roots in an array
                        safety_factor_array = [root_1,root_2,root_3,root_4,root_5,root_6];

                        % The smallest of all positive roots is the critical safety
                        % factor (if equal to 1, on the verge of failure, smaller
                        % than 1 = failure, bigger than 1 = safe).
                        critical_safety_factor = min(safety_factor_array(safety_factor_array>0));
                        % Saving the critical safety factor of this shell element in
                        % an array. Incorporating safety factor of this cross
                        % section within the failure index.
                        %                fail_idx_arr(i) = -(critical_safety_factor/safety_factor_input - 1);

                        obj.fail_idx_arr(i) = (1/(critical_safety_factor/safety_factor_input)^2 - 1);
                    end
                end
            elseif obj.failureEnvelopeType==2  
                for i=1:1:size(strain_in,2)
                    I1 =  strain_in(1,i) + strain_in(2,i); % volumetric strain invariant
                    I2 =  sqrt(((strain_in(1,i)-strain_in(2,i))/2)^2 + strain_in(3,i)^2 ); % maximum shear strain
                    if and((I1==0),(I2==0)) % no strains
                       critical_safety_factor=Inf;
                    else
                        % calculating failure indices for the provided strain array
                        % according to A. Khani 2011 (check also phd thesis
                        % in tudelft repository)
                        %constants
                        for iEnd=1:2

                            %strains
                            if iEnd==1
                                ex=strain_in_endA(1,i);
                                ey=strain_in_endA(2,i);
                                gamma=strain_in_endA(3,i);
                            elseif iEnd==2

                                ex=strain_in_endB(1,i);
                                ey=strain_in_endB(2,i);
                                gamma=strain_in_endB(3,i);
                            end

                            %mohr circle eq
                            eavg=(ex+ey)/2; %average
                            R=sqrt(((ex-ey)/2)^2+(gamma/2)^2); %radius

                            %principal strains
                            e1=eavg+R;
                            e2=eavg-R;

                            % this is from phd thesis of khani
                            a10=C0;
                            a11=C1*e1+C1*e2;
                            a12=C11*e1^2+C11*e2^2+2*C12*e1*e2;

                            % 1st and 2nd roots, from 2nd order polynomial
                            root_1 = (-a11 + sqrt(a11^2 - 4 * a12 * a10)) / (2 * a10);
                            root_2 = (-a11 - sqrt(a11^2 - 4 * a12 * a10)) / (2 * a10);

                            % storing all roots in an array
                            safety_factor_array = 1./[root_1,root_2];

                            % The smallest of all positive roots is the critical safety
                            % factor (if equal to 1, on the verge of failure, smaller
                            % than 1 = failure, bigger than 1 = safe).
                            critical_safety_factor_arr(iEnd,i) = min(safety_factor_array(safety_factor_array>0));


                        end
                    end
                end
            end
            obj.fail_idx_arr=(1./(critical_safety_factor_arr/safety_factor_input).^2 -1);

            % Saving the most critical safety factor AMONG ALL SHELL
            % ELEEMNTS in an array. Todo, KS-aggregation
            if obj.failureEnvelopeType==1
                obj.fail_idx_crit = max(obj.fail_idx_arr);
            elseif obj.failureEnvelopeType==2
                obj.fail_idx_crit = constAgg(obj.fail_idx_arr(:),500);
            end
        end
        
        function obj = calculate_safety_factor_from_princ(obj,e1,e2,sf, aggregationFactor)
            C0=-(1/4)*obj.u6k^2/obj.u4k-1;
            C1=-(1/2)*obj.u3k*obj.u6k/obj.u4k+obj.u5k;
            C11=-(1/4)*obj.u3k^2/obj.u4k+obj.u1k+obj.u2k;
            C12=obj.u1k-(1/4)*obj.u3k^2/obj.u4k;
            
            a10=C0;
            
            a11=C1*e1+C1*e2;
            a12=C11*e1.^2+C11*e2.^2+2*C12*e1.*e2;
            
            % 1st and 2nd roots, from 2nd order polynomial
            root_Part =sqrt(a11.^2 - 4 * a12 .* a10) ;
            rootArr=(([-a11-root_Part, -a11+root_Part])/(2*a10));
            % storing all roots in an array
            rootArr2=1./rootArr;
            obj.fail_idx_arr=(1./(rootArr2(rootArr2>0)/sf).^2 -1);
            obj.fail_idx_crit = constAgg(obj.fail_idx_arr(:),aggregationFactor);
        end
        
        function obj = update_laminate(obj)
            
            %first, all the inputs necessary for each laminate have to be
            %generated again
            
            % saving original layup of original laminates
            original_layup = obj.angle_lst_deg_original;
            
            % saving ply properties of original laminates
            ply_properties = [obj.E_x_arr(1), obj.E_y_arr(1), obj.mu_xy_arr(1), obj.mu_yx_arr(1), obj.G_xy_arr(1)];

            % Calculating new ply thickness arrays (inputs required: old
            % laminate ply thickness array, and updated laminate thickness
            % value)
            thickness_arr = obj.layup_generator_simple(obj.thick_arr, obj.t_lam);
            
            
            % Saving material invariants of old laminate (needed for
            % design envelope)
            u_arr = [obj.u1, obj.u2, obj.u3, obj.u4, obj.u5, obj.u6];
            u_arr_khani = [obj.u1k, obj.u2k, obj.u3k, obj.u4k, obj.u5k, obj.u6k];
            
            
            % Saving new offset angle (the new value should have been set
            % before this function is called)
            offset_angle_in = obj.offset_angle;
            
            
            obj = obj.execute_laminate_update(original_layup, ply_properties, thickness_arr, u_arr,u_arr_khani, offset_angle_in);

        end
        
        
        % function which updates the lamintes according to the
        % optimization output.
        % BEFORE CALLING THIS FUNCTION, THE t_lam AND offset_angle
        % PROPERTIES (BOTH BELONGING TO LAMINATE OBJECTS) NEED TO HAVE BEEN
        % SET TO THE UPDATED VALUES DETERMINED BY THE OPTIMIZATION PROCESS.
        % THE FUNCTION SHOULD BE CALLED FROM ITS CROSS SECTION PARENT.
        function obj = execute_laminate_update(obj,angle_lst_deg_in, ply_properties_input, ply_thick_array, design_envelope_u_arr,design_envelope_u_arr_khani, offset_angle)
           
            % saving input layup and flipping the sign. Sign inversion done
            % for internal purposes, the final result is the correct one
            % associated with the input layup.
            obj.angle_lst_deg_original = angle_lst_deg_in;
            
            obj.angle_lst_deg = -angle_lst_deg_in - offset_angle;
            % checks nr of plies in order to properly generate the arrays
            % below
            mat_dim = size(obj.angle_lst_deg);

            % generates arrays where the properties of each ply are stored.
            obj.E_x_arr = ones(mat_dim) * ply_properties_input(1);
            obj.E_y_arr= ones(mat_dim) * ply_properties_input(2);
            obj.mu_xy_arr= ones(mat_dim) * ply_properties_input(3);
            obj.mu_yx_arr= ones(mat_dim) * ply_properties_input(4);
            obj.G_xy_arr = ones(mat_dim) * ply_properties_input(5);
            obj.thick_arr = ply_thick_array;
            
            % assinging material invariants to laminate
            obj.u1 = design_envelope_u_arr(1);
            obj.u2 = design_envelope_u_arr(2);
            obj.u3 = design_envelope_u_arr(3);
            obj.u4 = design_envelope_u_arr(4);
            obj.u5 = design_envelope_u_arr(5);
            obj.u6 = design_envelope_u_arr(6);
            obj.u1k = design_envelope_u_arr_khani(1);
            obj.u2k = design_envelope_u_arr_khani(2);
            obj.u3k = design_envelope_u_arr_khani(3);
            obj.u4k = design_envelope_u_arr_khani(4);
            obj.u5k = design_envelope_u_arr_khani(5);
            obj.u6k = design_envelope_u_arr_khani(6);
            
            % calls next function
            obj = obj.laminate_maker(); 
        end
    
    end
end


