%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%
% ======================================================================
%> @brief class for a finite beam element
%> This class is the base class for a finite element beam, the follwing
%> coordinate system is used
%> Global Coordinate Definition (according to Airbus Wing Geometry Definitions ):
%>   y= along starboard wing
%>   x= top view pointing aftwards
%>   z= top view pointing upwards
%> Local Coordinate Definition (Cartesian Coordinate System )
%>   y= beam axis
%>   x= top view pointing aftwards
%>   z= top view pointing upwards
% ======================================================================



% This class is a re-adaptment of the class_beam_element class.
% The class was modified in order to properly handle anisotropic materials
classdef class_beamelement_anisotropic
    properties
        %% general element information
        
        % The mechanical properties given refer to a generic 0/90 QI CFRP
        %> element length                                [m]
        le=0.0;
        %> element rotation around z axis (sweep)        [rad]
        phi=0.0;
        %> element rotation around x axis (dihedral)     [rad]
        nu=0.0;
        %> element rotation around y axis (twist)        [rad]
        epsilon=0.0;
        %> Youngs modulus                                [N/m2]
        E1;
        E2;
        %> Shear modulus                                 [N/m2]
        G12;
        %> element cross-section area                    [m2]
        A=1;
        
        %> excentrictiy
        e=0;
        
        J=1; %torsional moment
        
        %% Forces and Moments (in element local coordinates)
        %> element transverse pressure load              [N/m]
        qz=0.0;
        %>
        dqz=0.0;
        
        %> element transverse pressure load              [N/m]
        qy=0.0;
        
        dqy=0.0;
        %> element transverse pressure load              [N/m]
        qx=0.0;
        dqx=0.0;
        %> element distributed torsion load              [Nm/m]
        mt=0.0;
        dmt=0.0;
        %> element distributed moment load               [Nm/m]
        mx=0.0;
        dmx=0.0;
        %> element distributed moment load               [Nm/m]
        mz=0.0;
        dmz=0.0;
        %> element distributed load due to eigenmass     [Nm/m]
        qm=0.0;
        
        %% accelerations defined in global coordinate system
        %> element distributed load due to x acceleration     [Nm/m]
        ax=0.0;
        %> element distributed load due to y acceleration       [Nm/m]
        ay=0.0;
        %> element distributed load due to z acceleration      [Nm/m]
        az=0.0;
        
        
        %> element distrubuted load due to fuel          [Nm/m]
        qf=0.0;
        
        %% additions for non linear element
        elp=0.0;
        elpglobal=0.0;              % element internal force vector
        nodal_deflections_loc;
        
        %% Masses
        %>  element structural mass per unit length       [kg/m]
        m=1;
        
        %> element fueled y/n                     []
        is_fueled=0;
        %> fuel volume per element                [m3]
        el_fuel_vol=0.0
        %> enclosed area
        A_enclosed=0;
        
        sigma_b=0.0;
        %> element degrees of freedom
        eldof;
        %> element stiffness matrix in local coordinates
        elK;
        %> element stiffness matrix in global system coordinates
        elKglobal;
        %> element mass matrix
        elM;
        %> element mass matrix in global coordinates
        elMglobal;
        %> element mass matrix
        elM_lumped;
        %> element mass matrix in global coordinates
        elM_lumped_global;
        %> element load vector in local coordinates
        elq;
        %> element load vector in global coordinates
        elqglobal;
        %> element mass
        el_m_s=0.0;
        %> element system mass
        el_m_sys=0.0;
        el_m_p=0.0;
        
        %> element fuel mass
        el_m_fuel=0.0;
        
        %> transformation matrix from global system to local system
        T;
        
        Nm=0;
        Vm;
        Mm;
        kappa;
        
        % The crosssection is where objects belonging to
        % class_crosssection_wingbox_anisotropic are stored.
        crosssection;
        %> reference to beam to which the element belongs to
        ref_parent_beam;
        %> flag for massless beam elements
        is_massless;
        
        % Arrays containing all shell strains and stresses recovered through the cross
        % sectional modeler. The elements are ordered anti-clockwise,
        % starting from the bottom point of the front spar, and ending with
        % the same point, but approached from the right-most point of the
        % bottom skin.
        % Each element has 6DoF, meaning that 6DoFs are present per element
        % (the arrays are 1D, so each element has 6 entries, apart from the von mises ones).
        element_strain_crossmod;         % strain array
        element_strain_von_mis_crossmod; % von mises strain array
        element_stress_crossmod;         % stress arrray
        element_stress_von_mis_crossmod; % von mises stress array
        
        % Flag indicating whether the current beam element is composed of
        % anisotropic materials. Used to quickly check whether beam element
        % objects belong to class_beamelement_anisotropic or
        % class_beamelement.
        anisotropic = 1;
        
        
    end
    
    methods
        function obj =class_beamelement_anisotropic(crosssection, parent_ref)
            % Saving the input cross section object as a class property.
            obj.crosssection=crosssection;
            % Saving the reference to the parent beam as a class property.
            obj.ref_parent_beam=parent_ref;
        end
        
        %compute element stiffness matrix for 6 DOF element
        obj=lin_elK_6dof_anisotropic(obj);
        
        % Compute element mass matrix for 6DOF element
        obj=elM_6dof_anisotropic(obj);
        
        obj=elM_lumped_6dof_anisotropic(obj);
        
        %compute element load vector for 6 DOF element
        obj=elq_6dof(obj,add_eigenmass,add_fuelmass);
        
        % Calculating cross sectional properies of the beam element.
        % PLEASE NOTE that all the inertia properties obtained from this function
        % are not actually used, as the cross sectional modeller calculates
        % its own. The calculated area is used.
        obj = f_calcCrossProp_anisotropic(obj)
        
        
        %set geometry of element
        function obj=setElementGeometry(obj,le,phi,nu,twist)
            obj.le=le;
            obj.phi=phi;
            obj.nu=nu;
            obj.epsilon=twist;
            obj.T=obj.f_rotM_6dof_anisotropic(real(obj.nu),real(obj.epsilon),real(obj.phi));
        end
        
        
        % Conpute ONLY the shell strains in the cross-section from the
        % Euler disp. of the beam
        function obj = f_calc_strain_crossmod(obj,nodal_deflection_input,ndof_node)
            euler_delf_node_1 = [nodal_deflection_input(2),-nodal_deflection_input(4), nodal_deflection_input(6), nodal_deflection_input(5)];
            euler_delf_node_2 = [nodal_deflection_input(8),-nodal_deflection_input(10), nodal_deflection_input(12), nodal_deflection_input(11)];
            
            euler_defl_rel = (euler_delf_node_2 - euler_delf_node_1)/obj.le;
            
            obj.element_strain_crossmod = obj.crosssection.Gamma_euler * euler_defl_rel';
            
            el_idx_arr = [sum(obj.crosssection.shell_id_arr==1), sum(obj.crosssection.shell_id_arr==2), sum(obj.crosssection.shell_id_arr==3), sum(obj.crosssection.shell_id_arr==4)];
            el_idx_arr = cumsum(el_idx_arr);
                       
            strain_fs_unformatted = obj.element_strain_crossmod(1:el_idx_arr(1)*ndof_node);
            strain_rs_unformatted = obj.element_strain_crossmod(el_idx_arr(2)*ndof_node+1:el_idx_arr(3)*ndof_node);
            strain_sk_up_unformatted = obj.element_strain_crossmod(el_idx_arr(1)*ndof_node+1:el_idx_arr(2)*ndof_node);
            strain_sk_lo_unformatted = obj.element_strain_crossmod(el_idx_arr(3)*ndof_node+1:el_idx_arr(4)*ndof_node);
                       
            obj.crosssection.strain_fs = reshape(strain_fs_unformatted, [ndof_node, length(strain_fs_unformatted)/ndof_node]);
            obj.crosssection.strain_sk_up = reshape(strain_sk_up_unformatted, [ndof_node, length(strain_sk_up_unformatted)/ndof_node]);
            obj.crosssection.strain_rs = reshape(strain_rs_unformatted, [ndof_node, length(strain_rs_unformatted)/ndof_node]);
            obj.crosssection.strain_sk_lo = reshape(strain_sk_lo_unformatted, [ndof_node, length(strain_sk_lo_unformatted)/ndof_node]);
        end
        
        
        %Compute ONLY shell strains in the CS from the strains computed in
        %f_calc_strain_crossmod
        %Cannot be called without the strain function
        function obj = f_calc_stress_crossmod(obj)
            
            % The stresses are obtained by taking the
            % strains and multipliying them with the ABD stiffness matrices of their
            % respective skin/spar
            
            el_idx_arr = [sum(obj.crosssection.shell_id_arr==1), sum(obj.crosssection.shell_id_arr==2), sum(obj.crosssection.shell_id_arr==3), sum(obj.crosssection.shell_id_arr==4)];
            el_idx_arr = cumsum(el_idx_arr);
            
            strain_von_mis_1 = obj.element_strain_crossmod(1:6:end).^2 + obj.element_strain_crossmod(2:6:end).^2;
            strain_von_mis_2 = obj.element_strain_crossmod(1:6:end) .* obj.element_strain_crossmod(2:6:end);
            strain_von_mis_3 = obj.element_strain_crossmod(3:6:end).^2;
            
            obj.element_strain_von_mis_crossmod = sqrt(4/9 * (strain_von_mis_1 - strain_von_mis_2) + 1/3 * strain_von_mis_3);
            
            obj.crosssection.stress_fs = (obj.crosssection.laminate_fs.ABD_stiff(1:3,1:3) * obj.crosssection.strain_fs(1:3,:))/obj.crosssection.t_sp_fr;
            obj.crosssection.stress_sk_up = (obj.crosssection.laminate_sk_up.ABD_stiff(1:3,1:3) * obj.crosssection.strain_sk_up(1:3,:))/obj.crosssection.t_sk_up;
            obj.crosssection.stress_rs = (obj.crosssection.laminate_rs.ABD_stiff(1:3,1:3) * obj.crosssection.strain_rs(1:3,:))/obj.crosssection.t_sp_re;
            obj.crosssection.stress_sk_lo = (obj.crosssection.laminate_sk_lo.ABD_stiff(1:3,1:3) * obj.crosssection.strain_sk_lo(1:3,:))/obj.crosssection.t_sk_lo;
            
            obj.element_stress_crossmod = [obj.crosssection.stress_fs(:);
                                           obj.crosssection.stress_sk_up(:);
                                           obj.crosssection.stress_rs(:);
                                           obj.crosssection.stress_sk_lo(:)];
            
           
            % GENERAL PLANE STRESS VON MISES
            von_mis_11 = obj.element_stress_crossmod(1:3:end);
            von_mis_22 = obj.element_stress_crossmod(2:3:end);
            von_mis_12 = obj.element_stress_crossmod(3:3:end);
            
            % GENERAL VON MISES
            %             obj.element_stress_von_mis_crossmod = sqrt(0.5 * (von_mis_11_22 + von_mis_22_33 + von_mis_11_33 + 6 * von_mis_shear));
            
            % GENERAL PLANE STRESS VON MISES
            obj.element_stress_von_mis_crossmod = sqrt(von_mis_11.^2 - von_mis_11 .* von_mis_22 + von_mis_22.^2 + 3*von_mis_12.^2);
            
            obj.crosssection.stress_von_mis_fs = obj.element_stress_von_mis_crossmod(1:el_idx_arr(1))';
            obj.crosssection.stress_von_mis_rs = obj.element_stress_von_mis_crossmod(el_idx_arr(2)+1:el_idx_arr(3))';
            obj.crosssection.stress_von_mis_sk_up = obj.element_stress_von_mis_crossmod(el_idx_arr(1)+1:el_idx_arr(2))';
            obj.crosssection.stress_von_mis_sk_lo = obj.element_stress_von_mis_crossmod(el_idx_arr(3)+1:el_idx_arr(4))';
            
        end
        
        
        %Compute safety factor when asked to 
        function obj = f_calc_safety_factor(obj)
                obj.crosssection.laminate_sk_up = obj.crosssection.laminate_sk_up.calculate_safety_factor(obj.crosssection.strain_sk_up, obj.crosssection.safety_factor);
                obj.crosssection.laminate_sk_lo = obj.crosssection.laminate_sk_lo.calculate_safety_factor(obj.crosssection.strain_sk_lo, obj.crosssection.safety_factor);
                obj.crosssection.laminate_fs = obj.crosssection.laminate_fs.calculate_safety_factor(obj.crosssection.strain_fs, obj.crosssection.safety_factor);
                obj.crosssection.laminate_rs = obj.crosssection.laminate_rs.calculate_safety_factor(obj.crosssection.strain_rs, obj.crosssection.safety_factor);            
        end
        
        % Thsi function uses the Euler Bernoulli beam's nodal displacements
        % found by solving the linear system and combines them with
        % matrices calculated in the cross sectional modeler in order to
        % recover the shell strains from the Euler displacements.
        function obj = f_calc_stress_strain_crossmod(obj,nodal_deflection_input,ndof_node)
            
            %check calling function (same name, within class_beam) to check
            %whether the euler_defl stuff are either the deformations or
            %loads
            
            % Using gamma euler
            % going to temporarily change the order of the deformations
            % below, to see if perhaps they should have been arranged as in
            % eq12 from abdalla's paper
            
            % Saving the Euler nodal deflections of node 1 (start of beam
            % element) and node 2 (end of beam element) in separate arrays.
            % The order of the DOFs is: [Axial, Flap, Lag, Torsion]
            euler_delf_node_1 = [nodal_deflection_input(2),-nodal_deflection_input(4), nodal_deflection_input(6), nodal_deflection_input(5)];
            euler_delf_node_2 = [nodal_deflection_input(8),-nodal_deflection_input(10), nodal_deflection_input(12), nodal_deflection_input(11)];
            
            %equation12
            %             euler_delf_node_1 = [nodal_deflection_input(2),nodal_deflection_input(5), -nodal_deflection_input(4), nodal_deflection_input(6)];
            %             euler_delf_node_2 = [nodal_deflection_input(8), nodal_deflection_input(11),-nodal_deflection_input(10), nodal_deflection_input(12)];
            
            % calculating the Euler strains. These are just the difference
            % between Euler displacements at nodes 2 and 1, divided by the
            % element length (in order to turn displacements into strains).
            euler_defl_rel = (euler_delf_node_2 - euler_delf_node_1)/obj.le;
            % Multipying the Gamma_euler matrix obtained in the cross
            % sectional modeller with the euler nodal strains yields the
            % shell strains. The first 3DoFs are planar x, y and xy
            % strains. The last 3DoFs are curvatures of the shell normal.
            
            obj.element_strain_crossmod = obj.crosssection.Gamma_euler * euler_defl_rel';
            
            %             obj.element_strain_crossmod = obj.crosssection.Gamma * euler_defl_rel';
            
            
            %             nr_shells = length(obj.crosssection.shell_id_arr);
            % The strains have to be assigned to their respective
            % skin/spar. The array is 1D, the elements are ordered from the
            % bottom-most one in the front spar to the right-most one in
            % the bottom skin, going in a counter-clockwise direction
            % (looking from the tip of the beam to the root).
            
            % The el_idx_arr array containts 4 entries. They are the
            % indices where entries element_strain_crossmod switch from fs
            % to sk_up, from sk_up to rs, from rs to sk_lo and the final
            % index of the array.
            
            el_idx_arr = [sum(obj.crosssection.shell_id_arr==1), sum(obj.crosssection.shell_id_arr==2), sum(obj.crosssection.shell_id_arr==3), sum(obj.crosssection.shell_id_arr==4)];
            el_idx_arr = cumsum(el_idx_arr);
            
            % The element_strain_crossmod is split into 4 separate 1D
            % array, 1 for each spar/skin.
            
            strain_fs_unformatted = obj.element_strain_crossmod(1:el_idx_arr(1)*ndof_node);
            strain_rs_unformatted = obj.element_strain_crossmod(el_idx_arr(2)*ndof_node+1:el_idx_arr(3)*ndof_node);
            strain_sk_up_unformatted = obj.element_strain_crossmod(el_idx_arr(1)*ndof_node+1:el_idx_arr(2)*ndof_node);
            strain_sk_lo_unformatted = obj.element_strain_crossmod(el_idx_arr(3)*ndof_node+1:el_idx_arr(4)*ndof_node);
            
            % The 1D strain arrays of each spar/skin are reshaped into 2D
            % arrays, where each row is a DOF (there are 6) and each column
            % a new element.
            
            obj.crosssection.strain_fs = reshape(strain_fs_unformatted, [ndof_node, length(strain_fs_unformatted)/ndof_node]);
            obj.crosssection.strain_sk_up = reshape(strain_sk_up_unformatted, [ndof_node, length(strain_sk_up_unformatted)/ndof_node]);
            obj.crosssection.strain_rs = reshape(strain_rs_unformatted, [ndof_node, length(strain_rs_unformatted)/ndof_node]);
            obj.crosssection.strain_sk_lo = reshape(strain_sk_lo_unformatted, [ndof_node, length(strain_sk_lo_unformatted)/ndof_node]);
            
            
            % Von mises equivalent strain for shell elements according to
            %https://knowledge.autodesk.com/support/nastran/learn-explore/caas/CloudHelp/cloudhelp/2017/ENU/NSTRN-Reference/files/GUID-5DE89858-947A-4B35-B244-03BEE2CA779C-htm.html
            
            strain_von_mis_1 = obj.element_strain_crossmod(1:6:end).^2 + obj.element_strain_crossmod(2:6:end).^2;
            strain_von_mis_2 = obj.element_strain_crossmod(1:6:end) .* obj.element_strain_crossmod(2:6:end);
            strain_von_mis_3 = obj.element_strain_crossmod(3:6:end).^2;
            
            obj.element_strain_von_mis_crossmod = sqrt(4/9 * (strain_von_mis_1 - strain_von_mis_2) + 1/3 * strain_von_mis_3);
            
            
            % The stresses are obtained by taking the aforementioned
            % strains and multipliying them with the ABD stiffness matrices of their
            % respective skin/spar. Note that these are probably wrong, as
            % the last 3 DOFs are probably not the curvatures of the
            % element. They need to be revised.
            
            % Using 6DOF of strain and full ABD
            %             obj.crosssection.stress_fs = (obj.crosssection.fs_ABD * obj.crosssection.strain_fs)/obj.crosssection.t_sp_fr;
            %             obj.crosssection.stress_sk_up = (obj.crosssection.sk_up_ABD * obj.crosssection.strain_sk_up)/obj.crosssection.t_sk_up;
            %             obj.crosssection.stress_rs = (obj.crosssection.rs_ABD * obj.crosssection.strain_rs)/obj.crosssection.t_sp_re;
            %             obj.crosssection.stress_sk_lo = (obj.crosssection.sk_lo_ABD * obj.crosssection.strain_sk_lo)/obj.crosssection.t_sk_lo;
            
            % Using only first 3DOF of strain (in plane) and only A part of
            % ABD matrix
            obj.crosssection.stress_fs = (obj.crosssection.laminate_fs.ABD_stiff(1:3,1:3) * obj.crosssection.strain_fs(1:3,:))/obj.crosssection.t_sp_fr;
            obj.crosssection.stress_sk_up = (obj.crosssection.laminate_sk_up.ABD_stiff(1:3,1:3) * obj.crosssection.strain_sk_up(1:3,:))/obj.crosssection.t_sk_up;
            obj.crosssection.stress_rs = (obj.crosssection.laminate_rs.ABD_stiff(1:3,1:3) * obj.crosssection.strain_rs(1:3,:))/obj.crosssection.t_sp_re;
            obj.crosssection.stress_sk_lo = (obj.crosssection.laminate_sk_lo.ABD_stiff(1:3,1:3) * obj.crosssection.strain_sk_lo(1:3,:))/obj.crosssection.t_sk_lo;
            
            
            % The stresses are assembled into a 1D array with the same
            % structure as element_strain_crossmod
            obj.element_stress_crossmod = [obj.crosssection.stress_fs(:);
                                           obj.crosssection.stress_sk_up(:);
                                           obj.crosssection.stress_rs(:);
                                           obj.crosssection.stress_sk_lo(:)];
            
            
            
            % GENERAL VON MISES
            %             von_mis_11_22 = (obj.element_stress_crossmod(1:6:end)-obj.element_stress_crossmod(2:6:end)).^2;
            %             von_mis_22_33 = (obj.element_stress_crossmod(2:6:end)-obj.element_stress_crossmod(3:6:end)).^2;
            %             von_mis_11_33 = (obj.element_stress_crossmod(1:6:end)-obj.element_stress_crossmod(3:6:end)).^2;
            %             von_mis_shear = (obj.element_stress_crossmod(4:6:end).^2 + obj.element_stress_crossmod(5:6:end).^2 + obj.element_stress_crossmod(6:6:end).^2);
            
            
            % GENERAL PLANE STRESS VON MISES
            von_mis_11 = obj.element_stress_crossmod(1:3:end);
            von_mis_22 = obj.element_stress_crossmod(2:3:end);
            von_mis_12 = obj.element_stress_crossmod(3:3:end);
            
            % GENERAL VON MISES
            %             obj.element_stress_von_mis_crossmod = sqrt(0.5 * (von_mis_11_22 + von_mis_22_33 + von_mis_11_33 + 6 * von_mis_shear));
            
            % GENERAL PLANE STRESS VON MISES
            obj.element_stress_von_mis_crossmod = sqrt(von_mis_11.^2 - von_mis_11 .* von_mis_22 + von_mis_22.^2 + 3*von_mis_12.^2);
            
            obj.crosssection.stress_von_mis_fs = obj.element_stress_von_mis_crossmod(1:el_idx_arr(1))';
            obj.crosssection.stress_von_mis_rs = obj.element_stress_von_mis_crossmod(el_idx_arr(2)+1:el_idx_arr(3))';
            obj.crosssection.stress_von_mis_sk_up = obj.element_stress_von_mis_crossmod(el_idx_arr(1)+1:el_idx_arr(2))';
            obj.crosssection.stress_von_mis_sk_lo = obj.element_stress_von_mis_crossmod(el_idx_arr(3)+1:el_idx_arr(4))';
            
            % Calculating the fail indices of each skin/spar.
            %
            % SAFETY FACTOR IS INCLUDED:
            % If fail index > 0, safety factor is not met.
            % If fail index = 0, safety factor perfectly met.
            % If fail index < 0, safety factor exceeded.
            %
            % fail index = -[(design envelope safety factor)/(cross section user safety factor) - 1]
            if obj.crosssection.calcul_failIdx == 1
                obj.crosssection.laminate_sk_up = obj.crosssection.laminate_sk_up.calculate_safety_factor(obj.crosssection.strain_sk_up, obj.crosssection.safety_factor);
                obj.crosssection.laminate_sk_lo = obj.crosssection.laminate_sk_lo.calculate_safety_factor(obj.crosssection.strain_sk_lo, obj.crosssection.safety_factor);
                obj.crosssection.laminate_fs = obj.crosssection.laminate_fs.calculate_safety_factor(obj.crosssection.strain_fs, obj.crosssection.safety_factor);
                obj.crosssection.laminate_rs = obj.crosssection.laminate_rs.calculate_safety_factor(obj.crosssection.strain_rs, obj.crosssection.safety_factor);
            end
        end
        
        
        % This function generates a 3D plot of the strains or stresses
        % within the whole wingbox. The plot is similar to a color-mapped
        % FEM plot.
        function obj = f_plot_crossmod_data(obj, node_coords, next_node_coords, stress_flag, dof_interest)
            
            % used to print on the command window the y axis location of
            % the current beam element being plotted.
            display(node_coords(2));
            
            crs = obj.crosssection;
            
            % Depending on the flag value, the value arrays are attributed
            % different data.
            
            % Used to plot strains
            if stress_flag==0
                
                value_fs = crs.strain_fs(dof_interest,:);
                value_rs = crs.strain_rs(dof_interest,:);
                value_sk_up = crs.strain_sk_up(dof_interest,:);
                value_sk_lo = crs.strain_sk_lo(dof_interest,:);
                
                % Used to plot stresses
            elseif stress_flag == 1
                
                value_fs = crs.stress_fs(dof_interest,:);
                value_rs = crs.stress_rs(dof_interest,:);
                value_sk_up = crs.stress_sk_up(dof_interest,:);
                value_sk_lo = crs.stress_sk_lo(dof_interest,:);
                
                % Used to plot von mises stresses
            elseif stress_flag == 2
                
                value_fs = crs.stress_von_mis_fs;
                value_rs = crs.stress_von_mis_rs;
                value_sk_up = crs.stress_von_mis_sk_up;
                value_sk_lo = crs.stress_von_mis_sk_lo;
                
            end
            
            % The surfaces of each spar/skin is created, composed of a
            % mesh whose elements correspond to the shell element used in
            % the discretized cross section. Each element is assigned a
            % color depending on the value of the data assigned to it.
            % Each block of code below corresponds to the front spar, then
            % the rear spar, then the upper skin, and finally the lower
            % skin.
            [Y,Z] = meshgrid([node_coords(2),next_node_coords(2)],[crs.fs_nodes_yz(:,2)+node_coords(3);crs.fs_nodes_yz(:,2)+next_node_coords(3)]);
            X =[node_coords(1)+crs.fs_nodes_yz(:,1),next_node_coords(1)+crs.fs_nodes_yz(:,1);node_coords(1)+crs.fs_nodes_yz(:,1),next_node_coords(1)+crs.fs_nodes_yz(:,1)];
            value_fs_formatted = [[value_fs(1)';value_fs'],[value_fs(1)';value_fs'];[value_fs(1)';value_fs'],[value_fs(1)';value_fs']];
            surf(X,Y,Z,value_fs_formatted,'EdgeColor','none','LineStyle','none');
            
            [Y,Z] = meshgrid([node_coords(2),next_node_coords(2)],[crs.rs_nodes_yz(:,2)+node_coords(3);crs.rs_nodes_yz(:,2)+next_node_coords(3)]);
            X =[node_coords(1)+crs.rs_nodes_yz(:,1),next_node_coords(1)+crs.rs_nodes_yz(:,1);node_coords(1)+crs.rs_nodes_yz(:,1),next_node_coords(1)+crs.rs_nodes_yz(:,1)];
            value_rs_formatted = [[value_rs(1)';value_rs'],[value_rs(1)';value_rs'];[value_rs(1)';value_rs'],[value_rs(1)';value_rs']];
            surf(X,Y,Z,value_rs_formatted,'EdgeColor','none','LineStyle','none');
            
            [Y,X] = meshgrid([node_coords(2),next_node_coords(2)],[crs.sk_up_nodes_yz(:,1)+node_coords(1);crs.sk_up_nodes_yz(:,1)+next_node_coords(1)]);
            Z =[node_coords(3)+crs.sk_up_nodes_yz(:,2),next_node_coords(3)+crs.sk_up_nodes_yz(:,2);node_coords(3)+crs.sk_up_nodes_yz(:,2),next_node_coords(3)+crs.sk_up_nodes_yz(:,2)];
            value_sk_up_formatted = [[value_sk_up(1)';value_sk_up'],[value_sk_up(1)';value_sk_up'];[value_sk_up(1)';value_sk_up'],[value_sk_up(1)';value_sk_up']];
            surf(X,Y,Z,value_sk_up_formatted,'EdgeColor','none','LineStyle','none');
            
            [Y,X] = meshgrid([node_coords(2),next_node_coords(2)],[crs.sk_lo_nodes_yz(:,1)+node_coords(1);crs.sk_lo_nodes_yz(:,1)+next_node_coords(1)]);
            Z =[node_coords(3)+crs.sk_lo_nodes_yz(:,2),next_node_coords(3)+crs.sk_lo_nodes_yz(:,2);node_coords(3)+crs.sk_lo_nodes_yz(:,2),next_node_coords(3)+crs.sk_lo_nodes_yz(:,2)];
            value_sk_lo_formatted = [[value_sk_lo(1)';value_sk_lo'],[value_sk_lo(1)';value_sk_lo'];[value_sk_lo(1)';value_sk_lo'],[value_sk_lo(1)';value_sk_lo']];
            surf(X,Y,Z,value_sk_lo_formatted,'EdgeColor','none','LineStyle','none');
        end
        
    end
end