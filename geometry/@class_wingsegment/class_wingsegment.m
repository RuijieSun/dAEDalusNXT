%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%
classdef class_wingsegment
    %class wingsegment geometric definitions for a wing segment
    %   Detailed explanation goes here
    
    properties
        %> reference position (c4 position)
        pos;
        %> is this surface symmetric?
        symmetric;
        %> grid start index of this segment in the global grid
        grid_start_idx;
        %> panel start index of this segment in the global grid
        panel_start_idx;
        %> quarter chordline coordinates
        c4_coords;
        %> spanwise forces
        span_forces;
        %> spanwise moments
        span_moments;
        %> special case
        real_sweep=[];
        %> quarter chordline sweep
        c4_sweep;
        %> leading edge sweep
        le_sweep;
        %> taper ratio
        TR;
        %> span
        b;
        %> root chord
        c_r;
        %> tip chord
        c_t;
        %> root profile data points
        profile_r;
        %> name of root profile
        profile_name_r;
        %> tip profile data points
        profile_t;
        %> name of root profile
        profile_name_t;
        %> root twist
        Theta_r;
        %> tip twist
        Theta_t;
        %> dihedral at c4 line
        dihed;
        %> dihedral at leading edge line
        le_dihed;
        
        %> thickness ratio
        tc=0.2;
        %> panel edge coordinates
        xyz;
        %> part of wing fixed
        xyz_fixed;
        %> edge coordinates of trailing edge device
        xyz_te_device;
        %> edge coordinates of leading edge device
        xyz_le_device;
        %> edge coordinates of spoiler device
        xyz_sp_device;
        %> segment area
        S;
        %> segment wetted area
        S_wet;
        %> mean aerodynamic chord
        c_mac;
        
        is_laminar;
        
        skeleton_line_r;
        skeleton_line_t;
        
        has_te_cs=0; % information flag if there is a trailing edge control surface
        has_le_cs=0; % information flag if there is a leading edge control surface
        has_sp_cs=0; % information flag if there is a spoiler control surface
        
        le_device; % class containing the information about the leading edge control surface
        te_device; % class containing the information about the leading edge control surface
        sp_device; % class containing the information about the spoiler control surface
        
        c_le_device;
        c_te_device;
        c_sp_device;
        
        delta_te_device;
        delta_le_device;
        delta_sp_device;
        
        te_max_defl;
        te_min_defl;
        
        n_span;
        n_chord;
        n_le_panels;
        n_te_panels;
        n_sp_panels;
        n_ctr1_panels;
        n_ctr2_panels;
        n_te_panels_overlap;
        n_te_panels_free;
        
        D_f;
        CD_f;
        
        grid;
        grid_flat;
        %for potential flow
        is_te;
        grid_vol_upper;
        grid_vol_lower;
        % grid for wake
        grid_wake;
        % panels for wake
        panels_wake;
        
        
        
        panels;
        nxt_ds_pt;
        te_idx;
        
        grid_le;
        panels_le;
        grid_te;
        panels_te;
        
        wingbox_coords;
        wingbox_c4;
        wingbox_height;
        
        %> number of beam elements; will be determined by grid spacing if
        % zero
        nBeamelements = 0;
        
        % This struct will hold the following information about this wing
        % segment's structure:
        % fs_segments: [xsi_root, xsi_tip,  t_web, t_top, t_bottom]
        % rs_segments: [xsi_root, xsi_tip,  t_web, t_top, t_bottom]
        % where xsi is the normalized chordwise position of the spar.
        structural_properties;
    end
    
    methods (Static)
        function obj = create_from_cpacs(tixi, tigl, wingIndex, segmentIndex)
            obj = class_wingsegment();
            
            xr = zeros(3,1);
            xt = zeros(3,1);
            xle_r = zeros(3,1);
            xle_t = zeros(3,1);
            xte_r = zeros(3,1);
            xte_t = zeros(3,1);

            % Get the root and tip quarter chord points in global coordinates.
            [xr(1), xr(2), xr(3)] = tiglWingGetChordPoint(tigl, wingIndex, segmentIndex, 0, .25);
            [xt(1), xt(2), xt(3)] = tiglWingGetChordPoint(tigl, wingIndex, segmentIndex, 1, .25);

            % The dihedral angle of the quarter chord line can be found by taking the 
            % inverse tangent of the length along the global z-axis over the length 
            % along the global y-axis.
            Gamma = atand((xt(3)-xr(3))/(xt(2)-xr(2)));
            
            % Calculate the rotation matrix for the diherdral angle
            Rx = [1 0 0; 0 cosd(Gamma) -sind(Gamma); 0 sind(Gamma) cosd(Gamma)];

            % De-rotate the quarter chord line from the dihedral down onto the global
            % x-y-plane. Make sure to subtract the root coordinates first, because the
            % rotation is done w.r.t. the root quarter point. Then at the root
            % coordinates again.
            xt_ = Rx\(xt-xr) + xr;

            % The quarter chord sweep is then the angle between the global y-axis and
            % this de-rotated line.
            Lambda = atand((xt_(1)-xr(1))/(xt_(2)-xr(2)));

            % In dAEDalus the span is defined as the length between the two quarter
            % chord points.
            b = norm(xt - xr);

            % Now the span needs to be corrected for dAEDalus' 'fake' sweep.
            b = b*cosd(Lambda);

            % Get the root and tip chord lengths by simply calculating the distance
            % between the leading and trailing edge points of the root and tip in
            % global coordinates.
            [xle_r(1), xle_r(2), xle_r(3)] = tiglWingGetChordPoint(tigl, wingIndex, segmentIndex, 0, 0);
            [xte_r(1), xte_r(2), xte_r(3)] = tiglWingGetChordPoint(tigl, wingIndex, segmentIndex, 0, 1);
            [xle_t(1), xle_t(2), xle_t(3)] = tiglWingGetChordPoint(tigl, wingIndex, segmentIndex, 1, 0);
            [xte_t(1), xte_t(2), xte_t(3)] = tiglWingGetChordPoint(tigl, wingIndex, segmentIndex, 1, 1);

            c_r = norm(xte_r - xle_r);
            c_t = norm(xte_t - xle_t);

            % To calculate the twist angles de-rotate the leading edge points of the
            % root and tip sections first. Again, subtract and then add the quarter
            % point first to assure the rotation happens relative to the quarter chord
            % points.
            xle_r_ = Rx\(xle_r - xr) + xr;
            xle_t_ = Rx\(xle_t - xr) + xr;

            % Then the twists can be determined by calculating the angles between the
            % local x-axes and the line connecting the leading edges and quarter chord
            % points.
            twist_r = atand((xle_r_(3) - xr(3))/abs(xle_r_(1) - xr(1)));
            twist_t = atand((xle_t_(3) - xt_(3))/abs(xle_t_(1) - xt_(1)));
            
            % Set the object's properties from these calculated properties
            obj.pos = xr';
            obj.b = b;
            obj.c_r = c_r;
            obj.c_t = c_t;
            obj.dihed = Gamma;
            obj.c4_sweep = Lambda;
            obj.Theta_r = twist_r;
            obj.Theta_t = twist_t;
            obj.has_te_cs = 0;
            obj.has_le_cs = 0;
            %obj.c4_sweep=atan2d(obj.b*tand(obj.le_sweep)+obj.c_t/4-obj.c_r/4, obj.b);
            obj=obj.complete_params_from_stdinit();
                
            [secUID_r, elemUID_r] = tiglWingGetInnerSectionAndElementUID(tigl, wingIndex, segmentIndex);        
            [secUID_t, elemUID_t] = tiglWingGetOuterSectionAndElementUID(tigl, wingIndex, segmentIndex);
            
            x_sec_r = tixiUIDGetXPath(tixi, secUID_r);
            x_sec_t = tixiUIDGetXPath(tixi, secUID_t);       
            
            [c_r, ~, t_r] = tixiGetPoint(tixi, [x_sec_r, '/transformation/scaling']);
            [c_t, ~, t_t] = tixiGetPoint(tixi, [x_sec_t, '/transformation/scaling']);
                   
            x_elem_r = tixiUIDGetXPath(tixi, elemUID_r);
            x_elem_t = tixiUIDGetXPath(tixi, elemUID_t);
            
            obj.profile_name_r = tixiGetTextElement(tixi, [x_elem_r, '/airfoilUID']);
            obj.profile_name_t = tixiGetTextElement(tixi, [x_elem_t, '/airfoilUID']);
                                         
            x_af_r = tixiUIDGetXPath(tixi, obj.profile_name_r);
            x_af_t = tixiUIDGetXPath(tixi, obj.profile_name_t);
            
            try
                n_af_r = tixiGetVectorSize(tixi, [x_af_r, '/pointList/x']);
                n_af_t = tixiGetVectorSize(tixi, [x_af_t, '/pointList/x']);

                af_r_x = tixiGetFloatVector(tixi, [x_af_r, '/pointList/x'], n_af_r);
                af_r_z = tixiGetFloatVector(tixi, [x_af_r, '/pointList/z'], n_af_r);

                af_t_x = tixiGetFloatVector(tixi, [x_af_t, '/pointList/x'], n_af_t);
                af_t_z = tixiGetFloatVector(tixi, [x_af_t, '/pointList/z'], n_af_t);

                if mod(n_af_r, 2) ~= 0
                    ids = [1:ceil(n_af_r/2), ceil(n_af_r/2):n_af_r];
                    af_r_x = af_r_x(ids);
                    af_r_z = af_r_z(ids);
                end
                
                if mod(n_af_t, 2) ~= 0
                    ids = [1:ceil(n_af_t/2), ceil(n_af_t/2):n_af_t];
                    af_t_x = af_t_x(ids);
                    af_t_z = af_t_z(ids); 
                end
                
                n_af_r = ceil(n_af_r/2);
                n_af_t = ceil(n_af_t/2);
                
                af_r_z_upper = af_r_z(1:n_af_r)';
                af_r_z_lower = af_r_z(n_af_r+1:end)';
                
                af_t_z_upper = af_t_z(1:n_af_t)';
                af_t_z_lower = af_t_z(n_af_t+1:end)';
                
                t_max_r = max(abs(af_r_z_upper - af_r_z_lower));
                t_max_t = max(abs(af_t_z_upper - af_t_z_lower));

                obj.tc = mean([t_max_r*t_r/c_r, t_max_t*t_t/c_t]);

                obj.profile_r = [n_af_r, n_af_r; ...
                    flipud([af_r_x(1:n_af_r)', af_r_z_upper]); ...
                    af_r_x(n_af_r+1:end)', af_r_z_lower];

                obj.profile_t = [n_af_t, n_af_t; ...
                    flipud([af_t_x(1:n_af_t)', af_t_z_upper]); ...
                    af_t_x(n_af_t+1:end)', af_t_z_lower];               
            catch
            end

            obj.xyz = [xle_r, xle_t, xte_t, xte_r];
            obj.le_dihed=asind((xle_t(3)-xle_r(3))/obj.b);
            
            obj=obj.compute_controlsurface_coordinates();
            obj=obj.compute_xyz_fixed();
        end
    end
    
    methods
        
        function obj = class_wingsegment(varargin)
            
            % input parameter is edge coordinates
            if nargin == 0
                % pass
            elseif nargin==1
                obj=obj.read_xml_definition(varargin{1});
            elseif nargin==4
                obj.profile_name_r='flat_plate';
                obj.profile_name_t='flat_plate';
                obj=obj.load_airfoils();
                sz=size(varargin{1});
                if sz(1)==3 && sz(2)==4
                    obj.xyz=varargin{1};
                end
                
                obj=obj.compute_parameters_from_coords();
                obj.symmetric=varargin{2};
                obj.profile_r=varargin{3};
                obj.profile_t=varargin{4};
                % input parameter is edge coordinates and control surfaces
            elseif nargin==6
                obj.profile_name_r='flat_plate';
                obj.profile_name_t='flat_plate';
                obj=obj.load_airfoils();
                sz=size(varargin{1});
                if sz(1)==3 && sz(2)==4
                    obj.xyz=varargin{1};
                end
                
                obj=obj.compute_parameters_from_coords();
                obj.symmetric=varargin{2};
                obj.profile_r=varargin{5};
                obj.profile_t=varargin{6};
                
                if strcmp(varargin{3},'cs_te')
                    cs_te=varargin{4};
                    for i=1:length(cs_te)
                        if i<length(cs_te)
                            obj.c_te_device(i)=obj.c_r*(cs_te(i)-cs_te(i+1));
                        else
                            obj.c_te_device(i)=obj.c_r*cs_te(i);
                        end
                        obj.delta_te_device(i)=25*i;
                    end
                    obj.has_te_cs=1;
                elseif strcmp(varargin{3},'cs_le')
                    cs_le=varargin{4};
                    for i=1:length(cs_le)
                        if i<length(cs_le)
                            obj.c_le_device(i)=obj.c_r*(cs_le(i)-cs_le(i+1));
                        else
                            obj.c_le_device(i)=obj.c_r*cs_le(i);
                        end
                        obj.delta_le_device(i)=-10*i;
                    end
                    obj.has_le_cs=1;
                end
                
                % standard initialization
            else
                obj.pos=varargin{1};
                obj.symmetric=varargin{2};
                obj.dihed=varargin{3};
                
                LambdaSpec=varargin{4};
                if strcmp(LambdaSpec,'c4-sweep')
                    obj.c4_sweep=varargin{5};
                elseif  strcmp(LambdaSpec,'real-sweep')  
                    obj.real_sweep=varargin{5};
                elseif strcmp(LambdaSpec,'le-sweep')
                    obj.le_sweep=varargin{5};
                end
                
                if nargin==5
                    obj.profile_name_r='flat_plate';
                    obj.profile_name_t='flat_plate';
                    obj=obj.load_airfoils();
                    obj.Theta_r=0;
                    obj.Theta_t=0;
                    obj.b=5;
                    obj.TR=0.5;
                    obj.c_r=3;
                    obj.S=obj.b/2*(obj.c_r*(1+obj.TR));
                    obj.S_wet=2*obj.S;
                    obj.c_t=obj.TR*obj.c_r;
                    obj.c_mac=8*obj.S/((1+obj.TR)^2*obj.b^2)*(obj.b/2-(1-obj.TR)*obj.b/2+(1-obj.TR)^2*obj.b/6);
                    
                    if strcmp(LambdaSpec,'c4-sweep')
                        obj.le_sweep=atan((obj.b*tan(obj.c4_sweep*pi/180)+obj.c_r/4-obj.c_t/4)/obj.b)*180/pi;
                    elseif strcmp(LambdaSpec,'le-sweep')
                        obj.c4_sweep=atan((obj.b*tan(obj.le_sweep*pi/180)+obj.c_t/4-obj.c_r/4)/obj.b)*180/pi;
                    elseif strcmp(LambdaSpec,'real-sweep')
                        obj.le_sweep=[];
                        obj.c4_sweep=[];
                    end
                elseif nargin==19
                    property1=varargin{6};
                    property2=varargin{8};
                    property3=varargin{10};
                    property4=varargin{12};
                    property5=varargin{14};
                    property6=varargin{16};
                    property7=varargin{18};
                    
                    if(strcmp(property1,'b'))&&(strcmp(property2,'c_r'))&&(strcmp(property3,'c_t'))&&(strcmp(property4,'Theta_r'))&&(strcmp(property5,'Theta_t'))&&(strcmp(property6,'Profile_r'))&&(strcmp(property7,'Profile_t'))
                        
                        obj.profile_name_r=varargin{17};
                        obj.profile_name_t=varargin{19};
                        obj=obj.load_airfoils();
                        
                        obj.b=varargin{7};
                        obj.c_r=varargin{9};
                        obj.c_t=varargin{11};
                        obj.Theta_r=varargin{13};
                        obj.Theta_t=varargin{15};
                        
                        obj=obj.complete_params_from_stdinit();
                        obj.has_le_cs=0;
                        obj.has_te_cs=0;
                    end
                    
                elseif nargin==21
                    property1=varargin{6};
                    property2=varargin{8};
                    property3=varargin{10};
                    property4=varargin{12};
                    property5=varargin{14};
                    property6=varargin{16};
                    property7=varargin{18};
                    property8=varargin{20};
                    
                    if(strcmp(property1,'b'))&&(strcmp(property2,'c_r'))&&(strcmp(property3,'c_t'))&&(strcmp(property4,'Theta_r'))&&(strcmp(property5,'Theta_t'))&&(strcmp(property6,'Profile_r'))&&(strcmp(property7,'Profile_t'))&&(strcmp(property8,'cs_te'))
                        obj.profile_name_r=varargin{17};
                        obj.profile_name_t=varargin{19};
                        
                        obj=obj.load_airfoils();
                        
                        obj.b=varargin{7};
                        obj.c_r=varargin{9};
                        obj.c_t=varargin{11};
                        obj.Theta_r=varargin{13};
                        obj.Theta_t=varargin{15};
                        cs=varargin{21};
                        
                        obj=obj.complete_params_from_stdinit();
                        
                        for i=1:length(cs)
                            if i<length(cs)
                                obj.c_te_device(i)=obj.c_r*(cs(i)-cs(i+1));
                                
                            else
                                obj.c_te_device(i)=obj.c_r*cs(i);
                            end
                            obj.delta_te_device(i)=-10*i;
                        end
                        
                        obj.has_te_cs=1;
                    end
                    
                elseif nargin==23
                    property1=varargin{6};
                    property2=varargin{8};
                    property3=varargin{10};
                    property4=varargin{12};
                    property5=varargin{14};
                    property6=varargin{16};
                    property7=varargin{18};
                    property8=varargin{20};
                    property9=varargin{22};
                    if(strcmp(property1,'b'))&&(strcmp(property2,'c_r'))&&(strcmp(property3,'c_t'))&&(strcmp(property4,'Theta_r'))&&(strcmp(property5,'Theta_t'))&&(strcmp(property6,'Profile_r'))&&(strcmp(property7,'Profile_t'))&&(strcmp(property8,'cs_te'))&&(strcmp(property9,'cs_le'))
                        
                        obj.profile_name_r=varargin{17};
                        obj.profile_name_t=varargin{19};
                        obj=obj.load_airfoils();
                        
                        obj.b=varargin{7};
                        obj.c_r=varargin{9};
                        obj.c_t=varargin{11};
                        obj.Theta_r=varargin{13};
                        obj.Theta_t=varargin{15};
                        cs_te=varargin{21};
                        cs_le=varargin{23};
                        
                        obj=obj.complete_params_from_stdinit();
                        
                        for i=1:length(cs_te)
                            if i<length(cs_te)
                                obj.c_te_device(i)=obj.c_r*(cs_te(i)-cs_te(i+1));
                            else
                                obj.c_te_device(i)=obj.c_r*cs_te(i);
                            end
                            obj.delta_te_device(i)=-10*i;
                        end
                        
                        for i=1:length(cs_le)
                            if i==1
                                obj.c_le_device(i)=obj.c_r*cs_le(i);
                            else
                                obj.c_le_device(i)=obj.c_r*(cs_le(i)-cs_le(i-1));
                            end
                            obj.delta_le_device(i)=0;
                        end
                        obj.has_te_cs=1;
                        obj.has_le_cs=1;
                    elseif(strcmp(property1,'b'))&&(strcmp(property2,'c_r'))&&(strcmp(property3,'c_t'))&&(strcmp(property4,'Theta_r'))&&(strcmp(property5,'Theta_t'))&&(strcmp(property6,'cs_te'))&&(strcmp(property7,'delta_te'))
                        
                        obj.profile_name_r=varargin{17};
                        obj.profile_name_t=varargin{19};
                        obj=obj.load_airfoils();
                        
                        obj.b=varargin{7};
                        obj.c_r=varargin{9};
                        obj.c_t=varargin{11};
                        obj.Theta_r=varargin{13};
                        obj.Theta_t=varargin{15};
                        cs_te=varargin{21};
                        delta_te=varargin{23};
                        
                        obj=obj.complete_params_from_stdinit();
                        
                        for i=1:length(cs_te)
                            if i<length(cs_te)
                                obj.c_te_device(i)=obj.c_r*(cs_te(i)-cs_te(i+1));
                                
                            else
                                obj.c_te_device(i)=obj.c_r*cs_te(i);
                            end
                            obj.delta_te_device(i)=0;
                        end
                        
                        for i=1:length(delta_te)
                            obj.delta_te_device(i)=delta_te(i);
                        end
                        obj.has_te_cs=1;
                    end
                end
                obj=obj.compute_segment_coordinates();
            end
            
            if nargin ~= 0
                % compute flap edge points
                obj=obj.compute_controlsurface_coordinates();
                obj=obj.compute_xyz_fixed();
            end
        end
        
        function obj=add_control_surface(obj,control_surface,varargin)
            
            % In case the control surfaces are untapered, the chord of the
            % CS should be constant throughout the whole CS, meaning that
            % the "hinge" property of the CS is multiplied with the root
            % chord of the parent segment. In case of a CS spread over
            % multiple segments, the children would lose the root chord
            % info of the parent, and their chord would be wrongly
            % sized using the root chords of the children, rather than the parent's. 
            % This is fixed by thenargin=3 case below, where varargin{1} is the root chord of
            % the segment where the original CS surface starts. The hinge
            % of the CS is later multiplied with the hinge_mult, so that
            % when multiplying the hinge of the child with the root chord
            % of the child, the chord of the parent CS is obtained.
            
            % Obviously if the CS is tapered, none of this is applied, and
            % the hinge_mult is 1.
            if nargin == 3
                hinge_mult =  varargin{1}/obj.c_r;
            elseif nargin==2
                hinge_mult =  1;
            end
            
            switch(control_surface.pos)
                case 0
                    obj.te_device=control_surface;
                    obj.te_device.hinge = obj.te_device.hinge*hinge_mult;
                    cs_te=obj.te_device.hinge;
                    
                    for i=1:length(cs_te)
                        if i<length(cs_te)
                            obj.c_te_device(i)=obj.c_r*(cs_te(i)-cs_te(i+1));
                        else
                            obj.c_te_device(i)=obj.c_r*cs_te(i);
                        end
                        obj.delta_te_device(i)=0;
                    end
                    obj.has_te_cs=1;
                    
                case 1
                    obj.le_device=control_surface;
                    obj.le_device.hinge = obj.le_device.hinge*hinge_mult;
                    cs_le=obj.le_device.hinge;
                    for i=1:length(cs_le)
                        if i==1
                            obj.c_le_device(i)=obj.c_r*cs_le(i);
                        else
                            obj.c_le_device(i)=obj.c_r*(cs_le(i)-cs_le(i-1));
                        end
                        obj.delta_le_device(i)=0;
                    end
                    obj.has_le_cs=1;
                case 2
                    obj.sp_device=control_surface;
                    obj.sp_device.hinge = obj.sp_device.hinge*hinge_mult;
                    cs_sp=obj.sp_device.hinge;
                    obj.c_sp_device=obj.c_r*cs_sp;
                    obj.delta_sp_device=0;
                    obj.has_sp_cs=1;
            end
            
            obj=obj.compute_controlsurface_coordinates();
            obj=obj.compute_xyz_fixed();
        end
        
        
        function obj= complete_params_from_stdinit(obj)
            %                                     if strcmp(LambdaSpec,'c4-sweep')
            %                             obj.le_sweep=atan((obj.b*tan(obj.c4_sweep*pi/180)+obj.c_r/4-obj.c_t/4)/obj.b)*180/pi;
            %                         elseif strcmp(LambdaSpec,'le-sweep')
            %                             obj.c4_sweep=atan((obj.b*tan(obj.le_sweep*pi/180)+obj.c_t/4-obj.c_r/4)/obj.b)*180/pi;
            %                         end
            if isempty(obj.le_sweep)
                obj.le_sweep=atan((obj.b*tan(obj.c4_sweep*pi/180)+obj.c_r/4-obj.c_t/4)/obj.b)*180/pi;
            elseif isempty(obj.c4_sweep)
                obj.c4_sweep=atan((obj.b*tan(obj.le_sweep*pi/180)+obj.c_t/4-obj.c_r/4)/obj.b)*180/pi;
            end
            
            obj.TR=obj.c_t/obj.c_r;
            obj.S=obj.b/2*(obj.c_r*(1+obj.TR));
            obj.S_wet=2*obj.S;
            obj.c_mac =obj.c_r * 2/3 * (( 1 + obj.TR +obj.TR^2 )/ ( 1 + obj.TR ));
        end
        
        function obj = compute_friction_drag(obj,state,S_ref,varargin)
            obj.S_wet=2*obj.S;
            Re_l=state.rho_air*norm(state.V_inf)*obj.c_mac/state.mu;
            if nargin==3
                % laminar or turbulent
                if Re_l>5E5
                    cf=0.455/log10(Re_l)^2.58;
                else
                    cf=1.328/sqrt(Re_l);
                end
                
                if obj.is_laminar==1
                    cf=1.328/sqrt(Re_l);
                end
            end
            % controlsurface gap drag
            % http://adg.stanford.edu/aa241/drag/gapdrag.html
            CD_gaps=0;
            if obj.has_te_cs
                CD_gaps = 0.0002*cos(obj.c4_sweep*pi/180)^2*obj.S/S_ref;
            end
            % form factor formula
            % (adg.stanford.edu/aa241/drag/lsformfactor.html)
            if state.Ma*cos(obj.le_sweep*pi/180)>1
                C=0;
            else
                C=1.1;
            end
            
            cossw2=cos(obj.le_sweep*pi/180)^2;
            Ma2=state.Ma^2;
            form_factor=1+2*C*obj.tc*cossw2/sqrt(1-Ma2*cossw2)+C^2*cossw2*obj.tc^2*(1+5*cossw2)/(2*(1-Ma2*cossw2));
            
            obj.CD_f=form_factor*cf*obj.S_wet/S_ref+CD_gaps;
            obj.D_f=form_factor*1/2*state.rho_air*norm(state.V_inf)^2*cf*obj.S_wet;
        end
        
        function obj=load_airfoils(obj)
            airfoil_name_r=obj.profile_name_r;
            airfoil_name_t=obj.profile_name_t;
            if strcmp(airfoil_name_r,'flat_plate')
                obj.profile_name_r='flat_plate';
                obj.profile_r=[3 3;0 0;0.5 0;1 0;0 0;0.5 0;1 0];
            else
                obj.profile_name_r=['airfoil/',airfoil_name_r,'.DAT'];
                obj.profile_r=load(['airfoil/',airfoil_name_r,'.DAT']);
            end
            
            if strcmp(airfoil_name_r,'flat_plate')
                obj.profile_t=[3 3;0 0;0.5 0;1 0;0 0;0.5 0;1 0];
                obj.profile_name_t='flat_plate';
            else
                obj.profile_t=load(['airfoil/',airfoil_name_t,'.DAT']);
                obj.profile_name_t=['airfoil/',airfoil_name_t,'.DAT'];
            end
        end
        
        function obj=compute_parameters_from_coords(obj)
            dv=obj.xyz(2:3,2)-obj.xyz(2:3,1);
            obj.b=hypot(dv(1),dv(2));
            obj.c_r=sqrt(sum((obj.xyz(:,4)-obj.xyz(:,1)).^2));
            obj.c_t=sqrt(sum((obj.xyz(:,3)-obj.xyz(:,2)).^2));
            obj.S=(obj.c_t+obj.c_r)/2*obj.b;
            
            obj.pos=[obj.xyz(:,1)+0.25*(obj.xyz(:,4)-obj.xyz(:,1))]';
            pos_t=obj.xyz(:,2)+0.25*(obj.xyz(:,3)-obj.xyz(:,2));
            
            dx=pos_t(1)-obj.pos(1);
            dy=pos_t(2)-obj.pos(2);
            dz=pos_t(3)-obj.pos(3);
            obj.c4_sweep=atand(dx/obj.b);
            if dy>1e-9
                obj.dihed=asind(dz/obj.b);
            elseif dy<-1e-9
                obj.dihed=180-asind(dz/obj.b);
            else
                obj.dihed=90;
            end
            
            obj.le_sweep=atand((obj.b*tan(obj.c4_sweep*pi/180)+obj.c_r/4-obj.c_t/4)/obj.b);
            
            dz=obj.xyz(3,4)-obj.xyz(3,1);
            dy=obj.xyz(2,4)-obj.xyz(2,1);
            if obj.dihed==90
                dn=-sqrt(dz^2+dy^2)*sign(dy);
            else
                dn=sqrt(dz^2+dy^2)*sign(dz);
            end
            obj.Theta_r=-asind(dn/obj.c_r);
            
            dz=obj.xyz(3,3)-obj.xyz(3,2);
            dy=obj.xyz(2,3)-obj.xyz(2,2);
            
            if obj.dihed==90
                dn=-sqrt(dz^2+dy^2)*sign(dy);
            else
                dn=sqrt(dz^2+dy^2)*sign(dz);
            end
            
            obj.le_dihed=asind((obj.xyz(3,2)-obj.xyz(3,1))/obj.b);
            
            obj.Theta_t=-asind(dn/obj.c_t);
            obj.TR=obj.c_t/obj.c_r;
            
            %obj.c_mac=8*obj.S/((1+obj.TR)^2*obj.b^2)*(obj.b/2-(1-obj.TR)*obj.b/2+(1-obj.TR)^2*obj.b/6);
            obj.c_mac =obj.c_r * 2/3 * (( 1 + obj.TR +obj.TR^2 )/ ( 1 + obj.TR ));
            %c_mac= obj.c_r-(2*(obj.c_r-obj.c_t)*(0.5*obj.c_r+obj.c_t)/(3*(obj.c_r+obj.c_t)))
            %obj.c_mac=c_mac;
        end
        
        function obj=set_cs_deflections(obj,varargin)
            if nargin==3
                if strcmp(varargin{1},'te')
                    obj.delta_te_device=varargin{2};
                elseif strcmp(varargin{1},'le')
                    obj.delta_le_device=varargin{2};
                end
            elseif nargin==5
                if strcmp(varargin{1},'te') && strcmp(varargin{3},'le')
                    obj.delta_te_device=varargin{2};
                    obj.delta_le_device=varargin{4};
                elseif strcmp(varargin{3},'le') && strcmp(varargin{1},'te')
                    obj.delta_te_device=varargin{4};
                    obj.delta_le_device=varargin{2};
                end
            else
                
            end
            obj=obj.compute_controlsurface_coordinates();
        end
        
        function obj=f_deflect_control_surface(obj,name,deflection,varargin)
            
            side=[];
            
            if nargin==4
                side=varargin{1};
            end
            
            
            if  ~isempty(obj.te_device)
                if strcmp(name,obj.te_device.name)
                    if ~isempty(side)
                        if strcmp(side,'left')
                            obj.te_device.delta=deflection;
                            obj.te_device.delta_l_r(1)=deflection;
                        elseif strcmp(side,'right')
                            obj.te_device.delta_l_r(2)=deflection;
                        end
                    end
                    obj.te_device.delta=deflection;
                end
            end
            
            
            if  ~isempty(obj.le_device)
                if strcmp(name,obj.le_device.name)
                    obj.le_device.delta=deflection;
                end
            end
            
            
            if  ~isempty(obj.sp_device)
                if strcmp(name,obj.sp_device.name)
                    obj.sp_device.delta=deflection;
                end
            end
            
            
            obj=obj.compute_controlsurface_coordinates();
            
        end
        
        function obj=deflect_control_surface(obj,name,delta)
            if strcmp(obj.te_device.name,name)
                %% TODO: check if right format!
                obj.delta_te_cs=delta;
            elseif strcmp(obj.le_device.name,name)
                obj.delta_le_cs=delta;
            end
            obj=obj.compute_controlsurface_coordinates();
        end
        
        
        function obj=compute_segment_coordinates(obj)
            
            if isempty(obj.real_sweep)
                %the following code does not work properly, e.g. input span
                %10m, sweep 50, dihed 20 -> real span = 15m
            % compute segment edge points
            p1=obj.pos+[-obj.c_r/4*cos(obj.Theta_r*pi/180), -obj.c_r/4*sin(obj.Theta_r*pi/180)*sin(obj.dihed*pi/180), +obj.c_r/4*sin(obj.Theta_r*pi/180)*cos(obj.dihed*pi/180)];
            p2=obj.pos+[+obj.b*tan(obj.c4_sweep*pi/180), obj.b*cos(obj.dihed*pi/180), +obj.b*sin(obj.dihed*pi/180)]+[-obj.c_t/4*cos(obj.Theta_t*pi/180) -obj.c_t/4*sin(obj.Theta_t*pi/180)*sin(obj.dihed*pi/180) +obj.c_t/4*sin(obj.Theta_t*pi/180)*cos(obj.dihed*pi/180)];
            p3=obj.pos+[+obj.b*tan(obj.c4_sweep*pi/180), obj.b*cos(obj.dihed*pi/180), +obj.b*sin(obj.dihed*pi/180)]-[-obj.c_t*3/4*cos(obj.Theta_t*pi/180) -obj.c_t*3/4*sin(obj.Theta_t*pi/180)*sin(obj.dihed*pi/180) +obj.c_t*3/4*sin(obj.Theta_t*pi/180)*cos(obj.dihed*pi/180)];
            p4=obj.pos-[-obj.c_r*3/4*cos(obj.Theta_r*pi/180), -obj.c_r*3/4*sin(obj.Theta_r*pi/180)*sin(obj.dihed*pi/180), +obj.c_r*3/4*sin(obj.Theta_r*pi/180)*cos(obj.dihed*pi/180)];
            
            obj.le_dihed=asind((p2(3)-p1(3))/obj.b);
            obj.xyz=[p1' p2' p3' p4'];
            else
                %the following code does not work and is commented out
%                 p1=[0 0 0];
%                 p2=[0 obj.b 0];
%                 p3=[obj.c_r obj.b 0];
%                 p4=[obj.c_r 0 0];
%                 a=obj.dihed;
%                 b=obj.Theta_r;
%                 
%                 c=obj.real_sweep;
%                 Lx=[1       0       0
%                     0   cosd(a)  sind(a)
%                     0   -sind(a) cosd(a)];
%                 
%                 Ly=[cosd(b) 0 sind(b)
%                     0      1    0
%                     -sind(b)  0   cosd(b)];
%                 
%                 Lz=[cosd(c) sind(c)   0
%                     -sind(c) cosd(c)  0
%                     0           0   1];
%                 
%                 T=Lz*Ly*Lx;
%                 
%                 p1=obj.pos+(T*p1')';
%                 p2=obj.pos+(T*p2')';
%                 p3=obj.pos+(T*p3')';
%                 p4=obj.pos+(T*p4')';
%                 obj.xyz=[p1' p2' p3' p4'];
                    disp('warning, this does not work properly')
            end
            
        end
        
        function obj=compute_controlsurface_coordinates(obj)                
            if obj.has_te_cs
                obj.delta_te_device=obj.te_device.delta;
            end
            if obj.has_le_cs
                obj.delta_le_device=obj.le_device.delta;
            end
            if obj.has_sp_cs
                obj.delta_sp_device=obj.sp_device.delta;
            end 
            
            p1=obj.xyz(:,1)';
            p2=obj.xyz(:,2)';
            p3=obj.xyz(:,3)';
            p4=obj.xyz(:,4)';
            
            if obj.c_te_device~=0
                if obj.te_device.is_tapered==1
                      t1=norm(p4-p1);
                      t2=norm(p3-p2);
                      p1=p1*(sum(obj.c_te_device)/t1)+p4*(1-sum(obj.c_te_device)/t1);
                      p2=p2*(sum(obj.c_te_device)*obj.c_t/obj.c_r/t2)+p3*(1-sum(obj.c_te_device)*obj.c_t/obj.c_r/t2);
%                     p1=p4+[-sum(obj.c_te_device)*cos(obj.Theta_r*pi/180) -sum(obj.c_te_device)*sin(obj.Theta_r*pi/180)*sin(obj.dihed*pi/180) sum(obj.c_te_device)*sin(obj.Theta_r*pi/180)*cos(obj.dihed*pi/180)];
%                     p2=p3+[-sum(obj.c_te_device)*cos(obj.Theta_t*pi/180) -sum(obj.c_te_device)*sin(obj.Theta_t*pi/180)*sin(obj.dihed*pi/180) sum(obj.c_te_device)*sin(obj.Theta_t*pi/180)*cos(obj.dihed*pi/180)]*obj.c_t/obj.c_r;
               

                % in the else loop below t2 was introduced (absent in
                % previous versions of the code). This is fundamental to
                % calculate coordinates of an untapered TE CS (it was
                % always calculated tapered before, regardless of the
                % is_tapered attribute value)
                else
                      t1=norm(p4-p1);
                      t2=norm(p3-p2);
                      p1 = p1*(sum(obj.c_te_device)/t1)+p4*(1-sum(obj.c_te_device)/t1);
                      p2 = p2*(sum(obj.c_te_device)/t2)+p3*(1-sum(obj.c_te_device)/t2);
 
%                     p1=p4+[-sum(obj.c_te_device)*cos(obj.Theta_r*pi/180) -sum(obj.c_te_device)*sin(obj.Theta_r*pi/180)*sin(obj.dihed*pi/180) sum(obj.c_te_device)*sin(obj.Theta_r*pi/180)*cos(obj.dihed*pi/180)];
%                     p2=p3+[-sum(obj.c_te_device)*cos(obj.Theta_t*pi/180) -sum(obj.c_te_device)*sin(obj.Theta_t*pi/180)*sin(obj.dihed*pi/180) sum(obj.c_te_device)*sin(obj.Theta_t*pi/180)*cos(obj.dihed*pi/180)];
                end
                for i=1:length(obj.c_te_device)
                    %if obj.te_device.is_tapered==1
                    
                    hinge_vec=p2-p1;
                    
                    u=hinge_vec(1);
                    v=hinge_vec(2);
                    w=hinge_vec(3);
                    
                    x=p3(1);
                    y=p3(2);
                    z=p3(3);
                    
                    a=p1(1);
                    b=p1(2);
                    c=p1(3);
                    
                    Theta=obj.delta_te_device(1)*pi/180;
                    L=u^2+v^2+w^2;
                    
                    p3=[(a*(v^2+w^2)-u*(b*v+c*w-u*x-v*y-w*z))*(1-cos(Theta))+L*x*cos(Theta)+sqrt(L)*(-c*v+b*w-w*y+v*z)*sin(Theta);
                          (b*(u^2+w^2)-v*(a*u+c*w-u*x-v*y-w*z))*(1-cos(Theta))+L*y*cos(Theta)+sqrt(L)*(c*u-a*w+w*x-u*z)*sin(Theta);
                          (c*(u^2+v^2)-w*(a*u+b*v-u*x-v*y-w*z))*(1-cos(Theta))+L*z*cos(Theta)+sqrt(L)*(-b*u+a*v-v*x+u*y)*sin(Theta);]'/L;
                      
                    x=p4(1);
                    y=p4(2);
                    z=p4(3);  
                      
                      
                     p4=[(a*(v^2+w^2)-u*(b*v+c*w-u*x-v*y-w*z))*(1-cos(Theta))+L*x*cos(Theta)+sqrt(L)*(-c*v+b*w-w*y+v*z)*sin(Theta);
                          (b*(u^2+w^2)-v*(a*u+c*w-u*x-v*y-w*z))*(1-cos(Theta))+L*y*cos(Theta)+sqrt(L)*(c*u-a*w+w*x-u*z)*sin(Theta);
                          (c*(u^2+v^2)-w*(a*u+b*v-u*x-v*y-w*z))*(1-cos(Theta))+L*z*cos(Theta)+sqrt(L)*(-b*u+a*v-v*x+u*y)*sin(Theta);]'/L;
                      
                    obj.xyz_te_device(i,:,:)=[p1' p2' p3' p4'];
                    p1=p4;
                    p2=p3;
                end
            end
            
            if obj.c_le_device~=0
                
                obj.xyz_le_device=zeros(length(obj.c_le_device),3,4);
                
                % Before the taper_correction was introduced below, the LE
                % CS was always considered untapered, regardless of its
                % is_tapered attribute value.
                if obj.le_device.is_tapered == 0
                    taper_correction = 1;
                else
                    taper_correction = obj.c_r/obj.c_t;
                end
                    
                % Location where taper_correction was introduced (p3)
                p4=obj.xyz(:,1)'-[-sum(obj.c_le_device)*cos(obj.Theta_r*pi/180) -(sum(obj.c_le_device))*sin(obj.Theta_r*pi/180)*sin(obj.dihed*pi/180) +(sum(obj.c_le_device))*sin(obj.Theta_r*pi/180)*cos(obj.dihed*pi/180)];
                p3=obj.xyz(:,2)'-[-sum(obj.c_le_device)*cos(obj.Theta_t*pi/180) -(sum(obj.c_le_device))*sin(obj.Theta_t*pi/180)*sin(obj.dihed*pi/180) +(sum(obj.c_le_device))*sin(obj.Theta_t*pi/180)*cos(obj.dihed*pi/180)] / taper_correction;
                
                for i=length(obj.c_le_device):-1:1
                    dvec1=[(cos(obj.delta_le_device(i)*pi/180)*obj.c_le_device(i))*cos(obj.Theta_r*pi/180) (cos(obj.delta_le_device(i)*pi/180)*obj.c_le_device(i))*sin(obj.Theta_r*pi/180)*sin(obj.dihed*pi/180) -(cos(obj.delta_le_device(i)*pi/180)*obj.c_le_device(i))*sin(obj.Theta_r*pi/180)*cos(obj.dihed*pi/180)];
                    dvec3=[(cos(obj.delta_le_device(i)*pi/180)*obj.c_le_device(i))*cos(obj.Theta_t*pi/180) (cos(obj.delta_le_device(i)*pi/180)*obj.c_le_device(i))*sin(obj.Theta_t*pi/180)*sin(obj.dihed*pi/180) -(cos(obj.delta_le_device(i)*pi/180)*obj.c_le_device(i))*sin(obj.Theta_t*pi/180)*cos(obj.dihed*pi/180)];
                    dvec2=[obj.c_le_device(i)*sin(obj.delta_le_device(i)*pi/180)*cos(obj.dihed*pi/180)*sin(obj.Theta_r*pi/180) -obj.c_le_device(i)*sin(obj.delta_le_device(i)*pi/180)*sin(obj.dihed*pi/180) obj.c_le_device(i)*sin(obj.delta_le_device(i)*pi/180)*cos(obj.dihed*pi/180)*cos(obj.Theta_r*pi/180)];
                    dvec4=[obj.c_le_device(i)*sin(obj.delta_le_device(i)*pi/180)*cos(obj.dihed*pi/180)*sin(obj.Theta_t*pi/180) -obj.c_le_device(i)*sin(obj.delta_le_device(i)*pi/180)*sin(obj.dihed*pi/180) obj.c_le_device(i)*sin(obj.delta_le_device(i)*pi/180)*cos(obj.dihed*pi/180)*cos(obj.Theta_t*pi/180)];
                    
                    % Location where taper_correction was introduced (p2)
                    p2=p3-(dvec3+dvec4) / taper_correction;
                    p1=p4- (dvec1+dvec2);
                    
                    obj.xyz_le_device(i,:,:)=[p1' p2' p3' p4'];
                    
                    p4=p1;
                    p3=p2;
                end
            end
            
            if obj.c_sp_device~=0
                
                obj.xyz_sp_device=zeros(length(obj.c_sp_device),3,4);
                
                % equivalents of vectors t1 and t2 used in TE surface
                vect_root = obj.xyz(:,4) - obj.xyz(:,1);
                vect_tip = obj.xyz(:,3) - obj.xyz(:,2);
                
                % if the CS are constant (untapered), point nr 3 is sized
                % by relating the chord of the Spoiler to the segment ROOT
                % CHORD, rather than the tip chord, which is used if
                % tapered.
                if obj.sp_device.is_tapered==0
                    
                    % if innerBorder_coords is empty, we are reading from
                    % our in-house XML, so our variables for defining the
                    % Spoiler coordinats are used. If it is not empty, then
                    % CPACS is being used, and the CPACS variables are
                    % called.
                    if isempty(obj.sp_device.innerBorder_coords)

                        p1 = obj.xyz(:,1) + vect_root * obj.sp_device.xsi_LE_inner;
                        p4 = obj.xyz(:,1) + vect_root *(obj.sp_device.xsi_LE_inner + obj.sp_device.hinge);
                        p2 = obj.xyz(:,2) + vect_tip * obj.sp_device.xsi_LE_outer;
                        p3 = obj.xyz(:,2) + vect_tip * obj.sp_device.xsi_LE_outer + obj.sp_device.hinge * vect_root;

                    else
                        p1 = obj.xyz(:,1) + vect_root * obj.sp_device.innerBorder_coords(3);
                        p4 = obj.xyz(:,1) + vect_root * obj.sp_device.innerBorder_coords(4);
                        p2 = obj.xyz(:,2) + vect_tip * obj.sp_device.outerBorder_coords(3);
                        p3 = obj.xyz(:,2) + vect_tip * obj.sp_device.outerBorder_coords(3) + vect_root * (obj.sp_device.innerBorder_coords(4) - obj.sp_device.innerBorder_coords(3));
                    end
                
                % same comments as the "parent if statement". Only
                % difference is that this is for tapered surfaces, so point
                % 3 is sized using the CS hinge and tip chord.
                else
                    if isempty(obj.sp_device.innerBorder_coords)

                        p1 = obj.xyz(:,1) + vect_root * obj.sp_device.xsi_LE_inner;
                        p4 = obj.xyz(:,1) + vect_root *(obj.sp_device.xsi_LE_inner + obj.sp_device.hinge);
                        p2 = obj.xyz(:,2) + vect_tip * obj.sp_device.xsi_LE_outer;
                        p3 = obj.xyz(:,2) + vect_tip * (obj.sp_device.xsi_LE_outer + obj.sp_device.hinge);

                    else
                        p1 = obj.xyz(:,1) + vect_root * obj.sp_device.innerBorder_coords(3);
                        p4 = obj.xyz(:,1) + vect_root * obj.sp_device.innerBorder_coords(4);
                        p2 = obj.xyz(:,2) + vect_tip * obj.sp_device.outerBorder_coords(3);
                        p3 = obj.xyz(:,2) + vect_tip * obj.sp_device.outerBorder_coords(4);
                    end
                end
                
                obj.xyz_sp_device(1,:,:)=[p1, p2, p3, p4];
            end
        end
        
        function obj=compute_xyz_fixed(obj)
            if isempty(obj.c_te_device)&& isempty(obj.c_le_device)
                obj.xyz_fixed=obj.xyz;
            elseif isempty(obj.c_le_device)
                obj.xyz_fixed=[obj.xyz(:,1) obj.xyz(:,2) squeeze(obj.xyz_te_device(1,:,2))' squeeze(obj.xyz_te_device(1,:,1))'];
            elseif isempty(obj.c_te_device)
                obj.xyz_fixed=[squeeze(obj.xyz_le_device(end,:,4))' squeeze(obj.xyz_le_device(end,:,3))' obj.xyz(:,3) obj.xyz(:,4)];
            elseif ~isempty(obj.c_te_device)&&~isempty(obj.c_le_device)
                obj.xyz_fixed=[squeeze(obj.xyz_le_device(end,:,4))' squeeze(obj.xyz_le_device(end,:,3))'  squeeze(obj.xyz_te_device(1,:,2))' squeeze(obj.xyz_te_device(1,:,1))'];
            end
        end
        
        function obj = read_wingbox_coords(obj,pathNodeCoords,n,iNodes)
            span_grid_aero=0:1/obj.n_span:1;
            front_coords=[];
            rear_coords=[];
            c4_coords=[];
            
            nodeCoords = importdata(pathNodeCoords);
            segmentNodeCoords = nodeCoords(iNodes:iNodes+n,:);
            %to find the c4 coords, we need to project each nodeCoord on
            %the segments one quarter line
            %vector of c4 line
            c4Vec = ((obj.xyz(:,2)*0.75+obj.xyz(:,3)*0.25)-(obj.xyz(:,1)*0.75+obj.xyz(:,4)*0.25));
            c4Vec2=dot(c4Vec,c4Vec);
            for i=1:n+1
                %vector from inboard c4 point to nodeCoord
                c4inboardnC = (segmentNodeCoords(i,:)'-(obj.xyz(:,1)*0.75+obj.xyz(:,4)*0.25));
                %portion of the c4 coord where the shortest distance happens 
                p = dot(c4inboardnC,c4Vec)/c4Vec2;
                %check if p is less than 0
                if p<0
                    p=0;
                elseif p>1
                    p=1;
                end
                c4_coords(:,i) =  ( obj.xyz(:,1)*0.75 + obj.xyz(:,4) * 0.25 ) + p * c4Vec;    
            end
            
            obj.wingbox_coords=zeros(3,n+1,2);
            obj.wingbox_coords(:,:,1) = segmentNodeCoords';
            obj.wingbox_coords(:,:,2) = segmentNodeCoords';
            obj.wingbox_c4 = c4_coords;
%             
            for i=1:obj.n_span
                obj.c4_coords(:,i)=obj.xyz(:,4)+(obj.xyz(:,1)-obj.xyz(:,4))*0.75+span_grid_aero(i)*(obj.xyz(:,3)+(obj.xyz(:,2)-obj.xyz(:,3))*0.75-obj.xyz(:,4)-(obj.xyz(:,1)-obj.xyz(:,4))*0.75);
            end
        end
        
        function obj=compute_wingbox_coords(obj,frontspar,rearspar,n)
            %% TODO: generalize for n number of spars
            span_grid=0:1/n:1;
            span_grid_aero=0:1/obj.n_span:1;
            front_coords=[];
            rear_coords=[];
            c4_coords=[];
            
            if isempty(frontspar) && isempty(rearspar)
                eta_input = [0, 1];
                xsi_fs_input = [obj.structural_properties.fs_segments(:, 1)', obj.structural_properties.fs_segments(end, 2)];
                xsi_rs_input = [obj.structural_properties.rs_segments(:, 1)', obj.structural_properties.rs_segments(end, 2)];

                front_sp = interp1(eta_input, xsi_fs_input, span_grid, 'lin');
                rear_sp = interp1(eta_input, xsi_rs_input, span_grid, 'lin');
                
                eta_input = eta_input(1:end-1);
                t_fs_input = obj.structural_properties.fs_segments(:, 3)';
                t_rs_input = obj.structural_properties.rs_segments(:, 3)';
                t_ts_input = (obj.structural_properties.fs_segments(:, 4)' + obj.structural_properties.rs_segments(:, 4)')/2;
                t_bs_input = (obj.structural_properties.fs_segments(:, 5)' + obj.structural_properties.rs_segments(:, 5)')/2;
                
                if length(eta_input) == 1
                    eta_input = [0, 1];
                    t_fs_input = [t_fs_input, t_fs_input];
                    t_rs_input = [t_rs_input, t_rs_input];
                    t_ts_input = [t_ts_input, t_ts_input];
                    t_bs_input = [t_bs_input, t_bs_input];
                end
                
                t_fs = interp1(eta_input, t_fs_input, span_grid, 'lin', 'extrap');
                t_rs = interp1(eta_input, t_rs_input, span_grid, 'lin', 'extrap');
                t_ts = interp1(eta_input, t_ts_input, span_grid, 'lin', 'extrap');
                t_bs = interp1(eta_input, t_bs_input, span_grid, 'lin', 'extrap');
                
                obj.structural_properties.t_fs = t_fs;
                obj.structural_properties.t_rs = t_rs;
                obj.structural_properties.t_ts = t_ts;
                obj.structural_properties.t_bs = t_bs;
            else
                front_sp=frontspar(1)*(1-span_grid)+frontspar(2)*span_grid;
                rear_sp=rearspar(1)*(1-span_grid)+rearspar(2)*span_grid;
            end
            
            for i=1:n+1
                front_coords(:,i)=(obj.xyz(:,1)*(1-front_sp(i))+obj.xyz(:,4)*front_sp(i))+span_grid(i)*(obj.xyz(:,2)*(1-front_sp(i))+obj.xyz(:,3)*front_sp(i)-obj.xyz(:,1)*(1-front_sp(i))-obj.xyz(:,4)*front_sp(i));
                rear_coords(:,i)=(obj.xyz(:,1)*(1-rear_sp(i))+obj.xyz(:,4)*rear_sp(i))+span_grid(i)*(obj.xyz(:,2)*(1-rear_sp(i))+obj.xyz(:,3)*rear_sp(i)-obj.xyz(:,1)*(1-rear_sp(i))-obj.xyz(:,4)*rear_sp(i));
                c4_coords(:,i)=(obj.xyz(:,1)*0.75+obj.xyz(:,4)*0.25)+span_grid(i)*(obj.xyz(:,2)*0.75+obj.xyz(:,3)*0.25-obj.xyz(:,1)*0.75-obj.xyz(:,4)*0.25);
            end
            
            obj.wingbox_coords=zeros(3,n+1,2);
            obj.wingbox_coords(:,:,1)=front_coords(:,:);
            obj.wingbox_coords(:,:,2)=rear_coords(:,:);
            obj.wingbox_c4=c4_coords;
            
            nprof_upper_t=obj.profile_t(1,1);
            nprof_lower_t=obj.profile_t(1,2);
            coords_upper_t=obj.profile_t(2:1+nprof_upper_t,1);
            profile_upper_t=obj.profile_t(2:1+nprof_upper_t,2);
            coords_lower_t=obj.profile_t(2+nprof_upper_t:1+nprof_upper_t+nprof_lower_t,1);
            profile_lower_t=obj.profile_t(2+nprof_upper_t:1+nprof_upper_t+nprof_lower_t,2);
            
            nprof_upper_r=obj.profile_r(1,1);
            nprof_lower_r=obj.profile_r(1,2);
            coords_upper_r=obj.profile_r(2:1+nprof_upper_r,1);
            profile_upper_r=obj.profile_r(2:1+nprof_upper_r,2);
            coords_lower_r=obj.profile_r(2+nprof_upper_r:1+nprof_upper_r+nprof_lower_r,1);
            profile_lower_r=obj.profile_r(2+nprof_upper_r:1+nprof_upper_r+nprof_lower_r,2);
            
            h_fr=zeros(1,length(front_sp));
            h_re=zeros(1,length(front_sp));
            
            for j=1:length(front_sp)
                profile_height_t=abs(interp1(coords_upper_t,profile_upper_t,front_sp(j),'lin','extrap')-interp1(coords_lower_t,profile_lower_t,front_sp(j),'lin','extrap'))*obj.c_t;
                profile_height_r=abs(interp1(coords_lower_r,profile_lower_r,front_sp(j),'lin','extrap')-interp1(coords_upper_r,profile_upper_r,front_sp(j),'lin','extrap'))*obj.c_r;
                h_fr(j)=profile_height_r*(1-span_grid(j))+profile_height_t*span_grid(j);
                
                profile_height_t=abs(interp1(coords_upper_t,profile_upper_t,rear_sp(j),'lin','extrap')-interp1(coords_lower_t,profile_lower_t,rear_sp(j),'lin','extrap'))*obj.c_t;
                profile_height_r=abs(interp1(coords_lower_r,profile_lower_r,rear_sp(j),'lin','extrap')-interp1(coords_upper_r,profile_upper_r,rear_sp(j),'lin','extrap'))*obj.c_r;
                h_re(j)=profile_height_r*(1-span_grid(j))+profile_height_t*span_grid(j);
            end
            obj.wingbox_height(:,1)=h_fr;
            obj.wingbox_height(:,2)=h_re;
            
            for i=1:obj.n_span
                obj.c4_coords(:,i)=obj.xyz(:,4)+(obj.xyz(:,1)-obj.xyz(:,4))*0.75+span_grid_aero(i)*(obj.xyz(:,3)+(obj.xyz(:,2)-obj.xyz(:,3))*0.75-obj.xyz(:,4)-(obj.xyz(:,1)-obj.xyz(:,4))*0.75);
            end
            
            
        end
        
        function obj=plot_segment(obj)
            
            if ~isempty(obj.te_device)
                for i=1:length(obj.c_te_device)
                    if ~obj.te_device.is_sym_defl
                        obj.te_device.delta=obj.te_device.delta_l_r(1);
                    end
                    obj=obj.compute_controlsurface_coordinates();
                    handle=fill3(squeeze(obj.xyz_te_device(i,1,:))',squeeze(obj.xyz_te_device(i,2,:))',squeeze(obj.xyz_te_device(i,3,:))','b');
                    alpha(handle,0.4)
                    if obj.symmetric==1
                        if obj.te_device.is_sym_defl
                            handle=fill3(squeeze(obj.xyz_te_device(i,1,:))',-squeeze(obj.xyz_te_device(i,2,:))',squeeze(obj.xyz_te_device(i,3,:))','b');
                            alpha(handle,0.4)
                        else
%                             obj.te_device.delta
%                             obj.te_device.name
%                             zw=obj.te_device.delta(1);
                            obj.te_device.delta=obj.te_device.delta_l_r(2)*-1;
                            obj=obj.compute_controlsurface_coordinates();
                            handle=fill3(squeeze(obj.xyz_te_device(i,1,:))',-squeeze(obj.xyz_te_device(i,2,:))',squeeze(obj.xyz_te_device(i,3,:))','b');
                            alpha(handle,0.4)
                            obj.te_device.delta=obj.te_device.delta_l_r(1);
                            obj=obj.compute_controlsurface_coordinates();
                        end
                    end
                end
            end
            
            if ~isempty(obj.le_device)
                for i=1:length(obj.c_le_device)
                    handle=fill3(squeeze(obj.xyz_le_device(i,1,:))',squeeze(obj.xyz_le_device(i,2,:))',squeeze(obj.xyz_le_device(i,3,:))','g');
                    alpha(handle,0.4)
                    if obj.symmetric==1
                        if obj.le_device.is_sym_defl
                            handle=fill3(squeeze(obj.xyz_le_device(i,1,:))',-squeeze(obj.xyz_le_device(i,2,:))',squeeze(obj.xyz_le_device(i,3,:))','g');
                            alpha(handle,0.4);
                        else
                            obj.le_device.delta=obj.le_device.delta*-1;
                            obj=obj.compute_controlsurface_coordinates();
                            handle=fill3(squeeze(obj.xyz_le_device(i,1,:))',-squeeze(obj.xyz_le_device(i,2,:))',squeeze(obj.xyz_le_device(i,3,:))','g');
                            alpha(handle,0.4);
                            obj.le_device.delta=obj.le_device.delta*-1;
                            obj=obj.compute_controlsurface_coordinates();
                        end
                    end
                end
            end
            
            if ~isempty(obj.sp_device)
                for i=1:length(obj.c_sp_device)
                    handle=fill3(squeeze(obj.xyz_sp_device(i,1,:))',squeeze(obj.xyz_sp_device(i,2,:))',squeeze(obj.xyz_sp_device(i,3,:)+0.05)','y');
                    alpha(handle,0.4)
                    if obj.symmetric==1
                        if obj.sp_device.is_sym_defl
                            handle=fill3(squeeze(obj.xyz_sp_device(i,1,:))',-squeeze(obj.xyz_sp_device(i,2,:))',squeeze(obj.xyz_sp_device(i,3,:)+0.05)','y');
                            alpha(handle,0.4);
                        else
                            obj.sp_device.delta=obj.sp_device.delta*-1;
                            obj=obj.compute_controlsurface_coordinates();
                            handle=fill3(squeeze(obj.xyz_sp_device(i,1,:))',-squeeze(obj.xyz_sp_device(i,2,:))',squeeze(obj.xyz_sp_device(i,3,:)+0.05)','y');
                            alpha(handle,0.4);
                            obj.sp_device.delta=obj.sp_device.delta*-1;
                            obj=obj.compute_controlsurface_coordinates();
                        end
                    end
                end
            end
            
            handle=fill3(obj.xyz_fixed(1,:),obj.xyz_fixed(2,:),obj.xyz_fixed(3,:),'r');
            hold on
            alpha(handle,0.4)
            if obj.symmetric==1
                handle=fill3(obj.xyz_fixed(1,:),-obj.xyz_fixed(2,:),obj.xyz_fixed(3,:),'r');
                hold on
                alpha(handle,0.4)
            end
            axis equal
        end
        
        function pos=get_tip_c4point(obj)
            pos=obj.pos+[obj.b*tan(obj.c4_sweep*pi/180) obj.b*cos(obj.dihed*pi/180) +obj.b*sin(obj.dihed*pi/180)];
        end
        
        function c=get_tip_chord(obj)
            c=obj.c_t;
        end
        
        function obj=compute_forces(obj,panel_forces,panels,grid)
            n_pan=sum([obj.n_le_panels obj.n_ctr1_panels obj.n_sp_panels obj.n_ctr2_panels  obj.n_te_panels]);
            obj.span_forces=zeros(3,obj.n_span);
            obj.span_moments=zeros(3,obj.n_span);
            span_grid=0.5/obj.n_span:1/obj.n_span:1-0.5/obj.n_span;
            for i=1:obj.n_span

                %zw=obj.c4_coords(:,i);
%                 hold on
%                                 plot3(zw(1,:),zw(2,:),zw(3,:),'o-');
%                 %                pause
                for j=1:n_pan
                    idx=obj.panel_start_idx+n_pan*(i-1)-1+j;

                    obj.span_forces(:,i)=obj.span_forces(:,i)+panel_forces(:,idx)/(obj.b/obj.n_span);
                    r=grid(:,panels(4,idx))+(grid(:,panels(1,idx))-grid(:,panels(4,idx)))*0.75+0.5*(grid(:,panels(3,idx))+(grid(:,panels(2,idx))-grid(:,panels(3,idx)))*0.75-grid(:,panels(4,idx))-(grid(:,panels(1,idx))-grid(:,panels(4,idx)))*0.75);
%                     hold on
%                      fill3(grid(1,panels(:,idx)),grid(2,panels(:,idx)),grid(3,panels(:,idx)),'r');
                    obj.span_moments(:,i)=obj.span_moments(:,i)+cross(r-obj.c4_coords(:,i),panel_forces(:,idx))/(obj.b/obj.n_span);
                    % plot3(zw(1,:),zw(2,:),zw(3,:),'x-');
                end
            end
        end
        
        function grid=compute_deflected_grid_left(obj,panels,grid,deflections_structmesh,offset)
            span_grid=0:1/obj.n_span:1;
            n_pan=sum([obj.n_le_panels obj.n_ctr1_panels obj.n_sp_panels obj.n_ctr2_panels  obj.n_te_panels]);
            
            le_ctr=1;
            le_loc_ctr=1;
            te_ctr=1;
            te_loc_ctr=1;
% %             
                deflections_structmesh(4,:)=deflections_structmesh(4,:);
% %              % deflections_structmesh(5,:)=-deflections_structmesh(5,:);
          deflections_structmesh(6,:)=deflections_structmesh(6,:);
    
%                            deflections_structmesh(4,:)=deflections_structmesh(4,:);
%              % deflections_structmesh(5,:)=-deflections_structmesh(5,:);
%                deflections_structmesh(6,:)=deflections_structmesh(6,:);   
%                
            for i=1:obj.n_span
                Theta=(obj.Theta_r*(1-span_grid(i))+obj.Theta_t*span_grid(i))*pi/180*0;
                
                a=-obj.dihed*pi/180;
                %a=obj.dihed*pi/180;
                c=obj.c4_sweep;
                
                Lx=[1       0       0
                    0   cos(a)  sin(a)
                    0   -sin(a) cos(a)];
                
                
                Lz=[cos(c) sin(c)   0
                    -sin(c) cos(c)  0
                    0           0   1];
                
                R=Lz*Lx*[0;Theta;0];

                c4_coords_edge(:,i)=(obj.xyz(:,1)*0.75+obj.xyz(:,4)*0.25)+span_grid(i)*(obj.xyz(:,2)*0.75+obj.xyz(:,3)*0.25-obj.xyz(:,1)*0.75-obj.xyz(:,4)*0.25);
                c4_coords_edge(2,i)=-c4_coords_edge(2,i);
%                  hold on
%                  plot3(c4_coords_edge(1,i),c4_coords_edge(2,i),c4_coords_edge(3,i),'x')
                c4_coords_edge(:,i)=c4_coords_edge(:,i)+deflections_structmesh(1:3,i);
                
            %     plot3(c4_coords_edge(1,i),c4_coords_edge(2,i),c4_coords_edge(3,i),'ro')
                for j=1:n_pan
                    idx=obj.panel_start_idx+n_pan*(i-1)+j-1+offset;
                    
                    grid(:,panels(2,idx))=grid(:,panels(2,idx))+deflections_structmesh(1:3,i);
                    
                    dist_x=norm(grid(1:3,panels(2,idx))-c4_coords_edge(1:3,i));
                    sgnx=sign(grid(1,panels(2,idx))-c4_coords_edge(1,i));
                    
                    dTheta=deflections_structmesh(5,i);
                    Theta=R(2);
                    dx2=-sgnx*(dist_x*cos(Theta)-(dist_x*cos(Theta+dTheta)));
                    dz2=-sgnx*(dist_x*sin(Theta+dTheta)-dist_x*sin(Theta));
                    Theta=R(3);
                    dTheta=deflections_structmesh(6,i);
                    
                    dx3=-sgnx*(dist_x*cos(Theta)-(dist_x*cos(Theta+dTheta)));
                    dy3=-sgnx*(dist_x*sin(Theta+dTheta)-dist_x*sin(Theta));
                    delta_twist_1=[dx2+dx3;-dy3;dz2];
                    
                    grid(:,panels(2,idx))=grid(:,panels(2,idx))+delta_twist_1;
                    
                   % plot3(grid(1,panels(2,idx)),grid(2,panels(2,idx)),grid(3,panels(2,idx)),'go')
                    if obj.has_le_cs
                        if le_ctr<=length(obj.n_le_panels)
                            le_loc_ctr=le_loc_ctr+1;
                            if le_loc_ctr>obj.n_le_panels(le_ctr)
                                le_loc_ctr=1;
                                le_ctr=le_ctr+1;
                                grid(:,panels(3,idx))=grid(:,panels(3,idx))+deflections_structmesh(1:3,i);
                                dist_x=abs(grid(1,panels(3,idx))-c4_coords_edge(1,i));
                                sgnx=sign(grid(1,panels(3,idx))-c4_coords_edge(1,i));
                               % plot3(c4_coords_edge(1,i)-grid(1,panels(4,idx)),c4_coords_edge(2,i),c4_coords_edge(3,i),'cx')
                                dTheta=deflections_structmesh(5,i);
                                Theta=R(2);
                                dx2=-sgnx*(dist_x-(dist_x/cos(Theta)*cos(Theta+dTheta)));
                                dz2=-sgnx*(dist_x/cos(Theta)*sin(Theta+dTheta)-dist_x/cos(Theta)*sin(Theta));
                                dTheta=deflections_structmesh(6,i);
                                Theta=R(3);
                                dx3=-sgnx*(dist_x-(dist_x/cos(Theta)*cos(Theta+dTheta)));
                                dy3=-sgnx*(dist_x/cos(Theta)*sin(Theta+dTheta)-dist_x/cos(Theta)*sin(Theta));
                                delta_twist_4=[dx2+dx3;-dy3;dz2];
                                grid(:,panels(3,idx))=grid(:,panels(3,idx))+delta_twist_4;
                            end
                            
                        end
                    end
                    
                    if obj.has_te_cs
                        if j>n_pan-sum(obj.n_te_panels)
                            if te_ctr<=length(obj.n_te_panels)
                                te_loc_ctr=te_loc_ctr+1;
                                if te_loc_ctr>obj.n_te_panels(te_ctr)
                                    te_loc_ctr=1;
                                    te_ctr=te_ctr+1;
                                    grid(:,panels(3,idx))=grid(:,panels(3,idx))+deflections_structmesh(1:3,i);
                                    dist_x=abs(grid(1,panels(3,idx))-c4_coords_edge(1,i));
                                    sgnx=sign(grid(1,panels(3,idx))-c4_coords_edge(1,i));
                                    dTheta=deflections_structmesh(5,i);
                                    Theta=R(2);
                                   %  plot3(c4_coords_edge(1,i)-grid(1,panels(4,idx)),c4_coords_edge(2,i),c4_coords_edge(3,i),'cx')
                                    dx2=-sgnx*(dist_x-(dist_x/cos(Theta)*cos(Theta+dTheta)));
                                    dz2=-sgnx*(dist_x/cos(Theta)*sin(Theta+dTheta)-dist_x/cos(Theta)*sin(Theta));
                                    dTheta=deflections_structmesh(6,i);
                                    Theta=R(3);
                                    dx3=-sgnx*(dist_x-(dist_x/cos(Theta)*cos(Theta+dTheta)));
                                    dy3=-sgnx*(dist_x/cos(Theta)*sin(Theta+dTheta)-dist_x/cos(Theta)*sin(Theta));
                                    delta_twist_4=[dx2+dx3;-dy3;dz2];
                                    grid(:,panels(3,idx))=grid(:,panels(3,idx))+delta_twist_4;
                                end
                            end
                        elseif j==n_pan-sum(obj.n_te_panels)
                            grid(:,panels(3,idx))=grid(:,panels(3,idx))+deflections_structmesh(1:3,i);
                            dist_x=abs(grid(1,panels(3,idx))-c4_coords_edge(1,i));
                            sgnx=sign(grid(1,panels(3,idx))-c4_coords_edge(1,i));
                            dTheta=deflections_structmesh(5,i);
                            Theta=R(2);
                           %  plot3(c4_coords_edge(1,i)-grid(1,panels(4,idx)),c4_coords_edge(2,i),c4_coords_edge(3,i),'cx')
                            dx2=-sgnx*(dist_x-(dist_x/cos(Theta)*cos(Theta+dTheta)));
                            dz2=-sgnx*(dist_x/cos(Theta)*sin(Theta+dTheta)-dist_x/cos(Theta)*sin(Theta));
                            dTheta=deflections_structmesh(6,i);
                            Theta=R(3);
                            dx3=-sgnx*(dist_x-(dist_x/cos(Theta)*cos(Theta+dTheta)));
                            dy3=-sgnx*(dist_x/cos(Theta)*sin(Theta+dTheta)-dist_x/cos(Theta)*sin(Theta));
                            delta_twist_4=[dx2+dx3;-dy3;dz2];
                            grid(:,panels(3,idx))=grid(:,panels(3,idx))+delta_twist_4;
                        end
                        
                    end
                    if obj.has_sp_cs
                        if j==obj.n_ctr1_panels+sum(obj.n_le_panels)
                            grid(:,panels(3,idx))=grid(:,panels(3,idx))+deflections_structmesh(1:3,i);
                                dist_x=abs(grid(1,panels(3,idx))-c4_coords_edge(1,i));
                            sgnx=sign(grid(1,panels(3,idx))-c4_coords_edge(1,i));
                            dTheta=deflections_structmesh(5,i);
                            Theta=R(2);
                           %  plot3(c4_coords_edge(1,i)-grid(1,panels(4,idx)),c4_coords_edge(2,i),c4_coords_edge(3,i),'cx')
                            dx2=-sgnx*(dist_x-(dist_x/cos(Theta)*cos(Theta+dTheta)));
                            dz2=-sgnx*(dist_x/cos(Theta)*sin(Theta+dTheta)-dist_x/cos(Theta)*sin(Theta));
                            dTheta=deflections_structmesh(6,i);
                            Theta=R(3);
                            dx3=-sgnx*(dist_x-(dist_x/cos(Theta)*cos(Theta+dTheta)));
                            dy3=-sgnx*(dist_x/cos(Theta)*sin(Theta+dTheta)-dist_x/cos(Theta)*sin(Theta));
                            delta_twist_4=[dx2+dx3;-dy3;dz2];
                            grid(:,panels(3,idx))=grid(:,panels(3,idx))+delta_twist_4;
                        end
                    end
                    if and(obj.has_sp_cs,~obj.has_te_cs)
                        if j==obj.n_ctr1_panels+sum(obj.n_le_panels)+obj.n_sp_panels
                            grid(:,panels(3,idx))=grid(:,panels(3,idx))+deflections_structmesh(1:3,i);
                            dist_x=abs(grid(1,panels(3,idx))-c4_coords_edge(1,i));
                            sgnx=sign(grid(1,panels(3,idx))-c4_coords_edge(1,i));
                            dTheta=deflections_structmesh(5,i);
                            Theta=R(2);
                           %  plot3(c4_coords_edge(1,i)-grid(1,panels(4,idx)),c4_coords_edge(2,i),c4_coords_edge(3,i),'cx')
                            dx2=-sgnx*(dist_x-(dist_x/cos(Theta)*cos(Theta+dTheta)));
                            dz2=-sgnx*(dist_x/cos(Theta)*sin(Theta+dTheta)-dist_x/cos(Theta)*sin(Theta));
                            dTheta=deflections_structmesh(6,i);
                            Theta=R(3);
                            dx3=-sgnx*(dist_x-(dist_x/cos(Theta)*cos(Theta+dTheta)));
                            dy3=-sgnx*(dist_x/cos(Theta)*sin(Theta+dTheta)-dist_x/cos(Theta)*sin(Theta));
                            delta_twist_4=[dx2+dx3;-dy3;dz2];
                            grid(:,panels(3,idx))=grid(:,panels(3,idx))+delta_twist_4;
                        end
                    end
                end
                if  (isempty(obj.has_te_cs))||(obj.has_te_cs==0)
                    grid(:,panels(3,idx))=grid(:,panels(3,idx))+deflections_structmesh(1:3,i);
                    dist_x=abs(grid(1,panels(3,idx))-c4_coords_edge(1,i));
                    sgnx=sign(grid(1,panels(3,idx))-c4_coords_edge(1,i));
                    dTheta=deflections_structmesh(5,i);
                    Theta=R(2);
                   %  plot3(c4_coords_edge(1,i)-grid(1,panels(3,idx)),c4_coords_edge(2,i),c4_coords_edge(3,i),'cx')
                    dx2=-sgnx*(dist_x-(dist_x/cos(Theta)*cos(Theta+dTheta)));
                    dz2=-sgnx*(dist_x/cos(Theta)*sin(Theta+dTheta)-dist_x/cos(Theta)*sin(Theta));
                    dTheta=deflections_structmesh(6,i);
                    Theta=R(3);
                    dx3=-sgnx*(dist_x-(dist_x/cos(Theta)*cos(Theta+dTheta)));
                    dy3=-sgnx*(dist_x/cos(Theta)*sin(Theta+dTheta)-dist_x/cos(Theta)*sin(Theta));
                    delta_twist_4=[dx2+dx3;-dy3;dz2];
                    grid(:,panels(3,idx))=grid(:,panels(3,idx))+delta_twist_4;
                end
                le_ctr=1;
                te_ctr=1;
            end
            le_loc_ctr=1;
            te_loc_ctr=1;
            
            Theta=(obj.Theta_r*(1-span_grid(i+1))+obj.Theta_t*span_grid(i+1))*pi/180*0;
            
            a=-obj.dihed*pi/180;
            
            b=Theta;
            
            c=obj.c4_sweep;
            
            
            Lx=[1       0       0
                0   cos(a)  sin(a)
                0   -sin(a) cos(a)];
            
            
            
            Lz=[cos(c) sin(c)   0
                -sin(c) cos(c)  0
                0           0   1];
            
            R=Lz*Lx*[0;Theta;0];
            
            c4_coords_edge(:,i+1)=(obj.xyz(:,1)*0.75+obj.xyz(:,4)*0.25)+span_grid(i+1)*(obj.xyz(:,2)*0.75+obj.xyz(:,3)*0.25-obj.xyz(:,1)*0.75-obj.xyz(:,4)*0.25);
            %c4_coords_edge(2,i)=-c4_coords_edge(2,i);
           %  plot3(c4_coords_edge(1,i+1),c4_coords_edge(2,i+1),c4_coords_edge(3,i+1),'x')
            c4_coords_edge(:,i+1)=c4_coords_edge(:,i+1)+deflections_structmesh(1:3,i+1);
            
           %  plot3(c4_coords_edge(1,i+1),c4_coords_edge(2,i+1),c4_coords_edge(3,i+1),'ro')
            for j=1:n_pan
                idx=obj.panel_start_idx+n_pan*(i-1)+j-1+offset;
                
                grid(:,panels(1,idx))=grid(:,panels(1,idx))+deflections_structmesh(1:3,i+1);
                
                dist_x=abs(grid(1,panels(1,idx))-c4_coords_edge(1,i+1));
                sgnx=sign(grid(1,panels(1,idx))-c4_coords_edge(1,i+1));
                dTheta=deflections_structmesh(5,i+1);
                Theta=R(2);
                dx2=-sgnx*(dist_x-(dist_x/cos(Theta)*cos(Theta+dTheta)));
                dz2=-sgnx*(dist_x/cos(Theta)*sin(Theta+dTheta)-dist_x/cos(Theta)*sin(Theta));
                Theta=R(3);
                dTheta=deflections_structmesh(6,i+1);
                
                dx3=-sgnx*(dist_x-(dist_x/cos(Theta)*cos(Theta+dTheta)));
                dy3=-sgnx*(dist_x/cos(Theta)*sin(Theta+dTheta)-dist_x/cos(Theta)*sin(Theta));
                delta_twist_2=[dx2+dx3;-dy3;dz2];
                
               % plot3(grid(1,panels(1,idx))-c4_coords_edge(1,i),c4_coords_edge(2,i),c4_coords_edge(3,i),'gx')
                
                grid(:,panels(1,idx))=grid(:,panels(1,idx))+delta_twist_2;
                if obj.has_le_cs
                    if le_ctr<=length(obj.n_le_panels)
                        le_loc_ctr=le_loc_ctr+1;
                        if le_loc_ctr>obj.n_le_panels(le_ctr)
                            le_loc_ctr=1;
                            le_ctr=le_ctr+1;
                            grid(:,panels(4,idx))=grid(:,panels(4,idx))+deflections_structmesh(1:3,i+1);
                            dist_x=abs(grid(1,panels(4,idx))-c4_coords_edge(1,i+1));
                            sgnx=sign(grid(1,panels(4,idx))-c4_coords_edge(1,i+1));
                         %     plot3(c4_coords_edge(1,i)-grid(1,panels(4,idx)),c4_coords_edge(2,i),c4_coords_edge(3,i),'cx')
                            dTheta=deflections_structmesh(5,i+1);
                            Theta=R(2);
                            dx2=-sgnx*(dist_x-(dist_x/cos(Theta)*cos(Theta+dTheta)));
                            dz2=-sgnx*(dist_x/cos(Theta)*sin(Theta+dTheta)-dist_x/cos(Theta)*sin(Theta));
                            Theta=R(3);
                            dTheta=deflections_structmesh(6,i+1);
                            
                            dx3=-sgnx*(dist_x-(dist_x/cos(Theta)*cos(Theta+dTheta)));
                            dy3=-sgnx*(dist_x/cos(Theta)*sin(Theta+dTheta)-dist_x/cos(Theta)*sin(Theta));
                            delta_twist_3=[dx2+dx3;-dy3;dz2];
                            
                            grid(:,panels(4,idx))=grid(:,panels(4,idx))+delta_twist_3;
                        end
                    end
                end
                
                if obj.has_te_cs
                    if j>n_pan-sum(obj.n_te_panels)
                        if te_ctr<=length(obj.n_te_panels)
                            te_loc_ctr=te_loc_ctr+1;
                            if te_loc_ctr>obj.n_te_panels(te_ctr)
                                te_loc_ctr=1;
                                te_ctr=te_ctr+1;
                                grid(:,panels(4,idx))=grid(:,panels(4,idx))+deflections_structmesh(1:3,i+1);
                                dist_x=abs(grid(1,panels(4,idx))-c4_coords_edge(1,i+1));
                                sgnx=sign(grid(1,panels(4,idx))-c4_coords_edge(1,i+1));
                            %     plot3(c4_coords_edge(1,i)-grid(1,panels(4,idx)),c4_coords_edge(2,i),c4_coords_edge(3,i),'cx')
                                dTheta=deflections_structmesh(5,i+1);
                                Theta=R(2);
                                dx2=-sgnx*(dist_x-(dist_x/cos(Theta)*cos(Theta+dTheta)));
                                dz2=-sgnx*(dist_x/cos(Theta)*sin(Theta+dTheta)-dist_x/cos(Theta)*sin(Theta));
                                Theta=R(3);
                                dTheta=deflections_structmesh(6,i+1);
                                
                                dx3=-sgnx*(dist_x-(dist_x/cos(Theta)*cos(Theta+dTheta)));
                                dy3=-sgnx*(dist_x/cos(Theta)*sin(Theta+dTheta)-dist_x/cos(Theta)*sin(Theta));
                                delta_twist_3=[dx2+dx3;-dy3;dz2];
                                
                                grid(:,panels(4,idx))=grid(:,panels(4,idx))+delta_twist_3;
                            end
                        end
                    elseif j==n_pan-sum(obj.n_te_panels)
                        grid(:,panels(4,idx))=grid(:,panels(4,idx))+deflections_structmesh(1:3,i+1);
                        dist_x=abs(grid(1,panels(4,idx))-c4_coords_edge(1,i+1));
                        sgnx=sign(grid(1,panels(4,idx))-c4_coords_edge(1,i+1));
                        %   plot3(c4_coords_edge(1,i)-grid(1,panels(4,idx)),c4_coords_edge(2,i),c4_coords_edge(3,i),'cx')
                        dTheta=deflections_structmesh(5,i+1);
                        Theta=R(2);
                        dx2=-sgnx*(dist_x-(dist_x/cos(Theta)*cos(Theta+dTheta)));
                        dz2=-sgnx*(dist_x/cos(Theta)*sin(Theta+dTheta)-dist_x/cos(Theta)*sin(Theta));
                        Theta=R(3);
                        dTheta=deflections_structmesh(6,i+1);
                        
                        dx3=-sgnx*(dist_x-(dist_x/cos(Theta)*cos(Theta+dTheta)));
                        dy3=-sgnx*(dist_x/cos(Theta)*sin(Theta+dTheta)-dist_x/cos(Theta)*sin(Theta));
                        delta_twist_3=[dx2+dx3;-dy3;dz2];
                        
                        grid(:,panels(4,idx))=grid(:,panels(4,idx))+delta_twist_3;
                    end
                end
                if obj.has_sp_cs
                    if j==obj.n_ctr1_panels+sum(obj.n_le_panels)
                         grid(:,panels(4,idx))=grid(:,panels(4,idx))+deflections_structmesh(1:3,i+1);
                        dist_x=abs(grid(1,panels(4,idx))-c4_coords_edge(1,i+1));
                        sgnx=sign(grid(1,panels(4,idx))-c4_coords_edge(1,i+1));
                        %   plot3(c4_coords_edge(1,i)-grid(1,panels(4,idx)),c4_coords_edge(2,i),c4_coords_edge(3,i),'cx')
                        dTheta=deflections_structmesh(5,i+1);
                        Theta=R(2);
                        dx2=-sgnx*(dist_x-(dist_x/cos(Theta)*cos(Theta+dTheta)));
                        dz2=-sgnx*(dist_x/cos(Theta)*sin(Theta+dTheta)-dist_x/cos(Theta)*sin(Theta));
                        Theta=R(3);
                        dTheta=deflections_structmesh(6,i+1);
                        
                        dx3=-sgnx*(dist_x-(dist_x/cos(Theta)*cos(Theta+dTheta)));
                        dy3=-sgnx*(dist_x/cos(Theta)*sin(Theta+dTheta)-dist_x/cos(Theta)*sin(Theta));
                        delta_twist_3=[dx2+dx3;-dy3;dz2];
                        
                        grid(:,panels(4,idx))=grid(:,panels(4,idx))+delta_twist_3;
                    end
                end
                if and(obj.has_sp_cs,~obj.has_te_cs)
                    if j==obj.n_ctr1_panels+sum(obj.n_le_panels)+obj.n_sp_panels
                         grid(:,panels(4,idx))=grid(:,panels(4,idx))+deflections_structmesh(1:3,i+1);
                        dist_x=abs(grid(1,panels(4,idx))-c4_coords_edge(1,i+1));
                        sgnx=sign(grid(1,panels(4,idx))-c4_coords_edge(1,i+1));
                        %   plot3(c4_coords_edge(1,i)-grid(1,panels(4,idx)),c4_coords_edge(2,i),c4_coords_edge(3,i),'cx')
                        dTheta=deflections_structmesh(5,i+1);
                        Theta=R(2);
                        dx2=-sgnx*(dist_x-(dist_x/cos(Theta)*cos(Theta+dTheta)));
                        dz2=-sgnx*(dist_x/cos(Theta)*sin(Theta+dTheta)-dist_x/cos(Theta)*sin(Theta));
                        Theta=R(3);
                        dTheta=deflections_structmesh(6,i+1);
                        
                        dx3=-sgnx*(dist_x-(dist_x/cos(Theta)*cos(Theta+dTheta)));
                        dy3=-sgnx*(dist_x/cos(Theta)*sin(Theta+dTheta)-dist_x/cos(Theta)*sin(Theta));
                        delta_twist_3=[dx2+dx3;-dy3;dz2];
                        
                        grid(:,panels(4,idx))=grid(:,panels(4,idx))+delta_twist_3;
                    end
                end
              
            end
            if (isempty(obj.has_te_cs)) || (obj.has_te_cs==0)
                grid(:,panels(4,idx))=grid(:,panels(4,idx))+deflections_structmesh(1:3,i+1);
                dist_x=abs(grid(1,panels(4,idx))-c4_coords_edge(1,i+1));
                sgnx=sign(grid(1,panels(4,idx))-c4_coords_edge(1,i+1));
                % plot3(c4_coords_edge(1,i)-grid(1,panels(4,idx)),c4_coords_edge(2,i),c4_coords_edge(3,i),'cx')
                dTheta=deflections_structmesh(5,i+1);
                Theta=R(2);
                dx2=-sgnx*(dist_x-(dist_x/cos(Theta)*cos(Theta+dTheta)));
                dz2=-sgnx*(dist_x/cos(Theta)*sin(Theta+dTheta)-dist_x/cos(Theta)*sin(Theta));
                Theta=R(3);
                dTheta=deflections_structmesh(6,i+1);
                
                dx3=-sgnx*(dist_x-(dist_x/cos(Theta)*cos(Theta+dTheta)));
                dy3=-sgnx*(dist_x/cos(Theta)*sin(Theta+dTheta)-dist_x/cos(Theta)*sin(Theta));
                delta_twist_3=[dx2+dx3;-dy3;dz2];
                
                grid(:,panels(4,idx))=grid(:,panels(4,idx))+delta_twist_3;
            end
            %grid(2,:)=-grid(2,:);
        end
        
        function grid=compute_deflected_grid(obj,panels,grid,deflections_structmesh)

            span_grid=0:1/obj.n_span:1;
            n_pan=sum([obj.n_le_panels obj.n_ctr1_panels obj.n_sp_panels obj.n_ctr2_panels  obj.n_te_panels]);
            
            le_ctr=1;
            le_loc_ctr=1;
            sp_ctr=1;
            sp_loc_ctr=1;
            te_ctr=1;
            te_loc_ctr=1;
            
            for i=1:obj.n_span
                %twist set to zerO?
                Theta=(obj.Theta_r*(1-span_grid(i))+obj.Theta_t*span_grid(i))*pi/180*0;
                
                a=obj.dihed*pi/180;
                c=obj.c4_sweep*pi/180;
                
                Lx=[1       0       0
                    0   cos(a)  sin(a)
                    0   -sin(a) cos(a)];
                
                
                Lz=[cos(c) sin(c)   0
                    -sin(c) cos(c)  0
                    0           0   1];
                
                %R is always zero beacuse twist is set to zero
                R=Lz*Lx*[0;Theta;0];
                
                c4_coords_edge(:,i)=(obj.xyz(:,1)*0.75+obj.xyz(:,4)*0.25)+span_grid(i)*(obj.xyz(:,2)*0.75+obj.xyz(:,3)*0.25-obj.xyz(:,1)*0.75-obj.xyz(:,4)*0.25);
              %    hold on
              %    plot3(c4_coords_edge(1,i),c4_coords_edge(2,i),c4_coords_edge(3,i),'x')
                c4_coords_edge(:,i)=c4_coords_edge(:,i)+deflections_structmesh(1:3,i);
                
         %        plot3(c4_coords_edge(1,i),c4_coords_edge(2,i),c4_coords_edge(3,i),'ro')
                for j=1:n_pan
                    idx=obj.panel_start_idx+n_pan*(i-1)+j-1;
                    
                    grid(:,panels(1,idx))=grid(:,panels(1,idx))+deflections_structmesh(1:3,i);
                    
                    dist_x=norm(grid(1:3,panels(1,idx))-c4_coords_edge(1:3,i));
                    sgnx=sign(grid(1,panels(1,idx))-c4_coords_edge(1,i));
                    
                    dTheta=deflections_structmesh(5,i);
                    Theta=R(2);
                    dx2=-sgnx*(dist_x*cos(Theta)-(dist_x*cos(Theta+dTheta)));
                    dz2=-sgnx*(dist_x*sin(Theta+dTheta)-dist_x*sin(Theta));
                    Theta=R(3);
                    
                    dTheta=deflections_structmesh(6,i);
                    dx3=-sgnx*(dist_x*cos(Theta)-(dist_x*cos(Theta+dTheta)));
                    dy3=sgnx*(dist_x*sin(Theta+dTheta)-dist_x*sin(Theta));
                    
                    
                    delta_twist_1=[dx2+dx3;dy3;dz2];
                    grid(:,panels(1,idx))=grid(:,panels(1,idx))+delta_twist_1;
                    
                  %    plot3(grid(1,panels(1,idx)),grid(2,panels(1,idx)),grid(3,panels(1,idx)),'go')
                    if obj.has_le_cs
                        if le_ctr<=length(obj.n_le_panels)
                            le_loc_ctr=le_loc_ctr+1;
                            if le_loc_ctr>obj.n_le_panels(le_ctr)
                                le_loc_ctr=1;
                                le_ctr=le_ctr+1;
                                grid(:,panels(4,idx))=grid(:,panels(4,idx))+deflections_structmesh(1:3,i);
                                dist_x=abs(grid(1,panels(4,idx))-c4_coords_edge(1,i));
                                sgnx=sign(grid(1,panels(4,idx))-c4_coords_edge(1,i));
                            %     plot3(c4_coords_edge(1,i)-grid(1,panels(4,idx)),c4_coords_edge(2,i),c4_coords_edge(3,i),'cx')
                                dTheta=deflections_structmesh(5,i);
                                Theta=R(2);
                                dx2=-sgnx*(dist_x-(dist_x/cos(Theta)*cos(Theta+dTheta)));
                                dz2=-sgnx*(dist_x/cos(Theta)*sin(Theta+dTheta)-dist_x/cos(Theta)*sin(Theta));
                                dTheta=deflections_structmesh(6,i);
                                Theta=R(3);
                                dx3=-sgnx*(dist_x-(dist_x/cos(Theta)*cos(Theta+dTheta)));
                                dy3=sgnx*(dist_x/cos(Theta)*sin(Theta+dTheta)-dist_x/cos(Theta)*sin(Theta));
                                delta_twist_4=[dx2+dx3;dy3;dz2];
                                grid(:,panels(4,idx))=grid(:,panels(4,idx))+delta_twist_4;
                            end
                            
                        end
                    end
                    
                    if obj.has_te_cs
                        if j>n_pan-sum(obj.n_te_panels)
                            if te_ctr<=length(obj.n_te_panels)
                                te_loc_ctr=te_loc_ctr+1;
                                if te_loc_ctr>obj.n_te_panels(te_ctr)
                                    te_loc_ctr=1;
                                    te_ctr=te_ctr+1;
                                    grid(:,panels(4,idx))=grid(:,panels(4,idx))+deflections_structmesh(1:3,i);
                                    dist_x=abs(grid(1,panels(4,idx))-c4_coords_edge(1,i));
                                    sgnx=sign(grid(1,panels(4,idx))-c4_coords_edge(1,i));
                                    dTheta=deflections_structmesh(5,i);
                                    Theta=R(2);
                                   %  plot3(c4_coords_edge(1,i)-grid(1,panels(4,idx)),c4_coords_edge(2,i),c4_coords_edge(3,i),'cx')
                                    dx2=-sgnx*(dist_x-(dist_x/cos(Theta)*cos(Theta+dTheta)));
                                    dz2=-sgnx*(dist_x/cos(Theta)*sin(Theta+dTheta)-dist_x/cos(Theta)*sin(Theta));
                                    dTheta=deflections_structmesh(6,i);
                                    Theta=R(3);
                                    dx3=-sgnx*(dist_x-(dist_x/cos(Theta)*cos(Theta+dTheta)));
                                    dy3=sgnx*(dist_x/cos(Theta)*sin(Theta+dTheta)-dist_x/cos(Theta)*sin(Theta));
                                    delta_twist_4=[dx2+dx3;dy3;dz2];
                                    grid(:,panels(4,idx))=grid(:,panels(4,idx))+delta_twist_4;
                                end
                            end
                        elseif j==n_pan-sum(obj.n_te_panels)
                            grid(:,panels(4,idx))=grid(:,panels(4,idx))+deflections_structmesh(1:3,i);
                            dist_x=abs(grid(1,panels(4,idx))-c4_coords_edge(1,i));
                            sgnx=sign(grid(1,panels(4,idx))-c4_coords_edge(1,i));
                            dTheta=deflections_structmesh(5,i);
                            Theta=R(2);
                            % plot3(c4_coords_edge(1,i)-grid(1,panels(4,idx)),c4_coords_edge(2,i),c4_coords_edge(3,i),'cx')
                            dx2=-sgnx*(dist_x-(dist_x/cos(Theta)*cos(Theta+dTheta)));
                            dz2=-sgnx*(dist_x/cos(Theta)*sin(Theta+dTheta)-dist_x/cos(Theta)*sin(Theta));
                            dTheta=deflections_structmesh(6,i);
                            Theta=R(3);
                            dx3=-sgnx*(dist_x-(dist_x/cos(Theta)*cos(Theta+dTheta)));
                            dy3=sgnx*(dist_x/cos(Theta)*sin(Theta+dTheta)-dist_x/cos(Theta)*sin(Theta));
                            delta_twist_4=[dx2+dx3;dy3;dz2];
                            grid(:,panels(4,idx))=grid(:,panels(4,idx))+delta_twist_4;
                        end
                    end
                    if obj.has_sp_cs
                      if j==obj.n_ctr1_panels+sum(obj.n_le_panels)
                            grid(:,panels(4,idx))=grid(:,panels(4,idx))+deflections_structmesh(1:3,i);
                            dist_x=abs(grid(1,panels(4,idx))-c4_coords_edge(1,i));
                            sgnx=sign(grid(1,panels(4,idx))-c4_coords_edge(1,i));
                            dTheta=deflections_structmesh(5,i);
                            Theta=R(2);
                            % plot3(c4_coords_edge(1,i)-grid(1,panels(4,idx)),c4_coords_edge(2,i),c4_coords_edge(3,i),'cx')
                            dx2=-sgnx*(dist_x-(dist_x/cos(Theta)*cos(Theta+dTheta)));
                            dz2=-sgnx*(dist_x/cos(Theta)*sin(Theta+dTheta)-dist_x/cos(Theta)*sin(Theta));
                            dTheta=deflections_structmesh(6,i);
                            Theta=R(3);
                            dx3=-sgnx*(dist_x-(dist_x/cos(Theta)*cos(Theta+dTheta)));
                            dy3=sgnx*(dist_x/cos(Theta)*sin(Theta+dTheta)-dist_x/cos(Theta)*sin(Theta));
                            delta_twist_4=[dx2+dx3;dy3;dz2];
                            grid(:,panels(4,idx))=grid(:,panels(4,idx))+delta_twist_4;
                       end
                    end
                    if and(obj.has_sp_cs,~obj.has_te_cs)
                      if j==obj.n_ctr1_panels+sum(obj.n_le_panels)+obj.n_sp_panels
                            grid(:,panels(4,idx))=grid(:,panels(4,idx))+deflections_structmesh(1:3,i);
                            dist_x=abs(grid(1,panels(4,idx))-c4_coords_edge(1,i));
                            sgnx=sign(grid(1,panels(4,idx))-c4_coords_edge(1,i));
                            dTheta=deflections_structmesh(5,i);
                            Theta=R(2);
                            % plot3(c4_coords_edge(1,i)-grid(1,panels(4,idx)),c4_coords_edge(2,i),c4_coords_edge(3,i),'cx')
                            dx2=-sgnx*(dist_x-(dist_x/cos(Theta)*cos(Theta+dTheta)));
                            dz2=-sgnx*(dist_x/cos(Theta)*sin(Theta+dTheta)-dist_x/cos(Theta)*sin(Theta));
                            dTheta=deflections_structmesh(6,i);
                            Theta=R(3);
                            dx3=-sgnx*(dist_x-(dist_x/cos(Theta)*cos(Theta+dTheta)));
                            dy3=sgnx*(dist_x/cos(Theta)*sin(Theta+dTheta)-dist_x/cos(Theta)*sin(Theta));
                            delta_twist_4=[dx2+dx3;dy3;dz2];
                            grid(:,panels(4,idx))=grid(:,panels(4,idx))+delta_twist_4;
                       end
                    end
                end
                if  (isempty(obj.has_te_cs))||(obj.has_te_cs==0)
                    grid(:,panels(4,idx))=grid(:,panels(4,idx))+deflections_structmesh(1:3,i);
                    dist_x=abs(grid(1,panels(4,idx))-c4_coords_edge(1,i));
                    sgnx=sign(grid(1,panels(4,idx))-c4_coords_edge(1,i));
                    dTheta=deflections_structmesh(5,i);
                    Theta=R(2);
                   %  plot3(c4_coords_edge(1,i)-grid(1,panels(4,idx)),c4_coords_edge(2,i),c4_coords_edge(3,i),'cx')
                    dx2=-sgnx*(dist_x-(dist_x/cos(Theta)*cos(Theta+dTheta)));
                    dz2=-sgnx*(dist_x/cos(Theta)*sin(Theta+dTheta)-dist_x/cos(Theta)*sin(Theta));
                    dTheta=deflections_structmesh(6,i);
                    Theta=R(3);
                    dx3=-sgnx*(dist_x-(dist_x/cos(Theta)*cos(Theta+dTheta)));
                    dy3=sgnx*(dist_x/cos(Theta)*sin(Theta+dTheta)-dist_x/cos(Theta)*sin(Theta));
                    delta_twist_4=[dx2+dx3;dy3;dz2];
                    grid(:,panels(4,idx))=grid(:,panels(4,idx))+delta_twist_4;
                end
                le_ctr=1;
                te_ctr=1;
            end
            % last spanwise row
            le_loc_ctr=1;
            te_loc_ctr=1;
            
            Theta=(obj.Theta_r*(1-span_grid(i+1))+obj.Theta_t*span_grid(i+1))*pi/180*0;
            
            a=obj.dihed*pi/180;
            
            b=Theta;
            
            c=obj.c4_sweep;
            
            
            Lx=[1       0       0
                0   cos(a)  sin(a)
                0   -sin(a) cos(a)];
            
            
            
            Lz=[cos(c) sin(c)   0
                -sin(c) cos(c)  0
                0           0   1];
            
            R=Lz*Lx*[0;Theta;0];
            
            c4_coords_edge(:,i+1)=(obj.xyz(:,1)*0.75+obj.xyz(:,4)*0.25)+span_grid(i+1)*(obj.xyz(:,2)*0.75+obj.xyz(:,3)*0.25-obj.xyz(:,1)*0.75-obj.xyz(:,4)*0.25);
            
           %  plot3(c4_coords_edge(1,i+1),c4_coords_edge(2,i+1),c4_coords_edge(3,i+1),'x')
            c4_coords_edge(:,i+1)=c4_coords_edge(:,i+1)+deflections_structmesh(1:3,i+1);
           %  plot3(c4_coords_edge(1,i+1),c4_coords_edge(2,i+1),c4_coords_edge(3,i+1),'ro')
            for j=1:n_pan
                idx=obj.panel_start_idx+n_pan*(i-1)+j-1;
                
                grid(:,panels(2,idx))=grid(:,panels(2,idx))+deflections_structmesh(1:3,i+1);
                
                dist_x=abs(grid(1,panels(2,idx))-c4_coords_edge(1,i+1));
                sgnx=sign(grid(1,panels(2,idx))-c4_coords_edge(1,i+1));
                dTheta=deflections_structmesh(5,i+1);
                Theta=R(2);
                dx2=-sgnx*(dist_x-(dist_x/cos(Theta)*cos(Theta+dTheta)));
                dz2=-sgnx*(dist_x/cos(Theta)*sin(Theta+dTheta)-dist_x/cos(Theta)*sin(Theta));
                Theta=R(3);
                dTheta=deflections_structmesh(6,i+1);
                
                dx3=-sgnx*(dist_x-(dist_x/cos(Theta)*cos(Theta+dTheta)));
                dy3=sgnx*(dist_x/cos(Theta)*sin(Theta+dTheta)-dist_x/cos(Theta)*sin(Theta));
                delta_twist_2=[dx2+dx3;dy3;dz2];
                
                
                %plot3(grid(1,panels(2,idx))-c4_coords_edge(1,i),c4_coords_edge(2,i),c4_coords_edge(3,i),'gx')
                
                grid(:,panels(2,idx))=grid(:,panels(2,idx))+delta_twist_2;
                if obj.has_le_cs
                    if le_ctr<=length(obj.n_le_panels)
                        le_loc_ctr=le_loc_ctr+1;
                        if le_loc_ctr>obj.n_le_panels(le_ctr)
                            le_loc_ctr=1;
                            le_ctr=le_ctr+1;
                            grid(:,panels(3,idx))=grid(:,panels(3,idx))+deflections_structmesh(1:3,i+1);
                            dist_x=abs(grid(1,panels(3,idx))-c4_coords_edge(1,i+1));
                            sgnx=sign(grid(1,panels(3,idx))-c4_coords_edge(1,i+1));
                          %     plot3(c4_coords_edge(1,i)-grid(1,panels(3,idx)),c4_coords_edge(2,i),c4_coords_edge(3,i),'cx')
                            dTheta=deflections_structmesh(5,i+1);
                            Theta=R(2);
                            dx2=-sgnx*(dist_x-(dist_x/cos(Theta)*cos(Theta+dTheta)));
                            dz2=-sgnx*(dist_x/cos(Theta)*sin(Theta+dTheta)-dist_x/cos(Theta)*sin(Theta));
                            Theta=R(3);
                            dTheta=deflections_structmesh(6,i+1);
                            
                            dx3=-sgnx*(dist_x-(dist_x/cos(Theta)*cos(Theta+dTheta)));
                            dy3=sgnx*(dist_x/cos(Theta)*sin(Theta+dTheta)-dist_x/cos(Theta)*sin(Theta));
                            delta_twist_3=[dx2+dx3;dy3;dz2];
                            
                            grid(:,panels(3,idx))=grid(:,panels(3,idx))+delta_twist_3;
                        end
                    end
                end
                
                if obj.has_te_cs
                    if j>n_pan-sum(obj.n_te_panels)
                        if te_ctr<=length(obj.n_te_panels)
                            te_loc_ctr=te_loc_ctr+1;
                            if te_loc_ctr>obj.n_te_panels(te_ctr)
                                te_loc_ctr=1;
                                te_ctr=te_ctr+1;
                                grid(:,panels(3,idx))=grid(:,panels(3,idx))+deflections_structmesh(1:3,i+1);
                                dist_x=abs(grid(1,panels(3,idx))-c4_coords_edge(1,i+1));
                                sgnx=sign(grid(1,panels(3,idx))-c4_coords_edge(1,i+1));
                             %    plot3(c4_coords_edge(1,i)-grid(1,panels(3,idx)),c4_coords_edge(2,i),c4_coords_edge(3,i),'cx')
                                dTheta=deflections_structmesh(5,i+1);
                                Theta=R(2);
                                dx2=-sgnx*(dist_x-(dist_x/cos(Theta)*cos(Theta+dTheta)));
                                dz2=-sgnx*(dist_x/cos(Theta)*sin(Theta+dTheta)-dist_x/cos(Theta)*sin(Theta));
                                Theta=R(3);
                                dTheta=deflections_structmesh(6,i+1);
                                
                                dx3=-sgnx*(dist_x-(dist_x/cos(Theta)*cos(Theta+dTheta)));
                                dy3=sgnx*(dist_x/cos(Theta)*sin(Theta+dTheta)-dist_x/cos(Theta)*sin(Theta));
                                delta_twist_3=[dx2+dx3;dy3;dz2];
                                
                                grid(:,panels(3,idx))=grid(:,panels(3,idx))+delta_twist_3;
                            end
                        end
                    elseif j==n_pan-sum(obj.n_te_panels)
                        grid(:,panels(3,idx))=grid(:,panels(3,idx))+deflections_structmesh(1:3,i+1);
                        dist_x=abs(grid(1,panels(3,idx))-c4_coords_edge(1,i+1));
                        sgnx=sign(grid(1,panels(3,idx))-c4_coords_edge(1,i+1));
                        %   plot3(c4_coords_edge(1,i)-grid(1,panels(3,idx)),c4_coords_edge(2,i),c4_coords_edge(3,i),'cx')
                        dTheta=deflections_structmesh(5,i+1);
                        Theta=R(2);
                        dx2=-sgnx*(dist_x-(dist_x/cos(Theta)*cos(Theta+dTheta)));
                        dz2=-sgnx*(dist_x/cos(Theta)*sin(Theta+dTheta)-dist_x/cos(Theta)*sin(Theta));
                        Theta=R(3);
                        dTheta=deflections_structmesh(6,i+1);
                        
                        dx3=-sgnx*(dist_x-(dist_x/cos(Theta)*cos(Theta+dTheta)));
                        dy3=sgnx*(dist_x/cos(Theta)*sin(Theta+dTheta)-dist_x/cos(Theta)*sin(Theta));
                        delta_twist_3=[dx2+dx3;dy3;dz2];
                        
                        grid(:,panels(3,idx))=grid(:,panels(3,idx))+delta_twist_3;
                    end
                end
                if obj.has_sp_cs
                    if j==obj.n_ctr1_panels+sum(obj.n_le_panels)
                        grid(:,panels(3,idx))=grid(:,panels(3,idx))+deflections_structmesh(1:3,i+1);
                        dist_x=abs(grid(1,panels(3,idx))-c4_coords_edge(1,i+1));
                        sgnx=sign(grid(1,panels(3,idx))-c4_coords_edge(1,i+1));
                        %   plot3(c4_coords_edge(1,i)-grid(1,panels(3,idx)),c4_coords_edge(2,i),c4_coords_edge(3,i),'cx')
                        dTheta=deflections_structmesh(5,i+1);
                        Theta=R(2);
                        dx2=-sgnx*(dist_x-(dist_x/cos(Theta)*cos(Theta+dTheta)));
                        dz2=-sgnx*(dist_x/cos(Theta)*sin(Theta+dTheta)-dist_x/cos(Theta)*sin(Theta));
                        Theta=R(3);
                        dTheta=deflections_structmesh(6,i+1);
                        
                        dx3=-sgnx*(dist_x-(dist_x/cos(Theta)*cos(Theta+dTheta)));
                        dy3=sgnx*(dist_x/cos(Theta)*sin(Theta+dTheta)-dist_x/cos(Theta)*sin(Theta));
                        delta_twist_3=[dx2+dx3;dy3;dz2];
                        
                        grid(:,panels(3,idx))=grid(:,panels(3,idx))+delta_twist_3;
                    end
                end
                if and(obj.has_sp_cs,~obj.has_te_cs)
                    if j==obj.n_ctr1_panels+sum(obj.n_le_panels)+obj.n_sp_panels
                        grid(:,panels(3,idx))=grid(:,panels(3,idx))+deflections_structmesh(1:3,i+1);
                        dist_x=abs(grid(1,panels(3,idx))-c4_coords_edge(1,i+1));
                        sgnx=sign(grid(1,panels(3,idx))-c4_coords_edge(1,i+1));
                        %   plot3(c4_coords_edge(1,i)-grid(1,panels(3,idx)),c4_coords_edge(2,i),c4_coords_edge(3,i),'cx')
                        dTheta=deflections_structmesh(5,i+1);
                        Theta=R(2);
                        dx2=-sgnx*(dist_x-(dist_x/cos(Theta)*cos(Theta+dTheta)));
                        dz2=-sgnx*(dist_x/cos(Theta)*sin(Theta+dTheta)-dist_x/cos(Theta)*sin(Theta));
                        Theta=R(3);
                        dTheta=deflections_structmesh(6,i+1);

                        dx3=-sgnx*(dist_x-(dist_x/cos(Theta)*cos(Theta+dTheta)));
                        dy3=sgnx*(dist_x/cos(Theta)*sin(Theta+dTheta)-dist_x/cos(Theta)*sin(Theta));
                        delta_twist_3=[dx2+dx3;dy3;dz2];
                        grid(:,panels(3,idx))=grid(:,panels(3,idx))+delta_twist_3;
                    end
                end
            end
            if (isempty(obj.has_te_cs)) || (obj.has_te_cs==0)
                grid(:,panels(3,idx))=grid(:,panels(3,idx))+deflections_structmesh(1:3,i+1);
                dist_x=abs(grid(1,panels(3,idx))-c4_coords_edge(1,i+1));
                sgnx=sign(grid(1,panels(3,idx))-c4_coords_edge(1,i+1));
                 % plot3(c4_coords_edge(1,i)-grid(1,panels(3,idx)),c4_coords_edge(2,i),c4_coords_edge(3,i),'cx')
                dTheta=deflections_structmesh(5,i+1);
                Theta=R(2);
                dx2=-sgnx*(dist_x-(dist_x/cos(Theta)*cos(Theta+dTheta)));
                dz2=-sgnx*(dist_x/cos(Theta)*sin(Theta+dTheta)-dist_x/cos(Theta)*sin(Theta));
                Theta=R(3);
                dTheta=deflections_structmesh(6,i+1);
                
                dx3=-sgnx*(dist_x-(dist_x/cos(Theta)*cos(Theta+dTheta)));
                dy3=sgnx*(dist_x/cos(Theta)*sin(Theta+dTheta)-dist_x/cos(Theta)*sin(Theta));
                delta_twist_3=[dx2+dx3;dy3;dz2];
                
                grid(:,panels(3,idx))=grid(:,panels(3,idx))+delta_twist_3;
            end
        end
        
        function obj=mirror_control_surfaces(obj)
            if ~isempty(obj.te_device)
                if obj.te_device.is_sym_defl==0
                    
                    obj.te_device.delta=obj.te_device.delta*-1;
                end
            elseif ~isempty(obj.le_device)
                if obj.le_device.is_sym_defl==0
                    obj.le_device.delta=obj.le_device.delta*-1;
                end
            end
            obj=obj.compute_controlsurface_coordinates();
        end
        
        function obj=right_control_surfaces(obj)
            if ~isempty(obj.te_device)
                if obj.te_device.is_sym_defl==0
                    obj.te_device.delta=obj.te_device.delta_l_r(2)*-1;
                end
                obj=obj.compute_controlsurface_coordinates();
            elseif ~isempty(obj.le_device)
                if obj.le_device.is_sym_defl==0
                    obj.le_device.delta=obj.le_device.delta_l_r(2)*-1;
                end
                obj=obj.compute_controlsurface_coordinates();
            end
            
        end
        
        function obj=left_control_surfaces(obj)
            if ~isempty(obj.te_device)
                if obj.te_device.is_sym_defl==0
                    obj.te_device.delta(1)=obj.te_device.delta_l_r(1);
                end
                obj=obj.compute_controlsurface_coordinates();
            elseif ~isempty(obj.le_device)
                if obj.le_device.is_sym_defl==0
                    obj.le_device.delta(1)=obj.le_device.delta_l_r(1);
                end
                obj=obj.compute_controlsurface_coordinates();
            end
        end
        
        function obj=compute_grid(obj,x_max,y_max,n_x_min,wake)
            
            grid3D=1;
            
            % runs if nr of spanwise panels is NOT specified
             if isempty(obj.n_span)
            % number of spanwise panels
                obj.n_span=ceil(obj.b/y_max);
             end
            % number of chordwise panels (accounting only for wing
            % segment geometry, no control surfaces yet)
            if isempty(obj.n_chord)
                obj.n_chord=ceil(norm((obj.xyz_fixed(:,4)-obj.xyz_fixed(:,1))*0.5+(obj.xyz_fixed(:,3)-obj.xyz_fixed(:,2))*0.5)/x_max);
                x_max=norm((obj.xyz_fixed(:,4)-obj.xyz_fixed(:,1))*0.5+(obj.xyz_fixed(:,3)-obj.xyz_fixed(:,2))*0.5)/obj.n_chord;
                 if obj.n_chord<n_x_min
                    obj.n_chord = n_x_min;
                    x_max=norm((obj.xyz_fixed(:,4)-obj.xyz_fixed(:,1))*0.5+(obj.xyz_fixed(:,3)-obj.xyz_fixed(:,2))*0.5)/obj.n_chord;
                 end

                % runs if nr of spanwise panels is specified
            else
                y_max = obj.b/obj.n_span;
                x_max = norm((obj.xyz_fixed(:,4)-obj.xyz_fixed(:,1))*0.5+(obj.xyz_fixed(:,3)-obj.xyz_fixed(:,2))*0.5)/obj.n_chord;
            end
            % Changing x_max to comply with min requirement of
            % chordwise panels (if req is not met)
          
            
            % compute span spacing vector
            %% TODO: optional use nonlinear cosine distribution of panels
            
            % vector from 0 to 1 of points which will serve as "span
            % stations"
            span_spacing=0:1/obj.n_span:1;
            
            % vector from 0 to 1 of points which will serve as "chord
            % stations"
            chord_spacing_ctr=0:1/obj.n_chord:1;
            
            % Initializing a whole bunch of array that will be useful later
            c_c_arr = zeros(1,length(span_spacing));
            rel_start_sp = zeros(1,length(span_spacing));
            rel_end_sp = zeros(1,length(span_spacing));
            xsi_sp_TE_ctr2_arr = ones(1,length(span_spacing));
            taper_correction = zeros(1,length(span_spacing));
            xsi_sp_LE_arr = zeros(1,length(span_spacing));
            te_sp_flap_overlap = zeros(1,length(span_spacing));
            le_sp_flap_overlap = zeros(1,length(span_spacing));
            overlap_flag = 0;
            
            % ctr2 zone is assumed to exist. If there is any overlap
            % between SP and TE surfaces at any span station, it will be
            % set to zero.
            ctr_2_exists = 1;
                
            
            % looping once through the whole span to calculate the current
            % chord at each span location and the taper correction at each
            % span location
            for i=1:length(span_spacing)
                % computes chord of current spanwise location
                c_c_arr(i)=obj.c_r*(1-span_spacing(i))+obj.c_t*span_spacing(i);
                % taper correction is used whenever the CS are tapered.
                taper_correction(i) = obj.c_r / c_c_arr(i);
            end
            
            
            % if spoiler is present, the area between LE and TE CS needs to
            % be split into 2 sections, one before the SP, and one after.
            if (obj.has_sp_cs)
                
                % if loop below allows to calculate arrays containing the
                % LE and TE xsi coordinates of the spoilers. They are
                % stored in arrays since will vary at each span station
                % (Spoiler have variable orientation = 2 different xsi for
                % inner and outer LE, which means different intermediate
                % LE xsi's at the intermediate span stations.
                    
                % runs if we are reading in-house XML file

                % saves inner and outer LE xsi coords of the spoiler
                % (as specified in the XML)
                xsi_sp_LE_arr(1) = obj.sp_device.xsi_LE_inner;
                xsi_sp_LE_arr(end) = obj.sp_device.xsi_LE_outer;

                % calculates chord between LE of wing-segment and LE of
                % spoiler for innermost and outermost span stations
                c_sp_LE_in = obj.c_r * obj.sp_device.xsi_LE_inner;
                c_sp_LE_out = obj.c_t * obj.sp_device.xsi_LE_outer;

                %calculates total difference between aformetioned
                %chords (gradient of change).
                c_sp_LE_delta = abs(c_sp_LE_in-c_sp_LE_out);

                % calculates all chords between LE of wing-segment and LE of
                % spoiler for intermediate span stations
                c_sp_LE_intermediate = c_sp_LE_in - (c_sp_LE_delta * span_spacing(2:end-1));

                % calculates spoiler LE xsi coords for all intermediate
                % span stations
                xsi_sp_LE_arr(2:end-1) = rdivide(c_sp_LE_intermediate,c_c_arr(2:end-1));

                % calculates spoiler TE xsi coords for all span
                % stations. It assumes a tapered spoiler, but
                % corrections are introduced later on to obtain the
                % correct value in case of untapered spoiler.
                xsi_sp_TE_arr = xsi_sp_LE_arr + obj.sp_device.hinge;
                
                % checks whether corrections due to untapered spoiler are
                % necessary
                
                % if tapered, the code below does not alter relevant values
                if obj.sp_device.is_tapered == 1
                    sp_taper_correction = ones(1,length(span_spacing));
                    xsi_sp_TE_ctr2_arr = times(xsi_sp_TE_ctr2_arr, xsi_sp_TE_arr);
                    
                % if untapered, the actual spoiler TE xsi coords are
                % calculated by applying a correction factor
                % (taper_correction)
                else
                    sp_taper_correction = taper_correction;
                    xsi_sp_TE_ctr2_arr = xsi_sp_LE_arr + times((xsi_sp_TE_arr-xsi_sp_LE_arr), sp_taper_correction); 
                end
                
                % looping over all span sections to calculate dimensions
                % and chords of areas which will be divided into panels/
                
                % rel_start_sp is the relative chord (min 0, max 1) between
                % the TE of the LE surface and the LE of the spoilers. It is
                % calculated at each span station (array)
                for i=1:length(span_spacing)
                    
                    % if LE surfaces are present, rel_start_sp is affected.
                    % Corrections are applied in case the LE surface is
                    % untapered.
                    if(obj.has_le_cs)
                        if obj.le_device.is_tapered==1
                            LE_correction = sum(obj.le_device.hinge);
                        elseif obj.le_device.is_tapered==0
                            LE_correction = sum(obj.le_device.hinge)*obj.c_r/c_c_arr(i);
                        end
                        rel_start_sp(i) = (xsi_sp_LE_arr(i) - LE_correction);
                    else
                        rel_start_sp(i) = xsi_sp_LE_arr(i);
                    end
                
                    
                    % same as above, but for rel_end_sp, which is
                    % influenced by TE surfaces.
                    
                    % rel_end_sp is the relative chord (min 0, max 1) between
                    % the TE of the spoiler and the LE of the TE surface. It is
                    % calculated at each span station (array)
                    if (obj.has_te_cs)
                        
                        % if LE surfaces are present, rel_end_sp is affected.
                        % Corrections are applied in case the TE surface is
                        % untapered.
                        if obj.te_device.is_tapered==1
                            TE_correction = sum(obj.te_device.hinge);
                        elseif obj.te_device.is_tapered==0
                            TE_correction = sum(obj.te_device.hinge)*obj.c_r/c_c_arr(i);
                        end
                        
                        % corrections are also necessary in case the
                        % spoiler is untapered.
                        if obj.sp_device.is_tapered==1
                            rel_end_sp(i) = (1 - xsi_sp_TE_arr(i) - TE_correction);
                        elseif obj.sp_device.is_tapered==0
                            rel_end_sp(i) = (1 - (xsi_sp_LE_arr(i) + obj.sp_device.hinge*obj.c_r/c_c_arr(i)) - TE_correction);
                        end
                        
                    % in case TE surfaces are not present. Corrections in
                    % case of untaepered spoilers still need to be applied
                    else
                        if obj.sp_device.is_tapered==1
                            rel_end_sp(i) = (1 - xsi_sp_TE_arr(i));
                        elseif obj.sp_device.is_tapered==0
                            rel_end_sp(i) = 1 - (xsi_sp_LE_arr(i) + obj.sp_device.hinge*obj.c_r/c_c_arr(i));
                        end
                    end
                    
                    if (obj.has_te_cs)
                        % checks whether TE of spoiler overlaps LE of TE
                        % control surface (at current span station)
                        if 1 - xsi_sp_TE_ctr2_arr(i) - TE_correction < 0
                            %if overlap occurs, overlap flag of current span
                            %station "i" is switched from 0 to 1.
                            te_sp_flap_overlap(i) = 1;

                            % overlap will either span the entire surfaces, or
                            % it will span only a portion. In the first case,
                            % ctr2 does not exist. In the second, we end up
                            % with triangular elements, so an error message is
                            % thrown. Hence ctr2 has no reason to exist in case
                            % of ANY overlap at ANY span station.
                            ctr_2_exists = 0;

                            % LE of spoiler can only overlap TE surface if the
                            % TE of the spoiler does so as well

                            if 1 - xsi_sp_LE_arr(i) - TE_correction < 0
                                le_sp_flap_overlap(i) = 1;
                            end
                        end
                    end
                end
                
                
                %ctr stands for center. Center 1 and 2 are the regions
                %respectively between LE surface and spoiler, and between
                %spoiler and TE surface. In case spoilers are absent, only
                %ctr1 exists.
                
                % since rel_start_sp and rel_end_sp are realative chords
                % (between 0 and 1), they represent the percentage of
                % segment taken up by ctr1 and ctr2 respectively.
                % Multiplying them by the nr of chord panels we would have
                % assigned to the whole segment surface gives us the panels
                % to be assigned to ctr1 and ctr2.
                n_chord_panels_ctr_1 = ceil(rel_start_sp(1) * obj.n_chord);
                n_chord_panels_ctr_2 = ceil((rel_end_sp(1)) * obj.n_chord * ctr_2_exists);
                
                obj.n_ctr1_panels = n_chord_panels_ctr_1;
                obj.n_ctr2_panels = n_chord_panels_ctr_2;
                
                % relative spacing of chord-wise points WITHIN REFERENCE
                % FRAME OF chord_spacing_ctr_1
                chord_spacing_ctr_1 = linspace(0,1,n_chord_panels_ctr_1+1);
                
                % calculates various grid size related quantities,
                % depending on whether ctr2 exists.
                if ctr_2_exists
                    chord_spacing_ctr_2 = linspace(0,1,n_chord_panels_ctr_2+1);
                    n_chordwise_points = length(chord_spacing_ctr_1) + length(chord_spacing_ctr_2);
                    n_chordwise_panels = n_chordwise_points -2;
                    size_grid_ctr_1=length(chord_spacing_ctr_1)*length(span_spacing);
                    size_grid_ctr_2=length(chord_spacing_ctr_2)*length(span_spacing);
                    
                else
                    n_chordwise_points = length(chord_spacing_ctr_1);
                    n_chordwise_panels = n_chordwise_points -1;
                    size_grid_ctr_1=length(chord_spacing_ctr_1)*length(span_spacing);
                    size_grid_ctr_2= 0;
                end
                
                
                % all the code below (until the else which runs if no
                % spoielrs are present) is related to the detection and
                % handling of overlaps between the spoiler and trailing
                % edge surface.
                
                % creating a new region (like ctr1 or ctr2) to handle
                % overlap would have resulted in a lot of extra code being
                % added to handle exceptions. It was simpler to treat the
                % ovrerlap region as a "special region" of the flaps. While
                % less elegant, this solution is simpler and easier to
                % implement.
                
                if sum(te_sp_flap_overlap) > 0
                    
                    % checks whether overlap occurs at all span stations.
                    % If it doesn't, we will end up with undesireable
                    % triangular elements. A warning message is raised, and
                    % the user is prompted to modify the XML/CPACS in order
                    % to guarantee a full overlap.
                    if sum(te_sp_flap_overlap) ~= length(te_sp_flap_overlap)
                        msg = ['Trailing edge of the Spoiler does not fully overlap with the leading edge of the flap'...
                        'This would result in undesireable triangular elements. Please modify the XML or CPACS file to ensure a full overlap'];
                        error(msg);
                    end
                    
                    % if the leading of the spoiler is also partially
                    % overlapping the flap, we might have triangular
                    % elements
                    if sum(le_sp_flap_overlap) > 0
                        
                        % if the overlap between the LE of the spoiler and
                        % the flap is only partial, triangular elements
                        % will result. Like above, a similar error message
                        % is raised
                        if sum(le_sp_flap_overlap) < length(le_sp_flap_overlap)
                            msg = ['Leading edge of the Spoiler does not fully overlap with the leading edge of the flap'...
                        'This would result in undesireable triangular elements. Please modify the XML or CPACS file to either ensure a full'...
                        'overlap, or no overlap at all'];
                        error(msg);
                        
                        % if the overlap between the LE of the spoiler and
                        % the flap is full, a warning message is raised,
                        % since this configuration is unusual and the user
                        % might have set it up by accident
                        else
                            msg = ['The spoiler is completely overlapping the TE surface (both the LE and TE of the spoiler are behind the LE of the flap'...
                                'Are you sure that this is what you were trying to design?'];
                            warning(msg);
                        end
                    end
                    
                    % if the code has run this far, no errors were raised,
                    % hence the overlap can be considered "proper"
                    
                    overlap_flag = 1;
                    
                    
                end
                
                
            % if no spoilers are present chord_spacing_ctr_1 is used as the
            % default spacing, and chord_spacing_ctr_2 is set to zero.
            
            % code below is basically snippets of the code that ran in case
            % spoilers were present, only it accounts for the fact that
            % none are present. All the comments made to code sections
            % above apply to the equivalent code snippets below.
            else
                rel_start_sp = ones(1,length(span_spacing));
                
                for i=1:length(span_spacing)
                    if(obj.has_le_cs)
                        if obj.le_device.is_tapered==1
                            LE_correction = sum(obj.le_device.hinge);
                        elseif obj.le_device.is_tapered==0
                            LE_correction = sum(obj.le_device.hinge)*obj.c_r/c_c_arr(i);
                        end
                        rel_start_sp(i) = rel_start_sp(i) - LE_correction;
                    end

                    if (obj.has_te_cs)
                        if obj.te_device.is_tapered==1
                            TE_correction = sum(obj.te_device.hinge);
                        elseif obj.te_device.is_tapered==0
                            TE_correction = sum(obj.te_device.hinge)*obj.c_r/c_c_arr(i);
                        end
                        rel_start_sp(i) = rel_start_sp(i) - TE_correction;
                    end
                    
                    n_chord_panels_ctr_1 = ceil(rel_start_sp(1) * obj.n_chord);

                    chord_spacing_ctr_1 = linspace(0,1,n_chord_panels_ctr_1+1);

                    n_chordwise_points = length(chord_spacing_ctr_1);
                    n_chordwise_panels = n_chord_panels_ctr_1;
                    obj.n_ctr1_panels = n_chord_panels_ctr_1;

                    size_grid_ctr_1=length(chord_spacing_ctr_1)*length(span_spacing);

                    ctr_2_exists = 0;
                    size_grid_ctr_2= 0;
                    obj.n_ctr2_panels=[];
                end
            end

            % initialize spacing vectors specifically for control surfaces.
            % functionalities are identical to span and chord spacing of
            % main wing, but these are in the reference system of their
            % respective control surface
            chord_spacing_le = [];
            chord_spacing_te = [];
            chord_spacing_sp = [];
            
            
            % initialize grid size of control surfaces
            size_grid_le=0;
            size_grid_te=0;
            size_grid_sp=0;
            
            
            % initialize nr of panels for control surfaces
            obj.n_le_panels=[];
            obj.n_te_panels=[];
            obj.n_sp_panels=[];
            
            
            if(obj.has_sp_cs)
                
                n_sp_devices=length(obj.xyz_sp_device(:,1,1));
                
                % loops over all SP CS and uses same criteria as in main
                % wing to find nr of panels (chordwise) necessary to
                % represent the spoiler
                for j=1:n_sp_devices
                    obj.n_sp_panels(j)=ceil(norm((obj.xyz_sp_device(j,:,4)-obj.xyz_sp_device(j,:,1))*0.5+(obj.xyz_sp_device(j,:,3)-obj.xyz_sp_device(j,:,2))*0.5)/x_max);
                    
                    % resizing of nr of panels needed for spoiler is done
                    % looking at the root of the segment.
                    if overlap_flag == 1
                        % ratio_sp_overlap is the ratio of the
                        % non-overlapping spoiler root chord over the whole
                        % spoiler root chord.
                        original_n_sp_panels = obj.n_sp_panels(j);
                        ratio_sp_overlap = (obj.sp_device(j).hinge -(sum(obj.te_device.hinge)-(1-xsi_sp_TE_ctr2_arr(1))))/obj.sp_device(j).hinge;
                        obj.n_sp_panels(j) = ceil(obj.n_sp_panels(j) * ratio_sp_overlap);
                    end
                end
                
                % updates grid sizes of SP surface and main wing to account
                % for presence of panels within the CS.
                size_grid_sp=length(span_spacing)*(sum(obj.n_sp_panels)+length(obj.n_sp_panels));
                n_chordwise_points=n_chordwise_points+sum(obj.n_sp_panels)+length(obj.n_sp_panels);
                n_chordwise_panels=n_chordwise_panels+sum(obj.n_sp_panels);
            end
            
            
            
            % runs if segment has a LE CS
            if(obj.has_le_cs)
                
                % calculates taper correction due to presence of leading
                % edge surfaces (to be used later). This one is actually
                % applied when the surface IS TAPERED
                if obj.le_device.is_tapered == 0
                    le_taper_correction = ones(1,length(span_spacing));
                else
                    le_taper_correction = taper_correction;
                end
                
                n_le_devices=length(obj.xyz_le_device(:,1,1));
                
                % loops over all LE CS and uses same criteria as in main
                % wing to find nr of panels (chordwise) necessary to
                % represent the leading edge
                for j=1:n_le_devices
                    obj.n_le_panels(j)=ceil(norm((obj.xyz_le_device(j,:,4)-obj.xyz_le_device(j,:,1))*0.5+(obj.xyz_le_device(j,:,3)-obj.xyz_le_device(j,:,2))*0.5)/x_max);
                end
                
                % updates grid sizes of LE surface and main wing to account
                % for presence of panels within the CS.
                size_grid_le=length(span_spacing)*(sum(obj.n_le_panels)+length(obj.n_le_panels));
                n_chordwise_points=n_chordwise_points+sum(obj.n_le_panels)+length(obj.n_le_panels);
                n_chordwise_panels=n_chordwise_panels+sum(obj.n_le_panels);  
                
            % still defines a LE taper correction even if there are no LE surfaces.
            % Done in order to avoid errors later on when the variable is called.
            else
                le_taper_correction = ones(1,length(span_spacing));
            end
            
            % runs if segment has a TE CS
            if(obj.has_te_cs)
                
                % taper corrections to be applied in case the TE device is UNTAPERED
                if obj.te_device.is_tapered == 1
                    te_taper_correction = ones(1,length(span_spacing));
                else
                    te_taper_correction = taper_correction;
                end
                
                n_te_devices=length(obj.xyz_te_device(:,1,1));
                n_sp_overlap_panels = zeros(1,n_te_devices);
                n_te_free_panels = zeros(1,n_te_devices);
                
                % loops over all TE CS and uses same criteria as in main
                % wing to find nr of panels (chordwise) necessary to
                % represent the trailing edge
                for j=1:n_te_devices
                    %% Problem: ceil is unstable 4.000000 will be 5
                    %% set proper limit
                    n=norm((obj.xyz_te_device(j,:,4)-obj.xyz_te_device(j,:,1))*0.5+(obj.xyz_te_device(j,:,3)-obj.xyz_te_device(j,:,2))*0.5)/x_max;
                    if ceil(n)-n>0.9
                        obj.n_te_panels(j)=ceil(norm((obj.xyz_te_device(j,:,4)-obj.xyz_te_device(j,:,1))*0.5+(obj.xyz_te_device(j,:,3)-obj.xyz_te_device(j,:,2))*0.5)/x_max)-1;
                    else
                        obj.n_te_panels(j)=ceil(norm((obj.xyz_te_device(j,:,4)-obj.xyz_te_device(j,:,1))*0.5+(obj.xyz_te_device(j,:,3)-obj.xyz_te_device(j,:,2))*0.5)/x_max);
                    end
                    
                    
                    if overlap_flag == 1
                        % ratio_sp_overlap is the ratio of the non-overlapping spoiler root chord over the whole
                        % spoiler root chord.
                        n_sp_overlap_panels(j) = ceil(original_n_sp_panels * (1-ratio_sp_overlap));
                        
                        % ratio_flap_freeis the ratio of the non-overlapping spoiler root chord over the whole
                        % spoiler root chord.
                        ratio_flap_free = (1-xsi_sp_TE_ctr2_arr(1))/obj.te_device(j).hinge;
                        n_te_free_panels(j) = ceil(obj.n_te_panels(j)*ratio_flap_free);
                        obj.n_te_panels(j) = n_te_free_panels(j) + n_sp_overlap_panels(j);
                        obj.n_te_panels_overlap = n_sp_overlap_panels;
                        obj.n_te_panels_free = n_te_free_panels;
                    end
                    
                end
                
                % updates grid sizes of TE surface and main wing to account
                % for presence of panels within the CS.
                size_grid_te=length(span_spacing)*(sum(obj.n_te_panels)+length(obj.n_te_panels));
                n_chordwise_points=n_chordwise_points+sum(obj.n_te_panels)+length(obj.n_te_panels);
                n_chordwise_panels=n_chordwise_panels+sum(obj.n_te_panels);
                
            end
            
            % initializes all grid variables
            grid_len = size_grid_ctr_1 + size_grid_ctr_2 + size_grid_le + size_grid_te + size_grid_sp;
            grid=zeros(3,grid_len);
            grid_flat=zeros(3,grid_len);
            grid_upper=zeros(3,grid_len);
            grid_lower=zeros(3,grid_len);
            te_idx=zeros(1,grid_len);
            
            % initializes counter of grid points
            k=1;
            
            % iterates over each spanwise station
            for i=1:length(span_spacing)
                
                % computes chord of current spanwise location
                c_c=obj.c_r*(1-span_spacing(i))+obj.c_t*span_spacing(i);
                
                % computes chord of LE CS relative to ROOT CHORD, not the
                % current one
                c_le_device_cs=[0 cumsum(obj.c_le_device)]/obj.c_r;
                
                % initializes current chord position (0 to 1) wrt LE CS
                % reference frame
                relative_chord_pos_le=0;
                
                % initializes previous chord position within LE CS and
                % previous skeleton point.
                relative_chord_pos_prv=0;
                skeleton_point_prv=0;
                
                % runs if wing segment has LE CS
                if(obj.has_le_cs)
                    
                    % iterate over all LE CS
                    for j=1:n_le_devices
                        
                        % compute chord spacing vector (define coords from
                        % 0 to 1 of each chord station)
                        chord_spacing_le=0:1/ obj.n_le_panels(j):1;
                        
                        % compute xyz coords of r1 and r2. They are the
                        % points resulting from the intersection between
                        % the Leading and Trailing edges of the LEADING
                        % EDGE CONTROL SURFACE and the current spanwise
                        % station we are looking at.
                        r1=obj.xyz_le_device(j,:,1)+span_spacing(i)*(obj.xyz_le_device(j,:,2)-obj.xyz_le_device(j,:,1));
                        r2=obj.xyz_le_device(j,:,4)+span_spacing(i)*(obj.xyz_le_device(j,:,3)-obj.xyz_le_device(j,:,4));
                        
                        % iterate over all chord stations within current
                        % span stations within LE CS.
                        for ii=1:length(chord_spacing_le)
                            
                            % if current point is not the starting point,
                            % previous point is set, else both prv and
                            % currents points are zero.
                            if ii>1
                                relative_chord_pos_prv=relative_chord_pos_le;
                                skeleton_point_prv=skeleton_point2;
                            end
                            
                            
                            % if CS LE is tapered, the relative chord position of the current point is calulcated
                            % by taking the starting point of the current LE CS (if no slotted LE CS are present it is zero),
                            % then adding the chord spacing of the current point (0 to 1 in CS frame) multiplied by the
                            % chord of the current LE device and then scaled with the root chord of the wing segment.
                            
                            % if CS LE is NOT tapered, the relative chord position of the current point is calulcated
                            % by taking the starting point of the current LE CS (if no slotted LE CS are present it is zero) multiplied with the ratio between the
                            % root chord and current chord, then adding the chord spacing of the current point (0 to 1 in CS frame) multiplied by the
                            % chord of the current LE device and then scaled with the current chord of the wing segment
                            if obj.le_device.is_tapered==1
                                relative_chord_pos_le=c_le_device_cs(j)+chord_spacing_le(ii)*obj.c_le_device(j)/obj.c_r;
                            else
                                relative_chord_pos_le=c_le_device_cs(j)*obj.c_r/c_c+chord_spacing_le(ii)*obj.c_le_device(j)/c_c;
                            end
                            
                            % computes skeleton points
                            skeleton_point=obj.compute_skeleton_point(relative_chord_pos_le,span_spacing(i));
                            skeleton_point2=obj.compute_skeleton_point_new(relative_chord_pos_le,relative_chord_pos_prv,skeleton_point_prv, span_spacing(i));

                            skeleton_point=skeleton_point2;


                            % calculates xyz of current grid point (current
                            % chord station and current span station on
                            % current LE CS)
                            grid(:,k)=r1+chord_spacing_le(ii)*(r2-r1);
                            grid_flat(:,k)=grid(:,k);
                            
                            if grid3D==1
                                [lower_point,upper_point]=obj.compute_thickness_point(relative_chord_pos_le,span_spacing(i));
                                grid_upper(:,k)=grid(:,k);
                                grid_lower(:,k)=grid(:,k);
                                grid_upper(2,k)=grid(2,k)-upper_point*sind(obj.dihed);
                                grid_upper(3,k)=grid(3,k)+upper_point*cosd(obj.dihed);
                                grid_lower(2,k)=grid(2,k)-lower_point*sind(obj.dihed);
                                grid_lower(3,k)=grid(3,k)+lower_point*cosd(obj.dihed);
                            end
                            grid(2,k)=grid(2,k)-skeleton_point*sind(obj.dihed);
                            grid(3,k)=grid(3,k)+skeleton_point*cosd(obj.dihed);
                            
                            % increases grid point counter by 1
                            k=k+1;
                        end
                    end
                end
                
                % Once all grid points at current span station have been
                % added all along the LE CS chord, we move to the main wing

                % same r1 and r2 concept as the one described in LE CS.
                % These however refer to the main wing
                
                r1_ext=obj.xyz(:,1)+span_spacing(i)*(obj.xyz(:,2)-obj.xyz(:,1));
                r2_ext=obj.xyz(:,4)+span_spacing(i)*(obj.xyz(:,3)-obj.xyz(:,4));
                
                % iterate over all chord stations
                for j=1:length(chord_spacing_ctr_1)
                    
                    % calculates relative chord position of chord-wise
                    % points at current span station within ctr1 zone.
                    % Basically defines points according to
                    % chord_spacing_ctr_1 within ctr1, so between the TE of
                    % the LE device and the LE of the spoiler. Taper
                    % corrections are applied both explicitly
                    % (le_taper_correction below) and implicitly (within
                    % rel_start_sp)
                    relative_chord_pos=sum(obj.c_le_device)/(le_taper_correction(i)*c_c) +chord_spacing_ctr_1(j)*rel_start_sp(i);
                    
                    % Stores previous chord position and accounts for taper
                    % of TE CS when looking at chord position on the main
                    % wing
                    if j>1
                        relative_chord_pos_prv=sum(obj.c_le_device)/(le_taper_correction(i)*c_c)+chord_spacing_ctr_1(j-1)*rel_start_sp(i);
                        skeleton_point_prv=skeleton_point2;
                    end
                           
                    skeleton_point=obj.compute_skeleton_point(relative_chord_pos,span_spacing(i));
                    skeleton_point2=obj.compute_skeleton_point_new(relative_chord_pos,relative_chord_pos_prv,skeleton_point_prv, span_spacing(i));

                    skeleton_point=skeleton_point2;


                    % calculates xyz of current grid point (current
                    % chord station and current span station on
                    % main wing)
                    
                    grid(:,k)=r1_ext+(relative_chord_pos) * (r2_ext - r1_ext);
                    
                    grid_flat(:,k)=grid(:,k);
                    if grid3D==1
                       [lower_point,upper_point]=obj.compute_thickness_point(relative_chord_pos,span_spacing(i));
                       grid_upper(:,k)=grid(:,k);
                       grid_lower(:,k)=grid(:,k);
                       grid_upper(2,k)=grid(2,k)-upper_point*sind(obj.dihed);
                       grid_upper(3,k)=grid(3,k)+upper_point*cosd(obj.dihed);
                       grid_lower(2,k)=grid(2,k)-lower_point*sind(obj.dihed);
                       grid_lower(3,k)=grid(3,k)+lower_point*cosd(obj.dihed);
                    end
                    grid(2,k)=grid(2,k)-skeleton_point*sind(obj.dihed);
                    grid(3,k)=grid(3,k)+skeleton_point*cosd(obj.dihed);
                    
                    % increases grid point counter by 1
                    k=k+1;
                end
                
                % runs if wing segment has spoiler
                if(obj.has_sp_cs)
                    
                    % iterate over all spoilers
                    for j=1:n_sp_devices
                        
                        % compute chord spaciing vector (define coords from
                        % 0 to 1 of each chord station)
                        chord_spacing_sp=0:1/obj.n_sp_panels(j):1;
                        
                        % same as r1 and r2 coords previous illustrated,
                        % only for spoilers
                        r1=obj.xyz_sp_device(j,:,1)+span_spacing(i)*(obj.xyz_sp_device(j,:,2)-obj.xyz_sp_device(j,:,1));
                        r2=obj.xyz_sp_device(j,:,4)+span_spacing(i)*(obj.xyz_sp_device(j,:,3)-obj.xyz_sp_device(j,:,4));
                        
                        % iterate over all chord stations 
                        for ii=1:length(chord_spacing_sp)
                            
                            % saves prev relative chord pos and skeleton point
                            if ii>1
                                relative_chord_pos_prv=relative_chord_pos;
                                skeleton_point_prv=skeleton_point2;
                            end
                            
                            % relative chord pos of chord-wise points at
                            % current span station within of spoiler
                            % surface. Contained between LE of spoiler and
                            % TE of spoiler.
                            
                            if overlap_flag == 0
                                relative_chord_pos = xsi_sp_LE_arr(i) + chord_spacing_sp(ii) * obj.sp_device.hinge * sp_taper_correction(i);
                                grid(:,k)=r1+chord_spacing_sp(ii)*(r2-r1);
                            else
                                % te_surface_LE is the flap leading edge xsi coordinate at the current span
                                % station (code snippet taken from TE surface if loop below)
                                te_surface_LE = (1-obj.te_device.hinge*te_taper_correction(i));
                                
                                % spoiler is limited to region between LE of spoiler and LE of flap
                                relative_chord_pos = xsi_sp_LE_arr(i) + (te_surface_LE - xsi_sp_LE_arr(i))*chord_spacing_sp(ii);
                                
                                grid(:,k)=r1+chord_spacing_sp(ii)*(te_surface_LE - xsi_sp_LE_arr(i))/(obj.sp_device.hinge * sp_taper_correction(i))*(r2-r1);
                            end
                            
                            skeleton_point=obj.compute_skeleton_point(relative_chord_pos,span_spacing(i));  
                            skeleton_point2=obj.compute_skeleton_point_new(relative_chord_pos,relative_chord_pos_prv,skeleton_point_prv,span_spacing(i));  

                            skeleton_point=skeleton_point2;
                                
                            
                            grid_flat(:,k)=grid(:,k);
                            if grid3D==1
                                [lower_point,upper_point]=obj.compute_thickness_point(relative_chord_pos,span_spacing(i));
                                grid_upper(:,k)=grid(:,k);
                                grid_lower(:,k)=grid(:,k);
                                grid_upper(2,k)=grid(2,k)-upper_point*sind(obj.dihed);
                                grid_upper(3,k)=grid(3,k)+upper_point*cosd(obj.dihed);
                                grid_lower(2,k)=grid(2,k)-lower_point*sind(obj.dihed);
                                grid_lower(3,k)=grid(3,k)+lower_point*cosd(obj.dihed);
                            end
                            grid(2,k)=grid(2,k)-skeleton_point*sind(obj.dihed);
                            grid(3,k)=grid(3,k)+skeleton_point*cosd(obj.dihed);
                            k=k+1;
                        end
                    end
                end
        
                % Once all grid in the spoilers have been
                % added all along the spoiler chord, we move to the section
                % of the wing between the spoiler and the control surface,
                % but only if this section is actually present. If the
                % spoiler overlaps a TE CS, or if the spoiler is absent,
                % we go directly to the TE
                
                if ctr_2_exists

                    % iterate over all chord stations
                    for j=1:length(chord_spacing_ctr_2)
                        
                        % starting from spoiler TE, calculates relative
                        % chord position of chord-wise points. Applies
                        % chord_spacing_ctr_2 spacing to rel_end_sp 
                        relative_chord_pos = xsi_sp_TE_ctr2_arr(i) + chord_spacing_ctr_2(j)*(rel_end_sp(i));
                        
                        % Stores previous chord position and accounts for taper
                        % of TE CS when looking at chord position on the main
                        % wing
                        if j>1
                            relative_chord_pos_prv=xsi_sp_TE_ctr2_arr(i) + chord_spacing_ctr_2(j-1)*(rel_end_sp(i));
                            skeleton_point_prv=skeleton_point2;
                        end

                        skeleton_point=obj.compute_skeleton_point(relative_chord_pos,span_spacing(i));

                        skeleton_point2=obj.compute_skeleton_point_new(relative_chord_pos,relative_chord_pos_prv,skeleton_point_prv, span_spacing(i));

                        skeleton_point=skeleton_point2;
                        
                        % the r1_ext and r2_ext are calculated within the
                        % spoiler section. This does not cause issues as
                        % ctr2 can only exist if the spoilers exist as
                        % well. Moreover, the current span station is
                        % always the same for both.
                        grid(:,k)=r1_ext + relative_chord_pos * (r2_ext - r1_ext);
                        
                        grid_flat(:,k)=grid(:,k);
                        if grid3D==1
                           [lower_point,upper_point]=obj.compute_thickness_point(relative_chord_pos,span_spacing(i));
                           grid_upper(:,k)=grid(:,k);
                           grid_lower(:,k)=grid(:,k);
                           grid_upper(2,k)=grid(2,k)-upper_point*sind(obj.dihed);
                           grid_upper(3,k)=grid(3,k)+upper_point*cosd(obj.dihed);
                           grid_lower(2,k)=grid(2,k)-lower_point*sind(obj.dihed);
                           grid_lower(3,k)=grid(3,k)+lower_point*cosd(obj.dihed);
                        end
                        grid(2,k)=grid(2,k)-skeleton_point*sind(obj.dihed);
                        grid(3,k)=grid(3,k)+skeleton_point*cosd(obj.dihed);

                        % increases grid point counter by 1
                        k=k+1;
                    end
                end


                % runs if TE CS is present on wing segment
                if(obj.has_te_cs)
                    
                    % loops over all TE devices
                    for j=1:n_te_devices
                        
                        % defines chord spacing (0 to 1, in CS TE reference
                        % frame) of all chord stations
                        chord_spacing_te=0:1/obj.n_te_panels(j):1;
                        chord_spacing_te_overlap = 0:1/n_sp_overlap_panels(j):1;
                        chord_spacing_te_free = 0:1/n_te_free_panels(j):1;
                        % the first element of chord_spacing_te_free is
                        % excluded since [1] of chord_spacing_te_overlap
                        % and [0] of chord_spacing_te_free represent the
                        % same point.
                        chord_spacing_te_free = chord_spacing_te_free(2:end);
                        
                        % same as r1 and r2 coords previous illustrated,
                        % only for CS TE.
                        r1=obj.xyz_te_device(j,:,1)+span_spacing(i)*(obj.xyz_te_device(j,:,2)-obj.xyz_te_device(j,:,1));
                        r2=obj.xyz_te_device(j,:,4)+span_spacing(i)*(obj.xyz_te_device(j,:,3)-obj.xyz_te_device(j,:,4));
                        
                        if overlap_flag == 1
                            r1_overlap = r1;
                            
                            rel_chord_overlap =  (xsi_sp_TE_ctr2_arr(i)-(1-obj.te_device.hinge*te_taper_correction(i)))/(obj.te_device.hinge*te_taper_correction(i));
                            r2_overlap = r1 + rel_chord_overlap*(r2-r1);
%                             r2_overlap = obj.xyz_sp_device(1,:,4)+span_spacing(i)*(obj.xyz_sp_device(1,:,3)-obj.xyz_sp_device(1,:,4));
                            
                            r1_free = r2_overlap;
                            r2_free = r2;
                        end
                        
                        % loops over all chord stations within current span
                        % station (j)
                        for ii=1:length(chord_spacing_te)
                            
                            % saves prev chord pos and prev skeleton point
                            if ii>1
                                relative_chord_pos_prv=relative_chord_pos;
                                skeleton_point_prv=skeleton_point2;
                            end
                                   
                            if overlap_flag == 0
                                relative_chord_pos = (1-obj.te_device.hinge*te_taper_correction(i)) + chord_spacing_te(ii) * obj.te_device.hinge * te_taper_correction(i);
                                grid(:,k)=r1+chord_spacing_te(ii)*(r2-r1);
                            else
                                % te_surface_LE is the flap leading edge xsi coordinate at the current span
                                % station (code snippet taken from TE surface if loop below)
                                te_surface_LE = (1-obj.te_device.hinge*te_taper_correction(i));

                                % checking whether we are looking at the
                                % overlapping region or the free one
                                if ii<= length(chord_spacing_te_overlap)
                                    relative_chord_pos = te_surface_LE + (xsi_sp_TE_ctr2_arr(i)-te_surface_LE)*chord_spacing_te_overlap(ii);
                                    grid(:,k)=r1_overlap+chord_spacing_te_overlap(ii)*(r2_overlap-r1_overlap);
                                else
                                    ii_loc = ii - length(chord_spacing_te_overlap);
                                    relative_chord_pos = xsi_sp_TE_ctr2_arr(i) + (1-xsi_sp_TE_ctr2_arr(i)) * chord_spacing_te_free(ii_loc);
                                    grid(:,k)=r1_free+chord_spacing_te_free(ii_loc)*(r2_free-r1_free);
                                end
                            end
                            
                            
                            
                            
                            skeleton_point=obj.compute_skeleton_point(relative_chord_pos,span_spacing(i));  
                            skeleton_point2=obj.compute_skeleton_point_new(relative_chord_pos,relative_chord_pos_prv,skeleton_point_prv,span_spacing(i));  

                            skeleton_point=skeleton_point2;
                            
                            % calculates xyz of current grid point (current
                            % chord station and current span station on
                            % current TE device)

                            
                            grid_flat(:,k)=grid(:,k);
                            if grid3D==1
                                [lower_point,upper_point]=obj.compute_thickness_point(relative_chord_pos,span_spacing(i));
                                grid_upper(:,k)=grid(:,k);
                                grid_lower(:,k)=grid(:,k);
                                grid_upper(2,k)=grid(2,k)-upper_point*sind(obj.dihed);
                                grid_upper(3,k)=grid(3,k)+upper_point*cosd(obj.dihed);
                                grid_lower(2,k)=grid(2,k)-lower_point*sind(obj.dihed);
                                grid_lower(3,k)=grid(3,k)+lower_point*cosd(obj.dihed);
                            end
                            grid(2,k)=grid(2,k)-skeleton_point*sind(obj.dihed);
                            grid(3,k)=grid(3,k)+skeleton_point*cosd(obj.dihed);
                            
                            % increases grid point counter by 1
                            k=k+1;
                        end
                    end
                end
                
                te_idx(k-n_chordwise_points:k-1)=k-1;
                

                
                % if wake grid is desired save wake front line
                if wake==1
                    grid_wake(:,i)=grid(:,k-1);
                elseif wake==2
                    grid_wake(:,i)=grid(:,k-1);
                    %grid_wake(:,i)=grid(:,k-1)+0.20*231*2.7056e-04*(grid(:,k-1)-grid(:,k-2))/norm(grid(:,k-1)-grid(:,k-2));
                    %grid_wake(:,i)=grid(:,k-1)+0.25*(grid(:,k-1)-grid(:,k-2));
                end
            end
            
            x=[];
            
            % compute panel indices
            
            % chordwise_pan will change depending on whether ctr2 exists or not
            if ctr_2_exists
                chordwise_pan=[obj.n_le_panels n_chord_panels_ctr_1 obj.n_sp_panels n_chord_panels_ctr_2 obj.n_te_panels];
            else
                chordwise_pan=[obj.n_le_panels n_chord_panels_ctr_1 obj.n_sp_panels obj.n_te_panels];
            end
            end_points=cumsum(chordwise_pan+1);
            offset=0;
            for ii=1:obj.n_span
                for jj=1:(length(end_points))
                    if jj==1
                        x=[x 1+offset:(end_points(jj)-1)+offset];
                    else
                        x=[x (end_points(jj-1)+1)+offset:(end_points(jj)-1)+offset];
                    end
                end
                offset=offset+sum(chordwise_pan+1);
            end
            panels=zeros(4,obj.n_span*n_chordwise_panels);
           
            % for potential flow solution
            obj.is_te=zeros(1,size(panels,2));
            
            for i=1:obj.n_span*n_chordwise_panels
                panels(:,i)=[x(i);x(i)+sum(chordwise_pan+1);x(i)+sum(chordwise_pan+1)+1;x(i)+1];
                if mod(i,n_chordwise_panels)==0
                    obj.is_te(i)=1;
                end
            end
            
            
            % if wake grid is desired
            if (wake==1)||(wake==2)
               len=length(span_spacing);
               for i=1:len-1
                   panels_wake(1,i)=i;
                   panels_wake(2,i)=i+1;
                   panels_wake(3,i)=i+1+len;
                   panels_wake(4,i)=i+len;
               end 
               obj.grid_wake=[grid_wake grid_wake];
               obj.panels_wake=panels_wake;
            end
            
            obj.grid=grid;
            obj.grid_vol_upper=grid_upper;
            obj.grid_vol_lower=grid_lower;
            obj.grid_flat=grid_flat;
            obj.panels=panels;
            obj.te_idx=te_idx;
        end
        
        function skeleton_point=compute_skeleton_point(obj,relative_chord_pos,relative_span_pos)
            nprof_upper_t=obj.profile_t(1,1);
            nprof_lower_t=obj.profile_t(1,2);
            coords_upper_t=obj.profile_t(2:1+nprof_upper_t,1);
            profile_upper_t=obj.profile_t(2:1+nprof_upper_t,2);
            coords_lower_t=obj.profile_t(2+nprof_upper_t:1+nprof_upper_t+nprof_lower_t,1);
            profile_lower_t=obj.profile_t(2+nprof_upper_t:1+nprof_upper_t+nprof_lower_t,2);
           % skeleton_line_t=obj.c_t*0.5*(interp1(coords_lower_t,profile_lower_t,relative_chord_pos,'lin','extrap')+interp1(coords_upper_t,profile_upper_t,relative_chord_pos,'lin','extrap'));
            
            skeleton_line_t=obj.c_t*0.5*(nakeinterp1(coords_lower_t,profile_lower_t,relative_chord_pos)+nakeinterp1(coords_upper_t,profile_upper_t,relative_chord_pos));
            nprof_upper_r=obj.profile_r(1,1);
            nprof_lower_r=obj.profile_r(1,2);
            coords_upper_r=obj.profile_r(2:1+nprof_upper_r,1);
            profile_upper_r=obj.profile_r(2:1+nprof_upper_r,2);
            coords_lower_r=obj.profile_r(2+nprof_upper_r:1+nprof_upper_r+nprof_lower_r,1);
            profile_lower_r=obj.profile_r(2+nprof_upper_r:1+nprof_upper_r+nprof_lower_r,2);
           % skeleton_line_r=obj.c_r*0.5*(interp1(coords_lower_r,profile_lower_r,relative_chord_pos,'lin','extrap')+interp1(coords_upper_r,profile_upper_r,relative_chord_pos,'lin','extrap'));
            skeleton_line_r=obj.c_r*0.5*(nakeinterp1(coords_lower_r,profile_lower_r,relative_chord_pos)+nakeinterp1(coords_upper_r,profile_upper_r,relative_chord_pos));
            skeleton_point=skeleton_line_r*(1-relative_span_pos)+skeleton_line_t*relative_span_pos;
        end
        function skeletonPoint=compute_skeleton_point_new(obj,relative_chord_pos,relative_chord_pos_prv,skeletonPointPrv,relative_span_pos)
            nprof_upper_t=obj.profile_t(1,1);
            nprof_lower_t=obj.profile_t(1,2);
            coords_upper_t=obj.profile_t(2:1+nprof_upper_t,1);
            profile_upper_t=obj.profile_t(2:1+nprof_upper_t,2);
            coords_lower_t=obj.profile_t(2+nprof_upper_t:1+nprof_upper_t+nprof_lower_t,1);
            profile_lower_t=obj.profile_t(2+nprof_upper_t:1+nprof_upper_t+nprof_lower_t,2);
            if ~isequal(coords_upper_t,coords_lower_t)
                lower_xOld=coords_lower_t;
                lower_yOld=profile_lower_t;
                profile_lower_t=[];
                coords_lower_t=coords_upper_t;
                for iPoint=1:length(coords_upper_t)
                    profile_lower_t(iPoint,1)=nakeinterp1(lower_xOld,lower_yOld,coords_lower_t(iPoint));
                end
            end
                camberline_t=(profile_upper_t+profile_lower_t)/2;
                coords_camberline_t=coords_lower_t;
                slope_t=diff(camberline_t)./diff(coords_camberline_t);
                coords_slope_t=(coords_camberline_t(1:end-1)+coords_camberline_t(2:end))/2;
            slope_pos=3/4*(relative_chord_pos-relative_chord_pos_prv)+relative_chord_pos_prv;
            slope_p_t=nakeinterp1(coords_slope_t,slope_t,slope_pos);
            skeletonPointT=slope_p_t*(relative_chord_pos-relative_chord_pos_prv);
%            skeleton_line_t=obj.c_t*0.5*(interp1(coords_lower_t,profile_lower_t,relative_chord_pos,'lin','extrap')+interp1(coords_upper_t,profile_upper_t,relative_chord_pos,'lin','extrap'));
            
            nprof_upper_r=obj.profile_r(1,1);
            nprof_lower_r=obj.profile_r(1,2);
            coords_upper_r=obj.profile_r(2:1+nprof_upper_r,1);
            profile_upper_r=obj.profile_r(2:1+nprof_upper_r,2);
            coords_lower_r=obj.profile_r(2+nprof_upper_r:1+nprof_upper_r+nprof_lower_r,1);
            profile_lower_r=obj.profile_r(2+nprof_upper_r:1+nprof_upper_r+nprof_lower_r,2);
            if ~isequal(coords_upper_r,coords_lower_r)
                lower_xOld=coords_lower_r;
                lower_yOld=profile_lower_r;
                coords_lower_r=coords_upper_r;
                profile_lower_r=[];
                for iPoint=1:length(coords_upper_r)
                    profile_lower_r(iPoint,1)=nakeinterp1(lower_xOld,lower_yOld,coords_lower_r(iPoint));
                end
            end
                camberline_r=(profile_upper_r+profile_lower_r)/2;
                coords_camberline_r=coords_lower_r;
                slope_r=diff(camberline_r)./diff(coords_camberline_r);
                coords_slope_r=(coords_camberline_r(1:end-1)+coords_camberline_r(2:end))/2;
            slope_p_r=nakeinterp1(coords_slope_r,slope_r,slope_pos);
            skeletonPointR=slope_p_r*(relative_chord_pos-relative_chord_pos_prv);
%             skeleton_line_r=obj.c_r*0.5*(interp1(coords_lower_r,profile_lower_r,relative_chord_pos,'lin','extrap')+interp1(coords_upper_r,profile_upper_r,relative_chord_pos,'lin','extrap'));
%             skeletonPoint=skeleton_line_r*(1-relative_span_pos)+skeleton_line_t*relative_span_pos;
             skeletonPoint=(obj.c_r*skeletonPointR*(1-relative_span_pos)+obj.c_t*skeletonPointT*relative_span_pos)+skeletonPointPrv;
        end
        function skeletonAngle=compute_skeleton_angle(obj,relative_chord_pos,relative_chord_pos_0,relative_chord_pos_1,relative_span_pos)
            %relative_chord_pos_0: leading edge of the panel
            %relative_chord_pos_1: trailing edge of the panel
            
            % function to determine which panel of the discretized
            % camberline is crossed by the perpendicular segment of the chosen
            % panel (depending on approach) - detailed report has been
            % writen
            
            %tip initialization
            nPointsUpperTip=obj.profile_t(1,1);
            nPointsLowerTip=obj.profile_t(1,2);
            coordsUpperTip=obj.profile_t(2:1+nPointsUpperTip,1);
            profileUpperTip=obj.profile_t(2:1+nPointsUpperTip,2);
            coordsLowerTip=obj.profile_t(2+nPointsUpperTip:1+nPointsUpperTip+nPointsLowerTip,1);
            profileLowerTip=obj.profile_t(2+nPointsUpperTip:1+nPointsUpperTip+nPointsLowerTip,2);
            %root initialization
            nPointsUpperRoot=obj.profile_r(1,1);
            nPointsLowerRoot=obj.profile_r(1,2);
            coordsUpperRoot=obj.profile_r(2:1+nPointsUpperRoot,1);
            profileUpperRoot=obj.profile_r(2:1+nPointsUpperRoot,2);
            coordsLowerRoot=obj.profile_r(2+nPointsUpperRoot:1+nPointsUpperRoot+nPointsLowerRoot,1);
            profileLowerRoot=obj.profile_r(2+nPointsUpperRoot:1+nPointsUpperRoot+nPointsLowerRoot,2);
            %settings
            nSteps=100;
            coordsInterp=0:1/nSteps:1;
            %tip interpolation and skeleton line
            coordsUpperTipNorm=coordsUpperTip/(coordsUpperTip(end) - coordsUpperTip(1));
            profileUpperTipNorm=profileUpperTip/(coordsUpperTip(end) - coordsUpperTip(1));
            profileUpperTipInterp=interp1(coordsUpperTipNorm,profileUpperTipNorm,coordsInterp);
            
            coordsLowerTipNorm=coordsLowerTip/(coordsLowerTip(end) - coordsLowerTip(1));
            profileLowerTipNorm=profileLowerTip/(coordsLowerTip(end) - coordsLowerTip(1));
            profileLowerTipInterp=interp1(coordsLowerTipNorm,profileLowerTipNorm,coordsInterp);
                        
            skeletonLineTipInterp=0.5*(profileUpperTipInterp+profileLowerTipInterp);
            % skeletonAngleTip=atan(diff(skeletonLineTipInterp)/(1/nSteps));
            
            % Position of leading edge, trailing edge and 3/4 chord of the
            % panel
            y_relative_chord_pos_0_Tip=interp1(coordsInterp,skeletonLineTipInterp,relative_chord_pos_0);
            y_relative_chord_pos_1_Tip=interp1(coordsInterp,skeletonLineTipInterp,relative_chord_pos_1);
            y_relative_chord_pos=interp1([relative_chord_pos_0 relative_chord_pos_1],[y_relative_chord_pos_0_Tip y_relative_chord_pos_1_Tip],relative_chord_pos);
            Point_ref=[relative_chord_pos,y_relative_chord_pos];
            CamberLineTip=[coordsInterp',skeletonLineTipInterp']; % Coordinates of the camberline
            
            %Slope of the perpendicular segment (m_ref) on 3/4 chord of the panel
            if relative_chord_pos_0==relative_chord_pos_1
                m_ref=0;
            elseif y_relative_chord_pos_0_Tip==y_relative_chord_pos_1_Tip
                m_ref=1E10;
            else
                m=(y_relative_chord_pos_0_Tip - y_relative_chord_pos_1_Tip)/(relative_chord_pos_0 - relative_chord_pos_1);
                m_ref=(-1)/m;
            end
            
            % Find the crossed segment
            for iPointCamberLine=1:(size(CamberLineTip,1)-1)
                P_0=CamberLineTip(iPointCamberLine,:); % Point of the camberline
                P_1=CamberLineTip(iPointCamberLine+1,:); % Following point of the camberline
                
                yaux=y_relative_chord_pos+20; 
                
                Point_aux=[Point_ref(1)-((Point_ref(2)-yaux)/m_ref) yaux]; % Given point which has m_ref slope, and y coordinate higher than the highest y of the segment
                segment_1=Point_aux-Point_ref; % Segment perpendicular to the chosen panel
                segment_2=P_1-Point_ref; % Segment of the discretized camberline (from 3/4 chord to given point of the discretized camberline)
                segment_3=P_0-Point_ref; % Segment of the discretized camberline (from 3/4 chord to  a next given point of the discretized camberline)
                
                VecProd1=(segment_1(1,1)*segment_2(1,2))-(segment_1(1,2)*segment_2(1,1));
                VecProd2=(segment_1(1,1)*segment_3(1,2))-(segment_1(1,2)*segment_3(1,1));
                proof=VecProd1*VecProd2; % To find the panel crossed by the perpendicular segment, the third components of the cross product must have opposite directions (negative result) 
                
                if proof<=0
                    %                     if P_0(1,1)==P_1(1,1)
                    %                         relativeChordSkeletonAngleTip=0;
                    %                     elseif P_0(1,2)==P_1(1,2)
                    %                         relativeChordSkeletonAngleTip=pi/2;
                    %                     else
                    %                     m_Panel=(P_1(2)-P_0(2)) / (P_1(1)-P_0(1));
                    %                     relativeChordSkeletonAngleTip=atan((-1) / m_Panel);
                    %                     end
                    
                    relativeChordSkeletonAngleTip=(P_1(2)-P_0(2)) / (P_1(1)-P_0(1)); % Slope of the segment
                end
            end
            
            % The same implementation has been done to the root
            % camberline
            
            %root interpolation and skeleton line
            coordsUpperRootNorm=coordsUpperRoot/(coordsUpperRoot(end) - coordsUpperRoot(1));
            profileUpperRootNorm=profileUpperRoot/(coordsUpperRoot(end) - coordsUpperRoot(1));
            profileUpperRootInterp=interp1(coordsUpperRootNorm,profileUpperRootNorm,coordsInterp);
            
            coordsLowerRootNorm=coordsLowerRoot/(coordsLowerRoot(end) - coordsLowerRoot(1));
            profileLowerRootNorm=profileLowerRoot/(coordsLowerRoot(end) - coordsLowerRoot(1));
            profileLowerRootInterp=interp1(coordsLowerRootNorm,profileLowerRootNorm,coordsInterp);
                        
            skeletonLineRootInterp=0.5*(profileUpperRootInterp+profileLowerRootInterp);
            %             skeletonAngleRoot=atan(diff(skeletonLineRootInterp)/(1/nSteps));
            
            y_relative_chord_pos_0_Root=interp1(coordsInterp,skeletonLineRootInterp,relative_chord_pos_0);
            y_relative_chord_pos_1_Root=interp1(coordsInterp,skeletonLineRootInterp,relative_chord_pos_1);
            y_relative_chord_pos=interp1([relative_chord_pos_0 relative_chord_pos_1],[y_relative_chord_pos_0_Root y_relative_chord_pos_1_Root],relative_chord_pos);
            Point_ref=[relative_chord_pos,y_relative_chord_pos];
            CamberLineRoot=[coordsInterp',skeletonLineRootInterp'];
            
            if relative_chord_pos_0==relative_chord_pos_1
                m_ref=0;
            elseif y_relative_chord_pos_0_Root==y_relative_chord_pos_1_Root
                m_ref=1E10;
            else
                m=(y_relative_chord_pos_0_Root - y_relative_chord_pos_1_Root)/(relative_chord_pos_0 - relative_chord_pos_1);
                m_ref=(-1)/m;
            end
            
            for iPointCamberLine=1:(size(CamberLineRoot,1)-1)
                P_0=CamberLineRoot(iPointCamberLine,:);
                P_1=CamberLineRoot(iPointCamberLine+1,:);
                
                yaux=y_relative_chord_pos+20;
                
                Point_aux=[Point_ref(1)-((Point_ref(2)-yaux)/m_ref) yaux];
                segment_1=Point_aux-Point_ref;
                segment_2=P_1-Point_ref;
                segment_3=P_0-Point_ref;
                
                VecProd1=(segment_1(1,1)*segment_2(1,2))-(segment_1(1,2)*segment_2(1,1));
                VecProd2=(segment_1(1,1)*segment_3(1,2))-(segment_1(1,2)*segment_3(1,1));
                proof=VecProd1*VecProd2;
                
                if proof<=0
                    %                     if P_0(1)==P_1(1)
                    %                         relativeChordSkeletonAngleRoot=0;
                    %                     elseif P_0(2)==P_1(2)
                    %                         relativeChordSkeletonAngleRoot=pi/2;
                    %                     else
                    %                     m_Panel=(P_1(2)-P_0(2)) / (P_1(1)-P_0(1));
                    %                     relativeChordSkeletonAngleRoot=atan((-1) / m_Panel);
                    %                     end
                    relativeChordSkeletonAngleRoot=(P_1(2)-P_0(2)) / (P_1(1)-P_0(1)); % Slope of the segment
                end
            end
            
            %             relativeChordSkeletonAngleTip=interp1(coordsInterp(2:end),skeletonAngleTip,relative_chord_pos);
            %             relativeChordSkeletonAngleRoot=interp1(coordsInterp(2:end),skeletonAngleRoot,relative_chord_pos);
            
            skeletonAngle=relative_span_pos*relativeChordSkeletonAngleTip+(1-relative_span_pos)*relativeChordSkeletonAngleRoot;
            % plot profile:
%             figure
%             hold on
%             grid on
%             plot(coordsInterp,skeletonLineRootInterp)
%             plot(coordsLowerRootNorm,profileLowerRootNorm)
%             plot(coordsUpperRootNorm,profileUpperRootNorm)
        end
        
        function [lower_point,upper_point]=compute_thickness_point(obj,relative_chord_pos,relative_span_pos)
            nprof_upper_t=obj.profile_t(1,1);
            nprof_lower_t=obj.profile_t(1,2);
            coords_upper_t=obj.profile_t(2:1+nprof_upper_t,1);
            profile_upper_t=obj.profile_t(2:1+nprof_upper_t,2);
            coords_lower_t=obj.profile_t(2+nprof_upper_t:1+nprof_upper_t+nprof_lower_t,1);
            profile_lower_t=obj.profile_t(2+nprof_upper_t:1+nprof_upper_t+nprof_lower_t,2);
            
            lower_line_t=obj.c_t*nakeinterp1(coords_lower_t,profile_lower_t,relative_chord_pos);
            %lower_line_t=obj.c_t*interp1(coords_lower_t,profile_lower_t,relative_chord_pos,'lin','extrap');
            
            %upper_line_t=obj.c_t*interp1(coords_upper_t,profile_upper_t,relative_chord_pos,'lin','extrap');
            upper_line_t=obj.c_t*nakeinterp1(coords_upper_t,profile_upper_t,relative_chord_pos);
            
            nprof_upper_r=obj.profile_r(1,1);
            nprof_lower_r=obj.profile_r(1,2);
            coords_upper_r=obj.profile_r(2:1+nprof_upper_r,1);
            profile_upper_r=obj.profile_r(2:1+nprof_upper_r,2);
            coords_lower_r=obj.profile_r(2+nprof_upper_r:1+nprof_upper_r+nprof_lower_r,1);
            profile_lower_r=obj.profile_r(2+nprof_upper_r:1+nprof_upper_r+nprof_lower_r,2);
            %lower_line_r=obj.c_r*interp1(coords_lower_r,profile_lower_r,relative_chord_pos,'lin','extrap');
            lower_line_r=obj.c_r*nakeinterp1(coords_lower_r,profile_lower_r,relative_chord_pos);
            
            %upper_line_r=obj.c_r*interp1(coords_upper_r,profile_upper_r,relative_chord_pos,'lin','extrap');
            upper_line_r=obj.c_r*nakeinterp1(coords_upper_r,profile_upper_r,relative_chord_pos);
            
            lower_point=lower_line_r*(1-relative_span_pos)+lower_line_t*relative_span_pos;
            upper_point=upper_line_r*(1-relative_span_pos)+upper_line_t*relative_span_pos; 
        end
        
        function profile=get_profile(obj,relative_span_pos)
            nprof_upper_t=obj.profile_t(1,1);
            nprof_lower_t=obj.profile_t(1,2);
            coords_upper_t=obj.profile_t(2:1+nprof_upper_t,1);
            profile_upper_t=obj.profile_t(2:1+nprof_upper_t,2);
            coords_lower_t=obj.profile_t(2+nprof_upper_t:1+nprof_upper_t+nprof_lower_t,1);
            profile_lower_t=obj.profile_t(2+nprof_upper_t:1+nprof_upper_t+nprof_lower_t,2);
            
            nprof_upper_r=obj.profile_r(1,1);
            nprof_lower_r=obj.profile_r(1,2);
            coords_upper_r=obj.profile_r(2:1+nprof_upper_r,1);
            profile_upper_r=obj.profile_r(2:1+nprof_upper_r,2);
            coords_lower_r=obj.profile_r(2+nprof_upper_r:1+nprof_upper_r+nprof_lower_r,1);
            profile_lower_r=obj.profile_r(2+nprof_upper_r:1+nprof_upper_r+nprof_lower_r,2);
            %% TODO: profiles don#t always have same size
            profile_upper=profile_upper_r*(1-relative_span_pos)+profile_upper_t*relative_span_pos;
            profile_lower=profile_lower_r*(1-relative_span_pos)+profile_lower_t*relative_span_pos;
            
            coords_upper=coords_upper_r*(1-relative_span_pos)+coords_upper_t*relative_span_pos;
            coords_lower=coords_lower_r*(1-relative_span_pos)+coords_lower_t*relative_span_pos;
            
            profile=[[length(profile_upper) coords_upper' coords_lower']; [length(profile_lower) profile_upper' profile_lower']]';
        end
        
        function obj=write_tecplot_grid(obj,fileID,name,id)
            nxt=0;
            fprintf(fileID,'zone t="%sSEGMENT%i_UPPER", i=%i, j=1, k=1, f=point \n',name,id,obj.te_idx(1));
            if id==1
                i=1;
            else
                i=obj.te_idx(1)+1;
            end
            while i<=size(obj.grid_vol_upper,2)
                if nxt==1
                  i=obj.te_idx(end)-obj.te_idx(1)+1;  
                  fprintf(fileID,'zone t="%sSEGMENT%i_UPPER", i=%i, j=1, k=1, f=point \n',name,id+i,obj.te_idx(i)-i+1); 
                  nxt=0;
                end
                if (i==obj.te_idx(i))
                  nxt=1;
                end
              fprintf(fileID,'%f %f %f\n',obj.grid_vol_upper(1,i),obj.grid_vol_upper(2,i),obj.grid_vol_upper(3,i));         
              i=i+1;
            end
            nxt=0;
            fprintf(fileID,'zone t="%sSEGMENT%i_LOWER", i=%i, j=1, k=1, f=point \n',name,id,obj.te_idx(1)); 
            if id==1
                i=1;
            else
                i=obj.te_idx(1)+1;
            end
            while i<=length(obj.grid_vol_lower)
              if nxt==1
                i=obj.te_idx(end)-obj.te_idx(1)+1; 
                fprintf(fileID,'zone t="%sSEGMENT%i_LOWER", i=%i, j=1, k=1, f=point \n',name,id+i,obj.te_idx(i)-i+1); 
                nxt=0;
              end   
              if (i==obj.te_idx(i))
                  nxt=1;
              end
              fprintf(fileID,'%f %f %f\n',obj.grid_vol_lower(1,i),obj.grid_vol_lower(2,i),obj.grid_vol_lower(3,i));         
              i=i+1;
            end
        end
        
        function obj=read_xml_definition(obj,xmlstruct)
            if strcmp(xmlstruct.tag,'SEGMENT')
                obj.pos(1)=str2double(xmlstruct.child(1).child(1).value);
                obj.pos(2)=str2double(xmlstruct.child(1).child(2).value);
                obj.pos(3)=str2double(xmlstruct.child(1).child(3).value);
                obj.dihed=str2double(xmlstruct.child(2).value);
                if strcmp(xmlstruct.child(3).tag,'REAL_SWEEP')
                    obj.le_sweep=[];
                    obj.c4_sweep=[];
                    obj.real_sweep=str2double(xmlstruct.child(3).value);
                end
                obj.le_sweep=str2double(xmlstruct.child(3).value);
                obj.b=str2double(xmlstruct.child(4).value);
                obj.c_r=str2double(xmlstruct.child(5).value);
                obj.c_t=str2double(xmlstruct.child(6).value);
                obj.Theta_r=str2double(xmlstruct.child(7).value);
                obj.Theta_t=str2double(xmlstruct.child(8).value);
                obj.profile_name_r=xmlstruct.child(9).value;
                obj.profile_name_t=xmlstruct.child(10).value;
                
                obj.has_te_cs=0;
                obj.has_le_cs=0;
                if ~isempty(obj.le_sweep)
                obj.c4_sweep=atan((obj.b*tan(obj.le_sweep*pi/180)+obj.c_t/4-obj.c_r/4)/obj.b)*180/pi;
                end
                obj=obj.complete_params_from_stdinit();
                obj=obj.load_airfoils();
                obj=obj.compute_segment_coordinates();
                
            else
                fprintf('Unknown Data Format: %s \n', xmlstruct.tag);
            end
        end
        
        % Finds and assigns repsective panel Ids to their corresponding
        % control surface (within wing-segment obj).
        function obj=  computeControlSurfacePanelIds(obj)
            
            % Total number of chordwise panels, from LE to TE (including
            % all surfaces and regions)
            n_panels = sum([obj.n_le_panels, obj.n_ctr1_panels, obj.n_sp_panels, obj.n_ctr2_panels, obj.n_te_panels]);
            
            % sum is used to turn empty elements ([]) into zeros.
            n_le_panels = sum(obj.n_le_panels);
            n_ctr1_panels = sum(obj.n_ctr1_panels);
            n_sp_panels = sum(obj.n_sp_panels);
            n_ctr2_panels = sum(obj.n_ctr2_panels);
            n_te_panels = sum(obj.n_te_panels);
            
            
            % Looks at panels belonging to TE
            if obj.has_te_cs      
                
                % In case of NO overlap with spoilers
                if isempty(obj.n_te_panels_overlap)
                    obj.te_device.panelIds=zeros(1,obj.n_span*obj.n_te_panels);
                    for iSpan=1:obj.n_span
                        obj.te_device.panelIds((iSpan-1)*obj.n_te_panels+1:(iSpan)*obj.n_te_panels)= (obj.panel_start_idx+iSpan*(n_panels-obj.n_te_panels)+(iSpan-1)*obj.n_te_panels:obj.panel_start_idx-1+iSpan*(n_panels-obj.n_te_panels)+(iSpan)*obj.n_te_panels);
                    end
                    
                % In case of overlap with spoilers
                else
                    % Panels belonging to free flap region are stored in panelIds
                    obj.te_device.panelIds=zeros(1,obj.n_span*obj.n_te_panels_free);
                    
                    % Panels belonging to overlapped flap region are stored in panelIds_special
                    obj.te_device.panelIds_special=zeros(1,obj.n_span*obj.n_te_panels_overlap);
                    
                    for iSpan=1:obj.n_span
                        obj.te_device.panelIds_special((iSpan-1)*obj.n_te_panels_overlap+1:(iSpan)*obj.n_te_panels_overlap)= ...
                            (obj.panel_start_idx + n_le_panels + n_ctr1_panels + n_sp_panels+(iSpan-1)*n_panels:obj.panel_start_idx-1 + n_le_panels + n_ctr1_panels + n_sp_panels + obj.n_te_panels_overlap + (iSpan-1)*n_panels);
                        obj.te_device.panelIds((iSpan-1)*obj.n_te_panels_free+1:(iSpan)*obj.n_te_panels_free)= ...
                            (obj.panel_start_idx + n_le_panels + n_ctr1_panels + n_sp_panels + obj.n_te_panels_overlap + (iSpan-1)*n_panels:obj.panel_start_idx-1 + n_le_panels + n_ctr1_panels + n_sp_panels + obj.n_te_panels_overlap + obj.n_te_panels_free + (iSpan-1)*n_panels);
                    end
                end
            end
            
            % Looks at panels belonging to LE
            if obj.has_le_cs
                obj.le_device.panelIds=zeros(1,obj.n_span*n_le_panels);
                for iSpan=1:obj.n_span
                    obj.le_device.panelIds((iSpan-1)*n_le_panels+1:(iSpan)*n_le_panels)= (obj.panel_start_idx+(iSpan-1)*n_panels:obj.panel_start_idx -1 + n_le_panels +(iSpan-1)*n_panels);
                end
            end
            
            % Looks at panels belonging to spoilers
            if obj.has_sp_cs
                obj.sp_device.panelIds=zeros(1,obj.n_span*n_sp_panels);
                for iSpan=1:obj.n_span
                    obj.sp_device.panelIds((iSpan-1)*n_sp_panels+1:(iSpan)*n_sp_panels)= (obj.panel_start_idx + n_le_panels + n_ctr1_panels+(iSpan-1)*n_panels:obj.panel_start_idx-1 + n_le_panels + n_ctr1_panels + n_sp_panels + (iSpan-1)*n_panels);
                end
            end
            
            % Code below makes sure that panelIds are arranged according to
            % our convention. Convention specified that:
                 % - PanelIds array shall contain all panels belonging to
                 %   that surface. This means that in case of spoiler/flap
                 %   overlap, the panel Ids arrays of both spoilers and
                 %   flaps will contain panels in both the free region and
                 %   overlap region
                 %
                 % - PanelIds_special array shall contain all panels
                 %   belonging to special areas (such as overlap). Both CS
                 %   involved will have a PanelIds_special array, which
                 %   will of course contain the same IDs.
                 %
                 % - PanelIds_standard array shall contain all panels
                 %   belonging to standard areas (only used when
                 %   PanelIds_special is not empty). The union of
                 %   PanelIds_special and PanelIds_standard shall yield
                 %   PanelIds.
                 
           if obj.has_te_cs
               % creating panelIds_standard array
               obj.te_device.panelIds_standard = obj.te_device.panelIds;
               if ~isempty(obj.te_device.panelIds_special)
                  
                  % adding special IDs to panelIDs
                  obj.te_device.panelIds = [obj.te_device.panelIds, obj.te_device.panelIds_special];
                  
                  % rearranging panelIds arrays so that they are in
                  % ascending order
                  obj.te_device.panelIds = sort(obj.te_device.panelIds,'ascend');
              end
           end
           
           if obj.has_sp_cs
              % creating panelIds_standard array
              obj.sp_device.panelIds_standard = obj.sp_device.panelIds;
                  
               %the check here is performed within the TE, since at this
               %point of the code special panels are only assigned to the
               %flap (in case of overlap)
               if obj.has_te_cs
                  if ~isempty(obj.te_device.panelIds_special)

                      % creating panelIds_special array (IDs from flap are
                      % used. Overlap panels will be the same)
                      obj.sp_device.panelIds_special = obj.te_device.panelIds_special;

                      % adding special IDs to panelIDs
                      obj.sp_device.panelIds = [obj.sp_device.panelIds, obj.sp_device.panelIds_special];

                      % rearranging panelIds arrays so that they are in
                      % ascending order
                      obj.sp_device.panelIds = sort(obj.sp_device.panelIds,'ascend');
                  end
               end
           end
           
        end
    end
end
