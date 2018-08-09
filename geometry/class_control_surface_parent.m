%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%
classdef class_control_surface_parent < handle
    %CONTROLSURFACE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        %> control surface name
        name;
        %>  1 leading then edge devic, 0 trailing edge device
        pos = 0; 
        %> tapered or constant
        is_tapered = 0;
        %> number of hinge lines
        n_hinge;
        %> chord
        hinge;
        %> deflection
        delta;
        
        %> Aerodynamic Panel Ids
        panelIds=[];
        %> Aerodynamic Panel Ids of Left Wing (so symmetric)
        panelIdsL=[];
        % Aerodynamic Panel Ids of "special regions" (ig overlap between
        % flap and spoiler)
        panelIds_special=[];
        % Same as above, but for Left Wing (symmetric)
        panelIds_special_L=[];
        % Aerodynamic Panel Ids of "NON-special regions", used to denote
        % standard areas of CS. Only present when special panels are also
        % present
        panelIds_standard = [];
        % Same as above, but for Left Wing (symmetric)
        panelIds_standard_L=[];
        
        delta_l_r;
        %> is symmetric
        is_sym;
        %> deflection
        is_sym_defl;
        % CPACS uID
        uID;
        % Parent uID
        parent_uID;
        
        % x-path to the specfici control surface instance
        xpath_cs;
        
        % eta TE, eta LE, xsi LE coordinates of the inner (root) and outer
        % (tip) borders of the control surface
        innerBorder_coords;
        outerBorder_coords;
        
        % xsi coordinate of the hinge line
        hinge_xsi;
        
        % max and min deflections as specfied in the path section (if
        % multiple paths are present, it might be a problem)
        max_defl;
        min_defl;
        
        
        % Info concerning the wing segment where the specific control
        % surface starts
        start_segmentUID; %UID of wing segment
        start_segment_index; %index of wing segment
        
        %Eta and Xsi coords (in wing segment reference frame) of point 1 of
        %CS if the CS is at the TE, or point 4 if the CS is at the LE.
        start_segmentEta; % Eta of CS start in wing segment reference frame
        start_segmentXsi; % Xsi of CS start in wing segment reference frame
        
        
        % Info concerning the wing segment where the specific control
        % surface ends
        end_segmentUID; %UID of wing segment
        end_segment_index; %index of wing segment

        
        %Eta and Xsi coords (in wing segment reference frame) of point 2 of
        %CS if the CS is at the TE, or point 3 if the CS is at the LE.
        end_segmentEta; % Eta of CS end in wing segment reference frame
        end_segmentXsi; % Xsi of CS end in wing segment reference frame

        
        %UID and index of wing containing CS
        wing_UID;
        wing_index;
        
        %tixi and tigl handles
        tixi_handle;
        tigl_handle;
        
        % xyz coord as in CPACS from eta & xsi
        xyz_CPACS;
        
        % attribute specific to xml spoilers, inboard and
        % outboard xsi coords of the spoiler leading edge
        xsi_LE_inner;
        xsi_LE_outer;
        
        % array containing all the child instances of this parent
        children;
        
        % 3D vector in xyz coorindates representing hinge line. VLM normal
        % vector are rotated around this one whenver they are deflected
        rotation_vect;
    end

    methods
        function obj=class_control_surface_parent(name,pos,hinge,varargin)
            if nargin == 0
                %pass
            elseif nargin>0
                    obj.name=name;
                    obj.pos=pos;
                    obj.hinge=str2num(hinge);
                    obj.n_hinge=length(obj.hinge);
                    obj.delta=zeros(1,obj.n_hinge);
                    
                %spoilers always have 2 extra varargin compared to
                %other CS, due to presence of hinge_start. Hence the
                %TE and LE equivalent of nargin==4 runs in nargin==5
                %for spoilers, and so on.
                if nargin==4
                    obj.is_sym=1;
                    obj.is_sym_defl=varargin{1};

                    if obj.is_sym_defl==0
                        obj.delta_l_r=[0 0];
                    end
                end
                
                if nargin==5
                    
                    if pos ~= 2
                        obj.is_sym=1;
                        obj.is_sym_defl=varargin{1};
                        obj.is_tapered=varargin{2};
                        if obj.is_sym_defl==0
                            obj.delta_l_r=[0 0];
                        end
                        
                    else
                        obj.xsi_LE_inner = str2num(varargin{1});
                        obj.xsi_LE_outer = str2num(varargin{2});
                    end
                end
                
                % equivalent of nargin==5 but for spoilers. Pos is not
                % checked since only spoilers can have nargin==6. 
                if nargin==6
                    
                    obj.is_sym=1;
                    obj.is_sym_defl=varargin{1};

                    if obj.is_sym_defl==0
                        obj.delta_l_r=[0 0];
                    end
                    obj.xsi_LE_inner = str2num(varargin{2});
                    obj.xsi_LE_outer = str2num(varargin{3});
                end
                 
               % read comment for nargin==6
                if nargin ==7
                    obj.is_sym=1;
                    obj.is_sym_defl=varargin{1};
                    obj.xsi_LE_inner = str2num(varargin{2});
                    obj.xsi_LE_outer = str2num(varargin{3});
                    obj.is_tapered=varargin{4};
                    if obj.is_sym_defl==0
                        obj.delta_l_r=[0 0];
                    end
                end
            end
        end
    end

    
    methods (Static)
        function obj=create_from_cpacs(tixi, tigl, xpath_cs, pos)
            
             %% ============== Looking at basic info ==============
             
             % creates an empty instance of the control surface class
             obj = class_control_surface_parent();
             
             % saves tixi and tigl handles
                          
             obj.tixi_handle = tixi;
             obj.tigl_handle = tigl;
             
             % xml path to specfic control surface instance
             obj.xpath_cs = xpath_cs;
             
             % saves uID of control surface and parent component segment
             obj.uID = tixiGetTextAttribute(tixi, xpath_cs, 'uID');
             obj.parent_uID = tixiGetTextElement(tixi, [xpath_cs, '/parentUID']);
             
             % saves name of control surface
             obj.name = tixiGetTextElement(tixi, [xpath_cs, '/name']);
             
             % saves whether CS is LE(1) or TE(0)
             obj.pos=pos;
                         
             
             %% ============== Looking at control surface inner and outer border ==============
             
             % String array to store the paths of all coordinates of the
             % control surface. The first row is the inner border paths, the 2nd
             % is the outer.
             x_path_inner = [xpath_cs,'/outerShape','/innerBorder'];
             x_path_outer = [xpath_cs,'/outerShape','/outerBorder'];
             x_path_arr = strings(2,3);
             x_path_arr(1,:) = [[x_path_inner,'/etaLE']; [x_path_inner,'/etaTE']; [x_path_inner,'/xsiLE']];
             x_path_arr(2,:) = [[x_path_outer,'/etaLE']; [x_path_outer,'/etaTE']; [x_path_outer,'/xsiLE']];
             

             % Fetches eta LE of inner and outer border 
             obj.innerBorder_coords(1) = tixiGetDoubleElement(tixi,char(x_path_arr(1,1)));
             obj.outerBorder_coords(1) = tixiGetDoubleElement(tixi,char(x_path_arr(2,1)));
             
             
             % Tries to fetch eta TE. Since it is not always defined, if an error is
             % given, the eta TE coord is set as NaN.
             try
                obj.innerBorder_coords(2) = tixiGetDoubleElement(tixi,char(x_path_arr(1,2)));
                obj.outerBorder_coords(2) = tixiGetDoubleElement(tixi,char(x_path_arr(2,2)));
             catch
                obj.innerBorder_coords(2) = NaN;
                obj.outerBorder_coords(2) = NaN;
             end

             
             
             
             % If the CS is at the TE, the xsi coordinate will be defined
             % as xsiLE, which is already present in the x_path_arr
             if obj.pos==0
                 
                 % Fetches xsi LE
                 obj.innerBorder_coords(3) = tixiGetDoubleElement(tixi,char(x_path_arr(1,3)));
                 obj.outerBorder_coords(3) = tixiGetDoubleElement(tixi,char(x_path_arr(2,3)));

             
             % If the CS is at the LE, the xsi coordinate is defined as
             % either xsiTE or xsiTEUpper & xsiTELower.
             elseif obj.pos == 1
                 
                 % Tries to fetch xsi TE. If undefined, it tries to fetch
                 % xsiTEUpper and xsiTELower
                 try
                     x_path_arr(:,3) = [x_path_inner,'/xsiTE'];
                     obj.innerBorder_coords(3) = tixiGetDoubleElement(tixi,char(x_path_arr(1,3)));
                     obj.outerBorder_coords(3) = tixiGetDoubleElement(tixi,char(x_path_arr(2,3)));
                     
                 catch
                     x_path_arr(:,3) = [x_path_inner,'/xsiTEUpper'];
                     x_path_arr(:,4) = [x_path_inner,'/xsiTELower'];
                     
                     obj.innerBorder_coords(3) = tixiGetDoubleElement(tixi,char(x_path_arr(1,3)));
                     obj.innerBorder_coords(4) = tixiGetDoubleElement(tixi,char(x_path_arr(1,4)));
                     
                     obj.outerBorder_coords(3) = tixiGetDoubleElement(tixi,char(x_path_arr(2,3)));
                     obj.outerBorder_coords(4) = tixiGetDoubleElement(tixi,char(x_path_arr(2,4)));
                 end
                 
             % If the CS is a spoiler, both xsiLE and xsiTE will be present
             elseif obj.pos == 2
                 x_path_arr(:,3) = [x_path_inner,'/xsiLE'];
                 x_path_arr(:,4) = [x_path_inner,'/xsiTE'];

                 obj.innerBorder_coords(3) = tixiGetDoubleElement(tixi,char(x_path_arr(1,3)));
                 obj.innerBorder_coords(4) = tixiGetDoubleElement(tixi,char(x_path_arr(1,4)));

                 obj.outerBorder_coords(3) = tixiGetDoubleElement(tixi,char(x_path_arr(2,3)));
                 obj.outerBorder_coords(4) = tixiGetDoubleElement(tixi,char(x_path_arr(2,4)));
                 
                 obj.xsi_LE_inner = obj.innerBorder_coords(3);
                 obj.xsi_LE_outer = obj.outerBorder_coords(3);
             end
             
             
             % Checks whether the taper ratio of the control surface is
             % exactly the same as the wing. If not, throws a warning
             % message, and then forces the xsi coordinates of the CS to be the same.
             % This is because daedalus cannot handle control surfaces whose
             % distance from the LE has to be described with more than one
             % parameter (property "hinge" in this class)
             if obj.innerBorder_coords(3) ~=  obj.outerBorder_coords(3)
                 
                msg = 'The taper ratio of the CS has to be the same as the wing segment housing it';
                warning(msg);
                
                msg = 'Due to Daedalus allowing for only a single xsi coordinate, the xsi coordinates have been forced to be the same';
                warning(msg);
                
                % Determining taper
                if obj.innerBorder_coords(3) <= obj.outerBorder_coords(3)
                    obj.is_tapered = 1;
                elseif obj.innerBorder_coords(3) > obj.outerBorder_coords(3)
                    obj.is_tapered = 0;
                end
                
                % Forcing outer border xsi to be equal to inner border xsi
                if obj.pos ~= 2
                    obj.outerBorder_coords(3) = obj.innerBorder_coords(3);
                end
                
                % Hinge is basically the chord of the CS. This will be
                % calculated differently depending on whether the CS is at
                % the LE or TE.
                if obj.pos == 1
                    obj.hinge = obj.innerBorder_coords(3);
                    
                elseif obj.pos==0
                    obj.hinge = 1. - obj.innerBorder_coords(3);
                    
                elseif obj.pos==2
                    obj.hinge = obj.innerBorder_coords(4) - obj.innerBorder_coords(3);
                end
                
                % n_hinge states how many hinge lines are present within
                % the same CS. More than one hinge line means that we are
                % dealing with slotted CS, which the current
                % read_from_CPACS cannot handle. Therefore, a warning is
                % thrown, reminding that slotted CS are currently not a
                % feature.
                obj.n_hinge=length(obj.hinge); 
                msg = 'warning: not actually checking amount of hinge lines and presence of slotted flaps. Assuming that no slotted CS are present';
                warning(msg);

             elseif obj.innerBorder_coords(3) ==  obj.outerBorder_coords(3)
                
                % is_tapered attribute set = 1 if xsi coords of inner and
                % outer border are the same, since it means that taper
                % ratio of CS is following exactly taper ratio of wing.
                obj.is_tapered = 1;
                
                % calculating CS chord. See same code within "if statement"
                % above
                if obj.pos == 1
                    obj.hinge = obj.innerBorder_coords(3);
                    
                elseif obj.pos==0
                    obj.hinge = 1. - obj.innerBorder_coords(3);
                    
                elseif obj.pos==2
                    obj.hinge = obj.innerBorder_coords(4) - obj.innerBorder_coords(3);
                end
                
                % n_hinge states how many hinge lines are present within
                % the same CS. More than one hinge line means that we are
                % dealing with slotted CS, which the current
                % read_from_CPACS cannot handle. Therefore, a warning is
                % thrown, reminding that slotted CS are currently not a
                % feature.
                obj.n_hinge=length(obj.hinge); 
                msg = 'warning: not actually checking amount of hinge lines and presence of slotted flaps. Assuming that no slotted CS are present';
                warning(msg);
                 
             % Throws generic error message if anything unexpected occurs
             % with the dimensions (xsiLE) of the control surface.
             else
                msg = 'Something unexpected occurred with the dimensions of the control surface (xsiLE coordinates of inner and outer border)';
                error(msg);
             end
             

             
             %% ============== Looking at hinge line ==============
             
             % Sets CS deflections to zero (1 element per hinge line)
             obj.delta=zeros(1,obj.n_hinge);
             
             % Saves xsi coords of hinge line (1st element inner coords,
             % 2nd element outer coords). These are the coordinates of the
             % ACTUAL HINGE LINE, unlike the hinge property in the code
             % above, which is just the chord of the CS.
             
             % Daedalus for now does not use them, but storing them might be
             % useful nonetheless (and it is effortless)
             obj.hinge_xsi(1) = tixiGetDoubleElement(tixi,[xpath_cs,'/path/innerHingePoint/hingeXsi']);
             obj.hinge_xsi(2) = tixiGetDoubleElement(tixi,[xpath_cs,'/path/outerHingePoint/hingeXsi']);
             

             
             
             %% ============== Looking at control surface displacement path ==============
             
             % fidning path to steps, initializing variables
             x_path_steps = [xpath_cs,'/path/steps'];
             nr_steps = tixiGetNamedChildrenCount(tixi, x_path_steps, 'step');
             rel_defl = zeros(1,nr_steps);
             hinge_line_rotation = zeros(1,nr_steps);
             
             % Saving relative deflection and hinge line rotation
             % (deflection)
             for step_index = 1:nr_steps
                string_in = [x_path_steps,'/step[',num2str(step_index),']'];
                rel_defl(1,step_index) = tixiGetDoubleElement(tixi, [string_in,'/relDeflection']);
                hinge_line_rotation(1,step_index) = tixiGetDoubleElement(tixi, [string_in,'/hingeLineRotation']); 
             end
             
             % Comparing rotation of all steps to find the minimum and
             % maximum deflections, as defined in CPACS
             
             [~,min_index] = min(rel_defl);
             [~,max_index] = max(rel_defl);
             
             obj.min_defl = hinge_line_rotation(1,min_index);
             obj.max_defl = hinge_line_rotation(1,max_index);   
             
             
             %% ============== Checking for intersection of control surface with wing segments ==============
             
             % gets wing UID, segment UID, and eta and xsi coordinate (in
             % the segment reference frame) of point 1 of the conrtol
             % surface
             [obj.wing_UID, obj.start_segmentUID, obj.start_segmentEta, obj.start_segmentXsi, ~]...
                 = tiglWingComponentSegmentPointGetSegmentEtaXsi(tigl, obj.parent_uID, obj.innerBorder_coords(1), obj.innerBorder_coords(3));
                     
             % gets wing UID, segment UID, and eta and xsi coordinate (in
             % the segment reference frame) of point 2 of the conrtol
             % surface
             [obj.wing_UID, obj.end_segmentUID, obj.end_segmentEta, obj.end_segmentXsi, ~]...
                 = tiglWingComponentSegmentPointGetSegmentEtaXsi(tigl, obj.parent_uID, obj.outerBorder_coords(1), obj.outerBorder_coords(3));
             
             % Gets segment indeces of start and end segments
             [obj.start_segment_index,~] = tiglWingGetSegmentIndex(tigl, obj.start_segmentUID);
             [obj.end_segment_index,~] = tiglWingGetSegmentIndex(tigl, obj.end_segmentUID);
             
             % Gets wing index
             obj.wing_index = tiglWingGetIndex(tigl, obj.wing_UID);
             
             
                          %% ============== Forcing variables which are defined in Daedalus but not in CPACS ==============
             
             % The two variables "is_sym" and "is_sym_defl" below are necessary to run Daedalus
             % properly, but are not defined in CPACS specifically for
             % control surfaces
             
             % If this is_sym refers the same control surface being
             % present on the opposite wing, then it is indeed defined in
             % CPACS. Symmetry in CPACS is defined on a wing-level, so the
             % wing attributes should be checked for
             % <symmetry>true</symmetry>. 
             wing_xpath = tixiUIDGetXPath(tixi, obj.wing_UID);
             
             % tries to get the symmetry attribute from the wing where the
             % CS is located. If it does not exist, a warning message is
             % thrown, and lack of symmetry is assumed.
             try
                wing_sym = tixiGetTextAttribute(tixi, wing_xpath, 'symmetry');
                
                % Checks whether <symmetry> is true or false. If it is
                % neither of them, a warning message is thrown, and lack of
                % symmetry is assumed
                
                if wing_sym == 'true'
                    obj.is_sym=1;
                    
                elseif wing_sym == 'false'
                    obj.is_sym=0;
                    
                else
                    obj.is_sym=0;
                    msg = 'The symmetry attribute of the wing where the CS is located has unexpected content';
                    warning(msg);
                end
                
             % if <symmetry> is not defined at all, lack of symmetry is
             % assumed
             catch
                 obj.is_sym=0;
             end

             % The variable below is not defined in CPACS. It is assumed to
             % be true, and a warning is thrown to point out this
             % occurrence.
             
             obj.is_sym_defl=1;
             msg = 'The variable is_sym_defl is not defined in CPACS, therefore we are assuming it is true (=1)';
                    warning(msg);
        end
    end    
    
    methods
        
        % gets xyz coordinates of points 1,2,3,4 as defined in CPACS
        function xyz_CPACS = get_xyz(obj)
        
        [x,y,z] = tiglWingComponentSegmentGetPoint(obj.tigl_handle, obj.parent_uID, obj.innerBorder_coords(1),obj.innerBorder_coords(3));
        p1 = [x,y,z];
        
        [x,y,z] = tiglWingComponentSegmentGetPoint(obj.tigl_handle, obj.parent_uID, obj.outerBorder_coords(1),obj.outerBorder_coords(3));
        p2 = [x,y,z];
        
        [x,y,z] = tiglWingComponentSegmentGetPoint(obj.tigl_handle, obj.parent_uID, obj.outerBorder_coords(2),obj.outerBorder_coords(4));
        p3 = [x,y,z];
        
        [x,y,z] = tiglWingComponentSegmentGetPoint(obj.tigl_handle, obj.parent_uID, obj.innerBorder_coords(2),obj.innerBorder_coords(4));
        p4 = [x,y,z];
        
        xyz_CPACS(1,:,:) = [p1' p2' p3' p4'];
        
        end        

        function obj=set.delta(obj, def)
            obj.delta=def;
            for i=1:length(obj.children) %#ok<MCSUP>
                obj.children{i}.delta=def; %#ok<MCSUP>
            end
        end
        
        function obj=set.delta_l_r(obj, def_l_r)
            obj.delta_l_r=def_l_r;
            for i=1:length(obj.children) %#ok<MCSUP>
                obj.children{i}.delta_l_r=def_l_r; %#ok<MCSUP>
            end
        end
        
        function dummy=get.panelIds(obj)
            dummy = [];
            for i=1:length(obj.children)
                dummy = [dummy, obj.children{i}.panelIds]; %#ok<*AGROW>
            end
        end
        
        function dummy=get.panelIdsL(obj)
            dummy = [];
            for i=1:length(obj.children)
                dummy = [dummy, obj.children{i}.panelIdsL];
            end
        end
        
        function dummy=get.panelIds_special(obj)
            dummy = [];
            for i=1:length(obj.children)
                dummy = [dummy, obj.children{i}.panelIds_special];
            end
        end
        
        function dummy=get.panelIds_special_L(obj)
            dummy = [];
            for i=1:length(obj.children)
                dummy = [dummy, obj.children{i}.panelIds_special_L];
            end
        end
        
        function dummy=get.panelIds_standard(obj)
            dummy = [];
            for i=1:length(obj.children)
                dummy = [dummy, obj.children{i}.panelIds_standard];
            end
        end
        
        function dummy=get.panelIds_standard_L(obj)
            dummy = [];
            for i=1:length(obj.children)
                dummy = [dummy, obj.children{i}.panelIds_standard_L];
            end
        end
        
        function obj=set.is_sym(obj, sym_cnd)
            obj.is_sym=sym_cnd;
            for i=1:length(obj.children) %#ok<MCSUP>
                obj.children{i}.is_sym=sym_cnd; %#ok<MCSUP>
            end
        end
    end
end



                
                
