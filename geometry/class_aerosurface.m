%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%
classdef class_aerosurface
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
        name;
        isExternalFEM = 0;  % 0 if wing should be selfdesigned in Daedalus, 
                            % 1 if wing is from external FEM model
        pathNodeCoords;     % if wing is from external FEM model, gives path
                            % of node coordinates matlab file. Nodes of 
                            % right half wing from root to tip
        pathMassMatrix;         % if wing is from external FEM model, gives path
                                % of mass matrix matlab file. Mass matrix
                                % of whole wing
        pathStiffnessMatrix;    % if wing is from external FEM model, gives path
                                % of stiffness matrix matlab file. Stiffness matrix
                                % of whole wing
        
        wing_segments;  % array of wing segments
        
        % This struct will hold the following information about this wing
        % segment's structure:
        % fs_segments: [eta_root, xsi_root, eta_tip, xsi_tip,  t_web, t_top, t_bottom]
        % rs_segments: [eta_root, xsi_root, eta_tip, xsi_tip,  t_web, t_top, t_bottom]
        wings_structural_properties;
        
        % for 2D grid
        grid_start_idx;
        panel_start_idx;
        
        % for 3D grid
        grid_start_idx_vol;
        panel_start_idx_vol;
        
        
        c4_coords;
        c4_rl_coords;
        c4_forces;
        c4_moments;
        
        wingbox_coords;
        wingbox_rl_coords;
        wingbox_c4;
        wingbox_height;
        air;
        
        wingbox_rl_coords_mid
        c4_forces_structmesh;
        c4_moments_structmesh;
        
        
        beam_forces_structmesh;
        beam_moments_structmesh;
        
        pos;
        symmetric;
        
        S_ref;      %   wing reference area
        S_wet;      %   wetted area
        S;          %   wing area
        AR;         %   aspect ratio
        TR;         %   taper ratio
        c4sweep;    %   sweep of quarter chord line
        Theta;      %   twist
        b;          %   segment span
        Lamda;      %   quarter chord sweep
        
        c_r;
        c_t;
        c_mac;
        
        CD_0;
        D_f;
        CD_f;
        
        Re;
        grid;
        grid_flat;
        is_te;
        
        grid_vol_lo;
        grid_vol_up;
        grid_deflected;
        panels;
        te_idx;
        
        % volume grid
        grid_vol;
        panels_vol;
        is_te_vol;
        opposite_te_vol;
        te_idx_vol;
        
        grid_wake;
        panels_wake;
        
        T;
        T_sort;
        db;
        db_sort;
        
        panel_to_beam_element;
        
        grid_3D=1;

        %array of control surface class instance
        control_surfaces;
        geom_arr;
        curr_geom_cell;
        curr_etas_cell;
        final_geom_cell;
        final_etas_cell;
        etas_cuts;
        % parametric , wing spanning control surfaces
        parCS
    end
    
    methods (Static)
        
        function [obj, b_ref] = create_from_cpacs(tixi, tigl, wingIndex)
            obj = class_aerosurface();
            
            uID = tiglWingGetUID(tigl, wingIndex);
            
            b_ref = tiglWingGetSpan(tigl, uID);

            x_wing = tixiUIDGetXPath(tixi, uID);
            try
                obj.name = tixiGetTextElement(tixi, [x_wing, '/name']);
            catch
                obj.name = '';
            end

            obj.symmetric = 0;
            if tiglWingGetSymmetry(tigl, wingIndex) ~= 0
                obj.symmetric = 1;
            end
            
            % First we can create the segments. We can split these later if
            % necessary.
            n_segments = tiglWingGetSegmentCount(tigl, wingIndex);
            if n_segments == 1
                obj.wing_segments = ...
                    class_wingsegment.create_from_cpacs(tixi, tigl, wingIndex, 1);
                obj.wing_segments.symmetric = obj.symmetric;
            else
                for i = 1:n_segments
                    segment = ...
                        class_wingsegment.create_from_cpacs(tixi, tigl, wingIndex, i);
                    segment.symmetric = obj.symmetric;

                    if isempty(obj.wing_segments)
                        obj.wing_segments = segment;
                    else
                        obj.wing_segments(end+1) = segment;
                    end
                end
            end

            n_compSegments = tiglWingGetComponentSegmentCount(tigl, wingIndex);
            obj.wings_structural_properties.fs_segments = [];
            obj.wings_structural_properties.rs_segments = [];
            obj.wings_structural_properties.frontspar = [];
            obj.wings_structural_properties.rearspar = [];
            obj.wings_structural_properties.is_fueled = [];
            obj.wings_structural_properties.material = {};
            
            % Saves original amount of segments, before daedalus starts
            % adding them. Used later on
            obj.geom_arr = ones(1,length(obj.wing_segments));
            
            % initializes array which will be later use to store segment
            % indices without losing info when the for loop below changes
            % cycle
            segmentIndices_updated = [];
            
            % loop over all component segments
            for i = 1:n_compSegments
                compSegUID = tiglWingGetComponentSegmentUID(tigl, wingIndex, i);
                x_compSeg = tixiUIDGetXPath(tixi, compSegUID);
                n_segs_in_compSeg = tiglWingComponentSegmentGetNumberOfSegments(tigl, compSegUID);

                % check if a structure is defined for this segment
                if tixiGetNamedChildrenCount(tixi, x_compSeg, 'structure')
                    x_struct = [x_compSeg, '/structure'];

                    % check if there are spars for this segment
                    if tixiGetNamedChildrenCount(tixi, x_struct, 'spars')
                        x_spSegs = [x_struct, '/spars/sparSegments'];
                        n_spSegs = tixiGetNamedChildrenCount(tixi, x_spSegs, 'sparSegment');

                        segmentUIDs = cell(1, n_segs_in_compSeg);
                        innerElemUIDs = cell(1, n_segs_in_compSeg);
                        outerElemUIDs = cell(1, n_segs_in_compSeg);
                        sectionEtas = nan(2*(n_segs_in_compSeg + 1), 1);
                        segmentIndices = zeros(1, n_segs_in_compSeg);
                        


                        for j = 1:n_segs_in_compSeg
                            segmentUIDs{j} = tiglWingComponentSegmentGetSegmentUID(tigl, compSegUID, j);
                            [segmentIndices(j), wingIndex] = tiglWingGetSegmentIndex(tigl, segmentUIDs{j});

                            [~, innerElemUIDs{j}] = tiglWingGetInnerSectionAndElementUID(tigl, wingIndex, segmentIndices(j));
                            [~, outerElemUIDs{j}] = tiglWingGetOuterSectionAndElementUID(tigl, wingIndex, segmentIndices(j));
                            [eta_inner, ~] = tiglWingSegmentPointGetComponentSegmentEtaXsi(tigl, segmentUIDs{j}, compSegUID, 0, 0);  
                            [eta_outer, ~] = tiglWingSegmentPointGetComponentSegmentEtaXsi(tigl, segmentUIDs{j}, compSegUID, 1, 0);

                            sectionEtas(2*j - 1) = eta_inner;
                            sectionEtas(2*j - 0) = eta_outer;
                        end
                        
                        % If this is the first time the loop is running
                        % (only 1 component segment) then
                        % segmentIndices_updated is set equal to
                        % segmentIndices, as there cannot have been any
                        % loss of info yet due to the for loop advancing 1
                        % cycle
                        if isempty(segmentIndices_updated) == 1
                            segmentIndices_updated = segmentIndices;
                            prev_index = 0;
                            
                        % If this is not the 1st time the loop is running,
                        % segmentIndices is not up to date anymore, hence segmentIndices_updated
                        % references itself
                        else
                            prev_index = segmentIndices_updated(end);
                            segmentIndices_updated = 1:1:segmentIndices_updated(end)+n_segs_in_compSeg;
                        end
                        
                        
                        sectionEtas = unique(sectionEtas(~isnan(sectionEtas)));

                        % spSegs will encode 2 point spar segments in its rows
                        % as follows: [eta_1, xsi_1, eta_2, xsi_2, t_web, t_top, t_bottom]
                        spSegs = zeros(0, 7);
                        points = zeros(0, 2);

                        for j = 1:n_spSegs
                            x_spSeg = sprintf('%s/sparSegment[%i]', x_spSegs, j);
                            n_spPos = tixiGetNamedChildrenCount(tixi, [x_spSeg, '/sparPositionUIDs'], 'sparPositionUID');

                            pos_ = zeros(n_spPos, 2);
                            % loop over all positions and store in pos array
                            for k = 1:n_spPos
                                x_spPosUID = sprintf('%s/sparPositionUIDs/sparPositionUID[%i]', x_spSeg, k);
                                spPosUID = tixiGetTextElement(tixi, x_spPosUID);
                                x_spPos = tixiUIDGetXPath(tixi, spPosUID);

                                xsi = tixiGetDoubleElement(tixi, [x_spPos, '/xsi']);

                                % Obtaining the eta is a little more
                                % complicated
                                if tixiGetNamedChildrenCount(tixi, x_spPos, 'eta') == 0
                                    elemUID = tixiGetTextElement(tixi, [x_spPos, '/elementUID']);

                                    index = find(strcmp(innerElemUIDs, elemUID));
                                    if ~isempty(index)
                                        eta = 0;
                                    else
                                        index = find(strcmp(outerElemUIDs, elemUID));
                                        eta = 1;
                                    end

                                    [eta, ~] = tiglWingSegmentPointGetComponentSegmentEtaXsi(tigl, segmentUIDs{index}, compSegUID, eta, xsi);
                                    [~, id_sec] = min(abs(sectionEtas - eta));
                                    eta = sectionEtas(id_sec);
                                else
                                    eta = tixiGetDoubleElement(tixi, [x_spPos, '/eta']);
                                    d_eta = abs(sectionEtas - eta);
                                    id_sec = find(round(d_eta, 2) == 0);
                                    if ~isempty(id_sec)
                                        warning('A spar was positioned very close to a section (eta difference = %.6f). Spar will be put exactly at this section.', d_eta);
                                        eta = sectionEtas(id_sec(1));
                                    end
                                end

                                point = [eta, xsi];
                                pos_(k, :) = point;
                                points = [points; point];
                            end

                            ts = ones(1,3);            
                            % Obtain the thicknesses of the wingbox
                            x_spCS = [x_spSeg, '/sparCrossSection'];
                            if tixiCheckElement(tixi, x_spCS)
                                x_t_web = [x_spCS, '/web1/material/thickness'];
                                if tixiCheckElement(tixi, x_t_web)
                                    ts(1) = tixiGetDoubleElement(tixi, x_t_web);
                                end

                                x_t_top = [x_spCS, '/upperCap/material/thickness'];
                                if tixiCheckElement(tixi, x_t_top)
                                    ts(2) =  tixiGetDoubleElement(tixi, x_t_top);
                                end

                                x_t_bot = [x_spCS, '/lowerCap/material/thickness'];
                                if tixiCheckElement(tixi, x_t_bot)
                                    ts(3) = tixiGetDoubleElement(tixi, x_t_bot);
                                end
                            end

                            pos = [pos_(1:end-1, :), pos_(2:end, :)];
                            spSeg = [pos, repmat(ts, size(pos, 1), 1)];
                            spSegs = [spSegs; spSeg];
                        end

                        % To construct front and rear spars we loop
                        % over all points and find the most forward and
                        % aft points at their eta. The most forward
                        % point becomes a node in the front spar, the
                        % most aft becomes a node in the rear spar. At
                        % these nodes we note the eta, xsi, and the
                        % thicknesses of the corresponding segment. We
                        % store these nodes in two matrices, one for
                        % the front and another for the rear spar.
                        % These matrices will be have size N x 5 such
                        % that their rows are:
                        % [eta, xsi, t_web, t_top, t_bottom]

                        % We should present the user with a warning if
                        % there are more than 2 points on a line,
                        % because this means we are removing one or
                        % possibly more middle spars.

                        % Furthermore, we should check for each node of
                        % the newly constructed front and rear spars if
                        % it lies on a boundary between two wing
                        % segments. If it doesn't we should split the
                        % wing segment at this node.

                        etas = unique(points(:, 1));
                        n_points = size(etas, 1);

                        nodes_fs_ = zeros(n_points, 5);
                        nodes_fs_(:, 1) = etas;
                        nodes_rs_ = nodes_fs_;

                        % We do this loop in reverse to aid the process
                        % of splitting segments.
                        for j = n_points:-1:1
                            % Only consider the spar segments that
                            % start at, end at, or cross the eta of
                            % the current point.
                            segs_ = spSegs(spSegs(:, 1) <= etas(j) & spSegs(:, 3) >= etas(j), :);

                            % Calculate the xsi values at the locations
                            % at which a vertical line through the
                            % current eta intersects each of the
                            % segments under consideration.
                            xsis = segs_(:, 2) + (segs_(:, 4) - segs_(:, 2))./(segs_(:, 3) - segs_(:, 1)).*(etas(j) - segs_(:, 1));

                            % Check if there are more than 2 spar
                            % segments that cross this eta. Give a
                            % warning to let the user know we are
                            % removing these segments here.
                            if length(xsis) > 2
                                warning('There are more than two spars specified at a spanwise location. Only the most forward and aft will be kept.');
                            end

                            % The minimum xsi in this vector
                            % corresponds to the front, and the maximum
                            % to the rear spar.
                            [nodes_fs_(j, 2), id_fs] = min(xsis);
                            [nodes_rs_(j, 2), id_rs] = max(xsis);

                            % Then, store the thicknesses of the
                            % corresponding spar segments alongside the
                            % new nodes.
                            nodes_fs_(j, 3:end) = segs_(id_fs, 5:end);
                            nodes_rs_(j, 3:end) = segs_(id_rs, 5:end);

                            % Finally, check if this node lies
                            % precisely on the boundary between two
                            % wing segments. If it doesn't, split the
                            % corresponding wing segment at this node.
                            if all(etas(j) ~= sectionEtas)
                                % First get the id of the segment we're
                                % gonna split.
                                seg_num = find(etas(j) < sectionEtas, 1) - 1;
                                seg_id = segmentIndices_updated(seg_num);

                                % Now we need to find the eta of the
                                % split w.r.t. the segment.
                                inner_section_eta = sectionEtas(seg_num);
                                outer_section_eta = sectionEtas(seg_num + 1);                          
                                eta_seg = (etas(j) - inner_section_eta)/(outer_section_eta - inner_section_eta);

                                % Then we can split the segment there
                                % and add the current node eta to the
                                % list of section etas.
                                obj = obj.split_segment(seg_id, eta_seg, 1);
                                sectionEtas = [sectionEtas(1:seg_num); etas(j); sectionEtas((seg_num+1):end)];
                                n_segs_in_compSeg = n_segs_in_compSeg + 1;
                                n_segments = n_segments + 1;
                                segmentIndices_updated = [segmentIndices_updated,segmentIndices_updated(end) + 1];
                                segmentIndices = [segmentIndices(1:seg_num), segmentIndices(seg_num:end)+1];
                            end
                        end

                        % Now we need to make sure that there are spar
                        % positions at each section. 
                        nodes_fs = zeros(n_segs_in_compSeg + 1, 5);
                        nodes_fs(2:end-1, 1) = sectionEtas(2:end-1);
                        nodes_fs(end, 1) = 1;
                        nodes_rs = nodes_fs;
                        for j = 1:4
                            nodes_fs(:, j+1) = interp1(nodes_fs_(:, 1), nodes_fs_(:, j+1), nodes_fs(:, 1), 'linear', 'extrap');
                            nodes_rs(:, j+1) = interp1(nodes_rs_(:, 1), nodes_rs_(:, j+1), nodes_rs(:, 1), 'linear', 'extrap');
                        end

                        % Now we have to turn the lists of nodes back
                        % into lists of spar segments with two end
                        % points.
                        fs_segments = [nodes_fs(1:end-1, 1:2), nodes_fs(2:end, 1:2), (nodes_fs(1:end-1, 3:end) + nodes_fs(2:end, 3:end))/2];
                        rs_segments = [nodes_rs(1:end-1, 1:2), nodes_rs(2:end, 1:2), (nodes_rs(1:end-1, 3:end) + nodes_rs(2:end, 3:end))/2];

                        % Finally, we can store the segments where they
                        % belong.
                        for j = 1:n_segs_in_compSeg
                            obj.wing_segments(segmentIndices_updated(j+prev_index)).structural_properties.fs_segments = fs_segments(j, [2, 4:end]);
                            obj.wing_segments(segmentIndices_updated(j+prev_index)).structural_properties.rs_segments = rs_segments(j, [2, 4:end]);
                        end
                        
                        % Note that prev_index is used above in order to give
                        % the correct starting index in case the for loop
                        % is at its 2nd or higher cycle

                        obj.wings_structural_properties.fs_segments = [obj.wings_structural_properties.fs_segments; fs_segments];
                        obj.wings_structural_properties.rs_segments = [obj.wings_structural_properties.rs_segments; rs_segments];

                        obj.wings_structural_properties.frontspar = [obj.wings_structural_properties.frontspar, nodes_fs(:, 2)'];
                        obj.wings_structural_properties.rearspar = [obj.wings_structural_properties.rearspar, nodes_rs(:, 2)'];
                    end
                end
            end        

            % loop over all component segments
            control_counter = 1;
            for i = 1:n_compSegments
                
                % Set component segment UID and xml path
                compSegUID = tiglWingGetComponentSegmentUID(tigl, wingIndex, i);
                x_compSeg = tixiUIDGetXPath(tixi, compSegUID);

                % checks if control surfaces are defined within this segment. If yes, saves X-path
                if tixiGetNamedChildrenCount(tixi, x_compSeg, 'controlSurfaces') ~= 0
                    x_CS = [x_compSeg, '/controlSurfaces'];


                    % checks if trailing edge control surfaces exist. If yes, saves X-path
                    if tixiGetNamedChildrenCount(tixi, x_CS, 'trailingEdgeDevices') ~= 0
                        x_CS_TE = [x_CS, '/trailingEdgeDevices'];

                        % checks how many TE_CS there are
                        nr_CS_TE = tixiGetNamedChildrenCount(tixi, x_CS_TE, 'trailingEdgeDevice');

                        % loop over every single TE_CS within TE_CS block
                        for CS_TE_index = 1:nr_CS_TE

                            %writes temp xpath to single TE_CS
                            string_in = ['/trailingEdgeDevice[',num2str(CS_TE_index),']'];
                            x_CS_TE_temp = [x_CS_TE,string_in ];

                            %creates CS instance in CS class (pos=0 since TE)
                            control_surfaces_temp = class_control_surface_parent.create_from_cpacs(tixi,tigl,x_CS_TE_temp,0);
                            control_surfaces(control_counter) = control_surfaces_temp;
                            control_counter = control_counter + 1;
                        end  
                    end


                    % checks if leading edge control surfaces exist. If yes, saves X-path
                    if tixiGetNamedChildrenCount(tixi, x_CS, 'leadingEdgeDevices') ~= 0
                        x_CS_LE = [x_CS, '/leadingEdgeDevices'];

                        % checks how many LE_CS there are
                        nr_CS_LE = tixiGetNamedChildrenCount(tixi, x_CS_LE, 'leadingEdgeDevice');

                        % loop over every single LE_CS within LE_CS block
                        for CS_LE_index = 1:nr_CS_LE

                            %writes temp xpath to single LE_CS
                            string_in = ['/leadingEdgeDevice[',num2str(CS_LE_index),']'];
                            x_CS_LE_temp = [x_CS_LE,string_in];

                            %creates CS instance in CS class (pos=1 since LE)
                            control_surfaces_temp = class_control_surface_parent.create_from_cpacs(tixi,tigl,x_CS_LE_temp,1);
                            control_surfaces(control_counter) = control_surfaces_temp;
                            control_counter = control_counter + 1;

                        end 
                    end
                    
                    
                    % checks if spoiler control surfaces exist. If yes, saves X-path
                    if tixiGetNamedChildrenCount(tixi, x_CS, 'spoilers') ~= 0
                        x_CS_SP = [x_CS, '/spoilers'];

                        % checks how many TE_CS there are
                        nr_CS_SP = tixiGetNamedChildrenCount(tixi, x_CS_SP, 'spoiler');

                        % loop over every single TE_CS within TE_CS block
                        for CS_SP_index = 1:nr_CS_SP

                            %writes temp xpath to single TE_CS
                            string_in = ['/spoiler[',num2str(CS_SP_index),']'];
                            x_CS_SP_temp = [x_CS_SP,string_in ];

                            %creates CS instance in CS class (pos=0 since TE)
                            control_surfaces_temp = class_control_surface_parent.create_from_cpacs(tixi,tigl,x_CS_SP_temp,2);
                            control_surfaces(control_counter) = control_surfaces_temp;
                            control_counter = control_counter + 1;
                        end  
                    end

                    
                    % saves all the control surfaces we saw so far as a
                    % class property (class_aerosurface)
                    obj.control_surfaces = control_surfaces;
                    
                    % looks at the location of all control surfaces
                    % (belonging to the same component_segment) on this
                    % aerosurface, then calculates exactly where the wing
                    % will have to be cut in order to generate all the
                    % segments necessary to perfectly contain all control
                    % surfaces
                    obj = calculate_final_segmentation(obj);
                    
                    % looks at the output of calculate_final_segmentation,
                    % and actually does all the splitting, resulting in
                    % wing_segments being created through split_segment.
                    % All of them are of course associated to this specific
                    % instance (obj) of aerosurface
                    obj = wing_segmentation(obj);  
                    
                    % looks at the output of wing_segmentation (so the
                    % current wing_segment disposition) and assigns each
                    % control surface to each spefic wing_segment
                    obj = assign_control_surfaces(obj);
                    
                end
            end
            n_segments = length(obj.wing_segments);
            obj.wings_structural_properties.is_fueled = ones(1, n_segments);
            obj.wings_structural_properties.material = repmat({'aluminum'}, 1, n_segments);        
        end
    end
    
    methods
        % constructor
        function obj = class_aerosurface(varargin)
            if nargin==1
                obj=obj.read_xml_definition(varargin{1});
            elseif nargin==19
                pos=varargin{1};
                symmetric=varargin{2};
                dihed=varargin{3};
                LambdaSpec=varargin{4};
                Lambda=varargin{5};
                property1=varargin{6};
                property2=varargin{8};
                property3=varargin{10};
                property4=varargin{12};
                property5=varargin{14};
                property6=varargin{16};
                property7=varargin{18};
                value1=varargin{7};
                value2=varargin{9};
                value3=varargin{11};
                value4=varargin{13};
                value5=varargin{15};
                value6=varargin{17};
                value7=varargin{19};
                obj.wing_segments=class_wingsegment(pos,symmetric, dihed,LambdaSpec,Lambda,property1,value1,property2,value2,property3,value3,property4,value4,property5,value5,property6,value6,property7,value7);
            elseif nargin==21
                pos=varargin{1};
                symmetric=varargin{2};
                dihed=varargin{3};
                LambdaSpec=varargin{4};
                Lambda=varargin{5};
                property1=varargin{6};
                property2=varargin{8};
                property3=varargin{10};
                property4=varargin{12};
                property5=varargin{14};
                property6=varargin{16};
                property7=varargin{18};
                property8=varargin{20};
                value1=varargin{7};
                value2=varargin{9};
                value3=varargin{11};
                value4=varargin{13};
                value5=varargin{15};
                value6=varargin{17};
                value7=varargin{19};
                value8=varargin{21};
                obj.wing_segments=class_wingsegment(pos,symmetric, dihed,LambdaSpec,Lambda,property1,value1,property2,value2,property3,value3,property4,value4,property5,value5,property6,value6,property7,value7,property8,value8);
            elseif nargin ~= 0
                obj.wing_segments=class_wingsegment(pos,symmetric, dihed,LambdaSpec,Lambda);
            end
            
            if nargin ~= 0
            
            %             obj.c_mac=cell2mat(obj.wing_segments.c_mac);
            %             obj.S=cell2mat(obj.wing_segments.S);
            obj.symmetric=obj.wing_segments(1).symmetric;
            %             property1=varargin{1};
            %             property2=varargin{3};
            %             property3=varargin{5};
            %
            %
            %             if(strcmp(property1,'SurfaceArea'))&&(strcmp(property2,'AR'))&&(strcmp(property3,'TR'))
            %
            %                 obj.S=varargin{2};
            %                 obj.S_ref=obj.S;
            %                 obj.AR=varargin{4};
            %                 obj.TR=varargin{6};
            %                 obj.S_wet=2*obj.S;
            %                 obj.b=sqrt(obj.S*obj.AR);
            %                 obj.c_r=2*obj.S/((1+obj.TR)*obj.b);
            %                 obj.c_t=obj.TR*obj.c_r;
            %
            %                 obj.c_mac=8*obj.S/((1+obj.TR)^2*obj.b^2)*(obj.b/2-(1-obj.TR)*obj.b/2+(1-obj.TR)^2*obj.b/6);
            %             end
            end
            
        end
        
        
        function obj =add_segment(obj,symmetric,dihed,LambdaSpec,Lambda,varargin)
            if nargin==19
                property1=varargin{1};
                property2=varargin{3};
                property3=varargin{5};
                property4=varargin{7};
                property5=varargin{9};
                property6=varargin{11};
                property7=varargin{13};
                
                value1=varargin{2};
                value2=varargin{4};
                value3=varargin{6};
                value4=varargin{8};
                value5=varargin{10};
                value6=varargin{12};
                value7=varargin{14};
                
                new_segment=class_wingsegment(obj.wing_segments(end).get_tip_c4point,symmetric, dihed,LambdaSpec, Lambda,property1,value1,property2,value2,property3,value3,property4,value4,property5,value5,property6,value6,property7,value7);
            elseif nargin==20
                
                property1=varargin{1};
                property2=varargin{3};
                property3=varargin{5};
                property4=varargin{7};
                property5=varargin{9};
                property6=varargin{11};
                property7=varargin{13};
                property8=varargin{15};
                
                value1=varargin{2};
                value2=varargin{4};
                value3=varargin{6};
                value4=varargin{8};
                value5=varargin{10};
                value6=varargin{12};
                value7=varargin{14};
                
                new_segment=class_wingsegment(property8,symmetric, dihed,LambdaSpec, Lambda,property1,value1,property2,value2,property3,value3,property4,value4,property5,value5,property6,value6,property7,value7);
                
            elseif nargin==21
                property1=varargin{1};
                property2=varargin{3};
                property3=varargin{5};
                property4=varargin{7};
                property5=varargin{9};
                property6=varargin{11};
                property7=varargin{13};
                property8=varargin{15};
                
                value1=varargin{2};
                value2=varargin{4};
                value3=varargin{6};
                value4=varargin{8};
                value5=varargin{10};
                value6=varargin{12};
                value7=varargin{14};
                value8=varargin{16};
                
                new_segment=class_wingsegment(obj.wing_segments(end).get_tip_c4point,symmetric, dihed,LambdaSpec,Lambda,property1,value1,property2,value2,property3,value3,property4,value4,property5,value5,property6,value6,property7,value7,property8,value8);
            elseif nargin==23
                property1=varargin{1};
                property2=varargin{3};
                property3=varargin{5};
                property4=varargin{7};
                property5=varargin{9};
                property6=varargin{11};
                property7=varargin{13};
                property8=varargin{15};
                property9=varargin{17};
                
                value1=varargin{2};
                value2=varargin{4};
                value3=varargin{6};
                value4=varargin{8};
                value5=varargin{10};
                value6=varargin{12};
                value7=varargin{14};
                value8=varargin{16};
                value9=varargin{18};
                
                new_segment=class_wingsegment(obj.wing_segments(end).get_tip_c4point,symmetric, dihed,LambdaSpec,Lambda,property1,value1,property2,value2,property3,value3,property4,value4,property5,value5,property6,value6,property7,value7,property8,value8,property9,value9);
            end
            obj.wing_segments=[obj.wing_segments new_segment];
            obj.c_mac=0;
            obj.b=0;
            for i=1:length(obj.wing_segments)
                obj.c_mac=obj.c_mac+obj.wing_segments(i).c_mac*obj.wing_segments(i).b;
                obj.b=obj.b+obj.wing_segments(i).b;
            end
            obj.c_mac=obj.c_mac/obj.b;
            obj.S=obj.S+obj.wing_segments(end).S;
        end
        
        function obj =split_segment(obj,segment_idx,f_1,f_2,varargin)
            
            xyz=obj.wing_segments(segment_idx).xyz;
            p1n=xyz(:,1)+f_1*(xyz(:,2)-xyz(:,1));
            p4n=xyz(:,4)+f_1*(xyz(:,3)-xyz(:,4));
            
            p2n=xyz(:,1)+f_2*(xyz(:,2)-xyz(:,1));
            p3n=xyz(:,4)+f_2*(xyz(:,3)-xyz(:,4));
            symm=obj.symmetric;
            
            if nargin==4
                if f_2==1
                    profile_r=obj.wing_segments(segment_idx).get_profile(f_1);
                    profile_t=obj.wing_segments(segment_idx).profile_t;
                    s1=class_wingsegment([xyz(:,1) p1n p4n xyz(:,4)],symm,profile_r,profile_t);
                    s2=class_wingsegment([p1n p2n p3n p4n],symm,profile_r,profile_t);
                elseif f_1==0
                    profile_r1=obj.wing_segments(segment_idx).profile_r;
                    profile_t1=obj.wing_segments(segment_idx).get_profile(f_2);
                    profile_r2=profile_t1;
                    profile_t2=obj.wing_segments(segment_idx).profile_t;
                    s2=class_wingsegment([xyz(:,1) p2n p3n xyz(:,4)],symm,profile_r1,profile_t1);
                    s3=class_wingsegment([p2n xyz(:,2) xyz(:,3) p3n],symm,profile_r2,profile_t2);
                else
                    profile_r=obj.wing_segments(segment_idx).get_profile(f_1);
                    profile_t=obj.wing_segments(segment_idx).get_profile(f_2);
                    s1=class_wingsegment([xyz(:,1) p1n p4n xyz(:,4)],symm,profile_r,profile_t);
                    s2=class_wingsegment([p1n p2n p3n p4n],symm,profile_r,profile_t);
                    s3=class_wingsegment([p2n xyz(:,2) xyz(:,3) p3n],symm,profile_r,profile_t);          
                    
                end
            elseif nargin==5
                if f_2==1
                    profile_r=obj.wing_segments(segment_idx).get_profile(f_1);
                    profile_t=obj.wing_segments(segment_idx).profile_t;
                    s1=class_wingsegment([xyz(:,1) p1n p4n xyz(:,4)],symm,profile_r,profile_t);
                    s2=class_wingsegment([p1n p2n p3n p4n],symm,profile_r,profile_t);
                    s2=s2.add_control_surface(varargin{1});
                else
                    profile_r=obj.wing_segments(segment_idx).get_profile(f_1);
                    profile_t=obj.wing_segments(segment_idx).get_profile(f_2);
                    s1=class_wingsegment([xyz(:,1) p1n p4n xyz(:,4)],symm,profile_r,profile_t);
                    s2=class_wingsegment([p1n p2n p3n p4n],symm,profile_r,profile_t);
                    s3=class_wingsegment([p2n xyz(:,2) xyz(:,3) p3n],symm,profile_r,profile_t);
                    s2=s2.add_control_surface(varargin{1});
                end
            elseif nargin==6
                x=varargin{1};
                y=varargin{2};
                if f_2==1
                    profile_r=obj.wing_segments(segment_idx).get_profile(f_1);
                    profile_t=obj.wing_segments(segment_idx).profile_t;
                    s1=class_wingsegment([xyz(:,1) p1n p4n xyz(:,4)],symm,profile_r,profile_t);
                    s2=class_wingsegment([p1n p2n p3n p4n],symm,x,y,profile_r,profile_t);
                else
                    profile_r=obj.wing_segments(segment_idx).get_profile(f_1);
                    profile_t=obj.wing_segments(segment_idx).get_profile(f_2);
                    s1=class_wingsegment([xyz(:,1) p1n p4n xyz(:,4)],symm,profile_r,profile_t);
                    s2=class_wingsegment([p1n p2n p3n p4n],symm,x,y,profile_r,profile_t);
                    s3=class_wingsegment([p2n xyz(:,2) xyz(:,3) p3n],symm,profile_r,profile_t);
                    
                end
                
            % If nargin==7, the split+segment function takes the structural
            % properties (front and rear spar) of the segment, and
            % interpolates them so as to find the equivalent ones for all
            % children segments.(nargin variabes 5 and 6 and 7 are dummy
            % variables necessary to trigger nargin==7)
            
            elseif nargin==7
                
                if f_2==1
                    
                    %Code below same as nargin==4
                    profile_r=obj.wing_segments(segment_idx).get_profile(f_1);
                    profile_t=obj.wing_segments(segment_idx).profile_t;
                    s1=class_wingsegment([xyz(:,1) p1n p4n xyz(:,4)],symm,profile_r,profile_t);
                    s2=class_wingsegment([p1n p2n p3n p4n],symm,profile_r,profile_t);
                    %Code above same as nargin==4
                    
                    parent_properties = obj.wing_segments(segment_idx).structural_properties;
                    
                    % xsi coordinate of front spar and rear spar at the eta of split are calculated through interpolation 
                    xsi_interp_fs = interp1([0,1],parent_properties.fs_segments(1:2),f_1);
                    xsi_interp_rs = interp1([0,1],parent_properties.rs_segments(1:2),f_1);
                    
                    
                    % structural properties of children segments are
                    % initialized by copying parent
                    s1.structural_properties = parent_properties;
                    s2.structural_properties = parent_properties;
                    
                    % structural_properties of children segments are
                    % modified to incorporate interpolated xsi.
                    s1.structural_properties.fs_segments(2) = xsi_interp_fs;
                    s1.structural_properties.rs_segments(2) = xsi_interp_rs;
                    
                    s2.structural_properties.fs_segments(1) = xsi_interp_fs;
                    s2.structural_properties.rs_segments(1) = xsi_interp_rs;
                    
                    
                else
                    
                    %Code below same as nargin==4
                    profile_r=obj.wing_segments(segment_idx).get_profile(f_1);
                    profile_t=obj.wing_segments(segment_idx).get_profile(f_2);
                    s1=class_wingsegment([xyz(:,1) p1n p4n xyz(:,4)],symm,profile_r,profile_t);
                    s2=class_wingsegment([p1n p2n p3n p4n],symm,profile_r,profile_t);
                    s3=class_wingsegment([p2n xyz(:,2) xyz(:,3) p3n],symm,profile_r,profile_t);      
                    %Code above same as nargin==4
                    
                    
                    % Code below does exactly the same of its equivalent
                    % shown within the if statement above. The only
                    % difference is that now each spar is divided into 3
                    % segments, so a total of 4 xsi need to be
                    % interpolated
                    parent_properties = obj.wing_segments(segment_idx).structural_properties;
                    
                    xsi_interp_fs_1 = interp1([0,1],parent_properties.fs_segments(1:2),f_1);
                    xsi_interp_rs_1 = interp1([0,1],parent_properties.rs_segments(1:2),f_1);
                    xsi_interp_fs_2 = interp1([0,1],parent_properties.fs_segments(1:2),f_2);
                    xsi_interp_rs_2 = interp1([0,1],parent_properties.rs_segments(1:2),f_2);
                    
                    s1.structural_properties = parent_properties;
                    s2.structural_properties = parent_properties;
                    s3.structural_properties = parent_properties;
                    
                    s1.structural_properties.fs_segments(2) = xsi_interp_fs_1;
                    s1.structural_properties.rs_segments(2) = xsi_interp_rs_1;
                    
                    s2.structural_properties.fs_segments(1) = xsi_interp_fs_1;
                    s2.structural_properties.rs_segments(1) = xsi_interp_rs_1;
                    s2.structural_properties.fs_segments(2) = xsi_interp_fs_2;
                    s2.structural_properties.rs_segments(2) = xsi_interp_rs_2;
                    
                    
                    s3.structural_properties.fs_segments(1) = xsi_interp_fs_2;
                    s3.structural_properties.rs_segments(1) = xsi_interp_rs_2;
                    
                end
                
                
            end
            if f_2==1
                if segment_idx+1<=length(obj.wing_segments)
                    obj.wing_segments=[obj.wing_segments(1:segment_idx-1) s1 s2 obj.wing_segments(segment_idx+1:end)];
                    obj = save_geometry_change(obj,segment_idx,f_1,f_2, 1);
                else
                    obj.wing_segments=[obj.wing_segments(1:segment_idx-1) s1 s2 ];
                    obj = save_geometry_change(obj,segment_idx,f_1,f_2, 1);
                end
            elseif f_1==0
                obj.wing_segments=[obj.wing_segments(1:segment_idx-1) s2 s3 obj.wing_segments(segment_idx+1:end)];
                obj = save_geometry_change(obj,segment_idx,f_1,f_2, 1);
            else
                obj.wing_segments=[obj.wing_segments(1:segment_idx-1) s1 s2 s3 obj.wing_segments(segment_idx+1:end)];
                obj = save_geometry_change(obj,segment_idx,f_1,f_2, 2);
            end
            obj.c_mac=0;
            for i=1:length(obj.wing_segments)
                obj.c_mac=obj.c_mac+obj.wing_segments(i).c_mac/length(obj.wing_segments);
            end
        end
        
        function obj = compute_friction_drag(obj,state,S_ref)
            obj.CD_f=0;
            obj.D_f=0;
            for i=1:1:length(obj.wing_segments)
                obj.wing_segments(i)=obj.wing_segments(i).compute_friction_drag(state,S_ref);
                obj.D_f=obj.D_f+obj.wing_segments(i).D_f;
                obj.CD_f=obj.CD_f+obj.wing_segments(i).CD_f;
            end
        end
        
        function obj= set_cs_deflections(obj,varargin)
            if nargin==4
                idx=varargin{2};
                delta=varargin{3};
                
                for i=1:length(idx)
                    if strcmp(varargin{1},'te')
                        obj.wing_segments(idx(i))=obj.wing_segments(idx(i)).set_cs_deflections('te',delta(i,:));
                    elseif strcmp(varargin{1},'le')
                        obj.wing_segments(idx(i))=obj.wing_segments(idx(i)).set_cs_deflections('le',delta(i,:));
                    end
                end
                
            elseif nargin==7
                idx1=varargin{2};
                delta1=varargin{3};
                
                idx2=varargin{5};
                delta2=varargin{6};
                
                if strcmp(varargin{1},'te') && strcmp(varargin{4},'le')
                    for i=1:length(idx1)
                        obj.wing_segments(idx1(i))=obj.wing_segments(idx1(i)).set_cs_deflections('te',delta1(i,:));
                    end
                    
                    for i=1:length(idx2)
                        obj.wing_segments(idx2(i))=obj.wing_segments(idx2(i)).set_cs_deflections('le',delta2(i,:));
                    end
                    
                elseif strcmp(varargin{1},'le') && strcmp(varargin{4},'te')
                    for i=1:length(idx1)
                        obj.wing_segments(idx1(i))=obj.wing_segments(idx1(i)).set_cs_deflections('le',delta1(i,:));
                    end
                    
                    for i=1:length(idx2)
                        obj.wing_segments(idx2(i))=obj.wing_segments(idx2(i)).set_cs_deflections('te',delta2(i,:));
                    end
                end
            end
        end
        
        function obj = plot(obj)
            for i=1:length(obj.wing_segments)
                obj.wing_segments(i).plot_segment();
                %obj.wing_segments(i).compute_grid();
                hold on
            end
        end
        
        function obj=update_idx(obj)
            for i=1:length(obj.wing_segments)
                obj.wing_segments(i).grid_start_idx=obj.wing_segments(i).grid_start_idx+(obj.grid_start_idx-1);
                obj.wing_segments(i).panel_start_idx=obj.wing_segments(i).panel_start_idx+(obj.panel_start_idx-1);
            end
        end
        
        function obj=compute_bspline_grid(obj,x_max,y_max)
            
            
        end
        
        function obj=compute_grid(obj,x_max,y_max,n_x_min,wake)
            
            grid=[];
            grid_flat=[];
            panels=[];
            grid_len_b4=[];
            grid_vol_up=[];
            grid_vol_lo=[];
            is_te=[];
            
            te_idx=[];
            
            if (wake==1)||(wake==2)
                grid_wake=[];
                panels_wake=[];
            end
            
            for i=1:length(obj.wing_segments)
                obj.wing_segments(i)=obj.wing_segments(i).left_control_surfaces();
            end
            
            for i=1:length(obj.wing_segments)
                obj.wing_segments(i)=obj.wing_segments(i).compute_grid(x_max,y_max,n_x_min,wake);
                grid_len_b4=length(grid);
                grid=[grid obj.wing_segments(i).grid];
                grid_flat=[grid_flat obj.wing_segments(i).grid_flat];
                
                grid_vol_up=[grid_vol_up obj.wing_segments(i).grid_vol_upper];
                grid_vol_lo=[grid_vol_lo obj.wing_segments(i).grid_vol_lower];
                
                
                is_te=[is_te obj.wing_segments(i).is_te];
                
                panel_len_b4=size(panels,2);
                te_idx=[te_idx, obj.wing_segments(i).te_idx+grid_len_b4];
                panels=[panels, obj.wing_segments(i).panels+grid_len_b4];
                panel_len_aft=length(panels);
                obj.wing_segments(i).grid_start_idx=grid_len_b4+1;
                obj.wing_segments(i).panel_start_idx=panel_len_b4+1;
                
                if (wake==1)||(wake==2)
                    grid_wake_len_b4=length(grid_wake);
                    grid_wake=[grid_wake obj.wing_segments(i).grid_wake];
                    panels_wake=[panels_wake obj.wing_segments(i).panels_wake+grid_wake_len_b4];
                end 
            end
            
            if obj.symmetric==1
                grid_mirror=[];
                grid_flat_mirror=[];
                
                grid_vol_mirror_up=[];
                grid_vol_mirror_lo=[];
                is_te_mirror=[];
                
                if (wake==1)||(wake==2)
                    grid_wake_mirror=[];
                    panels_wake_mirror=[];
                end
                
                for i=1:length(obj.wing_segments)
                    obj.wing_segments(i)=obj.wing_segments(i).right_control_surfaces();
                    obj.wing_segments(i)=obj.wing_segments(i).compute_grid(x_max,y_max,n_x_min,wake);
                    grid_mirror=[grid_mirror obj.wing_segments(i).grid];
                    grid_flat_mirror=[grid_flat_mirror obj.wing_segments(i).grid_flat];
                    
                    is_te_mirror=[is_te_mirror obj.wing_segments(i).is_te];
                    grid_vol_mirror_up=[grid_vol_mirror_up obj.wing_segments(i).grid_vol_upper];
                    grid_vol_mirror_lo=[grid_vol_mirror_lo obj.wing_segments(i).grid_vol_lower];
                    
                    if (wake==1)||(wake==2)
                        grid_wake_mirror=[grid_wake_mirror obj.wing_segments(i).grid_wake];
                    end
                end
                
                grid_len_b4=length(grid_mirror);
                obj.grid=[grid [grid_mirror(1,:);-grid_mirror(2,:);grid_mirror(3,:)]];
                obj.grid_flat=[grid_flat [grid_flat_mirror(1,:);-grid_flat_mirror(2,:);grid_flat_mirror(3,:)]];

                obj.is_te=[is_te is_te_mirror];
                obj.grid_vol_up=[grid_vol_up [grid_vol_mirror_up(1,:);-grid_vol_mirror_up(2,:);grid_vol_mirror_up(3,:)]];
                obj.grid_vol_lo=[grid_vol_lo [grid_vol_mirror_lo(1,:);-grid_vol_mirror_lo(2,:);grid_vol_mirror_lo(3,:)]];
                if (wake==1)||(wake==2)
                    grid_wake_len_b4=length(grid_wake_mirror);
                    obj.grid_wake=[grid_wake [grid_wake_mirror(1,:);-grid_wake_mirror(2,:);grid_wake_mirror(3,:)]];
                end
                
                obj.panels=[panels,[panels(2,:);panels(1,:);panels(4,:);panels(3,:)]+grid_len_b4];
                
                if (wake==1)||(wake==2)
                   obj.panels_wake=[panels_wake,[panels_wake(2,:);panels_wake(1,:);panels_wake(4,:);panels_wake(3,:)]+grid_wake_len_b4]; 
                end
                
                obj.te_idx=[te_idx,te_idx+grid_len_b4];
%                 for i=1:length(obj.wing_segments)
%                     obj.wing_segments(i)=obj.wing_segments(i).mirror_control_surfaces();
%                 end
                for i=1:length(obj.wing_segments)
                    obj.wing_segments(i)=obj.wing_segments(i).left_control_surfaces();
                end
            else
                
                obj.grid=grid;
                obj.grid_flat=grid_flat;
                obj.is_te=is_te;
                obj.grid_vol_up=grid_vol_up;
                obj.grid_vol_lo=grid_vol_lo;
                
                obj.panels=panels;
                obj.te_idx=te_idx;
                if (wake==1)||(wake==2)
                    obj.grid_wake=grid_wake;
                    obj.panels_wake=panels_wake;
                end
            end
            
            obj.grid_vol=[obj.grid_vol_up obj.grid_vol_lo];
            inv_pan(1,:)=obj.panels(4,:);
            inv_pan(2,:)=obj.panels(3,:);
            inv_pan(4,:)=obj.panels(1,:);
            inv_pan(3,:)=obj.panels(2,:);
            obj.panels_vol=[obj.panels inv_pan+length(obj.grid_vol_up)]; 
            
            obj.is_te_vol=[obj.is_te obj.is_te*0];
            obj.te_idx_vol=[obj.te_idx obj.te_idx+length(obj.grid)];
            op_idx=1:1:length(obj.is_te);
            op_idx=op_idx+length(obj.is_te);
            obj.opposite_te_vol=[op_idx obj.is_te*0];
        end
        
        function panel_to_beam_element=compute_force_interpolation_matrix(obj,panel_to_beam_element)
            %% clear data in panel_to_beam_element 
            panel_to_beam_element(obj.panel_start_idx:obj.panel_start_idx+size(obj.panels,2)-1,:)=panel_to_beam_element(obj.panel_start_idx:obj.panel_start_idx+size(obj.panels,2)-1,:)*0;
            %% str grid points
            if obj.symmetric
                aerGrid=obj.grid(:,1:end/2);
                aerPan=obj.panels(:,1:end/2);
            else
                aerGrid=obj.grid;
                aerPan=obj.panels;
            end
            strGrid=(obj.wingbox_coords(:,:,1)+obj.wingbox_coords(:,:,2))/2;
%             figure(1);
%             hold on; scatter3(strGrid(1,:),strGrid(2,:),strGrid(3,:),'ro','filled')

            %% loop over elements
            for iEl=1:size(strGrid,2)-1
                if iEl==1
                    distIB=1*ones(1,size(aerGrid,2));
                    pIB=[];
                    nIB=[];
                else
                    pIB=pOB;
                    nIB=nOB;
                    distIB=distOB;
                end
                if iEl==length(strGrid)-1
                    distOB=-1*ones(1,size(aerGrid,2));
                else
                    if iEl==1
                        %compute plane normal vector Outboard
                        direction=(((strGrid(:,iEl+1)-strGrid(:,iEl))+(strGrid(:,iEl+2)-strGrid(:,iEl+1)))/2);
                    else
                        %compute plane normal vector Outboard
                        direction=(((strGrid(:,iEl)-strGrid(:,iEl-1))+(strGrid(:,iEl+1)-strGrid(:,iEl))+(strGrid(:,iEl+2)-strGrid(:,iEl+1)))/2);
                    end
                    nOB=direction/norm(direction);
                    pOB=strGrid(:,iEl+1);
                    %claculate distance of aerGrid to plane between beam element iEl and beam element iEl+1
                    pqOB=aerGrid-pOB;
                    distOB=dot(pqOB,repmat(nOB,1,size(pqOB,2)));
                    % for all beam elements 
                end
                insideFlag=and(distOB<=0,distIB>0);

                %% check plot
%                 figure;
%                 scatter3(obj.grid(1,1:end/2),obj.grid(2,1:end/2),obj.grid(3,1:end/2),'o')
%                 hold on;
%                 scatter3(aerGrid(1,insideFlag),aerGrid(2,insideFlag),aerGrid(3,insideFlag),'+')
                %% loop over all panels and fill table
                b1=strGrid(:,iEl);
                b2=strGrid(:,iEl+1);
                bAx = b2 - b1;
                for iPan=1:size(aerPan,2)
                    if ~any(insideFlag(aerPan(:,iPan)))
                        % special case of small beam elements and large aero panels:
                        % panels which are both, not completely outboard nor completely inboard
                        % of this beam element need to be split as well; 
                        if and(all(distIB(aerPan([1,4],iPan))<0), all(distOB(aerPan([2 3],iPan))>0))
                            insideFlag(aerPan(1,iPan))=1; %set insideFlag active for one vertex so that it will be split in the next step
                        end
                    end

                    if any(insideFlag(aerPan(:,iPan))) %panel iPan fully or partially belonging to this iEl
                        % calc normal distance from 1/4 point of panel to beam element
                        % axis
                        % quarter point of panel
                        qp=.75*mean(aerGrid(:,aerPan(1:2,iPan)),2)+.25*mean(aerGrid(:,aerPan(3:4,iPan)),2);
                        % point of intersection of a plane normal to the beam axis
                        % through qp
                        s=b2 + ((bAx / norm(bAx)) * (qp - (b2))') * (bAx / norm(bAx));
                        % distance of qp  from s
                        dist=qp-s;
                        if all(insideFlag(aerPan(:,iPan))) %fully belonging to this iEl
                            % fill ptb matrix for fully contained panel
                            panel_to_beam_element(obj.panel_start_idx-1+iPan,[1:5])=[iEl 1 dist'];
                        else %partially belonging to this iEl
                            % estimate part of area of panel positive of IB
                            if iEl==1
                                pAreaIB=1;
                            else
                                p1distIB=distIB(aerPan(1,iPan));
                                p2distIB=distIB(aerPan(2,iPan));
                                p3distIB=distIB(aerPan(3,iPan));
                                p4distIB=distIB(aerPan(4,iPan));
                                if and(p1distIB<=0,p2distIB<=0)
                                    pLeIB=0;
                                elseif and(p1distIB>=0,p2distIB>=0)
                                    pLeIB=1;
                                else
                                    pLeIB=p2distIB/(abs(p1distIB)+abs(p2distIB));
                                end
                                if and(p4distIB<=0,p3distIB<=0)
                                    pTeIB=0;
                                elseif and(p4distIB>=0,p3distIB>=0)
                                    pTeIB=1;
                                else
                                    pTeIB=p3distIB/(abs(p4distIB)+abs(p3distIB));
                                end
                                pAreaIB=0.5*pTeIB+0.5*pLeIB;
                            end
                            % estimate part of area of panel positive of OB
                            if iEl==size(strGrid,2)-1
                                pAreaOB=0;
                            else
                                p1distOB=distOB(aerPan(1,iPan));
                                p2distOB=distOB(aerPan(2,iPan));
                                p3distOB=distOB(aerPan(3,iPan));
                                p4distOB=distOB(aerPan(4,iPan));
                                if and(p1distOB<=0,p2distOB<=0)
                                    pLeOB=0;
                                elseif and(p1distOB>=0,p2distOB>=0)
                                    pLeOB=1;
                                else
                                    pLeOB=p2distOB/(abs(p1distOB)+abs(p2distOB));
                                end
                                if and(p4distOB<=0,p3distOB<=0)
                                    pTeOB=0;
                                elseif and(p4distOB>=0,p3distOB>=0)
                                    pTeOB=1;
                                else
                                    pTeOB=p3distOB/(abs(p4distOB)+abs(p3distOB));
                                end
                                pAreaOB=0.5*pTeOB+0.5*pLeOB;
                            end
                            % get part inbetween
                            pArea=pAreaIB-pAreaOB;
                            % fill ptb matrix for panel partially belonging to this iEl
                            if all(panel_to_beam_element(obj.panel_start_idx-1+iPan,[1:5])==0)
                                %panel is first associated to this beam element
                                panel_to_beam_element(obj.panel_start_idx-1+iPan,[1:5])=[iEl pArea dist'];
                            elseif all(panel_to_beam_element(obj.panel_start_idx-1+iPan,[6:10])==0)
                                %panel is already associated to other beam element
                                panel_to_beam_element(obj.panel_start_idx-1+iPan,[6:10])=[iEl pArea dist'];
                            elseif all(panel_to_beam_element(obj.panel_start_idx-1+iPan,[11:15])==0)
                                %panel is already associated to two other beam elements
                                panel_to_beam_element(obj.panel_start_idx-1+iPan,[11:15])=[iEl pArea dist'];
                            else
                                %panel is associated to three other beam elements
                                disp('warning: panel can not be associated to more than three beam elements')
                            end
                        end
                    end
                end
            end
            nBeamEl=size(strGrid,2)-1;
            if obj.symmetric==1
                offset=length(obj.panels)/2;
                for j=1:size(obj.panels,2)/2
                    panel_to_beam_element(j+offset+obj.panel_start_idx-1,:)=panel_to_beam_element(j+obj.panel_start_idx-1,:);
                    panel_to_beam_element(j+offset+obj.panel_start_idx-1,[4 9 14])=-panel_to_beam_element(j+offset+obj.panel_start_idx-1,[4 9 14]);
                    for i=1:5:size(panel_to_beam_element,2)
                        if not(panel_to_beam_element(j+offset+obj.panel_start_idx-1,i)==0)
                            panel_to_beam_element(j+offset+obj.panel_start_idx-1,i)=-panel_to_beam_element(j+offset+obj.panel_start_idx-1,i)+nBeamEl+1;
%                             panel_to_beam_element(j+offset+obj.panel_start_idx-1,i+2)=-panel_to_beam_element(j+offset+obj.panel_start_idx-1,i+2);
%                             panel_to_beam_element(j+obj.panel_start_idx-1,i+2)=-panel_to_beam_element(j+obj.panel_start_idx-1,i+2);
                        end
                    end
                end   
                for j=1:size(obj.panels,2)/2
                    for i=1:5:size(panel_to_beam_element,2)
                        if not(panel_to_beam_element(j+obj.panel_start_idx-1,i)==0)
                            com_idx=panel_to_beam_element(j+obj.panel_start_idx-1,i)+nBeamEl;%+(n_beam-2*(panel_to_beam_element(j+obj.panel_start_idx-1,i)-1));
                            panel_to_beam_element(j+obj.panel_start_idx-1,i)=com_idx; 
                        end
                    end
                end
            end 
            
        end
        
        function obj=T_matrix(obj,panel_to_beam_element,grid,panels,beam)
            
            if obj.symmetric==1
                f=2;
            else
                f=1;
            end
           %  obj.T=zeros(length(obj.wingbox_coords(1,:,1))*6*f-6,size(panel_to_beam_element,1));
           %  obj.db=zeros(length(obj.wingbox_coords(1,:,1))*6*f-6,size(panel_to_beam_element,1));
            obj.T=zeros(length(beam.beamelement)*6+6,size(panel_to_beam_element,1));
            obj.db=zeros(length(beam.beamelement)*6+6,size(panel_to_beam_element,1));
            for j=obj.panel_start_idx:1:(obj.panel_start_idx+size(obj.panels,2))-1
                for i=1:5:size(panel_to_beam_element,2)
                    if not(panel_to_beam_element(j,i)==0)
                        diff=norm(panel_to_beam_element(j,i+2:i+4));
                        %                       idx=j-obj.panel_start_idx+1;
                        %                       r=(0.75*obj.grid(:,obj.panels(1,idx))+0.25*obj.grid(:,obj.panels(4,idx)))-(0.75*obj.grid(:,obj.panels(2,idx))+0.25*obj.grid(:,obj.panels(3,idx)));
                        %phi=180/pi*beam.beamelement(panel_to_beam_element(j,i)).phi;
                        r=(0.75*grid(:,panels(1,j))+0.25*grid(:,panels(4,j)))-(0.75*grid(:,panels(2,j))+0.25*grid(:,panels(3,j)));
                        
                        s=abs(r(2));
                        Lz=eye(3,3);%beam.beamelement(panel_to_beam_element(j,i)).T(1:3,1:3);
                        
                        F=[0;
                           0;
                           -1/2*panel_to_beam_element(j,i+1)]*s;
                       
                        Ft=Lz'*F;
                        
                        
                        Lz=beam.beamelement(panel_to_beam_element(j,i)).T(1:3,1:3);
                        M=[-1/24*panel_to_beam_element(j,i+1)*s^2;
                            -diff*1/2*panel_to_beam_element(j,i+1)*s;
                            0];

                        Mt1=Lz'*M;
                        
                        M=[+1/24*panel_to_beam_element(j,i+1)*s^2;
                            -diff*1/2*panel_to_beam_element(j,i+1)*s;
                            0];
                        
                        Mt2=Lz'*M;

                        %beamlength=norm(obj.wingbox_coords(:,panel_to_beam_element(j,i)+1)-obj.wingbox_coords(:,panel_to_beam_element(j,i)));
                        obj.T((panel_to_beam_element(j,i)-1)*6+1,j)=obj.T((panel_to_beam_element(j,i)-1)*6+1,j)+Ft(1);
                        obj.T((panel_to_beam_element(j,i)-1)*6+2,j)=obj.T((panel_to_beam_element(j,i)-1)*6+2,j)+Ft(2);
                        obj.T((panel_to_beam_element(j,i)-1)*6+3,j)=obj.T((panel_to_beam_element(j,i)-1)*6+3,j)+Ft(3);
                        obj.T((panel_to_beam_element(j,i)-1)*6+4,j)=obj.T((panel_to_beam_element(j,i)-1)*6+4,j)+Mt1(1);
                        obj.T((panel_to_beam_element(j,i)-1)*6+5,j)=obj.T((panel_to_beam_element(j,i)-1)*6+5,j)+Mt1(2);
                        obj.T((panel_to_beam_element(j,i)-1)*6+6,j)=obj.T((panel_to_beam_element(j,i)-1)*6+6,j)+Mt1(3);    

                        obj.T((panel_to_beam_element(j,i))*6+1,j)=obj.T((panel_to_beam_element(j,i))*6+1,j)+Ft(1);
                        obj.T((panel_to_beam_element(j,i))*6+2,j)=obj.T((panel_to_beam_element(j,i))*6+2,j)+Ft(2);
                        obj.T((panel_to_beam_element(j,i))*6+3,j)=obj.T((panel_to_beam_element(j,i))*6+3,j)+Ft(3);
                        obj.T((panel_to_beam_element(j,i))*6+4,j)=obj.T((panel_to_beam_element(j,i))*6+4,j)+Mt2(1);
                        obj.T((panel_to_beam_element(j,i))*6+5,j)=obj.T((panel_to_beam_element(j,i))*6+5,j)+Mt2(2);
                        obj.T((panel_to_beam_element(j,i))*6+6,j)=obj.T((panel_to_beam_element(j,i))*6+6,j)+Mt2(3);
                        
                        obj.db((panel_to_beam_element(j,i)-1)*6+1,j)=0;
                        obj.db((panel_to_beam_element(j,i)-1)*6+2,j)=0;
                        obj.db((panel_to_beam_element(j,i)-1)*6+3,j)=0;

                        dind=Lz'*[0  -1/2*panel_to_beam_element(j,i+1)*4*pi  0]';
                        %  obj.db((panel_to_beam_element(j,i)-1)*6+4,j)=obj.db((panel_to_beam_element(j,i)-1)*6+4,j)-1/7*panel_to_beam_element(j,i+1)*4*pi*sind(phi)*cosd(phi);
                        %  obj.db((panel_to_beam_element(j,i)-1)*6+5,j)=obj.db((panel_to_beam_element(j,i)-1)*6+5,j)+1/2*panel_to_beam_element(j,i+1)*4*pi*cosd(phi)*cosd(phi)*cosd(dihed);
                        obj.db((panel_to_beam_element(j,i)-1)*6+4,j)=obj.db((panel_to_beam_element(j,i)-1)*6+4,j)+dind(1);%-1/2*panel_to_beam_element(j,i+1)*4*pi*sind(phi);
                        obj.db((panel_to_beam_element(j,i)-1)*6+5,j)=obj.db((panel_to_beam_element(j,i)-1)*6+5,j)+dind(2);%-1/4*panel_to_beam_element(j,i+1)*4*pi*cosd(phi);
                        obj.db((panel_to_beam_element(j,i)-1)*6+6,j)=obj.db((panel_to_beam_element(j,i)-1)*6+6,j)+dind(3);%-1/2*panel_to_beam_element(j,i+1)*4*pi*cosd(phi)*sin(dihed);
                        
                        obj.db((panel_to_beam_element(j,i))*6+1,j)=0;
                        obj.db((panel_to_beam_element(j,i))*6+2,j)=0;
                        obj.db((panel_to_beam_element(j,i))*6+3,j)=0;
                        %  obj.db((panel_to_beam_element(j,i))*6+4,j)=obj.db((panel_to_beam_element(j,i))*6+4,j)-1/7*panel_to_beam_element(j,i+1)*4*pi*sind(phi)*cosd(phi);
                        %  obj.db((panel_to_beam_element(j,i))*6+5,j)=obj.db((panel_to_beam_element(j,i))*6+5,j)+1/2*panel_to_beam_element(j,i+1)*4*pi*cosd(phi)*cosd(phi)*cosd(dihed);
                        obj.db((panel_to_beam_element(j,i))*6+4,j)=obj.db((panel_to_beam_element(j,i))*6+4,j)+dind(1);%-1/2*panel_to_beam_element(j,i+1)*4*pi*sind(phi);
                        obj.db((panel_to_beam_element(j,i))*6+5,j)=obj.db((panel_to_beam_element(j,i))*6+5,j)+dind(2);%-1/4*panel_to_beam_element(j,i+1)*4*pi*cosd(phi);
                        obj.db((panel_to_beam_element(j,i))*6+6,j)=obj.db((panel_to_beam_element(j,i))*6+6,j)+dind(3);%-1/2*panel_to_beam_element(j,i+1)*4*pi*cosd(phi)*sind(dihed);
                    end
                end
            end
        end
        
        function panel_to_beam_element=initialize_beam_to_panel_for_block(obj,beam_idx,k,block,panels,grid,panel_to_beam_element,mirr,comp_idx)
            observer=0;
            plotBtp=0; % debug plot by simon 
            if mirr==1
                grid_offset=length(obj.panels)/2;
            else
                grid_offset=0;
            end
           
            for i=1:obj.wing_segments(k).n_span
                
                n_pan=sum([obj.wing_segments(k).n_le_panels obj.wing_segments(k).n_ctr1_panels obj.wing_segments(k).n_sp_panels obj.wing_segments(k).n_ctr2_panels  obj.wing_segments(k).n_te_panels]);

                for j=1:n_pan
                    
                    idx=obj.wing_segments(k).panel_start_idx+n_pan*(i-1)-1+j+grid_offset;
                    % check if current panel lies within block
                    % if ymin of panel is smaller or equal to ymax of block  OR  ymax of panel is smaller or equal to min 
                    if ~(min(grid(comp_idx,panels(:,idx)))> 100*eps+max(block(comp_idx,:))) || ~(max(grid(comp_idx,panels(:,idx)))> min(block(comp_idx,:)-100*eps))
                        
                        if size(block,2)==3
                            in=points_in_triangle(block,grid(:,panels(:,idx)));
                        else
                            in=points_in_quad(block,grid(:,panels(:,idx))); %which points of current panel lie in block
                        end
                        
                        if in(1)==1 && in(2)==1 && in(3)==1 && in(4)==1     %if panel lies completely in block
                            if plotBtp==1
                                observer=0;
                                figure;
                                h=fill3(block(1,:),block(2,:),block(3,:),'b');
                                set(h,'facealpha',.25);
                                hold on;
                                axis equal;
                                plot3(obj.wingbox_coords(1,[beam_idx, beam_idx+1],1),obj.wingbox_coords(2,[beam_idx, beam_idx+1],1),obj.wingbox_coords(3,[beam_idx, beam_idx+1],1),'--ks')
                                plot3(obj.wingbox_coords(1,[beam_idx, beam_idx+1],2),obj.wingbox_coords(2,[beam_idx, beam_idx+1],2),obj.wingbox_coords(3,[beam_idx, beam_idx+1],2),'--ks')
                                plot3((obj.wingbox_coords(1,[beam_idx, beam_idx+1],1)+obj.wingbox_coords(1,[beam_idx, beam_idx+1],2))/2,(obj.wingbox_coords(2,[beam_idx, beam_idx+1],1)+obj.wingbox_coords(2,[beam_idx, beam_idx+1],2))/2,(obj.wingbox_coords(3,[beam_idx, beam_idx+1],1)+obj.wingbox_coords(3,[beam_idx, beam_idx+1],2))/2,'-ks')
                                h=fill3(grid(1,panels(:,idx)),grid(2,panels(:,idx)),grid(3,panels(:,idx)),'g');
                                set(h,'facealpha',.25);
                            end
                            
                            
                            rb=0.25*obj.wingbox_coords(:,beam_idx+1,1)+0.25*obj.wingbox_coords(:,beam_idx+1,2)+0.25*obj.wingbox_coords(:,beam_idx,1)+0.25*obj.wingbox_coords(:,beam_idx,2);
                            r=0.5*(0.75*grid(:,panels(1,idx))+0.25*grid(:,panels(4,idx)))+0.5*(0.75*grid(:,panels(2,idx))+0.25*grid(:,panels(3,idx)));
                            vecb=(0.5*obj.wingbox_coords(:,beam_idx+1,1)+0.5*obj.wingbox_coords(:,beam_idx+1,2))-(0.5*obj.wingbox_coords(:,beam_idx,1)+0.5*obj.wingbox_coords(:,beam_idx,2));
                            theta=cross(vecb,r-rb)/(norm(r-rb)*norm(vecb));
                            s=      (0.5*obj.wingbox_coords(:,beam_idx+1,1)+0.5*obj.wingbox_coords(:,beam_idx+1,2))  + ((vecb / norm(vecb)) * (r - (0.5*obj.wingbox_coords(:,beam_idx+1,1)+0.5*obj.wingbox_coords(:,beam_idx+1,2)))'       ) * (vecb / norm(vecb));
%                                     if theta(3)<0
%                                         dist=-norm(r-rb);
%                                     else
%                                         dist=norm(r-rb);
%                                     end
                            dist=r-s;
                            if observer==1
                                fill3(grid(1,panels(:,idx)),grid(2,panels(:,idx)),grid(3,panels(:,idx)),[0.2 0.7 0.4])
                            end
                            panel_to_beam_element(idx,1)=beam_idx;
                            panel_to_beam_element(idx,2)=1.000;
                            panel_to_beam_element(idx,3:5)=dist;
                        else   %panel needs to be cut and only portion within block goes to panel to beam element matrix
                                [newPoints,faceNewPoint]=cut_polygons(block,grid(:,panels(:,idx)));
                            if sum(faceNewPoint)>=1 
                                
                                %plot3(new_polygon(1,:),new_polygon(2,:),new_polygon(3,:),'bo');
                                if size(block,2)==3
                                    idOrigPanelPointsInBlock=points_in_triangle_poly(block,grid(:,panels(:,idx)));
                                else
                                    idOrigPanelPointsInBlock=points_in_quad_poly(block,grid(:,panels(:,idx))); 
                                end

                                %red_polygon is the polygon inside the
                                %block consisting of newPoints and
                                %origPanelPointsInBlock
                                nRedPolygonPoints=size(newPoints,2)+sum(idOrigPanelPointsInBlock);
                                red_polygon=zeros(3,nRedPolygonPoints);
                                iRedPolygonPoints=1;
                                for iFace=1:4  %for every face of orig Panel
                                    if idOrigPanelPointsInBlock(iFace)==1
                                        red_polygon(:,iRedPolygonPoints)=grid(:,panels(iFace,idx));
                                        iRedPolygonPoints=iRedPolygonPoints+1;
                                    end
                                    for jNewPoints=1:length(faceNewPoint)
                                        if iFace==faceNewPoint(jNewPoints)
                                            red_polygon(:,iRedPolygonPoints)=newPoints(:,jNewPoints);
                                            iRedPolygonPoints=iRedPolygonPoints+1;
                                        end
                                    end
                                end
                                %fill3(red_polygon(1,:),red_polygon(2,:),red_polygon(3,:),0.1)
                                if size(red_polygon,1)>=3 && size(red_polygon,2)>=3
                                    if observer==1
                                        handle=fill3(red_polygon(1,:),red_polygon(2,:),red_polygon(3,:),[0.2 0.7 0.4]);
                                        alpha(handle,0.4);
                                    end
                                    if plotBtp==1
                                        observer=0;
                                        figure;
                                        h=fill3(block(1,:),block(2,:),block(3,:),'b');
                                        set(h,'facealpha',.25);
                                        hold on;
                                        axis equal;
                                        plot3(obj.wingbox_coords(1,[beam_idx, beam_idx+1],1),obj.wingbox_coords(2,[beam_idx, beam_idx+1],1),obj.wingbox_coords(3,[beam_idx, beam_idx+1],1),'--ks')
                                        plot3(obj.wingbox_coords(1,[beam_idx, beam_idx+1],2),obj.wingbox_coords(2,[beam_idx, beam_idx+1],2),obj.wingbox_coords(3,[beam_idx, beam_idx+1],2),'--ks')
                                        plot3((obj.wingbox_coords(1,[beam_idx, beam_idx+1],1)+obj.wingbox_coords(1,[beam_idx, beam_idx+1],2))/2,(obj.wingbox_coords(2,[beam_idx, beam_idx+1],1)+obj.wingbox_coords(2,[beam_idx, beam_idx+1],2))/2,(obj.wingbox_coords(3,[beam_idx, beam_idx+1],1)+obj.wingbox_coords(3,[beam_idx, beam_idx+1],2))/2,'-ks')
                                        handle=fill3(red_polygon(1,:),red_polygon(2,:),red_polygon(3,:),'y');
                                        scatter3(grid(1,panels(:,idx)),grid(2,panels(:,idx)),grid(3,panels(:,idx)),'dr')
                                        alpha(handle,0.25);
                                    end
                                        
                                    panel_area=1/2*norm(cross(grid(:,panels(3,idx))-grid(:,panels(1,idx)),grid(:,panels(4,idx))-grid(:,panels(2,idx))));
                                    moment_line=[(0.75*grid(:,panels(1,idx))+0.25*grid(:,panels(4,idx))) (0.75*grid(:,panels(2,idx))+0.25*grid(:,panels(3,idx)))];
                                    red_moment_line=cut_line_polygon(red_polygon,moment_line);
                                    %plot3(red_moment_line(1,:),red_moment_line(2,:),red_moment_line(3,:),'-ko');
                                
                                    if isnan(panel_area)
                                        panel_area;
                                    end
                                    if isnan(abs(planar_polygon_area(red_polygon)))
                                        red_polygon;
                                        abs(planar_polygon_area(red_polygon));
                                    end
                                    
                                    rb=0.25*obj.wingbox_coords(:,beam_idx+1,1)+0.25*obj.wingbox_coords(:,beam_idx+1,2)+0.25*obj.wingbox_coords(:,beam_idx,1)+0.25*obj.wingbox_coords(:,beam_idx,2);
                                    %r= sum(red_polygon,2)/size(red_polygon,2);
                                    r=0.5*red_moment_line(:,1)+0.5*red_moment_line(:,end); 
                                    vecb=(0.5*obj.wingbox_coords(:,beam_idx+1,1)+0.5*obj.wingbox_coords(:,beam_idx+1,2))-(0.5*obj.wingbox_coords(:,beam_idx,1)+0.5*obj.wingbox_coords(:,beam_idx,2));
                                    theta=cross(vecb,r-rb)/(norm(r-rb)*norm(vecb));
                                    s=(0.5*obj.wingbox_coords(:,beam_idx+1,1)+0.5*obj.wingbox_coords(:,beam_idx+1,2)) + ((vecb / norm(vecb)) * (r - (0.5*obj.wingbox_coords(:,beam_idx+1,1)+0.5*obj.wingbox_coords(:,beam_idx+1,2)))') * (vecb / norm(vecb));
%                                     if theta(3)<0
%                                         dist=-norm(r-rb);
%                                     else
%                                         dist=norm(r-rb);
%                                     end
                                    dist=r-s;
                                    
                                    if plotBtp==1
                                        plot3(moment_line(1,:),moment_line(2,:),moment_line(3,:),'-k');
                                        plot3(red_moment_line(1,:),red_moment_line(2,:),red_moment_line(3,:),'-ro');
                                        scatter3(rb(1),rb(2),rb(3),'r','filled');
                                        scatter3(r(1),r(2),r(3),'rd','filled');
                                        scatter3(s(1),s(2),s(3),'bd','filled');
                                        line([r(1), s(1)],[r(2), s(2)],[r(3), s(3)])
                                    end
                                    if panel_to_beam_element(idx,2)==0
                                        panel_to_beam_element(idx,1)=beam_idx;
                                        if abs(planar_polygon_area(red_polygon))/panel_area<1
                                            panel_to_beam_element(idx,2)=abs(planar_polygon_area(red_polygon))/panel_area;
                                        else
                                            panel_to_beam_element(idx,2)=1;
                                        end
                                        if isnan(panel_to_beam_element(idx,2))
                                            panel_to_beam_element(idx,2)=0;
                                        end
                                        panel_to_beam_element(idx,3:5)=dist;
                                        
                                    elseif panel_to_beam_element(idx,7)==0
                                        panel_to_beam_element(idx,6)=beam_idx;
                                        if abs(planar_polygon_area(red_polygon))/panel_area<1
                                            panel_to_beam_element(idx,7)=abs(planar_polygon_area(red_polygon))/panel_area;
                                        else
                                           panel_to_beam_element(idx,7)=1; 
                                        end
                                        panel_to_beam_element(idx,8:10)=dist;
                                        
                                    else
                                        panel_to_beam_element(idx,12);
                                        
                                        panel_to_beam_element(idx,11)=beam_idx;
                                        share=abs(planar_polygon_area(red_polygon));
                                        if ~(isinf(share)) && ~(isnan(share))
                                            panel_to_beam_element(idx,12)=share/panel_area;
                                            panel_to_beam_element(idx,13:15)=dist;
                                        else
                                            panel_to_beam_element(idx,12)=0;
                                            panel_to_beam_element(idx,7)=0;
                                        end
                                    end
                                else
                                    if panel_to_beam_element(idx,1)==0
                                        panel_to_beam_element(idx,1)=beam_idx;
                                        panel_to_beam_element(idx,2)=0.000;
                                        panel_to_beam_element(idx,3:5)=0.000;
                                    end
                                end
                                
                            end
                        end
                    end
                end
            end
        end
        
        function obj=compute_beam_forces(obj,panel_to_beam_element,panel_forces,panel_moments, beam)
            if obj.symmetric==1
                f=2;
            else
                f=1;
            end
            k=obj.panel_start_idx;
            obj.beam_forces_structmesh=zeros(3,size(obj.wingbox_coords,2)*f-1);
            obj.beam_moments_structmesh=zeros(3,size(obj.wingbox_coords,2)*f-1);
            end_idx=(obj.panel_start_idx+size(obj.panels,2)-1);
            ctr=1;
            overall_force=0;
            overall_force_loc=0;
            for j=obj.panel_start_idx:1:end_idx                             %for every panel i
                for i=1:5:size(panel_to_beam_element,2)-1                   %for all 3 corresponding beamelements j of the panel 
                    if not(panel_to_beam_element(j,i)==0)                   %only if there is parts of the panels assigned to the current beamelement j
                        beamlength=beam.beamelement(panel_to_beam_element(k,i)).le;%norm(obj.wingbox_coords(:,panel_to_beam_element(k,i)+1)-obj.wingbox_coords(:,panel_to_beam_element(k,i)));
                        %Trot=beam.beamelement(panel_to_beam_element(j,i)).f_rotM_6dof(real(beam.beamelement(panel_to_beam_element(j,i)).nu),real(beam.beamelement(panel_to_beam_element(j,i)).epsilon),0);
                        if obj.symmetric==1
                            tf=panel_forces(:,j);
                            %tf(2)=-tf(2);                                   %mirror global y component
                            force_loc=beam.beamelement(panel_to_beam_element(j,i)).T(1:3,1:3)*tf;%[0 0 dot(crp,tf)]';
                            moment_loc=beam.beamelement(panel_to_beam_element(j,i)).T(1:3,1:3)*panel_moments(:,j);
                            %beam.beamelement(panel_to_beam_element(j,i)).T(1:3,1:3)
                            force_glob=tf;
                            % force_loc=beam.beamelement(panel_to_beam_element(j,i)).T(1:3,1:3)*panel_forces(:,j);
                            % force_loc(2)=-force_loc(2);
                        else
                            
                           tf=panel_forces(:,j);
                          %  tf(2)=tf(2);
                           force_loc=beam.beamelement(panel_to_beam_element(j,i)).T(1:3,1:3)*tf;%[0 0 dot(crp,tf)]';
                            moment_loc=beam.beamelement(panel_to_beam_element(j,i)).T(1:3,1:3)*panel_moments(:,j);
                            %beam.beamelement(panel_to_beam_element(j,i)).T(1:3,1:3)
                            force_glob=tf;
                            % force_loc=panel_forces(:,j);%Trot(1:3,1:3)*
                            %  force_loc=beam.beamelement(panel_to_beam_element(j,i)).T(1:3,1:3)*panel_forces(:,j);
                          %   force_loc=panel_forces(:,j);%[0 0 dot(crp,panel_forces(:,j))]';
%                             
%                             
                         %    force_glob=panel_forces(:,j);
                        end
                        overall_force=overall_force+force_glob(3)*panel_to_beam_element(j,i+1);
                        overall_force_loc=overall_force_loc+force_loc(3)*panel_to_beam_element(j,i+1);
                        obj.beam_forces_structmesh(:,panel_to_beam_element(j,i))=obj.beam_forces_structmesh(:,panel_to_beam_element(j,i))+force_loc*panel_to_beam_element(j,i+1)/beamlength;
                        %obj.beam_moments_structmesh(2,panel_to_beam_element(j,i))=obj.beam_moments_structmesh(2,panel_to_beam_element(j,i))+sign(force_loc(3))*norm(force_loc(1:2:3))*panel_to_beam_element(j,i+1)*norm(panel_to_beam_element(j,i+2:i+4))/beamlength;
                        % transform distance to local coordinate system
                        dist = beam.beamelement(panel_to_beam_element(j,i)).T(1:3,1:3)*panel_to_beam_element(j,i+2:i+4)';
                        obj.beam_moments_structmesh(2,panel_to_beam_element(j,i))=obj.beam_moments_structmesh(2,panel_to_beam_element(j,i))+panel_to_beam_element(j,i+1)*(-force_loc(3)*dist(1)+force_loc(1)*dist(3)+moment_loc(2))/beamlength;
                    end
                end
                ctr=ctr+1;
                k=k+1;
                if k==(size(obj.panels,2)/2)+1
                   k=obj.panel_start_idx; 
                end
            end
%             overall_force
%             overall_force_loc
%             if obj.symmetric==1
%                 mm=obj.beam_forces_structmesh;
%                 obj.beam_forces_structmesh=mm(:,end-1:-1:1);
%                 mm=obj.beam_moments_structmesh;
%                 obj.beam_moments_structmesh=mm(:,end-1:-1:1);   
%             end          
        end
        
        function obj=compute_c4_forces(obj,panel_forces,panels,grid)
            
            obj.c4_forces=[];
            obj.c4_moments=[];
            
            
            for i=1:length(obj.wing_segments)
                obj.wing_segments(i)=obj.wing_segments(i).compute_forces(panel_forces,panels,grid);
                obj.c4_forces=[obj.c4_forces obj.wing_segments(i).span_forces];
                obj.c4_moments=[obj.c4_moments obj.wing_segments(i).span_moments];
            end
            
            obj.c4_forces_structmesh=zeros(3,length(obj.wingbox_rl_coords_mid));
            obj.c4_moments_structmesh=zeros(3,length(obj.wingbox_rl_coords_mid));
            obj.c4_forces_structmesh(1,:)=interp1(obj.c4_rl_coords,obj.c4_forces(1,:), obj.wingbox_rl_coords_mid,'linear','extrap');
            obj.c4_forces_structmesh(2,:)=interp1(obj.c4_rl_coords,obj.c4_forces(2,:), obj.wingbox_rl_coords_mid,'linear','extrap');
            obj.c4_forces_structmesh(3,:)=interp1(obj.c4_rl_coords,obj.c4_forces(3,:), obj.wingbox_rl_coords_mid,'linear','extrap');
            obj.c4_moments_structmesh(1,:)=interp1(obj.c4_rl_coords,obj.c4_moments(1,:),obj.wingbox_rl_coords_mid,'linear','extrap');
            obj.c4_moments_structmesh(2,:)=interp1(obj.c4_rl_coords,obj.c4_moments(2,:), obj.wingbox_rl_coords_mid,'linear','extrap');
            obj.c4_moments_structmesh(3,:)=interp1(obj.c4_rl_coords,obj.c4_moments(3,:), obj.wingbox_rl_coords_mid,'linear','extrap');
            
            %             figure
            %             hold on
            %             plot(c4_rl_coord,obj.c4_forces(3,:));
            %             plot(c4_rl_coord_sm,obj.c4_forces_structmesh(3,:));
        end
        
        function grid=compute_deflected_grid(obj,panels,grid,deflections_structmesh)
            % computation of  station of panels along the c4 line 
            c4_rl_coords_edge=0;
            for i=1:length(obj.wing_segments)
                x=0:1/obj.wing_segments(i).n_span:1;
                %seg_rl : stations times lenght of the c4 line of the segment
                seg_rl=obj.wing_segments(i).b*x/cosd(obj.wing_segments(i).c4_sweep);
                c4_rl_coords_edge=[c4_rl_coords_edge seg_rl(2:end)+c4_rl_coords_edge(end)];
            end
            
            for j=1:6
                def(j,:)=squeeze(deflections_structmesh(j:6:end));
                %leftdef(j,:)=interp1(squeeze(obj.mesh_struct(i).y),squeeze(obj.deflect(i).def(j:6:end)),obj.AC_geo.LiftingSurfaceList(i).Panels.s_R,'spline','extrap');
            end
            %             figure
            %             plot(obj.wingbox_rl_coords)
            %             figure
            %             plot(c4_rl_coords_edge)
            %deflections_aeromesh are the deflections at the c/4 line. for
            %determination, the deflections are interpolated onto the
            %aero stations on the c/4 line
            if obj.symmetric==1
                n_def_hw=(length(def(1,:))-1)/2; %halfwing number deflections structmesh
                %left wing
                deflections_aeromesh_left(1,:)=interp1(obj.wingbox_rl_coords, def(1,(n_def_hw+1):-1:1), c4_rl_coords_edge,'lin','extrap');
                deflections_aeromesh_left(2,:)=interp1(obj.wingbox_rl_coords, def(2,(n_def_hw+1):-1:1), c4_rl_coords_edge,'lin','extrap');
                deflections_aeromesh_left(3,:)=interp1(obj.wingbox_rl_coords, def(3,(n_def_hw+1):-1:1), c4_rl_coords_edge,'lin','extrap');
                deflections_aeromesh_left(4,:)=interp1(obj.wingbox_rl_coords, def(4,(n_def_hw+1):-1:1), c4_rl_coords_edge,'lin','extrap');
                deflections_aeromesh_left(5,:)=interp1(obj.wingbox_rl_coords,def(5,(n_def_hw+1):-1:1), c4_rl_coords_edge,'lin','extrap');
                deflections_aeromesh_left(6,:)=interp1(obj.wingbox_rl_coords, def(6,(n_def_hw+1):-1:1), c4_rl_coords_edge,'lin','extrap');
                %right wing
                deflections_aeromesh(1,:)=interp1(obj.wingbox_rl_coords,[ def(1,(n_def_hw+1):1:end)], c4_rl_coords_edge,'lin','extrap');
                deflections_aeromesh(2,:)=interp1(obj.wingbox_rl_coords,[ def(2,(n_def_hw+1):1:end)], c4_rl_coords_edge,'lin','extrap');
                deflections_aeromesh(3,:)=interp1(obj.wingbox_rl_coords,[ def(3,(n_def_hw+1):1:end)], c4_rl_coords_edge,'lin','extrap');
                deflections_aeromesh(4,:)=interp1(obj.wingbox_rl_coords,[ def(4,(n_def_hw+1):1:end)], c4_rl_coords_edge,'lin','extrap');
                deflections_aeromesh(5,:)=interp1(obj.wingbox_rl_coords,[ def(5,(n_def_hw+1):1:end)], c4_rl_coords_edge,'lin','extrap');
                deflections_aeromesh(6,:)=interp1(obj.wingbox_rl_coords,[ def(6,(n_def_hw+1):1:end)], c4_rl_coords_edge,'lin','extrap');

            else
                deflections_aeromesh(1,:)=interp1(obj.wingbox_rl_coords,def(1,:), c4_rl_coords_edge,'lin','extrap');
                deflections_aeromesh(2,:)=interp1(obj.wingbox_rl_coords,def(2,:), c4_rl_coords_edge,'lin','extrap');
                deflections_aeromesh(3,:)=interp1(obj.wingbox_rl_coords,def(3,:), c4_rl_coords_edge,'lin','extrap');
                deflections_aeromesh(4,:)=interp1(obj.wingbox_rl_coords,def(4,:), c4_rl_coords_edge,'lin','extrap');
                deflections_aeromesh(5,:)=interp1(obj.wingbox_rl_coords,def(5,:), c4_rl_coords_edge,'lin','extrap');
                deflections_aeromesh(6,:)=interp1(obj.wingbox_rl_coords,def(6,:), c4_rl_coords_edge,'lin','extrap');
            end
            
            start_idx=1;
%                         hold on
%                         plot3(obj.wingbox_coords(1,:,1),obj.wingbox_coords(2,:,1),obj.wingbox_coords(3,:,1),'cx');
%             plot3(obj.wingbox_coords(1,:,1)+def(1,24:end),obj.wingbox_coords(2,:,1)+def(2,24:end),obj.wingbox_coords(3,:,1)+def(3,24:end),'ro');
% %             
            for i=1:length(obj.wing_segments)
                end_idx=start_idx+size(obj.wing_segments(i).c4_coords,2);
                grid=obj.wing_segments(i).compute_deflected_grid(panels,grid,deflections_aeromesh(:,start_idx:end_idx));
                start_idx=end_idx;
            end
            
            
%             if obj.grid_3D==1
%                 for i=1:length(obj.wing_segments)
%                         obj.grid_vol
%                 end  
%             end
            
%                         hold on
%                         for i=1:length(obj.panels)
%                             handle= fill3(grid(1,obj.panels(:,i)), grid(2,obj.panels(:,i)),grid(3,obj.panels(:,i)),'b');
%                             alpha(handle,0.4);
%                         end
%                         axis equal
            start_idx=1;
            if obj.symmetric==1
                for i=1:length(obj.wing_segments)
                    end_idx=start_idx+size(obj.wing_segments(i).c4_coords,2);
                    grid=obj.wing_segments(i).compute_deflected_grid_left(panels,grid,deflections_aeromesh_left(:,start_idx:end_idx),length(obj.panels)/2);
                    start_idx=end_idx;
                end
               % grid(:,obj.grid_start_idx:obj.grid_start_idx-1+size(obj.grid,2))=[grid(:,obj.grid_start_idx:obj.grid_start_idx-1+size(obj.grid,2)/2) [grid(1,obj.grid_start_idx:obj.grid_start_idx-1+size(obj.grid,2)/2);-grid(2,obj.grid_start_idx:obj.grid_start_idx-1+size(obj.grid,2)/2);grid(3,obj.grid_start_idx:obj.grid_start_idx-1+size(obj.grid,2)/2)]];
            end
          %  obj.grid_deflected=grid;
%             hold on
%                      plot3(obj.wingbox_coords(1,:,1),obj.wingbox_coords(2,:,1),obj.wingbox_coords(3,:,1),'x-');
%                      plot3(obj.wingbox_coords(1,:,2),obj.wingbox_coords(2,:,2),obj.wingbox_coords(3,:,2),'x-');
%                      plot3(0.5*obj.wingbox_coords(1,:,1)+0.5*obj.wingbox_coords(1,:,2),0.5*obj.wingbox_coords(2,:,1)+0.5*obj.wingbox_coords(2,:,2),0.5*obj.wingbox_coords(3,:,1)+0.5*obj.wingbox_coords(3,:,2),'o-');
%                      figure
%                      plot(c4_rl_coords_edge,deflections_aeromesh(1:6,:))
        end
        
        function obj=compute_wingbox_coords(obj,le_max,varargin)
            frontspar=0.25*ones(1,length(obj.wing_segments)+1);
            rearspar=0.75*ones(1,length(obj.wing_segments)+1);
            
            if nargin==4
                frontspar=varargin{1};
                rearspar=varargin{2};
            end
            
            wingbox_coords=[];
            wingbox_c4=[];
            wingbox_height=[];
            
            if(length(frontspar)==2*length(obj.wing_segments))
                incr=2;
            else
                incr=1;
            end
            ctr=1;
            iNodes = 1;
           % incr
            for i=1:length(obj.wing_segments)
      %          i
                n=ceil(obj.wing_segments(i).b/le_max);
                if obj.wing_segments(i).nBeamelements>0
                    if obj.wing_segments(i).nBeamelements<n
                        disp('WARNING: The number of beamelements selected for this wingsegment results in a lower resolution than defined in the grid spacing!')
                    end
                    n = obj.wing_segments(i).nBeamelements;
                else
%                     if n<3
%                         n=3;
%                     end
                end
                %ctr;
%                 frontspar
%                 rearspar
                if obj.isExternalFEM==0
                    if nargin == 3
                        obj.wing_segments(i) = obj.wing_segments(i).compute_wingbox_coords([], [], n);
                    else
                        obj.wing_segments(i)=obj.wing_segments(i).compute_wingbox_coords(frontspar(ctr:ctr+1),rearspar(ctr:ctr+1),n);
                    end
                elseif obj.isExternalFEM==1
                    obj.wing_segments(i) = obj.wing_segments(i).read_wingbox_coords(obj.pathNodeCoords, n, iNodes);
                else
                    disp('ERROR: wing.isExternalFEM not correctly defined!')
                end
                
                if isempty(wingbox_coords)
                    wingbox_coords=[wingbox_coords obj.wing_segments(i).wingbox_coords];
                    wingbox_c4=[wingbox_c4 obj.wing_segments(i).wingbox_c4];
                    wingbox_height=[wingbox_height obj.wing_segments(i).wingbox_height'];
                else
                    wingbox_coords=[wingbox_coords obj.wing_segments(i).wingbox_coords(:,2:end,:)];
                    wingbox_c4=[wingbox_c4 obj.wing_segments(i).wingbox_c4(:,2:end)];
                    wingbox_height=[wingbox_height obj.wing_segments(i).wingbox_height(2:end,:)'];
                end
                ctr=ctr+incr;
                iNodes = iNodes+n;
            end
            obj.wingbox_coords=wingbox_coords;
            obj.wingbox_c4=wingbox_c4;
            obj.wingbox_height=wingbox_height';
            
            c4_rl_coord_sm(1)=0;
            for i=1:size(obj.wingbox_c4,2)-1
                c4_rl_coord_sm(i+1)=norm(obj.wingbox_c4(1:3,i+1)-obj.wingbox_c4(1:3,i));
            end
            
            c4_rl_coord_sm=cumsum(c4_rl_coord_sm);
            
            obj.wingbox_rl_coords=c4_rl_coord_sm; 
            
            for i=1:length(obj.wingbox_rl_coords)-1
                c4_rl_coord_sm_mid(i)=0.5* obj.wingbox_rl_coords(i)+0.5*(obj.wingbox_rl_coords(i+1));
            end
            

            obj.wingbox_rl_coords_mid=c4_rl_coord_sm_mid;
            
            obj.c4_coords=[];
            for i=1:length(obj.wing_segments)
                obj.c4_coords=[obj.c4_coords obj.wing_segments(i).c4_coords];
            end
            
            c4_rl_coord(1)=norm(obj.c4_coords(1:3,1))-norm(obj.wingbox_c4(1:3,1));
            for i=1:size(obj.c4_coords,2)-1
                c4_rl_coord(i+1)= norm(obj.c4_coords(1:3,i+1)-obj.c4_coords(1:3,i));
            end
            c4_rl_coord=cumsum(c4_rl_coord);
            obj.c4_rl_coords=c4_rl_coord;
            
            %             plot3(wingbox_coords(1,:,1),wingbox_coords(2,:,1),wingbox_coords(3,:,1),'x-');
            %             plot3(wingbox_coords(1,:,2),wingbox_coords(2,:,2),wingbox_coords(3,:,2),'x-');
            
        end
        
        function obj=plot_grid(obj,varargin)
            hold on
            if  ~isempty(varargin)
                if size(varargin{1},1)==size(obj.panels,2)
                    for i=1:length(obj.panels)
                        handle= fill3(obj.grid(1,obj.panels(:,i)), obj.grid(2,obj.panels(:,i)),obj.grid(3,obj.panels(:,i)),'b');
                        alpha(handle,varargin{1}(i));
                    end     
                end
            else
                for i=1:length(obj.panels)
                    handle= fill3(obj.grid(1,obj.panels(:,i)), obj.grid(2,obj.panels(:,i)),obj.grid(3,obj.panels(:,i)),'b');
                    alpha(handle,0.4);
                end
            end
        end
        
        function obj=plot_grid_vol(obj)
            hold on
            for i=1:length(obj.panels)
                if obj.is_te(i)==0
                    handle= fill3(obj.grid_vol_up(1,obj.panels(:,i)), obj.grid_vol_up(2,obj.panels(:,i)),obj.grid_vol_up(3,obj.panels(:,i)),'b');
                else
                    handle= fill3(obj.grid_vol_up(1,obj.panels(:,i)), obj.grid_vol_up(2,obj.panels(:,i)),obj.grid_vol_up(3,obj.panels(:,i)),'r');   
                end
                alpha(handle,1);
            end
            
            for i=1:length(obj.panels)
                 if obj.is_te(i)==0
                    handle= fill3(obj.grid_vol_lo(1,obj.panels(:,i)), obj.grid_vol_lo(2,obj.panels(:,i)),obj.grid_vol_lo(3,obj.panels(:,i)),'b');
                 else
                     handle= fill3(obj.grid_vol_lo(1,obj.panels(:,i)), obj.grid_vol_lo(2,obj.panels(:,i)),obj.grid_vol_lo(3,obj.panels(:,i)),'r');
                 end 
                alpha(handle,1);
            end    
        end
        
        
        function obj=plot_grid_flat(obj)
            hold on
            %fprintf(fileID,'VARIABLES = "X", "Y", "Z"\n');
            for i=1:length(obj.panels)
                handle= fill3(obj.grid_flat(1,obj.panels(:,i)), obj.grid_flat(2,obj.panels(:,i)),obj.grid_flat(3,obj.panels(:,i)),'r');
                alpha(handle,0.4);
            end
        end
        
        function obj=write_tecplot_grid(obj,fileID)
            for i=1:length(obj.wing_segments)
                obj.wing_segments(i).write_tecplot_grid(fileID,obj.name,i);
            end
        end
        

        
        function obj=read_xml_definition(obj,xmlstruct)
            di=0;
            if strcmp(xmlstruct.tag,'WING')
                obj.name=xmlstruct.attribs(1).value;
                obj.symmetric=str2double(xmlstruct.attribs(2).value);
                
                % Initializing the array which contains all control surface
                % class instances refering to the aerosurface (parent
                % control surfaces). Each column is a wingsegment, and each
                % row a different control surface (row1=TE,row2=LE,row3=SP)
                obj.control_surfaces = cell(3,1);
                obj.control_surfaces(:) = {0};
                %initializing the counter which will keep track of the
                %wingsegments
                seg_count = 0;
                
                for i=1:length(xmlstruct.child)
                    if strcmp(xmlstruct.child(i).tag,'EXTERNALFEM')
                        obj.isExternalFEM = 1;
                        for j=1:length(xmlstruct.child(i).child)
                            if strcmp(xmlstruct.child(i).child(j).tag,'PATHNODECOORDS')
                                obj.pathNodeCoords = xmlstruct.child(i).child(j).value;
                            elseif strcmp(xmlstruct.child(i).child(j).tag,'PATHMASSMATRIX')
                                obj.pathMassMatrix = xmlstruct.child(i).child(j).value;
                            elseif strcmp(xmlstruct.child(i).child(j).tag,'PATHSTIFFNESSMATRIX')
                                obj.pathStiffnessMatrix = xmlstruct.child(i).child(j).value;
                            end
                        end
                    end
                    if strcmp(xmlstruct.child(i).tag,'SEGMENT')
                        seg_count = seg_count + 1;
                        obj.control_surfaces(:,seg_count) = {0};
                        if (length(xmlstruct.child(i).child(1).child)==3)
                            n=length(obj.wing_segments);
                            if n~=0
                                obj.wing_segments(n+1)=class_wingsegment(xmlstruct.child(i));
                                obj.wing_segments(n+1).symmetric=obj.symmetric;
                            else
                                obj.wing_segments=class_wingsegment(xmlstruct.child(i));
                                obj.wing_segments.symmetric=obj.symmetric;
                            end
                        else
                            if  strcmp(xmlstruct.child(i).child(1).child(1).tag,'ATTACH_TO')
                                if strcmp(xmlstruct.child(i).child(1).child(1).value,'previous')
                                    attach_pos=obj.wing_segments(end).get_tip_c4point();
                                    xmlstruct.child(i).child(1).child(1).tag='X';
                                    xmlstruct.child(i).child(1).child(2).tag='Y';
                                    xmlstruct.child(i).child(1).child(3).tag='Z';
                                    xmlstruct.child(i).child(1).child(1).value=num2str(attach_pos(1));
                                    xmlstruct.child(i).child(1).child(2).value=num2str(attach_pos(2));
                                    xmlstruct.child(i).child(1).child(3).value=num2str(attach_pos(3));
                                    n=length(obj.wing_segments);
                                    obj.wing_segments(n+1)=class_wingsegment(xmlstruct.child(i));
                                    obj.wing_segments(n+1).symmetric=obj.symmetric;
                                end
                            end
                        end
                        
                        % Initializing variables which will be later
                        % necessary to split segments (to accomodate for
                        % all control surfaces)
                        CS_etas_arr = [0,1];
                        te_device_start = [];
                        te_device_end = [];
                        le_device_start = [];
                        le_device_end = [];
                        sp_device_start = [];
                        sp_device_end = [];
                        
                        
                        % Adding comments to each line of code would be redundant. Basically code below tries to read the
                        % control surfaces and beam structure from 3 different indices of the xml segment (indices 11,
                        % 12 and 13, since at most we will have 3 surfaces, one for each index).
                        
                        % Within each index try (11, 12 and 13) we check wether a TE, LE or spoiler control surface is
                        % present. If one is, it proceeds to check whether the attributes denoting its deflection symmetry
                        % and taper are defined. Once all data is collected from the XML, a local instance of the
                        % class_control_surface is created (nargin varies depending on how much stuff was defined within
                        % XML). For each surface, the start and end etas are saved, to be later used in order to
                        % segment the whole wing simultaneously (thus avoiding assignment of control surface to the
                        % wrong segment).
                        
                        try
                            cs_idx=11;
                            if strcmp(xmlstruct.child(i).child(cs_idx).tag,'CONTROL_SURFACE')
                                if strcmp(xmlstruct.child(i).child(cs_idx).attribs(2).value,'trailing_edge')
                                    try
                                        if strcmp(xmlstruct.child(i).child(cs_idx).attribs(3).value,'symmetric')
                                            symdefl=1;
                                        elseif strcmp(xmlstruct.child(i).child(cs_idx).attribs(3).value,'asymmetric')
                                            symdefl=0;
                                        end
                                        tapered=0;
                                        try
                                            if strcmp(xmlstruct.child(i).child(cs_idx).child(1).attribs(1).value,'tapered')
                                                tapered=1;
                                            elseif  strcmp(xmlstruct.child(i).child(cs_idx).child(1).attribs(1).value,'constant')
                                                tapered=0;
                                            end
                                        end
                                        te_device=class_control_surface_parent(xmlstruct.child(i).child(cs_idx).attribs(1).value,0,xmlstruct.child(i).child(cs_idx).child(1).value,symdefl,tapered);
                                        obj.control_surfaces(1,seg_count) = {te_device};
                                        
                                    catch
                                        tapered=0;
                                        try
                                            if strcmp(xmlstruct.child(i).child(cs_idx).child(1).attribs(1).value,'tapered')
                                                tapered=1;
                                            elseif  strcmp(xmlstruct.child(i).child(cs_idx).child(1).attribs(1).value,'constant')
                                                tapered=0;
                                            end
                                        end
                                        te_device=class_control_surface_parent(xmlstruct.child(i).child(cs_idx).attribs(1).value,0,xmlstruct.child(i).child(cs_idx).child(1).value,1,tapered);
                                        obj.control_surfaces(1,seg_count) = {te_device};
                                        
                                    end
                                    
                                    te_device_start = str2double(xmlstruct.child(i).child(cs_idx).child(2).value);
                                    te_device_end = str2double(xmlstruct.child(i).child(cs_idx).child(3).value);
                                    CS_etas_arr = [CS_etas_arr, te_device_start, te_device_end];
                                
                                    
                                elseif strcmp(xmlstruct.child(i).child(cs_idx).attribs(2).value,'leading_edge')
                                    try
                                        if strcmp(xmlstruct.child(i).child(cs_idx).attribs(3).value,'symmetric')
                                            symdefl=1;
                                        elseif strcmp(xmlstruct.child(i).child(cs_idx).attribs(3).value,'asymmetric')
                                            symdefl=0;
                                        end
                                        
                                        tapered=0;
                                        try
                                            if strcmp(xmlstruct.child(i).child(cs_idx).child(1).attribs(1).value,'tapered')
                                                tapered=1;
                                            elseif  strcmp(xmlstruct.child(i).child(cs_idx).child(1).attribs(1).value,'constant')
                                                tapered=0;
                                            end
                                        end
                                        
                                        le_device=class_control_surface_parent(xmlstruct.child(i).child(cs_idx).attribs(1).value,1,xmlstruct.child(i).child(cs_idx).child(1).value,symdefl,tapered);
                                        obj.control_surfaces(2,seg_count) = {le_device};
                                        
                                    catch
                                        
                                        tapered=0;
                                        try
                                            if strcmp(xmlstruct.child(i).child(cs_idx).child(1).attribs(1).value,'tapered')
                                                tapered=1;
                                            elseif  strcmp(xmlstruct.child(i).child(cs_idx).child(1).attribs(1).value,'constant')
                                                tapered=0;
                                            end
                                        end
                                        
                                        le_device=class_control_surface_parent(xmlstruct.child(i).child(cs_idx).attribs(1).value,1,xmlstruct.child(i).child(cs_idx).child(1).value,1,tapered);
                                        obj.control_surfaces(2,seg_count) = {le_device};
                                        
                                    end
                                    
                                    le_device_start = str2double(xmlstruct.child(i).child(cs_idx).child(2).value);
                                    le_device_end = str2double(xmlstruct.child(i).child(cs_idx).child(3).value);
                                    CS_etas_arr = [CS_etas_arr, le_device_start, le_device_end];
                                    
                                elseif strcmp(xmlstruct.child(i).child(cs_idx).attribs(2).value,'center_wing')
                                    try
                                        if strcmp(xmlstruct.child(i).child(cs_idx).attribs(3).value,'symmetric')
                                            symdefl=1;
                                        elseif strcmp(xmlstruct.child(i).child(cs_idx).attribs(3).value,'asymmetric')
                                            symdefl=0;
                                        end
                                        
                                        tapered=0;
                                        try
                                            if strcmp(xmlstruct.child(i).child(cs_idx).child(1).attribs(1).value,'tapered')
                                                tapered=1;
                                            elseif  strcmp(xmlstruct.child(i).child(cs_idx).child(1).attribs(1).value,'constant')
                                                tapered=0;
                                            end
                                        end
                                        
                                        sp_device=class_control_surface_parent(xmlstruct.child(i).child(cs_idx).attribs(1).value,2,xmlstruct.child(i).child(cs_idx).child(1).value,symdefl,xmlstruct.child(i).child(cs_idx).child(4).value,xmlstruct.child(i).child(cs_idx).child(5).value,tapered);
                                        obj.control_surfaces(3,seg_count) = {sp_device};
                                        
                                    catch
                                        sp_device=class_control_surface_parent(xmlstruct.child(i).child(cs_idx).attribs(1).value,2,xmlstruct.child(i).child(cs_idx).child(1).value,xmlstruct.child(i).child(cs_idx).child(4).value,xmlstruct.child(i).child(cs_idx).child(5).value);
                                        obj.control_surfaces(3,seg_count) = {sp_device};
                                        
                                    end
                                    
                                    sp_device_start = str2double(xmlstruct.child(i).child(cs_idx).child(2).value);
                                    sp_device_end = str2double(xmlstruct.child(i).child(cs_idx).child(3).value);
                                    CS_etas_arr = [CS_etas_arr, sp_device_start, sp_device_end];
                                end
                                
                            % NBEAMELEMENTS defines the number of the wingsegment's beamelements. This number (if defined) will be used instead of the grid spacing
                            elseif strcmp(xmlstruct.child(i).child(cs_idx).tag,'NBEAMELEMENTS')  
                                obj.wing_segments(i+di).nBeamelements = str2num(xmlstruct.child(i).child(cs_idx).value);
                            end
                        end
                        
                        try
                            cs_idx=12;
                            if strcmp(xmlstruct.child(i).child(cs_idx).tag,'CONTROL_SURFACE')
                                if strcmp(xmlstruct.child(i).child(cs_idx).attribs(2).value,'trailing_edge')
                                    try
                                        if strcmp(xmlstruct.child(i).child(cs_idx).attribs(3).value,'symmetric')
                                            symdefl=1;
                                        elseif strcmp(xmlstruct.child(i).child(cs_idx).attribs(3).value,'asymmetric')
                                            symdefl=0;
                                        end
                                        tapered=0;
                                        try
                                            if strcmp(xmlstruct.child(i).child(cs_idx).child(1).attribs(1).value,'tapered')
                                                tapered=1;
                                            elseif  strcmp(xmlstruct.child(i).child(cs_idx).child(1).attribs(1).value,'constant')
                                                tapered=0;
                                            end
                                        end
                                        te_device=class_control_surface_parent(xmlstruct.child(i).child(cs_idx).attribs(1).value,0,xmlstruct.child(i).child(cs_idx).child(1).value,symdefl,tapered);
                                        obj.control_surfaces(1,seg_count) = {te_device};
                                        
                                    catch
                                        tapered=0;
                                        try
                                            if strcmp(xmlstruct.child(i).child(cs_idx).child(1).attribs(1).value,'tapered')
                                                tapered=1;
                                            elseif  strcmp(xmlstruct.child(i).child(cs_idx).child(1).attribs(1).value,'constant')
                                                tapered=0;
                                            end
                                        end
                                        te_device=class_control_surface_parent(xmlstruct.child(i).child(cs_idx).attribs(1).value,0,xmlstruct.child(i).child(cs_idx).child(1).value,1,tapered);
                                        obj.control_surfaces(1,seg_count) = {te_device};
                                       
                                    end
                                    
                                    te_device_start = str2double(xmlstruct.child(i).child(cs_idx).child(2).value);
                                    te_device_end = str2double(xmlstruct.child(i).child(cs_idx).child(3).value);
                                    CS_etas_arr = [CS_etas_arr, te_device_start, te_device_end];
                                    
                                elseif strcmp(xmlstruct.child(i).child(cs_idx).attribs(2).value,'leading_edge')
                                    try
                                        if strcmp(xmlstruct.child(i).child(cs_idx).attribs(3).value,'symmetric')
                                            symdefl=1;
                                        elseif strcmp(xmlstruct.child(i).child(cs_idx).attribs(3).value,'asymmetric')
                                            symdefl=0;
                                        end
                                        
                                        tapered=0;
                                        try
                                            if strcmp(xmlstruct.child(i).child(cs_idx).child(1).attribs(1).value,'tapered')
                                                tapered=1;
                                            elseif  strcmp(xmlstruct.child(i).child(cs_idx).child(1).attribs(1).value,'constant')
                                                tapered=0;
                                            end
                                        end
                                        le_device=class_control_surface_parent(xmlstruct.child(i).child(cs_idx).attribs(1).value,1,xmlstruct.child(i).child(cs_idx).child(1).value,symdefl,tapered);
                                        obj.control_surfaces(2,seg_count) = {le_device};
                                        
                                    catch
                                        
                                        tapered=0;
                                        try
                                            if strcmp(xmlstruct.child(i).child(cs_idx).child(1).attribs(1).value,'tapered')
                                                tapered=1;
                                            elseif  strcmp(xmlstruct.child(i).child(cs_idx).child(1).attribs(1).value,'constant')
                                                tapered=0;
                                            end
                                        end
                                        le_device=class_control_surface_parent(xmlstruct.child(i).child(cs_idx).attribs(1).value,1,xmlstruct.child(i).child(cs_idx).child(1).value,1,tapered);
                                        obj.control_surfaces(2,seg_count) = {le_device};
                                       
                                    end
                                    
                                    le_device_start = str2double(xmlstruct.child(i).child(cs_idx).child(2).value);
                                    le_device_end = str2double(xmlstruct.child(i).child(cs_idx).child(3).value);
                                    CS_etas_arr = [CS_etas_arr, le_device_start, le_device_end];
                                    
                                elseif strcmp(xmlstruct.child(i).child(cs_idx).attribs(2).value,'center_wing')
                                    try
                                        if strcmp(xmlstruct.child(i).child(cs_idx).attribs(3).value,'symmetric')
                                            symdefl=1;
                                        elseif strcmp(xmlstruct.child(i).child(cs_idx).attribs(3).value,'asymmetric')
                                            symdefl=0;
                                        end
                                        
                                        tapered=0;
                                        try
                                            if strcmp(xmlstruct.child(i).child(cs_idx).child(1).attribs(1).value,'tapered')
                                                tapered=1;
                                            elseif  strcmp(xmlstruct.child(i).child(cs_idx).child(1).attribs(1).value,'constant')
                                                tapered=0;
                                            end
                                        end
                                        
                                        sp_device=class_control_surface_parent(xmlstruct.child(i).child(cs_idx).attribs(1).value,2,xmlstruct.child(i).child(cs_idx).child(1).value,symdefl,xmlstruct.child(i).child(cs_idx).child(4).value,xmlstruct.child(i).child(cs_idx).child(5).value,tapered);
                                        obj.control_surfaces(3,seg_count) = {sp_device};
                                        
                                    catch
                                        sp_device=class_control_surface_parent(xmlstruct.child(i).child(cs_idx).attribs(1).value,2,xmlstruct.child(i).child(cs_idx).child(1).value,xmlstruct.child(i).child(cs_idx).child(4).value,xmlstruct.child(i).child(cs_idx).child(5).value);
                                        obj.control_surfaces(3,seg_count) = {sp_device};
                                        
                                    end
                                    
                                    sp_device_start = str2double(xmlstruct.child(i).child(cs_idx).child(2).value);
                                    sp_device_end = str2double(xmlstruct.child(i).child(cs_idx).child(3).value);
                                    CS_etas_arr = [CS_etas_arr, sp_device_start, sp_device_end];
                                end
                            elseif strcmp(xmlstruct.child(i).child(cs_idx).tag,'NBEAMELEMENTS')  
                                obj.wing_segments(i+di).nBeamelements = str2num(xmlstruct.child(i).child(cs_idx).value);
                            end  
                        end
                        
                        
                        try
                            cs_idx=13;
                            if strcmp(xmlstruct.child(i).child(cs_idx).tag,'CONTROL_SURFACE')
                                if strcmp(xmlstruct.child(i).child(cs_idx).attribs(2).value,'trailing_edge')
                                    try
                                        if strcmp(xmlstruct.child(i).child(cs_idx).attribs(3).value,'symmetric')
                                            symdefl=1;
                                        elseif strcmp(xmlstruct.child(i).child(cs_idx).attribs(3).value,'asymmetric')
                                            symdefl=0;
                                        end
                                        tapered=0;
                                        try
                                            if strcmp(xmlstruct.child(i).child(cs_idx).child(1).attribs(1).value,'tapered')
                                                tapered=1;
                                            elseif  strcmp(xmlstruct.child(i).child(cs_idx).child(1).attribs(1).value,'constant')
                                                tapered=0;
                                            end
                                        end
                                        te_device=class_control_surface_parent(xmlstruct.child(i).child(cs_idx).attribs(1).value,0,xmlstruct.child(i).child(cs_idx).child(1).value,symdefl,tapered);
                                        obj.control_surfaces(1,seg_count) = {te_device};
                                        
                                    catch
                                        tapered=0;
                                        try
                                            if strcmp(xmlstruct.child(i).child(cs_idx).child(1).attribs(1).value,'tapered')
                                                tapered=1;
                                            elseif  strcmp(xmlstruct.child(i).child(cs_idx).child(1).attribs(1).value,'constant')
                                                tapered=0;
                                            end
                                        end
                                        te_device=class_control_surface_parent(xmlstruct.child(i).child(cs_idx).attribs(1).value,0,xmlstruct.child(i).child(cs_idx).child(1).value,1,tapered);
                                        obj.control_surfaces(1,seg_count) = {te_device};
                                       
                                    end
                                    
                                    te_device_start = str2double(xmlstruct.child(i).child(cs_idx).child(2).value);
                                    te_device_end = str2double(xmlstruct.child(i).child(cs_idx).child(3).value);
                                    CS_etas_arr = [CS_etas_arr, te_device_start, te_device_end];
                                    
                                elseif strcmp(xmlstruct.child(i).child(cs_idx).attribs(2).value,'leading_edge')
                                    try
                                        if strcmp(xmlstruct.child(i).child(cs_idx).attribs(3).value,'symmetric')
                                            symdefl=1;
                                        elseif strcmp(xmlstruct.child(i).child(cs_idx).attribs(3).value,'asymmetric')
                                            symdefl=0;
                                        end
                                        
                                        tapered=0;
                                        try
                                            if strcmp(xmlstruct.child(i).child(cs_idx).child(1).attribs(1).value,'tapered')
                                                tapered=1;
                                            elseif  strcmp(xmlstruct.child(i).child(cs_idx).child(1).attribs(1).value,'constant')
                                                tapered=0;
                                            end
                                        end
                                        
                                        le_device=class_control_surface_parent(xmlstruct.child(i).child(cs_idx).attribs(1).value,1,xmlstruct.child(i).child(cs_idx).child(1).value,symdefl,tapered);
                                        obj.control_surfaces(2,seg_count) = {le_device};
                                        
                                        
                                    catch
                                        tapered=0;
                                        try
                                            if strcmp(xmlstruct.child(i).child(cs_idx).child(1).attribs(1).value,'tapered')
                                                tapered=1;
                                            elseif  strcmp(xmlstruct.child(i).child(cs_idx).child(1).attribs(1).value,'constant')
                                                tapered=0;
                                            end
                                        end
                                        le_device=class_control_surface_parent(xmlstruct.child(i).child(cs_idx).attribs(1).value,1,xmlstruct.child(i).child(cs_idx).child(1).value,1,tapered);
                                        obj.control_surfaces(2,seg_count) = {le_device};
                                        
                                    end
                                    
                                    le_device_start = str2double(xmlstruct.child(i).child(cs_idx).child(2).value);
                                    le_device_end = str2double(xmlstruct.child(i).child(cs_idx).child(3).value);
                                    CS_etas_arr = [CS_etas_arr, le_device_start, le_device_end];
                                    
                                elseif strcmp(xmlstruct.child(i).child(cs_idx).attribs(2).value,'center_wing')
                                    try
                                        if strcmp(xmlstruct.child(i).child(cs_idx).attribs(3).value,'symmetric')
                                            symdefl=1;
                                        elseif strcmp(xmlstruct.child(i).child(cs_idx).attribs(3).value,'asymmetric')
                                            symdefl=0;
                                        end
                                        tapered=0;
                                        try
                                            if strcmp(xmlstruct.child(i).child(cs_idx).child(1).attribs(1).value,'tapered')
                                                tapered=1;
                                            elseif  strcmp(xmlstruct.child(i).child(cs_idx).child(1).attribs(1).value,'constant')
                                                tapered=0;
                                            end
                                        end
                                        
                                        sp_device=class_control_surface_parent(xmlstruct.child(i).child(cs_idx).attribs(1).value,2,xmlstruct.child(i).child(cs_idx).child(1).value,symdefl,xmlstruct.child(i).child(cs_idx).child(4).value,xmlstruct.child(i).child(cs_idx).child(5).value,tapered);
                                        obj.control_surfaces(3,seg_count) = {sp_device};
                                        
                                    catch
                                        sp_device=class_control_surface_parent(xmlstruct.child(i).child(cs_idx).attribs(1).value,2,xmlstruct.child(i).child(cs_idx).child(1).value,xmlstruct.child(i).child(cs_idx).child(4).value,xmlstruct.child(i).child(cs_idx).child(5).value);
                                        obj.control_surfaces(3,seg_count) = {sp_device};
                                        
                                    end
                                    
                                    sp_device_start = str2double(xmlstruct.child(i).child(cs_idx).child(2).value);
                                    sp_device_end = str2double(xmlstruct.child(i).child(cs_idx).child(3).value);
                                    CS_etas_arr = [CS_etas_arr, sp_device_start, sp_device_end];
                                end
                                
                            elseif strcmp(xmlstruct.child(i).child(cs_idx).tag,'NBEAMELEMENTS')  
                                obj.wing_segments(i+di).nBeamelements = str2num(xmlstruct.child(i).child(cs_idx).value);
                            end  
                        end
                        
                        try
                            cs_idx=14;
                            if strcmp(xmlstruct.child(i).child(cs_idx).tag,'NBEAMELEMENTS')  
                                obj.wing_segments(i+di).nBeamelements = str2num(xmlstruct.child(i).child(cs_idx).value);
                            end  
                        end
                        
                        % puts all etas within the array of cuts to be
                        % applied in ascending order. Removes all duplicate
                        % elements (only 1 cut will be needed in same
                        % location). Removes 0 and 1 coordinates since
                        % obviously no cuts are needed there.
                        CS_etas_arr_ordered = sort(CS_etas_arr,'ascend');
                        CS_etas_arr_ordered = unique(CS_etas_arr_ordered);
                        CS_etas_cuts = CS_etas_arr_ordered(CS_etas_arr_ordered~=0);
                        CS_etas_cuts = CS_etas_cuts(CS_etas_cuts~=1);
                        
                        % resize_factor is used to still calculate coorect
                        % eta coords after change in dimensions due to
                        % segment splitting.
                        resize_factor = 1;
                        
                        % di is the current segment counter
                        % within global reference frame of aerosurface, not
                        % just local frame.
                        
                        % di_base is the di before any segments are split
                        
                        di_base = di;
                        prev_eta=0;
                        
                        % loop that cuts the segment at CS_etas_cuts
                        % locations by splitting them.
                        for cut_count = 1:length(CS_etas_cuts)
                            %splits the segment. prev_eta and resize-factor
                            %are used to ensure that CS_etas_cuts is
                            %translated into the new coordinate system
                            %(segment splits change geometry).
                            obj=obj.split_segment(i+di, 0, (CS_etas_cuts(cut_count)-prev_eta)/resize_factor);
                            prev_eta = CS_etas_cuts(cut_count);
                            resize_factor = (1-CS_etas_cuts(cut_count));
                            di = di+1;
                        end
                        
                        
                        % Reading the names array, to check
                        % whether the newly added surfaces are simply a
                        % continuation of previously added ones, in which
                        % case they will be removed from the parent array,
                        % but still integrated as children within parents,
                        % and of course added to their corresponding wing
                        % segment
                        if seg_count>1  %if seg_count==1, current surface must be parent
                            
                            %looping through array to check TE, LE and SP
                            for cs_count = 1:size(obj.control_surfaces,1)
                                for seg_iter = 1:seg_count-1
                                    % if the names are the same, it means it is  the same surface, in which case the
                                    % corresponding entry in the control_surfaces array is turned into
                                    % the index of the segment where the parent is stored
                                    if ~isnumeric(obj.control_surfaces{cs_count,seg_count})
                                        if ~isnumeric(obj.control_surfaces{cs_count,seg_iter})
                                            if strcmp(obj.control_surfaces{cs_count,seg_count}.name,obj.control_surfaces{cs_count,seg_iter}.name)
                                                obj.control_surfaces(cs_count,seg_count) = {seg_iter};
                                            end
                                        end
                                    end
                                end
                            end
                        end                      
                        
                        % runs if TE surfaces are present
                        if isempty(te_device_start)==0
                            
                            % loops over CS_etas_arr_ordered array to find
                            % start and end segment indices of CS
                            for temp_count = 1:length(CS_etas_arr_ordered)
                                if CS_etas_arr_ordered(temp_count) == te_device_start
                                    te_device_start_idx = (temp_count-1) + di_base + i;

                                elseif CS_etas_arr_ordered(temp_count) == te_device_end
                                    te_device_end_idx = (temp_count-2) + di_base + i;
                                end
                            end
                            
                            % if untapered, the root chord of the parent CS is passed to all children CS instances
                            % created within multiple segments. This is to ensure that the chord is effectively constant
                            % over multiple segments. See add_control_surface function for more detail
                            if te_device.is_tapered==0
                                
                                c_r_parent = obj.wing_segments(te_device_start_idx).c_r;
                                
                                % Loops over all segments containing parent CS, in order to add the children CS to each of them
                                for loc_seg_count = te_device_start_idx:te_device_end_idx
                                    
                                    % if the current entry in control_surface array is a number, it means that this surface cannot be
                                    % considered the parent, and the parent is stored at the index indicated by the entry.
                                    if isnumeric(obj.control_surfaces{1,seg_count})
                                        parent_idx = obj.control_surfaces{1,seg_count};
                                        child_idx = length(obj.control_surfaces{1,parent_idx}.children);
                                        obj.control_surfaces{1,parent_idx}.children{1,child_idx+1} = class_control_surface(te_device,child_idx+1,loc_seg_count);
                                        obj.wing_segments(loc_seg_count) = obj.wing_segments(loc_seg_count).add_control_surface(obj.control_surfaces{1,parent_idx}.children{1,child_idx+1},c_r_parent);
                                    
                                        % if the entry is not numeric, then the current entry is also the parent of all children applied to this segment
                                    else
                                        child_idx = length(obj.control_surfaces{1,seg_count}.children);
                                        obj.control_surfaces{1,seg_count}.children{1,child_idx+1} = class_control_surface(te_device,child_idx+1,loc_seg_count);
                                        obj.wing_segments(loc_seg_count) = obj.wing_segments(loc_seg_count).add_control_surface(obj.control_surfaces{1,seg_count}.children{1,child_idx+1},c_r_parent);
                                    end
                                end
                                
                            % if the device is tapered, there is no need to pass the parent's root chord to the children.
                            elseif te_device.is_tapered==1
                                for loc_seg_count = te_device_start_idx:te_device_end_idx
                                    % if the current entry in control_surface array is a number, it means that this surface cannot be
                                    % considered the parent, and the parent is stored at the index indicated by the entry.
                                    if isnumeric(obj.control_surfaces{1,seg_count})
                                        parent_idx = obj.control_surfaces{1,seg_count};
                                        child_idx = length(obj.control_surfaces{1,parent_idx}.children);
                                        obj.control_surfaces{1,parent_idx}.children{1,child_idx+1} = class_control_surface(te_device,child_idx+1,loc_seg_count);
                                        obj.wing_segments(loc_seg_count) = obj.wing_segments(loc_seg_count).add_control_surface(obj.control_surfaces{1,parent_idx}.children{1,child_idx+1});
                                    
                                        % if the entry is not numeric, then the current entry is also the parent of all children applied to this segment
                                    else
                                        child_idx = length(obj.control_surfaces{1,seg_count}.children);
                                        obj.control_surfaces{1,seg_count}.children{1,child_idx+1} = class_control_surface(te_device,child_idx+1,loc_seg_count);
                                        obj.wing_segments(loc_seg_count) = obj.wing_segments(loc_seg_count).add_control_surface(obj.control_surfaces{1,seg_count}.children{1,child_idx+1});
                                    end
                                end
                            end
                        end
                        
                        % runs if LE surface is present. Code is identical
                        % to TE surface case above, see comments of the TE
                        % case for more details.
                        if isempty(le_device_start)==0
                            for temp_count = 1:length(CS_etas_arr_ordered)
                                if CS_etas_arr_ordered(temp_count) == le_device_start
                                    le_device_start_idx = (temp_count-1) + di_base + i;

                                elseif CS_etas_arr_ordered(temp_count) == le_device_end
                                    le_device_end_idx = (temp_count-2) + di_base + i;
                                end
                            end

                            if le_device.is_tapered==0
                                
                                c_r_parent = obj.wing_segments(le_device_start_idx).c_r;
                                
                                for loc_seg_count = le_device_start_idx:le_device_end_idx
                                    
                                    % if the current entry in control_surface array is a number, it means that this surface cannot be
                                    % considered the parent, and the parent is stored at the index indicated by the entry.
                                    if isnumeric(obj.control_surfaces{2,seg_count})
                                        parent_idx = obj.control_surfaces{2,seg_count};
                                        child_idx = length(obj.control_surfaces{2,parent_idx}.children);
                                        obj.control_surfaces{2,parent_idx}.children{1,child_idx+1} = class_control_surface(le_device,child_idx+1,loc_seg_count);
                                        obj.wing_segments(loc_seg_count) = obj.wing_segments(loc_seg_count).add_control_surface(obj.control_surfaces{2,parent_idx}.children{1,child_idx+1},c_r_parent);
                                    
                                        % if the entry is not numeric, then the current entry is also the parent of all children applied to this segment
                                    else
                                        child_idx = length(obj.control_surfaces{2,seg_count}.children);
                                        obj.control_surfaces{2,seg_count}.children{1,child_idx+1} = class_control_surface(le_device,child_idx+1,loc_seg_count);
                                        obj.wing_segments(loc_seg_count) = obj.wing_segments(loc_seg_count).add_control_surface(obj.control_surfaces{2,seg_count}.children{1,child_idx+1},c_r_parent);
                                    end
                                end
                                
                            elseif le_device.is_tapered==1
                                for loc_seg_count = le_device_start_idx:le_device_end_idx
                                    
                                    % if the current entry in control_surface array is a number, it means that this surface cannot be
                                    % considered the parent, and the parent is stored at the index indicated by the entry.
                                    if isnumeric(obj.control_surfaces{2,seg_count})
                                        parent_idx = obj.control_surfaces{2,seg_count};
                                        child_idx = length(obj.control_surfaces{2,parent_idx}.children);
                                        obj.control_surfaces{2,parent_idx}.children{1,child_idx+1} = class_control_surface(le_device,child_idx+1,loc_seg_count);
                                        obj.wing_segments(loc_seg_count) = obj.wing_segments(loc_seg_count).add_control_surface(obj.control_surfaces{2,parent_idx}.children{1,child_idx+1});
                                    
                                        % if the entry is not numeric, then the current entry is also the parent of all children applied to this segment
                                    else
                                        child_idx = length(obj.control_surfaces{2,seg_count}.children);
                                        obj.control_surfaces{2,seg_count}.children{1,child_idx+1} = class_control_surface(le_device,child_idx+1,loc_seg_count);
                                        obj.wing_segments(loc_seg_count) = obj.wing_segments(loc_seg_count).add_control_surface(obj.control_surfaces{2,seg_count}.children{1,child_idx+1});
                                    end
                                end
                            end
                        end
                        
                        
                        % runs if spoiler surface is present. Code is identical
                        % to TE surface case above, see comments of the TE
                        % case for more details.
                        if isempty(sp_device_start)==0
                            for temp_count = 1:length(CS_etas_arr_ordered)
                                if CS_etas_arr_ordered(temp_count) == sp_device_start
                                    sp_device_start_idx = (temp_count-1) + di_base + i;

                                elseif CS_etas_arr_ordered(temp_count) == sp_device_end
                                    sp_device_end_idx = (temp_count-2) + di_base + i;
                                end
                            end

                            if sp_device.is_tapered==0
                                
                                c_r_parent = obj.wing_segments(sp_device_start_idx).c_r;
                                
                                for loc_seg_count = sp_device_start_idx:sp_device_end_idx
                                    
                                    % if the current entry in control_surface array is a number, it means that this surface cannot be
                                    % considered the parent, and the parent is stored at the index indicated by the entry.
                                    if isnumeric(obj.control_surfaces{3,seg_count})
                                        parent_idx = obj.control_surfaces{3,seg_count};
                                        child_idx = length(obj.control_surfaces{3,parent_idx}.children);
                                        obj.control_surfaces{3,parent_idx}.children{1,child_idx+1} = class_control_surface(sp_device,child_idx+1,loc_seg_count);
                                        obj.wing_segments(loc_seg_count) = obj.wing_segments(loc_seg_count).add_control_surface(obj.control_surfaces{3,parent_idx}.children{1,child_idx+1},c_r_parent);
                                    
                                        % if the entry is not numeric, then the current entry is also the parent of all children applied to this segment
                                    else
                                        child_idx = length(obj.control_surfaces{3,seg_count}.children);
                                        obj.control_surfaces{3,seg_count}.children{1,child_idx+1} = class_control_surface(sp_device,child_idx+1,loc_seg_count);
                                        obj.wing_segments(loc_seg_count) = obj.wing_segments(loc_seg_count).add_control_surface(obj.control_surfaces{3,seg_count}.children{1,child_idx+1},c_r_parent);
                                    end
                                end
                                
                            elseif sp_device.is_tapered==1
                                for loc_seg_count = sp_device_start_idx:sp_device_end_idx
                                    
                                    % if the current entry in control_surface array is a number, it means that this surface cannot be
                                    % considered the parent, and the parent is stored at the index indicated by the entry.
                                    if isnumeric(obj.control_surfaces{3,seg_count})
                                        parent_idx = obj.control_surfaces{3,seg_count};
                                        child_idx = length(obj.control_surfaces{3,parent_idx}.children);
                                        obj.control_surfaces{3,parent_idx}.children{1,child_idx+1} = class_control_surface(sp_device,child_idx+1,loc_seg_count);
                                        obj.wing_segments(loc_seg_count) = obj.wing_segments(loc_seg_count).add_control_surface(obj.control_surfaces{3,parent_idx}.children{1,child_idx+1});
                                    
                                        % if the entry is not numeric, then the current entry is also the parent of all children applied to this segment
                                    else
                                        child_idx = length(obj.control_surfaces{3,seg_count}.children);
                                        obj.control_surfaces{3,seg_count}.children{1,child_idx+1} = class_control_surface(sp_device,child_idx+1,loc_seg_count);
                                        obj.wing_segments(loc_seg_count) = obj.wing_segments(loc_seg_count).add_control_surface(obj.control_surfaces{3,seg_count}.children{1,child_idx+1});
                                    end
                                end
                            end
                        end
                    end
                end
                
                
                
                % once looping over all wing-segments is finisehd, we need
                % to unravel the control_surfaces array and remove all
                % numeric elements
                
                % unraveling
                obj.control_surfaces = reshape(obj.control_surfaces,[1,size(obj.control_surfaces,2)*size(obj.control_surfaces,1)]);
                
                % removing numeric elements, thus leaving only parent
                % surfaces within the array
                numeric_logic = cellfun(@isnumeric,obj.control_surfaces);
                obj.control_surfaces = obj.control_surfaces(~numeric_logic);
                
                
                % converting obj.control_surfaces back into an array for
                % consistency with the rest of the code
                
                if ~isempty(obj.control_surfaces)
                    
                    for temp_count = 1:length(obj.control_surfaces)
                        
                        if obj.symmetric
                            obj.control_surfaces{temp_count}.is_sym = 1;
                        else
                            obj.control_surfaces{temp_count}.is_sym = 0;
                        end
                        temp_array(temp_count) = obj.control_surfaces{temp_count};
                    end
                    
                    obj.control_surfaces = temp_array;
                    
                else
                    obj.control_surfaces = [];
                end
                
            else
                fprintf('Unknown Data Format: %s \n', xmlstruct.tag);
            end
        end
        
        function obj=  computeControlSurfacePanelIds(obj)
            
            for iSeg=1:length(obj.wing_segments)
                obj.wing_segments(iSeg)=obj.wing_segments(iSeg).computeControlSurfacePanelIds();
            end
            
            % preparation of spanwise and chordwise stations for full wing
            if ~isempty(obj.parCS)
                nPanSeg=[];
                for iSeg=1:length(obj.wing_segments)
                    nPanSeg=[nPanSeg size(obj.wing_segments(iSeg).panels,2)];
                end
                
                normPrvMidlines=0;
                spanWisePointsI=[];
                spanWisePointsO=[];
                chordWisePointsF=[];
                chordWisePointsR=[];
                
                for iSeg=1:length(obj.wing_segments)
                    
                    nPrvPanels=0;
                    for jSeg=1:iSeg-1
                        nPrvPanels=nPrvPanels+(nPanSeg(jSeg));
                    end

                    nPanels=nPanSeg(iSeg);
                    nSpan=sum(obj.is_te(nPrvPanels+1:nPrvPanels+nPanels));
                    nChord=nPanels/nSpan;
                    innerLEPanel=1+nPrvPanels;
                    outerLEPanel=(nSpan-1)*nChord+1+nPrvPanels;
                    innerTEPanel=nChord+nPrvPanels;
                    outerTEPanel=nPanels+nPrvPanels;

                    midLine=(obj.grid(:,obj.panels(2,outerLEPanel))+obj.grid(:,obj.panels(3,outerTEPanel)))./2-(obj.grid(:,obj.panels(1,innerLEPanel))+obj.grid(:,obj.panels(4,innerTEPanel)))./2;
                    linPoints=normPrvMidlines+linspace(0,1,nSpan+1)*norm(midLine);
                    spanWisePointsI=[spanWisePointsI reshape(repmat(linPoints(1:end-1),nChord,1),1,nPanels)];
                    spanWisePointsO=[spanWisePointsO reshape(repmat(linPoints(2:end),nChord,1),1,nPanels)];
                    normPrvMidlines=normPrvMidlines+norm(midLine);
                    for iChord=1:nSpan
                        %midpoints of panels in this strip
                        panelIdx=((iChord-1)*nChord+nPrvPanels+1:(iChord)*nChord+nPrvPanels);
                        midpoints=(obj.grid(:,obj.panels(1,panelIdx))+obj.grid(:,obj.panels(2,panelIdx))+obj.grid(:,obj.panels(3,panelIdx))+obj.grid(:,obj.panels(4,panelIdx)))./4;
                        Fmidpoints=(obj.grid(:,obj.panels(1,panelIdx))+obj.grid(:,obj.panels(2,panelIdx)))./2;
                        Rmidpoints=(obj.grid(:,obj.panels(3,panelIdx))+obj.grid(:,obj.panels(4,panelIdx)))./2;
                        dist=midpoints(:,2:end)-midpoints(:,1);
                        distF=Fmidpoints(:,2:end)-Fmidpoints(:,1);
                        distR=Rmidpoints(:,2:end)-Fmidpoints(:,1);
                        totDist=norm(Rmidpoints(:,end)-Fmidpoints(:,1));
                        chordWisePointsF=[chordWisePointsF [0 sqrt(sum(distF.^2))./totDist]];
                        chordWisePointsR=[chordWisePointsR [sqrt(sum(distF(:,1).^2))/totDist sqrt(sum(distR.^2))./totDist]];
                    end

                %     hold on; scatter3(obj.grid(1,obj.panels(:,innerLEPanel)),obj.grid(2,obj.panels(:,innerLEPanel)),obj.grid(3,obj.panels(:,innerLEPanel)),'o')
                %     hold on; scatter3(obj.grid(1,obj.panels(:,outerLEPanel)),obj.grid(2,obj.panels(:,outerLEPanel)),obj.grid(3,obj.panels(:,outerLEPanel)),'+')
                %     hold on; scatter3(obj.grid(1,obj.panels(:,innerTEPanel)),obj.grid(2,obj.panels(:,innerTEPanel)),obj.grid(3,obj.panels(:,innerTEPanel)),'d')
                %     hold on; scatter3(obj.grid(1,obj.panels(:,outerTEPanel)),obj.grid(2,obj.panels(:,outerTEPanel)),obj.grid(3,obj.panels(:,outerTEPanel)),'s')

                end
                    spanWisePointsI=spanWisePointsI/normPrvMidlines;
                    spanWisePointsO=spanWisePointsO/normPrvMidlines;
                    spanWiseLength=normPrvMidlines;
                    
            end
            
            
            for iParCS=1:length(obj.parCS)
                %check if panels are partially inside
                %panels are partially inside when one of the F/R points or
                %I/O points are inside
                pointsR= and(chordWisePointsR<=obj.parCS{iParCS}.posChordR, chordWisePointsR>=obj.parCS{iParCS}.posChordF);
                pointsF= and(chordWisePointsF<=obj.parCS{iParCS}.posChordR, chordWisePointsF>=obj.parCS{iParCS}.posChordF);
                pointsI= and(spanWisePointsI<=obj.parCS{iParCS}.posSpanOB,spanWisePointsI>=obj.parCS{iParCS}.posSpanIB);
                pointsO= and(spanWisePointsO<=obj.parCS{iParCS}.posSpanOB,spanWisePointsO>=obj.parCS{iParCS}.posSpanIB);
                panelids=and(or(and(pointsF,pointsR),xor(pointsF,pointsR)),or(and(pointsI,pointsO),xor(pointsI,pointsO)));
                
                fullyInside=and(and(pointsF,pointsR),and(pointsI,pointsO));
                %compute fractions of panels which are partially inside
                rect1=[chordWisePointsF' spanWisePointsI' (chordWisePointsR-chordWisePointsF)' (spanWisePointsO-spanWisePointsI)'];
                rect2=[obj.parCS{iParCS}.posChordF, obj.parCS{iParCS}.posSpanIB, obj.parCS{iParCS}.posChordR-obj.parCS{iParCS}.posChordF , obj.parCS{iParCS}.posSpanOB-obj.parCS{iParCS}.posSpanIB];
                allFract=rectint(rect1,rect2)./[rect1(:,3).*rect1(:,4)];
                obj.parCS{iParCS}.fraction=allFract(find(panelids));
             
                
                
                %calculate inner point
                %find id of panel closest to inner front hinge axis which
                %is the first panel in panelids
                id1=find(panelids,1);
                panelPoints=obj.grid(:,obj.panels(:,id1));
                %find fraction of span of this panel at which the control
                %surface starts
                fractChord=interp1([chordWisePointsF(id1) chordWisePointsR(id1)],[0 1],obj.parCS{iParCS}.posChordF);
                fractSpan=interp1([spanWisePointsI(id1) spanWisePointsO(id1)],[0 1],obj.parCS{iParCS}.posSpanIB);
                %calculate inner point
                innerHingePoint=(panelPoints(:,1)*(1-fractChord)+panelPoints(:,4)*fractChord)*(1-fractSpan)+(panelPoints(:,2)*(1-fractChord)+panelPoints(:,3)*fractChord)*(fractSpan);
                
                %calculate outer hinge point for axis
                %first take last panel id, take this spanwise position and
                %find first panel id in panelids which has this spanwise
                %position
                idLast=find(panelids,1,'last');
                id2=find(and(panelids,spanWisePointsI==spanWisePointsI(idLast)),1);
                panelPoints2=obj.grid(:,obj.panels(:,id2));
                %find fraction of span of this panel at which the control
                %surface leading edge ends
                fractChord2=interp1([chordWisePointsF(id2) chordWisePointsR(id2)],[0 1],obj.parCS{iParCS}.posChordF);
                fractSpan2=interp1([spanWisePointsI(id2) spanWisePointsO(id2)],[0 1],obj.parCS{iParCS}.posSpanOB);
                outerHingePoint=(panelPoints2(:,1)*(1-fractChord2)+panelPoints2(:,4)*fractChord2)*(1-fractSpan2)+(panelPoints2(:,2)*(1-fractChord)+panelPoints2(:,3)*fractChord)*(fractSpan2);
                
                %calculate hinge direction
                obj.parCS{iParCS}.panelIds=find(panelids)+obj.panel_start_idx-1;
                obj.parCS{iParCS}.hingeAxis=(outerHingePoint-innerHingePoint)/norm(outerHingePoint-innerHingePoint);
                obj.parCS{iParCS}.hingePoint=innerHingePoint;
            end
            
            if obj.symmetric
                
                
                for iParCS=1:length(obj.parCS)
                    obj.parCS{iParCS}.panelIdsL=obj.parCS{iParCS}.panelIds+size(obj.panels,2)/2;
                end
                
                
                
                
                for iSeg=1:length(obj.wing_segments)
                    
                    if obj.wing_segments(iSeg).has_te_cs
                        
                        obj.wing_segments(iSeg).te_device.panelIdsL = obj.wing_segments(iSeg).te_device.panelIds+size(obj.panels,2)/2;
                        obj.wing_segments(iSeg).te_device.panelIds_standard_L = obj.wing_segments(iSeg).te_device.panelIds_standard+size(obj.panels,2)/2;
                        
                        if ~isempty(obj.wing_segments(iSeg).te_device.panelIds_special)
                            obj.wing_segments(iSeg).te_device.panelIds_special_L = obj.wing_segments(iSeg).te_device.panelIds_special+size(obj.panels,2)/2;
                        end
                        
                    end
                    
                    if obj.wing_segments(iSeg).has_le_cs
                        obj.wing_segments(iSeg).le_device.panelIdsL=obj.wing_segments(iSeg).le_device.panelIds+size(obj.panels,2)/2;
                    end
                    
                    if obj.wing_segments(iSeg).has_sp_cs
                        
                        obj.wing_segments(iSeg).sp_device.panelIdsL=obj.wing_segments(iSeg).sp_device.panelIds+size(obj.panels,2)/2;
                        obj.wing_segments(iSeg).sp_device.panelIds_standard_L=obj.wing_segments(iSeg).sp_device.panelIds_standard+size(obj.panels,2)/2;
                        
                        if ~isempty(obj.wing_segments(iSeg).sp_device.panelIds_special)
                            obj.wing_segments(iSeg).sp_device.panelIds_special_L = obj.wing_segments(iSeg).sp_device.panelIds_special+size(obj.panels,2)/2;
                        end
                    end
                    
                end
            end
        end
        

        
        % this function is called whenever split_segment is used. It saves
        % the changes in geometry in a way that makes it possible to
        % related the current geometry (and its segment reference system)
        % with the original one that was described in CPACS.
        function obj = save_geometry_change(obj,segment_idx,f_1,f_2, split_type)
            
            if ~(isempty(obj.geom_arr))
            
            % determines what type of split segment was called
                switch split_type

                    % split_segment type resulting in the parent segment being
                    % split into 2 children segments (s1 and s2 or s2 and s3). Techinically
                    % s1 should be to the left of the original segment, and s3
                    % to the right, but this is not entirely certain.
                    case 1

                        if f_2 == 1

                            sub_eta_1 = f_1 * obj.geom_arr(segment_idx);
                            sub_eta_2 = obj.geom_arr(segment_idx) - sub_eta_1;

                            obj.geom_arr = [obj.geom_arr(1:segment_idx-1), sub_eta_1, sub_eta_2,obj.geom_arr(segment_idx+1:end)];

                        elseif f_1 == 0
                            sub_eta_3 = f_2 * obj.geom_arr(segment_idx);
                            sub_eta_2 = obj.geom_arr(segment_idx) - sub_eta_3;
                            obj.geom_arr = [obj.geom_arr(1:segment_idx-1), sub_eta_2, sub_eta_3,obj.geom_arr(segment_idx+1:end)];

                        else
                            msg = 'split_segment method behavior is unexpected';
                            error(msg);
                        end
                    % in both instances of geom_arr above, the etas are saved
                    % in such a way that they still represent the eta
                    % coordinate of the cut taken from the beginning of the
                    % original CPACS segment which was being cut. The reason
                    % behind this will become clear in the gen_cell_arrays
                    % function


                    % split_segment type resulting in the parent segment being
                    % split into 3 children segments (s1, s2 and s3). Techinically
                    % s1 should be to the left of the original segment, and s3
                    % to the right, but this is not entirely certain.
                    case 2

                        sub_eta_1 = f_1 * obj.geom_arr(segment_idx);
                        sub_eta_3 = f_2 * obj.geom_arr(segment_idx);
                        sub_eta_2 = obj.geom_arr(segment_idx) - sub_eta_1 - sub_eta_3;

                        obj.geom_arr = [obj.geom_arr(1:segment_idx-1), sub_eta_1, sub_eta_2, sub_eta_3,obj.geom_arr(segment_idx+1:end)];
                    % geom_arr within this case follows the same logic as in
                    % case 1.
                end
            end
        end
        
        % This function takes the geom_arr of the aerosurface (generated
        % through save_geometry_change, and actually converts it in a
        % format which can be easily used to relate the new wing
        % segmentation with the original CPACS one.
        function [obj, geom_cell, etas_cell] = gen_cell_arrays(obj)
            
            % The section below takes the geom array and re-arranges it
            % into an easier to read cell array. Each element of the cell
            % array is one of the original CPACS wing segments. The 
            geom_cell = cell(0,0); %will contain the geometry of the wing segment cuts
            etas_cell = cell(0,0); %will contain the etas of the wing segments in the new coord sys

            seg_count = 1;
            
            % the while loop does the following: it begins scanning the
            % geometry array from the start. Each iteration, the index
            % advances by one. Once the elements within seg_count and
            % loc_count sum up to 1, it means that all the segments within
            % those 2 indices used to be part of a single CPACS segment.
            % This is due to the way that the geom array is created.
            while seg_count < length(obj.geom_arr) + 1
                
                loc_sum = 0;
                loc_count = seg_count - 1;
                
                while loc_sum<1
                    
                    loc_count = loc_count + 1;
                    loc_sum = sum(obj.geom_arr(seg_count:loc_count));                   
                end
                
                geom_cell{end+1} = obj.geom_arr(seg_count:loc_count);
                etas_cell{end+1} = cumsum(obj.geom_arr(seg_count:loc_count));
                
                seg_count = loc_count + 1;
                
                % the geom_cell output is simply the cell version of the geom arr
                % the etas_cell is a nested (2D) cell array. The first dimention
                % contains the "parent CPACS" segments. Within each parent
                % segment, the eta coordinates of the cuts necessary to form
                % the new additional segments can be found. The coordinates
                % within each parent elements must of course be within 0 and 1.
                
                if sum(geom_cell{end}) ~= 1.
                    msg = 'The sum of the elements within each cell is not 1';
                    error(msg);
                end
            end
        end


        % looks at the location of all control surfaces
        % (belonging to the same component_segment) on this
        % aerosurface, then calculates exactly where the wing
        % will have to be cut in order to generate all the
        % segments necessary to perfectly contain all control
        % surfaces
        function obj = calculate_final_segmentation(obj)
            
            % generates the geometry and eta cells of the wing segmentation
            % AFTER that extra segments have been added due to structure,
            % but BEFORE that extra segments are added due to Control
            % Surfaces.
            [obj, obj.curr_geom_cell, obj.curr_etas_cell] = gen_cell_arrays(obj);
            
            % initializes "final" version of eats cells (which will contain
            % the segmentation which accounts for control surfaces)
            obj.final_etas_cell = obj.curr_etas_cell;
            
            % initializes "cuts" version of etas_cell. The reference system
            % is the same as etas_cell, but the etas stored in etas_cuts
            % are only the ones which will be later fed to split_segment in
            % order to actually split the segments
            obj.etas_cuts = cell(1,length(obj.final_etas_cell));

            % loops over all control surfaces associated with this
            % aerosurface object.
            for CS_counter = 1:length(obj.control_surfaces)

                CS_current = obj.control_surfaces(CS_counter);

                % saves the "start" eta and "end" eta of this control
                % surface instance. The eta coords are in the reference
                % frames belonging respectively to the "start" and "end"
                % segments. The aforementioned segments are the CPACS
                % segments, NOT THE CURRENT ONES DEFINED IN DAEDALUS
                eta_start = CS_current.start_segmentEta;
                eta_end = CS_current.end_segmentEta;
                
                % same as above, but regarding the start and end segments
                % indices. Again, these are CPACS INDICES, NOT DAEDALUS
                start_seg_idx = CS_current.start_segment_index;
                end_seg_idx = CS_current.end_segment_index;
                
                % adding eta coordinate of cuts due to start and end of
                % control surfaces to their respective segments.
                obj.etas_cuts{start_seg_idx}(end+1) = eta_start;
                obj.etas_cuts{end_seg_idx}(end+1) = eta_end;
                
                % all the if statements below are necessary in case no cuts
                % are required within the parent CPACS segment, in which
                % case the output needs to be specified in a certain way,
                % so that later code knows that that parent segment is to
                % be left alone.
                if obj.final_etas_cell{start_seg_idx} == 1
                    obj.final_etas_cell{start_seg_idx} = [eta_start,1];
                else
                    obj.final_etas_cell{start_seg_idx}(end+1) = eta_start;
                end
                
                if obj.final_etas_cell{end_seg_idx} == 1
                    obj.final_etas_cell{end_seg_idx} = [eta_end,1];
                    
                elseif isequal(obj.final_etas_cell{end_seg_idx},[0,1])==1 && eta_end == 1
                    
                else
                    obj.final_etas_cell{end_seg_idx}(end+1) = eta_end;
                end
            end
            
            % within this for loops, the elements of the final_etas_cell
            % and etas_cuts arrays are re-ordered and if necessary removed
            for cell_counter = 1:length(obj.final_etas_cell)
                
                % rounds elements of both cells to 4 decimal digits. This
                % allows for comparison later on, in order to determine
                % which control surface belongs to which segment
                obj.final_etas_cell{cell_counter} = round(obj.final_etas_cell{cell_counter},4);
                obj.etas_cuts{cell_counter} = round(obj.etas_cuts{cell_counter},4);
                
                % arranges etas within each parent segment in ascending order
                obj.final_etas_cell{cell_counter} = sort(obj.final_etas_cell{cell_counter},'ascend');
                
                % Removes etas within same parent segment which are not
                % unique (THIS OCCURS AFTER ROUNDING, MEANING THAT
                % SEGMENTS LESS THAN 1e-4 APART WILL BE MERGED INTO 1.
                % POSSIBLE LOSS OF DATA.
                obj.final_etas_cell{cell_counter} = unique(obj.final_etas_cell{cell_counter});
                
                % etas_cuts might contain some etas which are the same as
                % in curr_etas_cell. Since in etas_cuts we only want the
                % etas where we have to apply a cut, the block below first
                % removes all elements from etas_cuts which are in common
                % with curr_etas_cell (within the same CPACS parent segment
                % of course). Then the elements within etas_cuts are
                % soreted in ascending order. Finally, all repeated
                % elements are removed (in case 2 control surfaces have
                % extremeties on the same eta).
                obj.etas_cuts{cell_counter} = setdiff(obj.etas_cuts{cell_counter}, obj.curr_etas_cell{cell_counter});
                obj.etas_cuts{cell_counter} = sort(obj.etas_cuts{cell_counter},'ascend');
                obj.etas_cuts{cell_counter} = unique(obj.etas_cuts{cell_counter});
            end
        end
        
        % looks at the output of calculate_final_segmentation,
        % and actually does all the splitting, resulting in
        % wing_segments being created through split_segment.
        % All of them are of course associated to this specific
        % instance (obj) of aerosurface.
        
        % get ready for a headache, the loops below are pretty nested.
        function obj = wing_segmentation(obj)
            
            idx_seg = length(obj.final_etas_cell); 
            
            %idx_seg represents the CPACS parent segment index that we are
            %currently looking at
            while idx_seg>0
                
                % the if statements below are necessary in case a parent
                % segment is not supposed to be cut
                if length(obj.final_etas_cell{idx_seg}) == 2
                    
                    if obj.final_etas_cell{idx_seg}(1) == 0 && obj.final_etas_cell{idx_seg}(2) == 1
                        cnd = 1;
                    else
                        cnd = 2;
                    end
                    
                elseif obj.final_etas_cell{idx_seg} == 1
                    cnd = [0,0];
                    
                else
                    cnd = [1,1];
                end
                
                % the if statement below runs if the parent segment is
                % supposed to be cut
                if  sum(cnd)== 2
                    
                    % idx2_cut is the segment index (WITHIN THE CPACS
                    % PARENT) of the cut that we want to apply
                    for idx2_cut = length(obj.etas_cuts{idx_seg}):-1:1

                        % idx2_curr is the segment index (WITHIN THE CPACS
                        % PARENT) of the eta of an already existing segment
                        % within the CPACS parent
                        idx2_curr = length(obj.curr_etas_cell{idx_seg});
                        
                        % the loops below shows the logic behind this function. Basically the idx_seg and idx2_cut are
                        % fixed, so we are looking at the eta of a specific cut that we want to apply to a specific (parent)
                        % segment. The idx2_curr is decreased, until the eta of an already existing segment border is to
                        % the left (smaller) of the cut that we want to apply. This means that we pick a cut within a
                        % parent segment, then we swipe, from right to left, the already existing children segment
                        % edges, in order to find a segmernt TO THE LEFT OF THE CUT WE WANT TO CREATE. Once we find one,
                        % split_segment is used on that segment, thus generaing 2 segments (hence 1 new one). The "new"
                        % segment will be to the right of the "original one". The advantage of this is that we do not
                        % have to keep track within the index system of the new segments we are creating, as they will all be
                        % to the RIGHT of the ones we are interested in, meaning that the indices that we will work with
                        % in the furure are untouched.
                        while obj.etas_cuts{idx_seg}(idx2_cut) < obj.curr_etas_cell{idx_seg}(idx2_curr)
                            idx2_curr = idx2_curr - 1;
                        
                        
                            % In case there is no pre-existing segment edge
                            % to the left of the new cut, we have to split
                            % the segment to the right, thus affecting the
                            % indices of the segments before it. If that is
                            % the case, the if statement below is executed,
                            % the segments are split, and the change in
                            % indices is taken into account.
                            if idx2_curr == 0
                                
                                % whenever you see code with cellfun, it
                                % means that the segment index IN THE
                                % GLOBAL DAEDALUS frame of reference is
                                % being calculated. This is found by going
                                % through each CPACS parent before the
                                % current segment, summing the segments
                                % they contain, and finally adding the
                                % segments within the current CPACS parent,
                                % but preceding (to the left) the specific
                                % segment we are looking at
                                seg_idx_cut = sum(cellfun(@(x) numel(x),obj.curr_etas_cell(1:idx_seg)))...
                                - sum(cellfun(@(x) numel(x),obj.curr_etas_cell(idx_seg))) + 1;
                                
                                % does the cutting
                                f_2 = (obj.etas_cuts{idx_seg}(idx2_cut))/ (obj.curr_etas_cell{idx_seg}(1));
                                obj = obj.split_segment(seg_idx_cut,0,f_2,0,0,0);
                                
                                % accounts for segment being added to the
                                % left and indices changing
                                obj.curr_etas_cell{idx_seg} = [obj.etas_cuts{idx_seg}(idx2_cut),obj.curr_etas_cell{idx_seg}];
                                break;
                            end 
                        end
                    
                        
                        % this if statement runs if it is possible to
                        % create a new segment to the right of a
                        % pre-existing one
                        if idx2_curr ~= 0
                            
                            % see comment above for similar code (using
                            % cellfun)
                            seg_idx_cut = sum(cellfun(@(x) numel(x),obj.curr_etas_cell(1:idx_seg)))...
                                - sum(cellfun(@(x) numel(x),obj.curr_etas_cell(idx_seg)))...
                                + idx2_curr;


                            % Always split creating segments on the outer section
                            f_2 = (obj.etas_cuts{idx_seg}(idx2_cut) - obj.curr_etas_cell{idx_seg}(idx2_curr)) / (1-obj.curr_etas_cell{idx_seg}(idx2_curr));
                            obj = obj.split_segment(seg_idx_cut,0,f_2,0,0,0);
                            
                            obj.curr_etas_cell{idx_seg} = [obj.curr_etas_cell{idx_seg}(1:idx2_curr), obj.etas_cuts{idx_seg}(idx2_cut),obj.curr_etas_cell{idx_seg}(idx2_curr+1:end)];
                            
                        end
                    end
                end
                
                % once we looked at all cuts within a CPACS parent segment,
                % we can move to the next one (moving from right to left,
                % so from tip to root).
                idx_seg = idx_seg - 1;

            end
        end
    
        % looks at the output of wing_segmentation (so the
        % current wing_segment disposition) and assigns each
        % control surface to each spefic wing_segment
        function obj = assign_control_surfaces(obj)
            
            % the obj.final_etas_cell cell array is both used to calculate
            % segment indices and to associate the starting and ending
            % points of control surfaces on the wing. In special cases (when a
            % control surface is perfectly contained within a parent segment) the 
            % final_etas_cell of that parent segment is written as [0,1].
            % This should be read as 1 segment, but is instead read as 2. 
            
            % This knockdown_arr
            % fixes that. Each parent segment has a certain knockdown
            % value, which needs to be applied to all indices of all
            % wingsegments belonging to that parent.
            knockdown_arr = zeros(1,length(obj.final_etas_cell));
            
            for knock_index = 1:length(obj.final_etas_cell)
                
                if isequal(obj.final_etas_cell{knock_index}, [0,1])
                    
                    knockdown_arr(knock_index) = 1;
                end
            end
            
            knockdown_arr = cumsum(knockdown_arr) - knockdown_arr;
            
            % loops through all control surfaces belonging to this
            % aerosurface object
            for CS_counter = 1:length(obj.control_surfaces)

                CS_current = obj.control_surfaces(CS_counter);
                
                % the loop below takes the starting eta of a control
                % surface (within the parent CPACS segment frame of
                % reference) and compares it with the eta of the final
                % segmentation of the wing, within the same CPACS parent
                % elements, and within of course the same reference frame.
                % Once the same eta coordinate is found, the inde of the
                % start child segment within its parent is found (for that
                % specific control surface)
                start_sub_idx = 1;
                while round(obj.final_etas_cell{CS_current.start_segment_index}(start_sub_idx),4) ~= round(CS_current.start_segmentEta,4)
                    start_sub_idx = start_sub_idx + 1;
                end
                
                if round(CS_current.start_segmentEta,4) == 1.0
                    start_sub_idx = start_sub_idx + 1;
                end
                
                % the loop below does the same as the one above, only for
                % the end of the control surface (which means that WE MAY
                % OR MAY NOT be looking within the same CPACS parent
                % segment. It depends on the CS_current.end_segment_index).
                
                end_sub_idx = 1;
                while round(obj.final_etas_cell{CS_current.end_segment_index}(end_sub_idx),4) ~= round(CS_current.end_segmentEta,4)
                    end_sub_idx = end_sub_idx + 1;
                end
                
                % calculates the start segment index within the GLOBAL
                % DAEDALUS frame of reference. Functioning of the cellfun
                % code below is explained above similar code within
                % wing_segmentation method.
                final_start_seg_idx = sum(cellfun(@(x) numel(x),obj.final_etas_cell(1:CS_current.start_segment_index)))...
                - sum(cellfun(@(x) numel(x),obj.final_etas_cell(CS_current.start_segment_index)))...
                + start_sub_idx + 1;
                
                % same as above, but for the end segment index in the GLOBAL
                % DAEDALUS frame of reference.
                final_end_seg_idx = sum(cellfun(@(x) numel(x),obj.final_etas_cell(1:CS_current.end_segment_index)))...
                - sum(cellfun(@(x) numel(x),obj.final_etas_cell(CS_current.end_segment_index)))...
                + end_sub_idx ;
            
                % applying index knockdown
                final_start_seg_idx = final_start_seg_idx - knockdown_arr(CS_current.start_segment_index);
                final_end_seg_idx = final_end_seg_idx - knockdown_arr(CS_current.end_segment_index);
                
                % makes sure that right index is saved in case eta_end is
                % right on the edge between the current and next segment,
                % and the segment was never split (original parent).
                if obj.final_etas_cell{CS_current.end_segment_index}(end_sub_idx) == 1 && isequal(obj.final_etas_cell{CS_current.end_segment_index}, [0,1])
                    final_end_seg_idx = final_end_seg_idx - 1;
                end
                
                % necessary to avoid error due to exceeding matrix index
                % when running elseif statement below this one while the CS
                % "starts" at eta=1
                if round(CS_current.start_segmentEta,4) == 1
                    
                % makes sure that right index is saved in case eta_end is
                % right on the edge between the current and next segment,
                % and the segment was never split (original parent).
                elseif obj.final_etas_cell{CS_current.start_segment_index}(start_sub_idx) == 0 && isequal(obj.final_etas_cell{CS_current.start_segment_index}, [0,1])
                    final_start_seg_idx = final_start_seg_idx - 1;                    
                end
                
                % adds current control surface to its start wing segment
                nr_segments_covered = abs(final_end_seg_idx - final_start_seg_idx)+1;
                obj.control_surfaces(CS_counter).children = cell(1,nr_segments_covered);
                obj.control_surfaces(CS_counter).children{1,1} = class_control_surface(CS_current,1,final_start_seg_idx);
                obj.wing_segments(final_start_seg_idx) = obj.wing_segments(final_start_seg_idx).add_control_surface(obj.control_surfaces(CS_counter).children{1,1});
                
                
                % if the difference beetween the end and start wing segment
                % indices is more than 1, it means that there are additonal
                % wing segments present between them. As the control
                % surface begins at the start wing segment and finishes at
                % the end wing segment, it must necessarily span all wing
                % segments in between. The loop below applies the control
                % surface to the "in-between" segments.
                if abs(final_end_seg_idx - final_start_seg_idx) > 1
                    
                    for mid_idx = final_start_seg_idx+1 : final_end_seg_idx-1
                        obj.control_surfaces(CS_counter).children{1,mid_idx+1-final_start_seg_idx} = class_control_surface(CS_current,mid_idx+1-final_start_seg_idx,mid_idx);
                        obj.wing_segments(mid_idx) = obj.wing_segments(mid_idx).add_control_surface(obj.control_surfaces(CS_counter).children{1,mid_idx+1-final_start_seg_idx});
                    end
                end
                
                
                if abs(final_end_seg_idx - final_start_seg_idx) > 0
                    % adds current control surface to its end wing segment
                    obj.control_surfaces(CS_counter).children{1,end} = class_control_surface(CS_current,nr_segments_covered,final_end_seg_idx);
                    obj.wing_segments(final_end_seg_idx) = obj.wing_segments(final_end_seg_idx).add_control_surface(obj.control_surfaces(CS_counter).children{1,end});
                end
            end
        end
    end
end
