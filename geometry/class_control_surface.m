%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%
classdef class_control_surface < handle% &  matlab.mixin.SetGetExactNames
    %UNTITLED Summary of this class goes here
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
        panelIds;
        %> Aerodynamic Panel Ids of Left Wing (so symmetric)
        panelIdsL;
        % Aerodynamic Panel Ids of "special regions" (ig overlap between
        % flap and spoiler)
        panelIds_special;
        % Same as above, but for Left Wing (symmetric)
        panelIds_special_L;
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
        
        % index of CS instance within parent instance
        child_ixd;
        % index wing-segment where instance is contained
        seg_idx;
    end
    
    methods
        function obj=class_control_surface(cs_parent, child_ixd, seg_idx)
            if nargin == 0
                %pass
            elseif nargin>0
                obj.name=[cs_parent.name, '_part_',num2str(child_ixd)];
                obj.child_ixd = child_ixd;
                obj.seg_idx = seg_idx;
                obj.pos=cs_parent.pos;
                obj.hinge=cs_parent.hinge;
                obj.n_hinge=cs_parent.n_hinge;
                obj.delta=cs_parent.delta;

                %spoilers always have 2 extra varargin compared to
                %other CS, due to presence of hinge_start. Hence the
                %TE and LE equivalent of nargin==4 runs in nargin==5
                %for spoilers, and so on.
                obj.is_sym = cs_parent.is_sym;
                obj.is_sym_defl = cs_parent.is_sym_defl;
                obj.delta_l_r = cs_parent.delta_l_r;
                obj.is_tapered = cs_parent.is_tapered;
                    
                if obj.pos == 2
                    obj.xsi_LE_inner = cs_parent.xsi_LE_inner;
                    obj.xsi_LE_outer = cs_parent.xsi_LE_outer;
                end
            end
        end
        
        
    end
end
