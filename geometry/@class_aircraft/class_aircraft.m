%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%
classdef class_aircraft
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        %> aircraft name    
        name;
        
        %> fuselage object(s)
        fuselages;
        fuselage_structure_properties;
        %> wing object(s)
        wings;
        wings_structural_properties;
        %> additional settings
        addSettings;
        
        nacelles;
        pylons;
        
        %
        engines;
        
        aero_surfaces;
        
        control_surfaces;
        control_deflections;
        controlSurfaceMoments;

        S_wet_fus;
        S_wet_wing;
        S_wet_stabilizer;
        S_wet_engine;
        
        W_TO;       % Gross Take-off weight
        OWE;
        W_E;        % Empty weight
        W_F;        % Mission fuel weight
        P_TO;       % Take-off Power
        
        CL_max_L;
        CL_max_TO;
        CL_max;
        
        CD;
        CD_f;
        grid;
        grid_flat;
        grid_deflected;
        
        panels;
        te_idx;
        
        % volume grid
        grid_3D;  %flag if 3D grid required
        grid_vol;
        panels_vol;
        is_te;
        is_te_vol;
        te_idx_vol;
        opposite_te_vol;
        
        grid_wake;
        panels_wake;
        
        reference;

        weights;
        panel_to_beam_element;
        grid_settings
        J;
        
        control;
        acc;
        boundary_conditions;
        %> clamping conditions (restraining) from xml
        clamping_conditions;        
        %{nGustZones x 1} cell array with panelIds belonging to GustZones
        gustZones;
        
        % array containing handles to parent control surfaces (originals
        % are stored in their respective class_aerosurface object)
        control_surfaces_parents;
        
        %maximum operating altitude
        Zmo=13500;
        
        % struct containing the points for cheb type distributions 
        dataForCheb
    end
    
    methods (Static)
        
        function obj = create_from_cpacs(filename)
            obj = class_aircraft();
            %% check if tixi is installed
            if or(~exist('tixiOpenDocument'), ~exist('tigl_matlab'))
                if ~exist('tixiOpenDocument')
                    disp('Tixi v2.2.4 is required. https://github.com/DLR-SC/tixi/releases/download/v2.2.4/TIXI-2.2.4-win64.exe')
                end
                if ~exist('tigl_matlab')
                    disp('Tigl v2.2.0 is required. https://github.com/DLR-SC/tigl/releases/download/v2.2.0/TIGL-2.2.0-win64.exe')
                end
                return
            end
            tixi = tixiOpenDocument(filename);
     
            x_ac = '/cpacs/vehicles/aircraft';
            if tixiGetNamedChildrenCount(tixi, '/cpacs/vehicles', 'aircraft') > 1
                x_ac = [a_ac, '[1]'];
            end

            x_mod = [x_ac, '/model'];
            if tixiGetNamedChildrenCount(tixi, x_ac, 'model') > 1
                x_mod = [x_mod, '[1]'];
            end

            try
                obj.name = tixiGetTextElement(tixi, [x_mod, '/name']);
            catch
                obj.name = 'unspeficied';
            end

            ac_uID = tixiGetTextAttribute(tixi, x_mod, 'uID');
            tigl = tiglOpenCPACSConfiguration(tixi, ac_uID);
            
            x_ref = [x_mod, '/reference'];
            obj.reference.S_ref = tixiGetDoubleElement(tixi, [x_ref, '/area']);
            obj.reference.c_ref = tixiGetDoubleElement(tixi, [x_ref, '/length']);
            
            [x,y,z] = tixiGetPoint(tixi, [x_ref, '/point']);
            obj.reference.p_ref = [x,y,z];

            n_wings = tiglGetWingCount(tigl);
            if n_wings == 1
                [wing, b_ref] = ...
                    class_aerosurface.create_from_cpacs(tixi, tigl, 1);
                obj.wings = wing;
                obj.wings_structural_properties = obj.wings.wings_structural_properties;
                
                obj.reference.b_ref = b_ref;
            else
                for i = n_wings:-1:1
                    [wing, b_ref] = ...
                        class_aerosurface.create_from_cpacs(tixi, tigl, i);
                    wings_structural_properties(i) = wing.wings_structural_properties;
                    wings(i) = wing;
                    
                    if i == 1
                        obj.reference.b_ref = b_ref;
                    end
                end
                
                obj.wings = wings;
                obj.wings_structural_properties = wings_structural_properties;
            end
            
            obj.weights = class_weights;
            obj.weights.WingSystemsEstimate = [];
            obj.weights.WingSkinEstimate = [];
            obj.weights.FuselageNonStructuralEstimate = [];
            obj.weights.FuselageSystemsEstimate = [];
            obj.weights.WingInitialGuess = [];
            
            if tixiGetNamedChildrenCount(tixi, x_mod, 'analyses')
                if tixiGetNamedChildrenCount(tixi, [x_mod, '/analyses'], 'massBreakdown')
                    x_mbd = [x_mod, '/analyses/massBreakdown'];
                    x_desMass = [x_mbd, '/designMasses'];
                    
                    obj.weights.MTOW = tixiGetDoubleElement(tixi, [x_desMass, '/mTOM/mass']);
                    obj.weights.MZFW = tixiGetDoubleElement(tixi, [x_desMass, '/mZFM/mass']);
                    obj.weights.MLW = tixiGetDoubleElement(tixi, [x_desMass, '/mMLM/mass']);
                    obj.weights.OWE = tixiGetDoubleElement(tixi, [x_mbd, '/mOEM/massDescription/mass']);
                    
                    f_oew_systs = 0.3;
                    f_wingbox_skin = 0.3;
                    f_wingstruct_skin = 0.2;
                    
                    if tixiGetNamedChildrenCount(tixi, [x_mbd, '/mOEM'], 'mEM')
                        if tixiGetNamedChildrenCount(tixi, [x_mbd, '/mOEM/mEM'], 'mSystems')
                            W_systs = tixiGetDoubleElement(tixi, [x_mbd, '/mOEM/mEM/mSystems/massDescription/mass']);
                        else
                            W_systs = f_oew_systs*obj.weights.OWE;
                        end
                        
                        if tixiGetNamedChildrenCount(tixi, [x_mbd, '/mOEM/mEM'], 'mStructure')
                            W_struct = tixiGetDoubleElement(tixi, [x_mbd, '/mOEM/mEM/mStructure/massDescription/mass']);
                            
                            if tixiGetNamedChildrenCount(tixi, [x_mbd, '/mOEM/mEM/mStructure'], 'mWingsStructure')
                                n_mWingStructure = tixiGetNamedChildrenCount(tixi, [x_mbd, '/mOEM/mEM/mStructure/mWingsStructure'], 'mWingStructure');
                                if n_mWingStructure ~= n_wings
                                    error('Number of mWingStructures does not match number of wings!');             
                                else
                                    for i = 1:n_mWingStructure
                                        x_mWingStructure = sprintf('%s/mOEM/mEM/mStructure/mWingsStructure/mWingStructure[%i]', x_mbd, i);
                                        W_wingstruct = tixiGetDoubleElement(tixi, [x_mWingStructure, '/massDescription/mass']);
                                        
                                        WingSkinEstimate = 0;
                                        n_mComponentSegment = tixiGetNamedChildrenCount(tixi, x_mWingStructure, 'mComponentSegment');
                                        for j = 1:n_mComponentSegment
                                            x_mComponentSegment = sprintf('%s/mComponentSegment[%i]', x_mWingStructure, j);
                                            if tixiGetNamedChildrenCount(tixi, x_mComponentSegment, 'mWingBox')
                                                if tixiGetNamedChildrenCount(tixi, [x_mComponentSegment, '/mWingBox'], 'mSkins')
                                                    WingSkinEstimate = WingSkinEstimate + tixiGetDoubleElement(tixi, [x_mComponentSegment, '/mWingBox/mSkins/massDescription/mass']);
                                                else
                                                    WingSkinEstimate = WingSkinEstimate + f_wingbox_skin*tixiGetDoubleElement(tixi, [x_mComponentSegment, '/mWingBox/massDescription/mass']);
                                                end
                                            else
                                                WingSkinEstimate = WingSkinEstimate + f_wingstruct_skin*tixiGetDoubleElement(tixi, [x_mComponentSegment, '/massDescription/mass']);
                                            end
                                        end
                                        
                                        obj.weights.WingSkinEstimate(i) = WingSkinEstimate;
                                        obj.weights.WingInitialGuess(i) = tixiGetDoubleElement(tixi, [x_mWingStructure, '/massDescription/mass']);
                                        obj.weights.WingSystemsEstimate(i) = W_systs*W_wingstruct/W_struct;
                                    end
                                end
                            end
                        end
                    end
                end
            end
            
            % moving all parent control surfaces from wing objects to
            % aircraft object
            
            for counter = 1:length(obj.wings)
                
                for counter_cs = 1:length(obj.wings(counter).control_surfaces)
                    obj.wings(counter).control_surfaces(counter_cs).wing_index = counter;
                end
                
                obj.control_surfaces_parents(counter) = obj.wings(counter).control_surfaces;
            end
            
            obj.grid_settings=class_grid_settings;    
            mkdir('results', obj.name);
        end
        
    end
    
    methods
        %> function to read xml file
        obj=read_xml_definition(obj,filename);
        
        function obj = class_aircraft(wing,varargin)
            if nargin == 0
                % pass
            elseif nargin==1
                obj.CD_f=0;
                aero_surface(1)=wing;
                obj.wings=aero_surface;
                obj.name='blank';
            else
                obj=read_xml_definition(obj,wing);
                obj.runInputChecks();
            end
            
            if nargin ~= 0
                obj.grid_settings=class_grid_settings;
            
                mkdir('results',obj.name);
            end
        end
        
        function obj = add_aerosurface(obj,aero_surface)
            obj.aero_surfaces=[obj.aero_surfaces aero_surface];
        end
        
        function obj = add_wing(obj,wing)
            obj.wings=[obj.wings wing];
        end
        
        function obj = add_fuselage(obj,fuselage)
            obj.fuselages=[obj.fuselages fuselage];
        end
        
        function obj = add_nacelle(obj,nacelle)
            obj.nacelles=[obj.nacelles nacelle];
        end
        
        function obj=plot(obj)
            hold on
            for i=1:length(obj.aero_surfaces)
                obj.aero_surfaces(i).plot(); 
            end
            for i=1:length(obj.wings)
                obj.wings(i).plot();
            end

            for i=1:length(obj.fuselages)
                obj.fuselages(i).plot_fuselage_surface();
            end
            
            for i=1:length(obj.nacelles)
                obj.nacelles(i).plot_fuselage_surface();
            end
        end
        
        function obj = compute_CD_f(obj,state,S_ref)
            obj.CD_f=0;
            for i=1:length(obj.aero_surfaces)
                obj.aero_surfaces(i)=obj.aero_surfaces(i).compute_friction_drag(state,S_ref);
                obj.CD_f=obj.CD_f+obj.aero_surfaces(i).CD_f*(1+obj.wings(i).symmetric);
            end
            
            for i=1:length(obj.wings)
                obj.wings(i)=obj.wings(i).compute_friction_drag(state,S_ref);
                if obj.wings(i).symmetric==1
                    obj.CD_f=obj.CD_f+obj.wings(i).CD_f*2;
                else
                    obj.CD_f=obj.CD_f+obj.wings(i).CD_f;
                end
            end
            
            for i=1:length(obj.fuselages)
                obj.fuselages(i)=obj.fuselages(i).compute_wetted_area();
                obj.fuselages(i)=obj.fuselages(i).compute_friction_drag(state,S_ref);
                obj.CD_f=obj.CD_f+obj.fuselages(i).CD_f;
            end
            
            for i=1:length(obj.nacelles)
                obj.nacelles(i)=obj.nacelles(i).compute_wetted_area();
                obj.nacelles(i)=obj.nacelles(i).compute_friction_drag(state,S_ref);
                obj.CD_f=obj.CD_f+obj.nacelles(i).CD_f;
            end
        end
        
        function obj=compute_force_interpolation_matrix(obj,aircraft_structure)
            obj.panel_to_beam_element=zeros(size(obj.panels,2),15);
            for i=1:length(obj.wings)
                obj.panel_to_beam_element=obj.wings(i).compute_force_interpolation_matrix(obj.panel_to_beam_element);
            end
            
            for i=1:length(obj.wings)
                obj.wings(i)=obj.wings(i).T_matrix(obj.panel_to_beam_element,obj.grid,obj.panels,aircraft_structure.beam(i));
            end
        end

        function obj=compute_beam_forces(obj,F_body, M_body ,aircraft_structure,varargin)
            if nargin==5
                nw=varargin{1};
            elseif nargin==3
                aircraft_structure=M_body;
                M_body=F_body*0;
                nw=length(obj.wings);
            else
                nw=length(obj.wings);
            end
            for i=1:nw
                obj.wings(i)=obj.wings(i).compute_beam_forces(obj.panel_to_beam_element,F_body, M_body ,aircraft_structure.beam(i));
            end
        end
        
        function obj=compute_forces(obj,F_aero,panels,grid) 
            for i=1:length(obj.wings)
                obj.wings(i)=obj.wings(i).compute_c4_forces(F_aero,panels,grid);
            end
            
            for i=1:length(obj.aero_surfaces)
                obj.aero_surfaces(i)=obj.aero_surfaces(i).compute_c4_forces(F_aero,panels,grid);
            end
        end
        
        function obj=compute_deflected_grid(obj,deflections_structmesh)
            x_max=obj.grid_settings.y_max_grid_size;
            obj.grid_deflected=obj.grid;
            for i=1:length(obj.wings)
                obj.grid_deflected=obj.wings(i).compute_deflected_grid(obj.panels,obj.grid_deflected,deflections_structmesh(i).def);
                 %obj.grid_deflected=obj.wings(i).grid_deflected;
            end
            if obj.grid_settings.aerodynamic_fuselage==1
                for j=i+1:i+length(obj.fuselages)
                    obj.grid_deflected=obj.fuselages(j-i).compute_deflected_grid_flat(obj.panels,obj.grid_deflected,deflections_structmesh(j).def,obj.fuselages(j-i).center_coords(3,1),x_max);
                end
                for k=j+1:j+length(obj.nacelles)
                    obj.grid_deflected=obj.nacelles(k-j).compute_deflected_grid_flat(obj.panels,obj.grid_deflected,deflections_structmesh(j+length(obj.nacelles)-(k-j)+1).def,obj.nacelles(k-j).center_coords(3,1),x_max);
                end
            end
        end
        
        function obj=compute_wingbox_coords(obj)
            if ~isempty(obj.wings_structural_properties)
                if isfield(obj.wings_structural_properties, 'fs_segments')
                    for i = 1:length(obj.wings)
                        obj.wings(i) = obj.wings(i).compute_wingbox_coords(obj.grid_settings.dy_max_struct_grid, true);
                    end
                else                
                    for i=1:length(obj.wings)
                        obj.wings(i)=obj.wings(i).compute_wingbox_coords(obj.grid_settings.dy_max_struct_grid,obj.wings_structural_properties(i).frontspar,obj.wings_structural_properties(i).rearspar);
                    end
                end
            else
                for i=1:length(obj.wings)
                    obj.wings(i)=obj.wings(i).compute_wingbox_coords(obj.grid_settings.dy_max_struct_grid);
                end
            end
        end
        
        function obj=compute_shell_coords(obj)
            for i=1:length(obj.fuselages)
                obj.fuselages(i)=obj.fuselages(i).compute_shell_coords(obj.grid_settings.y_max_grid_size);
            end
        end
        
        function obj=compute_acceleration(obj,wingaero,wingstructure,state)
            
            CX=wingaero.CX;
            CZ=wingaero.CZ;
            CY=wingaero.CY;
            CL=wingaero.CL;
            CM=wingaero.CM;
            CN=wingaero.CN;
            rho=wingaero.state.rho_air;
            V=norm(wingaero.Uinf);
            Sref=wingaero.reference.S_ref;
            bref=wingaero.reference.b_ref;
            cref=wingaero.reference.c_ref;
            
            m=wingstructure.f_compute_totalMass;
            
            J=wingstructure.f_compute_moment_of_inertia(wingaero.reference.p_ref);
            J1=J(1,1);%state.aircraft_state.I_xyz(1,1);
            J2=J(2,2);%state.aircraft_state.I_xyz(2,2);
            J3=J(3,3);%state.aircraft_state.I_xyz(3,3);
            
            obj.J=[J1 J2 J3];
            ax=1/m*(1/2*rho*V^2*CX*Sref);
            ay=1/m*(1/2*rho*V^2*CY*Sref);
            az=1/m*(1/2*rho*V^2*CZ*Sref);
            ap=1/J1*(1/2*rho*V^2*CL*Sref*bref);
            aq=1/J2*(1/2*rho*V^2*CM*Sref*cref);
            ar=1/J3*(1/2*rho*V^2*CN*Sref*bref);
            
            obj.acc=[ax ay az ap aq ar]; 
        end
        
        
        function obj=compute_grid(obj)
            % read desired grid size from grid_settings structure
            x_max=obj.grid_settings.x_max_grid_size;
            y_max=obj.grid_settings.y_max_grid_size;
            n_x_min=obj.grid_settings.n_x_min;
            wake=obj.grid_settings.wake;

            % initialize variables for 2D grid
            grid=[];
            grid_flat=[];
            te_idx=[];
            panels=[];
            
            % initialize variables for 3D grid
            obj.grid_vol=[];
            obj.panels_vol=[];
            obj.is_te_vol=[];
            obj.is_te=[];
            obj.te_idx_vol=[];
            obj.opposite_te_vol=[];
            
            obj.grid_wake=[];
            obj.panels_wake=[]; 
            
            % assemble grid from all wings
            for i=1:length(obj.wings)
                % first compute all sub grids
                obj.wings(i)=obj.wings(i).compute_grid(x_max,y_max,n_x_min,wake);
                
                % all operations for 2D grid
                grid_len_b4=length(grid);
                grid=[grid obj.wings(i).grid];
                grid_flat=[grid_flat obj.wings(i).grid_flat];
                panel_len_b4=length(panels);
                panels=[panels,obj.wings(i).panels+grid_len_b4];
                te_idx=[te_idx,obj.wings(i).te_idx+grid_len_b4];
                obj.wings(i).grid_start_idx=grid_len_b4+1;
                obj.wings(i).panel_start_idx=panel_len_b4+1;
                
                obj.wings(i)=obj.wings(i).update_idx();

                obj.is_te=[obj.is_te obj.wings(i).is_te_vol(1:end/2)];
                % all operations for 3D grid
                grid_len_b4_vol=length(obj.grid_vol);
                panel_len_b4_vol=length(obj.panels_vol);
                obj.grid_vol=[obj.grid_vol obj.wings(i).grid_vol];
                obj.panels_vol=[obj.panels_vol obj.wings(i).panels_vol+grid_len_b4_vol];
                obj.is_te_vol=[obj.is_te_vol obj.wings(i).is_te_vol];
                obj.te_idx_vol=[obj.te_idx_vol obj.wings(i).te_idx_vol+grid_len_b4_vol];
                obj.opposite_te_vol=[obj.opposite_te_vol obj.wings(i).opposite_te_vol+panel_len_b4_vol];
                obj.wings(i).grid_start_idx_vol=grid_len_b4_vol+1;
                obj.wings(i).panel_start_idx_vol=panel_len_b4_vol+1;
                
                if (wake==1)||(wake==2)
                    grid_len_wake_b4_vol=length(obj.grid_wake);
                    obj.grid_wake=[obj.grid_wake obj.wings(i).grid_wake];
                    obj.panels_wake=[obj.panels_wake obj.wings(i).panels_wake+grid_len_wake_b4_vol];
                end
            end
            
            obj.grid=grid;
            obj.grid_flat=grid_flat;
            obj.panels=panels;
            obj.te_idx=te_idx;
            if obj.grid_settings.aerodynamic_fuselage==1
                obj=compute_fuselage_grid(obj);
            end
            obj.gustZones={1:size(obj.panels,2)};
            
            obj = obj.computeControlSurfacePanelIds();
        end    
        

        function obj=compute_fuselage_grid(obj)
            
            % read desired grid size from grid_settings structure
            x_max=obj.grid_settings.y_max_grid_size;
            y_max=obj.grid_settings.x_max_grid_size;
            wake=obj.grid_settings.wake;
            
            % initialize variables for 2D grid
           % grid=obj.grid_vol;
            grid_flat=obj.grid;
            te_idx=obj.te_idx_vol;
           % panels=obj.panels_vol;
            te_idx_flat=obj.te_idx;
            panels_flat=obj.panels;
            len=[];
            
            for i=1:length(obj.fuselages)
                
                obj.fuselages(i)=obj.fuselages(i).compute_shell_coords(x_max);
                obj.fuselages(i)=obj.fuselages(i).compute_grid();
                obj.fuselages(i)=obj.fuselages(i).compute_grid_flat(obj.fuselages(i).center_coords(3,1),x_max);
                
                % all operations for 2D grid
                grid_len_b4=length(grid_flat);
                panel_len_b4=length(panels_flat);
                grid_flat=[grid_flat obj.fuselages(i).grid_flat];
                panels_flat=[panels_flat,obj.fuselages(i).panels_flat+grid_len_b4];
                te_idx_flat=[te_idx_flat,obj.fuselages(i).te_idx_flat+grid_len_b4];
                obj.fuselages(i).grid_start_idx_flat=grid_len_b4+1;
                obj.fuselages(i).panel_start_idx_flat=panel_len_b4+1; 
                
                
                obj.fuselages(i)=obj.fuselages(i).update_idx();
%                 grid_len_b4=length(grid);
%                 grid=[grid obj.fuselages(i).grid_flat];
%                 panel_len_b4=length(panels);
%                 panels=[panels,obj.fuselages(i).panels+grid_len_b4];
%                 te_idx=[te_idx,obj.fuselages(i).te_idx+grid_len_b4];
%                 obj.fuselages(i).grid_start_idx_vol=grid_len_b4+1;
%                 obj.fuselages(i).panel_start_idx=panel_len_b4+1; 
                %obj.wings(i)=obj.wings(i).update_idx();

                grid_len_b4_vol=length(obj.grid_vol);
                panel_len_b4_vol=length(obj.panels_vol);
                obj.grid_vol=[obj.grid_vol obj.fuselages(i).grid_vol];
                obj.panels_vol=[obj.panels_vol obj.fuselages(i).panels_vol+grid_len_b4_vol];
                obj.te_idx_vol=[obj.te_idx_vol obj.fuselages(i).te_idx_vol+grid_len_b4_vol];
                obj.is_te_vol=[obj.is_te_vol obj.fuselages(i).is_te_vol];
                obj.opposite_te_vol=[obj.opposite_te_vol obj.fuselages(i).opposite_te_vol+panel_len_b4_vol];
                obj.fuselages(i).grid_start_idx_vol=grid_len_b4_vol+1;
                obj.fuselages(i).panel_start_idx_vol=panel_len_b4_vol+1;
            end
             
            for i=1:length(obj.nacelles)
                obj.nacelles(i)=obj.nacelles(i).compute_shell_coords(x_max);
                obj.nacelles(i)=obj.nacelles(i).compute_grid();
                obj.nacelles(i)=obj.nacelles(i).compute_grid_flat(obj.nacelles(i).center_coords(3,1),x_max);
                
                grid_len_b4=length(grid_flat);
                grid_flat=[grid_flat obj.nacelles(i).grid_flat];
                panels_flat=[panels_flat,obj.nacelles(i).panels_flat+grid_len_b4];
                te_idx_flat=[te_idx_flat,obj.nacelles(i).te_idx_flat+grid_len_b4];
                obj.nacelles(i).grid_start_idx_flat=grid_len_b4+1;
                obj.nacelles(i).panel_start_idx_flat=panel_len_b4+1; 
                
                obj.nacelles(i)=obj.nacelles(i).update_idx();

                grid_len_b4_vol=length(obj.grid_vol);
                panel_len_b4_vol=length(obj.panels_vol);
                obj.grid_vol=[obj.grid_vol obj.nacelles(i).grid_vol];
                obj.panels_vol=[obj.panels_vol obj.nacelles(i).panels_vol+grid_len_b4_vol];
                obj.te_idx_vol=[obj.te_idx_vol obj.nacelles(i).te_idx_vol+grid_len_b4_vol];
                obj.is_te_vol=[obj.is_te_vol obj.nacelles(i).is_te_vol];
                obj.opposite_te_vol=[obj.opposite_te_vol obj.nacelles(i).opposite_te_vol+panel_len_b4_vol]; 
                obj.nacelles(i).grid_start_idx_vol=grid_len_b4_vol+1;
                obj.nacelles(i).panel_start_idx_vol=panel_len_b4_vol+1;
                
            end 
            obj.grid=grid_flat;
            obj.grid_flat=grid_flat;
            obj.panels=panels_flat;
            obj.te_idx=te_idx_flat;
        end
        
        function obj=write_modes_in_tecplot(obj,aircraft_structure,nmodes,exaturation_factor)
            for i=1:length(obj.control_surfaces)
                obj=obj.f_set_control_surface(obj.control_surfaces{i},0);
            end
            obj=obj.compute_grid();
            for mod=1:nmodes
                for omega=[0 pi/2 3*pi/2]
                    aircraft_structure.nodal_deflections=aircraft_structure.modeshapes(:,mod)*exaturation_factor*sin(omega);
                    aircraft_structure=aircraft_structure.f_postprocess();
                    obj=obj.compute_deflected_grid(aircraft_structure.f_get_deflections);
                    obj.write_grid_deflected(['results/' obj.name '/mode_' num2str(mod) '_' num2str(omega)],[0 0 0]',[0 0 0]',[0 0 0]')
                    %obj.plot_grid_deflected(['results/' obj.name '/mode_' num2str(mod) '_' num2str(omega)]);
                end
            end
        end
        
        function obj=write_modes_in_paraview(obj, aircraft_structure,Modes,exagFactor,folder, varargin)
            if ~isempty(varargin)
                shape=varargin{1};
            else
                shape=[];
            end
                mkdir(folder);
            for i=Modes
                aircraft_structure.nodal_deflections=shape+exagFactor*aircraft_structure.modeshapes(:,i);
                aircraft_structure=aircraft_structure.f_postprocess();
                aircraft=obj.compute_deflected_grid(aircraft_structure.f_get_deflections_c4);
                aircraft.write_grid_deflected([folder, '\modeshape_', sprintf('%.3i',i)],[0 0 0]',[0 0 0]',0);
            end                
        end
        
        function obj=write_mode_animation_in_paraview(obj,aircraft_structure,Modes,exagFactor,folder,frames,varargin)
            if ~isempty(varargin)
                shape=varargin{1};
            else
                shape=aircraft_structure.modeshapes(:,1)*0;
            end
            inputExagFactor=exagFactor;
            exagFactor=ones(1,max(Modes));
            if length(inputExagFactor)==length(Modes)
                exagFactor(Modes)=inputExagFactor;
            end
            if length(inputExagFactor)==1
                exagFactor=ones(1,max(Modes))*inputExagFactor;
            end
            mkdir(folder);
            for i=Modes
                for j=1:frames
                    aircraft_structure.nodal_deflections=shape+cos((j/frames)*2*pi)*exagFactor(i)*aircraft_structure.modeshapes(:,i);
                    aircraft_structure=aircraft_structure.f_postprocess();
                    aircraft=obj.compute_deflected_grid(aircraft_structure.f_get_deflections);
                    aircraft.write_grid_deflected([folder, '\modeshape_', sprintf('%.3i',i),'_frame_', int2str(j)],[0 0 0]',[0 0 0]',0);
                end
            end
        end
        
        function obj=write_modal_time_history_in_paraview(obj,aircraft_structure,history,folder,name)
            mkdir(folder);
            for iFrame=1:size(history,1)
                    aircraft_structure.nodal_deflections=aircraft_structure.modeshapes(:,1:size(history,2))*history(iFrame,:)';
                    aircraft_structure=aircraft_structure.f_postprocess();
                    aircraft=obj.compute_deflected_grid(aircraft_structure.f_get_deflections);
                    aircraft.write_grid_deflected([folder, '\' name '_frame_', int2str(iFrame)],[0 0 0]',[0 0 0]',0);
            end
        end
                
        function obj=plot_grid(obj)
            
            hold on
            
          
            allCSPanels=[];
            
            for wing_count = 1:length(obj.wings)
            
                % loops through all wing-segments in order to highlight panels
                % belonging to control surfaces.
                for seg_count = 1:length(obj.wings(wing_count).wing_segments)

                    % Plots panels belonging to leading edge control surfaces in GREEN
                    if ~isempty(obj.wings(wing_count).wing_segments(seg_count).n_le_panels) && obj.wings(wing_count).wing_segments(seg_count).n_le_panels > 0

                        panels_loc = obj.wings(wing_count).wing_segments(seg_count).le_device.panelIds;
                        panels_locL = obj.wings(wing_count).wing_segments(seg_count).le_device.panelIdsL;

                        for i=1:length(panels_loc)
                            handle= fill3(obj.grid(1,obj.panels(:,panels_loc(i))), obj.grid(2,obj.panels(:,panels_loc(i))),obj.grid(3,obj.panels(:,panels_loc(i))),'g');
                            alpha(handle,1);
                        end
                        
                        if obj.wings(wing_count).wing_segments(seg_count).symmetric
                           for i=1:length(panels_locL)
                               handle= fill3(obj.grid(1,obj.panels(:,panels_locL(i))), obj.grid(2,obj.panels(:,panels_locL(i))),obj.grid(3,obj.panels(:,panels_locL(i))),'g');
                               alpha(handle,1);
                           end
                        end
                        allCSPanels=[allCSPanels panels_loc panels_locL];
                    end

                    % Plots panels belonging to spoiler control surfaces in YELLOW
                    if ~isempty(obj.wings(wing_count).wing_segments(seg_count).n_sp_panels) && obj.wings(wing_count).wing_segments(seg_count).n_sp_panels > 0

                        panels_loc = obj.wings(wing_count).wing_segments(seg_count).sp_device.panelIds_standard;
                        panels_locL = obj.wings(wing_count).wing_segments(seg_count).sp_device.panelIds_standard_L;

                        for i=1:length(panels_loc)
                            handle= fill3(obj.grid(1,obj.panels(:,panels_loc(i))), obj.grid(2,obj.panels(:,panels_loc(i))),obj.grid(3,obj.panels(:,panels_loc(i))),'y');
                            alpha(handle,1);
                        end
                        
                        if obj.wings(wing_count).wing_segments(seg_count).symmetric
                            for i=1:length(panels_locL)
                                handle= fill3(obj.grid(1,obj.panels(:,panels_locL(i))), obj.grid(2,obj.panels(:,panels_locL(i))),obj.grid(3,obj.panels(:,panels_locL(i))),'y');
                                alpha(handle,1);
                            end
                        end
                        allCSPanels=[allCSPanels panels_loc panels_locL];
                    end

                    % Plots panels belonging to trailing edge surfaces in CYAN
                    if ~isempty(obj.wings(wing_count).wing_segments(seg_count).n_te_panels) && obj.wings(wing_count).wing_segments(seg_count).n_te_panels > 0

                        panels_loc = obj.wings(wing_count).wing_segments(seg_count).te_device.panelIds_standard;
                        panels_locL = obj.wings(wing_count).wing_segments(seg_count).te_device.panelIds_standard_L;

                        for i=1:length(panels_loc)
                            handle= fill3(obj.grid(1,obj.panels(:,panels_loc(i))), obj.grid(2,obj.panels(:,panels_loc(i))),obj.grid(3,obj.panels(:,panels_loc(i))),'c');
                            alpha(handle,1);
                        end
                        
                        if obj.wings(wing_count).wing_segments(seg_count).symmetric
                            for i=1:length(panels_locL)
                                handle= fill3(obj.grid(1,obj.panels(:,panels_locL(i))), obj.grid(2,obj.panels(:,panels_locL(i))),obj.grid(3,obj.panels(:,panels_locL(i))),'c');
                                alpha(handle,1);
                            end
                        end
                        allCSPanels=[allCSPanels panels_loc panels_locL];
                        % If there is an overlap between spoilers and flap, the
                        % overlapping region is plotted in RED, and the free
                        % region in CYAN
                        if ~isempty(obj.wings(wing_count).wing_segments(seg_count).te_device.panelIds_special)

                            panels_loc_special = obj.wings(wing_count).wing_segments(seg_count).te_device.panelIds_special;
                            panels_loc_special_L = obj.wings(wing_count).wing_segments(seg_count).te_device.panelIds_special_L;

                            for i=1:length(panels_loc_special)
                                handle= fill3(obj.grid(1,obj.panels(:,panels_loc_special(i))), obj.grid(2,obj.panels(:,panels_loc_special(i))),obj.grid(3,obj.panels(:,panels_loc_special(i))),'r');
                                alpha(handle,1);
                            end
                            
                            if obj.wings(wing_count).wing_segments(seg_count).symmetric
                                for i=1:length(panels_loc_special_L)
                                    handle= fill3(obj.grid(1,obj.panels(:,panels_loc_special_L(i))), obj.grid(2,obj.panels(:,panels_loc_special_L(i))),obj.grid(3,obj.panels(:,panels_loc_special_L(i))),'r');
                                    alpha(handle,1);
                                end
                            end                            
                        end 
                    end
                end
            end
            % plot parametric cs
            for wing_count = 1:length(obj.wings)
                if ~isempty(obj.wings(wing_count).parCS)
                    for iCs=1:length(obj.wings(wing_count).parCS)
                        panels_loc_special = obj.wings(wing_count).parCS{iCs}.panelIds;
                        for i=1:length(panels_loc_special)
                            handle= fill3(obj.grid(1,obj.panels(:,panels_loc_special(i))), obj.grid(2,obj.panels(:,panels_loc_special(i))),obj.grid(3,obj.panels(:,panels_loc_special(i))),'m');
                            alpha(handle,obj.wings(wing_count).parCS{iCs}.fraction(i));
                        end

                        if obj.wings(wing_count).wing_segments(seg_count).symmetric
                            panels_loc_special_L = obj.wings(wing_count).parCS{iCs}.panelIdsL;

                            for i=1:length(panels_loc_special_L)
                                handle= fill3(obj.grid(1,obj.panels(:,panels_loc_special_L(i))), obj.grid(2,obj.panels(:,panels_loc_special_L(i))),obj.grid(3,obj.panels(:,panels_loc_special_L(i))),'m');
                                alpha(handle,obj.wings(wing_count).parCS{iCs}.fraction(i));
                            end
                        end
                        allCSPanels=[allCSPanels panels_loc_special panels_loc_special_L];
                    end
                end
                
            end
             % plots all the grid panels in LIGHT BLUE
             restPanels=setdiff(1:length(obj.panels),allCSPanels);
            for i=1:length(restPanels)
                handle= fill3(obj.grid(1,obj.panels(:,restPanels(i))), obj.grid(2,obj.panels(:,restPanels(i))),obj.grid(3,obj.panels(:,restPanels(i))),'b');
                alpha(handle,0.4);
            end
            axis equal
            axis tight
            grid on
            
        end
        
        function obj=plot_gustZones(obj)
            hold on
            for iZone=1:size(obj.gustZones,1)
                r=rand(1,3);
                for i=obj.gustZones{iZone}
                    handle= fill3(obj.grid(1,obj.panels(:,i)), obj.grid(2,obj.panels(:,i)),obj.grid(3,obj.panels(:,i)),'b');
                    alpha(handle,0.4);
                    handle.FaceColor=r;
                end
            end
            axis equal
            axis tight
            grid on
        end
        
        
        function obj=write_tecplot_volgrid(obj,filename)
            mode='W';
            fileID = fopen([filename '.tp'],mode);
            fprintf(fileID,'TITLE = "Example: Multi-Zone 2D Plot"\n');
            fprintf(fileID,'VARIABLES = "X", "Y", "Z","Cp"\n');
            fprintf(fileID,'zone n=%i, e=%i, f=fepoint, et=quadrilateral \n',length(obj.grid_vol),length(obj.panels_vol));
            
            for i =1:length(obj.grid_vol)
                fprintf(fileID,'%f %f %f %f\n',obj.grid_vol(1,i),obj.grid_vol(2,i),obj.grid_vol(3,i),1);
            end
            
            for i =1:length(obj.panels_vol)
                fprintf(fileID,'%i %i %i %i \n',obj.panels_vol(1,i),obj.panels_vol(2,i),obj.panels_vol(3,i),obj.panels_vol(4,i));
            end

            fclose(fileID);
        end
        
        function obj=write_tecplot_wingbox(obj,filename)

            mode='W';

            fileID = fopen([filename '.tp'],mode);
            fprintf(fileID,'TITLE = "Example: Multi-Zone 2D Plot"\n');
            fprintf(fileID,'VARIABLES = "X", "Y", "Z","Cp"\n');
            for nwing=1:length(obj.wings)
                fprintf(fileID,'zone n=%i, e=%i, f=fepoint, et=lineseg \n',size(obj.wings(nwing).wingbox_coords,2),size(obj.wings(nwing).wingbox_coords,2)-1);
                for i =1:size(obj.wings(nwing).wingbox_coords,2)
                    fprintf(fileID,'%f %f %f %f\n',obj.wings(nwing).wingbox_coords(1,i,1),obj.wings(nwing).wingbox_coords(2,i,1),obj.wings(nwing).wingbox_coords(3,i,1)+obj.wings(nwing).wingbox_height(i,1)/2,1);
                end
                for i=1:(size(obj.wings(nwing).wingbox_coords,2))
                    fprintf(fileID,'%i %i \n',i,i+1);
                end
                fprintf(fileID,'zone n=%i, e=%i, f=fepoint, et=lineseg \n',size(obj.wings(nwing).wingbox_coords,2),size(obj.wings(nwing).wingbox_coords,2)-1);
                for i =1:size(obj.wings(nwing).wingbox_coords,2)
                    fprintf(fileID,'%f %f %f %f\n',obj.wings(nwing).wingbox_coords(1,i,1),obj.wings(nwing).wingbox_coords(2,i,1),obj.wings(nwing).wingbox_coords(3,i,1)-obj.wings(nwing).wingbox_height(i,1)/2,1);
                end
                for i=1:(size(obj.wings(nwing).wingbox_coords,2))
                    fprintf(fileID,'%i %i \n',i,i+1);
                end
                fprintf(fileID,'zone n=%i, e=%i, f=fepoint, et=lineseg \n',size(obj.wings(nwing).wingbox_coords,2),size(obj.wings(nwing).wingbox_coords,2)-1);
                for i =1:size(obj.wings(nwing).wingbox_coords,2)
                    fprintf(fileID,'%f %f %f %f\n',obj.wings(nwing).wingbox_coords(1,i,2),obj.wings(nwing).wingbox_coords(2,i,2),obj.wings(nwing).wingbox_coords(3,i,2)+obj.wings(nwing).wingbox_height(i,2)/2,1);
                end
                for i=1:(size(obj.wings(nwing).wingbox_coords,2))
                    fprintf(fileID,'%i %i \n',i,i+1);
                end
                fprintf(fileID,'zone n=%i, e=%i, f=fepoint, et=lineseg \n',size(obj.wings(nwing).wingbox_coords,2),size(obj.wings(nwing).wingbox_coords,2)-1);
                for i =1:size(obj.wings(nwing).wingbox_coords,2)
                    fprintf(fileID,'%f %f %f %f\n',obj.wings(nwing).wingbox_coords(1,i,2),obj.wings(nwing).wingbox_coords(2,i,2),obj.wings(nwing).wingbox_coords(3,i,2)-obj.wings(nwing).wingbox_height(i,2)/2,1);
                end
                for i=1:(size(obj.wings(nwing).wingbox_coords,2))
                    fprintf(fileID,'%i %i \n',i,i+1);
                end
                
            end
            fclose(fileID);
        end
        
        
        function obj=write_tecplot_grid(obj,filename)
            fileID = fopen(filename,'W');
            fprintf(fileID,'TITLE = "Aircraft"\n');
            fprintf(fileID,'VARIABLES = "X", "Y", "Z"\n');
            
            for i=1:length(obj.wings)
                obj.wings(i).write_tecplot_grid(fileID);
            end
            fclose(fileID);
        end
        
        function obj=plot_grid_vol(obj)
            hold on
            for i=1:length(obj.wings)
                obj.wings(i).plot_grid_vol();
            end
            for i=1:length(obj.fuselages)
                obj.fuselages(i).plot_grid_vol();
            end
            for i=1:length(obj.nacelles)
                obj.nacelles(i).plot_grid_vol();
            end
            axis equal
            axis tight
            grid on
        end
        
        function obj=plot_wingbox_coordinates(obj)
           figure 
           hold on
           for i=1:length(obj.wings)
              plot3(obj.wings(i).wingbox_coords(1,:,1),obj.wings(i).wingbox_coords(2,:,1) ,obj.wings(i).wingbox_coords(3,:,1),'-rx')
              plot3(obj.wings(i).wingbox_coords(1,:,2),obj.wings(i).wingbox_coords(2,:,2) ,obj.wings(i).wingbox_coords(3,:,2),'-bx')
           end
        end
               
        function obj=plot_grid_deflected(obj,filename)
            hold on
            nel_fus=0;
            
%             if ~isempty(obj.fuselages)
%                 for i=1:length(obj.fuselages)
%                     nel_fus=nel_fus+length(obj.fuselages(i).panels);
%                 end
%             end
            for i=1:length(obj.panels)
                handle= fill3(obj.grid_deflected(1,obj.panels(:,i)), obj.grid_deflected(2,obj.panels(:,i)),obj.grid_deflected(3,obj.panels(:,i)),'b');
                alpha(handle,0.4);
            end
            axis equal
            axis tight
            grid on
    
        end
        
        function obj=write_grid_deflected(obj,filename,pos,Euler,ref)
             mode='W';
             pos(1)=0;
                             Lx=[1       0       0
                    0   cos(Euler(1))  -sin(Euler(1))
                    0   sin(Euler(1)) cos(Euler(1))];
                
                Ly=[cos(Euler(2)) 0 -sin(Euler(2))
                    0      1    0
                    sin(Euler(2))  0   cos(Euler(2))];
                
                Lz=[cos(Euler(3)) sin(Euler(3))   0
                    -sin(Euler(3)) cos(Euler(3))  0
                    0           0   1];
                
                
                M_BI=Lz*Ly*Lx;
             
            fileID = fopen([filename '.tp'],mode);
            append=0;
            if append~=1
                fprintf(fileID,'TITLE = "Example: Multi-Zone 2D Plot"\n');
                fprintf(fileID,'VARIABLES = "X", "Y", "Z","Cp"\n');
            end
            
            fprintf(fileID,'zone n=%i, e=%i, f=fepoint, et=quadrilateral \n',length(obj.grid),length(obj.panels));
            
            for i =1:length(obj.grid)
                ri=obj.grid_deflected(:,i)-ref;
                pt=M_BI'*ri+pos(1:3);
                fprintf(fileID,'%f %f %f %f\n',pt(1),pt(2),pt(3),0);
            end
            
            for i =1:length(obj.panels)
                fprintf(fileID,'%i %i %i %i \n',obj.panels(1,i),obj.panels(2,i),obj.panels(3,i),obj.panels(4,i));
            end
            fclose(fileID);
        end

        function obj=f_set_control_surface(obj,name,deflection) 
            
            % cs_idx indicates whether a surface was deflected. If yes, its value
            % will be the index of deflected surface within control_surfaces_parents
            cs_idx = 0;
            
            for i=1:length(obj.control_surfaces_parents)
                if strcmp(obj.control_surfaces_parents(i).name,name)
                    obj.control_surfaces_parents(i).delta = deflection;
                    if obj.control_surfaces_parents(i).is_sym
                        if obj.control_surfaces_parents(i).is_sym_defl
                            obj.control_surfaces_parents(i).delta_l_r=[deflection -deflection];
                        else
                            obj.control_surfaces_parents(i).delta_l_r=[deflection deflection];
                        end
                    else
                        obj.control_surfaces_parents(i).delta_l_r=[deflection];
                    end
                    cs_idx = i;
                end
                
                if strcmp([obj.control_surfaces_parents(i).name '_left'],name)
                    obj.control_surfaces_parents(i).delta_l_r(2)= deflection;
                    cs_idx = i;
                end
                
                if strcmp([ obj.control_surfaces_parents(i).name '_right'],name)
                    obj.control_surfaces_parents(i).delta_l_r(1)= deflection;
                    cs_idx = i;
                end
                
                % _asymm is a new deflection type that can be used for any
                % surface. It deflects surfaces asymmetrically, meaning
                % that left and right displacement are totally independent.
                % As such, the deflection input has to be an array [left deflection, right deflection]
                if strcmp([obj.control_surfaces_parents(i).name '_asymm'],name)
                    obj.control_surfaces_parents(i).delta = deflection;
                    obj.control_surfaces_parents(i).delta_l_r(1) = deflection;
                    obj.control_surfaces_parents(i).delta_l_r(2) = deflection;
                    cs_idx = i;
                end
            end
            
            if cs_idx>0
                wing_index = obj.control_surfaces_parents(cs_idx).wing_index;
                start_seg_idx = obj.control_surfaces_parents(cs_idx).children{1}.seg_idx;
                end_seg_idx = obj.control_surfaces_parents(cs_idx).children{end}.seg_idx;
                
                for seg_count = start_seg_idx:end_seg_idx
                    obj.wings(wing_index).wing_segments(seg_count) = obj.wings(wing_index).wing_segments(seg_count).compute_controlsurface_coordinates();        
                end
            end
            
            for i=1:length(obj.control_surfaces)
                if(strcmp(obj.control_surfaces{i},name))
                   obj.control_deflections{i}=deflection; 
                end
            end
        end
        
        function [] = check_panel_to_beam_element(obj, aircraft_structure,iPanel)
            %% panel to beam element checks
            figure
            hold on
            grid on
            axis equal
            % which wing
            for i=1:size(obj.wings,2) 
                if iPanel<=obj.wings(i).panel_start_idx+size(obj.wings(i).panels,2)
                    wingId=i;
                    break;
                end
            end
            
            %aircraft_ow.plot_grid;
            if iPanel-30<obj.wings(wingId).panel_start_idx;
                startPanel=obj.wings(wingId).panel_start_idx;
            else
                startPanel=iPanel-30;
            end
            if iPanel+30>obj.wings(wingId).panel_start_idx+size(obj.wings(i).panels,2)
                endPanel=obj.wings(wingId).panel_start_idx+size(obj.wings(i).panels,2)-1;
            else
                endPanel=iPanel+30;
            end
            for i=startPanel:endPanel
                h=fill3(obj.grid(1,obj.panels(:,i)),obj.grid(2,obj.panels(:,i)),obj.grid(3,obj.panels(:,i)),'b');
                set(h,'facealpha',.25)
            end
            h=fill3(obj.grid(1,obj.panels(:,iPanel)),obj.grid(2,obj.panels(:,iPanel)),obj.grid(3,obj.panels(:,iPanel)),'r');

            if ~obj.panel_to_beam_element(iPanel,1)==0
                P1=aircraft_structure.beam(wingId).node_coords(obj.panel_to_beam_element(iPanel,1),:);
                P2=aircraft_structure.beam(wingId).node_coords(obj.panel_to_beam_element(iPanel,1)+1,:);
                scatter3(P1(1),P1(2),P1(3),'bs');
                scatter3(P2(1),P2(2),P2(3),'bs');
                Pm=(P1+P2)/2;
                text(Pm(1),Pm(2),Pm(3), num2str(obj.panel_to_beam_element(iPanel,2),3))
            end

            if ~obj.panel_to_beam_element(iPanel,6)==0
                P1=aircraft_structure.beam(wingId).node_coords(obj.panel_to_beam_element(iPanel,6),:);
                P2=aircraft_structure.beam(wingId).node_coords(obj.panel_to_beam_element(iPanel,6)+1,:);
                scatter3(P1(1),P1(2),P1(3),'bs');
                scatter3(P2(1),P2(2),P2(3),'bs');
                Pm=(P1+P2)/2;
                text(Pm(1),Pm(2),Pm(3), num2str(obj.panel_to_beam_element(iPanel,7),3))
            end

            if ~obj.panel_to_beam_element(iPanel,11)==0
                P1=aircraft_structure.beam(wingId).node_coords(obj.panel_to_beam_element(iPanel,11),:);
                P2=aircraft_structure.beam(wingId).node_coords(obj.panel_to_beam_element(iPanel,11)+1,:);
                scatter3(P1(1),P1(2),P1(3),'bs');
                scatter3(P2(1),P2(2),P2(3),'bs');
                Pm=(P1+P2)/2;
                text(Pm(1),Pm(2),Pm(3), num2str(obj.panel_to_beam_element(iPanel,12),3))
            end
        end
        function [] = check_panel_to_beam_element2(obj, aircraft_structure,iBeam,iBe)

            %determine panel range for  wingidx for iBe
            panelRange=[obj.wings(iBeam).panel_start_idx:obj.wings(iBeam).panel_start_idx+size(obj.wings(iBeam).panels,2)-1]';

            %create matrix with panelIds belonging to iBe and Arearatio of Panel
            ptba=[panelRange obj.panel_to_beam_element(panelRange,1:2); panelRange obj.panel_to_beam_element(panelRange,6:7);panelRange obj.panel_to_beam_element(panelRange,11:12)];
            ptba2=ptba(ptba(:,2)==iBe,1:3);

            %plot beamelement
            figure; 
            hold on;
            for i=iBe-5:iBe+5
                 if ~(i>aircraft_structure.beam(iBeam).nel) && i>0
                     P1=aircraft_structure.beam(iBeam).node_coords(i,:);
                     P2=aircraft_structure.beam(iBeam).node_coords(i+1,:);
                     scatter3(P1(1),P1(2),P1(3),'bs');
                     scatter3(P2(1),P2(2),P2(3),'bs');
                     Pm=(P1+P2)/2;
                     text(Pm(1),Pm(2),Pm(3), num2str(i))
                 end
            end
            %check if panels are missing
            missing=setdiff(min(ptba2(:,1)):max(ptba2(:,1)), ptba2(:,1))';
            ptba2=[ptba2; [missing zeros(length(missing),1),  zeros(length(missing),1)]];
            %plotting also couple of panels before minNo and after maxNo
            n=20;
            ptba2=[ptba2; [max(ptba2(:,1))+1:max(ptba2(:,1))+n ; zeros(1,n); zeros(1,n)]'];
            ptba2=[ptba2; [min(ptba2(:,1))-n:min(ptba2(:,1)-1) ; zeros(1,n); zeros(1,n)]'];

            for i=1:size(ptba2,1)
                 if any(panelRange==ptba2(i,1))
                     %plot panel and fill according to area ratio
                     h=fill3(obj.grid(1,obj.panels(:,ptba2(i,1))),obj.grid(2,obj.panels(:,ptba2(i,1))),obj.grid(3,obj.panels(:,ptba2(i,1))),'r');
                     set(h,'facealpha',ptba2(i,3))
                 end
            end
            axis equal;
        end
        function obj = computeControlSurfacePanelIds(obj)
            for iWing=1:length(obj.wings)
                obj.wings(iWing)=obj.wings(iWing).computeControlSurfacePanelIds();
            end
            % correct spoilers piecewise linear
            
%             obj=obj.correctSpoilers('PL');
        end
        function obj = computeHingeMoments(obj, F)
            for iWing=1:length(obj.wings)
                for iSeg=1:length(obj.wings(iWing).wing_segments)
                    if obj.wings(iWing).wing_segments(iSeg).has_te_cs
                        CSHingeLine=[obj.wings(iWing).wing_segments(iSeg).xyz_te_device(:,:,1)' obj.wings(iWing).wing_segments(iSeg).xyz_te_device(:,:,2)'];
                        
                        CSPanels=obj.wings(iWing).wing_segments(iSeg).te_device.panelIds;
                        
                        hingeOrigin=CSHingeLine(:,1);
                        hingeAxis=CSHingeLine(:,2)-CSHingeLine(:,1)./norm(CSHingeLine(:,2)-CSHingeLine(:,1));
                        originMoment=zeros(3,1);
                        for iCSPanel=1:length(CSPanels)
                            fVap=.75*sum(obj.grid(:,obj.panels(1:2,CSPanels(iCSPanel))),2)/2+.25*sum(obj.grid(:,obj.panels(3:4,CSPanels(iCSPanel))),2)/2;
                            originMoment=originMoment+cross(fVap-hingeOrigin,F(:,CSPanels(iCSPanel)));
                        end
                        hingeMoment=dot(originMoment,hingeAxis);
                        if  ~obj.wings(iWing).wing_segments(iSeg).te_device.is_sym
                            [~, iCs] = ismember([obj.wings(iWing).wing_segments(iSeg).te_device.name], obj.control_surfaces);
                            obj.controlSurfaceMoments(iCs)=hingeMoment;
                            %here todo set hingemoment
                        else
                            [~, iCs] = ismember( [obj.wings(iWing).wing_segments(iSeg).te_device.name '_right'], obj.control_surfaces);
                            obj.controlSurfaceMoments(iCs)=hingeMoment;
                            %here todo set hingemoment
                            %mirror hingeline
                            CSHingeLine=[CSHingeLine(1,:); -CSHingeLine(2,:); CSHingeLine(3,:)];
                            %get panels of mirrored CS
                            CSPanels=obj.wings(iWing).wing_segments(iSeg).te_device.panelIdsL;
                            %calculate hingemoment
                            hingeOrigin=CSHingeLine(:,1);
                            hingeAxis=CSHingeLine(:,2)-CSHingeLine(:,1)./norm(CSHingeLine(:,2)-CSHingeLine(:,1));
                            originMoment=zeros(3,1);
                            for iCSPanel=1:length(CSPanels)
                                fVap=.75*sum(obj.grid(:,obj.panels(1:2,CSPanels(iCSPanel))),2)/2+.25*sum(obj.grid(:,obj.panels(3:4,CSPanels(iCSPanel))),2)/2;
                                originMoment=originMoment+cross(fVap-hingeOrigin,F(:,CSPanels(iCSPanel)));
                            end
                            hingeMoment=dot(originMoment,hingeAxis);
                            %set mirrored cs hingemoment
                            [~, iCs] = ismember( [obj.wings(iWing).wing_segments(iSeg).te_device.name '_left'], obj.control_surfaces);
                            obj.controlSurfaceMoments(iCs)=hingeMoment;
                            
                        end
                    end
                end
            end
        end
        function [pass]= runInputChecks(obj)
            % check if segmentboundaries are smooth (twist and profile)
            pass=1;
            for iWing=1:length(obj.wings)
                for iSeg=2:length(obj.wings(iWing).wing_segments)
                    if ~(abs(obj.wings(iWing).wing_segments(iSeg-1).Theta_t-obj.wings(iWing).wing_segments(iSeg).Theta_r) < 1e-8)
                        fprintf('Warning: nonsmooth twist transition between segment %i and segment %i on wing %i \n',iSeg-1,iSeg, iWing)
                        pass=0;
                    end
                    if ~isequal(obj.wings(iWing).wing_segments(iSeg-1).profile_name_t,obj.wings(iWing).wing_segments(iSeg).profile_name_r)
                        if ~isequal(round(obj.wings(iWing).wing_segments(iSeg-1).profile_t,10),round(obj.wings(iWing).wing_segments(iSeg).profile_r,10))
                            fprintf('Warning: nonsmooth profile transition between segment %i and segment %i on wing %i \n',iSeg-1,iSeg, iWing)
                            pass=0;
                        end
                    end
                    if  ~isequal(round(obj.wings(iWing).wing_segments(iSeg-1).c_t,10),round(obj.wings(iWing).wing_segments(iSeg).c_r,10))
                        fprintf('Warning: nonsmooth chord transition between segment %i and segment %i on wing %i \n',iSeg-1,iSeg, iWing)
                        pass=0;
                    end
                end
            end
        end
        function obj = correctSpoilers(obj,type)
            for iCs=1:length(obj.control_surfaces_parents)
                if obj.control_surfaces_parents(iCs).pos==2 %only spoilers
                    for iChild=1:length(obj.control_surfaces_parents(iCs).children)
                        spoilerPanelIds=obj.control_surfaces_parents(iCs).children{iChild}.panelIds;
                        panelsOnSpoiler=find(diff(spoilerPanelIds)~=1,1);
                        panelsPerStrip=diff(spoilerPanelIds(panelsOnSpoiler:panelsOnSpoiler+1))-1+panelsOnSpoiler;
                        panelsBehindSpoiler=find(obj.is_te(spoilerPanelIds(panelsOnSpoiler)+1:spoilerPanelIds(panelsOnSpoiler+1)));
                        panelsBeforeSpoiler=panelsPerStrip-panelsOnSpoiler-panelsBehindSpoiler;

                        newPanelIds=spoilerPanelIds(1)-panelsBeforeSpoiler:spoilerPanelIds(end)+panelsBehindSpoiler;

                        spoilerPanelIdsL=obj.control_surfaces_parents(iCs).children{iChild}.panelIdsL;
                        newPanelIdsL=spoilerPanelIdsL(1)-panelsBeforeSpoiler:spoilerPanelIdsL(end)+panelsBehindSpoiler;
                        nStrips=length(spoilerPanelIds)/panelsOnSpoiler;
                        nPanelsPerStrip=length(newPanelIds)/nStrips;
                        if strcmp(type,'PC') % piecewise constant correction
                            factors=[-0.075,1.1,.5];
                            fractionPerStrip=[factors(1)*ones(panelsBeforeSpoiler,1); factors(2)*ones(panelsOnSpoiler,1); factors(3)*ones(panelsBehindSpoiler,1)];
                            obj.control_surfaces_parents(iCs).children{iChild}.fractions=repmat(fractionPerStrip,nStrips,1);
                        elseif strcmp(type,'PL') % picewise linear correction
                            factors=[-0.075 -0.225 1.1 0.5];
                            fractions=zeros(length(newPanelIds),1);
                            for iStrip=1:nStrips
                                % panel ids for this strip
                                panIdx=newPanelIds((iStrip-1)*nPanelsPerStrip+1:iStrip*nPanelsPerStrip);
                                % get local positions of leading edge, spoiler
                                % leading edge, spoiler trailing edge and
                                % trailing edge
                                xPos=[  (obj.grid(1,obj.panels(1,panIdx(1)))+obj.grid(1,obj.panels(2,panIdx(1))))/2 ,...
                                        (obj.grid(1,obj.panels(1,panIdx(panelsBeforeSpoiler+1)))+obj.grid(1,obj.panels(2,panIdx(panelsBeforeSpoiler+1))))/2 ,...
                                        (obj.grid(1,obj.panels(1,panIdx(panelsBeforeSpoiler+panelsOnSpoiler+1)))+obj.grid(1,obj.panels(2,panIdx(panelsBeforeSpoiler+panelsOnSpoiler+1))))/2,...
                                        (obj.grid(1,obj.panels(3,panIdx(end)))+obj.grid(1,obj.panels(4,panIdx(end))))/2];
                                % x values for picewise linear function
                                xApprox=[xPos(1) xPos(2)-(xPos(3)-xPos(2))/2 (xPos(2)+xPos(3))/2 xPos(4)];
                                xApprox=[xPos(1) xPos(2)-(xPos(3)-xPos(2)) (xPos(2)+xPos(3))/2 xPos(4)];

                                % x coordinates of boundary condition
                                % application points of this strip
                                rdwnX=0.5*(obj.grid(1,obj.panels(2,panIdx))+obj.grid(1,obj.panels(1,panIdx)))*0.25+0.5*(obj.grid(1,obj.panels(3,panIdx))+obj.grid(1,obj.panels(4,panIdx)))*0.75;
                                fractions((iStrip-1)*nPanelsPerStrip+1:iStrip*nPanelsPerStrip)=interp1(xApprox,factors,rdwnX,'linear','extrap');
                            end
                            obj.control_surfaces_parents(iCs).children{iChild}.fractions=fractions;

                        end
                        obj.control_surfaces_parents(iCs).children{iChild}.panelIdsUncorr=obj.control_surfaces_parents(iCs).children{iChild}.panelIds;
                        obj.control_surfaces_parents(iCs).children{iChild}.panelIdsLUncorr= obj.control_surfaces_parents(iCs).children{iChild}.panelIdsL;
                        obj.control_surfaces_parents(iCs).children{iChild}.panelIds=newPanelIds;
                        obj.control_surfaces_parents(iCs).children{iChild}.panelIdsL=newPanelIdsL;
                    end
                end
            end
        end
        function obj = computeChebPoints(obj)
            nPanSeg={};
                for iWing=1:length(obj.wings)
                    temp=[];
                    for iSeg=1:length(obj.wings(iWing).wing_segments)
                        temp=[temp size(obj.wings(iWing).wing_segments(iSeg).panels,2)];
                    end
                    nPanSeg{iWing}=temp;
                end
            for iWing=1:length(obj.wings)
                iSeg=1;
                normPrvMidlines=0;
                spanWisePoints{iWing}=[];
                chordWisePoints{iWing}=[];

                for iSeg=1:length(obj.wings(iWing).wing_segments)
                    nPrvPanels=0;
                    for jWing=1:iWing-1
                        nPrvPanels=nPrvPanels+sum(nPanSeg{jWing});
                        if obj.wings(jWing).symmetric
                            nPrvPanels=nPrvPanels+sum(nPanSeg{jWing});
                        end
                    end
                    for jSeg=1:iSeg-1
                        nPrvPanels=nPrvPanels+(nPanSeg{iWing}(jSeg));
                    end

                    nPanels=nPanSeg{iWing}(iSeg);
                    nSpan=sum(obj.is_te(nPrvPanels+1:nPrvPanels+nPanels));
                    nChord=nPanels/nSpan;
                    innerLEPanel=1+nPrvPanels;
                    outerLEPanel=(nSpan-1)*nChord+1+nPrvPanels;
                    innerTEPanel=nChord+nPrvPanels;
                    outerTEPanel=nPanels+nPrvPanels;

                    midLine=(obj.grid(:,obj.panels(2,outerLEPanel))+obj.grid(:,obj.panels(3,outerTEPanel)))./2-(obj.grid(:,obj.panels(1,innerLEPanel))+obj.grid(:,obj.panels(4,innerTEPanel)))./2;
                    linPoints=normPrvMidlines+linspace(0,1,nSpan*2+1)*norm(midLine);
                    spanWisePoints{iWing}=[spanWisePoints{iWing} reshape(repmat(linPoints(2:2:end-1),nChord,1),1,nPanels)];
                    normPrvMidlines=normPrvMidlines+norm(midLine);
                    for iChord=1:nSpan
                        %midpoints of panels in this strip
                        panelIdx=((iChord-1)*nChord+nPrvPanels+1:(iChord)*nChord+nPrvPanels);
                        midpoints=(obj.grid(:,obj.panels(1,panelIdx))+obj.grid(:,obj.panels(2,panelIdx))+obj.grid(:,obj.panels(3,panelIdx))+obj.grid(:,obj.panels(4,panelIdx)))./4;
                        dist=midpoints(:,2:end)-midpoints(:,1);
                        chordWisePoints{iWing}=[chordWisePoints{iWing} [0 sqrt(sum(dist.^2))./max(sqrt(sum(dist.^2)))]];
                    end

                %     hold on; scatter3(obj.grid(1,obj.panels(:,innerLEPanel)),obj.grid(2,obj.panels(:,innerLEPanel)),obj.grid(3,obj.panels(:,innerLEPanel)),'o')
                %     hold on; scatter3(obj.grid(1,obj.panels(:,outerLEPanel)),obj.grid(2,obj.panels(:,outerLEPanel)),obj.grid(3,obj.panels(:,outerLEPanel)),'+')
                %     hold on; scatter3(obj.grid(1,obj.panels(:,innerTEPanel)),obj.grid(2,obj.panels(:,innerTEPanel)),obj.grid(3,obj.panels(:,innerTEPanel)),'d')
                %     hold on; scatter3(obj.grid(1,obj.panels(:,outerTEPanel)),obj.grid(2,obj.panels(:,outerTEPanel)),obj.grid(3,obj.panels(:,outerTEPanel)),'s')

                end
                if obj.wings(iWing).symmetric
                    spanWisePoints{iWing}=([normPrvMidlines+spanWisePoints{iWing} normPrvMidlines-spanWisePoints{iWing}])/(2*normPrvMidlines)*2-1;
                    spanWiseLength{iWing}=2*normPrvMidlines;
                    chordWisePoints{iWing}=[chordWisePoints{iWing} chordWisePoints{iWing}]*2-1;
                else
                    spanWisePoints{iWing}=spanWisePoints{iWing}/normPrvMidlines*2-1;
                    spanWiseLength{iWing}=normPrvMidlines;
                    chordWisePoints{iWing}=chordWisePoints{iWing}*2-1;
                    
                end
            end

            obj.dataForCheb.spanWisePoints=spanWisePoints;
            obj.dataForCheb.spanWiseLength=spanWiseLength;
            obj.dataForCheb.chordWisePoints=chordWisePoints;
            
        end
    end
end

