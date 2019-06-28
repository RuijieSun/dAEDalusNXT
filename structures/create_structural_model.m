%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%

% Due to implementation of anisotropic materials, structure_solver_settings
% is now part of the input. The object contains information on the
% discretization of the beam-element cross sections through the cross
% sectional modeller
function [aircraft,aircraft_structure]=create_structural_model(aircraft,structure_solver_settings)

for nwings=1:length(aircraft.wings_structural_properties)
    wing_settings=class_structural_settings_wing;
    
    % The nr of nodes and maximum shell width for the cross sectional
    % discretization are saved within wing_settings
    wing_settings.nr_nodes_crosssection_element=structure_solver_settings.nr_nodes_crosssection_element;
    wing_settings.ds_nodes_crosssection_element=structure_solver_settings.ds_nodes_crosssection_element;
    
    % reads the name of the material saved within the XML file. If the
    % material type is listed as "anisotropic" it runs
    % class_material_anisotropic and sets  anisotropic_flag = 1.
    % Otherwise class_material is run and the anisotropic_flag = 0.
    if strcmp(aircraft.wings_structural_properties(nwings).material(1),'anisotropic')
        anisotropic_flag = 1;
        if isfield(aircraft.addSettings,'kdf')
            material=class_material_anisotropic(aircraft.addSettings.skins_ply_material,aircraft.addSettings.kdf);
        else
            material=class_material_anisotropic(aircraft.addSettings.skins_ply_material); %assumes that spars use same ply material as skin
        end
    else
        material=class_material(aircraft.wings_structural_properties(nwings).material(1));
        anisotropic_flag = 0;
    end
    wing_settings=wing_settings.f_set_material(material);
    %modify minThicknesses if specified in XML
    if isfield(aircraft.addSettings, 'wingTminSK')
        wing_settings.t_min_sk=aircraft.addSettings.wingTminSK;
    end
    if isfield(aircraft.addSettings, 'wingTminSP')
        wing_settings.t_min_sp=aircraft.addSettings.wingTminSP;
    end
    aircraft=aircraft.compute_wingbox_coords();
    
    n_elems = (size(aircraft.wings(nwings).wingbox_coords,2)-1)*(aircraft.wings(nwings).symmetric+1);
    
    % Checks if the wing is made of anisotropic material. If so,
    % class_crosssection_wingbox_anisotropic is run. Otherwise
    % class_crosssection_wingbox is run.
    if anisotropic_flag == 1
        wingbox_crosssections = repmat([class_crosssection_wingbox_anisotropic()], 1, n_elems);
    else
        wingbox_crosssections = repmat([class_crosssection_wingbox()], 1, n_elems);
    end
    
    if isfield(aircraft.wings(nwings).wing_segments(1).structural_properties, 'fs_segments')
        c_elems = 1;
        for k = 1:(1 + aircraft.wings(nwings).symmetric)
            for i = 1:length(aircraft.wings(nwings).wing_segments)
                wing_segment = aircraft.wings(nwings).wing_segments(i);
                for j = 1:length(wing_segment.structural_properties.t_fs)-1
                    wingbox_crosssections(c_elems).t_sp_fr = wing_segment.structural_properties.t_fs(j);
                    wingbox_crosssections(c_elems).t_sp_re = wing_segment.structural_properties.t_rs(j);
                    wingbox_crosssections(c_elems).t_sk_up = wing_segment.structural_properties.t_ts(j);
                    wingbox_crosssections(c_elems).t_sk_lo = wing_segment.structural_properties.t_bs(j);
                    wingbox_crosssections(c_elems).segment_index = i;
                    c_elems = c_elems + 1;
                end
            end
        end
    end
    
    % Implementation of anisotropic materials requires to include
    % anisotropic_flag within the inputs to class_wing
    wingstructure=class_wing(n_elems, wingbox_crosssections, anisotropic_flag, aircraft.wings(nwings).name);
    %modify fuelingfactor if specified in XML
    if isfield(aircraft.addSettings, 'wingboxFuelingFactor')
        wingstructure=wingstructure.f_set_fuelingFactor(aircraft.addSettings.wingboxFuelingFactor);
    end
    wingstructure=wingstructure.f_init_structure(aircraft.wings(nwings),aircraft.wings_structural_properties(nwings).is_fueled);
    
    % set rib bays
    if isfield(aircraft.wings_structural_properties(nwings),'ribLocations')
        if ~isempty(aircraft.wings_structural_properties(nwings).ribLocations)
            wingstructure=wingstructure.f_init_ribLocations(aircraft.wings_structural_properties(nwings).ribLocations);
        end
    end  
    if wingstructure.isExternalFEM==0
        wingstructure=wingstructure.f_init_material_properties(wing_settings, anisotropic_flag, aircraft.addSettings);
        
        % if the wing is anisotropic, the cross sectional modeller is run
        % automatically
        if anisotropic_flag == 1
            wingstructure = wingstructure.cross_sectional_modeler(wing_settings);
        end
        
        wingstructure.update_M=1;
        wingstructure.update_K=1;
        wingstructure.update_Q=1;
    elseif wingstructure.isExternalFEM==1
        % TODO: do material properties and allowed stresses need to be given?
        
        % Load external FEM mass and stiffness matrices
        iWing = [];
        nameTemp = wingstructure.identifier;
        for j=1:length(aircraft.wings)
            if strcmp(nameTemp, aircraft.wings(j).name)
                iWing = j;
            end
        end
        
        wingstructure.externalFEM.Mext = importdata(aircraft.wings(iWing).pathMassMatrix);
        wingstructure.externalFEM.Kext = importdata(aircraft.wings(iWing).pathStiffnessMatrix);
        
        wingstructure.update_M=1;
        wingstructure.update_K=1;
        wingstructure.update_Q=1;
        
        for j=1:length(wingstructure.beamelement)
            wingstructure.beamelement(j).m = 0;
        end
    end
    for nengines=1:length(aircraft.engines)
        if strcmp(aircraft.engines(nengines).mounting,aircraft.wings(nwings).name)
            wingstructure=wingstructure.f_add_engine(aircraft.engines(nengines));
            wingstructure.wingmountedengines=1;
        end
    end 
    if nwings==1
        aircraft_structure=class_beam_collection(wingstructure);
    else
        aircraft_structure=aircraft_structure.f_add_beam(wingstructure);
    end
end

for nfuse=1:length(aircraft.fuselages)
    fuselage_settings=class_structural_settings_fuselage;  
    fuselage_settings.delta_pressure=8.9632e+04;
    fuselage_settings=fuselage_settings.f_set_material(class_material('aluminum'));
    aircraft=aircraft.compute_shell_coords();
    aircraft.fuselages(nfuse)=aircraft.fuselages(nfuse).compute_shell_coords(aircraft.grid_settings.dy_max_struct_grid);
    
    n_elems = size(aircraft.fuselages(nfuse).center_coords,2)-1;
    fuselage_crosssections = repmat([class_crosssection_fuselage()], 1, n_elems);
    
    fuselage_structure=class_fuselage(n_elems, fuselage_crosssections, aircraft.fuselages(nfuse).name);
    
    fuselage_structure=fuselage_structure.f_init_structure(aircraft.fuselages(nfuse));
    for i=1:size(fuselage_structure.beamelement,2)
        fuselage_structure.beamelement(i).crosssection.delta_pressure=fuselage_settings.delta_pressure;
    end
    
    if fuselage_structure.isExternalFEM==0
        fuselage_structure=fuselage_structure.f_init_material_properties(fuselage_settings);
    elseif fuselage_structure.isExternalFEM==1
        % TODO: do material properties and allowed stresses need to be given?
        
        % Load external FEM mass and stiffness matrices
        iFuselage = [];
        nameTemp = fuselage_structure.identifier;
        for j=1:length(aircraft.fuselages)
            if strcmp(nameTemp, aircraft.fuselages(j).name)
                iFuselage = j;
            end
        end
        
        fuselage_structure.externalFEM.Mext = importdata(aircraft.fuselages(iFuselage).pathMassMatrix);
        fuselage_structure.externalFEM.Kext = importdata(aircraft.fuselages(iFuselage).pathStiffnessMatrix);
        
        fuselage_structure.update_M=1;
        fuselage_structure.update_K=1;
        fuselage_structure.update_Q=1;
        
        for j=1:length(fuselage_structure.beamelement)
            fuselage_structure.beamelement(j).m = 0;
        end
    end
    
   % fuselage_structure=fuselage_structure.f_add_boundary_condition(class_boundary_condition(1,[1 1 1 1 1 1],[0 0 0 0 0 0]));
   % fuselage_structure=fuselage_structure.f_add_boundary_condition(class_boundary_condition(ceil((abs(aircraft.fuselages(1).fuselage_segments(1,1).pos(1,1))+aircraft.wings(1,1).wing_segments(1,1).pos(1,1))/aircraft.grid_settings.dy_max_struct_grid),[1 1 1 1 1 1],[0 0 0 0 0 0]));
   desiredClampingPoint = [aircraft.clamping_conditions{2},aircraft.clamping_conditions{3},aircraft.clamping_conditions{4}]';
   distanceToWing = fuselage_structure.node_coords(:,1)-desiredClampingPoint(1);
   [~,boundaryIdx] = min(abs(distanceToWing));
   fuselage_structure = fuselage_structure.f_add_boundary_condition(class_boundary_condition(boundaryIdx,[1 1 1 1 1 1],[0 0 0 0 0 0]));

   aircraft_structure=aircraft_structure.f_add_beam(fuselage_structure);
end

n_cc=length(aircraft.boundary_conditions);
for ne=1:length(aircraft.engines)
    %find closest node:
    for nw=1:length(aircraft.wings)
        if(strcmp(aircraft.engines(ne).mounting,aircraft.wings(nw).name))
            dist=zeros(length(aircraft_structure.beam(nw).node_coords(:,2)),1);
            for i=1:length(aircraft_structure.beam(nw).node_coords(:,2))
                dist(i)=sqrt(sum((aircraft_structure.beam(nw).node_coords(i,2)-aircraft.engines(ne).cg_pos(2)').^2));
            end
            [Y,node_idx_beam] = min(dist);
            %add forces and moments
            pyl_coords=[aircraft_structure.beam(nw).node_coords(node_idx_beam,:)' aircraft.engines(ne).cg_pos'];
            pylon=class_pylon(2,'none',['Engine' num2str(ne) 'Pylon' ]);
            geof.center_coords=pyl_coords;
            pylon=pylon.f_init_structure(geof);
            pylon.beamelement(1).el_m_sys=pylon.beamelement(1).m*pylon.beamelement(1).le;
            pylon.beamelement(2).el_m_sys=aircraft.engines(ne).m-pylon.beamelement(1).m;
            connection{1}=aircraft.engines(ne).mounting;
            connection{2}=['Engine' num2str(ne) 'Pylon' ];
            aircraft_structure=aircraft_structure.f_add_beam(pylon);
            aircraft.boundary_conditions{n_cc+ne}=connection;
        end
    end
    for nf=1:length(aircraft.fuselages)
        if(strcmp(aircraft.engines(ne).mounting,aircraft.fuselages(nf).name))
            dist=zeros(length(aircraft_structure.beam(nf+nw).node_coords(:,1)),1);
            for i=1:length(aircraft_structure.beam(nf+nw).node_coords(:,1))
                dist(i)=sqrt(sum((aircraft_structure.beam(nf+nw).node_coords(i,1)-aircraft.engines(ne).cg_pos(1)').^2));
            end
            [Y,node_idx_beam] = min(dist);
            %add forces and moments
            pyl_coords=[aircraft_structure.beam(nf+nw).node_coords(node_idx_beam,:)' aircraft.engines(ne).cg_pos'];
            pylon=class_pylon(2,'none',['Engine' num2str(ne) 'Pylon' ]);
            geof.center_coords=pyl_coords;
            pylon=pylon.f_init_structure(geof);
            pylon.beamelement(1).el_m_sys=pylon.beamelement(1).m*pylon.beamelement(1).le;
            pylon.beamelement(2).el_m_sys=aircraft.engines(ne).m-pylon.beamelement(1).m;
            connection{1}=aircraft.engines(ne).mounting;
            connection{2}=['Engine' num2str(ne) 'Pylon' ];
            aircraft_structure=aircraft_structure.f_add_beam(pylon);
            aircraft.boundary_conditions{n_cc+ne}=connection;
        end
    end
end

for nbc=1:length(aircraft.boundary_conditions)
    balljoint=0;
    hingejoint=0;
    zjoint=0;
    for nbeam=1:length(aircraft_structure.beam)
        if strcmp(aircraft.boundary_conditions{nbc}{1},aircraft_structure.beam(nbeam).identifier)
           idx_1=nbeam; 
        end
        if strcmp(aircraft.boundary_conditions{nbc}{2},aircraft_structure.beam(nbeam).identifier)
           idx_2=nbeam; 
        end
        if length(aircraft.boundary_conditions{nbc})==3
            if strcmp(aircraft.boundary_conditions{nbc}{3},'BallJoint')
                 balljoint=1;
            elseif strcmp(aircraft.boundary_conditions{nbc}{3},'ZJoint')
                 zjoint=1;
            elseif strcmp(aircraft.boundary_conditions{nbc}{3},'KardanJoint')
                 hingejoint=1;
            end
        end
    end
    dis=1E6;
    dis_min=dis;
    node_idx_1=0;
    node_idx_2=0;
    for nc_ctr1=1:size(aircraft_structure.beam(idx_1).node_coords,1)
         for nc_ctr2=1:size(aircraft_structure.beam(idx_2).node_coords,1)
            dis=norm(aircraft_structure.beam(idx_1).node_coords(nc_ctr1,:)'-aircraft_structure.beam(idx_2).node_coords(nc_ctr2,:)');
            if dis<dis_min
               dis_min=dis;
               node_idx_1=nc_ctr1;
               node_idx_2=nc_ctr2;
            elseif abs(dis-dis_min)<1E-5
               node_idx_1=[node_idx_1 nc_ctr1];
               node_idx_2=[node_idx_2 nc_ctr2];  
            end
         end
    end
    for nnodes=1:length(node_idx_1)
        if balljoint==0 && zjoint==0
            aircraft_structure=aircraft_structure.f_add_coupling_condition([idx_1 idx_2],[node_idx_1(nnodes) node_idx_2(nnodes)],[1 1 1 1 1 1]);
        elseif balljoint==1
            aircraft_structure=aircraft_structure.f_add_coupling_condition([idx_1 idx_2],[node_idx_1(nnodes) node_idx_2(nnodes)],[1 1 1 0 1 1]); 
        elseif zjoint==1
            aircraft_structure=aircraft_structure.f_add_coupling_condition([idx_1 idx_2],[node_idx_1(nnodes) node_idx_2(nnodes)],[0 0 1 0 0 0]); 
        elseif hingejoint==1
            aircraft_structure=aircraft_structure.f_add_coupling_condition([idx_1 idx_2],[node_idx_1(nnodes) node_idx_2(nnodes)],[1 1 1 0 1 1]); 
        end
    end
end

% pre-merge coupling conditions
del_idx=[];
for i=1:length(aircraft_structure.coupling_condition)-1
   for j=i+1:length(aircraft_structure.coupling_condition)
       for kk=1:min([length(aircraft_structure.coupling_condition(i).beam_idx) length(aircraft_structure.coupling_condition(j).beam_idx)])
            if (aircraft_structure.coupling_condition(i).beam_idx(kk)==aircraft_structure.coupling_condition(j).beam_idx(kk))&&(aircraft_structure.coupling_condition(i).node_idx(kk)==aircraft_structure.coupling_condition(j).node_idx(kk))
                if kk==1
                    aircraft_structure.coupling_condition(i).beam_idx=[aircraft_structure.coupling_condition(i).beam_idx aircraft_structure.coupling_condition(j).beam_idx(2)];
                    aircraft_structure.coupling_condition(i).node_idx=[aircraft_structure.coupling_condition(i).node_idx aircraft_structure.coupling_condition(j).node_idx(2)];
                 %   aircraft_structure.coupling_condition(i).dof=[aircraft_structure.coupling_condition(i).dof; aircraft_structure.coupling_condition(j).dof];
                elseif kk==2
                    aircraft_structure.coupling_condition(i).beam_idx=[aircraft_structure.coupling_condition(i).beam_idx aircraft_structure.coupling_condition(j).beam_idx(1)];
                    aircraft_structure.coupling_condition(i).node_idx=[aircraft_structure.coupling_condition(i).node_idx aircraft_structure.coupling_condition(j).node_idx(1)]; 
                   % aircraft_structure.coupling_condition(i).dof=[aircraft_structure.coupling_condition(i).dof; aircraft_structure.coupling_condition(j).dof];
                end
                del_idx=[del_idx j];
            end
       end
   end
end

del_idx=unique(del_idx);
if ~isempty(del_idx)
    temp_cc=[];
    jj=1;
    for i=1:length(aircraft_structure.coupling_condition)
        if i==del_idx(jj)
            if jj<length(del_idx)
                jj=jj+1;
            end
        else
            temp_cc=[temp_cc aircraft_structure.coupling_condition(i)];
        end
    end

    aircraft_structure.coupling_condition=temp_cc;
end
aircraft_structure=aircraft_structure.f_calc_mass(aircraft.weights);
aircraft_structure.identifier=aircraft.name;

% pre-merge coupling conditions
del_idx=[];
for i=1:length(aircraft_structure.coupling_condition)-1
   for j=i+1:length(aircraft_structure.coupling_condition)
       for kk=1:min([length(aircraft_structure.coupling_condition(i).beam_idx) length(aircraft_structure.coupling_condition(j).beam_idx)])
            if (aircraft_structure.coupling_condition(i).beam_idx(kk)==aircraft_structure.coupling_condition(j).beam_idx(kk))&&(aircraft_structure.coupling_condition(i).node_idx(kk)==aircraft_structure.coupling_condition(j).node_idx(kk))
                if kk==1
                    aircraft_structure.coupling_condition(i).beam_idx=[aircraft_structure.coupling_condition(i).beam_idx aircraft_structure.coupling_condition(j).beam_idx(2)];
                    aircraft_structure.coupling_condition(i).node_idx=[aircraft_structure.coupling_condition(i).node_idx aircraft_structure.coupling_condition(j).node_idx(2)];
                 %   aircraft_structure.coupling_condition(i).dof=[aircraft_structure.coupling_condition(i).dof; aircraft_structure.coupling_condition(j).dof];
                elseif kk==2
                    aircraft_structure.coupling_condition(i).beam_idx=[aircraft_structure.coupling_condition(i).beam_idx aircraft_structure.coupling_condition(j).beam_idx(1)];
                    aircraft_structure.coupling_condition(i).node_idx=[aircraft_structure.coupling_condition(i).node_idx aircraft_structure.coupling_condition(j).node_idx(1)]; 
                   % aircraft_structure.coupling_condition(i).dof=[aircraft_structure.coupling_condition(i).dof; aircraft_structure.coupling_condition(j).dof];
                end
                del_idx=[del_idx j];
            end
       end
   end
end

del_idx=unique(del_idx);
if ~isempty(del_idx)
    temp_cc=[];
    jj=1;
    for i=1:length(aircraft_structure.coupling_condition)
        if i==del_idx(jj)
            if jj<length(del_idx)
                jj=jj+1;
            end
        else
            temp_cc=[temp_cc aircraft_structure.coupling_condition(i)];
        end
    end

    aircraft_structure.coupling_condition=temp_cc;
end
aircraft_structure=aircraft_structure.f_calc_mass(aircraft.weights);
aircraft_structure=aircraft_structure.f_set_solver_settings(structure_solver_settings);
