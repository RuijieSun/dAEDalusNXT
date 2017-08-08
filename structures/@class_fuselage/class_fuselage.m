%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%
classdef class_fuselage<class_beam
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        %% Element positioning and type
        %> twist at each node position
        epsilon;   
        %> dihedral at each node position
        nu;
        %> sweep at each node position
        phi;                             
        %> flag if it is a rigid beam
        is_rigid=0;
        %% general values
        %> Fuselage volume
        Vfuselage=0;
        %> Stiffened shell total mass
        m_shell_total=0.0;
        %> fuel volume
        fuel_volume_total=0.0;

        %% Attached parts
        %> flag for wingmounted engines
        fuselagemountedengines=0;
        %> flag for wingmounted gears
        fuselagemountedgears=0;
        %> flag for adding gear forces
        landingimpact=0;
        %> number of engines
        n_engines=0;
        %> number of wing mounted gears
        n_gears=0; 
        %> array of engine classes    
        engine;
        %> array of engine classes
        gear;  
        
        n_PAX;
        
        
    end
    
    methods (Access=public)
        
         function obj = class_fuselage(nel, crosssections, varargin)
            %call class_beam constructor 
            obj=obj@class_beam(nel, crosssections, varargin);
         end
         
         function struct_deflections=f_get_deflections(obj)
            struct_deflections.def=obj.nodal_deflections'; 
         end
         
         function obj=f_add_engine(obj,engine)
             obj.n_engines=obj.n_engines+1;
             if(obj.n_engines==1)
                obj.engine=engine;
             else
                obj.engine(obj.n_engines)=engine;
             end
         end
         
         function obj=f_add_gear(obj,gear)
             obj.n_gears=obj.n_gears+1;
             if(obj.n_gears==1)
                obj.gear=gear;
             else
                obj.gear(obj.n_gears)=gear;
             end
         end
         
         function obj=f_set_state(obj,state)
             obj.load_factor=state.load_factor;
             obj.loadcase_index=state.loadcase_index;
             obj.landingimpact=0;
         end
         
         % plot_critical_case_idx(wing,varargin);
            
         % initialize beam with estimated distributed mass
         % obj = f_init_wingmass(obj,weights);
         
         % calculate estimated wing mass
         obj = f_calc_mass(obj,weights);
         
         % initialize beam manually 
         obj = f_init_stdBeam(obj,E,G,m);
         
%          obj = f_init_coords(le,phi,nu,epsilon);
%          
%          obj = f_init_geomProp(obj,r,t);
         
         obj = f_init_std_structure(obj,nu,phi,Ix,Iy,Iz,J,A,le);
         
         % assemble system matrices for solving
         obj = f_assemble(obj,add_eigenmass,add_fuelmass,add_engineforces,add_gearforces);
         
         % initialize elements with material properties
         obj=f_init_material_properties(obj,structure);

         
%          function structmesh = get_struct_mesh(obj)
%             structmesh.y=0:obj.wing_frontview_length/obj.nel:obj.wing_frontview_length;
%          end     
         
%          function plot_externalforces(wing,add_eigenmass,add_fuelmass,engine,gear,varargin)
%             hold on
% 
%            
%             midp=zeros(3,1);
%             
%             norm_def=wing.nodal_deflections;%/max(wing.nodal_deflections);
%             
%             x_dist=zeros(length(wing.beamelement)+1,1);
%             y_dist=zeros(length(wing.beamelement)+1,1);
%             z_dist=zeros(length(wing.beamelement)+1,1);
%             
%                 mid_def=zeros(3,1);
%             mid_rot=zeros(3,1);
%             
%             for i=1:1:length(wing.beamelement)
%                 for j=1:3
%                     midp(j)=(wing.node_coords(i+1,j)+wing.node_coords(i,j))/2;
%                     mid_def(j)=(norm_def(j+6*(i-1))+norm_def(6+j+6*(i-1)))/2;
%                     mid_rot(j)=(wing.nodal_deflections(3+j+6*(i-1))+wing.nodal_deflections(9+j+6*(i-1)))/2;
%                 end
%                 
%                 midp_defl=[midp(1)+mid_def(1),midp(2)+mid_def(2),midp(3)+mid_def(3)]';
%               
%                 midp=midp_defl;
%                 h=wing.beamelement(i).crosssection.h;
%                 w=wing.beamelement(i).crosssection.w;
%                 
%                 le=wing.beamelement(i).le;
%                 nu=wing.beamelement(i).nu;
%                 twist=wing.beamelement(i).epsilon;
%                 z_dist(i)=h/2;
%                 x_dist(i)=w/2;
%                 z_dist(i+1)=h/2;
%                 x_dist(i+1)=w/2;
%                 plotcube([midp(1),midp(2),midp(3)],[w,le,h],[nu+mid_rot(1),twist+mid_rot(2),mid_rot(3)],[1 1 1 1 1 1 1 1],0.3,1);
%             end
%             
%             axis equal
%             grid on
%             
%             qx=cell2mat({wing.beamelement(:).qx});
%             qy=cell2mat({wing.beamelement(:).qy});
%             qz=cell2mat({wing.beamelement(:).qz});
%             halfspan=sum(cell2mat({wing.beamelement(:).le}));
%             
%             max_q=max([abs(qx),abs(qy),abs(qz)]);
%             
%             scale=2*max_q/halfspan;
%             
%             optargin = size(varargin,2);
%             if optargin==1
%                 scale=varargin{1};
%             end
%             
%             qmrot=zeros(3,length(wing.beamelement));
%             qfrot=zeros(3,length(wing.beamelement));
%             
%             for i=1:length(wing.beamelement)    
%                 if add_eigenmass
%                     T=wing.beamelement(i).f_rotVec(0,wing.gamma,0);  % transformations matrix from NED to aerodynamic system
%                     qmrot(1:3,i)=T*[0,0,wing.beamelement(i).qm]';%wing.beamelement(i).qm];           % distributed loading due to eigenmass in x,y,z coordinates
%                 end
%                 
%                 if add_fuelmass
%                     T=wing.beamelement(i).f_rotVec(0,wing.gamma,0);  % transformations matrix from NED to aerodynamic system
%                     qfrot(:,i)=T*[0;0;wing.beamelement(i).qf];           % distributed loading due to eigenmass in x,y,z coordinates
%                 end
%             end
% 
%             xxx=subplot(1,1,1);
%             %plot3(wing.node_coords(:,1),wing.node_coords(:,2),wing.node_co
%             %ords(:,3),'-k','LineWidth',2)
%             wing.node_coords(:,1)=wing.node_coords(:,1)+wing.nodal_deflections(1:6:end);
%             wing.node_coords(:,2)=wing.node_coords(:,2)+wing.nodal_deflections(2:6:end);
%             wing.node_coords(:,3)=wing.node_coords(:,3)+wing.nodal_deflections(3:6:end);
%             
%             
%             wing.node_coords(:,3)=wing.node_coords(:,3)+z_dist;
%             plotglobalforce(xxx,qz,scale,wing.node_coords(:,:),'z',0,'b',0); 
%             wing.node_coords(:,3)=wing.node_coords(:,3)-2.*z_dist;
%             plotglobalforce(xxx,qmrot(3,:),scale,wing.node_coords(:,:),'z',0.45,'g',0); 
%             %wing.node_coords(:,3)=wing.node_coords(:,3)+2*z_dist;
%             wing.node_coords(1:end-1,3)=wing.node_coords(1:end-1,3)+qmrot(3,:)'./scale;
%             plotglobalforce(xxx,qfrot(3,:),scale,wing.node_coords(:,:),'z',0.9,'m',0); 
%             wing.node_coords(1:end-1,3)=wing.node_coords(1:end-1,3)-qmrot(3,:)'./scale;
%             wing.node_coords(:,3)=wing.node_coords(:,3)+z_dist;
%             
%             wing.node_coords(:,1)=wing.node_coords(:,1)+x_dist;
%             plotglobalforce(xxx,qx,scale,wing.node_coords(:,:),'x',0.3,'c',0); 
%             axis equal
%             grid on   
%             wing.node_coords(:,1)=wing.node_coords(:,1)-x_dist;
%             
%             wing.node_coords(:,1)=wing.node_coords(:,1)-wing.nodal_deflections(1:6:end);
%             wing.node_coords(:,2)=wing.node_coords(:,2)-wing.nodal_deflections(2:6:end);
%             wing.node_coords(:,3)=wing.node_coords(:,3)-wing.nodal_deflections(3:6:end);
%             %plot3(wing.node_coords(:,1),wing.node_coords(:,2),wing.node_coords(:,3),'-x','LineWidth',2)
%             %plot3(engine.cgpos(1),engine.cgpos(2),engine.cgpos(3),'c+');
%             %plot3(gear.pos(1),gear.pos(2),gear.pos(3),'go');
%           end
%          
         % ================================================================
         %> @brief plot deformations of the wing
         %>
         %> @param opt1 optional: color for plot 
         %> @param opt2 optional: transparency for plot 
         % ================================================================
         function plot_deformations(fuselage,varargin)
%             
%             col=0.1;
%             tr=0.5;
%             % Color and transparency
%             optargin = size(varargin,2);
%             if optargin==1
%                xcol=varargin{1};
%                xcol=cell2mat(xcol);
%                col=xcol(1);
%                tr=xcol(2);
%             end
            hold on
            
            caxis([0 1]);
            aux_points=linspace(0,2*pi,18)';
            for i=1:length(fuselage.beamelement)
                
                [r]=fuselage.beamelement(i).crosssection.get_dimensions();
              
                cy=r*cos(aux_points);
                cz=r*sin(aux_points);
                
                x_front=ones(18,1)*fuselage.node_coords(i,1);
                y_front=ones(18,1)*fuselage.node_coords(i,2)+cy;
                z_front=ones(18,1)*fuselage.node_coords(i,3)+cz;                
                
                x_rear=ones(18,1)*fuselage.node_coords(i+1,1);
                y_rear=ones(18,1)*fuselage.node_coords(i+1,2)+cy;
                z_rear=ones(18,1)*fuselage.node_coords(i+1,3)+cz;
                
                a_front=-fuselage.nodal_deflections(5+6*(i-1),1)*0;
                b_front=-fuselage.nodal_deflections(6+6*(i-1),1)*0;
                c_front=-fuselage.nodal_deflections(4+6*(i-1),1)*0;
                
                a_rear=-fuselage.nodal_deflections(11+6*(i-1),1)*0;
                b_rear=-fuselage.nodal_deflections(12+6*(i-1),1)*0;
                c_rear=-fuselage.nodal_deflections(10+6*(i-1),1)*0;
                
                T_front=T_plot(a_front,b_front,c_front);
                T_rear=T_plot(a_rear,b_rear,c_rear);
                
                x_plot=[];
                y_plot=[];
                z_plot=[];
                for j=1:18
                    coord_front=T_front*[x_front(j,1);y_front(j,1);z_front(j,1)];
                    coord_rear=T_rear*[x_rear(j,1);y_rear(j,1);z_rear(j,1)];
                    
                    x_plot=[x_plot; coord_front(1) coord_rear(1)];
                    y_plot=[y_plot; coord_front(2) coord_rear(2)];
                    z_plot=[z_plot; coord_front(3) coord_rear(3)];
                end
                
                x_defl=[ones(18,1)*fuselage.nodal_deflections(1+6*(i-1),1) ones(18,1)*fuselage.nodal_deflections(1+6*(i-1),1)];
                y_defl=[ones(18,1)*fuselage.nodal_deflections(2+6*(i-1),1) ones(18,1)*fuselage.nodal_deflections(2+6*(i-1),1)];
                z_defl=[ones(18,1)*fuselage.nodal_deflections(3+6*(i-1),1) ones(18,1)*fuselage.nodal_deflections(3+6*(i-1),1)];
                                
                x_plot=x_plot+x_defl;
                y_plot=y_plot+y_defl;
                z_plot=z_plot+z_defl;
                
                surface(x_plot,y_plot,z_plot,0.1*ones(18,2),'FaceAlpha',0.5);
            end
            
            axis equal
            grid on
         end         

          
         
        % ================================================================
         %> @brief plot geometry of the wing
         %>
         %> @param opt1 optional: color for plot 
         %> @param opt2 optional: transparency for plot 
         % ================================================================
          function plot_geometry(fuselage,varargin)
            hold on
            
            caxis([0 1]) 
            
%             col=1;
%             tr=0.5;
%             
%             optargin = size(varargin,2);
%             if optargin==1
%                xcol=varargin{1};
%                xcol=cell2mat(xcol);
%                col=xcol(1);
%                tr=xcol(2);
%             end

            aux_points=linspace(0,2*pi,18)';

            for i=1:length(fuselage.beamelement)
                
                [r]=fuselage.beamelement(i).crosssection.get_dimensions();
              
                cy=r*cos(aux_points);
                cz=r*sin(aux_points);
                
                x_front=ones(18,1)*fuselage.node_coords(i,1);
                y_front=ones(18,1)*fuselage.node_coords(i,2)+cy;
                z_front=ones(18,1)*fuselage.node_coords(i,3)+cz;                
                
                x_rear=ones(18,1)*fuselage.node_coords(i+1,1);
                y_rear=ones(18,1)*fuselage.node_coords(i+1,2)+cy;
                z_rear=ones(18,1)*fuselage.node_coords(i+1,3)+cz;
                
                surface([x_front x_rear],[y_front y_rear],[z_front z_rear],ones(18,2),'FaceAlpha',0.5);
                
            end
            axis equal
            grid on
          end  
         
          function plot_t(fuselage,varargin)
              
              if nargin==2
                  color=varargin{1};
              else
                  color='blue';
              end
              
              for i=1:length(fuselage.beamelement)
                  t=fuselage.beamelement(i).crosssection.t_sk_eq;
                  x=fuselage.node_coords(i,1);
                  plot(x,t,'*','Color',color)
                  hold on
              end
             grid on 
          end
          
          
          function write_structure_tecplot(fuselage,fileID,beam_nr)

            for i=1:1:length(fuselage.beamelement) 
              t_sk_eq(i)=cell2mat({fuselage.beamelement(i).crosssection.t_sk_eq});
            end
            
            aux_points=linspace(0,2*pi,18)';
            nodal_def=[fuselage.nodal_deflections(1:6:end) fuselage.nodal_deflections(2:6:end) fuselage.nodal_deflections(3:6:end)];
            for i=1:length(fuselage.beamelement)
                [r]=fuselage.beamelement(i).crosssection.get_dimensions();
              
                cy=r*cos(aux_points);
                cz=r*sin(aux_points);
                
                cy_inner=(r-t_sk_eq(i))*cos(aux_points);
                cz_inner=(r-t_sk_eq(i))*sin(aux_points);
                
                x_front=ones(18,1)*fuselage.node_coords(i,1)+nodal_def(i,1);
                y_front=ones(18,1)*fuselage.node_coords(i,2)+cy+nodal_def(i,2);
                z_front=ones(18,1)*fuselage.node_coords(i,3)+cz+nodal_def(i,3);     
                y_front_inner=ones(18,1)*fuselage.node_coords(i,2)+cy_inner+nodal_def(i,2);
                z_front_inner=ones(18,1)*fuselage.node_coords(i,3)+cz_inner+nodal_def(i,3);

                x_rear=ones(18,1)*fuselage.node_coords(i+1,1)+nodal_def(i+1,1);
                y_rear=ones(18,1)*fuselage.node_coords(i+1,2)+cy+nodal_def(i+1,2);
                z_rear=ones(18,1)*fuselage.node_coords(i+1,3)+cz+nodal_def(i+1,3);
                y_rear_inner=ones(18,1)*fuselage.node_coords(i+1,2)+cy_inner+nodal_def(i+1,2);
                z_rear_inner=ones(18,1)*fuselage.node_coords(i+1,3)+cz_inner+nodal_def(i+1,3);
                
                for kk=1:17
                    fprintf(fileID,'ZONE T="%s_%s_1" I=2, J=2, K=2 F=BLOCK \n',num2str(beam_nr),num2str(i));
                    fprintf(fileID,'%f %f %f %f %f %f %f %f \n\n',x_front(kk), x_front(kk+1), x_front(kk), x_front(kk+1), x_rear(kk),x_rear(kk+1), x_rear(kk), x_rear(kk+1));
                    fprintf(fileID,'%f %f %f %f %f %f %f %f \n\n',y_front_inner(kk), y_front_inner(kk+1), y_front(kk), y_front(kk+1), y_rear_inner(kk),y_rear_inner(kk+1), y_rear(kk), y_rear(kk+1));
                    fprintf(fileID,'%f %f %f %f %f %f %f %f \n\n',z_front_inner(kk), z_front_inner(kk+1), z_front(kk), z_front(kk+1), z_rear_inner(kk),z_rear_inner(kk+1), z_rear(kk), z_rear(kk+1));
                    fprintf(fileID,'%f %f %f %f %f %f %f %f \n\n',t_sk_eq(i),t_sk_eq(i),t_sk_eq(i),t_sk_eq(i),t_sk_eq(i),t_sk_eq(i),t_sk_eq(i),t_sk_eq(i));
                end
            end
            
        
%             if ~isempty(fuselage.engine)
%                 for ne=1:length(fuselage.engine)
%                     dist=zeros(fuselage.nel,1);
%                     for i=1:fuselage.nel
%                         dist(i)=sqrt(sum((fuselage.node_coords(i,1)-fuselage.engine(ne).cg_pos(1)').^2));
%                     end
%                     [Y,I] = min(dist);
%                     hold on
%                      
%                     plot3([fuselage.engine(ne).cg_pos(1) fuselage.node_coords(I,1)],[fuselage.engine(ne).cg_pos(2) fuselage.node_coords(I,2)],[fuselage.engine(ne).cg_pos(3) fuselage.node_coords(I,3)],'-k','LineWidth',2);
%                     plot3(fuselage.engine(ne).cg_pos(1),fuselage.engine(ne).cg_pos(2),fuselage.engine(ne).cg_pos(3),'-kx','MarkerSize',15,'LineWidth',2);
%                 end
%             end
%             
%            if ~isempty(fuselage.gear)
%                 for ne=1:length(fuselage.gear)
%                     dist=zeros(fuselage.nel,1);
%                     for i=1:fuselage.nel
%                         dist(i)=sqrt(sum((fuselage.node_coords(i,1)-fuselage.gear(ne).pos(1)').^2));
%                     end
%                     [Y,I] = min(dist);
%                     hold on
%                      
%                     plot3([fuselage.gear(ne).pos(1) fuselage.node_coords(I,1)],[fuselage.gear(ne).pos(2) fuselage.node_coords(I,2)],[fuselage.gear(ne).pos(3) fuselage.node_coords(I,3)],'-k','LineWidth',2);
%                     plot3(fuselage.gear(ne).pos(1),fuselage.gear(ne).pos(2),fuselage.gear(ne).pos(3),'-ko','MarkerSize',15,'LineWidth',2);
%                     
% 
%                 end
%            end
          end
         
         function plot_structure(fuselage,varargin)

            hold on
            
            for i=1:1:length(fuselage.beamelement) 
              t_sk_eq(i)=cell2mat({fuselage.beamelement(i).crosssection.t_sk_eq});
            end
            
            optargin = size(varargin,2);
            if optargin==2
                t_min=varargin{2};
                t_max=varargin{1};
            else
                t_max=0;
                t_min=100;
                t_max=max([t_sk_eq t_max])
                t_min=min([t_sk_eq t_min])
            end
            
            aux_points=linspace(0,2*pi,18)';

            for i=1:length(fuselage.beamelement)
                
                col_i=(t_sk_eq(i)-t_min)/(t_max-t_min);
                
                [r]=fuselage.beamelement(i).crosssection.get_dimensions();
              
                cy=r*cos(aux_points);
                cz=r*sin(aux_points);
                
                x_front=ones(18,1)*fuselage.node_coords(i,1);
                y_front=ones(18,1)*fuselage.node_coords(i,2)+cy;
                z_front=ones(18,1)*fuselage.node_coords(i,3)+cz;                
                
                x_rear=ones(18,1)*fuselage.node_coords(i+1,1);
                y_rear=ones(18,1)*fuselage.node_coords(i+1,2)+cy;
                z_rear=ones(18,1)*fuselage.node_coords(i+1,3)+cz;
                
                surface([x_front x_rear],[y_front y_rear],[z_front z_rear],ones(18,2)*col_i);
            end
            
            axis equal
            grid on 
            if ~isempty(fuselage.engine)
                for ne=1:length(fuselage.engine)
                    dist=zeros(fuselage.nel,1);
                    for i=1:fuselage.nel
                        dist(i)=sqrt(sum((fuselage.node_coords(i,1)-fuselage.engine(ne).cg_pos(1)').^2));
                    end
                    [Y,I] = min(dist);
                    hold on
                    plot3([fuselage.engine(ne).cg_pos(1) fuselage.node_coords(I,1)],[fuselage.engine(ne).cg_pos(2) fuselage.node_coords(I,2)],[fuselage.engine(ne).cg_pos(3) fuselage.node_coords(I,3)],'-k','LineWidth',2);
                    plot3(fuselage.engine(ne).cg_pos(1),fuselage.engine(ne).cg_pos(2),fuselage.engine(ne).cg_pos(3),'-kx','MarkerSize',15,'LineWidth',2);
                end
            end
            
            if ~isempty(fuselage.gear)
                for ne=1:length(fuselage.gear)
                    dist=zeros(fuselage.nel,1);
                    for i=1:fuselage.nel
                        dist(i)=sqrt(sum((fuselage.node_coords(i,1)-fuselage.gear(ne).pos(1)').^2));
                    end
                    [Y,I] = min(dist);
                    hold on
                    plot3([fuselage.gear(ne).pos(1) fuselage.node_coords(I,1)],[fuselage.gear(ne).pos(2) fuselage.node_coords(I,2)],[fuselage.gear(ne).pos(3) fuselage.node_coords(I,3)],'-k','LineWidth',2);
                    plot3(fuselage.gear(ne).pos(1),fuselage.gear(ne).pos(2),fuselage.gear(ne).pos(3),'-ko','MarkerSize',15,'LineWidth',2);
                end
            end
         end

% 
%          
%          function plot_internalforces(wing,varargin)
%             delta_x=0;
%                 delta_y=0;
%                 delta_z=0;
%                 med_x=0;
%                 med_y=0;
%                 med_z=0;
%                 scale=0;
%                 delta=0;
%                 deltaxy=0;
%                 
%             aQx=wing.node_loadings_loc(1:6:end);
%             aQy=wing.node_loadings_loc(2:6:end);
%             aQz=wing.node_loadings_loc(3:6:end);
%             aMx=wing.node_loadings_loc(4:6:end);
%             aMy=wing.node_loadings_loc(5:6:end);
%             aMz=wing.node_loadings_loc(6:6:end);
% 
%             optargin = size(varargin,2);
%             if optargin==15
%                 delta_x=varargin{1};
%                 delta_y=varargin{2};
%                 delta_z=varargin{3};
%                 med_x=varargin{4};
%                 med_y=varargin{5};
%                 med_z=varargin{6};
%                 scale=0.6*sqrt((delta_x)^2+(delta_z)^2+(delta_y)^2);
%                 delta=max([delta_y,delta_z]);
%                 deltaxy=max([delta_y,delta_x]);
%                 scale_bendingmoment=varargin{7};
%                 scale_torsionmoment=varargin{8};
%                 scale_force=varargin{9};
%                 maxMx=varargin{10};
%                 maxMt=varargin{11};
%                 maxMz=varargin{12};
%                 maxQx=varargin{13};
%                 maxQy=varargin{14};
%                 maxQz=varargin{15};
%             else 
%                 max_x=max(wing.node_coords(:,1));
%                 min_x=min(wing.node_coords(:,1));
%             
%                 max_y=max(wing.node_coords(:,2));
%                 min_y=min(wing.node_coords(:,2));
%             
%                 max_z=max(wing.node_coords(:,3));
%                 min_z=min(wing.node_coords(:,3)); 
%              
%                 delta_x=max_x-min_x;
%                 delta_z=max_z-min_z;
%                 delta_y=max_y-min_y;
% 
%                 med_y=(max_y+min_y)/2;
%                 med_z=(max_z+min_z)/2;
%                 med_x=(max_x+min_x)/2;
%                 
%                  scale=0.6*sqrt((delta_x)^2+(delta_z)^2+(delta_y)^2);
%                 delta=max([delta_y,delta_z]);
%                 deltaxy=max([delta_y,delta_x]);
%                 
%                      max_Mb=max([abs(aMx)' abs(aMz)' abs(aMy)']);
%                     %scale torsional moments
%                  max_Mt=max_Mb;
%                  max_Q =max([abs(aQx)' abs(aQy)' abs(aQz)']);
%             
%                     scale_bendingmoment=(2*max_Mb)/scale;
%                     scale_torsionmoment=(2*max_Mt)/scale;
%                                 scale_force=(2*max_Q)/scale;
%                                 
%                 maxMx=max(abs(aMx));
%                 maxMt=max(abs(aMy));
%                 maxMz=max(abs(aMz));
%                 maxQx=max(abs(aQx));
%                 maxQy=max(abs(aQy));
%                 maxQz=max(abs(aQz));
%             end
%              
%             % plot Mx 
%             handle1=subplot(3,2,2);
%             plotforce(handle1,aMx,scale_bendingmoment,wing.node_coords,'z',1,'r',1);
%            
%             title('Mx')
%             view(-90,0)
%             ylim([med_y-0.6*delta,med_y+0.6*delta]);
%             zlim([med_z-0.6*delta,med_z+0.6*delta]);
%             ylabel('Y')
%             zlabel('Z')
%             text(0,delta,delta/1.5,['Mx_{max}=' sprintf('%4.3gMNm',maxMx/10^6)],'BackgroundColor','white');
%             axis equal
%             
%             % plot Qz 
%             handle2=subplot(3,2,1);
%             plotforce(handle2,aQz,scale_force,wing.node_coords,'z',1,'r',1);
%             title('Qz')
%             view(-90,0)
%             ylim([med_y-0.6*delta,med_y+0.6*delta]);
%             zlim([med_z-0.6*delta,med_z+0.6*delta]);
%             ylabel('Y')
%             zlabel('Z')
%             text(0,delta,delta/1.5,['Qz_{max}=' sprintf('%4.3gMN',maxQz/10^6)],'BackgroundColor','white');
%             axis equal
%             
%             handle4=subplot(3,2,3);
%             plotforce(handle4,aQy,scale_force,wing.node_coords,'z',1,'r',1);
%             title('Qy')
%             view(-90,0)
%             ylim([med_y-0.6*delta,med_y+0.6*delta]);
%             zlim([med_z-0.6*delta,med_z+0.6*delta]);
%             ylabel('Y')
%             zlabel('Z')
%             text(0,delta,delta/1.5,['Qy_{max}=' sprintf('%4.3gMN',maxQy/10^6)],'BackgroundColor','white');
%             axis equal
%                 
%             handle5=subplot(3,2,4);
%             plotforce(handle5,aMy,scale_torsionmoment,wing.node_coords,'z',1,'r',1);
%             title('Mt')
%             view(-90,0)
%             ylabel('Y')
%             zlabel('Z')
%             ylim([med_y-0.6*delta,med_y+0.6*delta]);
%             zlim([med_z-0.6*delta,med_z+0.6*delta]);
%             text(0,delta,delta/1.5,['Mt_{max}=' sprintf('%4.3gMNm',maxMt/10^6)],'BackgroundColor','white');
%             axis equal
%             
%             handle3=subplot(3,2,6);
%             plotforce(handle3,aMz,scale_bendingmoment,wing.node_coords,'x',1,'r',1);
%             title('Mz')
%             view(-90,90)
%             ylim([med_y-0.6*deltaxy,med_y+0.6*deltaxy]);
%             xlim([med_x-0.3*deltaxy,med_x+0.3*deltaxy]);
%             ylabel('Y')
%             xlabel('X')
%               text(0,delta,-delta/2,['Mz_{max}=' sprintf('%4.3gMNm',maxMz/10^6)],'BackgroundColor','white');
%             axis equal
%             
%             handle6=subplot(3,2,5);
%             plotforce(handle6,aQx,scale_force,wing.node_coords,'x',1,'r',1);
%             title('Qx')
%             view(-90,90)
%             ylim([med_y-0.6*deltaxy,med_y+0.6*deltaxy]);
%             xlim([med_x-0.3*deltaxy,med_x+0.3*deltaxy]);
%             ylabel('Y')
%             xlabel('X')
%             text(0,delta,-delta/2,['Qx_{max}=' sprintf('%4.3gMN',maxQx/10^6)],'BackgroundColor','white');
%             axis equal
%          end   
    end
    
end

