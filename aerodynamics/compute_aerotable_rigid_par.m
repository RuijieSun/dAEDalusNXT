%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
% compute_aerodata_rigid file modified to use parallel processing
function [Rigid_Aerodata] = compute_aerotable_rigid_par(aircraft,ac_state,varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
state=ac_state.aerodynamic_state;
deflected=0;
control_surfs=aircraft.control_surfaces;
drag=1;
if nargin==2
    aircraft=aircraft.compute_grid();
    config=0;
else
    deflected=varargin{3};%1
    def=varargin{1};
    aircraft=aircraft.compute_deflected_grid(def);
    aircraft=aircraft.compute_grid();

    config=1;
    %standardGridding
    vecConf=length(ac_state);
    if nargin>=4
        vecDeltaR=varargin{2}.vecDeltaR;
        vecDeltaP=varargin{2}.vecDeltaP;
        vecDeltaQ=varargin{2}.vecDeltaQ;
        vecCsDefl=varargin{2}.vecCsDefl;
        vecAlpha=varargin{2}.vecAlpha;
        vecMa=varargin{2}.vecMa;
        vecBeta=varargin{2}.vecBeta;
        
    else        
        vecDeltaR=-30:10:30;
        vecDeltaP=-30:10:30;
        vecDeltaQ=-30:10:30;
        vecCsDefl=-30:5:30;
        vecAlpha=-10:2:15;
        vecMa=0.0:0.05:0.2;
        vecBeta=-10:5:10;
    end
end


conf=0;
Rigid_AerodataStatic_CoeffCX=zeros(length(vecConf),length(vecMa),length(vecAlpha),length(vecBeta));
Rigid_AerodataStatic_CoeffCY=zeros(length(vecConf),length(vecMa),length(vecAlpha),length(vecBeta));
Rigid_AerodataStatic_CoeffCZ=zeros(length(vecConf),length(vecMa),length(vecAlpha),length(vecBeta));
Rigid_AerodataStatic_CoeffCl=zeros(length(vecConf),length(vecMa),length(vecAlpha),length(vecBeta));
Rigid_AerodataStatic_CoeffCm=zeros(length(vecConf),length(vecMa),length(vecAlpha),length(vecBeta));
Rigid_AerodataStatic_CoeffCn=zeros(length(vecConf),length(vecMa),length(vecAlpha),length(vecBeta));

Rigid_AerodataDamping_CoeffCXp=zeros(length(vecConf),length(vecMa),length(vecAlpha),length(vecBeta));
Rigid_AerodataDamping_CoeffCYp=zeros(length(vecConf),length(vecMa),length(vecAlpha),length(vecBeta));
Rigid_AerodataDamping_CoeffCZp=zeros(length(vecConf),length(vecMa),length(vecAlpha),length(vecBeta));
Rigid_AerodataDamping_CoeffClp=zeros(length(vecConf),length(vecMa),length(vecAlpha),length(vecBeta));
Rigid_AerodataDamping_CoeffCmp=zeros(length(vecConf),length(vecMa),length(vecAlpha),length(vecBeta));
Rigid_AerodataDamping_CoeffCnp=zeros(length(vecConf),length(vecMa),length(vecAlpha),length(vecBeta));

Rigid_AerodataDamping_CoeffCXq=zeros(length(vecConf),length(vecMa),length(vecAlpha),length(vecBeta));
Rigid_AerodataDamping_CoeffCYq=zeros(length(vecConf),length(vecMa),length(vecAlpha),length(vecBeta));
Rigid_AerodataDamping_CoeffCZq=zeros(length(vecConf),length(vecMa),length(vecAlpha),length(vecBeta));
Rigid_AerodataDamping_CoeffClq=zeros(length(vecConf),length(vecMa),length(vecAlpha),length(vecBeta));
Rigid_AerodataDamping_CoeffCmq=zeros(length(vecConf),length(vecMa),length(vecAlpha),length(vecBeta));
Rigid_AerodataDamping_CoeffCnq=zeros(length(vecConf),length(vecMa),length(vecAlpha),length(vecBeta));

Rigid_AerodataDamping_CoeffCXr=zeros(length(vecConf),length(vecMa),length(vecAlpha),length(vecBeta));
Rigid_AerodataDamping_CoeffCYr=zeros(length(vecConf),length(vecMa),length(vecAlpha),length(vecBeta));
Rigid_AerodataDamping_CoeffCZr=zeros(length(vecConf),length(vecMa),length(vecAlpha),length(vecBeta));
Rigid_AerodataDamping_CoeffClr=zeros(length(vecConf),length(vecMa),length(vecAlpha),length(vecBeta));
Rigid_AerodataDamping_CoeffCmr=zeros(length(vecConf),length(vecMa),length(vecAlpha),length(vecBeta));
Rigid_AerodataDamping_CoeffCnr=zeros(length(vecConf),length(vecMa),length(vecAlpha),length(vecBeta));

%counterst
% i=1;
% j=1;
% k=1;
% l=1;
tic

if deflected
    wingaero=class_VLM_solver(aircraft.grid_deflected,aircraft.te_idx,aircraft.panels,state,aircraft.reference); 
else
    wingaero=class_VLM_solver(aircraft.grid,aircraft.te_idx,aircraft.panels,state,aircraft.reference); 
end
wingaero=wingaero.f_solve_std();
wingaero=wingaero.f_solve_full();
t_stop=toc

tic;
wingaero=wingaero.f_set_state(state);
wingaero=wingaero.set_grid(aircraft.grid,aircraft.panels);
%wingaero=class_VLM_solver(aircraft.grid,aircraft.te_idx,aircraft.panels,state,aircraft.reference);
wingaero=wingaero.f_solve_std();
wingaero=wingaero.f_solve_full();
t_stop2=toc;

%time estimation
t_loop=t_stop;
t_loop2=t_stop2;
t_total=0;
for conf=vecConf
    for Ma=vecMa
        for alpha=vecAlpha
            for beta=vecBeta
                t_total=t_total+t_loop;
            end
        end
    end
end

for conf=vecConf
    for im=1:length(vecMa)
        Ma=vecMa(im);
        for alpha=vecAlpha
            for beta=vecBeta
                  for cs=1:length(control_surfs)
                    if aircraft.control_surfaces{cs}(1:6) == 'ailero'
                        cs_defl_grid = vecDeltaP;
                    elseif aircraft.control_surfaces{cs}(1:6) == 'elevat'
                        cs_defl_grid = vecDeltaQ;
                    elseif aircraft.control_surfaces{cs}(1:6) == 'rudder'
                        cs_defl_grid = vecDeltaR;
                    else
                        cs_defl_grid =vecCsDefl;
                    end
                    for cs_defl=cs_defl_grid
                        t_total=t_total+t_loop;
                    end  
                  end
            end
        end
    end
end

fprintf('Estimated computation time: %f hours \n',t_total/3600);
%vecConf=1;
% static derivatives

% Uinf=200;                Correction
a=state.V_A/state.Ma;%324.5786;          %     correction
rho_air=state.rho_air;%0.7;
mu=state.mu;
for conf_cs=1:length(aircraft.control_surfaces)
    aircraft=aircraft.f_set_control_surface(aircraft.control_surfaces{conf_cs},0);
end
if deflected    
    aircraft=aircraft.compute_grid();
    aircraft=aircraft.compute_deflected_grid(def);
    wingaero=class_VLM_solver(aircraft.grid_deflected,aircraft.te_idx,aircraft.panels,state,aircraft.reference);
else
    aircraft=aircraft.compute_grid();
    wingaero=class_VLM_solver(aircraft.grid,aircraft.te_idx,aircraft.panels,state,aircraft.reference);
end
lengthAlpha=length(vecAlpha);
lengthBeta=length(vecBeta);
lengthConf=length(vecConf);
lengthMa=length(vecMa);

coeff_lim=1E-10;

for l=1:lengthConf
    conf=vecConf(l);
%     k=1;
    if length(ac_state)>1
        for conf_cs=1:length(ac_state(conf).aircraft_state.control_deflections)
            aircraft=aircraft.f_set_control_surface(ac_state(conf).aircraft_state.control_surfaces{conf_cs},ac_state(conf).aircraft_state.control_deflections{conf_cs}); %set all control surfaces to zero!!
        end
    end
    if deflected
        aircraft=aircraft.compute_grid();
        aircraft=aircraft.compute_deflected_grid(def);
    else
        aircraft=aircraft.compute_grid();
    end

    parfor k=1:lengthMa
        Ma=vecMa(k);
        alpha=0;
        beta=0;
        if Ma
        Uinf=Ma*a;
        else
            Uinf=200;
        end%correction
%         j=1;
        state=class_aero_state(Uinf,alpha,beta,Ma,rho_air,mu);
        aircraftma=aircraft;
        aircraftma=aircraftma.compute_CD_f(state,aircraftma.reference.S_ref);
        
        for j=1:lengthAlpha
            alpha=vecAlpha(j);
%             i=1;
            for i=1:lengthBeta
                beta=vecBeta(i);
                disp(['Ma=' num2str(Ma) ' aoa=' num2str(alpha) ' beta=' num2str(beta)])
                state=class_aero_state(Uinf,alpha,beta,Ma,rho_air,mu);
                if deflected
                    wingaero=class_VLM_solver(aircraftma.grid_deflected,aircraftma.te_idx,aircraftma.panels,state,aircraftma.reference);
                else
                    
                    wingaero=class_VLM_solver(aircraftma.grid,aircraftma.te_idx,aircraftma.panels,state,aircraftma.reference);
                end
                %  wingaero=wingaero.f_set_state(state);

                wingaero=wingaero.f_solve_std();
                wingaero=wingaero.f_solve_full(aircraftma.CD_f);
                if drag==0
                    wingaero.Cdi=0;
                    aircraftma.CD_f=0;
                else
                    %wingaero.Cdi=wingaero.compute_trefftz_drag_mex();
                end
%                 wingaerotest(i)=wingaero; %debug
%                 cdaic(i,:)=wingaero.cdaic;
                C= [wingaero.CX, wingaero.CY, wingaero.CZ, wingaero.CL, wingaero.CM, wingaero.CN];%correction
                Cp=[wingaero.CXp,wingaero.CYp,wingaero.CZp,wingaero.CLp,wingaero.CMp,wingaero.CNp];
                Cq=[wingaero.CXq,wingaero.CYq,wingaero.CZq,wingaero.CLq,wingaero.CMq,wingaero.CNq];
                Cr=[wingaero.CXr,wingaero.CYr,wingaero.CZr,wingaero.CLr,wingaero.CMr,wingaero.CNr];
                % static derivatives
                % static force derivatives
                Rigid_AerodataStatic_CoeffCX(l,k,j,i)=-C(1);
                Rigid_AerodataStatic_CoeffCY(l,k,j,i)=C(2);
                Rigid_AerodataStatic_CoeffCZ(l,k,j,i)=-C(3);
                % static moment derivatives
                Rigid_AerodataStatic_CoeffCl(l,k,j,i)=-C(4);
                Rigid_AerodataStatic_CoeffCm(l,k,j,i)=C(5);
                Rigid_AerodataStatic_CoeffCn(l,k,j,i)=-C(6);
                % damping derivatives
                % p derivatives
                Rigid_AerodataDamping_CoeffCXp(l,k,j,i)=-Cp(1);
                Rigid_AerodataDamping_CoeffCYp(l,k,j,i)=Cp(2);
                Rigid_AerodataDamping_CoeffCZp(l,k,j,i)=-Cp(3);
                Rigid_AerodataDamping_CoeffClp(l,k,j,i)=Cp(4);
                Rigid_AerodataDamping_CoeffCmp(l,k,j,i)=Cp(5);
                Rigid_AerodataDamping_CoeffCnp(l,k,j,i)=Cp(6);
                
                % q derivatives
                Rigid_AerodataDamping_CoeffCXq(l,k,j,i)=-Cq(1);%*cosd(alpha)+Cq(3)*sind(alpha);
                Rigid_AerodataDamping_CoeffCYq(l,k,j,i)=Cq(2);
                Rigid_AerodataDamping_CoeffCZq(l,k,j,i)=-Cq(3);%*cosd(alpha)-Cq(1)*sind(alpha);
                Rigid_AerodataDamping_CoeffClq(l,k,j,i)=-Cq(4);
                Rigid_AerodataDamping_CoeffCmq(l,k,j,i)=Cq(5);
                Rigid_AerodataDamping_CoeffCnq(l,k,j,i)=-Cq(6);
                
                % p derivatives
                Rigid_AerodataDamping_CoeffCXr(l,k,j,i)=-Cr(1);
                Rigid_AerodataDamping_CoeffCYr(l,k,j,i)=Cr(2);
                Rigid_AerodataDamping_CoeffCZr(l,k,j,i)=-Cr(3);
                Rigid_AerodataDamping_CoeffClr(l,k,j,i)=Cr(4);
                Rigid_AerodataDamping_CoeffCmr(l,k,j,i)=Cr(5);
                Rigid_AerodataDamping_CoeffCnr(l,k,j,i)=Cr(6);
             
                %                 C;
                %                 Ma;
                %                 alpha;
            end
%                j=j+1;
        end
%           k=k+1;
    end
%         l=l+1;
end
   Rigid_Aerodata.Static_Coeff.CX(1,:,:,:)=Rigid_AerodataStatic_CoeffCX(1,:,:,:);
    Rigid_Aerodata.Static_Coeff.CY(1,:,:,:)=Rigid_AerodataStatic_CoeffCY(1,:,:,:);
    Rigid_Aerodata.Static_Coeff.CZ(1,:,:,:)=Rigid_AerodataStatic_CoeffCZ(1,:,:,:);
    % static moment derivatives
    Rigid_Aerodata.Static_Coeff.Cl(1,:,:,:)=Rigid_AerodataStatic_CoeffCl(1,:,:,:);
    Rigid_Aerodata.Static_Coeff.Cm(1,:,:,:)=Rigid_AerodataStatic_CoeffCm(1,:,:,:);
    Rigid_Aerodata.Static_Coeff.Cn(1,:,:,:)=Rigid_AerodataStatic_CoeffCn(1,:,:,:);
    % damping derivatives
    % p derivatives
    Rigid_Aerodata.Damping_Coeff.CXp(1,:,:,:)=Rigid_AerodataDamping_CoeffCXp(1,:,:,:);
    Rigid_Aerodata.Damping_Coeff.CYp(1,:,:,:)=Rigid_AerodataDamping_CoeffCYp(1,:,:,:);
    Rigid_Aerodata.Damping_Coeff.CZp(1,:,:,:)=Rigid_AerodataDamping_CoeffCZp(1,:,:,:);
    Rigid_Aerodata.Damping_Coeff.Clp(1,:,:,:)=Rigid_AerodataDamping_CoeffClp(1,:,:,:);
    Rigid_Aerodata.Damping_Coeff.Cmp(1,:,:,:)=Rigid_AerodataDamping_CoeffCmp(1,:,:,:);
    Rigid_Aerodata.Damping_Coeff.Cnp(1,:,:,:)=Rigid_AerodataDamping_CoeffCnp(1,:,:,:);
    
    % q derivatives
    Rigid_Aerodata.Damping_Coeff.CXq(1,:,:,:)=Rigid_AerodataDamping_CoeffCXq(1,:,:,:);%*cosd(alpha)+Cq(3)*sind(alpha);
    Rigid_Aerodata.Damping_Coeff.CYq(1,:,:,:)=Rigid_AerodataDamping_CoeffCYq(1,:,:,:);
    Rigid_Aerodata.Damping_Coeff.CZq(1,:,:,:)=Rigid_AerodataDamping_CoeffCZq(1,:,:,:);%*cosd(alpha)-Cq(1)*sind(alpha);
    Rigid_Aerodata.Damping_Coeff.Clq(1,:,:,:)=Rigid_AerodataDamping_CoeffClq(1,:,:,:);
    Rigid_Aerodata.Damping_Coeff.Cmq(1,:,:,:)=Rigid_AerodataDamping_CoeffCmq(1,:,:,:);
    Rigid_Aerodata.Damping_Coeff.Cnq(1,:,:,:)=Rigid_AerodataDamping_CoeffCnq(1,:,:,:);
    
    % p derivatives
    Rigid_Aerodata.Damping_Coeff.CXr(1,:,:,:)=Rigid_AerodataDamping_CoeffCXr(1,:,:,:);
    Rigid_Aerodata.Damping_Coeff.CYr(1,:,:,:)=Rigid_AerodataDamping_CoeffCYr(1,:,:,:);
    Rigid_Aerodata.Damping_Coeff.CZr(1,:,:,:)=Rigid_AerodataDamping_CoeffCZr(1,:,:,:);
    Rigid_Aerodata.Damping_Coeff.Clr(1,:,:,:)=Rigid_AerodataDamping_CoeffClr(1,:,:,:);
    Rigid_Aerodata.Damping_Coeff.Cmr(1,:,:,:)=Rigid_AerodataDamping_CoeffCmr(1,:,:,:);
    Rigid_Aerodata.Damping_Coeff.Cnr(1,:,:,:)=Rigid_AerodataDamping_CoeffCnr(1,:,:,:);


if length(ac_state)==1
    Rigid_Aerodata.Static_Coeff.CX(2,:,:,:)=Rigid_Aerodata.Static_Coeff.CX(1,:,:,:);
    Rigid_Aerodata.Static_Coeff.CY(2,:,:,:)=Rigid_Aerodata.Static_Coeff.CY(1,:,:,:);
    Rigid_Aerodata.Static_Coeff.CZ(2,:,:,:)=Rigid_Aerodata.Static_Coeff.CZ(1,:,:,:);
    % static moment derivatives
    Rigid_Aerodata.Static_Coeff.Cl(2,:,:,:)=Rigid_Aerodata.Static_Coeff.Cl(1,:,:,:);
    Rigid_Aerodata.Static_Coeff.Cm(2,:,:,:)=Rigid_Aerodata.Static_Coeff.Cm(1,:,:,:);
    Rigid_Aerodata.Static_Coeff.Cn(2,:,:,:)=Rigid_Aerodata.Static_Coeff.Cn(1,:,:,:);
    % damping derivatives
    % p derivatives
    Rigid_Aerodata.Damping_Coeff.CXp(2,:,:,:)=Rigid_Aerodata.Damping_Coeff.CXp(1,:,:,:);
    Rigid_Aerodata.Damping_Coeff.CYp(2,:,:,:)=Rigid_Aerodata.Damping_Coeff.CYp(1,:,:,:);
    Rigid_Aerodata.Damping_Coeff.CZp(2,:,:,:)=Rigid_Aerodata.Damping_Coeff.CZp(1,:,:,:);
    Rigid_Aerodata.Damping_Coeff.Clp(2,:,:,:)=Rigid_Aerodata.Damping_Coeff.Clp(1,:,:,:);
    Rigid_Aerodata.Damping_Coeff.Cmp(2,:,:,:)=Rigid_Aerodata.Damping_Coeff.Cmp(1,:,:,:);
    Rigid_Aerodata.Damping_Coeff.Cnp(2,:,:,:)=Rigid_Aerodata.Damping_Coeff.Cnp(1,:,:,:);
    
    % q derivatives
    Rigid_Aerodata.Damping_Coeff.CXq(2,:,:,:)=Rigid_Aerodata.Damping_Coeff.CXq(1,:,:,:);%*cosd(alpha)+Cq(3)*sind(alpha);
    Rigid_Aerodata.Damping_Coeff.CYq(2,:,:,:)=Rigid_Aerodata.Damping_Coeff.CYq(1,:,:,:);
    Rigid_Aerodata.Damping_Coeff.CZq(2,:,:,:)=Rigid_Aerodata.Damping_Coeff.CZq(1,:,:,:);%*cosd(alpha)-Cq(1)*sind(alpha);
    Rigid_Aerodata.Damping_Coeff.Clq(2,:,:,:)=Rigid_Aerodata.Damping_Coeff.Clq(1,:,:,:);
    Rigid_Aerodata.Damping_Coeff.Cmq(2,:,:,:)=Rigid_Aerodata.Damping_Coeff.Cmq(1,:,:,:);
    Rigid_Aerodata.Damping_Coeff.Cnq(2,:,:,:)=Rigid_Aerodata.Damping_Coeff.Cnq(1,:,:,:);
    
    % p derivatives
    Rigid_Aerodata.Damping_Coeff.CXr(2,:,:,:)=Rigid_Aerodata.Damping_Coeff.CXr(1,:,:,:);
    Rigid_Aerodata.Damping_Coeff.CYr(2,:,:,:)=Rigid_Aerodata.Damping_Coeff.CYr(1,:,:,:);
    Rigid_Aerodata.Damping_Coeff.CZr(2,:,:,:)=Rigid_Aerodata.Damping_Coeff.CZr(1,:,:,:);
    Rigid_Aerodata.Damping_Coeff.Clr(2,:,:,:)=Rigid_Aerodata.Damping_Coeff.Clr(1,:,:,:);
    Rigid_Aerodata.Damping_Coeff.Cmr(2,:,:,:)=Rigid_Aerodata.Damping_Coeff.Cmr(1,:,:,:);
    Rigid_Aerodata.Damping_Coeff.Cnr(2,:,:,:)=Rigid_Aerodata.Damping_Coeff.Cnr(1,:,:,:);
end




for conf_cs=1:length(aircraft.control_surfaces)
    aircraft=aircraft.f_set_control_surface(aircraft.control_surfaces{conf_cs},0);
end

if deflected
    aircraft=aircraft.compute_grid();
    aircraft=aircraft.compute_deflected_grid(def);
else
    aircraft=aircraft.compute_grid();
end

%% Initializing the variables for parallel processing

                
for cs=1:length(control_surfs)
    if aircraft.control_surfaces{cs}(1:6) == 'ailero'
        cs_defl_grid = vecDeltaP;
    elseif aircraft.control_surfaces{cs}(1:6) == 'elevat'
        cs_defl_grid = vecDeltaQ;
    elseif aircraft.control_surfaces{cs}(1:6) == 'rudder'
        cs_defl_grid = vecDeltaR;
    else
        cs_defl_grid =vecCsDefl;
    end
    eval(['Rigid_AerodataControl_Surface_Coeff' control_surfs{cs} 'CX=zeros(length(vecConf),length(vecMa),length(vecAlpha),length(vecBeta),length(cs_defl_grid));']);
    eval(['Rigid_AerodataControl_Surface_Coeff' control_surfs{cs} 'CY=zeros(length(vecConf),length(vecMa),length(vecAlpha),length(vecBeta),length(cs_defl_grid));']);
    eval(['Rigid_AerodataControl_Surface_Coeff' control_surfs{cs} 'CZ=zeros(length(vecConf),length(vecMa),length(vecAlpha),length(vecBeta),length(cs_defl_grid));']);
    eval(['Rigid_AerodataControl_Surface_Coeff' control_surfs{cs} 'Cl=zeros(length(vecConf),length(vecMa),length(vecAlpha),length(vecBeta),length(cs_defl_grid));']);
    eval(['Rigid_AerodataControl_Surface_Coeff' control_surfs{cs} 'Cm=zeros(length(vecConf),length(vecMa),length(vecAlpha),length(vecBeta),length(cs_defl_grid));']);
    eval(['Rigid_AerodataControl_Surface_Coeff' control_surfs{cs} 'Cn=zeros(length(vecConf),length(vecMa),length(vecAlpha),length(vecBeta),length(cs_defl_grid));']);
end
                
% Calculation of control surface coefficients
lengthAlpha=length(vecAlpha);
lengthBeta=length(vecBeta);
lengthConf=length(vecConf);
lengthMa=length(vecMa);

for l=1:lengthConf
    conf=vecConf(l);
    if config~=1
        for conf_cs=1:length(config(conf).controls)
            aircraft=aircraft.f_set_control_surface(config(conf).controls{conf_cs},config(conf).deflections*0);
        end
    end
    if deflected
        aircraft=aircraft.compute_grid();
        aircraft=aircraft.compute_deflected_grid(def);
    else
        aircraft=aircraft.compute_grid();
    end
    
    parfor k=1:lengthMa
        Ma=vecMa(k)
        if Ma
        Uinf=Ma*a; 
        else
            Uinf=200;
        end%correction
        aircraftma=aircraft;
        for j=1:lengthAlpha
            alpha=vecAlpha(j);
            for i=1:lengthBeta
                beta=vecBeta(i)
                state=class_aero_state(Uinf,alpha,beta,Ma,rho_air);
                if deflected
                    wingaero=class_VLM_solver(aircraftma.grid_deflected,aircraftma.te_idx,aircraftma.panels,state,aircraftma.reference);
                else
                    wingaero=class_VLM_solver(aircraftma.grid,aircraftma.te_idx,aircraftma.panels,state,aircraftma.reference);
                end
                % wingaero=wingaero.f_set_state(state);
                %wingaero=wingaero.set_grid(aircraft.grid,aircraft.panels);
                
                wingaero=wingaero.f_solve_std();
                wingaero=wingaero.f_solve_full();
                if drag==0
                    wingaero.Cdi=0;
                end
                C0= [wingaero.CX, wingaero.CY, wingaero.CZ, wingaero.CL, wingaero.CM, wingaero.CN];
                lengthcs=length(control_surfs)
                for cs=1:lengthcs                 
                    %choose gridding for this control surface
                    if aircraftma.control_surfaces{cs}(1:6) == 'ailero'
                        cs_defl_grid = vecDeltaP;
                    elseif aircraftma.control_surfaces{cs}(1:6) == 'elevat'
                        cs_defl_grid = vecDeltaQ;
                    elseif aircraftma.control_surfaces{cs}(1:6) == 'rudder'
                        cs_defl_grid = vecDeltaR;
                    else
                        cs_defl_grid =vecCsDefl;
                    end
                    CX=zeros(length(vecConf),length(vecMa),length(vecAlpha),length(vecBeta),length(cs_defl_grid));
                    CY=zeros(length(vecConf),length(vecMa),length(vecAlpha),length(vecBeta),length(cs_defl_grid));
                    CZ=zeros(length(vecConf),length(vecMa),length(vecAlpha),length(vecBeta),length(cs_defl_grid));
                    Cl=zeros(length(vecConf),length(vecMa),length(vecAlpha),length(vecBeta),length(cs_defl_grid));
                    Cm=zeros(length(vecConf),length(vecMa),length(vecAlpha),length(vecBeta),length(cs_defl_grid));
                    Cn=zeros(length(vecConf),length(vecMa),length(vecAlpha),length(vecBeta),length(cs_defl_grid));
                    lengthcsdefl=length(cs_defl_grid)
                    for h=1:lengthcsdefl
                        cs_defl=cs_defl_grid(h);
                        disp(['Ma=' num2str(Ma) ' aoa=' num2str(alpha) ' beta=' num2str(beta) ' cs:' control_surfs{cs} ' csdefl=' num2str(cs_defl)])
                        aircraftma=aircraftma.f_set_control_surface(control_surfs{cs},cs_defl);
                        if deflected
                            aircraftma=aircraftma.compute_grid();
                            aircraftma=aircraftma.compute_deflected_grid(def);
                            wingaero=class_VLM_solver(aircraftma.grid_deflected,aircraftma.te_idx,aircraftma.panels,state,aircraftma.reference);
                        else
                            aircraftma=aircraftma.compute_grid();
                            wingaero=class_VLM_solver(aircraftma.grid,aircraftma.te_idx,aircraftma.panels,state,aircraftma.reference);
                        end
                        wingaero=wingaero.f_solve_std();
                        if drag==0
                            wingaero.Cdi=0;
                            aircraftma.CD_f=0;
                        end
                        %wingaero=wingaero.f_solve_full();
                        C= [wingaero.CX, wingaero.CY, wingaero.CZ, wingaero.CL, wingaero.CM, wingaero.CN]-C0;
                        % static force derivatives
                        CX(l,k,j,i,h)=-C(1);
                        CY(l,k,j,i,h)=C(2);
                        CZ(l,k,j,i,h)=-C(3);
                        % static moment derivatives
                        Cl(l,k,j,i,h)=-C(4);
                        Cm(l,k,j,i,h)=C(5);
                        Cn(l,k,j,i,h)=-C(6);                      
                    end
                    if aircraftma.control_surfaces{cs}(end-1:end) == 'ft'
                        Rigid_AerodataControl_Surface_Coeffaileron_leftCX(l,k,j,i,:)=CX(l,k,j,i,:);
                        Rigid_AerodataControl_Surface_Coeffaileron_leftCY(l,k,j,i,:)=CY(l,k,j,i,:);
                        Rigid_AerodataControl_Surface_Coeffaileron_leftCZ(l,k,j,i,:)=CZ(l,k,j,i,:);
                        Rigid_AerodataControl_Surface_Coeffaileron_leftCl(l,k,j,i,:)=Cl(l,k,j,i,:);
                        Rigid_AerodataControl_Surface_Coeffaileron_leftCm(l,k,j,i,:)=Cm(l,k,j,i,:);
                        Rigid_AerodataControl_Surface_Coeffaileron_leftCn(l,k,j,i,:)=Cn(l,k,j,i,:);
                    elseif aircraftma.control_surfaces{cs}(end-1:end) == 'or'
                        Rigid_AerodataControl_Surface_CoeffelevatorCX(l,k,j,i,:)=CX(l,k,j,i,:);
                        Rigid_AerodataControl_Surface_CoeffelevatorCY(l,k,j,i,:)=CY(l,k,j,i,:);
                        Rigid_AerodataControl_Surface_CoeffelevatorCZ(l,k,j,i,:)=CZ(l,k,j,i,:);
                        Rigid_AerodataControl_Surface_CoeffelevatorCl(l,k,j,i,:)=Cl(l,k,j,i,:);
                        Rigid_AerodataControl_Surface_CoeffelevatorCm(l,k,j,i,:)=Cm(l,k,j,i,:);
                        Rigid_AerodataControl_Surface_CoeffelevatorCn(l,k,j,i,:)=Cn(l,k,j,i,:);
                    elseif aircraftma.control_surfaces{cs}(end-1:end) == 'er'
                        Rigid_AerodataControl_Surface_CoeffrudderCX(l,k,j,i,:)=CX(l,k,j,i,:);
                        Rigid_AerodataControl_Surface_CoeffrudderCY(l,k,j,i,:)=CY(l,k,j,i,:);
                        Rigid_AerodataControl_Surface_CoeffrudderCZ(l,k,j,i,:)=CZ(l,k,j,i,:);
                        Rigid_AerodataControl_Surface_CoeffrudderCl(l,k,j,i,:)=Cl(l,k,j,i,:);
                        Rigid_AerodataControl_Surface_CoeffrudderCm(l,k,j,i,:)=Cm(l,k,j,i,:);
                        Rigid_AerodataControl_Surface_CoeffrudderCn(l,k,j,i,:)=Cn(l,k,j,i,:);
                    else
                        Rigid_AerodataControl_Surface_Coeffaileron_rightCX(l,k,j,i,:)=CX(l,k,j,i,:);
                        Rigid_AerodataControl_Surface_Coeffaileron_rightCY(l,k,j,i,:)=CY(l,k,j,i,:);
                        Rigid_AerodataControl_Surface_Coeffaileron_rightCZ(l,k,j,i,:)=CZ(l,k,j,i,:);
                        Rigid_AerodataControl_Surface_Coeffaileron_rightCl(l,k,j,i,:)=Cl(l,k,j,i,:);
                        Rigid_AerodataControl_Surface_Coeffaileron_rightCm(l,k,j,i,:)=Cm(l,k,j,i,:);
                        Rigid_AerodataControl_Surface_Coeffaileron_rightCn(l,k,j,i,:)=Cn(l,k,j,i,:);
                    end                                   
                    aircraftma=aircraftma.f_set_control_surface(control_surfs{cs},0);
                    if deflected
                        aircraftma=aircraftma.compute_grid();
                        aircraftma=aircraftma.compute_deflected_grid(def);
                    else
                        aircraftma=aircraftma.compute_grid();
                    end
                end
            end
        end
    end
end

%% Store in the table
for cs=1:length(control_surfs)
    eval(['Rigid_Aerodata.Control_Surface_Coeff.' control_surfs{cs} '.CX=Rigid_AerodataControl_Surface_Coeff' control_surfs{cs} 'CX;']);
    eval(['Rigid_Aerodata.Control_Surface_Coeff.' control_surfs{cs} '.CY=Rigid_AerodataControl_Surface_Coeff' control_surfs{cs} 'CY;']);
    eval(['Rigid_Aerodata.Control_Surface_Coeff.' control_surfs{cs} '.CZ=Rigid_AerodataControl_Surface_Coeff' control_surfs{cs} 'CZ;']);
    eval(['Rigid_Aerodata.Control_Surface_Coeff.' control_surfs{cs} '.Cl=Rigid_AerodataControl_Surface_Coeff' control_surfs{cs} 'Cl;']);
    eval(['Rigid_Aerodata.Control_Surface_Coeff.' control_surfs{cs} '.Cm=Rigid_AerodataControl_Surface_Coeff' control_surfs{cs} 'Cm;']);
    eval(['Rigid_Aerodata.Control_Surface_Coeff.' control_surfs{cs} '.Cn=Rigid_AerodataControl_Surface_Coeff' control_surfs{cs} 'Cn;']);
end
%%

if config==1
    k=1;
    l=2;
    for Ma=vecMa
        j=1; 
        for alpha=vecAlpha
            i=1;
            for beta=vecBeta
                h=1;
                
                for cs=1:length(control_surfs)
                    
                    eval(['Rigid_Aerodata.Control_Surface_Coeff.' control_surfs{cs} '.CX(l,k,j,i,:)=Rigid_Aerodata.Control_Surface_Coeff.' control_surfs{cs} '.CX(1,k,j,i,:);']);
                    eval(['Rigid_Aerodata.Control_Surface_Coeff.' control_surfs{cs} '.CY(l,k,j,i,:)=Rigid_Aerodata.Control_Surface_Coeff.' control_surfs{cs} '.CY(1,k,j,i,:);']);
                    eval(['Rigid_Aerodata.Control_Surface_Coeff.' control_surfs{cs} '.CZ(l,k,j,i,:)=Rigid_Aerodata.Control_Surface_Coeff.' control_surfs{cs} '.CZ(1,k,j,i,:);']);
                    eval(['Rigid_Aerodata.Control_Surface_Coeff.' control_surfs{cs} '.Cl(l,k,j,i,:)=Rigid_Aerodata.Control_Surface_Coeff.' control_surfs{cs} '.Cl(1,k,j,i,:);']);
                    eval(['Rigid_Aerodata.Control_Surface_Coeff.' control_surfs{cs} '.Cm(l,k,j,i,:)=Rigid_Aerodata.Control_Surface_Coeff.' control_surfs{cs} '.Cm(1,k,j,i,:);']);
                    eval(['Rigid_Aerodata.Control_Surface_Coeff.' control_surfs{cs} '.Cn(l,k,j,i,:)=Rigid_Aerodata.Control_Surface_Coeff.' control_surfs{cs} '.Cn(1,k,j,i,:);']);
                    h=1;
                end
                i=i+1;
            end
            j=j+1;
        end
        k=k+1;
        k
    end
end
time_taken=toc;
%% %% ground effect and landing gear dummy
Rigid_Aerodata.Landing_Gear.CX_gear=zeros(length(vecAlpha),1);
Rigid_Aerodata.Landing_Gear.CZ_gear=zeros(length(vecAlpha),1);

Rigid_Aerodata.Ground_Effect.CX_ground=zeros(length(vecConf),length(vecAlpha));
Rigid_Aerodata.Ground_Effect.CZ_ground=zeros(length(vecConf),length(vecAlpha));
Rigid_Aerodata.Ground_Effect.Cm_ground=zeros(length(vecConf),length(vecAlpha));

Xi_fact=zeros(40,1);
l=1;
landing_gear_down=1;
h=0;
for conf=vecConf
    i=1;
    for alpha=vecAlpha
        Rigid_Aerodata.Ground_Effect.CZ(l,i)=0;
        Rigid_Aerodata.Ground_Effect.CX(l,i)=0;
        Rigid_Aerodata.Ground_Effect.Cm(l,i)=0;
        i=i+1;
    end
    l=l+1;
end

if length(ac_state)==1
    i=1;
    for alpha=vecAlpha
        Rigid_Aerodata.Ground_Effect.CZ(2,i)= Rigid_Aerodata.Ground_Effect.CZ(1,i);
        Rigid_Aerodata.Ground_Effect.CX(2,i)=Rigid_Aerodata.Ground_Effect.CX(1,i);
        Rigid_Aerodata.Ground_Effect.Cm(2,i)=Rigid_Aerodata.Ground_Effect.Cm(1,i);
        i=i+1;
    end
end

i=1;
for alpha=vecAlpha
    Rigid_Aerodata.Landing_Gear.CZ(i)=0;
    Rigid_Aerodata.Landing_Gear.CX(i)=0;
    i=i+1;
end

for h=[0:1:40]
    Xi_fact(h+1)=0;
end
Rigid_Aerodata.Ground_Effect.Xi=Xi_fact;
Rigid_Aerodata.Ground_Effect.vecH=0:1:40;
Rigid_Aerodata.Landing_Gear.Xi=Xi_fact;
Rigid_Aerodata.Landing_Gear.vecH=0:1:40;

Rigid_Aerodata.Grid.vecAlpha=vecAlpha;
Rigid_Aerodata.Grid.vecMa=vecMa;
Rigid_Aerodata.Grid.vecBeta=vecBeta;
if config==1
    Rigid_Aerodata.Grid.vecConf=0:1;
else
    Rigid_Aerodata.Grid.vecConf=vecConf-1;
end
%% set cs grid (could be done within loop)
%%% TODO!!! automate!!! ? what
for cs=1:length(control_surfs)
    if aircraft.control_surfaces{cs}(1:6) == 'ailero'
        cs_defl_grid = vecDeltaP;
    elseif aircraft.control_surfaces{cs}(1:6) == 'elevat'
        cs_defl_grid = vecDeltaQ;
    elseif aircraft.control_surfaces{cs}(1:6) == 'rudder'
        cs_defl_grid = vecDeltaR;
    else
        cs_defl_grid =vecCsDefl;
    end
    eval(['Rigid_Aerodata.Control_Surface_Coeff.' control_surfs{cs} '.Grid= [' num2str(cs_defl_grid) '];']);
end
fprintf('Estimated computation time: %f hours \n',t_total/3600);
t_total=toc;
fprintf('Actual computation time: %f hours \n',t_total/3600);
% vecDeltaSP=vecDeltaR;
% %% save data
mkdir(aircraft.name)
save([aircraft.name '/Rigid_Aerodata'],'Rigid_Aerodata');

end

