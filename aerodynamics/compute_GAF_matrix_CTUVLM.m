%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%
function [Q,VGaf] = compute_GAF_matrix_CTUVLM(k,nmod,aircraft,aircraft_structure,state,name,UVLM_settings,exag,varargin)
% check if trim shape is given
if nargin==9
    shape=varargin{1};
else
    shape=0;
end
aeSSMSettings=AeSSMSettings;
if shape~=0
    aircraft_structure.nodal_deflections=shape;
    aircraft_structure=aircraft_structure.f_postprocess();
    def=aircraft_structure.f_get_deflections;
    aircraft=aircraft.compute_deflected_grid(def);
    %aircraft.grid=aircraft.grid_deflected;
else
    aeSSMSettings.deformedGrid=0;
end

aeSSMSettings.modes=7:nmod+6;
aeSSMSettings.loadsOut=0;
aeSSMSettings.mor='';
aeSSMSettings.poiDof=[];
aeSSMSettings.strRayleighDamping=0.000;
aeSSMSettings.strPropDamping=0.000;
aeSSMSettings.ctrInputs=1;

aeSSMF2=AeSSM(aircraft,aircraft_structure,state, aeSSMSettings);

aeSSMF2=aeSSMF2.generateTuvlm();
VGaf=state.V_A;
aerSSM=ss(   aeSSMF2.tUvlmSSM.a0+VGaf*aeSSMF2.tUvlmSSM.dadV,...
                        aeSSMF2.tUvlmSSM.b0+VGaf*aeSSMF2.tUvlmSSM.dbdV+VGaf^2*aeSSMF2.tUvlmSSM.dbdV2,...
                        aeSSMF2.tUvlmSSM.c0+VGaf*aeSSMF2.tUvlmSSM.dcdV,...
                        aeSSMF2.tUvlmSSM.d0+VGaf*aeSSMF2.tUvlmSSM.dddV+VGaf^2*aeSSMF2.tUvlmSSM.dddV2...
                        );
           aerSSM.InputName=aeSSMF2.tUvlmSSM.inputName;
           aerSSM.InputGroup=aeSSMF2.tUvlmSSM.inputGroup;
           aerSSM.OutputName=aeSSMF2.tUvlmSSM.outputName;
           aerSSM.OutputGroup=aeSSMF2.tUvlmSSM.outputGroup;
           

nModesGaf=nmod;
inputs=[aeSSMF2.uvlmSSM.inputName([ aeSSMF2.uvlmSSM.inputGroup.vBDot...
                                    aeSSMF2.uvlmSSM.inputGroup.vB...
                                    aeSSMF2.uvlmSSM.inputGroup.wBDot...
                                    aeSSMF2.uvlmSSM.inputGroup.wB...
                                    aeSSMF2.uvlmSSM.inputGroup.dCsDotDot...
                                    aeSSMF2.uvlmSSM.inputGroup.dCsDot...
                                    aeSSMF2.uvlmSSM.inputGroup.dCs; ]);...
        aeSSMF2.strSSM.outputName([ aeSSMF2.strSSM.outputGroup.modeDdot...
                                    aeSSMF2.strSSM.outputGroup.modeDot...
                                    aeSSMF2.strSSM.outputGroup.mode]);];
outputs=[{'X';'Y';'Z';'L';'M';'N'};cellstr([repmat('modalForce_',nModesGaf,1) num2str([1:nModesGaf]','%04d')]); cellstr([repmat('Mhinge_',length(aeSSMF2.uvlm.csNames),1) num2str([1:length(aeSSMF2.uvlm.csNames)]','%04d')])];

sparseSS=aerSSM(outputs,inputs);
% integrators for VbDot to vB
intSSM_T=ss(zeros(3),eye(3),[eye(3); zeros(3)],[zeros(3);eye(3)]);
intSSM_T.inputName={'TxDD';'TyDDT';'TzDDT';};
intSSM_T.OutputName={'vB_x';'vB_y';'vB_z';'vBDot_x';'vBDot_y';'vBDot_z';};
% wBdot to wB
intSSM_R1=ss(zeros(3),eye(3),[eye(3); zeros(3)],[zeros(3);eye(3)]);
intSSM_R1.inputName={'RxDD';'RyDD';'RzDD';};
intSSM_R1.OutputName={'wB_x';'wB_y';'wB_z';'wBDot_x';'wBDot_y';'wBDot_z';};
% wB_y to vBDot_z % pitch creates alpha which is an overlaying vertical body velocity
intSSM_R2=ss(VGaf);
intSSM_R2.InputName='wB_y';
intSSM_R2.OutputName='TzDDR';
% wB_y to vBDot_z % yaw creates beta which is an overlaying lateral body velocity
intSSM_R3=ss(VGaf);
intSSM_R3.InputName='wB_z';
intSSM_R3.OutputName='TyDDR';

intSSM_R=connect(intSSM_R1,intSSM_R2,intSSM_R3,{'RxDD';'RyDD';'RzDD';},{'wB_x';'wB_y';'wB_z';'wBDot_x';'wBDot_y';'wBDot_z';'TzDDR';'TyDDR'});
%add that to vBDot_zT
intSSM_vBzSum=sumblk('TzDDT=TzDDR+TzDD');
intSSM_vBySum=sumblk('TyDDT=TyDDR+TyDD');
intSSMtot=connect(intSSM_T,intSSM_R,intSSM_vBzSum,intSSM_vBySum,{'TxDD';'TyDD';'TzDD';'RxDD';'RyDD';'RzDD';},{'vB_x';'vB_y';'vB_z';'vBDot_x';'vBDot_y';'vBDot_z';'wB_x';'wB_y';'wB_z';'wBDot_x';'wBDot_y';'wBDot_z'});
% intSSMtot=intSSMtot({'vBDot_z';'vB_z';'wBDot_y';'wB_y'},{'TzDD';'RyDD'});
%integrator ssm for nModesGaf modes
% intSSMtot=[];
for iMode=1:nModesGaf
    intSSM=ss([0 0; 1 0],[1;0], [0 0; 1 0; 0 1],[1;0;0]);
    intSSM.OutputName={['modeDdot_'  num2str(iMode,'%04d')];['modeDot_'  num2str(iMode,'%04d')];['mode_'  num2str(iMode,'%04d')];};
    intSSM.InputName=['modeDdot_'  num2str(iMode,'%04d')];
    intSSMtot=append(intSSMtot,intSSM);
end

% for iCs=1:length(aeSSMF2.uvlm.csNames)
%     intSSM=ss([0 0; 1 0],[1;0], [0 0; 1 0; 0 1],[1;0;0]);
%     intSSM.OutputName={['dCsDotDot_'  num2str(iCs,'%04d')];['dCsDot_'  num2str(iCs,'%04d')];['dCs_'  num2str(iCs,'%04d')];};
%     intSSM.InputName=['dCsDotDot_'  num2str(iCs,'%04d')];
%     intSSMtot=append(intSSMtot,intSSM);
% end
%% combine symmetrically deflected inputs
intSSMCS=[];
for iWing=1:length(aircraft.wings)
    for iSeg=1:length(aircraft.wings(iWing).wing_segments)
         if aircraft.wings(iWing).wing_segments(iSeg).has_te_cs
            if and(aircraft.wings(iWing).wing_segments(iSeg).te_device.is_sym, aircraft.wings(iWing).symmetric)
                    iD_r=find(strcmp(aeSSMF2.uvlm.csNames,aircraft.wings(iWing).wing_segments(iSeg).te_device.name));
                    iD_l=find(strcmp(aeSSMF2.uvlm.csNames,[aircraft.wings(iWing).wing_segments(iSeg).te_device.name, '_l']));
                if aircraft.wings(iWing).wing_segments(iSeg).te_device.is_sym_defl
                    %add combined integrator
                    intSSM=ss([0 0; 1 0],[1;0], [0 0; 1 0; 0 1;0 0; -1 0; 0 -1],[1;0;0;-1;0;0]);
                    intSSM.OutputName={['dCsDotDot_'  num2str(iD_r,'%04d')];['dCsDot_'  num2str(iD_r,'%04d')];['dCs_'  num2str(iD_r,'%04d')];...
                                        ['dCsDotDot_'  num2str(iD_l,'%04d')];['dCsDot_'  num2str(iD_l,'%04d')];['dCs_'  num2str(iD_l,'%04d')];};
                    
                    intSSM.InputName=['dCsDotDot_'  num2str(iD_r,'%04d') num2str(iD_l,'%04d')];
                    intSSMCS=append(intSSMCS,intSSM);
                    % combine output
                    outSSM=sumblk(['Mhinge_'  num2str(iD_r,'%04d') num2str(iD_l,'%04d') '=' 'Mhinge_' num2str(iD_r,'%04d') '-Mhinge_' num2str(iD_l,'%04d') ]);
                    opt = connectOptions('Simplify',false);
                    outputs=[outputs(1:find(strcmp(outputs,'Mhinge_0003'))-1); ['Mhinge_'  num2str(iD_r,'%04d') num2str(iD_l,'%04d')]; outputs(find(strcmp(outputs,'Mhinge_0003'))+2:end)];
                    sparseSS=connect(sparseSS,outSSM,inputs,outputs, opt);
                else
                    %add two integrators
                    intSSM=ss([0 0; 1 0],[1;0], [0 0; 1 0; 0 1],[1;0;0]);
                    intSSM.OutputName={['dCsDotDot_'  num2str(iD_r,'%04d')];['dCsDot_'  num2str(iD_r,'%04d')];['dCs_'  num2str(iD_r,'%04d')];};
                    intSSM.InputName=['dCsDotDot_'  num2str(iD_r,'%04d')];
                    intSSMCS=append(intSSMCS,intSSM);
                    intSSM=ss([0 0; 1 0],[1;0], [0 0; 1 0; 0 1],[1;0;0]);
                    intSSM.OutputName={['dCsDotDot_'  num2str(iD_l,'%04d')];['dCsDot_'  num2str(iD_l,'%04d')];['dCs_'  num2str(iD_l,'%04d')];};
                    intSSM.InputName=['dCsDotDot_'  num2str(iD_l,'%04d')];
                    intSSMCS=append(intSSMCS,intSSM);
                end
            else
                % add single integrator
                iD_r=find(strcmp(aeSSMF2.uvlm.csNames,aircraft.wings(iWing).wing_segments(iSeg).te_device.name));
                intSSM=ss([0 0; 1 0],[1;0], [0 0; 1 0; 0 1],[1;0;0]);
                intSSM.OutputName={['dCsDotDot_'  num2str(iD_r,'%04d')];['dCsDot_'  num2str(iD_r,'%04d')];['dCs_'  num2str(iD_r,'%04d')];};
                intSSM.InputName=['dCsDotDot_'  num2str(iD_r,'%04d')];
                intSSMCS=append(intSSMCS,intSSM);
            end
         end
    end
end
 intSSMtot=append(intSSMtot,intSSMCS);
%% assemble

sparseSS=sss(series(intSSMtot,sparseSS,'name'));           
           
%%           
tic;
FreqVec=k*state.V_A*2/aircraft.reference.c_ref;
[vecMagData,vecPhaseData]=bode(sparseSS,FreqVec);
toc;           
           

vecMag=vecMagData.*-reshape(repmat(FreqVec.^2,size(vecMagData,1)*size(vecMagData,2),1),size(vecMagData,1),size(vecMagData,2),length(k));
vecPhase=(vecPhaseData)*pi/180;
QDaed = vecMag.*exp(vecPhase*sqrt(-1));           
%% norm gafs 
Q=QDaed/(.5*state.rho_air*VGaf^2*aircraft.reference.S_ref);           


% save([aircraft.name '/' name],'Q','k','aircraft_structure','aircraft','-v7.3');



end

