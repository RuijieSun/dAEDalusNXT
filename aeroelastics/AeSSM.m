%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%
classdef AeSSM
    %AESSM Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        % UVLM data struct for Simulink Model
        uvlmSSM
        % transformed UVLM 
        tUvlmSSM
        % Str data struct for Simulink Model
        strSSM
        % TAS data struct for Simulink Model (transformation from str to
        % aero displacements/boundary condition)
        tAS
        % TSA data struct for Simulink Model (transformation from aero to
        % str froces (modal))
        tSA
        % Settings (class:AeSSMSettings)
        settings@AeSSMSettings
        % linearized system
        linSSM
        % linearized system Velocity
        linV
        % V-g Data (Flutter Analysis)
        flutterTable
        % uvlm for ssm generation
        uvlm@class_UVLM_solver
        % strModel data needed (struct)
        strModel
        % state for initialization
        state@class_aero_state
        % reduction for aeroStates
        redBase
    end
    
    methods
        function obj=AeSSM(aircraft,aircraft_structure,aero_state,opt)
            
            %initialize settings or use input settings
            if nargin==4
                obj.settings=opt;
            else
                obj.settings=AeSSMSettings();
            end
            
            % set state
            obj.state=aero_state;
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % prepare aerodynamic model
            if ~obj.settings.deformedGrid
                aircraft_structure.nodal_deflections=aircraft_structure.nodal_deflections*0;
                aircraft_structure=aircraft_structure.f_postprocess();
            end
            aircraft=aircraft.compute_deflected_grid(aircraft_structure.f_get_deflections);
            %init uvlm settings
            uvlmSettings=class_UVLM_computation_settings();
            uvlmSettings.debug=0;
            uvlmSettings.movie=0;
            uvlmSettings.wakelength_factor=obj.settings.wakeLengthFactor/sqrt(1-obj.state.Ma^2);
            uvlmSettings.wakeGrowthFactor=obj.settings.rW;
            uvlmSettings.nFixedWakePanels=obj.settings.nFixedWakePanels;
            uvlmSettings.wakePanelDensity=obj.settings.wakePanelDensity;
            uvlmSettings.addLengthFactorLastPanel=obj.settings.addLengthFactorLastPanel;
            % compute cs surf panels
            aircraft=aircraft.computeControlSurfacePanelIds;
            % init uvlm
            obj.uvlm=class_UVLM_solver(aircraft,obj.state,uvlmSettings);
            obj.uvlm=obj.uvlm.generateSSM();
            obj.uvlm.ssm.settings.sectPress=0;
            flowAbs=norm(obj.state.V_inf);
            flowDir=obj.state.V_inf./flowAbs;
            V1=flowAbs;
            V2=(flowAbs-10);
            obj.uvlm.ssm=obj.uvlm.ssm.linearize(obj.uvlm,flowDir*V1);
            aeroSSM1=obj.uvlm.ssm.linSSM;
            obj.uvlm.ssm=obj.uvlm.ssm.linearize(obj.uvlm,flowDir*V2);
            aeroSSM2=obj.uvlm.ssm.linSSM;
            usedInputs=[aeroSSM1.inputGroup.b aeroSSM1.inputGroup.bDot];
            usedOutputs=[aeroSSM1.OutputGroup.Fp];
            if obj.settings.rbmInputs
                usedInputs=[usedInputs aeroSSM1.inputGroup.vB aeroSSM1.inputGroup.vBDot aeroSSM1.inputGroup.wB aeroSSM1.inputGroup.wBDot];
                usedOutputs=[usedOutputs aeroSSM1.OutputGroup.RBMForces aeroSSM1.OutputGroup.RBMMoments];
            end
            if obj.settings.gustInputs
                usedInputs=[usedInputs aeroSSM1.inputGroup.vG  aeroSSM1.inputGroup.vGDot];
            end
            if obj.settings.ctrInputs
                usedInputs=[usedInputs aeroSSM1.inputGroup.dCs  aeroSSM1.inputGroup.dCsDot  aeroSSM1.inputGroup.dCsDotDot];
                usedOutputs=[usedOutputs aeroSSM1.OutputGroup.HingeMoments];
            end
            if obj.settings.loadsOut
                %here future load outputs
            end
                
            % proceed with smaller aeroSSMs
            aeroSSM1=aeroSSM1(usedOutputs,usedInputs);
            aeroSSM2=aeroSSM2(usedOutputs,usedInputs);
            % velocity derivatives (TODO: directly create velocity sensitivities within
            % uvlm.generatesSSM() instead of "finite diff"
            obj.uvlmSSM.dadV=(aeroSSM1.a-aeroSSM2.a)./(V1-V2);
            obj.uvlmSSM.dbdV=(aeroSSM1.b-aeroSSM2.b)./(V1-V2);
            obj.uvlmSSM.dcdV=(aeroSSM1.c-aeroSSM2.c)./(V1-V2);
            obj.uvlmSSM.dddV=(aeroSSM1.d-aeroSSM2.d)./(V1-V2);
            % zero values
            obj.uvlmSSM.a0=aeroSSM1.a-obj.uvlmSSM.dadV*V1;
            obj.uvlmSSM.b0=aeroSSM1.b-obj.uvlmSSM.dbdV*V1;
            obj.uvlmSSM.c0=aeroSSM1.c-obj.uvlmSSM.dcdV*V1;
            obj.uvlmSSM.d0=aeroSSM1.d-obj.uvlmSSM.dddV*V1;
            
            % specification of input/output and statenames as well as in
            % and outputGroups
            obj.uvlmSSM.inputName=aeroSSM1.InputName;
            obj.uvlmSSM.inputGroup=aeroSSM1.InputGroup;
            obj.uvlmSSM.stateName=aeroSSM1.StateName;
            obj.uvlmSSM.outputGroup=aeroSSM1.OutputGroup;
            obj.uvlmSSM.outputName=aeroSSM1.OutputName;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % prepare structural model
            obj.strModel.modalBase=aircraft_structure.modeshapes(:,obj.settings.modes);
            obj.strModel.M=diag(diag(obj.strModel.modalBase'*aircraft_structure.Mff*obj.strModel.modalBase));
            obj.strModel.K=diag(diag(obj.strModel.modalBase'*aircraft_structure.Kff*obj.strModel.modalBase));
            if obj.settings.strRayleighDamping~=0
                obj.strModel.C=diag((obj.settings.strRayleighDamping*2*pi*aircraft_structure.modefrequencies(obj.settings.modes)));
            elseif obj.settings.strPropDamping~=0
                obj.strModel.C=diag(obj.settings.strPropDamping*(2*pi*aircraft_structure.modefrequencies(obj.settings.modes)));
            else
                obj.strModel.C=diag(0*(2*pi*aircraft_structure.modefrequencies(obj.settings.modes)));
            end
            nModes=length(obj.settings.modes);

            % init structural state space
            obj.strSSM.a=[ -obj.strModel.M^-1*obj.strModel.C, -obj.strModel.M^-1*obj.strModel.K; eye(nModes) zeros(nModes)];
            obj.strSSM.b=[obj.strModel.M^-1; zeros(nModes)];
            Cmodes=[obj.strSSM.a; zeros(nModes) eye(nModes)];
            Dmodes=[obj.strSSM.b; zeros(nModes)];
            poiBase=aircraft_structure.modeshapes(obj.settings.poiDof,obj.settings.modes);
            poiBase=blkdiag(poiBase,poiBase,poiBase);
            Cpoi=poiBase*Cmodes;
            Dpoi=poiBase*Dmodes;
            obj.strSSM.c=[Cmodes; Cpoi];
            obj.strSSM.d=[Dmodes; Dpoi];
            nPoi=length(obj.settings.poiDof);
            % specification of input/output and statenames as well as in
            % and outputGroups
            obj.strSSM.stateName=[  cellstr([repmat('modeDot_',nModes,1) num2str([1:nModes]','%04d')]);...
                                    cellstr([repmat('mode_',nModes,1) num2str([1:nModes]','%04d')])];
            obj.strSSM.inputName=cellstr([repmat('modalForce_',nModes,1) num2str([1:nModes]','%04d')]);
            obj.strSSM.outputName=[  cellstr([repmat('modeDdot_',nModes,1) num2str([1:nModes]','%04d')]);...
                                     cellstr([repmat('modeDot_',nModes,1) num2str([1:nModes]','%04d')]);...
                                     cellstr([repmat('mode_',nModes,1) num2str([1:nModes]','%04d')]);...
                                     cellstr([repmat('poiDdot_',nPoi,1) num2str([1:nPoi]','%04d')]);...
                                     cellstr([repmat('poiDot_',nPoi,1) num2str([1:nPoi]','%04d')]);...
                                     cellstr([repmat('poi_',nPoi,1) num2str([1:nPoi]','%04d')]);];
            obj.strSSM.inputGroup.modalForce=1:nModes;
            obj.strSSM.outputGroup.modeDdot=1:nModes;
            obj.strSSM.outputGroup.modeDot=nModes+1:2*nModes;
            obj.strSSM.outputGroup.mode=2*nModes+1:3*nModes;
            obj.strSSM.outputGroup.poiDdot=3*nModes+1:3*nModes+nPoi;
            obj.strSSM.outputGroup.poiDot=3*nModes+nPoi+1:3*nModes+2*nPoi;
            obj.strSSM.outputGroup.poi=3*nModes+2*nPoi+1:3*nModes+3*nPoi;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % prepare coupling matrices TAS and TSA
            nB=size(obj.uvlm.ssm.linSSM.InputGroup.b,2);

            if obj.settings.deformedGrid
                initDef=aircraft_structure.nodal_deflections; %<- CAREFULL!
            else
                initDef=aircraft_structure.nodal_deflections*0;
            end
            
            uvlmDef=obj.uvlm;
            incr=aircraft.reference.b_ref/100;
            % dBdMode (is linear function of V) describes e.g. influence of mode on b and modeDot on bDot
            % because of the linear dependency on V, only ddBBdModedV is determined and
            % later multiplied by the velocity (see TAS assembly)
            ddBdModedV=zeros(nB,nModes);
            % dBdModeDot - describes the influence of modeDot on b and modeDotDot on bDot
            dBdModeDot=zeros(nB,nModes);
            dCollocdModeXYZ=zeros(nB*3,nModes);
            dFvapdModeXYZ=zeros(nB*3,nModes);
            collocNvec=obj.uvlm.colloc_nvec;
            colloc=obj.uvlm.colloc;
            fvap=obj.uvlm.fvap;
            bZero=flowDir*collocNvec;
            for iMode=1:nModes
                aircraft_structure.nodal_deflections=initDef+incr*aircraft_structure.modeshapes(:,obj.settings.modes(iMode));
            %     aircraft_structure_def.nodal_deflections=initDef+incr*aircraft_structure.modeshapes(:,iMode);
                aircraft_structure=aircraft_structure.f_postprocess();
                aircraftDef=aircraft.compute_deflected_grid(aircraft_structure.f_get_deflections);
                uvlmDef=uvlmDef.set_grid(aircraftDef.grid_deflected,aircraftDef.panels);
                collocNvecMode=uvlmDef.colloc_nvec;
                collocMode=uvlmDef.colloc;
                fvapMode=uvlmDef.fvap;
                bMode=flowDir*collocNvecMode;
                ddBdModedV(:,iMode)=(bMode-bZero)/incr;
                dCollocdMode=(collocMode-colloc)/incr;
                dFvapdMode=(fvapMode-fvap)/incr;
                dBdModeDot(:,iMode)=sum((dCollocdMode.*collocNvec),1);
                dCollocdModeXYZ(:,iMode)=reshape(dCollocdMode,nB*3,1);
                dFvapdModeXYZ(:,iMode)=reshape(dFvapdMode,nB*3,1);
            end
        %TAS data
            obj.tAS.d0=[  zeros(nB,nModes)   -dBdModeDot      zeros(nB,nModes);...
                        -dBdModeDot          zeros(nB,nModes) zeros(nB,nModes)];
            obj.tAS.dddV=[  zeros(nB,nModes)   zeros(nB,nModes)      ddBdModedV;...
                            zeros(nB,nModes)       ddBdModedV zeros(nB,nModes)];
            obj.tAS.inputName=[	cellstr([repmat('modeDdot_',nModes,1) num2str([1:nModes]','%04d')]);...
                                cellstr([repmat('modeDot_',nModes,1) num2str([1:nModes]','%04d')]);...
                                cellstr([repmat('mode_',nModes,1) num2str([1:nModes]','%04d')])];
            obj.tAS.inputGroup.modeDdot=1:nModes;
            obj.tAS.inputGroup.modeDot=nModes+1:2*nModes;
            obj.tAS.inputGroup.mode=2*nModes+1:3*nModes;
            obj.tAS.outputName=aeroSSM1.inputName([aeroSSM1.InputGroup.b aeroSSM1.InputGroup.bDot]); 
            obj.tAS.outputGroup.b=1:nB;
            obj.tAS.outputGroup.bDot=nB+1:2*nB;        

            %TSA data
            obj.tSA.d=dFvapdModeXYZ';
            obj.tSA.inputName=aeroSSM1.OutputName(aeroSSM1.OutputGroup.Fp);
            obj.tSA.inputGroup.Fp=1:length(aeroSSM1.OutputGroup.Fp);
            obj.tSA.outputName=cellstr([repmat('modalForce_',nModes,1) num2str([1:nModes]','%04d')]);
            obj.tSA.outputGroup.modalForce=1:nModes;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % reduction of TAS Aero TSA system with specified reduction
%             obj=obj.redAero(V1);
            

        end % constructor
        
        function obj=redAero(obj,V1,varargin) %<- this should be in aero ssm class
            if nargin==3
                obj.redBase=varargin{1};
            else
                aerSSM=ss(   obj.uvlmSSM.a0+obj.uvlmSSM.dadV*V1,...
                            obj.uvlmSSM.b0+obj.uvlmSSM.dbdV*V1,...
                            obj.uvlmSSM.c0+obj.uvlmSSM.dcdV*V1,...
                            obj.uvlmSSM.d0+obj.uvlmSSM.dddV*V1);
                aerSSM.inputName=obj.uvlmSSM.inputName;
                aerSSM.InputGroup=obj.uvlmSSM.inputGroup;
                aerSSM.outputName=obj.uvlmSSM.outputName;
                aerSSM.OutputGroup=obj.uvlmSSM.outputGroup;
                aerSSM.StateName=obj.uvlmSSM.stateName;

                tASSSM=ss(obj.tAS.d0+obj.tAS.dddV*V1);
                tASSSM.inputName=obj.tAS.inputName;
                tASSSM.InputGroup=obj.tAS.inputGroup;
                tASSSM.outputName=obj.tAS.outputName;
                tASSSM.OutputGroup=obj.tAS.outputGroup;
                aerSSM2=seriesPreserve(tASSSM,aerSSM,{'b','bDot'});
                tSASSM=ss(obj.tSA.d);
                tSASSM.inputName=obj.tSA.inputName;
                tSASSM.InputGroup=obj.tSA.inputGroup;
                tSASSM.outputName=obj.tSA.outputName;
                tSASSM.OutputGroup=obj.tSA.outputGroup;
                aerSSM3=seriesPreserve(aerSSM2,tSASSM,{'Fp'});
                
%                  aerSSM3=aerSSM3(:,[aerSSM3.InputGroup.modeDdot aerSSM3.InputGroup.dCsDotDot aerSSM3.InputGroup.vBDot aerSSM3.InputGroup.wBDot aerSSM3.InputGroup.vGDot]);
%                  aerSSM3=aerSSM3(:,[aerSSM3.InputGroup.modeDdot aerSSM3.InputGroup.modeDot aerSSM3.InputGroup.mode]);
                
                
%                 %assemble system for redbase creation
%                 % connect TAS to uvlm and by that reduce input vector by
%                 % eliminating b and bDot inputs and introducing modal
%                 % inputs; this is done by transforming B and D
%                 %init b
%                 bpodSSMb=obj.uvlmSSM.b0+obj.uvlmSSM.dbdV*V1;
%                 %without b and bDot
%                 bpodSSMb(:,[obj.uvlmSSM.inputGroup.b obj.uvlmSSM.inputGroup.bDot])=[];
%                 %append transformed columns
%                 bpodSSMb=[bpodSSMb  (obj.uvlmSSM.b0(:,[obj.uvlmSSM.inputGroup.b obj.uvlmSSM.inputGroup.bDot])+obj.uvlmSSM.dbdV(:,[obj.uvlmSSM.inputGroup.b obj.uvlmSSM.inputGroup.bDot])*V1)*(obj.tAS.d0+obj.tAS.dddV*V1)];
% 
%                 %init d
%                 bpodSSMd=obj.uvlmSSM.d0+obj.uvlmSSM.dddV*V1;
%                 %without b and bDot
%                 bpodSSMd(:,[obj.uvlmSSM.inputGroup.b obj.uvlmSSM.inputGroup.bDot])=[];
%                 %append transformed columns
%                 bpodSSMd=[bpodSSMd  (obj.uvlmSSM.d0(:,[obj.uvlmSSM.inputGroup.b obj.uvlmSSM.inputGroup.bDot])+obj.uvlmSSM.dddV(:,[obj.uvlmSSM.inputGroup.b obj.uvlmSSM.inputGroup.bDot])*V1)*(obj.tAS.d0+obj.tAS.dddV*V1)];
% 
%                 %to reduce the outputs, c and d are transformed with tSA.d
%                 bpodSSMc=obj.uvlmSSM.c0+obj.uvlmSSM.dcdV*V1;
%                 bpodSSMcPartToTransform=bpodSSMc(obj.uvlmSSM.outputGroup.Fp,:);
%                 %cut out part which has to be transformed
%                 bpodSSMc(obj.uvlmSSM.outputGroup.Fp,:)=[];
%                 %append transformed part
%                 bpodSSMc=[bpodSSMc; obj.tSA.d*bpodSSMcPartToTransform];
% 
%                 %d was already initialized
%                 bpodSSMdPartToTransform=bpodSSMd(obj.uvlmSSM.outputGroup.Fp,:);
%                 %cut out part which has to be transformed
%                 bpodSSMd(obj.uvlmSSM.outputGroup.Fp,:)=[];
%                 %append transformed part
%                 bpodSSMd=[bpodSSMd; obj.tSA.d*bpodSSMdPartToTransform];
% 
% 
%                 bpodSSM=ss(obj.uvlmSSM.a0+obj.uvlmSSM.dadV*V1,bpodSSMb,bpodSSMc,bpodSSMd);
%                 % only use modal in/outputs for red:
%                 bpodSSM_small=ss(obj.uvlmSSM.a0+obj.uvlmSSM.dadV*V1,bpodSSMb_small,bpodSSMc_small,bpodSSMd_small2);
%                 %add integrators to all dot inputs....

%                 bpodOrder=size(aerSSM3,2); %<- if chosen too large, assembled ssm may have near zero unstable poles ! should match the order of aerodynamic inputs used for bpod.
                if strcmp(obj.settings.mor,'bPod')
                    tic;
                    [~,obj.redBase]=bpod(aerSSM3,obj.settings.morOrder, obj.settings.bpodT, obj.settings.nPODModes);
                    toc;
                elseif strcmp(obj.settings.mor,'balRed')
                    disp('balRed not implemented yet');
                    obj.redBase=eye(obj.uvlmSSM.a0);
                else
                    obj.redBase=eye(size(obj.uvlmSSM.a0));
                end
            end
            %transform aero with redBase and save only transformed
            %values.
%             obj.uvlmSSM.a0=obj.redBase'*obj.uvlmSSM.a0*obj.redBase;
%             obj.uvlmSSM.dadV=obj.redBase'*obj.uvlmSSM.dadV*obj.redBase;
%             obj.uvlmSSM.b0=obj.redBase'*obj.uvlmSSM.b0;
%             obj.uvlmSSM.dbdV=obj.redBase'*obj.uvlmSSM.dbdV;
%             obj.uvlmSSM.c0=obj.uvlmSSM.c0*obj.redBase;
%             obj.uvlmSSM.dcdV=obj.uvlmSSM.dcdV*obj.redBase;
            obj.uvlmSSM.a0=pinv(obj.redBase)*obj.uvlmSSM.a0*obj.redBase;
            obj.uvlmSSM.dadV=pinv(obj.redBase)*obj.uvlmSSM.dadV*obj.redBase;
            obj.uvlmSSM.b0=pinv(obj.redBase)*obj.uvlmSSM.b0;
            obj.uvlmSSM.dbdV=pinv(obj.redBase)*obj.uvlmSSM.dbdV;
            obj.uvlmSSM.c0=obj.uvlmSSM.c0*obj.redBase;
            obj.uvlmSSM.dcdV=obj.uvlmSSM.dcdV*obj.redBase;
            obj.uvlmSSM.stateName=cellstr([repmat('LagRed_',size(obj.redBase,2),1) num2str([1:size(obj.redBase,2)]','%04d')]);
        end
        
        function obj=getLinSSM(obj, Vlin)
           aerSSM=ss(   obj.tUvlmSSM.a0+Vlin*obj.tUvlmSSM.dadV,...
                        obj.tUvlmSSM.b0+Vlin*obj.tUvlmSSM.dbdV+Vlin^2*obj.tUvlmSSM.dbdV2,...
                        obj.tUvlmSSM.c0+Vlin*obj.tUvlmSSM.dcdV,...
                        obj.tUvlmSSM.d0+Vlin*obj.tUvlmSSM.dddV+Vlin^2*obj.tUvlmSSM.dddV2...
                        );
           aerSSM.InputName=obj.tUvlmSSM.inputName;
           aerSSM.InputGroup=obj.tUvlmSSM.inputGroup;
           aerSSM.OutputName=obj.tUvlmSSM.outputName;
           aerSSM.OutputGroup=obj.tUvlmSSM.outputGroup;

           structSSM=ss(obj.strSSM.a,obj.strSSM.b,obj.strSSM.c,obj.strSSM.d);
           structSSM.InputName=obj.strSSM.inputName;
           structSSM.OutputName=obj.strSSM.outputName;
           structSSM.InputGroup=obj.strSSM.inputGroup;
           structSSM.OutputGroup=obj.strSSM.outputGroup;
           
           obj.linSSM=connect(aerSSM,structSSM,[aerSSM.InputName],[aerSSM.outputName; structSSM.outputName]);
            
        end %linearize
        
        function obj=generateTuvlm(obj)
            % add in and outputs to tAS.d and tSA.d in order to connect it to uvlmSSM with
            % preserving all in and outputs of the uvlmSSM 
            conIdsInUvlm=[obj.uvlmSSM.inputGroup.b obj.uvlmSSM.inputGroup.bDot];
            conIdsOutUvlm=obj.uvlmSSM.outputGroup.Fp;
            nInUvlm=length(obj.uvlmSSM.inputName);
            nInUvlmPrv=min(conIdsInUvlm)-1;
            nInUvlmPost=nInUvlm-nInUvlmPrv-length(conIdsInUvlm);
            nOutUvlm=length(obj.uvlmSSM.outputName);
            nOutUvlmPrv=min(conIdsOutUvlm)-1; %<-------
            nOutUvlmPost=nOutUvlm-nOutUvlmPrv-length(conIdsOutUvlm);
            % expanded tAS
            tASe0=[eye(nInUvlm) [zeros(nInUvlmPrv,size(obj.tAS.d0,2)); obj.tAS.d0 ; zeros(nInUvlmPost,size(obj.tAS.d0,2))]];
            dtASedV=[zeros(nInUvlm) [zeros(nInUvlmPrv,size(obj.tAS.d0,2)); obj.tAS.dddV ; zeros(nInUvlmPost,size(obj.tAS.d0,2))]];
            %resulting inputvector
            nameIn=[obj.uvlmSSM.inputName;obj.tAS.inputName];
            pairs = [fieldnames(obj.uvlmSSM.inputGroup), struct2cell(obj.uvlmSSM.inputGroup); fieldnames(obj.tAS.inputGroup), struct2cell(obj.tAS.inputGroup)].';
            groupIn=struct(pairs{:});
            
            %discard b and bDot as input
            tASe0(:,conIdsInUvlm)=[];
            dtASedV(:,conIdsInUvlm)=[];
            nameIn(conIdsInUvlm)=[];
            groupIn=rmfield(groupIn,'b');
            groupIn=rmfield(groupIn,'bDot');
            
            % expanded tSA
            tSAe=[eye(nOutUvlm); zeros(size(obj.tSA.d,1),nOutUvlmPrv) obj.tSA.d zeros(size(obj.tSA.d,1),nOutUvlmPost)];
            %resulting outputvector
            nameOut=[obj.uvlmSSM.outputName;obj.tSA.outputName];
            pairs = [fieldnames(obj.uvlmSSM.outputGroup), struct2cell(obj.uvlmSSM.outputGroup); fieldnames(obj.tSA.outputGroup), struct2cell(obj.tSA.outputGroup)].';
            groupOut= struct(pairs{:});
            
            %discard Fp as output
            tSAe(conIdsOutUvlm,:)=[];
            nameOut(conIdsOutUvlm)=[];
            groupOut=rmfield(groupOut,'Fp');

            % analytically connect in series with uvlmSSM to obtain one ssm with having
            % b and d as a quadratic function of V called transformed uvlm ssm
            % (tUvlmSSM)

            obj.tUvlmSSM.a0=obj.uvlmSSM.a0;
            obj.tUvlmSSM.dadV=obj.uvlmSSM.dadV;

            obj.tUvlmSSM.b0=obj.uvlmSSM.b0*tASe0;
            obj.tUvlmSSM.dbdV=obj.uvlmSSM.b0*dtASedV+obj.uvlmSSM.dbdV*tASe0;
            obj.tUvlmSSM.dbdV2=obj.uvlmSSM.dbdV*dtASedV;

            obj.tUvlmSSM.c0=tSAe*obj.uvlmSSM.c0;
            obj.tUvlmSSM.dcdV=tSAe*obj.uvlmSSM.dcdV;

            obj.tUvlmSSM.d0=tSAe*obj.uvlmSSM.d0*tASe0;
            obj.tUvlmSSM.dddV=tSAe*(obj.uvlmSSM.d0*dtASedV+obj.uvlmSSM.dddV*tASe0);
            obj.tUvlmSSM.dddV2=tSAe*obj.uvlmSSM.dddV*dtASedV;
            obj.tUvlmSSM.inputName=nameIn;
            obj.tUvlmSSM.outputName=nameOut;
            
            
            fieldNames=fieldnames(groupOut);
            for iField=1:size(fieldNames,1)
                if isfield(obj.uvlmSSM.outputGroup, fieldNames{iField})
                    groupOut.(fieldNames{iField})=find(ismember(obj.tUvlmSSM.outputName, obj.uvlmSSM.outputName(obj.uvlmSSM.outputGroup.(fieldNames{iField}))))';
                elseif isfield(obj.tSA.outputGroup, fieldNames{iField})
                    groupOut.(fieldNames{iField})=find(ismember(obj.tUvlmSSM.outputName, obj.tSA.outputName(obj.tSA.outputGroup.(fieldNames{iField}))))';
                else
                    disp('STOP');
                end
            end
            
            fieldNames=fieldnames(groupIn);
            for iField=1:size(fieldNames,1)
                if isfield(obj.uvlmSSM.inputGroup, fieldNames{iField})
                    groupIn.(fieldNames{iField})=find(ismember(obj.tUvlmSSM.inputName, obj.uvlmSSM.inputName(obj.uvlmSSM.inputGroup.(fieldNames{iField}))))';
                elseif isfield(obj.tAS.inputGroup, fieldNames{iField})
                    groupIn.(fieldNames{iField})=find(ismember(obj.tUvlmSSM.inputName, obj.tAS.inputName(obj.tAS.inputGroup.(fieldNames{iField}))))';
                else
                    disp('STOP');
                end
            end
            
             obj.tUvlmSSM.outputGroup=groupOut;
             obj.tUvlmSSM.inputGroup=groupIn;
            
            obj.tUvlmSSM.stateName=obj.uvlmSSM.stateName;
        end
        
        function obj=runFlutterAnalysis(obj,VloopVec)
            obj.flutterTable=[];
            fFlag=0;
%             VloopVec=[1 5:5:55 60:10:Vmax*1.2];
            D3=obj.tSA.d;
            A4=obj.strSSM.a;
            B4=obj.strSSM.b;
            C4=obj.strSSM.c([obj.strSSM.outputGroup.modeDdot obj.strSSM.outputGroup.modeDot obj.strSSM.outputGroup.mode],:);
            D4=obj.strSSM.d([obj.strSSM.outputGroup.modeDdot obj.strSSM.outputGroup.modeDot obj.strSSM.outputGroup.mode],:);
            for Vloop=VloopVec
                % create aeroSSMV and TASV for certain speed
                A2=obj.uvlmSSM.a0+obj.uvlmSSM.dadV*Vloop;
                B2=obj.uvlmSSM.b0(:,[obj.uvlmSSM.inputGroup.b obj.uvlmSSM.inputGroup.bDot])+obj.uvlmSSM.dbdV(:,[obj.uvlmSSM.inputGroup.b obj.uvlmSSM.inputGroup.bDot])*Vloop;
                C2=obj.uvlmSSM.c0(obj.uvlmSSM.outputGroup.Fp,:)+obj.uvlmSSM.dcdV(obj.uvlmSSM.outputGroup.Fp,:)*Vloop;
                D2=obj.uvlmSSM.d0(obj.uvlmSSM.outputGroup.Fp,[obj.uvlmSSM.inputGroup.b obj.uvlmSSM.inputGroup.bDot])+obj.uvlmSSM.dddV(obj.uvlmSSM.outputGroup.Fp,[obj.uvlmSSM.inputGroup.b obj.uvlmSSM.inputGroup.bDot])*Vloop;
            %     aeroSSMV=ss(aeroSSMa,aeroSSMb,aeroSSMc,aeroSSMd);
                D1=obj.tAS.d0+obj.tAS.dddV*Vloop;
                N=(eye(size(D4,1))+D4*D3*D2*D1)^-1;
            %     A=[A2-B2*D1*D4*D3*C2 -B2*D1*C4;B4*D3*C2-B4*D3*D2*D1*D4*D3*C2 A4-B4*D3*D2*D1*C4];
            %     A=[redBase*A2*pinv(redBase)+redBase*B2*D1*D4*D3*C2*pinv(redBase) +redBase*B2*D1*C4;B4*D3*C2*pinv(redBase)+B4*D3*D2*D1*D4*D3*C2*pinv(redBase) A4+B4*D3*D2*D1*C4];
                A=[A2+B2*D1*N*D4*D3*C2 +B2*D1*N*C4;B4*D3*C2+B4*D3*D2*D1*N*D4*D3*C2 A4+B4*D3*D2*D1*N*C4];
                [~, eigVal]=eig(A);
                obj.flutterTable=[obj.flutterTable diag(eigVal)];
                hold on;  plot(real(diag(eigVal)),imag(diag(eigVal)),'.','MarkerFaceColor',[Vloop/max(VloopVec)*0.5 1-Vloop/max(VloopVec) 0],'MarkerEdgeColor',[Vloop/max(VloopVec)*0.5  1-Vloop/max(VloopVec) 0])
%               hold on;  plot(real(diag(eigVal)),imag(diag(eigVal)),'b+')
                if(any(real(diag(eigVal))>10^-10))
                    idx=find(real(diag(eigVal))>10^-10);
                    Freq=abs(eigVal(idx,idx));
                    fprintf('unstable pole @ %.02f m/s & %.02f Hz\n',Vloop,Freq(1,1)/2/pi) 
                    for iUnstable=1:length(diag(Freq));
                        if fFlag==0
                            plot(real(diag(eigVal(idx(iUnstable),idx(iUnstable)))),imag(diag(eigVal(idx(iUnstable),idx(iUnstable)))),'ro')
                            text(real(diag(eigVal(idx(iUnstable),idx(iUnstable)))),imag(diag(eigVal(idx(iUnstable),idx(iUnstable)))), sprintf('first @  %.02f m/s & %.02f Hz\n',Vloop,Freq(iUnstable,iUnstable)/2/pi));
                        end
                    end
                    fFlag=1;
                end

            end
        end % flutter
        function obj=runDivergenceAnalysis(obj,VloopVec) 
            TSA=obj.tSA.d;
            
            A4=obj.strSSM.a;
            B4=obj.strSSM.b;
            C4=obj.strSSM.c([obj.strSSM.outputGroup.mode],:);
            D4=obj.strSSM.d([obj.strSSM.outputGroup.mode],:);
            Kdcs=D4-C4*A4^-1*B4;
            figure; hold on;
            plot(VloopVec,ones(1,length(VloopVec)),'red');
            for Vloop=VloopVec 
                A2=obj.uvlmSSM.a0+obj.uvlmSSM.dadV*Vloop;
                B2=obj.uvlmSSM.b0(:,[obj.uvlmSSM.inputGroup.b obj.uvlmSSM.inputGroup.bDot])+obj.uvlmSSM.dbdV(:,[obj.uvlmSSM.inputGroup.b obj.uvlmSSM.inputGroup.bDot])*Vloop;
                C2=obj.uvlmSSM.c0(obj.uvlmSSM.outputGroup.Fp,:)+obj.uvlmSSM.dcdV(obj.uvlmSSM.outputGroup.Fp,:)*Vloop;
                D2=obj.uvlmSSM.d0(obj.uvlmSSM.outputGroup.Fp,[obj.uvlmSSM.inputGroup.b obj.uvlmSSM.inputGroup.bDot])+obj.uvlmSSM.dddV(obj.uvlmSSM.outputGroup.Fp,[obj.uvlmSSM.inputGroup.b obj.uvlmSSM.inputGroup.bDot])*Vloop;
                TAS=obj.tAS.d0([obj.tAS.outputGroup.b obj.uvlmSSM.inputGroup.bDot],obj.tAS.inputGroup.mode) +obj.tAS.dddV([obj.tAS.outputGroup.b obj.uvlmSSM.inputGroup.bDot],obj.tAS.inputGroup.mode)*Vloop;
                Kdca=D2-C2*A2^-1*B2;
                Kae=Kdcs*TSA*Kdca*TAS;
               plot(repmat(Vloop,size(Kae,1)),diag(Kae),'+');
            end
        end % divergenceAnalysis
    end
    
end

