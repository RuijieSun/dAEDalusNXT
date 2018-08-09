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
        % controllers
        ctr
        % actuators
        act
        % Settings (class:AeSSMSettings)
        settings@AeSSMSettings
        % linearized system
        linSSM
        % linearized system Velocity
        linV
        % V-g Data (Flutter Analysis)
        flutterTable
        % Gain Data (Diverg Analysis)
        divergenceData
        % Gust Data
        gustData
        % uvlm for ssm generation
        uvlm@class_UVLM_solver
        % strModel data needed (struct)
        strModel
        % state for initialization
        state@class_aero_state
        % reduction basis for aeroStates
        redBase
        % inverse of reduction basis for aeroStates
        redBaseInv
        % rank of the product of controllability and observability of the
        % aero system
        rankAero
        % reduced aero system with in and outputs as used for the reduction
        % in redAero
        redAeroSSM
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
			
			% recompute grid if necessary
            if ~obj.settings.deformedGrid
                aircraft_structure.nodal_deflections=zeros(size(aircraft_structure.Kff, 1), 1);
                aircraft_structure=aircraft_structure.f_postprocess();
                aircraft=aircraft.compute_deflected_grid(aircraft_structure.f_get_deflections);
            end
            if ~isempty(aircraft_structure.nodal_deflections)
                aircraft=aircraft.compute_deflected_grid(aircraft_structure.f_get_deflections);
            end
			
			% init aero model
			obj=obj.initUVLM(aircraft);
            
			% init flexible body model
			obj=obj.initSTRModel(aircraft_structure);
            
						
			% init coupling
			obj=obj.initCouplingMatrices(aircraft,aircraft_structure);
			
            % init Actuator Models (use default values for now)
            obj=obj.initActuatorModels();
            
            % init Sensor Models (position of aoa probe only)
            obj=obj.initSensorModels(aircraft);
            
			
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % reduction of TAS Aero TSA system with specified reduction
%             obj=obj.redAero(V1);
		end

        function obj=initUVLM(obj,aircraft) 
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % prepare aerodynamic model
            
            %init uvlm settings
            uvlmSettings=class_UVLM_computation_settings();
            uvlmSettings.debug=0;
            uvlmSettings.indComp=obj.settings.indComp;
            uvlmSettings.movie=0;
            uvlmSettings.wakelength_factor=obj.settings.wakeLengthFactor/sqrt(1-obj.state.Ma^2);
            uvlmSettings.wakeGrowthFactor=obj.settings.rW;
            uvlmSettings.nFixedWakePanels=obj.settings.nFixedWakePanels;
            uvlmSettings.wakePanelDensity=obj.settings.wakePanelDensity;
            uvlmSettings.addLengthFactorLastPanel=obj.settings.addLengthFactorLastPanel;
            % compute cs surf panels
            aircraft=aircraft.computeControlSurfacePanelIds;
            % compute info for cheb polynomials
            aircraft=aircraft.computeChebPoints;
            % init uvlm
            obj.uvlm=class_UVLM_solver(aircraft,obj.state,uvlmSettings);
            obj.uvlm=obj.uvlm.generateSSM(aircraft);
            
		end
		
        function obj=initSTRModel(obj,aircraft_structure)
		
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Prepare FD Model
            % velocity of the FLOW seen from the aeroBodySystem (positive
            % x Velocity means forward motion, positive z Velocity means
            % downward motion); no convertion needed in x and z when
            % transforming to flight dynamic body system; only beta and y
            % component change sign. What about the inertia?
            flowAbs=norm(obj.state.V_inf);
            flowDir=obj.state.V_inf./flowAbs;
            V1=flowAbs;
            V2=(flowAbs-10);
            if length(aircraft_structure.nodal_deflections)~=length(aircraft_structure.Mglob_lumped)
                aircraft_structure.nodal_deflections=[];
            end
            cg=aircraft_structure.f_compute_CG';
            cg(2)=0;
            m=aircraft_structure.f_compute_totalMass;
            Inertia=aircraft_structure.f_compute_moment_of_inertia(cg);
            fdSSM1=createFDSSM(m,V1*flowDir.*[1 -1 1],[0 obj.state.alpha -obj.state.beta]*pi/180,[0 0 0]', Inertia, 9.81);
            fdSSM2=createFDSSM(m,V2*flowDir.*[1 -1 1],[0 obj.state.alpha -obj.state.beta]*pi/180,[0 0 0]', Inertia, 9.81);
            fd_dadV=(fdSSM1.a-fdSSM2.a)./(V1-V2);
            fd_dbdV=(fdSSM1.b-fdSSM2.b)./(V1-V2);
            fd_dcdV=(fdSSM1.c-fdSSM2.c)./(V1-V2);
            fd_dddV=(fdSSM1.d-fdSSM2.d)./(V1-V2);
            fd_a0=fdSSM1.a-fd_dadV*V1;
            fd_b0=fdSSM1.b-fd_dbdV*V1;
            fd_c0=fdSSM1.c-fd_dcdV*V1;
            fd_d0=fdSSM1.d-fd_dddV*V1;
            fd_inputName=fdSSM1.InputName;
            fd_outputName=fdSSM1.OutputName;
            fd_stateName=fdSSM1.StateName;
            fd_outputGroup=fdSSM1.OutputGroup;
            fd_inputGroup=fdSSM1.InputGroup;
            % Transformation of L M N from reference point used in
            % aerodynamics to the CG
            D=eye(6);
            D(4:6,1:3)=skewsym(obj.uvlm.reference.p_ref-cg);
            fd_b0=fd_b0*D;
            fd_dbdV=fd_dbdV*D;
            fd_d0=fd_d0*D;
            fd_dddV=fd_dddV*D;
            % Transformation of vB/wB vBDot/wBDot around the CG (output of fd ssm)
            % to vB/wB vBDot/wBDot around the reference point (input to aero)
            D2=blkdiag(D',D');
            fd_c0([fd_outputGroup.vB fd_outputGroup.wB fd_outputGroup.vBDot fd_outputGroup.wBDot], :)=D2*fd_c0([fd_outputGroup.vB fd_outputGroup.wB fd_outputGroup.vBDot fd_outputGroup.wBDot], :);
            fd_d0([fd_outputGroup.vB fd_outputGroup.wB fd_outputGroup.vBDot fd_outputGroup.wBDot], :)=D2*fd_d0([fd_outputGroup.vB fd_outputGroup.wB fd_outputGroup.vBDot fd_outputGroup.wBDot], :);
            fd_dcdV([fd_outputGroup.vB fd_outputGroup.wB fd_outputGroup.vBDot fd_outputGroup.wBDot], :)=D2*fd_dcdV([fd_outputGroup.vB fd_outputGroup.wB fd_outputGroup.vBDot fd_outputGroup.wBDot], :);
            fd_dddV([fd_outputGroup.vB fd_outputGroup.wB fd_outputGroup.vBDot fd_outputGroup.wBDot], :)=D2*fd_dddV([fd_outputGroup.vB fd_outputGroup.wB fd_outputGroup.vBDot fd_outputGroup.wBDot], :);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % prepare structural model
            
            obj.strModel.modalBase=aircraft_structure.modeshapes(:,obj.settings.modes);
            obj.strModel.M=diag(diag(obj.strModel.modalBase'*aircraft_structure.Mff*obj.strModel.modalBase));
            obj.strModel.K=diag(diag(obj.strModel.modalBase'*aircraft_structure.Kff*obj.strModel.modalBase));
            if obj.settings.strRayleighDamping~=0
                obj.strModel.C=diag((obj.settings.strRayleighDamping*2*pi*aircraft_structure.modefrequencies(obj.settings.modes)));
            elseif obj.settings.strPropDamping~=0
                obj.strModel.C=diag(obj.settings.strPropDamping./(2*pi*aircraft_structure.modefrequencies(obj.settings.modes))).*obj.strModel.K;
            else
                obj.strModel.C=diag(0*(2*pi*aircraft_structure.modefrequencies(obj.settings.modes)));
            end
            nModes=length(obj.settings.modes);

            % init structural state space
            str_a=[ -obj.strModel.M^-1*obj.strModel.C, -obj.strModel.M^-1*obj.strModel.K; eye(nModes) zeros(nModes)];
            str_b=[obj.strModel.M^-1; zeros(nModes)];
            Cmodes=[str_a; zeros(nModes) eye(nModes)];
            Dmodes=[str_b; zeros(nModes)];
            defSave=aircraft_structure.nodal_deflections;
            if ~isempty(obj.settings.strLoadStations)
                %make load out matrix
                AllNodeLoadings=[];
                for iModForce=1:length(obj.settings.modes)

                    aircraft_structure.nodal_deflections=aircraft_structure.modeshapes(:,obj.settings.modes(iModForce));

                    aircraft_structure=aircraft_structure.f_postprocess();

                    %gather all nodal loadings
                    node_loadings=[];
                    for iBeam=1:length(aircraft_structure.beam)
                        node_loadings=[node_loadings; aircraft_structure.beam(iBeam).node_loadings_loc];
                    end
                    AllNodeLoadings=[AllNodeLoadings node_loadings];
                end
                loadIdx=sort([obj.settings.strLoadStations*6-5 obj.settings.strLoadStations*6-4 obj.settings.strLoadStations*6-3 obj.settings.strLoadStations*6-2 obj.settings.strLoadStations*6-1 obj.settings.strLoadStations*6]);
                Cloads=[zeros(length(loadIdx),length(obj.settings.modes)) AllNodeLoadings(loadIdx,:)];
                Dloads=zeros(length(loadIdx),length(obj.settings.modes));
            else
                Cloads=[];
                Dloads=[];
            end
            aircraft_structure.nodal_deflections=defSave;
            poiBase=aircraft_structure.modeshapes(obj.settings.poiDof,obj.settings.modes);
            poiBase=blkdiag(poiBase,poiBase,poiBase);
            Cpoi=poiBase*Cmodes;
            Dpoi=poiBase*Dmodes;
            
            str_c=[Cmodes; Cpoi; Cloads;];
            str_d=[Dmodes; Dpoi; Dloads;];
            nPoi=length(obj.settings.poiDof);
            % specification of input/output and statenames as well as in
            % and outputGroups
            str_stateName=[  cellstr([repmat('modeDot_',nModes,1) num2str([1:nModes]','%04d')]);...
                                    cellstr([repmat('mode_',nModes,1) num2str([1:nModes]','%04d')])];
            str_inputName=cellstr([repmat('modalForce_',nModes,1) num2str([1:nModes]','%04d')]);
            str_outputName=[  cellstr([repmat('modeDdot_',nModes,1) num2str([1:nModes]','%04d')]);...
                                     cellstr([repmat('modeDot_',nModes,1) num2str([1:nModes]','%04d')]);...
                                     cellstr([repmat('mode_',nModes,1) num2str([1:nModes]','%04d')]);];
            if nPoi>0
                str_outputName=[  str_outputName;...
                                         cellstr([repmat('poiDdot_',nPoi,1) num2str([1:nPoi]','%04d')]);...
                                         cellstr([repmat('poiDot_',nPoi,1) num2str([1:nPoi]','%04d')]);...
                                         cellstr([repmat('poi_',nPoi,1) num2str([1:nPoi]','%04d')]);];
            end
            for iLoadStation=1:length(obj.settings.strLoadStations)
                str_outputName=[str_outputName;
                                        ['load_' num2str(obj.settings.strLoadStations(iLoadStation)) '_T1'];
                                        ['load_' num2str(obj.settings.strLoadStations(iLoadStation)) '_T2'];
                                        ['load_' num2str(obj.settings.strLoadStations(iLoadStation)) '_T3'];
                                        ['load_' num2str(obj.settings.strLoadStations(iLoadStation)) '_R1'];
                                        ['load_' num2str(obj.settings.strLoadStations(iLoadStation)) '_R2'];
                                        ['load_' num2str(obj.settings.strLoadStations(iLoadStation)) '_R3'];];                                        
            end
            
           
            %merge fd + str
            if obj.settings.restrained
                nRBMDof=0;        
                obj.strSSM.a0=str_a;
                obj.strSSM.b0=str_b;
                obj.strSSM.c0=str_c;
                obj.strSSM.d0=str_d;
                obj.strSSM.dadV=0*str_a;
                obj.strSSM.dbdV=0*str_b;
                obj.strSSM.dcdV=0*str_c;
                obj.strSSM.dddV=0*str_d;

                obj.strSSM.stateName=[str_stateName];
                obj.strSSM.inputName=[str_inputName];
                obj.strSSM.outputName=[str_outputName];
            else
                nRBMDof=6;
                
                obj.strSSM.a0=blkdiag(fd_a0,str_a);
                obj.strSSM.b0=blkdiag(fd_b0,str_b);
                obj.strSSM.c0=blkdiag(fd_c0,str_c);
                obj.strSSM.d0=blkdiag(fd_d0,str_d);
                obj.strSSM.dadV=blkdiag(fd_dadV,0*str_a);
                obj.strSSM.dbdV=blkdiag(fd_dbdV,0*str_b);
                obj.strSSM.dcdV=blkdiag(fd_dcdV,0*str_c);
                obj.strSSM.dddV=blkdiag(fd_dddV,0*str_d);


                obj.strSSM.stateName=[fd_stateName;str_stateName];
                obj.strSSM.inputName=[fd_inputName;str_inputName];
                obj.strSSM.outputName=[fd_outputName;str_outputName];
                obj.strSSM.inputGroup=fd_inputGroup;
                obj.strSSM.outputGroup=fd_outputGroup;
            end
            obj.strSSM.inputGroup.modalForce=nRBMDof+(1:nModes);
            obj.strSSM.outputGroup.modeDdot=nRBMDof*4+(1:nModes);
            obj.strSSM.outputGroup.modeDot=nRBMDof*4+(nModes+1:2*nModes);
            obj.strSSM.outputGroup.mode=nRBMDof*4+(2*nModes+1:3*nModes);
            obj.strSSM.outputGroup.poiDdot=nRBMDof*4+(3*nModes+1:3*nModes+nPoi);
            obj.strSSM.outputGroup.poiDot=nRBMDof*4+(3*nModes+nPoi+1:3*nModes+2*nPoi);
            obj.strSSM.outputGroup.poi=nRBMDof*4+(3*nModes+2*nPoi+1:3*nModes+3*nPoi);
            obj.strSSM.outputGroup.loads=nRBMDof*4+(3*nModes+3*nPoi+1:3*nModes+3*nPoi+6*length(obj.settings.strLoadStations)); 
        end
        
        function obj=initCouplingMatrices(obj,aircraft,aircraft_structure)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % prepare coupling matrices TAS and TSA
            nB=size(obj.uvlm.panels,2);
            nModes=length(obj.settings.modes);
            flowAbs=norm(obj.state.V_inf);
            flowDir=obj.state.V_inf./flowAbs;
            if obj.settings.deformedGrid
                initDef=aircraft_structure.nodal_deflections; %<- CAREFULL!
            else
                if ~isempty(aircraft_structure.nodal_deflections)
                    initDef=aircraft_structure.nodal_deflections*0;
                else
                    initDef=zeros(size(aircraft_structure.modeshapes(:,1)));
                end
            end
            uvlmUnDef=obj.uvlm;
            uvlmUnDef.Ma_corr=1;
            
            aircraft_structure.nodal_deflections=initDef;
            aircraft_structure=aircraft_structure.f_postprocess();
            aircraftDef=aircraft.compute_deflected_grid(aircraft_structure.f_get_deflections);
            uvlmUnDef=uvlmUnDef.set_grid(aircraftDef.grid_deflected,aircraftDef.panels);
            
            uvlmDef=obj.uvlm;
            uvlmDef.Ma_corr=1;
            incr=aircraft.reference.b_ref/100;
            % dBdMode (is linear function of V) describes e.g. influence of mode on b and modeDot on bDot
            % because of the linear dependency on V, only ddBBdModedV is determined and
            % later multiplied by the velocity (see TAS assembly)
            ddBdModedV=zeros(nB,nModes);
            % dBdModeDot - describes the influence of modeDot on b and modeDotDot on bDot
            dBdModeDot=zeros(nB,nModes);
            dCollocdModeXYZ=zeros(nB*3,nModes);
            dFvapdModeXYZ=zeros(nB*3,nModes);
            dSegmentVecdMode=zeros(nB*3,nModes);
            collocNvec=uvlmUnDef.colloc_nvec;
            colloc=uvlmUnDef.colloc;
            fvap=uvlmUnDef.fvap;
            segmentVec=((0.75*uvlmUnDef.grid(:,uvlmUnDef.panels(2,:))+0.25*uvlmUnDef.grid(:,uvlmUnDef.panels(3,:)))-(0.75*uvlmUnDef.grid(:,uvlmUnDef.panels(1,:))+0.25*uvlmUnDef.grid(:,uvlmUnDef.panels(4,:))))';
            segmentVec=reshape(segmentVec',3*nB,1);

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
                segmentVecMode=((0.75*uvlmDef.grid(:,uvlmDef.panels(2,:))+0.25*uvlmDef.grid(:,uvlmDef.panels(3,:)))-(0.75*uvlmDef.grid(:,uvlmDef.panels(1,:))+0.25*uvlmDef.grid(:,uvlmDef.panels(4,:))))';
                segmentVecMode=reshape(segmentVecMode',3*nB,1);
                bMode=flowDir*collocNvecMode;
                ddBdModedV(:,iMode)=(bMode-bZero)/incr;
                dCollocdMode=(collocMode-colloc)/incr;
                dFvapdMode=(fvapMode-fvap)/incr;
                dSegmentVecdMode(:,iMode)=(segmentVecMode-segmentVec)/incr;
                dBdModeDot(:,iMode)=sum((dCollocdMode.*collocNvec),1);
                dCollocdModeXYZ(:,iMode)=reshape(dCollocdMode,nB*3,1);
                dFvapdModeXYZ(:,iMode)=reshape(dFvapdMode,nB*3,1);
            end
        %TAS data
            obj.tAS.d0=[  zeros(nB,nModes)   -dBdModeDot      zeros(nB,nModes);...
                        -dBdModeDot          zeros(nB,nModes) zeros(nB,nModes);...
                        zeros(3*nB,nModes)    zeros(3*nB,nModes)  dSegmentVecdMode;...
                        zeros(nB,nModes)   -dFvapdModeXYZ(1:3:end,:) zeros(nB,nModes) ;...
                        zeros(nB,nModes)   -dFvapdModeXYZ(2:3:end,:) zeros(nB,nModes) ;...
                        zeros(nB,nModes)   -dFvapdModeXYZ(3:3:end,:) zeros(nB,nModes) ;];
                    
            obj.tAS.dddV=[  zeros(nB,nModes)   zeros(nB,nModes)      ddBdModedV;...
                            zeros(nB,nModes)       ddBdModedV zeros(nB,nModes);...
                        zeros(3*nB,nModes)    zeros(3*nB,nModes)   zeros(3*nB,nModes);...
                        zeros(3*nB,nModes)    zeros(3*nB,nModes)   zeros(3*nB,nModes);];
                    
            obj.tAS.inputName=[	cellstr([repmat('modeDdot_',nModes,1) num2str([1:nModes]','%04d')]);...
                                cellstr([repmat('modeDot_',nModes,1) num2str([1:nModes]','%04d')]);...
                                cellstr([repmat('mode_',nModes,1) num2str([1:nModes]','%04d')])];
            obj.tAS.inputGroup.modeDdot=1:nModes;
            obj.tAS.inputGroup.modeDot=nModes+1:2*nModes;
            obj.tAS.inputGroup.mode=2*nModes+1:3*nModes;
            obj.tAS.outputName=[obj.uvlm.ssm.LPVKernel.inputNames([obj.uvlm.ssm.LPVKernel.InputGroup.b obj.uvlm.ssm.LPVKernel.InputGroup.bDot]);...
                                obj.uvlm.ssm.linOutT.InputName([obj.uvlm.ssm.linOutT.InputGroup.s obj.uvlm.ssm.linOutT.InputGroup.vSInf_x  obj.uvlm.ssm.linOutT.InputGroup.vSInf_y  obj.uvlm.ssm.linOutT.InputGroup.vSInf_z])];
            obj.tAS.outputGroup.b=1:nB;
            obj.tAS.outputGroup.bDot=nB+1:2*nB;        
            obj.tAS.outputGroup.s=2*nB+1:5*nB;
            obj.tAS.outputGroup.vSInf_x=5*nB+1:6*nB;   
            obj.tAS.outputGroup.vSInf_y=6*nB+1:7*nB;
            obj.tAS.outputGroup.vSInf_z=7*nB+1:8*nB;  

            %TSA data
            obj.tSA.d=dFvapdModeXYZ';
            obj.tSA.inputName=obj.uvlm.ssm.linOutT.OutputName(obj.uvlm.ssm.linOutT.OutputGroup.Fp);
            obj.tSA.inputGroup.Fp=1:length(obj.uvlm.ssm.linOutT.OutputGroup.Fp);
            obj.tSA.outputName=cellstr([repmat('modalForce_',nModes,1) num2str([1:nModes]','%04d')]);
            obj.tSA.outputGroup.modalForce=1:nModes;
        end
        
        function obj=initActuatorModels(obj)
            %for each control surface pair (ailerons, spoilers, elevator,
            %movCS) 2 control paths are defined (sym and asym def)
            %for every single cs (rudder) one act path is defined
            % this function connects partial control surfaces if needed
            % fixed parameter
            w0 = 10;
            DampCoeff = 0.7;
            

            % find different cs names and store ids of cs in cell
            for iUvlmCS=1:length(obj.uvlm.csNames)
                currentSplit=strsplit(obj.uvlm.csNames{iUvlmCS},'_');
                allCsNames{iUvlmCS}=currentSplit{1};
                
            end
            [csNames, iA, iC]=unique(allCsNames);
            
            %create dummyActuator (as long as they have all the same fixed
            %parameters)
            a = [0, 1; -w0^2, -2*DampCoeff*w0];
            b = [0; w0^2];
            c = [1, 0; 0, 1; -w0^2, -2*DampCoeff*w0];
            d = [0; 0; w0^2];
            %%
            dummyAct=ss(a,b,c,d);
            allAct=ss([]);
            %
            for iCs=1:length(csNames)
                dummyAct.StateName={[csNames{iCs} '_ActState1']; [csNames{iCs} '_ActState2']};
                %is left and right?
                leftRightFlag=any(~(cellfun(@isempty,strfind(obj.uvlm.csNames(iC==iCs),'_l'))));
                outputStartId=length(allAct.outputName)/3+1;
                if leftRightFlag % create two actuators
                    symAct=dummyAct;
                    symAct.StateName={[csNames{iCs} '_ActState1sym']; [csNames{iCs} '_ActState2sym']};
                    symAct.InputName=[csNames{iCs} '_sym'];
                    %add left outputs to symAct
                    symAct=series(symAct,ss([eye(3);-eye(3)]));
                    symAct.OutputName={['dCs_' num2str([(outputStartId)]','%04d')]; ['dCsDot_' num2str([(outputStartId)]','%04d')]; ['dCsDotDot_' num2str([(outputStartId)]','%04d')]; ['dCs_' num2str([(outputStartId+1)]','%04d')]; ['dCsDot_' num2str([(outputStartId+1)]','%04d')]; ['dCsDotDot_' num2str([(outputStartId+1)]','%04d')]};
                    
                    aSymAct=dummyAct;
                    aSymAct.StateName={[csNames{iCs} '_ActState1asym']; [csNames{iCs} '_ActState2asym']};
                    aSymAct.InputName=[csNames{iCs} '_asym'];
                    %add left outputs to aSymAct
                    aSymAct=series(aSymAct,ss([eye(3);eye(3)]));
                    aSymAct.OutputName={['dCs_' num2str([(outputStartId)]','%04d')]; ['dCsDot_' num2str([(outputStartId)]','%04d')]; ['dCsDotDot_' num2str([(outputStartId)]','%04d')]; ['dCs_' num2str([(outputStartId+1)]','%04d')]; ['dCsDot_' num2str([(outputStartId+1)]','%04d')]; ['dCsDotDot_' num2str([(outputStartId+1)]','%04d')]};
                    %add outputs for parts
                    if length(obj.uvlm.csNames(iC==iCs))==4
                        disp('not implemented yet')
                        break
                    elseif length(obj.uvlm.csNames(iC==iCs))==6
                        disp('not implemented yet')
                        break
                    end
                    %name in and outputs
                    %append to allAct
                    allAct=append(allAct,parallel(symAct,aSymAct,'name'));
                    
                else %
                    %only add ouputs
                    oneAct=dummyAct;
                    %name in and ouputs
                    oneAct.InputName=[csNames{iCs}];
                    oneAct.OutputName={['dCs_' num2str([(outputStartId)]','%04d')]; ['dCsDot_' num2str([(outputStartId)]','%04d')]; ['dCsDotDot_' num2str([(outputStartId)]','%04d')];};
                    %add outputs for parts
                    if length(obj.uvlm.csNames(iC==iCs))==2
                        disp('not implemented yet')
                        break
                    elseif length(obj.uvlm.csNames(iC==iCs))==3
                        disp('not implemented yet')
                        break
                    end
                    
                    %append to allAct
                    allAct=append(allAct,oneAct);
                end
            end
            
            
            % store rate limits for each cs
            
            allAct.OutputGroup.dCs=[1:3:size(allAct,1)];
            allAct.OutputGroup.dCsDot=[2:3:size(allAct,1)];
            allAct.OutputGroup.dCsDotDot=[3:3:size(allAct,1)];
            allAct.InputGroup.controlCommand=[1:size(allAct,2)];
            obj.act.SSM=allAct;
            obj.act.rateLimit=50;
        end
        
        function obj=initSensorModels(obj, aircraft)
            %calculate the quarterpoint of the first gust zone (aerodynamic
            %center!?)
            quarterPoint=sum(obj.uvlm.colloc(:,obj.uvlm.gustZones{1}).*repmat(obj.uvlm.area(obj.uvlm.gustZones{1}),3,1),2)./sum(obj.uvlm.area(obj.uvlm.gustZones{1}));
            %calculate the position of the nose in the same coordinate frame
            if ~isempty(aircraft.fuselages)
                 nose=aircraft.fuselages.fuselage_segments(1).pos';
                %relative pos:
                obj.settings.aoaProbePos=nose(1)-quarterPoint(1);
            else
                obj.settings.aoaProbePos=-30;
            end
        end
        
        function obj=redAero(obj,V1,varargin) %<- this should be in aero ssm class
            if nargin==4
                obj.redBase=varargin{1};
                obj.redBaseInv=varargin{2};
            else
                if strcmp(obj.settings.mor,'bPod')
                    %bPodMod takes real modeshapes as inputs
                    inputsOfAeroUsedForReduction=[obj.uvlm.ssm.ssm.inputGroup.b,...
                                            obj.uvlm.ssm.ssm.inputGroup.bDot,...
                                            obj.uvlm.ssm.ssm.inputGroup.dCs,...
                                            obj.uvlm.ssm.ssm.inputGroup.dCsDot,...
                                            obj.uvlm.ssm.ssm.inputGroup.dCsDotDot,...
                                            obj.uvlm.ssm.ssm.inputGroup.s,...
                                            obj.uvlm.ssm.ssm.inputGroup.vSInf_x,...
                                            obj.uvlm.ssm.ssm.inputGroup.vSInf_y,...
                                            obj.uvlm.ssm.ssm.inputGroup.vSInf_z];
%                                             obj.uvlm.ssm.ssm.inputGroup.vG,...
%                                             obj.uvlm.ssm.ssm.inputGroup.vGDot,...
%                                             obj.uvlm.ssm.ssm.inputGroup.vB,...
%                                             obj.uvlm.ssm.ssm.inputGroup.vBDot,...
%                                             obj.uvlm.ssm.ssm.inputGroup.wB,...
%                                             obj.uvlm.ssm.ssm.inputGroup.wBDot,...
                    outputsOfAeroUsedForReduction=obj.uvlm.ssm.ssm.outputGroup.Fp;
                    aerSSM=ss(  full(obj.uvlm.ssm.ssm.a0+obj.uvlm.ssm.ssm.dadV*V1),...
                                full(obj.uvlm.ssm.ssm.b0+obj.uvlm.ssm.ssm.dbdV*V1+obj.uvlm.ssm.ssm.dbddV*V1^2),...
                                full(obj.uvlm.ssm.ssm.c0+obj.uvlm.ssm.ssm.dcdV*V1+obj.uvlm.ssm.ssm.dcddV*V1^2+obj.uvlm.ssm.ssm.dcdddV*V1^3),...
                                full(obj.uvlm.ssm.ssm.d0+obj.uvlm.ssm.ssm.dddV*V1+obj.uvlm.ssm.ssm.ddddV*V1^2+obj.uvlm.ssm.ssm.dddddV*V1^3+obj.uvlm.ssm.ssm.ddddddV*V1^4));
                    aerSSM.inputName=obj.uvlm.ssm.ssm.inputName;
                    aerSSM.InputGroup=obj.uvlm.ssm.ssm.inputGroup;
                    aerSSM.outputName=obj.uvlm.ssm.ssm.outputName;
                    aerSSM.OutputGroup=obj.uvlm.ssm.ssm.outputGroup;
                    aerSSM.StateName=obj.uvlm.ssm.ssm.stateName;
                    %remove inputs vP and vPDot
                    aerSSM=aerSSM(outputsOfAeroUsedForReduction,inputsOfAeroUsedForReduction);
                    tASSSM=ss(obj.tAS.d0+obj.tAS.dddV*V1);
                    tASSSM.inputName=obj.tAS.inputName;
                    tASSSM.InputGroup=obj.tAS.inputGroup;
                    tASSSM.outputName=obj.tAS.outputName;
                    tASSSM.OutputGroup=obj.tAS.outputGroup;
                    aerSSM2=seriesPreserve(tASSSM,aerSSM,{'b','bDot','s','vSInf_x','vSInf_y','vSInf_z'});
                    tSASSM=ss(obj.tSA.d);
                    tSASSM.inputName=obj.tSA.inputName;
                    tSASSM.InputGroup=obj.tSA.inputGroup;
                    tSASSM.outputName=obj.tSA.outputName;
                    tSASSM.OutputGroup=obj.tSA.outputGroup;
                    aerSSM3=seriesPreserve(aerSSM2,tSASSM,{'Fp'});                   
                    %% reduce    
                    
                    tic;
                    [obj.redBase, obj.redBaseInv,obj.rankAero, singVal]=bpod(aerSSM3,obj.settings.morOrder, obj.settings.bpodT, obj.settings.errPODProj);
                    toc;
                elseif or(strcmp(obj.settings.mor,'bPodCheb'),strcmp(obj.settings.mor,'rPodCheb'))
                    %% identify separate lifting surfaces
                    
                    %% generate cheb modes
                    nAll=size(obj.uvlm.panels,2);
                    nWings=size(obj.uvlm.dataForCheb.spanWisePoints,2);
                    
                    modes=cell(nWings,1);
                    for iWing=1:nWings
                        orderMaxX=obj.settings.morChebOrder(2);
                        orderMaxY=ceil((obj.settings.morChebOrder(1)+1)/obj.uvlm.reference.b_ref*obj.uvlm.dataForCheb.spanWiseLength{iWing})-1;
                        for iX=0:orderMaxX
                            x=obj.uvlm.dataForCheb.chordWisePoints{iWing};
                            Tx=((x+sqrt(x.^2-1)).^iX+(x-sqrt(x.^2-1)).^iX)/2;
                            for iY=0:orderMaxY
                                x=obj.uvlm.dataForCheb.spanWisePoints{iWing};
                                Ty=((x+sqrt(x.^2-1)).^iY+(x-sqrt(x.^2-1)).^iY)/2;
                                modes{iWing}(:,(iX)*(orderMaxY+1)+iY+1)=(Tx.*Ty)';
                            end
                        end
                    end
                    modes=blkdiag(modes{:});
                    %% gen sys 2 reduce
                    obj.uvlm.ssm=obj.uvlm.ssm.linearize(V1);
                    uvlmSSM=obj.uvlm.ssm.linSSM;
                    inTrafo=ss(modes);
                    inTrafo.InputName=[cellstr([repmat('mode_',size(modes,2),1) num2str([1:size(modes,2)]','%04d')]);];
                    inTrafo.InputGroup.mode=1:size(modes,2);

                    inTrafo.OutputName=uvlmSSM.InputName([uvlmSSM.InputGroup.b]) ;
                    inTrafo.OutputGroup.b=1:length(uvlmSSM.InputGroup.b);

                    sys2=seriesPreserve(inTrafo,uvlmSSM,{'b'});
                    A=sys2.a;
                    B=sys2.b(:,[sys2.InputGroup.mode]); 
                    C=blkdiag(modes,modes)'*(obj.uvlm.ssm.LPVKernel.c0([obj.uvlm.ssm.LPVKernel.OutputGroup.gbeff obj.uvlm.ssm.LPVKernel.OutputGroup.gbDoteff],:)+obj.uvlm.ssm.LPVKernel.dcdV([obj.uvlm.ssm.LPVKernel.OutputGroup.gbeff obj.uvlm.ssm.LPVKernel.OutputGroup.gbDoteff],:)*V1);
                    D=zeros(size(C,1),size(B,2));
                    inDelay=sys2.inputdelay([sys2.InputGroup.mode]);
                    sys4red=ss(A,B,C,D);
                    sys4red.InputDelay=inDelay;                    
                    %% reduce      
                    tic;
                    if strcmp(obj.settings.mor,'bPodCheb')
                        [obj.redBase, obj.redBaseInv,obj.rankAero, singVal]=bpod(sys4red,obj.settings.morOrder, obj.settings.bpodT, obj.settings.errPODProj);
                    elseif strcmp(obj.settings.mor,'rPodCheb')
                        [obj.redBase, obj.redBaseInv,obj.rankAero, singVal]=rpod(sys4red,obj.settings.morOrder, obj.settings.bpodT, obj.settings.errPODProj);
                    end
                    toc;
                elseif or(or(strcmp(obj.settings.mor,'bPodRBF'),strcmp(obj.settings.mor,'rPodRBF')),startsWith(obj.settings.mor,'bPodRBF'))
                    %% identify factor for rbfYRad and rbfYRad
                    if length(obj.settings.mor)==8
                        radFactS=str2double(obj.settings.mor(8));
                        radFactC=str2double(obj.settings.mor(8));
                    elseif length(obj.settings.mor)==7
                        radFactS=2; %default setting
                        radFactC=2; %default setting
                    elseif length(obj.settings.mor)==9
                        radFactS=str2double(obj.settings.mor(8:9));
                        radFactC=str2double(obj.settings.mor(8:9));
                    elseif length(obj.settings.mor)==11
                        radFactS=str2double(obj.settings.mor(8:9));
                        radFactC=str2double(obj.settings.mor(10:11));
                    end
                    %% generate rbf modes
                    nAll=size(obj.uvlm.panels,2);
                    nWings=size(obj.uvlm.dataForCheb.spanWisePoints,2);
                    
                    modes=cell(nWings,1);
                    for iWing=1:nWings
                        orderMaxX=obj.settings.morChebOrder(2)+1;
                        orderMaxY=ceil((obj.settings.morChebOrder(1)+1)/obj.uvlm.reference.b_ref*obj.uvlm.dataForCheb.spanWiseLength{iWing});
                       
                        if orderMaxX==1
                            rbfXPos=0;
                            rbfXRad=0.5*radFactC;
                        else
                            rbfXPos=linspace(-1,1,orderMaxX);
                            rbfXRad=radFactC/orderMaxX;
                        end
                        if orderMaxY==1
                            rbfYPos=0;
                            rbfYRad=radFactS*0.5;
                        else
                            rbfYPos=linspace(-1,1,orderMaxY);
                            rbfYRad=radFactS/orderMaxY;
                        end
                        for iX=1:orderMaxX
                            etaX=abs(obj.uvlm.dataForCheb.chordWisePoints{iWing}-rbfXPos(iX))./rbfXRad; 
                            etaX(etaX>1)=1;
                            rbfX=(1-etaX).^4.*(4.*etaX+1);
                            for iY=1:orderMaxY
                                etaY=abs(obj.uvlm.dataForCheb.spanWisePoints{iWing}-rbfYPos(iY))./rbfYRad;
                                etaY(etaY>1)=1;
                                rbfY=(1-etaY).^4.*(4.*etaY+1);
                                 modes{iWing}(:,(iX-1)*(orderMaxY)+iY)=(rbfX.*rbfY)';
                            end
                        end
                        
                    end
                    modes=blkdiag(modes{:});
                    %% gen sys 2 reduce
                    obj.uvlm.ssm=obj.uvlm.ssm.linearize(V1);
                    uvlmSSM=obj.uvlm.ssm.linSSM;
                    inTrafo=ss(modes);
                    inTrafo.InputName=[cellstr([repmat('mode_',size(modes,2),1) num2str([1:size(modes,2)]','%04d')]);];
                    inTrafo.InputGroup.mode=1:size(modes,2);

                    inTrafo.OutputName=uvlmSSM.InputName([uvlmSSM.InputGroup.b]) ;
                    inTrafo.OutputGroup.b=1:length(uvlmSSM.InputGroup.b);

                    sys2=seriesPreserve(inTrafo,uvlmSSM,{'b'});
                    A=sys2.a;
                    B=sys2.b(:,[sys2.InputGroup.mode]); 
                    C=blkdiag(modes,modes)'*(obj.uvlm.ssm.LPVKernel.c0([obj.uvlm.ssm.LPVKernel.OutputGroup.gbeff obj.uvlm.ssm.LPVKernel.OutputGroup.gbDoteff],:)+obj.uvlm.ssm.LPVKernel.dcdV([obj.uvlm.ssm.LPVKernel.OutputGroup.gbeff obj.uvlm.ssm.LPVKernel.OutputGroup.gbDoteff],:)*V1);
                    D=zeros(size(C,1),size(B,2));
                    inDelay=sys2.inputdelay([sys2.InputGroup.mode]);
                    sys4red=ss(A,B,C,D);
                    sys4red.InputDelay=inDelay;                    
                    %% reduce      
                    tic;
                    if startsWith(obj.settings.mor,'bPodRBF')
                        [obj.redBase, obj.redBaseInv,obj.rankAero, singVal]=bpod(sys4red,obj.settings.morOrder, obj.settings.bpodT, obj.settings.errPODProj);
                    elseif strcmp(obj.settings.mor,'rPodRBF')
                        [obj.redBase, obj.redBaseInv,obj.rankAero, singVal]=rpod(sys4red,obj.settings.morOrder, obj.settings.bpodT, obj.settings.errPODProj);
                    end
                    toc;
                        
                elseif or(strcmp(obj.settings.mor,'rPodZones'),strcmp(obj.settings.mor,'bPodZones'))
                    nAll=size(obj.uvlm.panels,2);
                    nWings=size(obj.uvlm.dataForCheb.spanWisePoints,2);
                    
                    modes=cell(nWings,1);
                    for iWing=1:nWings
                        nSectionsX=obj.settings.morChebOrder(2)+1;
                        nSectionsY=ceil((obj.settings.morChebOrder(1)+1)/obj.uvlm.reference.b_ref*obj.uvlm.dataForCheb.spanWiseLength{iWing});
                        dividersX=linspace(-1,1,nSectionsX+1);
                        dividersY=linspace(-1,1,nSectionsY+1);
                        for iX=1:nSectionsX
                            %get ids of this x section
                            if iX==nSectionsX
                                idxX=and(obj.uvlm.dataForCheb.chordWisePoints{iWing}<=dividersX(iX+1),obj.uvlm.dataForCheb.chordWisePoints{iWing}>=dividersX(iX));                            
                            else
                                idxX=and(obj.uvlm.dataForCheb.chordWisePoints{iWing}<dividersX(iX+1),obj.uvlm.dataForCheb.chordWisePoints{iWing}>=dividersX(iX));
                            end
                            
                            
                            for iY=1:nSectionsY
                                %get ids of this y section
                                if iY==nSectionsY
                                    idxY=and(obj.uvlm.dataForCheb.spanWisePoints{iWing}<=dividersY(iY+1),obj.uvlm.dataForCheb.spanWisePoints{iWing}>=dividersY(iY));                            
                                else
                                    idxY=and(obj.uvlm.dataForCheb.spanWisePoints{iWing}<dividersY(iY+1),obj.uvlm.dataForCheb.spanWisePoints{iWing}>=dividersY(iY));
                                end
                                % current mode is the panels which are in
                                % this x and y section
                                modes{iWing}(:,(iX-1)*nSectionsY+iY)=and(idxX,idxY)';
                            end
                        end
                    end
                    modes=blkdiag(modes{:});
                    %% gen sys 2 reduce
                    obj.uvlm.ssm=obj.uvlm.ssm.linearize(V1);
                    uvlmSSM=obj.uvlm.ssm.linSSM;
                    inTrafo=ss(modes);
                    inTrafo.InputName=[cellstr([repmat('mode_',size(modes,2),1) num2str([1:size(modes,2)]','%04d')]);];
                    inTrafo.InputGroup.mode=1:size(modes,2);

                    inTrafo.OutputName=uvlmSSM.InputName([uvlmSSM.InputGroup.b]) ;
                    inTrafo.OutputGroup.b=1:length(uvlmSSM.InputGroup.b);

                    sys2=seriesPreserve(inTrafo,uvlmSSM,{'b'});
                    A=sys2.a;
                    B=sys2.b(:,[sys2.InputGroup.mode]); 
                    C=blkdiag(modes,modes)'*(obj.uvlm.ssm.LPVKernel.c0([obj.uvlm.ssm.LPVKernel.OutputGroup.gbeff obj.uvlm.ssm.LPVKernel.OutputGroup.gbDoteff],:)+obj.uvlm.ssm.LPVKernel.dcdV([obj.uvlm.ssm.LPVKernel.OutputGroup.gbeff obj.uvlm.ssm.LPVKernel.OutputGroup.gbDoteff],:)*V1);
                    D=zeros(size(C,1),size(B,2));
                    inDelay=sys2.inputdelay([sys2.InputGroup.mode]);
                    sys4red=ss(A,B,C,D);
                    sys4red.InputDelay=inDelay;
                    %% reduce
                    
                    tic;
                    if strcmp(obj.settings.mor,'bPodZones')
                        [obj.redBase, obj.redBaseInv,obj.rankAero, singVal]=bpod(sys4red,obj.settings.morOrder, obj.settings.bpodT, obj.settings.errPODProj);
                        
                    elseif strcmp(obj.settings.mor,'rPodZones')
                        [obj.redBase, obj.redBaseInv,obj.rankAero, singVal]=rpod(sys4red,obj.settings.morOrder, obj.settings.bpodT, obj.settings.errPODProj);
                    end
                    toc;
                elseif strcmp(obj.settings.mor,'bPodStrip')
                    %% generate strip modes
                    testsum=cumsum(obj.uvlm.is_te(end:-1:1));
                    testsum=testsum(end:-1:1);
                    stripModes=zeros(length(testsum),max(testsum));
                    for iStripMode=1:max(testsum)
                        stripModes(find(testsum==(max(testsum)+1)-iStripMode),iStripMode)=1;
                    end
                    %% gen sys 2 reduce
                    obj.uvlm.ssm=obj.uvlm.ssm.linearize(V1);
                    uvlmSSM=obj.uvlm.ssm.linSSM;
                    inTrafo=ss(stripModes);
                    inTrafo.InputName=[cellstr([repmat('mode_',size(stripModes,2),1) num2str([1:size(stripModes,2)]','%04d')]);];
                    inTrafo.InputGroup.mode=1:size(stripModes,2);

                    inTrafo.OutputName=uvlmSSM.InputName([uvlmSSM.InputGroup.b]) ;
                    inTrafo.OutputGroup.b=1:length(uvlmSSM.InputGroup.b);

                    sys2=seriesPreserve(inTrafo,uvlmSSM,{'b'});
                    A=sys2.a;
                    B=sys2.b(:,[sys2.InputGroup.mode]); 
                    C=blkdiag(stripModes,stripModes)'*(obj.uvlm.ssm.LPVKernel.c0([obj.uvlm.ssm.LPVKernel.OutputGroup.gbeff obj.uvlm.ssm.LPVKernel.OutputGroup.gbDoteff],:)+obj.uvlm.ssm.LPVKernel.dcdV([obj.uvlm.ssm.LPVKernel.OutputGroup.gbeff obj.uvlm.ssm.LPVKernel.OutputGroup.gbDoteff],:)*V1);
                    D=zeros(size(C,1),size(B,2));
                    inDelay=sys2.inputdelay([sys2.InputGroup.mode]);
                    sys4red=ss(A,B,C,D);
                    sys4red.InputDelay=inDelay;
                    %% reduce
                    
                    tic;
                        [obj.redBase, obj.redBaseInv,obj.rankAero, singVal]=bpod(sys4red,obj.settings.morOrder, obj.settings.bpodT, obj.settings.errPODProj);
                    toc;
                    
                elseif strcmp(obj.settings.mor,'bPodPan')
                    %% generate pan modes
                    panModes=eye(size(obj.uvlm.panels,2));
                    %% gen sys 2 reduce
                    obj.uvlm.ssm=obj.uvlm.ssm.linearize(V1);
                    uvlmSSM=obj.uvlm.ssm.linSSM;
                    inTrafo=ss(panModes);
                    inTrafo.InputName=[cellstr([repmat('mode_',size(panModes,2),1) num2str([1:size(panModes,2)]','%04d')]);];
                    inTrafo.InputGroup.mode=1:size(panModes,2);

                    inTrafo.OutputName=uvlmSSM.InputName([uvlmSSM.InputGroup.b]) ;
                    inTrafo.OutputGroup.b=1:length(uvlmSSM.InputGroup.b);

                    sys2=seriesPreserve(inTrafo,uvlmSSM,{'b'});
                    A=sys2.a;
                    B=sys2.b(:,[sys2.InputGroup.mode]); 
                    C=blkdiag(panModes,panModes)'*(obj.uvlm.ssm.LPVKernel.c0([obj.uvlm.ssm.LPVKernel.OutputGroup.gbeff obj.uvlm.ssm.LPVKernel.OutputGroup.gbDoteff],:)+obj.uvlm.ssm.LPVKernel.dcdV([obj.uvlm.ssm.LPVKernel.OutputGroup.gbeff obj.uvlm.ssm.LPVKernel.OutputGroup.gbDoteff],:)*V1);
                    D=zeros(size(C,1),size(B,2));
                    inDelay=sys2.inputdelay([sys2.InputGroup.mode]);
                    sys4red=ss(A,B,C,D);
                    sys4red.InputDelay=inDelay;
                    %% reduce
                    
                    tic;
                        [obj.redBase, obj.redBaseInv,obj.rankAero, singVal]=bpod(sys4red,obj.settings.morOrder, obj.settings.bpodT, obj.settings.errPODProj);
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
            obj.uvlm.ssm=obj.uvlm.ssm.reduce(obj.redBase, obj.redBaseInv);
        end
          
        function modes=genSynthModes(obj,V1,varargin) %<- this should be in aero ssm class
            if nargin==4
                obj.redBase=varargin{1};
                obj.redBaseInv=varargin{2};
            else
                if strcmp(obj.settings.mor,'bPod')
                  
                elseif or(strcmp(obj.settings.mor,'bPodCheb'),strcmp(obj.settings.mor,'rPodCheb'))
                    %% identify separate lifting surfaces
                    
                    %% generate cheb modes
                    nAll=size(obj.uvlm.panels,2);
                    nWings=size(obj.uvlm.dataForCheb.spanWisePoints,2);
                    
                    modes=cell(nWings,1);
                    for iWing=1:nWings
                        orderMaxX=obj.settings.morChebOrder(2);
                        orderMaxY=ceil((obj.settings.morChebOrder(1)+1)/obj.uvlm.reference.b_ref*obj.uvlm.dataForCheb.spanWiseLength{iWing})-1;
                        for iX=0:orderMaxX
                            x=obj.uvlm.dataForCheb.chordWisePoints{iWing};
                            Tx=((x+sqrt(x.^2-1)).^iX+(x-sqrt(x.^2-1)).^iX)/2;
                            for iY=0:orderMaxY
                                x=obj.uvlm.dataForCheb.spanWisePoints{iWing};
                                Ty=((x+sqrt(x.^2-1)).^iY+(x-sqrt(x.^2-1)).^iY)/2;
                                modes{iWing}(:,(iX)*(orderMaxY+1)+iY+1)=(Tx.*Ty)';
                            end
                        end
                    end
                    modes=blkdiag(modes{:});
                    %% gen sys 2 reduce                    
                    %% reduce   
                    modes=modes;
                elseif or(or(strcmp(obj.settings.mor,'bPodRBF'),strcmp(obj.settings.mor,'rPodRBF')),startsWith(obj.settings.mor,'bPodRBF'))
                    %% identify factor for rbfYRad and rbfYRad
                    if length(obj.settings.mor)==8
                        radFactS=str2double(obj.settings.mor(8));
                        radFactC=str2double(obj.settings.mor(8));
                    elseif length(obj.settings.mor)==7
                        radFactS=2; %default setting
                        radFactC=2; %default setting
                    elseif length(obj.settings.mor)==9
                        radFactS=str2double(obj.settings.mor(8:9));
                        radFactC=str2double(obj.settings.mor(8:9));
                    elseif length(obj.settings.mor)==11
                        radFactS=str2double(obj.settings.mor(8:9));
                        radFactC=str2double(obj.settings.mor(10:11));
                    end
                    %% generate rbf modes
                    nAll=size(obj.uvlm.panels,2);
                    nWings=size(obj.uvlm.dataForCheb.spanWisePoints,2);
                    
                    modes=cell(nWings,1);
                    for iWing=1:nWings
                        orderMaxX=obj.settings.morChebOrder(2)+1;
                        orderMaxY=ceil((obj.settings.morChebOrder(1)+1)/obj.uvlm.reference.b_ref*obj.uvlm.dataForCheb.spanWiseLength{iWing});
                       
                        if orderMaxX==1
                            rbfXPos=0;
                            rbfXRad=0.5*radFactC;
                        else
                            rbfXPos=linspace(-1,1,orderMaxX);
                            rbfXRad=radFactC/orderMaxX;
                        end
                        if orderMaxY==1
                            rbfYPos=0;
                            rbfYRad=radFactS*0.5;
                        else
                            rbfYPos=linspace(-1,1,orderMaxY);
                            rbfYRad=radFactS/orderMaxY;
                        end
                        for iX=1:orderMaxX
                            etaX=abs(obj.uvlm.dataForCheb.chordWisePoints{iWing}-rbfXPos(iX))./rbfXRad; 
                            etaX(etaX>1)=1;
                            rbfX=(1-etaX).^4.*(4.*etaX+1);
                            for iY=1:orderMaxY
                                etaY=abs(obj.uvlm.dataForCheb.spanWisePoints{iWing}-rbfYPos(iY))./rbfYRad;
                                etaY(etaY>1)=1;
                                rbfY=(1-etaY).^4.*(4.*etaY+1);
                                 modes{iWing}(:,(iX-1)*(orderMaxY)+iY)=(rbfX.*rbfY)';
                            end
                        end
                        
                    end
                    modes=blkdiag(modes{:});
                    %% gen sys 2 reduce                  
                    %% reduce    
                    modes=modes;
                        
                elseif or(strcmp(obj.settings.mor,'rPodZones'),strcmp(obj.settings.mor,'bPodZones'))
                    nAll=size(obj.uvlm.panels,2);
                    nWings=size(obj.uvlm.dataForCheb.spanWisePoints,2);
                    
                    modes=cell(nWings,1);
                    for iWing=1:nWings
                        nSectionsX=obj.settings.morChebOrder(2)+1;
                        nSectionsY=ceil((obj.settings.morChebOrder(1)+1)/obj.uvlm.reference.b_ref*obj.uvlm.dataForCheb.spanWiseLength{iWing});
                        dividersX=linspace(-1,1,nSectionsX+1);
                        dividersY=linspace(-1,1,nSectionsY+1);
                        for iX=1:nSectionsX
                            %get ids of this x section
                            if iX==nSectionsX
                                idxX=and(obj.uvlm.dataForCheb.chordWisePoints{iWing}<=dividersX(iX+1),obj.uvlm.dataForCheb.chordWisePoints{iWing}>=dividersX(iX));                            
                            else
                                idxX=and(obj.uvlm.dataForCheb.chordWisePoints{iWing}<dividersX(iX+1),obj.uvlm.dataForCheb.chordWisePoints{iWing}>=dividersX(iX));
                            end
                            
                            
                            for iY=1:nSectionsY
                                %get ids of this y section
                                if iY==nSectionsY
                                    idxY=and(obj.uvlm.dataForCheb.spanWisePoints{iWing}<=dividersY(iY+1),obj.uvlm.dataForCheb.spanWisePoints{iWing}>=dividersY(iY));                            
                                else
                                    idxY=and(obj.uvlm.dataForCheb.spanWisePoints{iWing}<dividersY(iY+1),obj.uvlm.dataForCheb.spanWisePoints{iWing}>=dividersY(iY));
                                end
                                % current mode is the panels which are in
                                % this x and y section
                                modes{iWing}(:,(iX-1)*nSectionsY+iY)=and(idxX,idxY)';
                            end
                        end
                    end
                    modes=blkdiag(modes{:});
                    %% gen sys 2 reduce
                    %% reduce
                    modes=modes;
                elseif strcmp(obj.settings.mor,'bPodStrip')
                    %% generate strip modes
                    testsum=cumsum(obj.uvlm.is_te(end:-1:1));
                    testsum=testsum(end:-1:1);
                    stripModes=zeros(length(testsum),max(testsum));
                    for iStripMode=1:max(testsum)
                        stripModes(find(testsum==(max(testsum)+1)-iStripMode),iStripMode)=1;
                    end
                    %% gen sys 2 reduce
                    %% reduce
                    
                    modes=stripModes;
                    
                elseif strcmp(obj.settings.mor,'bPodPan')
                    %% generate pan modes
                    panModes=eye(size(obj.uvlm.panels,2));
                    %% gen sys 2 reduce
                    %% reduce
                    modes=panModes;
                    
                elseif strcmp(obj.settings.mor,'balRed')
                   disp('care')
                else
                   disp('care')
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
        end
        
        function obj=getLinSSM(obj, Vlin)
            obj.linV=Vlin;
            aerSSM=ss(  full(obj.tUvlmSSM.a0+Vlin*obj.tUvlmSSM.dadV),...
                        full(obj.tUvlmSSM.b0+Vlin*obj.tUvlmSSM.dbdV+Vlin^2*obj.tUvlmSSM.dbdV2),...
                        full(obj.tUvlmSSM.c0+Vlin*obj.tUvlmSSM.dcdV+Vlin^2*obj.tUvlmSSM.dcdV2+Vlin^3*obj.tUvlmSSM.dcdV3),...
                        full(obj.tUvlmSSM.d0+Vlin*obj.tUvlmSSM.dddV+Vlin^2*obj.tUvlmSSM.dddV2+Vlin^3*obj.tUvlmSSM.dddV3+Vlin^4*obj.tUvlmSSM.dddV4+Vlin^5*obj.tUvlmSSM.dddV5)...
                        );
            aerSSM.InputName=obj.tUvlmSSM.inputName;
            aerSSM.InputDelay=obj.tUvlmSSM.inputDelayFV/Vlin;
            aerSSM.InputGroup=obj.tUvlmSSM.inputGroup;
            aerSSM.OutputName=obj.tUvlmSSM.outputName;
            aerSSM.OutputGroup=obj.tUvlmSSM.outputGroup;

            %add actuators to aerSSM
            if ~isempty(obj.act.SSM)
                aerSSM=seriesPreserve(obj.act.SSM,aerSSM,{'dCs','dCsDot','dCsDotDot'});
            end
            %
            structSSM=ss(   obj.strSSM.a0+Vlin*obj.strSSM.dadV,...
                            obj.strSSM.b0+Vlin*obj.strSSM.dbdV,...
                            obj.strSSM.c0+Vlin*obj.strSSM.dcdV,...
                            obj.strSSM.d0+Vlin*obj.strSSM.dddV...
                            );
            structSSM.InputName=obj.strSSM.inputName;
            structSSM.statename=obj.strSSM.stateName;
            structSSM.OutputName=obj.strSSM.outputName;
            structSSM.InputGroup=obj.strSSM.inputGroup;
            structSSM.OutputGroup=obj.strSSM.outputGroup;

            obj.linSSM=connect(aerSSM,structSSM,[aerSSM.InputName],[aerSSM.outputName; structSSM.outputName]);
            %merge gust zones if necessary
            if obj.settings.gustInputs
                nGz= size(obj.uvlm.gustZones,1);
                combIn=ss(blkdiag(repmat(eye(3),nGz,1),repmat(eye(3),nGz,1)));
                combIn.InputName={'vGx';'vGy';'vGz';'vGxDot';'vGyDot';'vGzDot'};    
                combIn.OutputName=obj.linSSM.InputName([obj.linSSM.InputGroup.vG obj.linSSM.InputGroup.vGDot]);
                combIn.OutputGroup.vG=1:3*nGz;
                combIn.OutputGroup.vGDot=3*nGz+1:6*nGz;
                combIn.InputGroup.vG=1:3;
                combIn.InputGroup.vGDot=4:6;
                obj.linSSM=seriesPreserve(combIn,obj.linSSM,{'vG'; 'vGDot'});
            end
            %add controller if existent
            
            if ~isempty(obj.ctr)
            if ~isempty(obj.ctr.contSSM)
                ctrSSM=obj.ctr.contSSM;
                
                % add other gust inputs (x y)
%                 ctrSSM=append(ss([],zeros(0,2),[],zeros(0,2)),ctrSSM);
%                 ctrSSM.inputName=obj.linSSM.InputName(obj.linSSM.InputGroup.vG);
                % add other control outputs (the control commands which are not used by the ff controller
                nCCall=length(obj.linSSM.InputName(obj.linSSM.InputGroup.controlCommand));
                nCCctr=length(ctrSSM.OutputName);
                ctrSSM=append(ctrSSM,ss([],[],zeros(nCCall-nCCctr,0),[]));
                % find missing outputnames
                cell1=obj.linSSM.InputName(obj.linSSM.InputGroup.controlCommand);
                cell2= ctrSSM.OutputName;
                [n,m]=size(cell1);
                [~,idx]=ismember(cell1(:),cell2(:));
                [ii,~]=ind2sub([n m],idx);
                output=ii(1:n)';
                ctrSSM.OutputName(nCCctr+1:end)=[ obj.linSSM.InputName(obj.linSSM.InputGroup.controlCommand(output==0))];
                
                % resort outputs 
                cell1=obj.linSSM.InputName(obj.linSSM.InputGroup.controlCommand);
                cell2= ctrSSM.OutputName;
                [n,m]=size(cell1);
                [~,idx]=ismember(cell1(:),cell2(:));
                [ii,~]=ind2sub([n m],idx);
                output=ii(1:n)';
                ctrSSM=ctrSSM(output,:);
                ctrSSM.outputGroup.controlCommand=1:nCCall;
                
                %add additional inputs for control commands
                a=ctrSSM.a;
                b=[ctrSSM.b zeros(order(ctrSSM),size(ctrSSM,1))];
                c=[ctrSSM.c];
                d=[ctrSSM.d eye(size(ctrSSM,1))];

                ffSSM=ss(a,b,c,d);
                ffSSM.InputName=[ctrSSM.InputName; obj.linSSM.InputName(obj.linSSM.InputGroup.controlCommand)];
                ffSSM.InputGroup.controlCommand=2:1+nCCall;

                ffSSM.StateName=ctrSSM.StateName;
                ffSSM.OutputName=[ ctrSSM.OutputName; ];
                ffSSM.outputGroup.controlCommand=1:nCCall;
                %add vgdotfeedthrough
                vGft=ss(eye(6));
                vGft.InputName=[ obj.linSSM.InputName(obj.linSSM.InputGroup.vG); obj.linSSM.InputName(obj.linSSM.InputGroup.vGDot)];
                vGft.InputGroup.vG=1:3;
                vGft.InputGroup.vGDot=4:6;
                vGft.OutputName= [ obj.linSSM.InputName(obj.linSSM.InputGroup.vG); obj.linSSM.InputName(obj.linSSM.InputGroup.vGDot)];
                vGft.OutputGroup.vG=1:3;
                vGft.OutputGroup.vGDot=4:6;
                ffSSM=append(vGft,ffSSM);
                % add outputdelay to vGust output
                ffSSM.OutputDelay([ffSSM.OutputGroup.vG ffSSM.OutputGroup.vGDot])=abs(obj.settings.aoaProbePos/Vlin);
                ffSSM.OutputDelay([ffSSM.OutputGroup.vG ffSSM.OutputGroup.vGDot])=abs(30/Vlin);
                
             
                %connect to linSSM with seriesPreserve
                obj.linSSM=seriesPreserve(ffSSM,obj.linSSM,{'vG','vGDot','controlCommand'});
                
            end
            end
        end %linearize
        
        function obj=getLinSSMonlyRBMAero(obj, Vlin)
            aerSSM=ss(   full(obj.tUvlmSSM.a0+Vlin*obj.tUvlmSSM.dadV),...
                        full(obj.tUvlmSSM.b0+Vlin*obj.tUvlmSSM.dbdV+Vlin^2*obj.tUvlmSSM.dbdV2),...
                        full(obj.tUvlmSSM.c0+Vlin*obj.tUvlmSSM.dcdV+Vlin^2*obj.tUvlmSSM.dcdV2+Vlin^3*obj.tUvlmSSM.dcdV3),...
                        full(obj.tUvlmSSM.d0+Vlin*obj.tUvlmSSM.dddV+Vlin^2*obj.tUvlmSSM.dddV2+Vlin^3*obj.tUvlmSSM.dddV3+Vlin^4*obj.tUvlmSSM.dddV4+Vlin^5*obj.tUvlmSSM.dddV5)...
                        );
            aerSSM.InputName=obj.tUvlmSSM.inputName;
            aerSSM.InputGroup=obj.tUvlmSSM.inputGroup;
            aerSSM.OutputName=obj.tUvlmSSM.outputName;
            aerSSM.OutputGroup=obj.tUvlmSSM.outputGroup;

            
            obj.linSSM=aerSSM({'RBMForces','RBMmoments'},{'vB','vBDot','wB', 'wBDot','dCS_0009','dCS_0010','dCS_0011','dCS_0012'});
            
        end %linearize only aero
        
        function obj=generateTuvlm(obj)
            
            % transformed uvlm ssm restrained
            %                               __________
            %              _____    ...--->|          |--->...       _____                            
            %--modeDdot-->|     |--b/bDot->|          |--Fp-------->|     |--modalF--->
            %--modeDot--->| tAS |--s------>| uvlm.ssm |             | tSA |
            %--mode------>|_____|--vSInf-->|__________|             |_____|
           

            % add in and outputs to tAS.d and tSA.d in order to connect it to uvlmSSM with
            % preserving all other in and outputs of the uvlmSSM 
            conIdsInUvlm=[obj.uvlm.ssm.ssm.inputGroup.b obj.uvlm.ssm.ssm.inputGroup.bDot obj.uvlm.ssm.ssm.inputGroup.s obj.uvlm.ssm.ssm.inputGroup.vSInf_x obj.uvlm.ssm.ssm.inputGroup.vSInf_y obj.uvlm.ssm.ssm.inputGroup.vSInf_z];
            conIdsOutUvlm=obj.uvlm.ssm.ssm.outputGroup.Fp;

            nInUvlm=length(obj.uvlm.ssm.ssm.inputName);
            nInUvlmPrv=min(conIdsInUvlm)-1;
            nInUvlmPost=nInUvlm-nInUvlmPrv-length(conIdsInUvlm);
            nOutUvlm=length(obj.uvlm.ssm.ssm.outputName);
            nOutUvlmPrv=min(conIdsOutUvlm)-1; %<-------
            nOutUvlmPost=nOutUvlm-nOutUvlmPrv-length(conIdsOutUvlm);
            %short check if we miss something
            if ~or((nOutUvlmPost+nOutUvlmPrv+length(conIdsOutUvlm)==length(obj.uvlm.ssm.ssm.outputName)),(nInUvlmPost+nInUvlmPrv+length(conIdsInUvlm)==length(obj.uvlm.ssm.ssm.inputName)))
                disp('Warning, this does not work (generateTuvlm), some in/outputs are lost')
            end
            % expanded tAS
            tASe0=sparse([eye(nInUvlm) [zeros(nInUvlmPrv,size(obj.tAS.d0,2)); obj.tAS.d0 ; zeros(nInUvlmPost,size(obj.tAS.d0,2))]]);
            dtASedV=sparse([zeros(nInUvlm) [zeros(nInUvlmPrv,size(obj.tAS.d0,2)); obj.tAS.dddV ; zeros(nInUvlmPost,size(obj.tAS.d0,2))]]);
            %resulting inputvector
            nameIn=[obj.uvlm.ssm.ssm.inputName;obj.tAS.inputName];
            inputDelayFV=[obj.uvlm.ssm.ssm.inputDelayFV; zeros(length(obj.tAS.inputName),1)];
            pairs = [fieldnames(obj.uvlm.ssm.ssm.inputGroup), struct2cell(obj.uvlm.ssm.ssm.inputGroup); fieldnames(obj.tAS.inputGroup), struct2cell(obj.tAS.inputGroup)].';
            groupIn=struct(pairs{:});
            
            %discard inputs and groups which are substituted by the
            %structural inputs
            tASe0(:,conIdsInUvlm)=[];
            dtASedV(:,conIdsInUvlm)=[];
            nameIn(conIdsInUvlm)=[];
            groupIn=rmfield(groupIn,'b');
            groupIn=rmfield(groupIn,'bDot');
            groupIn=rmfield(groupIn,'s');
            groupIn=rmfield(groupIn,'vSInf_x');
            groupIn=rmfield(groupIn,'vSInf_y');
            groupIn=rmfield(groupIn,'vSInf_z');
            inputDelayFV(conIdsInUvlm)=[];
            
            % expanded tSA
            tSAe=sparse([eye(nOutUvlm); zeros(size(obj.tSA.d,1),nOutUvlmPrv) obj.tSA.d zeros(size(obj.tSA.d,1),nOutUvlmPost)]);
            %resulting outputvector
            nameOut=[obj.uvlm.ssm.ssm.outputName;obj.tSA.outputName];
            pairs = [fieldnames(obj.uvlm.ssm.ssm.outputGroup), struct2cell(obj.uvlm.ssm.ssm.outputGroup); fieldnames(obj.tSA.outputGroup), struct2cell(obj.tSA.outputGroup)].';
            groupOut= struct(pairs{:});
            
            %discard Fp as output
            tSAe(conIdsOutUvlm,:)=[];
            nameOut(conIdsOutUvlm)=[];
            groupOut=rmfield(groupOut,'Fp');

            % analytically connect in series with uvlmSSM to obtain one ssm with having
            % b and d as a quadratic function of V called transformed uvlm ssm
            % (tUvlmSSM)

            obj.tUvlmSSM.a0=obj.uvlm.ssm.ssm.a0;
            obj.tUvlmSSM.dadV=obj.uvlm.ssm.ssm.dadV;

            obj.tUvlmSSM.b0=obj.uvlm.ssm.ssm.b0*tASe0;
            obj.tUvlmSSM.dbdV=obj.uvlm.ssm.ssm.b0*dtASedV+obj.uvlm.ssm.ssm.dbdV*tASe0;
            obj.tUvlmSSM.dbdV2=obj.uvlm.ssm.ssm.dbdV*dtASedV;

            obj.tUvlmSSM.c0=tSAe*obj.uvlm.ssm.ssm.c0;
            obj.tUvlmSSM.dcdV=tSAe*obj.uvlm.ssm.ssm.dcdV;
            obj.tUvlmSSM.dcdV2=tSAe*obj.uvlm.ssm.ssm.dcddV;
            obj.tUvlmSSM.dcdV3=tSAe*obj.uvlm.ssm.ssm.dcdddV;

            obj.tUvlmSSM.d0=tSAe*obj.uvlm.ssm.ssm.d0*tASe0;
            obj.tUvlmSSM.dddV=tSAe*(obj.uvlm.ssm.ssm.d0*dtASedV+obj.uvlm.ssm.ssm.dddV*tASe0);
            obj.tUvlmSSM.dddV2=tSAe*(obj.uvlm.ssm.ssm.dddV*dtASedV+obj.uvlm.ssm.ssm.ddddV*tASe0);
            obj.tUvlmSSM.dddV3=tSAe*(obj.uvlm.ssm.ssm.ddddV*dtASedV+obj.uvlm.ssm.ssm.dddddV*tASe0);
            obj.tUvlmSSM.dddV4=tSAe*(obj.uvlm.ssm.ssm.dddddV*dtASedV+obj.uvlm.ssm.ssm.ddddddV*tASe0);
            obj.tUvlmSSM.dddV5=tSAe*(obj.uvlm.ssm.ssm.ddddddV*dtASedV);
            obj.tUvlmSSM.inputName=nameIn;
            obj.tUvlmSSM.inputName=nameIn;
            obj.tUvlmSSM.outputName=nameOut;
            obj.tUvlmSSM.inputDelayFV=inputDelayFV;
            
            
            fieldNames=fieldnames(groupOut);
            for iField=1:size(fieldNames,1)
                if isfield(obj.uvlm.ssm.ssm.outputGroup, fieldNames{iField})
                    groupOut.(fieldNames{iField})=find(ismember(obj.tUvlmSSM.outputName, obj.uvlm.ssm.ssm.outputName(obj.uvlm.ssm.ssm.outputGroup.(fieldNames{iField}))))';
                elseif isfield(obj.tSA.outputGroup, fieldNames{iField})
                    groupOut.(fieldNames{iField})=find(ismember(obj.tUvlmSSM.outputName, obj.tSA.outputName(obj.tSA.outputGroup.(fieldNames{iField}))))';
                else
                    disp('STOP');
                end
            end
            
            fieldNames=fieldnames(groupIn);
            for iField=1:size(fieldNames,1)
                if isfield(obj.uvlm.ssm.ssm.inputGroup, fieldNames{iField})
                    groupIn.(fieldNames{iField})=find(ismember(obj.tUvlmSSM.inputName, obj.uvlm.ssm.ssm.inputName(obj.uvlm.ssm.ssm.inputGroup.(fieldNames{iField}))))';
                elseif isfield(obj.tAS.inputGroup, fieldNames{iField})
                    groupIn.(fieldNames{iField})=find(ismember(obj.tUvlmSSM.inputName, obj.tAS.inputName(obj.tAS.inputGroup.(fieldNames{iField}))))';
                else
                    disp('STOP');
                end
            end
            
             obj.tUvlmSSM.outputGroup=groupOut;
             obj.tUvlmSSM.inputGroup=groupIn;
            
            obj.tUvlmSSM.stateName=obj.uvlm.ssm.ssm.stateName;
        end %get a parameter varying representation of the transformed uvlm (i.e. with coupling matrices)
        
        function obj=updateStrQuick(obj,aircraft,aircraft_structure)
            %update strssm and flight dynamics
            obj=obj.initSTRModel(aircraft_structure);
            %update coupling matrices by reprojection
            
            %%>>>>>>>>> not needed?
%             oldModalBase=obj.strModel.modalBase;
%             nModes=length(obj.settings.modes);
%             obj.tAS.d0=obj.tAS.d0*kron(eye(3),oldModalBase')*kron(eye(3),obj.strModel.modalBase);
%             obj.tAS.dddV=obj.tAS.dddV*kron(eye(3),oldModalBase')*kron(eye(3),obj.strModel.modalBase);
%             obj.tAS.inputName=[	cellstr([repmat('modeDdot_',nModes,1) num2str([1:nModes]','%04d')]);...
%                                 cellstr([repmat('modeDot_',nModes,1) num2str([1:nModes]','%04d')]);...
%                                 cellstr([repmat('mode_',nModes,1) num2str([1:nModes]','%04d')])];
%             obj.tAS.inputGroup.modeDdot=1:nModes;
%             obj.tAS.inputGroup.modeDot=nModes+1:2*nModes;
%             obj.tAS.inputGroup.mode=2*nModes+1:3*nModes;
%             
%             
%             obj.tSA.d=obj.strModel.modalBase'*oldModalBase*obj.tSA.d;
%             
%             obj.tSA.outputName=cellstr([repmat('modalForce_',nModes,1) num2str([1:nModes]','%04d')]);
%             obj.tSA.outputGroup.modalForce=1:nModes;
            %<<<<<<<<<<<<<<<
            
            obj=obj.initCouplingMatrices(aircraft,aircraft_structure);
            %remove Tuvlm
            obj.tUvlmSSM=[];
            %remove linSSM
            obj.linSSM=[];
            
        end
        
        function obj=runFlutterAnalysis(obj,VloopVec,plotFlag)
            obj.flutterTable.eigVal=[];
            obj.flutterTable.eigVec=[];
            obj.flutterTable.speedVec=VloopVec;
            fFlag=0;
            if obj.settings.restrained
                outputsSTR=[    obj.strSSM.outputGroup.modeDdot...
                                obj.strSSM.outputGroup.modeDot...
                                obj.strSSM.outputGroup.mode];
                inputsUVLM= [   obj.uvlm.ssm.ssm.inputGroup.b...
                                obj.uvlm.ssm.ssm.inputGroup.bDot...
                                obj.uvlm.ssm.ssm.inputGroup.s...
                                obj.uvlm.ssm.ssm.inputGroup.vSInf_x...
                                obj.uvlm.ssm.ssm.inputGroup.vSInf_y...
                                obj.uvlm.ssm.ssm.inputGroup.vSInf_z];
                outputsUVLM=[   obj.uvlm.ssm.ssm.outputGroup.Fp];
                D3=sparse((obj.tSA.d));
            else
                outputsSTR=[     obj.strSSM.outputGroup.vB...
                                obj.strSSM.outputGroup.vBDot...
                                obj.strSSM.outputGroup.wB...
                                obj.strSSM.outputGroup.wBDot...
                                obj.strSSM.outputGroup.modeDdot...
                                obj.strSSM.outputGroup.modeDot...
                                obj.strSSM.outputGroup.mode];
                inputsUVLM= [   obj.uvlm.ssm.ssm.inputGroup.vB...
                                obj.uvlm.ssm.ssm.inputGroup.vBDot...
                                obj.uvlm.ssm.ssm.inputGroup.wB...
                                obj.uvlm.ssm.ssm.inputGroup.wBDot...
                                obj.uvlm.ssm.ssm.inputGroup.b...
                                obj.uvlm.ssm.ssm.inputGroup.bDot...
                                obj.uvlm.ssm.ssm.inputGroup.s...
                                obj.uvlm.ssm.ssm.inputGroup.vSInf_x...
                                obj.uvlm.ssm.ssm.inputGroup.vSInf_y...
                                obj.uvlm.ssm.ssm.inputGroup.vSInf_z];
                outputsUVLM=[   obj.uvlm.ssm.ssm.outputGroup.RBMForces...
                                obj.uvlm.ssm.ssm.outputGroup.RBMMoments...
                                obj.uvlm.ssm.ssm.outputGroup.Fp];
                D3=sparse(blkdiag(eye(6),obj.tSA.d));
            end
            eigVec_prv=[];
            for Vloop=VloopVec
                % create matrices for certain speed
                A4=sparse(obj.strSSM.a0+obj.strSSM.dadV*Vloop);
                B4=sparse(obj.strSSM.b0+obj.strSSM.dbdV*Vloop);
                C4=sparse(obj.strSSM.c0(outputsSTR,:)+obj.strSSM.dcdV(outputsSTR,:)*Vloop);
                D4=sparse(obj.strSSM.d0(outputsSTR,:)+obj.strSSM.dddV(outputsSTR,:)*Vloop);
                
                A2=obj.uvlm.ssm.ssm.a0+obj.uvlm.ssm.ssm.dadV*Vloop;
                B2=obj.uvlm.ssm.ssm.b0(:,inputsUVLM)+obj.uvlm.ssm.ssm.dbdV(:,inputsUVLM)*Vloop+obj.uvlm.ssm.ssm.dbddV(:,inputsUVLM)*Vloop^2;
                C2=obj.uvlm.ssm.ssm.c0(outputsUVLM,:)+obj.uvlm.ssm.ssm.dcdV(outputsUVLM,:)*Vloop+obj.uvlm.ssm.ssm.dcddV(outputsUVLM,:)*Vloop^2+obj.uvlm.ssm.ssm.dcdddV(outputsUVLM,:)*Vloop^3;
                D2=obj.uvlm.ssm.ssm.d0(outputsUVLM,inputsUVLM)+obj.uvlm.ssm.ssm.dddV(outputsUVLM,inputsUVLM)*Vloop+obj.uvlm.ssm.ssm.ddddV(outputsUVLM,inputsUVLM)*Vloop^2+obj.uvlm.ssm.ssm.dddddV(outputsUVLM,inputsUVLM)*Vloop^3+obj.uvlm.ssm.ssm.ddddddV(outputsUVLM,inputsUVLM)*Vloop^4;
                if obj.settings.restrained
                    D1=sparse(obj.tAS.d0+obj.tAS.dddV*Vloop);
                else
                    D1=sparse(blkdiag(eye(12),obj.tAS.d0+obj.tAS.dddV*Vloop));
                end
                
                N=sparse(eye(size(D4,1))-D4*D3*D2*D1)^-1;
                A=[A2+B2*D1*N*D4*D3*C2 +B2*D1*N*C4;B4*D3*C2+B4*D3*D2*D1*N*D4*D3*C2 A4+B4*D3*D2*D1*N*C4];
                [eigVec, eigVal]=eig(full(A));
                %sort according to previous modeshapes
                nModeStates=size(obj.strModel.modalBase,2)*2;
                if ~isempty(eigVec_prv)
                    %calc MAC values for all eigenvector perturbations
                    MAC=zeros(size(eigVec,1));
                    for iVec=1:size(eigVec,1)
                         vec=eigVec(end-nModeStates+1:end,iVec)./norm(eigVec(end-nModeStates+1:end,iVec));
                        for iVecPrv=1:size(eigVec,1)
                            vecPrv=eigVec_prv(end-nModeStates+1:end,iVecPrv)./norm(eigVec_prv(end-nModeStates+1:end,iVecPrv));
                            MAC(iVec,iVecPrv)=(vec'*vecPrv)^2./((vec'*vec)*(vecPrv'*vecPrv));
                        end
                    end
                    %get the maximum values of every row and their idx
                    [~,idxMax]=max(abs(MAC),[],2);
                    %since some idx are multiple, find the unique idx
                    [~,idxUnique]=unique(idxMax);
                    %find the not unique ids in the idxMax 
                    idxNotUnique=setdiff( 1:size(eigVec,1),idxUnique);
                    %check which idx are missing in idxMax
                    missingIdx=setdiff(1:size(eigVec,1),idxMax);
                    % now fill the idxNotUnique in idxMax with the missingIdx
                    idxMax(idxNotUnique)=missingIdx;
                    %check if the sort vector idxMax is now fully unique
                    if length(unique(idxMax))~=unique(idxMax)
                        disp('error in sorting')
                    end
                    %now rearrange
                    [~,idxMax]=sort(idxMax);
                else
                    %initial sort by small real part
                    [~,idxMax]=sort(-real(diag(eigVal)));
                end
                obj.flutterTable.eigVal=[obj.flutterTable.eigVal diag(eigVal(idxMax,idxMax))];
                obj.flutterTable.eigVec=[obj.flutterTable.eigVec eigVec(:,idxMax)];
                if plotFlag
                    hold on;  plot(real(diag(eigVal)),imag(diag(eigVal)),'.','MarkerFaceColor',[Vloop/max(VloopVec)*0.5 1-Vloop/max(VloopVec) 0],'MarkerEdgeColor',[Vloop/max(VloopVec)*0.5  1-Vloop/max(VloopVec) 0])
                end
%               hold on;  plot(real(diag(eigVal)),imag(diag(eigVal)),'b+')

                if fFlag==0
                    if(any(real(diag(eigVal))>10^-10))
                        [~, eigVal]=eig(full(A));
                        idx=find(real(diag(eigVal))>10^-10);
                        if plotFlag
                            hold on;  plot(real(diag(eigVal)),imag(diag(eigVal)),'.','MarkerFaceColor',[Vloop/max(VloopVec)*0.5 1-Vloop/max(VloopVec) 0],'MarkerEdgeColor',[Vloop/max(VloopVec)*0.5  1-Vloop/max(VloopVec) 0])
                        end
                        Freq=abs(eigVal(idx,idx));
                        if plotFlag
                            fprintf('unstable pole @ %.02f m/s & %.02f Hz\n',Vloop,Freq(1,1)/2/pi) 
                        end
                        for iUnstable=1:length(diag(Freq))
                            if fFlag==0
                                if plotFlag
                                    plot(real(diag(eigVal(idx(iUnstable),idx(iUnstable)))),imag(diag(eigVal(idx(iUnstable),idx(iUnstable)))),'ro')
                                    text(real(diag(eigVal(idx(iUnstable),idx(iUnstable)))),imag(diag(eigVal(idx(iUnstable),idx(iUnstable)))), sprintf('first @  %.02f m/s & %.02f Hz\n',Vloop,Freq(iUnstable,iUnstable)/2/pi));
                                end
                            end
                        end
                        fFlag=1;
                    end
                end
                eigVec_prv=eigVec(:,idxMax);

            end
        end % flutter
        
        function obj=runDivergenceAnalysis(obj,VloopVec) 
            tASin=obj.tAS.inputGroup.mode;
            tASout=[obj.tAS.outputGroup.b obj.tAS.outputGroup.s obj.tAS.outputGroup.vSInf_x obj.tAS.outputGroup.vSInf_y obj.tAS.outputGroup.vSInf_z];
            
%             aeroIn=[obj.uvlm.ssm.ssm.inputGroup.vB obj.uvlm.ssm.ssm.inputGroup.wB obj.uvlm.ssm.ssm.inputGroup.b obj.uvlm.ssm.ssm.inputGroup.s obj.uvlm.ssm.ssm.inputGroup.vSInf_x obj.uvlm.ssm.ssm.inputGroup.vSInf_y obj.uvlm.ssm.ssm.inputGroup.vSInf_z];
            aeroIn=[obj.uvlm.ssm.ssm.inputGroup.b obj.uvlm.ssm.ssm.inputGroup.s obj.uvlm.ssm.ssm.inputGroup.vSInf_x obj.uvlm.ssm.ssm.inputGroup.vSInf_y obj.uvlm.ssm.ssm.inputGroup.vSInf_z];
%             aeroOut=[obj.uvlm.ssm.ssm.outputGroup.RBMForces obj.uvlm.ssm.ssm.outputGroup.RBMMoments obj.uvlm.ssm.ssm.outputGroup.Fp];
            aeroOut=[ obj.uvlm.ssm.ssm.outputGroup.Fp];
            tSAin=obj.tSA.inputGroup.Fp;
            tSAout=obj.tSA.outputGroup.modalForce;
            
%             strIn=[obj.strSSM.inputGroup.RBMForces obj.strSSM.inputGroup.RBMMoments obj.strSSM.inputGroup.modalForce];
            strIn=[obj.strSSM.inputGroup.modalForce];
%             strOut=[obj.strSSM.outputGroup.vB obj.strSSM.outputGroup.wB obj.strSSM.outputGroup.mode];  
            strOut=[obj.strSSM.outputGroup.mode];  
            strStates=1:length(obj.strSSM.stateName);
            
            
            TSA=blkdiag(obj.tSA.d);
            
            figure; hold on;
            plot(VloopVec,ones(1,length(VloopVec)),'red');
            allGains=zeros(length(strOut),length(VloopVec));
            j=1;
            for Vloop=VloopVec
                TAS=blkdiag(obj.tAS.d0(tASout,tASin)+Vloop*obj.tAS.dddV(tASout,tASin));
                aeroA=obj.uvlm.ssm.ssm.a0+obj.uvlm.ssm.ssm.dadV*Vloop;
                aeroB=obj.uvlm.ssm.ssm.b0(:,aeroIn)+obj.uvlm.ssm.ssm.dbdV(:,aeroIn)*Vloop+obj.uvlm.ssm.ssm.dbddV(:,aeroIn)*Vloop^2;
                aeroC=obj.uvlm.ssm.ssm.c0(aeroOut,:)+obj.uvlm.ssm.ssm.dcdV(aeroOut,:)*Vloop+obj.uvlm.ssm.ssm.dcddV(aeroOut,:)*Vloop^2+obj.uvlm.ssm.ssm.dcdddV(aeroOut,:)*Vloop^3;
                aeroD=obj.uvlm.ssm.ssm.d0(aeroOut,aeroIn)+obj.uvlm.ssm.ssm.dddV(aeroOut,aeroIn)*Vloop+obj.uvlm.ssm.ssm.ddddV(aeroOut,aeroIn)*Vloop^2+obj.uvlm.ssm.ssm.dddddV(aeroOut,aeroIn)*Vloop^3+obj.uvlm.ssm.ssm.ddddddV(aeroOut,aeroIn)*Vloop^4;
                kDcAero=aeroD-aeroC*aeroA^-1*aeroB;
                
                strA=obj.strSSM.a0(strStates,strStates)+obj.strSSM.dadV(strStates,strStates)*Vloop;
                strB=obj.strSSM.b0(strStates,strIn)+obj.strSSM.dbdV(strStates,strIn)*Vloop;
                strC=obj.strSSM.c0(strOut,strStates)+obj.strSSM.dcdV(strOut,strStates)*Vloop;
                strD=obj.strSSM.d0(strOut,strIn)+obj.strSSM.dddV(strOut,strIn)*Vloop;
                kDcStr=strD-strC*strA^-1*strB;
                Kae=kDcStr*TSA*kDcAero*TAS;
                allGains(:,j)=diag(Kae);
                j=j+1;
            end
            plot(VloopVec,allGains') 
            obj.divergenceData.speedVec=VloopVec;
            obj.divergenceData.gains=allGains;
        end % divergenceAnalysis
        
        function obj=runGustAnalysis(obj,gustState,Fg, nGustLengths)
            
            
            V=gustState.aerodynamic_state.V_A;
            if isempty(obj.linSSM)
                if isempty(obj.tUvlmSSM)
                    obj=obj.generateTuvlm();
                end
                obj=obj.getLinSSM(V);
            elseif obj.linV~=V
                obj=obj.getLinSSM(V);
            end
        
            if obj.settings.restrained
                outputs={'loads','mode'};
            else
                outputs={'rI','loads','mode'};
            end
            if ~isempty(obj.ctr)
                if ~isempty(obj.ctr.contSSM)
                    inputs={'vGz','vGzDot','aoaProbeIn'};            
                else
                    inputs={'vGz','vGzDot'};     
                end
            else
                inputs={'vGz','vGzDot'};                
            end
            simSSM=obj.linSSM(outputs,inputs);
            %connect vGz to aoaProbeIn when present
            if ~isempty(obj.ctr)
                sys=ss([1 0; 0 1; 1 0]);
                sys.inputName={'vGz','vGzDot'};
                sys.outputName={'vGz','vGzDot','aoaProbeIn'}; 
                
                simSSM=series(sys,simSSM);
            end
            
            gradientVec=linspace(9,107,nGustLengths); %shortest to longest gust
            timeStep=9/V/25; %s 25 time steps for shortest gust
            totalTime=107/V*3+0.01; %s
            timeVec=0:timeStep:totalTime;
            gustStartTime=timeVec(min(find(timeVec>0.01))); %s
            allLoads=zeros(length(timeVec),size(simSSM,1),length(gradientVec));
            %Reference Gust Velocity (see CS25)
            h=gustState.h;
%             if h<=4572
%                Uref= 17.07 - (17.07 - 13.41)*h/4572;
%             elseif and(h>4572,h<=18288)
%                Uref=13.41 -(13.41 - 6.36)*(h-4572)/(18288-4572);
%             elseif h>18288
%                Uref=6.36;
%             end
            
            DynPress_Dive = 26000; Mach_Dive = 0.84; % to be defined based on OAD
            
            [Uds_grid, U_sig] = ComputeGustProperties_EASA_CS25(Fg, h, gradientVec, V, DynPress_Dive, Mach_Dive); % !!!Altitude is in meters!!! % a/c speed in m/s TAS!!!
            
            for iLength=1:length(gradientVec)
                    % definition of 1-cosine gust
                    gustGradient  = gradientVec(iLength); %m
                    gustAmplitude = Uds_grid(iLength);    %m/s 
                    tGust=2*gustGradient/V;
                    omegaGust=2*pi/tGust;
                    timeVecGust=0:timeStep:tGust;
                    uGust=zeros(1,length(timeVec));
                    uGustDot=zeros(1,length(timeVec));
                    uGust(find(timeVec==gustStartTime):find(timeVec==gustStartTime)+length(timeVecGust)-1)=gustAmplitude/2*(1-cos(omegaGust*timeVecGust));
                    uGustDot(find(timeVec==gustStartTime):find(timeVec==gustStartTime)+length(timeVecGust)-1)=gustAmplitude/2*omegaGust*sin(omegaGust*timeVecGust);
                    % gust loads simulation
                    u=[uGust; uGustDot]';
                    allData(:,:,iLength)=lsim(simSSM,u,timeVec);
                    allInputs(:,:,iLength)=u;
            end
            if ~obj.settings.restrained
                obj.gustData.r=allData(:,simSSM.OutputGroup.rI,:);
            end
            obj.gustData.loads=allData(:,simSSM.OutputGroup.loads,:);
            obj.gustData.inputs=allInputs;
            obj.gustData.modes=allData(:,simSSM.OutputGroup.mode,:);
            obj.gustData.lengthVec=gradientVec;
            obj.gustData.timeVec=timeVec;
        end

        function obj=runGustAnalysisClosedLoop(obj,gustState,Fg,nGustLengths,nomLoads)
            %% controller synthesis for gla and analysis
            
            % prepare used state space model and gust velocity
            V=gustState.aerodynamic_state.V_A;
            if isempty(obj.linSSM)
                if isempty(obj.tUvlmSSM)
                    obj=obj.generateTuvlm();
                end
                obj=obj.getLinSSM(V);
            elseif obj.linV~=V
                obj=obj.getLinSSM(V);
            end
            
            % Control synthesis model generation: current implementation requires gust input to be at the end
            if obj.settings.restrained
                simSSM=obj.linSSM({'loads','mode'},{'vGz','vGzDot'});
                % simSSM_cntr = obj.linSSM({'loads'},{'dCs_0005','dCsDot_0005','dCsDotDot_0005','dCs_0006','dCsDot_0006','dCsDotDot_0006','vGz','vGzDot'});
                simSSM_cntr = obj.linSSM({'loads'},...
                                         {'dCs_0001','dCsDot_0001','dCsDotDot_0001',...
                                          'dCs_0002','dCsDot_0002','dCsDotDot_0002',...
                                          'dCs_0003','dCsDot_0003','dCsDotDot_0003',...
                                          'dCs_0004','dCsDot_0004','dCsDotDot_0004',...
                                          'dCs_0005','dCsDot_0005','dCsDotDot_0005',...
                                          'dCs_0006','dCsDot_0006','dCsDotDot_0006',...
                                          'vGz','vGzDot'});  % added by Andreas Wildschek 14.02.2018  !!!Note that gust input is at the end!!!
            else
                simSSM=obj.linSSM({'rI','loads','mode'},{'vGz','vGzDot'});
                simSSM_cntr = obj.linSSM({'loads'},...
                                         {'dCs_0001','dCsDot_0001','dCsDotDot_0001',...
                                          'dCs_0002','dCsDot_0002','dCsDotDot_0002',...
                                          'dCs_0003','dCsDot_0003','dCsDotDot_0003',...
                                          'dCs_0004','dCsDot_0004','dCsDotDot_0004',...
                                          'dCs_0005','dCsDot_0005','dCsDotDot_0005',...
                                          'dCs_0006','dCsDot_0006','dCsDotDot_0006',...
                                          'vGz','vGzDot'});  % added by Simon for free case ; no implementation of choice of cs yet: put in gustState.glaSurf
            end
            %todo: this should happen in the initialization of AeSSM
            
            %% define discrete gusts and continuous turbulence
 
            gradientVec=linspace(9,107,nGustLengths);   %shortest to longest gust (CS25: H)
            timeStep=9/V/25;                            %s 25 time steps for shortest gust
            totalTime=107/V*5+0.01;                     %s
            timeVec=0:timeStep:totalTime;
            gustStartTime=timeVec(min(find(timeVec>0.01))); %s
            allLoads=zeros(length(timeVec),size(simSSM,1),length(gradientVec));
            %Reference Gust Velocity (see CS25)
            h=gustState.h;

            DynPress_Dive = 25961; Mach_Dive = 0.85; % to be defined based on OAD
            
            [Uds_grid, U_sig] = ComputeGustProperties_EASA_CS25(Fg, h, gradientVec, V, DynPress_Dive, Mach_Dive); % !!!Altitude is in meters!!! % a/c speed in m/s TAS!!!
            
            %% generate sample continuous turbulence for optimization
            QuadraticCoherence = 0.75; % 1 is perfect coherence; 0 is no coherence; for first wing bending mode around 1 Hz, 0.75 is a realistic choice, see PhD Wildschek
            Ts = timeStep; % Controller sampling Time
            timevec_turbsim = 0:Ts:totalTime;
            SD = 1234;
            rng(SD,'twister');
            randomnumber = (1 -2*rand(length(timevec_turbsim),1));
            whitenoise1 = U_sig/sqrt(Ts)  *  (1/rms(randomnumber)*randomnumber);
            u = whitenoise1;  % vertical gust speed in m/s 
            VonKarman_Filter = sqrt(762/V)*tf([.3398*(762/V)^2 2.7478*(762/V) 1],[.1539*(762/V)^3 1.9752*(762/V)^2 2.9958*(762/V) 1]);
            TurbulenceGenerator(1:length(timevec_turbsim)) = lsim(VonKarman_Filter, u, timevec_turbsim);
            TurbulenceGeneratorDot = (diff(TurbulenceGenerator(1:end-1))+diff(TurbulenceGenerator(2:end)))/(2*Ts);
            
            u_known = sqrt(QuadraticCoherence)*whitenoise1;
                       
            SD = 5678;
            rng(SD,'twister');
            randomnumber = (1 -2*rand(length(timevec_turbsim),1));
            whitenoise2 = U_sig/sqrt(Ts)  *  (1/rms(randomnumber)*randomnumber);
            u_unknown = sqrt(1-QuadraticCoherence)*whitenoise2;  % vertical gust speed in m/s 
            TurbulenceGenerator_unknown(1:length(timevec_turbsim)) = lsim(VonKarman_Filter, u_unknown, timevec_turbsim);
            
            u_PCP = u_known+u_unknown;
            u_SCP = u;
            
            
            
            
            %% generate control synthesis model
            
            Generate_ActuatorModel; %<- TODO: the actuator models should be part of AeSSM and already connected in linSSM
            
            FullAC_Plant = simSSM_cntr*ActuatorModel_sym;   %  antisymmetric simSSM tbd
            FullAC_Plant.InputName = {'Command1','Command2','Command3','vGz','vGzDot'};
            % !!!!!!!!!!!!!!!!!!!!!!!!
           
            
            
            
            
            
            HighPass = tf([1 0],[1 4]);
            FF_delay = tf([1 -150 7500],[1 150 7500]); % t2his is the overall delay of the feed-forward controller
            
            SecondaryControlPath = HighPass*FF_delay;
            SecondaryControlPath.InputName = {'VertGustSpeed'};
            SecondaryControlPath.OutputName = {'RefSignal'};
            
            % delay between alpha probe position and gust reference point
            AlphaProbe_Position = 30; % in meters in front of first gust zone center 
            [AlphaProbe_LeadTime_NUM,AlphaProbe_LeadTime_DEN] = pade(AlphaProbe_Position/V,2);
            
            AlphaProbe_LeadTimeFilter = tf([AlphaProbe_LeadTime_NUM],[AlphaProbe_LeadTime_DEN]);
            
            DerivativeApprox = tf([1 0],[timeStep 1]);   % !!!function to take into account 'vGzDot'
            PrimaryControlPath = DerivativeApprox*AlphaProbe_LeadTimeFilter*blkdiag([tf([1],[1 0]); 1]);
            PrimaryControlPath.InputName  = {'VertGustSpeed'};
            PrimaryControlPath.OutputName = {'vGz','vGzDot'};
            
            PrimaryControlPathSimon= DerivativeApprox*blkdiag([tf([1],[1 0]); 1]); %modified PCP: no lead time filter but inputdelay 
            PrimaryControlPathSimon.InputName  = {'VertGustSpeed'};
            PrimaryControlPathSimon.OutputName = {'vGz','vGzDot'};
            PrimaryControlPathSimon.InputDelay(1)=AlphaProbe_Position/V;
%             
%             PrimaryControlPathSimon2= DerivativeApprox*blkdiag([tf([1],[1 0]); 1]); %modified PCP: no lead time filter
%             PrimaryControlPathSimon2.InputName  = {'VertGustSpeed'};
%             PrimaryControlPathSimon2.OutputName = {'vGz','vGzDot'};
%             PrimaryControlPathSimon2.OutputDelay(1)=AlphaProbe_Position/V;
%             
            PrimaryControlPath=PrimaryControlPathSimon;
            FullAC_Plant_OL = connect(PrimaryControlPath, FullAC_Plant, {'Command1','Command2','Command3','VertGustSpeed'}, FullAC_Plant.OutputName);
            
            %% Generate aircraft time responses for feed-forward control law synthesis
            % note that TimeResponses_ContTurb now offers consideration of quadratic coherence function, unlike TimeResponses_ContTurb2 which is not used any more
            [TimeResponses_ContTurb, TimeResponses_ContTurb2, turb_input] = Generate_TimeResponses_ContTurb(FullAC_Plant, FullAC_Plant_OL, AlphaProbe_LeadTimeFilter, SecondaryControlPath, u_PCP, u_SCP, VonKarman_Filter, TurbulenceGenerator, TurbulenceGeneratorDot, timevec_turbsim);

            
            for iLength=1:length(gradientVec)
                    % definition of 1-cosine gust
                    gustGradient  = gradientVec(iLength); %m
                    gustAmplitude = Uds_grid(iLength);    %m/s 
                    tGust=2*gustGradient/V;
                    omegaGust=2*pi/tGust;
                    timeVecGust=0:timeStep:tGust;
                    uGust=zeros(1,length(timeVec));
                    uGustDot=zeros(1,length(timeVec));
                    uGust(find(timeVec==gustStartTime):find(timeVec==gustStartTime)+length(timeVecGust)-1)=gustAmplitude/2*(1-cos(omegaGust*timeVecGust));
                    uGustDot(find(timeVec==gustStartTime):find(timeVec==gustStartTime)+length(timeVecGust)-1)=gustAmplitude/2*omegaGust*sin(omegaGust*timeVecGust);
                    % pure open loop gust loads simulation
                    allData(:,:,iLength)=lsim(simSSM,[uGust; uGustDot]',timeVec);
                    
                    
                    % time response data for L-inf GLA synthesis
                    % !!!Note: data already includes properties of the reference sensor
                    %[TimeResponses_1minusCosGust, gust_input1] = Generate_TimeResponses_1minusCosGust(FullAC_Plant, AlphaProbe_LeadTimeFilter, SecondaryControlPath, V, gustAmplitude, gustGradient, totalTime, timeStep, DerivativeApprox);
                    [TimeResponses_1minusCosGust, TimeResponses_1minusCosGust2, gust_input1] = Generate_TimeResponses_1minusCosGust(FullAC_Plant, FullAC_Plant_OL, AlphaProbe_LeadTimeFilter, SecondaryControlPath, V, gustAmplitude, gustGradient, totalTime, timeStep, DerivativeApprox);
                  
                    allData_cntr_OL(:,:,:,iLength) = TimeResponses_1minusCosGust;
                    allData_gust(:,iLength) = gust_input1;
                    
            end
            % simon: why is allData response generated with simSSM
            % different from allData_cntr_OL generated with FUllAC_Plant_OL
            % response at (:,:,4,:)? because FullAC_Plant_OL contains the
            % leadtimefilter
            symLineBendingMomentId = [((length(obj.settings.strLoadStations)-1)/2)*6+4]; %

            figure,plot((0:length(timeVec)-1)*timeStep,squeeze(allData(:,symLineBendingMomentId,:)));  % plot Bending Moment for a/c symmetry line for all gusts
            figure,bode(FullAC_Plant(symLineBendingMomentId,1:4)); % plot transfer function from inputs of simSSM_cntr to Bending Moment for a/c symmetry line
            
            figure,plot((0:length(timeVec)-1)*timeStep,squeeze(TimeResponses_ContTurb(:,18*6+4,:))); legend('Command1','Command2','Command3','VertGustSpeed'); xlabel('Time [s]'); %Simon: 18*6 -> wrbm?
            
%             figure,plot((0:length(timevec_turbsim)-1)*Ts,squeeze(TimeResponses_ContTurb2(:,18*6+4,:))); legend('dCs_0005','dCs_0006', 'vGz'); xlabel('Time [s]');

            %% GLA Optimization             
            % simon: why zero in the middle?
            MassCoeffBend = [linspace(0, 1, (length(obj.settings.strLoadStations)-1)/2), 0, linspace(1, 0, (length(obj.settings.strLoadStations)-1)/2)];
            figure, plot(MassCoeffBend);
            
            StaticShear   = nomLoads(3:6:end);
            StaticBending = nomLoads(4:6:end);
            StaticTorsion = nomLoads(5:6:end);
            figure,subplot(3,1,1),plot(StaticShear),title('1g Shear'), subplot(3,1,2),plot(StaticBending),title('1g Bending'),subplot(3,1,3),plot(StaticTorsion),title('1g Torsion');
            % !!!!!!!!!!!!!!!!!
%             StaticBending = zeros(1, 37);
%             MassCoeffBend = ones(1, 37);
            % !!!!!!!!!!!!!!!!!!!
            
            N = 2^4; %order of FIR filter
            Ts = timeStep; %sample time of controller must be the same as available data with current implementation
            
            %initialize controller
            h =  ones(1,(length(FullAC_Plant.InputName)-2)*N);
            h0 = zeros(1,(length(FullAC_Plant.InputName)-2)*N);
            
            
            display('*******   Optimization runnning for optimresults_GLA...');
            class optimresults_GLA;
            % Algorithm can be chosen to be 'SQP','interior-point','active-set','trust-region-reflective'
            %options = optimset('algorithm','SQP','display','iter'); 
            GLAS_ActuatorModel.LB=[GLAS_ActuatorModel.Cs_0001.min_defl GLAS_ActuatorModel.Cs_0003.min_defl GLAS_ActuatorModel.Cs_0005.min_defl];
            GLAS_ActuatorModel.UB=[GLAS_ActuatorModel.Cs_0001.max_defl GLAS_ActuatorModel.Cs_0003.max_defl GLAS_ActuatorModel.Cs_0005.max_defl];
            GLAS_ActuatorModel.maxRate=[GLAS_ActuatorModel.Cs_0001.max_rate];
            options = optimoptions('fmincon','algorithm','sqp-legacy','Display','iter-detailed');
            optimresults_GLA.x = fmincon(...
            @(h)GLA_Costfun(h, N, MassCoeffBend, StaticBending, TimeResponses_ContTurb), ...
            ...
            h0,[],[],[],[],[],[],...
            ...
            @(h)GLA_Constraintsfun(h,N,Ts,GLAS_ActuatorModel,allData_gust(:,:)),...
            ...
            options);

            display('*-*-*-*-*-*-*-* Optimization process completed !! *-*-*-*-*-*-*-*');
            
            figure, hold on;
            bode(tf([optimresults_GLA.x(1:N)],[1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0],Ts),'b');
            bode(tf([optimresults_GLA.x(N+1:2*N)],[1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0],Ts),'g');
            bode(tf([optimresults_GLA.x(2*N+1:3*N)],[1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0],Ts),'r');
            legend('Discrete transfer function of Command1 controller','Discrete transfer function of Command2 controller','Discrete transfer function of Command3 controller');
            
            FFController = [tf([optimresults_GLA.x(1:N)],[1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0],Ts); tf([optimresults_GLA.x(N+1:2*N)],[1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0],Ts); tf([optimresults_GLA.x(2*N+1:3*N)],[1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0],Ts)];
            
            opt = d2cOptions('Method','tustin','PrewarpFrequency',0.5);
            FFControllerContinuous = d2c(FFController, opt);
            FFControllerContinuous.InputName = {'RefSignal'};
            FFControllerContinuous.OutputName = {'Command1','Command2','Command3'};
            
            
            FullAC_Plant_CL = connect(PrimaryControlPath, SecondaryControlPath, FFControllerContinuous, FullAC_Plant, {'Command1','Command2','Command3','VertGustSpeed'}, FullAC_Plant.OutputName);
            
            %% quick check of result:
            ValidateSection = [(length(obj.settings.strLoadStations)-1)/2]; %
            
            BendingVector = [4:6:length(obj.settings.strLoadStations)*6];
            
            
            BendingVector_Cont_Command1 = TimeResponses_ContTurb(1:end,BendingVector,1); BendingVector_Cont_Command2  = TimeResponses_ContTurb(1:end,BendingVector,2); BendingVector_Cont_Command3 = TimeResponses_ContTurb(1:end,BendingVector,3); BendingVector_Cont_Turb = TimeResponses_ContTurb(1:end,BendingVector,end); 
           
            figure,bode(FullAC_Plant_OL(BendingVector(ValidateSection(1)),{'VertGustSpeed'})*VonKarman_Filter,'b'),hold on,bode(FullAC_Plant_CL(BendingVector(ValidateSection(1)),{'VertGustSpeed'})*VonKarman_Filter,'r');
            
            % Wing root bending moment due to continuous turbulence
            SecondaryControlPathOutput = conv(optimresults_GLA.x(1:N),BendingVector_Cont_Command1(:,ValidateSection(1)))  + conv(optimresults_GLA.x(N+1:2*N),BendingVector_Cont_Command2(:,ValidateSection(1))) + conv(optimresults_GLA.x(N+1:2*N),BendingVector_Cont_Command3(:,ValidateSection(1)));
            
            % Method 1 cutting off the convolution vector "SecondaryControlPathOutput"
            WRB_turb_CL = (...
                            (...
                                (...
                                    abs(StaticBending(ValidateSection(1))'  )   ...
                                    + sqrt(...
                                            sum(...
                                                    (  ... 
                                                        (1/size(squeeze([BendingVector_Cont_Turb(:,ValidateSection(1))]),1)) ...
                                                        * ...
                                                        ( ...
                                                            [BendingVector_Cont_Turb(:,ValidateSection(1))]   ...
                                                            + ...
                                                            SecondaryControlPathOutput(1:length(timeVec)) ...
                                                        ) .^2 ...
                                                    )...
                                                )...
                                        )   ...
                                )...
                            ) ...
                          );
            WRB_turb_OL = (((abs(StaticBending(ValidateSection(1))')   + sqrt(sum((   (1/size(squeeze([BendingVector_Cont_Turb(:,ValidateSection(1))]),1))  * ([BendingVector_Cont_Turb(:,ValidateSection(1))] ) .^2)))   )));
            % Method 2 with zero padding
            WRB_turb_CL = (...
                            (...
                                (...
                                    abs(StaticBending(ValidateSection(1))')...
                                    +...
                                    sqrt( ...
                                        sum(...
                                            (...
                                                (1/size(squeeze([BendingVector_Cont_Turb(:,ValidateSection(1));zeros(N-1,1)]),1))...
                                                * ...
                                                (...
                                                    [BendingVector_Cont_Turb(:,ValidateSection(1));zeros(N-1,1)]   ...
                                                    + SecondaryControlPathOutput ...
                                                 ) .^2 ...
                                             )...
                                            )...
                                        )...
                                 )...
                             )...
                         );
            WRB_turb_OL = (...
                            (...
                                (...
                                    abs(StaticBending(ValidateSection(1))')...
                                    + ...
                                    sqrt( ...
                                        sum( ...
                                            (...
                                                (1/size(squeeze([BendingVector_Cont_Turb(:,ValidateSection(1));zeros(N-1,1)]),1))...
                                                * ...
                                                (...
                                                    [BendingVector_Cont_Turb(:,ValidateSection(1));zeros(N-1,1)] ...
                                                ) .^2 ...
                                             ) ...
                                            ) ...
                                        )...
                                ) ...
                            ) ...
                          );
            %method 3 without the last 16 elements
            
            WRB_nom=(StaticBending(ValidateSection(1))');
            WRB_OL_t=[BendingVector_Cont_Turb(:,ValidateSection(1))] ;
            WRB_CL_t=WRB_OL_t + SecondaryControlPathOutput(1:end-N+1);
            WRB_turb_OL=max(abs(WRB_nom+rms(WRB_OL_t)),abs(WRB_nom-rms(WRB_OL_t)));
            WRB_turb_CL=max(abs(WRB_nom+rms(WRB_CL_t)),abs(WRB_nom-rms(WRB_CL_t)));

            display('*-*-*-*-*-*-*-* Percentage of continuous turbulence loads reduction at the section "ValidateSection" !! *-*-*-*-*-*-*-*');
            
            WRB_turb_CL
            WRB_turb_OL
            
            disp(['bending load at wing root reduced by ' num2str(abs(WRB_turb_OL - WRB_turb_CL)/abs(WRB_turb_OL)*100) '%'])
            
            %% generate closed loop time responses           
            [TimeResponses_ContTurb_CL, TimeResponses_ContTurb2_CL, turb_input_CL] = Generate_TimeResponses_ContTurb(FullAC_Plant, FullAC_Plant_CL, AlphaProbe_LeadTimeFilter, SecondaryControlPath, u_PCP, u_SCP, VonKarman_Filter, TurbulenceGenerator, TurbulenceGeneratorDot, timevec_turbsim);
            
            figure,plot((0:length(timeVec)-1)*timeStep,squeeze(TimeResponses_ContTurb(:,BendingVector(ValidateSection(1)),end))),hold on,plot((0:length(timeVec)-1)*timeStep,squeeze(TimeResponses_ContTurb_CL(:,BendingVector(ValidateSection(1)),end)),'r'); legend('open loop','close loop'); xlabel('Time [s]');
            
            clear OpenLoopLoads; clear ClosedLoopLoads;
            OpenLoopLoads   = rms(TimeResponses_ContTurb(:,:,4));
            ClosedLoopLoads = rms(TimeResponses_ContTurb_CL(:,:,4));
            figure,grid on;
            plot(abs(StaticBending)'+OpenLoopLoads(4:6:end),'b'),hold on,plot(abs(StaticBending)'+ClosedLoopLoads(4:6:end),'r'); xlabel('wing span section [-]'); ylabel('Total bending load [Nm]'); legend('GLA off', 'GLA on');
            
            for iLength=1:length(gradientVec)
                    % definition of 1-cosine gust
                    gustGradient  = gradientVec(iLength); %m
                    gustAmplitude = Uds_grid(iLength);    %m/s 
                    
                    % Note: "TimeResponses_1minusCosGust2_CL" is only for verification of code but otherwise don't use it!!!
                    [TimeResponses_1minusCosGust_CL, TimeResponses_1minusCosGust2_CL, gust_input1_CL] = Generate_TimeResponses_1minusCosGust(FullAC_Plant, FullAC_Plant_CL, AlphaProbe_LeadTimeFilter, SecondaryControlPath, V, gustAmplitude, gustGradient, totalTime, timeStep, DerivativeApprox);
                    allData_cntr_CL(:,:,:,iLength) = TimeResponses_1minusCosGust_CL;
                    allData_gust_CL(:,iLength) = gust_input1_CL;
                    
                    % plot the commands in 1-cosine gust
                    Command1 = conv(optimresults_GLA.x(1:N),allData_gust(:,iLength));
                    Command2 = conv(optimresults_GLA.x(N+1:2*N),allData_gust(:,iLength));
                    Command3 = conv(optimresults_GLA.x(2*N+1:3*N),allData_gust(:,iLength));
                    figure(99), hold on
                    subplot(3,1,1), hold on, plot(0:Ts:totalTime, Command1(1:length(timeVec))); xlabel('Time []'), ylabel('Command1 []');
                    subplot(3,1,2), hold on, plot(0:Ts:totalTime, Command2(1:length(timeVec))); xlabel('Time []'), ylabel('Command2 []');
                    subplot(3,1,3), hold on, plot(0:Ts:totalTime, Command3(1:length(timeVec))); xlabel('Time []'), ylabel('Command3 []');
                
                    figure(100),grid on;
                    hold on,plot(StaticBending(ValidateSection(1))'   + allData_cntr_OL(:,BendingVector(ValidateSection(1)),4,iLength)), hold on, plot(StaticBending(ValidateSection(1))'   + allData_cntr_CL(:,BendingVector(ValidateSection(1)),4,iLength),'--');
                    hold on,plot(StaticBending(ValidateSection(1))'   - allData_cntr_OL(:,BendingVector(ValidateSection(1)),4,iLength)), hold on, plot(StaticBending(ValidateSection(1))'   - allData_cntr_CL(:,BendingVector(ValidateSection(1)),4,iLength),'--');
            end
            
            % these plots are only to verify the code
            figure,plot(squeeze(TimeResponses_1minusCosGust_CL(:,BendingVector(ValidateSection(1)),1:4))), hold on,plot(squeeze(TimeResponses_1minusCosGust2_CL(:,BendingVector(ValidateSection(1)),1:4)),'--')
            figure,plot(squeeze(TimeResponses_1minusCosGust(:,BendingVector(ValidateSection(1)),1:4))), hold on,plot(squeeze(TimeResponses_1minusCosGust2(:,BendingVector(ValidateSection(1)),1:4)),'--')
            

            %% write results
            if ~obj.settings.restrained
                obj.gustData.r=allData(:,simSSM.OutputGroup.rI,:);
            end
            obj.gustData.loads                           = allData(:,simSSM.OutputGroup.loads,:);
            obj.gustData.modes                           = allData(:,simSSM.OutputGroup.mode,:);
            obj.gustData.lengthVec                       = gradientVec;
            obj.gustData.timeVec                         = timeVec;
            
            obj.gustData.FullAC_Plant                    = FullAC_Plant;
            obj.gustData.FullAC_Plant_OL                 = FullAC_Plant_OL;
            obj.gustData.FullAC_Plant_CL                 = FullAC_Plant_CL;
            obj.gustData.TimeResponses_ContTurb          = TimeResponses_ContTurb;
            obj.gustData.TimeResponses_ContTurb_CL       = TimeResponses_ContTurb_CL;
            obj.gustData.TimeResponses_1minusCosGust_CL  = TimeResponses_1minusCosGust_CL;
            
            obj.gustData.optimresults_GLA                = optimresults_GLA;
            
            obj.gustData.loads_OL                        = allData_cntr_OL;
            obj.gustData.discrete_gusts                  = allData_gust;
            obj.gustData.loads_cntr_CL                   = allData_cntr_CL;
            obj.gustData.discrete_gusts_CL               = allData_gust_CL;
            
        end
        
        function obj=runContTurbAnalysis(obj, contTurbState,Fg)
            V=contTurbState.aerodynamic_state.V_A;
            %time step is determined as in discretegust:
            Ts=9/V/25; %s 25 time steps for shortest gust
            totalTime=107/V*5+0.01; %s
            tVec= 0:Ts:totalTime;
            %gen signal or take from settings
            if or(isempty(obj.settings.contTurbSignal),~isequal(tVec,tVec))
                SD = 5678;
                rng(SD,'twister');
                randomnumber2 = (1 -2*rand(length(tVec),1));
                SD = 1234;
                rng(SD,'twister');
                randomnumber1 = (1 -2*rand(length(tVec),1));
                %store for further analyses
                obj.settings.contTurbSignal.rand1=randomnumber1;
                obj.settings.contTurbSignal.rand2=randomnumber2;
                obj.settings.contTurbSignal.tVec=tVec;
            else
                %get previously generated signals
                randomnumber1=obj.settings.contTurbSignal.rand1;
                randomnumber2=obj.settings.contTurbSignal.rand2; 
                if length(tVec)~=length(obj.settings.contTurbSignal.tVec)
                    disp('smth went wrong, different time vector in prev cont turb sim')
                end
            end
            %get Usig from CS25
            h=contTurbState.h;

            DynPress_Dive = 25961; Mach_Dive = 0.85; % to be defined based on OAD
            gradientVec=[];
            [~, U_sig] = ComputeGustProperties_EASA_CS25(Fg, h, gradientVec, V, DynPress_Dive, Mach_Dive); % !!!Altitude is in meters!!! % a/c speed in m/s TAS!!!
            
            %coherence between SCP and PCP
            QuadraticCoherence = 0.75; % 1 is perfect coherence; 0 is no coherence; for first wing bending mode around 1 Hz, 0.75 is a realistic choice, see PhD Wildschek

            whitenoise1 = U_sig/sqrt(Ts)  *  (1/rms(randomnumber1)*randomnumber1);
            u = whitenoise1;  % vertical gust speed in m/s 
            u_known = sqrt(QuadraticCoherence)*whitenoise1;

            whitenoise2 = U_sig/sqrt(Ts)  *  (1/rms(randomnumber2)*randomnumber2);
            u_unknown = sqrt(1-QuadraticCoherence)*whitenoise2;  % vertical gust speed in m/s 

            u_PCP = u_known+u_unknown;
            u_SCP = u;
                
            
            %get sim model
            if isempty(obj.linSSM)
                if isempty(obj.tUvlmSSM)
                    obj=obj.generateTuvlm();
                end
                obj=obj.getLinSSM(V);
            elseif obj.linV~=V
                obj=obj.getLinSSM(V);
            end
            if obj.settings.restrained
                outputs={'loads','mode'};
            else
                outputs={'rI','loads','mode'};
            end
            if ~isempty(obj.ctr)
                if ~isempty(obj.ctr.contSSM)
                    inputs={'vGz','vGzDot','aoaProbeIn'};            
                else
                    inputs={'vGz','vGzDot'};     
                end
            else
                inputs={'vGz','vGzDot'};                
            end
            simSSM=obj.linSSM(outputs,inputs);
            % derivative approx and integrator
            DerivativeApprox = tf([1 0],[Ts 1]);   % !!!function to take into account 'vGzDot'
            PCP = blkdiag(DerivativeApprox*blkdiag([tf([1],[1 0]); 1]),1);
            PCP.InputName  = {'VertGustSpeed', 'VertGustSpeedSCP'};
            PCP.OutputName = {'vGz','vGzDot', 'aoaProbeIn'};
            

            % add pcp to simSSM
            simSSM=series(PCP,simSSM,'name');
           
            
            %von Karman Filter
            VonKarman_Filter = sqrt(762/V)*tf([.3398*(762/V)^2 2.7478*(762/V) 1],[.1539*(762/V)^3 1.9752*(762/V)^2 2.9958*(762/V) 1]);
            
            %simulate            
            [allData] = lsim(simSSM*blkdiag(VonKarman_Filter,VonKarman_Filter), [u_PCP u_SCP], tVec);
            
            %store contTurbData 
            if ~obj.settings.restrained
                obj.gustData.contTurb.r=allData(:,simSSM.OutputGroup.rI);
            end
            obj.gustData.contTurb.inputs= [u_PCP u_SCP];
            obj.gustData.contTurb.loads=allData(:,simSSM.OutputGroup.loads);
            obj.gustData.contTurb.modes=allData(:,simSSM.OutputGroup.mode);
            obj.gustData.contTurb.timeVec=tVec;
        end
        
        function obj=optimFFController(obj, state,Fg,nGustLengths, nomLoads)
			% get required variables
            V=state.aerodynamic_state.V_A;
            %time step is determined as in discretegust:
            Ts=9/V/25; %s 25 time steps for shortest gust
            % clear ffcontroller
            obj.ctr.contSSM=[];
            obj.linSSM=[];
            %regen linSSM without ctr
            if isempty(obj.linSSM)
                if isempty(obj.tUvlmSSM)
                    obj=obj.generateTuvlm();
                end
                obj=obj.getLinSSM(V);
            elseif obj.linV~=V
                obj=obj.getLinSSM(V);
            end
            %get controlled cs
            nCs=length(state.glaSurf);
            outputName=cell(nCs,1);
            for iCs=1:nCs
                outputName{iCs}=[state.glaSurf{iCs} '_sym'];
            end
            % init ctr class
            obj.ctr=class_FF_Controller(outputName);
            N=obj.ctr.order;
            % prepare struct with min max def and rates
			GLAS_ActuatorModel.UB=state.glaUB;
            GLAS_ActuatorModel.LB=state.glaLB;
            GLAS_ActuatorModel.maxRate=obj.act.rateLimit;
            % coefficients for loads required for objective function (i.e.
            % what mass benefit can be achived by a load reduction of this
            % station)
            MassCoeffBend = [linspace(0, 1, (length(obj.settings.strLoadStations)-1)/2), 0, linspace(1, 0, (length(obj.settings.strLoadStations)-1)/2)];
            %Static Bending from nominal loads ; This is required so that
            %the optimizer knows if it is better to reduce or to maximize
            %the dynamic load at a certain station  
            StaticBending = nomLoads(4:6:end);            
            
            %% generate time responses (gust to load)
            obj=obj.runContTurbAnalysis(state,1);
            %% generate models for time responses (control to load)
            %scp model
            scp=series(obj.ctr.delay, obj.ctr.highPass);
            %von karman filter
            VonKarman_Filter = sqrt(762/V)*tf([.3398*(762/V)^2 2.7478*(762/V) 1],[.1539*(762/V)^3 1.9752*(762/V)^2 2.9958*(762/V) 1]);
            %% uSCP  gen
            randomnumber1=obj.settings.contTurbSignal.rand1;
            randomnumber2=obj.settings.contTurbSignal.rand2;
            DynPress_Dive = 25961; Mach_Dive = 0.85; % to be defined based on OAD
            gradientVec=linspace(9,107,nGustLengths);
            
            [Uds_grid, U_sig] = ComputeGustProperties_EASA_CS25(Fg, state.h, gradientVec, V, DynPress_Dive, Mach_Dive); % !!!Altitude is in meters!!! % a/c speed in m/s TAS!!!
            
            %coherence between SCP and PCP
            QuadraticCoherence = 0.75; % 1 is perfect coherence; 0 is no coherence; for first wing bending mode around 1 Hz, 0.75 is a realistic choice, see PhD Wildschek

            whitenoise1 = U_sig/sqrt(Ts)  *  (1/rms(randomnumber1)*randomnumber1);
            u = whitenoise1;  % vertical gust speed in m/s 
            u_known = sqrt(QuadraticCoherence)*whitenoise1;

            whitenoise2 = U_sig/sqrt(Ts)  *  (1/rms(randomnumber2)*randomnumber2);
            u_unknown = sqrt(1-QuadraticCoherence)*whitenoise2;  % vertical gust speed in m/s 

            u_PCP = u_known+u_unknown;
            u_SCP = u;
            %% contTurb  sim
            dataContTurb=zeros(length(obj.gustData.contTurb.timeVec),length(obj.linSSM.OutputGroup.loads),nCs);
            for iCs=1:nCs
                dataContTurb(:,:,iCs)=lsim(obj.linSSM({'loads'},outputName(iCs))*scp*VonKarman_Filter,u_SCP,obj.gustData.contTurb.timeVec);
            end
            %discrete gust inputs
            timeVec=obj.gustData.contTurb.timeVec;
            allData_gust=zeros(length(obj.gustData.contTurb.timeVec),nGustLengths);
            for iLength=1:length(gradientVec)
                % definition of 1-cosine gust
                gustGradient  = gradientVec(iLength); %m
                gustAmplitude = Uds_grid(iLength);    %m/s 
                tGust=2*gustGradient/V;
                omegaGust=2*pi/tGust;
                timeVecGust=0:Ts:tGust;
                uGust=zeros(1,length(timeVec));
                uGust(1:length(timeVecGust))=gustAmplitude/2*(1-cos(omegaGust*timeVecGust));
                % sim
                allData_gust(:,iLength)=lsim(scp,uGust,timeVec);
            end
            
            
            
            
            
            
            
            %% merge time responses 
            % include time shift for gust to load responses (aoaprobe)
            % TODO: check effect on result if this matters or not...
            delay=abs(obj.settings.aoaProbePos/V);
            delay=abs(30/V);
            [~,delayIdx]=min(abs(timeVec-delay));
            dataContTurb(delayIdx:end,:,nCs+1)=obj.gustData.contTurb.loads(1:end-delayIdx+1,:);
			%% optim
            h0=0*ones(1,nCs*N);
%             h0=linspace(1,2,nCs*N);
            %% test
%             GLA_Costfun(h0, N, MassCoeffBend, StaticBending, dataContTurb)
%             GLA_CostfunOld(h0, N, MassCoeffBend, StaticBending, dataContTurb)
%             GLA_Constraintsfun(h0,N,Ts,GLAS_ActuatorModel,allData_gust);
            %%
			display('*******   Optimization runnning for optimresults_GLA...');
            options =optimoptions('fmincon','algorithm','sqp-legacy','Display','iter-detailed','MaxFunctionEvaluations',10000, 'ConstraintTolerance',1e-4);
            
            x = fmincon(...
            @(h)GLA_Costfun(h, N, MassCoeffBend, StaticBending, dataContTurb), ...
            ...
            h0,[],[],[],[],[],[],...
            ...
            @(h)GLA_Constraintsfun(h,N,Ts,GLAS_ActuatorModel,allData_gust(:,:)),...
            ...
            options);
            display('*-*-*-*-*-*-*-* Optimization process completed !! *-*-*-*-*-*-*-*');
			%% set controller coefficients
			% set
			obj.ctr.coeff=x;
			% gen ssm
			obj.ctr=obj.ctr.genSSM(Ts);
            obj.linSSM=[];
			
            
        end
        function obj=plotActDeflectionAndRates(obj, inputsignal,timevec)
            commandData=obj.ctr.getCommands(inputsignal,timevec);
            for iSim=1:size(inputsignal,2)
                actData(:,:,iSim)=lsim(obj.act.SSM(:,obj.ctr.outputName),commandData(:,:,iSim),timevec);
            end
            activeCSid=find(sum(actData(:,1:3:end,1)~=0)~=0);
            figure;
            for iCs=1:length(activeCSid)
                subplot(length(activeCSid),2,((iCs)-1)*2+1)
                plot(timevec, squeeze(actData(:,(activeCSid(iCs)-1)*3+1,:)))
                xlabel('time [s]')
                ylabel(obj.act.SSM.outputname((activeCSid(iCs)-1)*3+1));
                subplot(length(activeCSid),2,((iCs)-1)*2+2)
                plot(timevec, squeeze(actData(:,(activeCSid(iCs)-1)*3+2,:)))
                xlabel('time [s]')
                ylabel(obj.act.SSM.outputname((activeCSid(iCs)-1)*3+2));
            end
            suptitle('Deflection and Rate of Control Surfaces')
        end
    end
    
end


 
