%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%
classdef class_UVLM_SSM
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        % class which contains the LPV Kernel
        LPVKernel
        
        %parameter varying ssm
        ssm
        
        %reduced parameter varying ssm
        redSSM
        
        %reduction basis
        redBase
        
        %inverse of reduction basis
        redBaseInv
        
        %linearized state space model
        linSSM
        
        %speed for linearized model
        Vlin
        
        %matrix which transforms vector of panel velocities and their
        %derivatives to b,bDot and Vs, stored as sys
        vpT
        
        %matrix which transforms gust inputs and their
        %derivatives to b,bDot and Vs, stored as sys
        vgT
        
        %matrix which transforms vector of rigid body motion and their
        %derivatives to b,bDot and Vs, stored as sys
        rbmT
        
        %matrix which transforms vector of cs deflections and their
        %derivatives to b,bDot and Vs, stored as sys
        csT
        
        % matrix which transforms output of LPV kernel to Forces, stored as
        % sys
        linOutT
        
        % sys which transforms output of linOutT to integrated Forces, 
        intOutT
        
        % matrix which transforms output of linOutputT to total Force
        % Outputs
        totDerOutT
        
        %matrix which transforms output of linOutputT to sectional pressure
        %outputs (for kier comparison)
        sectPressOutT
        
        %settings class
        settings
    end
    
    methods
        function obj=class_UVLM_SSM(solver)
            obj.settings=class_UVLM_SSM_settings();
            obj.rbmT=obj.generateRigidBodyMotionInputTransformation(solver);
            obj.vpT=obj.generatePanelVelocityInputTransformation(solver);
            obj.LPVKernel=class_UVLM_LPVKernel(solver);
            obj=obj.generateLinearOutputTransformation(solver);
            obj=obj.generateIntegratedOutputTransformation(solver);
            obj.csT=obj.generateCsInputTransformation(solver);
            obj.vgT=obj.generateGustInputTransformation(solver);
            if obj.settings.redOrder~=0
                obj.LPVKernel=obj.LPVKernel.reduce(obj.settings.redOrder, obj.settings.redMethod);
            end
            %%% first part: influences on forces due to changes in Boundary
            %%% condition
            %            ______          _______             ______         
            %--b/bDot--->|    |--b/bDot->|     |--gbeff----->|    |--Fp----->
            %--V/omega-->|    |          | LPV |--gbDoteff-->|    |--RBMFo-->
            %--Vp------->| Ti |          | Ker-|--vSind_x--->| To |--RBMMo-->
            %--Vg------->|    |          | nel |--vSind_y--->|    |--HingeM->
            %--dCS------>|____|          |_____|--vSind_z--->|____|  
            %
            %%% second part: influences on forces due to changes in
            %%% inflow condition -> this tilts the force vector
            %            ______                              ______                    _
            %--b/bDot--->|    |------------vSinf_x---------->|    |--Fp----->           |
            %--V/omega-->|    |------------vSinf_y---------->|    |--RBMFo-->           |
            %--Vp------->| Ti2|------------vSinf_z---------->| To2|--RBMMo-->           |- this gets D2
            %--Vg------->|    |                              |    |--HingeM->           |
            %--dCS------>|____|                              |____|                    _|
            %    
            %%% third part: influences on forces due to changes in
            %%% segment orientation > this tilts the force vector
            %                                                ______         
            %--s-------------------------------------------->|    |--Fp----->
            %--vSinf_x-------------------------------------->|    |--RBMFo-->
            %--vSinf_y-------------------------------------->| To3|--RBMMo-->
            %--vSinf_z-------------------------------------->|    |--HingeM->
            %                                                |____|          


            % assemble Ti
            Ti0=[eye(size(obj.LPVKernel.b0,2)), ...
                obj.rbmT.d([obj.rbmT.OutputGroup.b obj.rbmT.OutputGroup.bDot],:), ...
                obj.vpT.d([obj.rbmT.OutputGroup.b obj.rbmT.OutputGroup.bDot],:), ...
                obj.vgT.d([obj.rbmT.OutputGroup.b obj.rbmT.OutputGroup.bDot],:), ...
                obj.csT.d0([obj.rbmT.OutputGroup.b obj.rbmT.OutputGroup.bDot],:)];
            dTidV=[ zeros(size(obj.LPVKernel.b0,2)), ...
                    zeros(size(obj.rbmT.d([obj.rbmT.OutputGroup.b obj.rbmT.OutputGroup.bDot],:))), ...
                    zeros(size(obj.vpT.d([obj.rbmT.OutputGroup.b obj.rbmT.OutputGroup.bDot],:))), ...
                    zeros(size(obj.vgT.d([obj.rbmT.OutputGroup.b obj.rbmT.OutputGroup.bDot],:))), ...
                    obj.csT.dddV([obj.rbmT.OutputGroup.b obj.rbmT.OutputGroup.bDot],:)];
            % assemble Ti2
            Ti20=[zeros(size(obj.LPVKernel.b0,2)/2*3,size(obj.LPVKernel.b0,2)), ...
                obj.rbmT.d([obj.rbmT.OutputGroup.vSInf_x obj.rbmT.OutputGroup.vSInf_y  obj.rbmT.OutputGroup.vSInf_z ],:), ...
                obj.vpT.d([obj.rbmT.OutputGroup.vSInf_x obj.rbmT.OutputGroup.vSInf_y obj.rbmT.OutputGroup.vSInf_z],:), ...
                obj.vgT.d([obj.rbmT.OutputGroup.vSInf_x obj.rbmT.OutputGroup.vSInf_y obj.rbmT.OutputGroup.vSInf_z],:), ...
                obj.csT.d0([obj.rbmT.OutputGroup.vSInf_x obj.rbmT.OutputGroup.vSInf_y obj.rbmT.OutputGroup.vSInf_z],:)];
            dTi2dV=[ zeros(size(obj.LPVKernel.b0,2)/2*3,size(obj.LPVKernel.b0,2)), ...
                   zeros(size(obj.rbmT.d([obj.rbmT.OutputGroup.vSInf_x obj.rbmT.OutputGroup.vSInf_y  obj.rbmT.OutputGroup.vSInf_z ],:))), ...
                zeros(size(obj.vpT.d([obj.rbmT.OutputGroup.vSInf_x obj.rbmT.OutputGroup.vSInf_y obj.rbmT.OutputGroup.vSInf_z],:))), ...
                zeros(size(obj.vgT.d([obj.rbmT.OutputGroup.vSInf_x obj.rbmT.OutputGroup.vSInf_y obj.rbmT.OutputGroup.vSInf_z],:))), ...
                obj.csT.dddV([obj.rbmT.OutputGroup.vSInf_x obj.rbmT.OutputGroup.vSInf_y obj.rbmT.OutputGroup.vSInf_z],:)];
                        
            % assemble To
            To2in=[obj.linOutT.InputGroup.gbeff obj.linOutT.InputGroup.gbDoteff obj.linOutT.InputGroup.vSInd_x obj.linOutT.InputGroup.vSInd_y obj.linOutT.InputGroup.vSInd_z];
            To2out=obj.linOutT.OutputGroup.Fp;
            To0=[   obj.linOutT.d0(To2out,To2in);...
                    obj.intOutT.d*obj.linOutT.d0(To2out,To2in)];
            dTodV=[ obj.linOutT.dddV(To2out,To2in);...
                    obj.intOutT.d*obj.linOutT.dddV(To2out,To2in)];
            dToddV=[ obj.linOutT.ddddV(To2out,To2in);...
                     obj.intOutT.d*obj.linOutT.ddddV(To2out,To2in)];
            % assemble To2
            To2In=[obj.linOutT.InputGroup.vSInf_x obj.linOutT.InputGroup.vSInf_y obj.linOutT.InputGroup.vSInf_z];
            To2Out=[obj.linOutT.OutputGroup.Fp];
            To20=[   obj.linOutT.d0(To2Out,To2In);...
                    obj.intOutT.d*obj.linOutT.d0(To2Out,To2In)];
            dTo2dV=[ obj.linOutT.dddV(To2Out,To2In);...
                    obj.intOutT.d*obj.linOutT.dddV(To2Out,To2In)];
            dTo2ddV=[ obj.linOutT.ddddV(To2Out,To2In);...
                    obj.intOutT.d*obj.linOutT.ddddV(To2Out,To2In)];
            % assemble To3 (this one feeds through vsInf and s for
            % structural coupling)
            To3In=[obj.linOutT.InputGroup.s obj.linOutT.InputGroup.vSInf_x obj.linOutT.InputGroup.vSInf_y obj.linOutT.InputGroup.vSInf_z];
            To3Out=[obj.linOutT.OutputGroup.Fp];
            To30=[   obj.linOutT.d0(To3Out,To3In);...
                    obj.intOutT.d*obj.linOutT.d0(To3Out,To3In)];
            dTo3dV=[ obj.linOutT.dddV(To3Out,To3In);...
                    obj.intOutT.d*obj.linOutT.dddV(To3Out,To3In)];
            dTo3ddV=[ obj.linOutT.ddddV(To3Out,To3In);...
                    obj.intOutT.d*obj.linOutT.ddddV(To3Out,To3In)];
            %assemble D2
            D20=To20*Ti20;
            D2dV=To20*dTi2dV + dTo2dV*Ti20 ;
            D2ddV=dTo2dV*dTi2dV + dTo2ddV* Ti20;
            D2dddV= dTo2ddV* dTi2dV;
              
            nP=length(obj.LPVKernel.InputGroup.b);
            nStates=length(obj.LPVKernel.stateNames);
            % get parts of c and d of LPVkernel which give gbeff/gbdoteff
            lpvOut=[obj.LPVKernel.OutputGroup.gbeff obj.LPVKernel.OutputGroup.gbDoteff obj.LPVKernel.OutputGroup.vSInd_x obj.LPVKernel.OutputGroup.vSInd_y obj.LPVKernel.OutputGroup.vSInd_z];
            c0=obj.LPVKernel.c0(lpvOut,:);
            dcdV=obj.LPVKernel.dcdV(lpvOut,:);
            d0=obj.LPVKernel.d0(lpvOut,:);
            dddV=obj.LPVKernel.dddV(lpvOut,:);
            % assemble a
            obj.ssm.a0=obj.LPVKernel.a0;
            obj.ssm.dadV=obj.LPVKernel.dadV;
            % assemble b
            obj.ssm.b0=[obj.LPVKernel.b0*Ti0 zeros(nStates,nP*6)];
            obj.ssm.dbdV=[obj.LPVKernel.b0*dTidV + obj.LPVKernel.dbdV*Ti0 zeros(nStates,nP*6)];
            obj.ssm.dbddV=[obj.LPVKernel.dbdV*dTidV zeros(nStates,nP*6)];
            %assemble c
            obj.ssm.c0=To0*c0;
            obj.ssm.dcdV=To0*dcdV+dTodV*c0;
            obj.ssm.dcddV=dTodV*dcdV+dToddV*c0;
            obj.ssm.dcdddV=dToddV*dcdV;
            %assemble d
            obj.ssm.d0=[To0*d0*Ti0 + D20 To30];
            obj.ssm.dddV=[To0*d0*dTidV + To0*dddV*Ti0 + dTodV*d0*Ti0 + D2dV dTo3dV];
            obj.ssm.ddddV=[To0*dddV*dTidV + dTodV*dddV*Ti0 + dTodV*d0*dTidV + dToddV*d0*Ti0 + D2ddV dTo3ddV];
            obj.ssm.dddddV=[dTodV*dddV*dTidV + dToddV*d0*dTidV + dToddV*dddV*Ti0 + D2dddV zeros(size(To30))];
            obj.ssm.ddddddV=[dToddV*dddV*dTidV zeros(size(To30))];
            
            
            
            obj.ssm.inputName=[ obj.LPVKernel.inputNames;...
                                obj.rbmT.InputName;...
                                obj.vpT.InputName;...
                                obj.vgT.InputName;...
                                obj.csT.InputName;...
                                obj.linOutT.InputName( [obj.linOutT.InputGroup.s obj.linOutT.InputGroup.vSInf_x obj.linOutT.InputGroup.vSInf_y obj.linOutT.InputGroup.vSInf_z]);];
            nGin=length(obj.vgT.InputGroup.vG);
            nCsIn=length(obj.csT.InputGroup.dCs);
            obj.ssm.inputGroup.b=(1:nP);
            obj.ssm.inputGroup.bDot=(nP+1:2*nP);
            obj.ssm.inputGroup.vB=(2*nP+1:2*nP+3);
            obj.ssm.inputGroup.vBDot=(2*nP+3+1:2*nP+3+3);
            obj.ssm.inputGroup.wB=(2*nP+3+3+1:2*nP+3+3+3);
            obj.ssm.inputGroup.wBDot=(2*nP+3+3+3+1:2*nP+3+3+3+3);
            obj.ssm.inputGroup.vP=(2*nP+12+1:2*nP+12+3*nP);
            obj.ssm.inputGroup.vPDot=(2*nP+12+3*nP+1:2*nP+12+3*nP+3*nP);
            obj.ssm.inputGroup.vG=(8*nP+12+1:8*nP+12+nGin);
            obj.ssm.inputGroup.vGDot=(8*nP+12+nGin+1:8*nP+12+2*nGin);
            obj.ssm.inputGroup.dCs=(8*nP+12+2*nGin+1:8*nP+12+2*nGin+nCsIn);
            obj.ssm.inputGroup.dCsDot=(8*nP+12+2*nGin+nCsIn+1:8*nP+12+2*nGin+2*nCsIn);
            obj.ssm.inputGroup.dCsDotDot=(8*nP+12+2*nGin+2*nCsIn+1:8*nP+12+2*nGin+3*nCsIn);
            obj.ssm.inputGroup.s=(8*nP+12+2*nGin+3*nCsIn+1:8*nP+12+2*nGin+3*nCsIn+3*nP);
            obj.ssm.inputGroup.vSInf_x=(8*nP+12+2*nGin+3*nCsIn+3*nP+1:8*nP+12+2*nGin+3*nCsIn+4*nP);
            obj.ssm.inputGroup.vSInf_y=(8*nP+12+2*nGin+3*nCsIn+4*nP+1:8*nP+12+2*nGin+3*nCsIn+5*nP);
            obj.ssm.inputGroup.vSInf_z=(8*nP+12+2*nGin+3*nCsIn+5*nP+1:8*nP+12+2*nGin+3*nCsIn+6*nP);

            obj.ssm.outputName=[obj.linOutT.OutputName(obj.linOutT.OutputGroup.Fp); obj.intOutT.OutputName];
            obj.ssm.outputGroup.Fp=obj.linOutT.OutputGroup.Fp;
            obj.ssm.outputGroup.RBMForces=nP*3+1:nP*3+3;
            obj.ssm.outputGroup.RBMMoments=nP*3+3+1:nP*3+3+3;
            obj.ssm.outputGroup.HingeMoments=nP*3+6+1:nP*3+6+nCsIn;
            obj.ssm.stateName=obj.LPVKernel.stateNames;
            % inputdelays for gust zones
             %inputdelay has to be divided by the speed at which the model
            %is used!!
            obj.ssm.inputDelayFV=zeros(length(obj.ssm.inputName),1);
            obj.ssm.inputDelayFV([obj.ssm.inputGroup.vG obj.ssm.inputGroup.vGDot])=obj.vgT.InputDelayFV;
        end
        
        function obj=linearize(obj,V) 
            obj.linSSM=ss(  full(obj.ssm.a0+obj.ssm.dadV*V),...
                            full(obj.ssm.b0+obj.ssm.dbdV*V+ obj.ssm.dbddV*V^2),...
                            full(obj.ssm.c0+obj.ssm.dcdV*V+ obj.ssm.dcddV*V^2+ obj.ssm.dcdddV*V^3),...
                            full(obj.ssm.d0+obj.ssm.dddV*V+ obj.ssm.ddddV*V^2+obj.ssm.dddddV*V^3 +obj.ssm.ddddddV*V^4));
            obj.linSSM.InputName=obj.ssm.inputName;
            obj.linSSM.InputGroup=obj.ssm.inputGroup;
            obj.linSSM.InputDelay=obj.ssm.inputDelayFV/V;
            obj.linSSM.OutputName=obj.ssm.outputName;
            obj.linSSM.OutputGroup=obj.ssm.outputGroup;
            obj.linSSM.StateName=obj.ssm.stateName;
        end
        function obj=generateLinearOutputTransformation(obj, solver)
            %matrix which transforms gb and gbeff of LPV kernel to Forces,
            %stored as sys
            V=solver.Uinf;
            %number bound panels
            nB=size(solver.panels,2);
            %Mach Correction Factor
            if solver.Ma_corr<1
                betaInf=solver.Ma_corr;
                betaInf=1; %<- comment: if the grid is stretched, no further postprocessing is needed for the forces (at least do the values then match the ones from compressible Nastran solutions (steady))
            else
                betaInf=1;
            end
            
            Areas=solver.area';
            rho=solver.state.rho_air;
            %matrix to repeat 3x1 vector npanel times
            Ihat=repmat(eye(3),nB,1);
            %matrix containing the skew symmetric matrices of s
            segmentVec=((0.75*solver.grid(:,solver.panels(2,:))+0.25*solver.grid(:,solver.panels(3,:)))-(0.75*solver.grid(:,solver.panels(1,:))+0.25*solver.grid(:,solver.panels(4,:))))';
            ScrossCellskew=skewsymmatrixCell(num2cell(segmentVec,2));
            Scross=blkdiag(ScrossCellskew{:})';
            %vector with stacked normal vectors
            NzeroHat=reshape(solver.colloc_nvec,nB*3,1);
            %matrix to repeat gamma vector elements 3 times
            gbrep=kron(eye(nB),[1 1 1]');
            %factor for gbdotf
            gbdotF=gbrep*Areas.*NzeroHat*rho*1;
            if solver.settings.indComp
                %get induced downwash due to static loading
%                 solver=solver.f_solve;
%                 solver=solver.determine_boundary_conditions;
                %boundary condition at linearization point; for rbmT, V
                %needs to be in the flight dynamic coordinate system with
                %z pointing downwards, x forward and y to right wing
                bLin=obj.rbmT.d(obj.rbmT.OutputGroup.b ,obj.rbmT.InputGroup.vB)*([-1 1 -1].*V)';
                % get dcgain of this boundary condition onto gbeff and
                % vsind
                inLPV=obj.LPVKernel.InputGroup.b;
                outLPV=[obj.LPVKernel.OutputGroup.gbeff obj.LPVKernel.OutputGroup.vSInd_x obj.LPVKernel.OutputGroup.vSInd_y obj.LPVKernel.OutputGroup.vSInd_z];
                Gsteady=(obj.LPVKernel.d0(outLPV,inLPV)+obj.LPVKernel.dddV(outLPV,inLPV)...
                    -(...
                    (obj.LPVKernel.c0(outLPV,:)+obj.LPVKernel.dcdV(outLPV,:))*...
                    (obj.LPVKernel.a0(:,:)+obj.LPVKernel.dadV(:,:))^-1*...
                    (obj.LPVKernel.b0(:,inLPV)+obj.LPVKernel.dbdV(:,inLPV))...
                    ))*-bLin;
            
                %induced velocity at linearization point
                Vind=[Gsteady(nB+1:2*nB)'; Gsteady(2*nB+1:3*nB)'; Gsteady(3*nB+1:4*nB)'];
                %factor for effective vorticity input WITH INDUCED COMPONENT
                gbF=Scross*(Ihat*V'+reshape(Vind,size(Vind,1)*size(Vind,2),1))*rho*1;
                %factor for segment direction input
                gb0VsiCrossCellSkew=skewsymmatrixCell(num2cell((Vind'+V).*repmat(Gsteady(1:nB),1,3),2));
                gb0VsiCross=blkdiag(gb0VsiCrossCellSkew{:});
                sF=rho*gb0VsiCross;
                %factor for segment center velocity input
                gb0SCrossCellSkew=skewsymmatrixCell(num2cell(segmentVec.*repmat(Gsteady(1:nB),1,3),2));
                gb0SCross=blkdiag(gb0SCrossCellSkew{:})';
                vF=rho*gb0SCross;
            else
                %factor for gbf WITHOUT INDUCED COMPONENT
                gbF=Scross*Ihat*V'*rho*1;
                sF=sparse(3*nB,3*nB);
                vF=sparse(3*nB,3*nB);
            end
            
            dvFdV=vF./norm(V);
            dsFdV2=sF./norm(V)^2;
            dgbFdV=gbF./norm(V);
            if strcmp(obj.settings.force,'steady')
                gbdotF=0*gbdotF;
            elseif strcmp(obj.settings.force,'unsteady')
                gbF=0*gbF;
                dgbFdV=0*dgbFdV;
            end
                
            %overall transformation matrix
            if obj.settings.splitForce
            else
                
                obj.linOutT.d0=sparse([diag(gbF*0)*gbrep -diag(gbdotF)*gbrep 0*dsFdV2 0*vF(:,1:3:end) 0*vF(:,2:3:end) 0*vF(:,3:3:end) 0*vF(:,1:3:end) 0*vF(:,2:3:end) 0*vF(:,3:3:end)]);
                obj.linOutT.dddV=sparse([diag(dgbFdV)*gbrep -diag(gbdotF*0)*gbrep 0*dsFdV2 dvFdV(:,1:3:end) dvFdV(:,2:3:end) dvFdV(:,3:3:end)  dvFdV(:,1:3:end) dvFdV(:,2:3:end) dvFdV(:,3:3:end)]);
                obj.linOutT.ddddV=sparse([diag(gbF*0)*gbrep -diag(gbdotF*0)*gbrep dsFdV2 0*vF(:,1:3:end) 0*vF(:,2:3:end) 0*vF(:,3:3:end) 0*vF(:,1:3:end) 0*vF(:,2:3:end) 0*vF(:,3:3:end)]);
                obj.linOutT.InputName=[     cellstr([repmat('gbEff_',nB,1) num2str([1:nB]','%04d')]);...
                                                    cellstr([repmat('gbEffDot_',nB,1) num2str([1:nB;]','%04d')]);...
                                                    cellstr([repmat(['s_x_'; 's_y_';'s_z_'],nB,1) num2str(reshape(repmat(1:nB,3,1),3*nB,1),'%04d')]);...
                                                    cellstr([repmat('vSInf_x_',nB,1) num2str([1:nB]','%04d')]);...
                                                    cellstr([repmat('vSInf_y_',nB,1) num2str([1:nB]','%04d')]);...
                                                    cellstr([repmat('vSInf_z_',nB,1) num2str([1:nB]','%04d')]);
                                                    cellstr([repmat('vSInd_x_',nB,1) num2str([1:nB]','%04d')]);...
                                                    cellstr([repmat('vSInd_y_',nB,1) num2str([1:nB]','%04d')]);...
                                                    cellstr([repmat('vSInd_z_',nB,1) num2str([1:nB]','%04d')]);];
                obj.linOutT.OutputName=[    cellstr([repmat(['Fp_x_'; 'Fp_y_';'Fp_z_'],nB,1) num2str(reshape(repmat(1:nB,3,1),3*nB,1),'%04d')]);];
                obj.linOutT.OutputGroup.Fp=[1:3*nB];
                obj.linOutT.InputGroup.gbeff=[1:nB];
                obj.linOutT.InputGroup.gbDoteff=[nB+1:2*nB];
                obj.linOutT.InputGroup.s=[2*nB+1:5*nB];
                obj.linOutT.InputGroup.vSInf_x=[5*nB+1:6*nB];
                obj.linOutT.InputGroup.vSInf_y=[6*nB+1:7*nB];
                obj.linOutT.InputGroup.vSInf_z=[7*nB+1:8*nB];
                obj.linOutT.InputGroup.vSInd_x=[8*nB+1:9*nB];
                obj.linOutT.InputGroup.vSInd_y=[9*nB+1:10*nB];
                obj.linOutT.InputGroup.vSInd_z=[10*nB+1:11*nB];
            end
        end
        function obj=generateIntegratedOutputTransformation(obj,solver)
            %number bound panels
            nB=size(solver.panels,2);
            %matrix to repeat 3x1 vector npanel times
            Ihat=repmat(eye(3),nB,1);
            
            %Rxs contains the cross product matrices of the vectors from the
            %reference point to the panel fvap
            rSPCellSkew=skewsymmatrixCell(num2cell(solver.r',2)); 
            Rxs=[rSPCellSkew{:}]';
            T=[Ihat';Rxs';solver.Rcss];  
            % transform T and Td so that the integrated coefficients are
            % put out in the flight dynamic coordinate frame with x
            % pointing rear, z pointing down and y pointing to the right
            % wing
            T2=blkdiag(diag([-1 1 -1,-1 1 -1]),eye(size(solver.Rcss,1))); 
            T=T2*T;
            
            obj.intOutT.d=sparse(T);
            obj.intOutT.InputName=cellstr([repmat(['Fp_x_'; 'Fp_y_';'Fp_z_'],nB,1) num2str(reshape(repmat(1:nB,3,1),3*nB,1),'%04d')]);
            nCs=size(solver.Rcs,1);
            if nCs>0
                obj.intOutT.OutputName=['X'; 'Y'; 'Z'; 'L'; 'M'; 'N'; cellstr([repmat('Chinge_',nCs,1) num2str([1:nCs]','%04d')])];
            else
                obj.intOutT.OutputName={'X'; 'Y'; 'Z'; 'L'; 'M'; 'N';};
            end
            obj.intOutT.OutputGroup.RBMForces=1:3;
            obj.intOutT.OutputGroup.RBMMoments=4:6;
            obj.intOutT.OutputGroup.HingeMoments=6+1:6+nCs;
            obj.intOutT.InputGroup.Fp=1:3*nB;
        end
        
        function obj=generateDerivativeOutputTransformation(obj,solver,V)
            %in this function, the totForceOutT sys is generated
            %outputs are all summed forces and moments in body fixed frame
            %as well as the hinge moments normed to qinf and S_ref
            %inputs are all total forces on the panels
            
            %number bound panels
            nB=size(solver.panels,2);
            %matrix to repeat 3x1 vector npanel times
            Ihat=repmat(eye(3),nB,1);
            
            %Rxs contains the cross product matrices of the vectors from the
            %reference point to the panel fvap
            rSPCellSkew=skewsymmatrixCell(num2cell(solver.r',2)); 
            Rxs=[rSPCellSkew{:}]';
            T=[Ihat';Rxs';solver.Rcss];
            Td=T/(solver.state.rho_air/2*norm(V)^2*solver.reference.S_ref);
            % transform T and Td so that the integrated coefficients are
            % put out in the flight dynamic coordinate frame with x
            % pointing rear, z pointing down and y pointing to the right
            % wing
            T2=blkdiag(diag([-1 1 -1,-1 1 -1]),eye(size(solver.Rcss,1))); 
            T=T2*T;
            Td=T2*Td;
            if obj.settings.splitForce
                obj.totDerOutT=ss([blkdiag(Td,Td,Td);blkdiag(T,T,T)]);
                obj.totDerOutT.InputName=[    cellstr([repmat(['Fp_x_'; 'Fp_y_';'Fp_z_'],nB,1) num2str(reshape(repmat(1:nB,3,1),3*nB,1),'%04d')]);...
                                                    cellstr([repmat(['Fp_st_x_'; 'Fp_st_y_';'Fp_st_z_'],nB,1) num2str(reshape(repmat(1:nB,3,1),3*nB,1),'%04d')]);...
                                                    cellstr([repmat(['Fp_ust_x_'; 'Fp_ust_y_';'Fp_ust_z_'],nB,1) num2str(reshape(repmat(1:nB,3,1),3*nB,1),'%04d')]);];
                nCs=size(solver.Rcs,1);
                if nCs > 0
                    obj.totDerOutT.OutputName=['CX'; 'CY'; 'CZ'; 'CL'; 'CM'; 'CN'; cellstr([repmat('Chinge_',nCs,1) num2str([1:nCs]','%04d')]);...
                                            'CX_st'; 'CY_st'; 'CZ_st'; 'CL_st'; 'CM_st'; 'CN_st'; cellstr([repmat('Chinge_st_',nCs,1) num2str([1:nCs]','%04d')]);...
                                            'CX_ust'; 'CY_ust'; 'CZ_ust'; 'CL_ust'; 'CM_ust'; 'CN_ust'; cellstr([repmat('Chinge_ust_',nCs,1) num2str([1:nCs]','%04d')]);
                                            'X'; 'Y'; 'Z'; 'L'; 'M'; 'N'; cellstr([repmat('Mhinge_',nCs,1) num2str([1:nCs]','%04d')]);...
                                            'X_st'; 'Y_st'; 'Z_st'; 'L_st'; 'M_st'; 'N_st'; cellstr([repmat('Mhinge_st_',nCs,1) num2str([1:nCs]','%04d')]);...
                                            'X_ust'; 'Y_ust'; 'Z_ust'; 'L_ust'; 'M_ust'; 'N_ust'; cellstr([repmat('Mhinge_ust_',nCs,1) num2str([1:nCs]','%04d')]);];
                else
                    obj.totDerOutT.OutputName={ 'CX'; 'CY'; 'CZ'; 'CL'; 'CM'; 'CN'; ...
                                                'CX_st'; 'CY_st'; 'CZ_st'; 'CL_st'; 'CM_st'; 'CN_st';...
                                                'CX_ust'; 'CY_ust'; 'CZ_ust'; 'CL_ust'; 'CM_ust'; 'CN_ust';...
                                                'X'; 'Y'; 'Z'; 'L'; 'M'; 'N'; ...
                                                'X_st'; 'Y_st'; 'Z_st'; 'L_st'; 'M_st'; 'N_st';...
                                                'X_ust'; 'Y_ust'; 'Z_ust'; 'L_ust'; 'M_ust'; 'N_ust';};
                    
                end
                obj.totDerOutT.OutputGroup.Derivatives=[1:(6+nCs)];
                obj.totDerOutT.OutputGroup.Derivatives_st=[(6+nCs)+1:2*(6+nCs)];
                obj.totDerOutT.OutputGroup.Derivatives_ust=[2*(6+nCs)+1:3*(6+nCs)]; 
                
                obj.totDerOutT.OutputGroup.RBMForces=[(6+nCs)*3+1:(6+nCs)*3+3];
                obj.totDerOutT.OutputGroup.RBMMoments=[(6+nCs)*3+4:(6+nCs)*3+6];
                
                obj.totDerOutT.OutputGroup.RBMForces_st=[(6+nCs)*4+1:(6+nCs)*4+3];
                obj.totDerOutT.OutputGroup.RBMMoments_st=[(6+nCs)*4+4:(6+nCs)*4+6];
                
                obj.totDerOutT.OutputGroup.RBMForces_ust=[(6+nCs)*5+1:(6+nCs)*5+3];
                obj.totDerOutT.OutputGroup.RBMMoments_ust=[(6+nCs)*5+4:(6+nCs)*5+6];
                
                obj.totDerOutT.InputGroup.Fp=[1:3*nB];
                obj.totDerOutT.InputGroup.Fp_st=[3*nB+1:6*nB];
                obj.totDerOutT.InputGroup.Fp_ust=[6*nB+1:9*nB];
                if nCs > 0
                    obj.totDerOutT.OutputGroup.HingeMoments=[(6+nCs)*3+6+1:(6+nCs)*3+6+nCs];
                    obj.totDerOutT.OutputGroup.HingeMoments_st=[(6+nCs)*4+6+1:(6+nCs)*4+6+nCs];
                    obj.totDerOutT.OutputGroup.HingeMoments_ust=[(6+nCs)*5+6+1:(6+nCs)*5+6+nCs];
                end
            else
                obj.totDerOutT=ss(T);
                obj.totDerOutT.InputName=cellstr([repmat(['Fp_x_'; 'Fp_y_';'Fp_z_'],nB,1) num2str(reshape(repmat(1:nB,3,1),3*nB,1),'%04d')]);
                nCs=size(solver.Rcs,1);
                obj.totDerOutT.OutputName=['CX'; 'CY'; 'CZ'; 'CL'; 'CM'; 'CN'; cellstr([repmat('Chinge_',nCs,1) num2str([1:nCs]','%04d')])];
                obj.totDerOutT.OutputGroup.Derivatives=[1:(6+nCs)];
                
                obj.totDerOutT.OutputGroup.RBMForces=[(6+nCs)+1:(6+nCs)+3];
                obj.totDerOutT.OutputGroup.RBMMoments=[(6+nCs)+4:(6+nCs)+6];
                obj.totDerOutT.OutputGroup.HingeMoments=[(6+nCs)+6+1:(6+nCs)+6+nCs];
                obj.totDerOutT.InputGroup.Fp=[1:3*nB];
            end
        end
        function obj=generateSectionalPressureOutputTransformation(obj,solver,V)
            %in this function, the sectPressOutT sys is generated
            %outputs are all spanwised summed forces and normed to qinf 
            %inputs are all total forces on the panels
            
            %number bound panels
            nB=size(solver.panels,2);
            %number spanwise panels (both left/right)
            nSpan=sum(solver.is_te);
            %number chordwise panels 
            nChord=nB/nSpan;
            %areas of the panels
            Areas=solver.area';
            %matrix
            T=[repmat(kron(eye(nChord),[0 0 1]),1,nSpan)]*diag(reshape(repmat(1./Areas,1,3)',nB*3,1))/(solver.state.rho_air/2*norm(V)^2)/(nSpan);
            if obj.settings.splitForce
                obj.sectPressOutT=ss(blkdiag(T,T,T));
                obj.sectPressOutT.InputName=[    cellstr([repmat(['Fp_x_'; 'Fp_y_';'Fp_z_'],nB,1) num2str(reshape(repmat(1:nB,3,1),3*nB,1),'%04d')]);...
                                                    cellstr([repmat(['Fp_st_x_'; 'Fp_st_y_';'Fp_st_z_'],nB,1) num2str(reshape(repmat(1:nB,3,1),3*nB,1),'%04d')]);...
                                                    cellstr([repmat(['Fp_ust_x_'; 'Fp_ust_y_';'Fp_ust_z_'],nB,1) num2str(reshape(repmat(1:nB,3,1),3*nB,1),'%04d')]);];
                obj.sectPressOutT.OutputName=[cellstr([repmat('sectPress_',nChord,1) num2str([1:nChord]','%04d')]);...
                    cellstr([repmat('sectPress_st_',nChord,1) num2str([1:nChord]','%04d')]);...
                    cellstr([repmat('sectPress_ust_',nChord,1) num2str([1:nChord]','%04d')]);];
                obj.sectPressOutT.OutputGroup.SectPress=[1:nChord];
                obj.sectPressOutT.OutputGroup.SectPress_st=[nChord+1:2*nChord];
                obj.sectPressOutT.OutputGroup.SectPress_ust=[2*nChord+1:3*nChord];
                obj.sectPressOutT.InputGroup.Fp=[1:3*nB];
                obj.sectPressOutT.InputGroup.Fp_st=[3*nB+1:6*nB];
                obj.sectPressOutT.InputGroup.Fp_ust=[6*nB+1:9*nB];
            else
                obj.sectPressOutT=ss(T);
                obj.sectPressOutT.InputName=cellstr([repmat(['Fp_x_'; 'Fp_y_';'Fp_z_'],nB,1) num2str(reshape(repmat(1:nB,3,1),3*nB,1),'%04d')]);
                obj.sectPressOutT.OutputName=[cellstr([repmat('sectPress_',nChord,1) num2str([1:nChord]','%04d')])];
                obj.sectPressOutT.OutputGroup.SectPress=[1:nChord];
                obj.sectPressOutT.InputGroup.Fp=[1:3*nB];
            end
        end
        function obj=reduce(obj,redBase, redBaseInv)
            obj.redBase=redBase;
            obj.redBaseInv=redBaseInv;
            obj.ssm.a0=sparse(redBaseInv*obj.ssm.a0*obj.redBase);
            obj.ssm.dadV=sparse(redBaseInv*obj.ssm.dadV*obj.redBase);
            obj.ssm.b0=sparse(redBaseInv*obj.ssm.b0);
            obj.ssm.dbdV=sparse(redBaseInv*obj.ssm.dbdV);
            obj.ssm.dbddV=sparse(redBaseInv*obj.ssm.dbddV);
            obj.ssm.c0=sparse(obj.ssm.c0*obj.redBase);
            obj.ssm.dcdV=sparse(obj.ssm.dcdV*obj.redBase);
            obj.ssm.dcddV=sparse(obj.ssm.dcddV*obj.redBase);
            obj.ssm.dcdddV=sparse(obj.ssm.dcdddV*obj.redBase);
            obj.ssm.stateName=cellstr([repmat('LagRed_',size(obj.redBase,2),1) num2str([1:size(obj.redBase,2)]','%04d')]);
        end
    end
    methods(Static) % static
     
        function rbmT=generateRigidBodyMotionInputTransformation(solver)
            % this matrix/this system transforms a velocity and rotation of
            % the aircraft onto the boundary condition which is implied at
            % each panel and its time derivative. Example: the body moves
            % with a velocity vB=[vB_x vB_y vB_z] body axes. ; the
            % resulting boundary condition for a panel with a normal vector
            % n=[nx ny nz] defined in the body reference frame would be
            % b=-(vB_x*nx+vB_y*ny+vB_z*nz).  

            %number bound panels 
            nB=size(solver.panels,2);            
            %matrix to repeat 3x1 vector npanel times
            Ihat=repmat(eye(3),nB,1);
            %normal vector blockdiagonal matrix
            normalvectorCell=num2cell(solver.colloc_nvec,1);
            NZero=blkdiag(normalvectorCell{:});
            %matrix containing the cross product matrices of the vectors
            %from the reference point to the collocation points
            rDwnPCellSkew=skewsymmatrixCell(num2cell(solver.colloc',2));
            Rx=[rDwnPCellSkew{:}]';
            %matrix containing the cross product matrices of the vectors
            %from the reference point to the FVAP points
            rSPCellSkew=skewsymmatrixCell(num2cell(solver.r',2)); 
            Rxs=[rSPCellSkew{:}]';
            %matrices selecting x,y and z components
            Sx=kron(eye(nB),[1 0 0]);
            Sy=kron(eye(nB),[0 1 0]);
            Sz=kron(eye(nB),[0 0 1]);
            
            %here change the sign of the x and z components to have a
            %flight mechanic coordinate system with x pointing rear and z
            %pointing down while y still pointing twoards right wing as in
            %the daedalus body frame
            T=diag([-1 1 -1,-1 1 -1,-1 1 -1,-1 1 -1,]);
            
            rbmT.d=sparse([-NZero'*Ihat ...
                    zeros(nB,3)...
                    -NZero'*Rx...
                    zeros(nB,3);
                    %
                    zeros(nB,3)...
                    -NZero'*Ihat ...
                    zeros(nB,3)...
                    -NZero'*Rx;...
                    
                    -Sx*Ihat ...
                    zeros(nB,3)...
                    -Sx*Rxs...
                    zeros(nB,3)...
                    %
                    -Sy*Ihat ...
                    zeros(nB,3)...
                    -Sy*Rxs...
                    zeros(nB,3);...
                    %
                    -Sz*Ihat ...
                    zeros(nB,3)...
                    -Sz*Rxs...
                    zeros(nB,3);...
                ]*T);
            rbmT.InputName={'vB_x';'vB_y';'vB_z';'vBDot_x';'vBDot_y';'vBDot_z';'wB_x';'wB_y';'wB_z'; 'wBDot_x';'wBDot_y';'wBDot_z'};
            
            rbmT.OutputName=[   cellstr([repmat('b_',nB,1) num2str([1:nB]','%04d')]);...
                                cellstr([repmat('bDot_',nB,1) num2str([1:nB]','%04d')]);...
                                cellstr([repmat('vSInf_x_',nB,1) num2str([1:nB]','%04d')]);...
                                cellstr([repmat('vSInf_y_',nB,1) num2str([1:nB]','%04d')]);...
                                cellstr([repmat('vSInf_z_',nB,1) num2str([1:nB]','%04d')]);];
            rbmT.OutputGroup.b=[1:nB];
            rbmT.OutputGroup.bDot=[nB+1:2*nB];
            rbmT.OutputGroup.vSInf_x=[nB*2+1:3*nB];
            rbmT.OutputGroup.vSInf_y=[nB*3+1:4*nB];
            rbmT.OutputGroup.vSInf_z=[nB*4+1:5*nB];
            rbmT.InputGroup.vB=[1:3];
            rbmT.InputGroup.vBDot=[4:6];
            rbmT.InputGroup.wB=[7:9];
            rbmT.InputGroup.wBDot=[10:12];
        end
        function vpT=generatePanelVelocityInputTransformation(solver)
            
            %number bound panels
            nB=size(solver.panels,2);            
            %normal vector blockdiagonal matrix
            normalvectorCell=num2cell(solver.colloc_nvec,1);
            NZero=blkdiag(normalvectorCell{:});
           
            %matrices selecting x,y and z components
            Sx=kron(eye(nB),[1 0 0]);
            Sy=kron(eye(nB),[0 1 0]);
            Sz=kron(eye(nB),[0 0 1]);
            
            
            vpT.d=sparse([blkdiag(-NZero',-NZero');[[Sx ;Sy;Sz] zeros(nB*3,nB*3)]]);
            
            
            vpT.InputName=[ cellstr([repmat(['vP_x_'; 'vP_y_';'vP_z_'],nB,1) num2str(reshape(repmat(1:nB,3,1),3*nB,1),'%04d')]);...
                            cellstr([repmat(['vPDot_x_'; 'vPDot_y_';'vPDot_z_'],nB,1) num2str(reshape(repmat(1:nB,3,1),3*nB,1),'%04d')]);];
            
            vpT.OutputName=[   cellstr([repmat('b_',nB,1) num2str([1:nB]','%04d')]);...
                                cellstr([repmat('bDot_',nB,1) num2str([1:nB]','%04d')]);...
                                cellstr([repmat('vSInf_x_',nB,1) num2str([1:nB]','%04d')]);...
                                cellstr([repmat('vSInf_y_',nB,1) num2str([1:nB]','%04d')]);...
                                cellstr([repmat('vSInf_z_',nB,1) num2str([1:nB]','%04d')]);];
            vpT.OutputGroup.b=[1:nB];
            vpT.OutputGroup.bDot=[nB+1:2*nB];
            vpT.OutputGroup.vSInf_x=[nB*2+1:3*nB];
            vpT.OutputGroup.vSInf_y=[nB*3+1:4*nB];
            vpT.OutputGroup.vSInf_z=[nB*4+1:5*nB];
            vpT.InputGroup.vP=[1:3*nB];
            vpT.InputGroup.vPDot=[3*nB+1:6*nB];
        end
        function csT=generateCsInputTransformation(solver)
            Vlin=(solver.Uinf)/norm(solver.Uinf);
            %number bound panels
            nB=size(solver.panels,2);    
            %number of control surface inputs
            nCs=size(solver.r_cs,3);
            %normal vector blockdiagonal matrix
            normalvectorCell=num2cell(solver.colloc_nvec,1);
            NZero=blkdiag(normalvectorCell{:});
            %matrix to repeat 3x1 vector npanel times
            Ihat=repmat(eye(3),nB,1);
                       
            %vector containing the distances of the collocation points on 
            %the control surface panels totheir hingepoint
            Rcs2=repmat(sum(reshape(solver.r_cs,nB*3,1,nCs),3),1,nCs); %careful no overlapping allowed!!!
             Rcs2=squeeze(reshape(solver.r_cs,nB*3,1,nCs)); %better
            %same for segment center points
            Rcs2s=repmat(sum(reshape(solver.r_css,nB*3,1,nCs),3),1,nCs);
             Rcs2s=squeeze(reshape(solver.r_css,nB*3,1,nCs)); %better
            %matrices selecting x,y and z components
            Sx=kron(eye(nB),[1 0 0]);
            Sy=kron(eye(nB),[0 1 0]);
            Sz=kron(eye(nB),[0 0 1]);
            %prealloc
            csT0=sparse(nB*5,nCs*3);
            csTdV=sparse(nB*5,nCs*3);
            for iCs=1:nCs
                csT0(:,(iCs-1)*3+1:(iCs)*3)=[    (solver.Khat(:,:,iCs)*NZero)'*Ihat*Vlin'*0 ... %<- influence due to deflection
                            -NZero'*solver.Khat(:,:,iCs)*Rcs2(:,iCs),...                             %<- influence on b due to speed
                            zeros(nB,1); ...
                            %
                            zeros(nB,1)...
                            (solver.Khat(:,:,iCs)*NZero)'*Ihat*Vlin'*0 ...
                            -NZero'*solver.Khat(:,:,iCs)*Rcs2(:,iCs); ...
                            %
                            zeros(nB,1)...            <- missing part 
                            -Sx*solver.Khat(:,:,iCs)*Rcs2s(:,iCs)...
                            zeros(nB,1); ...
                            %
                            zeros(nB,1)...            <- missing part 
                            -Sy*solver.Khat(:,:,iCs)*Rcs2s(:,iCs)...
                            zeros(nB,1); ...
                            %
                            zeros(nB,1)...            <- missing part 
                            -Sz*solver.Khat(:,:,iCs)*Rcs2s(:,iCs)...
                            zeros(nB,1);];
                csTdV(:,(iCs-1)*3+1:(iCs)*3)=[    (solver.Khat(:,:,iCs)*NZero)'*Ihat*Vlin'... %<- influence due to deflection
                            -NZero'*solver.Khat(:,:,iCs)*Rcs2(:,iCs)*0 ...                             %<- influence on b due to speed
                            zeros(nB,1); ...
                            %
                            zeros(nB,1)...
                            (solver.Khat(:,:,iCs)*NZero)'*Ihat*Vlin'...
                            -NZero'*solver.Khat(:,:,iCs)*Rcs2(:,iCs)*0; ...
                            %
                            zeros(nB,1)...            <- missing part 
                            -Sx*solver.Khat(:,:,iCs)*Rcs2s(:,iCs)*0 ...
                            zeros(nB,1); ...
                            %
                            zeros(nB,1)...            <- missing part 
                            -Sy*solver.Khat(:,:,iCs)*Rcs2s(:,iCs)*0 ...
                            zeros(nB,1); ...
                            %
                            zeros(nB,1)...            <- missing part 
                            -Sz*solver.Khat(:,:,iCs)*Rcs2s(:,iCs)*0 ...
                            zeros(nB,1);];
            end
            %resort cST0 so that input is [dCs(1..nCs); dCsDot(1..nCs); dCsDdot(1..nCs);]
            csT0=[csT0(:,1:3:nCs*3) csT0(:,2:3:nCs*3) csT0(:,3:3:nCs*3)];
            %resort cSTdV so that input is [dCs(1..nCs); dCsDot(1..nCs); dCsDdot(1..nCs);]
            csTdV=[csTdV(:,1:3:nCs*3) csTdV(:,2:3:nCs*3) csTdV(:,3:3:nCs*3)];
            csT.d0=csT0;
            csT.dddV=csTdV;
            %old naming (cs by cs sorting, not all cs then all csdot then all csdotdot
%             csT.InputName=strcat(repmat({'dCs_'; 'dCsDot_';'dCsDotDot_'},nCs,1), cellstr(num2str(reshape(repmat(1:nCs,3,1),3*nCs,1),'%04d')));
            %new naming
            if ~isempty(csT.d0)
                csT.InputName=[     cellstr([repmat('dCs_',nCs,1) num2str([1:nCs]','%04d')]);...
                                    cellstr([repmat('dCsDot_',nCs,1) num2str([1:nCs]','%04d')]);...
                                    cellstr([repmat('dCsDotDot_',nCs,1) num2str([1:nCs]','%04d')]);];
                csT.OutputName=[   cellstr([repmat('b_',nB,1) num2str([1:nB]','%04d')]);...
                                    cellstr([repmat('bDot_',nB,1) num2str([1:nB]','%04d')]);...
                                    cellstr([repmat('vSInf_x_',nB,1) num2str([1:nB]','%04d')]);...
                                    cellstr([repmat('vSInf_y_',nB,1) num2str([1:nB]','%04d')]);...
                                    cellstr([repmat('vSInf_z_',nB,1) num2str([1:nB]','%04d')]);];

                csT.OutputGroup.b=[1:nB];
                csT.OutputGroup.bDot=[nB+1:2*nB];
                csT.OutputGroup.vSInf_x=[nB*2+1:3*nB];
                csT.OutputGroup.vSInf_y=[nB*3+1:4*nB];
                csT.OutputGroup.vSInf_z=[nB*4+1:5*nB];
                csT.InputGroup.dCs=[1:nCs];
                csT.InputGroup.dCsDot=[nCs+1:2*nCs];
                csT.InputGroup.dCsDotDot=[2*nCs+1:3*nCs];
            else
                csT.InputName=[];
                csT.OutputName=[];

                csT.OutputGroup.b=[1:nB];
                csT.OutputGroup.bDot=[nB+1:2*nB];
                csT.OutputGroup.vSInf_x=[nB*2+1:3*nB];
                csT.OutputGroup.vSInf_y=[nB*3+1:4*nB];
                csT.OutputGroup.vSInf_z=[nB*4+1:5*nB];
                csT.InputGroup.dCs=[1:nCs];
                csT.InputGroup.dCsDot=[nCs+1:2*nCs];
                csT.InputGroup.dCsDotDot=[2*nCs+1:3*nCs];
            end
        end
        function vgT=generateGustInputTransformation(solver)
            Vlin=solver.Uinf/norm(solver.Uinf);
            %number bound panels
            nB=size(solver.panels,2);            
            %number of gust Zones
            nGz=length(solver.gustZones);
            %Mapping matrix for gust zones
            M=zeros(nB,nGz);
            for iGz=1:nGz
                M(solver.gustZones{iGz},iGz)=1;
            end
            %when there is no gust zone specified, use all panels as one
            %gust Zone
            if isempty(M)
                nGz=1;
                M=ones(nB,1);
            end
            %normal vector blockdiagonal matrix
            normalvectorCell=num2cell(solver.colloc_nvec,1);
            NZero=blkdiag(normalvectorCell{:});
            
            %matrices selecting x,y and z components
            Sx=kron(eye(nB),[1 0 0]);
            Sy=kron(eye(nB),[0 1 0]);
            Sz=kron(eye(nB),[0 0 1]);
            %T matrix which transforms gust velocities or accelerations on
            %b or bDot
            T=NZero'*kron(M,eye(3));
            %Txyz matrix which transforms gust velocities on Vsinfxyz
            Tx=Sx*kron(M,eye(3));
            Ty=Sy*kron(M,eye(3));
            Tz=Sz*kron(M,eye(3));
            
            
            vgT.d=sparse([blkdiag(T,T);[[Tx;Ty;Tz] zeros(3*nB,3*nGz)]]);
            
            
            vgT.InputName=[ cellstr([repmat(['vG_x_'; 'vG_y_';'vG_z_'],nGz,1) num2str(reshape(repmat(1:nGz,3,1),3*nGz,1),'%04d')]);...
                            cellstr([repmat(['vGDot_x_'; 'vGDot_y_';'vGDot_z_'],nGz,1) num2str(reshape(repmat(1:nGz,3,1),3*nGz,1),'%04d')]);];
            
            vgT.OutputName=[   cellstr([repmat('b_',nB,1) num2str([1:nB]','%04d')]);...
                                cellstr([repmat('bDot_',nB,1) num2str([1:nB]','%04d')]);...
                                cellstr([repmat('vSInf_x_',nB,1) num2str([1:nB]','%04d')]);...
                                cellstr([repmat('vSInf_y_',nB,1) num2str([1:nB]','%04d')]);...
                                cellstr([repmat('vSInf_z_',nB,1) num2str([1:nB]','%04d')]);];
            vgT.OutputGroup.b=[1:nB];
            vgT.OutputGroup.bDot=[nB+1:2*nB];
            vgT.OutputGroup.vSInf_x=[nB*2+1:3*nB];
            vgT.OutputGroup.vSInf_y=[nB*3+1:4*nB];
            vgT.OutputGroup.vSInf_z=[nB*4+1:5*nB];
            vgT.InputGroup.vG=[1:3*nGz];
            vgT.InputGroup.vGDot=[3*nGz+1:6*nGz];
            %% calculate delays if there is multiple gust zones
            %delay is caluclated by the distance of 1/4 points of gust
            %zones in flight direction divided by the magnitude of the
            %flight speed
            if nGz>1
                % first calculate 1/4 points of gust zone
                quarterPoints=zeros(3,nGz);
                for iGz=1:nGz
                     quarterPoints(:,iGz)=sum(solver.colloc(:,solver.gustZones{iGz}).*repmat(solver.area(solver.gustZones{iGz}),3,1),2)./sum(solver.area(solver.gustZones{iGz}));
                end
                % then calculate distances to reference point (which is most
                % forward point of gust zone 1/4 points) and resulting
                % delay
                delay=Vlin*(quarterPoints(:,2:end)-quarterPoints(:,1));
                
                vgT.InputDelayFV=repmat(reshape(repmat([0 delay],3,1),nGz*3,1),2,1)*solver.Ma_corr;
                % merge gust inputs
%                 combIn=ss(blkdiag(repmat(eye(3),nGz,1),repmat(eye(3),nGz,1)));
%                 combIn.InputName={'vGx';'vGy';'vGz';'vGxDot';'vGyDot';'vGzDot'};    
%                 combIn.OutputName=vgT.InputName;
%                 combIn.InputGroup.vG=1:3;
%                 combIn.InputGroup.vGDot=4:6;
%                 vgT=series(combIn,vgT,'name');
            else
                 vgT.InputDelayFV=0;
            end
                
        end
        
        
    end
    
end

