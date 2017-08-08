%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%
classdef class_UVLM_LPVKernel
    %UNTITLED Summary of this class goes here
    %  here describe state equation
    % or store as function handle, maybe function handles can be used to
    % evaluate with properties?
    
    properties
        linSSM
        a
        b
        c
        d
        cF
        cF2
        dF
        dF2
        inputNames
        InputGroup
        outputNames
        OutputGroup
        stateNames
        %reduction matrix
        redT
    end
    
    methods
        function obj=class_UVLM_LPVKernel(solver) 
            disp('Generating LPV Kernel Matrices')
            %prepare for unit velocity
            Vref=1;
            
            nS=sum(solver.is_te);
            %number of chordwise panels
            nC=length(solver.Abb)/nS;
            %number of bound panels
            nB=size(solver.panels,2);
            %xW describes the spacing of the wake grid
            xW=solver.grid_wake(1,solver.panels_wake(1,1+2*nS))-solver.grid_wake(1,solver.panels_wake(1,nS+1));
            
            % number of ALL wake panels
            nW=size(solver.Abw,2);
            % number of REST wake panels
            nWR=nW-nS;
            
            %K1 is the influence coefficients of bound on bound panels
            K1=-solver.Abb;
            %K2 is the influence coefficients of small first wake row on bound
            K2=-solver.Abw(:,1:nS);
            %K3 is the influence coefficients of wake on bound
            K3=-solver.Abw(:,nS+1:end);
            % same for M1X M1Y M1Z M2X M2Y M2Z ...
            M1X=-solver.Abb_x;
            M2X=-solver.Abw_x(:,1:nS);
            M3X=-solver.Abw_x(:,nS+1:end);
            M1Y=-solver.Abb_y;
            M2Y=-solver.Abw_y(:,1:nS);
            M3Y=-solver.Abw_y(:,nS+1:end);
            M1Z=-solver.Abb_z;
            M2Z=-solver.Abw_z(:,1:nS);
            M3Z=-solver.Abw_z(:,nS+1:end);
            %K4 maps all body gamma to vector of gammas of last row
            K4=diag(solver.is_te)'*1;
            K4(all(K4==0,2),:)=[];
            %K5 is an eye matrix of the number of spanwise panels 
            K5=-eye(nS); 
            %K6 is transport within rest wake
%             K6=([zeros(nS,nWR) ;eye(nWR-nS) zeros(nWR-nS,nS)]-blkdiag(eye(nWR-nS),eye(nS)))*Vref/xW;
%             K6(end-nS+1:end,:)=K6(end-nS+1:end,:)*2701;
            dXVec=[ones(1,solver.settings.nFixedWakePanels) solver.settings.wakeGrowthFactor.^(1:solver.n_step-2-solver.settings.nFixedWakePanels) solver.settings.wakeGrowthFactor.^(solver.n_step+1-2-solver.settings.nFixedWakePanels)+solver.settings.addLengthFactorLastPanel/solver.settings.wakelength_factor*solver.n_step];
            dXVec2=(((dXVec(2:end)+dXVec(1:end-1))/2));
%             dXVec2(1:end-1)=1;
            K6=([zeros(nS,nWR) ;diag(repelem(1./dXVec2,nS)) zeros(nWR-nS,nS)]-blkdiag(eye(nS),diag(repelem(1./dXVec2,nS))))*Vref/xW;
%             K6=([zeros(nS,nWR) ;eye(nWR-nS) zeros(nWR-nS,nS)]-blkdiag(eye(nWR-nS),eye(nS)))*Vref/xW;
            %K7 is transport within first wake row
            K7=[eye(nS); zeros(nWR-nS,nS)]*Vref/xW;
            %% Assembly
            K8=+K6+K7*(K5-K4*K1^-1*K2)^-1*K4*K1^-1*K3;
            K9=K7*(K5-K4*K1^-1*K2)^-1*K4*K1^-1;

            L1=[zeros(1,nB); -diag(~solver.is_te(1:end-1)') zeros(nB-1,1)]+ eye(nB,nB);
            L2=1/2*([zeros(1,nB); diag(~solver.is_te(1:end-1)') zeros(nB-1,1)]+ eye(nB,nB));

            L3=(K2*K5^-1*K4-K1)^-1*K3;
            L4=(K2*K5^-1*K4-K1)^-1;
            L5=L3*K8;
            L6=L3*K9;

            L7=(M1X-M2X*K5^-1*K4)*L3+M3X;
            L8=(M1X-M2X*K5^-1*K4)*L4;
            L9=(M1Y-M2Y*K5^-1*K4)*L3+M3Y;
            L10=(M1Y-M2Y*K5^-1*K4)*L4;
            L11=(M1Z-M2Z*K5^-1*K4)*L3+M3Z;
            L12=(M1Z-M2Z*K5^-1*K4)*L4;
            %% Storage of Matrices
            obj.a=K8;
            obj.b=[K9 zeros(nWR,nB)];
            obj.c=[L1*L3; L2*L5; L7; L9; L11];
            obj.d=[L1*L4 zeros(nB,nB); L2*L6 L2*L4; L8 zeros(nB,nB); L10 zeros(nB,nB); L12 zeros(nB,nB)];
           
            nB=size(obj.c,1)/5;
            nStates=size(obj.a,1);
            %define lpv factors for C and D matrix
            obj.cF=[zeros(nB,nStates); ones(nB,nStates); zeros(3*nB,nStates)];
            obj.cF2=(obj.cF==0);
            obj.dF=[zeros(nB) zeros(nB); ones(nB) zeros(nB); zeros(3*nB,2*nB)];
            obj.dF2=(obj.dF==0);
            %% Storage of Names (Inputs, Outputs, States)
            obj.inputNames=[cellstr([repmat('b_',nB,1) num2str([1:nB]','%04d')]); cellstr([repmat('bDot_',nB,1) num2str([1:nB]','%04d')]);];
            obj.stateNames=[cellstr([repmat('wakeVort_',nWR,1) num2str([1:nWR]','%04d')]);];
            obj.outputNames=[   cellstr([repmat('gbEff_',nB,1) num2str([1:nB]','%04d')]);...
                                cellstr([repmat('gbEffDot_',nB,1) num2str([1:nB]','%04d')]);...
                                cellstr([repmat('vSInd_x_',nB,1) num2str([1:nB]','%04d')]);...
                                cellstr([repmat('vSInd_y_',nB,1) num2str([1:nB]','%04d')]);...
                                cellstr([repmat('vSInd_z_',nB,1) num2str([1:nB]','%04d')])];
            
            obj.OutputGroup.gbeff=[1:nB];
            obj.OutputGroup.gbDoteff=[nB+1:2*nB];
            obj.OutputGroup.vSInd_x=[2*nB+1:3*nB];
            obj.OutputGroup.vSInd_y=[3*nB+1:4*nB];
            obj.OutputGroup.vSInd_z=[4*nB+1:5*nB];
            
            obj.InputGroup.b=[1:nB];
            obj.InputGroup.bDot=[nB+1:2*nB];
        end
        function obj=linearize(obj,V)
            fprintf('linearizing kernel at V=%.3f m/s\n', V)
            
            obj.linSSM=ss(obj.a*V,obj.b*V,obj.c.*(obj.cF*V+obj.cF2),obj.d.*(obj.dF*V+obj.dF2));
            obj.linSSM.StateName=obj.stateNames;
            obj.linSSM.InputName=obj.inputNames;
            obj.linSSM.OutputName=obj.outputNames;
            obj.linSSM.InputGroup=obj.InputGroup;
            obj.linSSM.OutputGroup=obj.OutputGroup;
        end
        function obj=reduce(obj,nOrder,type)
            fprintf('reduce kernel to order %i using ', nOrder)
            fullSys=ss(obj.a,obj.b,obj.c,obj.d);
            if strcmp(type,'balred')
                disp('balanced reduction')
             [obj.a,obj.b,obj.c,obj.d]=ssdata(balred(fullSys,nOrder));
            end
            if strcmp(type,'bpod')
                disp('balanced proper orthogonal decomposition')
                VBpod=50;
                nPanels=size(obj.b,2)/2;
                fullSys=ss(obj.a*VBpod,obj.b*VBpod,obj.c.*(obj.cF*VBpod+obj.cF2),obj.d.*(obj.dF*VBpod+obj.dF2));
                %add integrators for boundary condition prior to
                %controllability estimation
                fullSysInt=series(ss(zeros(nPanels),eye(nPanels),[eye(nPanels);zeros(nPanels)],[zeros(nPanels);eye(nPanels)]),fullSys);
                nPODModes=50;
                stepsimtime=0:0.1:3; % only 6 Timesteps
                % empirical controllability
                inputs=find(any(fullSysInt.b));  %<- only use inputs which have an influence, not sure if that works
                [~, impT,impX]=impulse(fullSysInt(1,inputs),stepsimtime);  %<- improve by using sss toolbox?
                impTStep=impT(3)-impT(2);
                impXper=permute(impX,[2 1 3]);
                %remove integratorstatesdata)
                impXper=impXper(1:order(fullSys),:,:);
                xMat=reshape(impXper,order(fullSys),length(impT)*size(impX,3));
                xMat=xMat.*sqrt(impTStep);
                wcemp=xMat*xMat';
                % POD Modes

                 xMatPodAllt=reshape(impXper,order(fullSys),length(impT)*size(impX,3));
                 xMatPod=fullSys.c*xMatPodAllt(:,1:end);
                [Vr,~,~]=svd(xMatPod);
                % reduction base
                PodRed=Vr(:,1:nPODModes);      %add criteria for determination of the projection error by using the sum of eigenvalues (see paper)
                % PodRed=Vr*Vr';

                % empirical observability
                fullSys_adj=ss(fullSys.a',fullSys.c'*PodRed,fullSys.b',0);
                inputs=find(any(fullSys_adj.b));  %<- only use inputs which have an influence, not sure if this is correct
                [~, impAdjT,impAdjX]=impulse(fullSys_adj(1,inputs),stepsimtime);
                impAdjTStep=impAdjT(3)-impAdjT(2);
                impAdjXper=permute(impAdjX,[2 1 3]);
                yMatAllt=reshape(impAdjXper,order(fullSys),length(impAdjT)*size(impAdjX,3));
                yMat=yMatAllt(:,1:end);
                yMat=yMat.*sqrt(impAdjTStep);
                woemp=yMat*yMat';
                % build transformation matrix
                [redBaseF, ~]=eig(woemp*wcemp);
                %
                redBase=real(redBaseF(:,1:nOrder));     %<- workaround for complex eigenvalues -> also not sure if that works
                % actual reduction
                obj.redT=redBase;
                redSys=ss(pinv(redBase)*obj.a*redBase,  pinv(redBase)*obj.b, obj.c*redBase, obj.d);
                [obj.a,obj.b,obj.c,obj.d]=ssdata(redSys);
            end
            nB=size(obj.c,1)/5;
            obj.cF=[zeros(nB,nOrder); ones(nB,nOrder); zeros(3*nB,nOrder)];
            obj.cF2=(obj.cF==0);
            obj.dF=[zeros(nB) zeros(nB); ones(nB) zeros(nB); zeros(3*nB,2*nB)];
            obj.dF2=(obj.dF==0);
             obj.stateNames=[cellstr([repmat('aeroRedState_',nOrder,1) num2str([1:nOrder]','%04d')]);];
            disp('...done')
        end
    end
    
end

