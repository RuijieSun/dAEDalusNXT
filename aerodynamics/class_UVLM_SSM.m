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
            obj.vgT=obj.generateGustInputTransformation(solver);
            obj.LPVKernel=class_UVLM_LPVKernel(solver);
            if obj.settings.redOrder~=0
                obj.LPVKernel=obj.LPVKernel.reduce(obj.settings.redOrder, obj.settings.redMethod);
            end
                
            
        end
        
        function obj=linearize(obj,solver,V) 
            % prepare linearized parts
            obj.Vlin=V;
            %number bound panels
            nB=size(solver.panels,2);
            
            obj.LPVKernel=obj.LPVKernel.linearize(norm(obj.Vlin));
            obj=obj.generateLinearOutputTransformation(solver, obj.Vlin);
            obj=obj.generateDerivativeOutputTransformation(solver, obj.Vlin);
            if obj.settings.sectPress
                obj=obj.generateSectionalPressureOutputTransformation(solver, obj.Vlin);
            end
            obj.csT=obj.generateCsInputTransformation(solver,obj.Vlin);
            % parallel connection of various inputTransformations
            if ~isempty(obj.csT)
                inputSys=parallel(obj.rbmT,obj.csT,'name');
            else
                inputSys=obj.rbmT;
            end
            %add panel velocities as inputs
            inputSys=parallel(inputSys,obj.vpT,'name');
            %add gust Zone velocities as inputs
            inputSys=parallel(inputSys,obj.vgT,'name');
            bT=ss(eye(2*nB));
            bT.inputName=obj.LPVKernel.linSSM.InputName([obj.LPVKernel.linSSM.InputGroup.b obj.LPVKernel.linSSM.InputGroup.bDot]);
            bT.InputGroup.b=1:nB;
            bT.InputGroup.bDot=nB+1:2*nB;
            bT.OutputGroup=bT.InputGroup;
            bT.OutputName=bT.InputName;
            %add normal velocities as inputs
            inputSys=parallel(inputSys,bT,'name');
            
            %prepare output matrix (from f to out)
            if obj.settings.sectPress
                allOutputs=[obj.totDerOutT; obj.sectPressOutT];
            else
                allOutputs=[obj.totDerOutT;];
            end
            %prepare outputsys(from kernelOut to out)
            if obj.settings.splitForce
                outputSys=series(obj.linOutT,allOutputs,[obj.linOutT.OutputGroup.Fp, obj.linOutT.OutputGroup.Fp_st, obj.linOutT.OutputGroup.Fp_ust],[allOutputs.InputGroup.Fp allOutputs.InputGroup.Fp_st allOutputs.InputGroup.Fp_ust]);
                
            else
                outputSys=series(obj.linOutT,allOutputs,obj.linOutT.OutputGroup.Fp,obj.totDerOutT.InputGroup.Fp);
            end
            %add obj.linOutT to have forces as outputs
            outputSys=parallel(outputSys,obj.linOutT,'name');
            
            %connect input, kernel, output
             obj.linSSM=series(...
                 series(...
                 inputSys,obj.LPVKernel.linSSM,...
                 [inputSys.OutputGroup.b, inputSys.OutputGroup.bDot],...
                 [obj.LPVKernel.linSSM.InputGroup.b, obj.LPVKernel.linSSM.InputGroup.bDot])...
                    ,outputSys,...
                    [obj.LPVKernel.linSSM.OutputGroup.gbeff, obj.LPVKernel.linSSM.OutputGroup.gbDoteff],...
                    [outputSys.InputGroup.gbeff, outputSys.InputGroup.gbDoteff]);
        end
        function obj=generateLinearOutputTransformation(obj, solver,V)
            %matrix which transforms gb and gbeff of LPV kernel to Forces,
            %stored as sys
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
            gbdotF=gbrep*Areas.*NzeroHat*rho*1/betaInf;
            %factor for gbf
            gbF=Scross*Ihat*V'*rho*1/betaInf;
            %overall transformation matrix
            if obj.settings.splitForce
                obj.linOutT=ss([diag(gbF)*gbrep -diag(gbdotF)*gbrep;...
                            diag(gbF)*gbrep zeros(nB*3,nB*3)*gbrep;...
                            zeros(nB*3,nB*3)*gbrep -diag(gbdotF)*gbrep]);
                obj.linOutT.InputName=[     cellstr([repmat('gbEff_',nB,1) num2str([1:nB]','%04d')]);...
                                                    cellstr([repmat('gbEffDot_',nB,1) num2str([1:nB;]','%04d')]);];
                obj.linOutT.OutputName=[    cellstr([repmat(['Fp_x_'; 'Fp_y_';'Fp_z_'],nB,1) num2str(reshape(repmat(1:nB,3,1),3*nB,1),'%04d')]);...
                                                    cellstr([repmat(['Fp_st_x_'; 'Fp_st_y_';'Fp_st_z_'],nB,1) num2str(reshape(repmat(1:nB,3,1),3*nB,1),'%04d')]);...
                                                    cellstr([repmat(['Fp_ust_x_'; 'Fp_ust_y_';'Fp_ust_z_'],nB,1) num2str(reshape(repmat(1:nB,3,1),3*nB,1),'%04d')]);];
                obj.linOutT.OutputGroup.Fp=[1:3*nB];
                obj.linOutT.OutputGroup.Fp_st=[3*nB+1:6*nB];
                obj.linOutT.OutputGroup.Fp_ust=[6*nB+1:9*nB];
                obj.linOutT.InputGroup.gbeff=[1:nB];
                obj.linOutT.InputGroup.gbDoteff=[nB+1:2*nB];
            else
                
                obj.linOutT=ss([diag(gbF)*gbrep -diag(gbdotF)*gbrep]);
                obj.linOutT.InputName=[     cellstr([repmat('gbEff_',nB,1) num2str([1:nB]','%04d')]);...
                                                    cellstr([repmat('gbEffDot_',nB,1) num2str([1:nB;]','%04d')]);];
                obj.linOutT.OutputName=[    cellstr([repmat(['Fp_x_'; 'Fp_y_';'Fp_z_'],nB,1) num2str(reshape(repmat(1:nB,3,1),3*nB,1),'%04d')]);];
                obj.linOutT.OutputGroup.Fp=[1:3*nB];
                obj.linOutT.InputGroup.gbeff=[1:nB];
                obj.linOutT.InputGroup.gbDoteff=[nB+1:2*nB];
            end
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
    end
    methods(Static) % static
     
        function rbmT=generateRigidBodyMotionInputTransformation(solver)
            
            %number bound panels
            nB=size(solver.panels,2);            
            %matrix to repeat 3x1 vector npanel times
            Ihat=repmat(eye(3),nB,1);
            %normal vector blockdiagonal matrix
            normalvectorCell=num2cell(solver.colloc_nvec,1);
            NZero=blkdiag(normalvectorCell{:});
            %matrix containing the cross product matrices of the vectors
            %from the reference point to the collocation points
            rDwnPCellSkew=skewsymmatrixCell(num2cell(solver.r_dwn',2));
            Rx=[rDwnPCellSkew{:}]';
            %matrix containing the cross product matrices of the vectors
            %from the reference point to the FVAP points
            rSPCellSkew=skewsymmatrixCell(num2cell(solver.r',2)); 
            Rxs=[rSPCellSkew{:}]';
            %matrices selecting x,y and z components
            Sx=kron(eye(nB),[1 0 0]);
            Sy=kron(eye(nB),[0 1 0]);
            Sz=kron(eye(nB),[0 0 1]);
            
            
            rbmT=ss([NZero'*Ihat ...
                    zeros(nB,3)...
                    -NZero'*Rx...
                    zeros(nB,3);
                    %
                    zeros(nB,3)...
                    NZero'*Ihat ...
                    zeros(nB,3)...
                    -NZero'*Rx;...
                    
                    Sx*Ihat ...
                    zeros(nB,3)...
                    -Sx*Rxs...
                    zeros(nB,3)...
                    %
                    Sy*Ihat ...
                    zeros(nB,3)...
                    -Sy*Rxs...
                    zeros(nB,3);...
                    %
                    Sz*Ihat ...
                    zeros(nB,3)...
                    -Sz*Rxs...
                    zeros(nB,3);...
                ]);
            
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
            
            
            vpT=ss([blkdiag(NZero',NZero');[[Sx ;Sy;Sz] zeros(nB*3,nB*3)]]);
            
            
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
        function csT=generateCsInputTransformation(solver,Vlin)
                       
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
            Rcs2=sum(reshape(solver.r_cs,nB*3,1,nCs),3);
            %same for segment center points
            Rcs2s=sum(reshape(solver.r_css,nB*3,1,nCs),3);
            
            %matrices selecting x,y and z components
            Sx=kron(eye(nB),[1 0 0]);
            Sy=kron(eye(nB),[0 1 0]);
            Sz=kron(eye(nB),[0 0 1]);
            %prealloc
            csT=zeros(nB*5,nCs*3);
            for iCs=1:nCs
                csT(:,(iCs-1)*3+1:(iCs)*3)=[    (solver.Khat(:,:,iCs)*NZero)'*Ihat*Vlin'... %<- influence due to deflection
                            -NZero'*solver.Khat(:,:,iCs)*Rcs2...                             %<- influence on b due to speed
                            zeros(nB,1); ...
                            %
                            zeros(nB,1)...
                            (solver.Khat(:,:,iCs)*NZero)'*Ihat*Vlin'...
                            -NZero'*solver.Khat(:,:,iCs)*Rcs2; ...
                            %
                            zeros(nB,1)...            <- missing part 
                            -Sx*solver.Khat(:,:,iCs)*Rcs2s...
                            zeros(nB,1); ...
                            %
                            zeros(nB,1)...            <- missing part 
                            -Sy*solver.Khat(:,:,iCs)*Rcs2s...
                            zeros(nB,1); ...
                            %
                            zeros(nB,1)...            <- missing part 
                            -Sz*solver.Khat(:,:,iCs)*Rcs2s...
                            zeros(nB,1);];
            end
            %resort cST so that input is [dCs(1..nCs); dCsDot(1..nCs); dCsDdot(1..nCs);]
            csT=[csT(:,1:3:nCs*3) csT(:,2:3:nCs*3) csT(:,3:3:nCs*3)];
            csT=ss(csT);
            %old naming (cs by cs sorting, not all cs then all csdot then all csdotdot
%             csT.InputName=strcat(repmat({'dCs_'; 'dCsDot_';'dCsDotDot_'},nCs,1), cellstr(num2str(reshape(repmat(1:nCs,3,1),3*nCs,1),'%04d')));
            %new naming
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
                    
                    
        end
        function vgT=generateGustInputTransformation(solver)
            
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
            
            
            vgT=ss([blkdiag(T,T);[[Tx;Ty;Tz] zeros(3*nB,3*nGz)]]);
            
            
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
        end
        
        
    end
    
end

