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
            obj.rbmT=obj.generateRigidBodyMotionInputTransormation(solver);
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
            obj=obj.generateSectionalPressureOutputTransformation(solver, obj.Vlin);
            obj.csT=obj.generateCsInputTransormation(solver,obj.Vlin);
            % connect in a way that all outputs and inputs are preserved
            inputSys=parallel(obj.rbmT,obj.csT,'name');
            %prepare output matrix (from f to out)
            if obj.settings.sectPress
                allOutputs=[obj.totDerOutT; obj.sectPressOutT];
            else
                allOutputs=[obj.totDerOutT;];
            end
            %prepare from kernelOut to out
            if obj.settings.splitForce
                outputSys=series(obj.linOutT,allOutputs,[obj.linOutT.OutputGroup.Fp, obj.linOutT.OutputGroup.Fp_st, obj.linOutT.OutputGroup.Fp_ust],[allOutputs.InputGroup.Fp allOutputs.InputGroup.Fp_st allOutputs.InputGroup.Fp_ust]);
            else
                outputSys=series(obj.linOutT,allOutputs,obj.linOutT.OutputGroup.Fp,obj.totDerOutT.InputGroup.Fp);
            end
            
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
            gbdotF=gbrep*Areas.*NzeroHat*rho;
            %factor for gbf
            gbF=Scross*Ihat*V'*rho;
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
            T=[Ihat';Rxs';solver.Rcs]/(solver.state.rho_air/2*norm(V)^2*solver.reference.S_ref);
            if obj.settings.splitForce
                obj.totDerOutT=ss(blkdiag(T,T,T));
                obj.totDerOutT.InputName=[    cellstr([repmat(['Fp_x_'; 'Fp_y_';'Fp_z_'],nB,1) num2str(reshape(repmat(1:nB,3,1),3*nB,1),'%04d')]);...
                                                    cellstr([repmat(['Fp_st_x_'; 'Fp_st_y_';'Fp_st_z_'],nB,1) num2str(reshape(repmat(1:nB,3,1),3*nB,1),'%04d')]);...
                                                    cellstr([repmat(['Fp_ust_x_'; 'Fp_ust_y_';'Fp_ust_z_'],nB,1) num2str(reshape(repmat(1:nB,3,1),3*nB,1),'%04d')]);];
                nCs=size(solver.Rcs,1);
                obj.totDerOutT.OutputName=['CX'; 'CY'; 'CZ'; 'CL'; 'CM'; 'CN'; cellstr([repmat('Chinge_',nCs,1) num2str([1:nCs]','%04d')]);...
                                        'CX_st'; 'CY_st'; 'CZ_st'; 'CL_st'; 'CM_st'; 'CN_st'; cellstr([repmat('Chinge_st_',nCs,1) num2str([1:nCs]','%04d')]);...
                                        'CX_ust'; 'CY_ust'; 'CZ_ust'; 'CL_ust'; 'CM_ust'; 'CN_ust'; cellstr([repmat('Chinge_ust_',nCs,1) num2str([1:nCs]','%04d')]);];
                obj.totDerOutT.OutputGroup.Derivatives=[1:(6+nCs)];
                obj.totDerOutT.OutputGroup.Derivatives_st=[(6+nCs)+1:2*(6+nCs)];
                obj.totDerOutT.OutputGroup.Derivatives_ust=[2*(6+nCs)+1:3*(6+nCs)];
                obj.totDerOutT.InputGroup.Fp=[1:3*nB];
                obj.totDerOutT.InputGroup.Fp_st=[3*nB+1:6*nB];
                obj.totDerOutT.InputGroup.Fp_ust=[6*nB+1:9*nB];
            else
                obj.totDerOutT=ss(T);
                obj.totDerOutT.InputName=cellstr([repmat(['Fp_x_'; 'Fp_y_';'Fp_z_'],nB,1) num2str(reshape(repmat(1:nB,3,1),3*nB,1),'%04d')]);
                nCs=size(solver.Rcs,1);
                obj.totDerOutT.OutputName=['CX'; 'CY'; 'CZ'; 'CL'; 'CM'; 'CN'; cellstr([repmat('Chinge_',nCs,1) num2str([1:nCs]','%04d')])];
                obj.totDerOutT.OutputGroup.Derivatives=[1:(6+nCs)];
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
     
        function rbmT=generateRigidBodyMotionInputTransormation(solver)
            
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
            
            rbmT.InputName={'vB_x';'vB_y';'vB_z';'wB_x';'wB_y';'wB_z';'vBDot_x';'vBDot_y';'vBDot_z'; 'wBDot_x';'wBDot_y';'wBDot_z'};
            
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
        end
        function csT=generateCsInputTransormation(solver,Vlin)
            %notiz für freitag:
            % es fehlt der einfluss der
            % csdeflection auf die vsinf velocity - weshalb der drag
            % wahrscheinlich nicht stimmt bei csdef
            %nicht beides auf einmal, erst diese cst matrix machen dann mit
            %dem InputTrafoLin von UVLMData vergleichen.
            %außerdem ist diese Trafo ebenfalls abhängig von Vinf3D.... das
            %bedeutet ich muss es auch in linearize definieren...... damnit
            
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
                csT(:,(iCs-1)*3+1:(iCs)*3)=[    (solver.Khat(:,:,iCs)*NZero)'*Ihat*Vlin'...
                            -NZero'*solver.Khat(:,:,iCs)*Rcs2...
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
            csT=ss(csT);
            csT.InputName=strcat(repmat({'dCs_'; 'dCsDot_';'dCsDotDot_'},nCs,1), cellstr(num2str(reshape(repmat(1:nCs,3,1),3*nCs,1),'%04d')));
            csT.OutputName=[   cellstr([repmat('b_',nB,1) num2str([1:nB]','%04d')]);...
                                cellstr([repmat('bDot_',nB,1) num2str([1:nB]','%04d')]);...
                                cellstr([repmat('vSInf_x_',nB,1) num2str([1:nB]','%04d')]);...
                                cellstr([repmat('vSInf_y_',nB,1) num2str([1:nB]','%04d')]);...
                                cellstr([repmat('vSInf_z_',nB,1) num2str([1:nB]','%04d')]);];
           
            csT.OutputGroup.b=[1:nB];
            csT.OutputGroup.bDot=[nB+1:2*nB];
            rbmT.OutputGroup.vSInf_x=[nB*2+1:3*nB];
            rbmT.OutputGroup.vSInf_y=[nB*3+1:4*nB];
            rbmT.OutputGroup.vSInf_z=[nB*4+1:5*nB];
                    
                    
        end
        
        
    end
    
end

