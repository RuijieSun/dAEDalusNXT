classdef OptimizationCase
    %OPTIMIZATIONCASE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        %case settings
        settings@OptimizationCaseSettings;
        %aircraft 
        ac@class_aircraft;        
        %structure 
        str@class_beam_collection;
        %aeSSM
        aeSSM;
        %array of manevuer cases
        maneuverCases
        %array of gustCases
        gustCases
        % Struct containing load hulls, nodal displacements and info
        loadHullData
        %struct containing all designVariables
        desVar
        % types of design variables as vector 
        desVarTypes
        % scales of all desVar
        desVarScale
        %struct containing all constraints
        constraints
        %objective
        objective
        %lowerBounds for desVar
        lowerBounds
        %upperBounds for desVar
        upperBounds
        %VLM Solvers for all Maneuver and Gust Cases
        maneuverVLM@class_VLM_solver
        %VLM Solvers for all Maneuver and Gust Cases
        gustVLM@class_VLM_solver
        % finite difference type
        fdType
        % obj gradient
        objGrad
        % const gradient
        constGrad
        % vector of stepsizes used for gradient computation
        gradStepVec
        % history of optimization run
        history
        end
    
    methods
        function obj=OptimizationCase(aircraft, aircraft_structure, aeSSM, maneuverCases, gustCases, settings)
            %% initialize 
            obj.settings=settings; 
            obj.ac=aircraft;
            obj.str=aircraft_structure;
            obj.aeSSM=aeSSM;
            obj.maneuverCases=maneuverCases;
            obj.gustCases=gustCases;
            %% prepare lighter aessm
            obj.aeSSM.uvlm.Abb=[];
            obj.aeSSM.uvlm.Abw=[];
            obj.aeSSM.uvlm.Aww=[];
            obj.aeSSM.uvlm.Cbb=[];
            obj.aeSSM.uvlm.Awb=[];
            obj.aeSSM.uvlm.Abb_x=[];
            obj.aeSSM.uvlm.Abb_y=[];
            obj.aeSSM.uvlm.Abb_z=[];
            obj.aeSSM.uvlm.Abw_x=[];
            obj.aeSSM.uvlm.Abw_y=[];
            obj.aeSSM.uvlm.Abw_z=[];
            obj.aeSSM.uvlm.ssm.LPVKernel.a0=[];
            obj.aeSSM.uvlm.ssm.LPVKernel.b0=[];
            obj.aeSSM.uvlm.ssm.LPVKernel.c0=[];
            obj.aeSSM.uvlm.ssm.LPVKernel.d0=[];
            obj.aeSSM.uvlm.ssm.LPVKernel.dadV=[];
            obj.aeSSM.uvlm.ssm.LPVKernel.dbdV=[];
            obj.aeSSM.uvlm.ssm.LPVKernel.dcdV=[];
            obj.aeSSM.uvlm.ssm.LPVKernel.dddV=[];
            obj.aeSSM.uvlm.Khat=[];
            obj.aeSSM.uvlm.Rcs=[];
            obj.aeSSM.uvlm.Rcss=[];
            obj.aeSSM.uvlm.grid_wake=[];
            obj.aeSSM.uvlm.grid_wake_old=[];
            obj.aeSSM.uvlm.panels_wake=[];
            obj.aeSSM.uvlm.ssm.redBase=[];
            obj.aeSSM.uvlm.ssm.redBaseInv=[];
            obj.aeSSM.redBase=[];
            obj.aeSSM.redBaseInv=[];
            %% get init str design variables and bounds (todo: add stiffness direction as design variables)
            for iBeam=1:length(obj.settings.beamsToOptimize)
                
                % check symmetry of beam
                if obj.str.beam(obj.settings.beamsToOptimize(iBeam)).is_sym 
                    offset=obj.str.beam(obj.settings.beamsToOptimize(iBeam)).nel/2;
                else
                    offset=0;
                end
                % calculate thicknessScaleFactor for Skin and Spar based on
                % respective Area; 
                lengths=[obj.str.beam(iBeam).beamelement(obj.settings.thicknessControlStations{iBeam}+offset).le];
                cs=[obj.str.beam(iBeam).beamelement(obj.settings.thicknessControlStations{iBeam}+offset).crosssection];
                ASkin=[cs.w].*lengths;
                ASpar=[cs.h].*lengths;
                scaleFactorSkin=ASkin/mean(ASkin);
                scaleFactorSpar=ASpar/mean(ASpar);
                obj.settings.thicknessScaleFactorSkin{iBeam}=scaleFactorSkin*obj.settings.thicknessScaleFactor;
                obj.settings.thicknessScaleFactorSpar{iBeam}=scaleFactorSpar*mean(ASpar)/mean(ASkin)*obj.settings.thicknessScaleFactor;
                
                
                % init str desVar
                if obj.str.beam(obj.settings.beamsToOptimize(iBeam)).anisotropic
                    
                    obj.desVar.Str(iBeam).tSkLo=zeros(1,length(obj.settings.thicknessControlStations{iBeam}));
                    obj.desVar.Str(iBeam).tSkUp=zeros(1,length(obj.settings.thicknessControlStations{iBeam}));
                    obj.desVar.Str(iBeam).tSpFr=zeros(1,length(obj.settings.thicknessControlStations{iBeam}));
                    obj.desVar.Str(iBeam).tSpRe=zeros(1,length(obj.settings.thicknessControlStations{iBeam}));
                    obj.upperBounds.Str(length(obj.settings.beamsToOptimize)).tSkLo=[];
                    obj.upperBounds.Str(length(obj.settings.beamsToOptimize)).tSkUp=[];
                    obj.lowerBounds.Str(length(obj.settings.beamsToOptimize)).tSpFr=[];
                    obj.lowerBounds.Str(length(obj.settings.beamsToOptimize)).tSpRe=[];
                    
                   
                    for iSt=1:length(obj.settings.thicknessControlStations{iBeam})
                        % get current thicknesses as starting value
                        obj.desVar.Str(iBeam).tSkLo(iSt)=obj.str.beam(obj.settings.beamsToOptimize(iBeam)).beamelement(obj.settings.thicknessControlStations{iBeam}(iSt)+offset).crosssection.t_sk_lo;
                        obj.desVar.Str(iBeam).tSkUp(iSt)=obj.str.beam(obj.settings.beamsToOptimize(iBeam)).beamelement(obj.settings.thicknessControlStations{iBeam}(iSt)+offset).crosssection.t_sk_up;
                        obj.desVar.Str(iBeam).tSpFr(iSt)=obj.str.beam(obj.settings.beamsToOptimize(iBeam)).beamelement(obj.settings.thicknessControlStations{iBeam}(iSt)+offset).crosssection.t_sp_fr;
                        obj.desVar.Str(iBeam).tSpRe(iSt)=obj.str.beam(obj.settings.beamsToOptimize(iBeam)).beamelement(obj.settings.thicknessControlStations{iBeam}(iSt)+offset).crosssection.t_sp_re;
                        %set upper bounds to inf
                        obj.upperBounds.Str(iBeam).tSkLo(iSt)=0.1;
                        obj.upperBounds.Str(iBeam).tSkUp(iSt)=0.1;
                        obj.upperBounds.Str(iBeam).tSpFr(iSt)=0.1;
                        obj.upperBounds.Str(iBeam).tSpRe(iSt)=0.1;
                        % set lower bounds to min thickness
                        obj.lowerBounds.Str(iBeam).tSkLo(iSt)=obj.str.beam(obj.settings.beamsToOptimize(iBeam)).beamelement(obj.settings.thicknessControlStations{iBeam}(iSt)+offset).crosssection.t_min_sk;
                        obj.lowerBounds.Str(iBeam).tSkUp(iSt)=obj.str.beam(obj.settings.beamsToOptimize(iBeam)).beamelement(obj.settings.thicknessControlStations{iBeam}(iSt)+offset).crosssection.t_min_sk;
                        obj.lowerBounds.Str(iBeam).tSpFr(iSt)=obj.str.beam(obj.settings.beamsToOptimize(iBeam)).beamelement(obj.settings.thicknessControlStations{iBeam}(iSt)+offset).crosssection.t_min_sp;
                        obj.lowerBounds.Str(iBeam).tSpRe(iSt)=obj.str.beam(obj.settings.beamsToOptimize(iBeam)).beamelement(obj.settings.thicknessControlStations{iBeam}(iSt)+offset).crosssection.t_min_sp;
                    end
                    
                    if ~isempty(obj.settings.stiffnessDirectionControlStations)
                        if ~isempty(obj.settings.stiffnessDirectionControlStations{iBeam})
                            obj.desVar.Str(iBeam).phiSkLo=zeros(1,length(obj.settings.stiffnessDirectionControlStations{iBeam}));
                            obj.desVar.Str(iBeam).phiSkUp=zeros(1,length(obj.settings.stiffnessDirectionControlStations{iBeam}));
                            obj.lowerBounds.Str(iBeam).phiSkLo=obj.settings.stiffnessDirectionBounds(1)*ones(1,length(obj.settings.stiffnessDirectionControlStations{iBeam}));
                            obj.lowerBounds.Str(iBeam).phiSkUp=obj.settings.stiffnessDirectionBounds(1)*ones(1,length(obj.settings.stiffnessDirectionControlStations{iBeam}));
                            obj.upperBounds.Str(iBeam).phiSkLo=obj.settings.stiffnessDirectionBounds(2)*ones(1,length(obj.settings.stiffnessDirectionControlStations{iBeam}));
                            obj.upperBounds.Str(iBeam).phiSkUp=obj.settings.stiffnessDirectionBounds(2)*ones(1,length(obj.settings.stiffnessDirectionControlStations{iBeam}));
                            if obj.settings.sparTailoring
                                obj.desVar.Str(iBeam).phiSpFr=zeros(1,length(obj.settings.stiffnessDirectionControlStations{iBeam}));
                                obj.desVar.Str(iBeam).phiSpRe=zeros(1,length(obj.settings.stiffnessDirectionControlStations{iBeam}));
                                obj.lowerBounds.Str(iBeam).phiSpFr=obj.settings.stiffnessDirectionBounds(1)*ones(1,length(obj.settings.stiffnessDirectionControlStations{iBeam}));
                                obj.lowerBounds.Str(iBeam).phiSpRe=obj.settings.stiffnessDirectionBounds(1)*ones(1,length(obj.settings.stiffnessDirectionControlStations{iBeam}));
                                obj.upperBounds.Str(iBeam).phiSpFr=obj.settings.stiffnessDirectionBounds(2)*ones(1,length(obj.settings.stiffnessDirectionControlStations{iBeam}));
                                obj.upperBounds.Str(iBeam).phiSpRe=obj.settings.stiffnessDirectionBounds(2)*ones(1,length(obj.settings.stiffnessDirectionControlStations{iBeam}));
                            end
                        end
                    end
                    if ~isempty(obj.settings.stiffnessDirectionPolynomialOrder)
                        if ~isempty(obj.settings.stiffnessDirectionPolynomialOrder{iBeam})
                            obj.desVar.Str(iBeam).phiPolSkLo=zeros(1,obj.settings.stiffnessDirectionPolynomialOrder{iBeam});
                            obj.desVar.Str(iBeam).phiPolSkUp=zeros(1,obj.settings.stiffnessDirectionPolynomialOrder{iBeam});
                            obj.lowerBounds.Str(iBeam).phiPolSkLo=obj.settings.stiffnessDirectionPolBounds(1)*ones(1,obj.settings.stiffnessDirectionPolynomialOrder{iBeam});
                            obj.lowerBounds.Str(iBeam).phiPolSkUp=obj.settings.stiffnessDirectionPolBounds(1)*ones(1,obj.settings.stiffnessDirectionPolynomialOrder{iBeam});
                            obj.upperBounds.Str(iBeam).phiPolSkLo=obj.settings.stiffnessDirectionPolBounds(2)*ones(1,obj.settings.stiffnessDirectionPolynomialOrder{iBeam});
                            obj.upperBounds.Str(iBeam).phiPolSkUp=obj.settings.stiffnessDirectionPolBounds(2)*ones(1,obj.settings.stiffnessDirectionPolynomialOrder{iBeam});
                            if obj.settings.sparTailoring
                                obj.desVar.Str(iBeam).phiPolSpFr=zeros(1,obj.settings.stiffnessDirectionPolynomialOrder{iBeam});
                                obj.desVar.Str(iBeam).phiPolSpRe=zeros(1,obj.settings.stiffnessDirectionPolynomialOrder{iBeam});
                                obj.lowerBounds.Str(iBeam).phiPolSpFr=obj.settings.stiffnessDirectionPolBounds(1)*ones(1,obj.settings.stiffnessDirectionPolynomialOrder{iBeam});
                                obj.lowerBounds.Str(iBeam).phiPolSpRe=obj.settings.stiffnessDirectionPolBounds(1)*ones(1,obj.settings.stiffnessDirectionPolynomialOrder{iBeam});
                                obj.upperBounds.Str(iBeam).phiPolSpFr=obj.settings.stiffnessDirectionPolBounds(2)*ones(1,obj.settings.stiffnessDirectionPolynomialOrder{iBeam});
                                obj.upperBounds.Str(iBeam).phiPolSpRe=obj.settings.stiffnessDirectionPolBounds(2)*ones(1,obj.settings.stiffnessDirectionPolynomialOrder{iBeam});
                            end
                        end
                    end
                else
                    % isotropic model accepts only same thickness for
                    % top/bottom skin and front/rear spar
                    obj.desVar.Str(length(obj.settings.beamsToOptimize)).tSk=[];
                    obj.desVar.Str(length(obj.settings.beamsToOptimize)).tSp=[];
                    obj.upperBounds.Str(length(obj.settings.beamsToOptimize)).tSk=[];
                    obj.lowerBounds.Str(length(obj.settings.beamsToOptimize)).tSp=[];
                    obj.desVar.Str(iBeam).tSk=zeros(1,length(obj.settings.thicknessControlStations{iBeam}));
                    obj.desVar.Str(iBeam).tSp=zeros(1,length(obj.settings.thicknessControlStations{iBeam}));
                    % get current thicknesses as starting value
                    for iSt=1:length(obj.settings.thicknessControlStations{iBeam})
                        obj.desVar.Str(iBeam).tSk(iSt)=obj.str.beam(obj.settings.beamsToOptimize(iBeam)).beamelement(obj.settings.thicknessControlStations{iBeam}(iSt)+offset).crosssection.t_sk_lo;
                        obj.desVar.Str(iBeam).tSp(iSt)=obj.str.beam(obj.settings.beamsToOptimize(iBeam)).beamelement(obj.settings.thicknessControlStations{iBeam}(iSt)+offset).crosssection.t_sp_fr;
                        obj.upperBounds.Str(iBeam).tSk(iSt)=Inf;
                        obj.upperBounds.Str(iBeam).tSp(iSt)=Inf;
                        obj.lowerBounds.Str(iBeam).tSk(iSt)=obj.str.beam(obj.settings.beamsToOptimize(iBeam)).beamelement(obj.settings.thicknessControlStations{iBeam}(iSt)+offset).crosssection.t_min_sk;
                        obj.lowerBounds.Str(iBeam).tSp(iSt)=obj.str.beam(obj.settings.beamsToOptimize(iBeam)).beamelement(obj.settings.thicknessControlStations{iBeam}(iSt)+offset).crosssection.t_min_sp;
                    end
                end
            end

            
            %% initialize VLM Solvers for every ManeuverCase
            %get reference from AC
            ref=obj.ac.reference;
            for iCase=1:length(obj.maneuverCases)
                %set target CG to ref
                ref.p_ref=obj.maneuverCases(iCase).aircraft_state.CG_ref;
                % init VLM solver
                obj.maneuverVLM(iCase)=class_VLM_solver(obj.ac.grid,obj.ac.te_idx,obj.ac.panels,obj.maneuverCases(iCase).aerodynamic_state,ref, obj.ac.control_surfaces_parents);
                %initial solve for AIC matrix
                obj.maneuverVLM(iCase)=obj.maneuverVLM(iCase).f_solve_std();
            end
            %% initialize VLM Solvers for every gustCase
            for iCase=1:length(obj.gustCases)
                obj.gustVLM(iCase)=class_VLM_solver(obj.ac.grid,obj.ac.te_idx,obj.ac.panels,obj.gustCases(iCase).aerodynamic_state,obj.ac.reference, obj.ac.control_surfaces_parents);
                %initial solve for AIC matrix
                obj.gustVLM(iCase)=obj.gustVLM(iCase).f_solve_std();
            end
            
            %% initialize MLA design variables and bounds
            %per loadcase an MLA desVar for every control surface allowed
            %for MLA is added
            if ~isempty(obj.maneuverCases)
                obj.desVar.maneuverCases(length(obj.maneuverCases)).mla=[];
                obj.lowerBounds.maneuverCases(length(obj.maneuverCases)).mla=[];
                obj.upperBounds.maneuverCases(length(obj.maneuverCases)).mla=[];
                obj.desVar.maneuverCases(length(obj.maneuverCases)).rollAlloc=[];
                for iCase=1:length(obj.maneuverCases)
                    if ~isempty(obj.maneuverCases(iCase).mlaSurf) % when mla is present
                        if isempty(obj.maneuverCases(iCase).mlaVal) % when mla is optimized and not preset
                            obj.desVar.maneuverCases(iCase).mla=zeros(1,length(obj.maneuverCases(iCase).mlaSurf));
                            
                            if isempty(obj.maneuverCases(iCase).rollSurf)
                                obj.lowerBounds.maneuverCases(iCase).mla=obj.maneuverCases(iCase).mlaLB;
                                obj.upperBounds.maneuverCases(iCase).mla=obj.maneuverCases(iCase).mlaUB;
                            else
                                obj.lowerBounds.maneuverCases(iCase).mla=-Inf*ones(size(obj.maneuverCases(iCase).mlaSurf));
                                obj.upperBounds.maneuverCases(iCase).mla=Inf*ones(size(obj.maneuverCases(iCase).mlaSurf));
                            end
                            if ~isempty(obj.maneuverCases(iCase).rollSurf) % if this is a roll case add roll desvar and bounds
                                rollAlloc=zeros(1,length(obj.maneuverCases(iCase).rollSurf));
                                %init point: only ail is used for roll control
                                ailId=contains(maneuverCases(iCase).rollSurf,'ail');
                                rollAlloc(ailId)=1;
                                obj.desVar.maneuverCases(iCase).rollAlloc=rollAlloc;
                                obj.lowerBounds.maneuverCases(iCase).rollAlloc=-Inf*ones(size(rollAlloc));
                                obj.upperBounds.maneuverCases(iCase).rollAlloc=Inf*ones(size(rollAlloc));
                            end
                                
                        end
                    end
                end
                
                %% initialize constraints
                % for every maneuver loadcase stress- and strain
                % constraints for every beam which is optimized are added
                obj.constraints.maneuverCases(length(obj.maneuverCases)).strength=[];
                for iCase=1:length(obj.maneuverCases)
                    for iBeam=1:length(obj.settings.beamsToOptimize)
                        obj.constraints.maneuverCases(iCase).strength(iBeam).spFr=zeros(1,obj.str.beam(obj.settings.beamsToOptimize(iBeam)).nel);
                        obj.constraints.maneuverCases(iCase).strength(iBeam).spRe=zeros(1,obj.str.beam(obj.settings.beamsToOptimize(iBeam)).nel);
                        obj.constraints.maneuverCases(iCase).strength(iBeam).skLo=zeros(1,obj.str.beam(obj.settings.beamsToOptimize(iBeam)).nel);
                        obj.constraints.maneuverCases(iCase).strength(iBeam).skUp=zeros(1,obj.str.beam(obj.settings.beamsToOptimize(iBeam)).nel);
                        if ~isempty(obj.str.beam(obj.settings.beamsToOptimize(iBeam)).ribLocations)
                            obj.constraints.maneuverCases(iCase).buckl(iBeam).spFr=zeros(1,obj.str.beam(obj.settings.beamsToOptimize(iBeam)).nel);
                            obj.constraints.maneuverCases(iCase).buckl(iBeam).spRe=zeros(1,obj.str.beam(obj.settings.beamsToOptimize(iBeam)).nel);
                            obj.constraints.maneuverCases(iCase).buckl(iBeam).skLo=zeros(1,obj.str.beam(obj.settings.beamsToOptimize(iBeam)).nel);
                            obj.constraints.maneuverCases(iCase).buckl(iBeam).skUp=zeros(1,obj.str.beam(obj.settings.beamsToOptimize(iBeam)).nel);
                        end
                            
                    end
                    if and(~isempty(obj.maneuverCases(iCase).rollSurf), isempty(maneuverCases(iCase).rollAllocVal)) %if rollcase, mla cs bounds will be handled by constraints as cs has multiple functions
                        if and(~isempty(obj.maneuverCases(iCase).mlaSurf), isempty(obj.maneuverCases(iCase).mlaVal)) % if mla is present and not preset
                            if length(obj.maneuverCases(iCase).mlaSurf)~=length(obj.maneuverCases(iCase).rollSurf)
                                disp('warning, vector of mla surfaces different from roll surfaces, not implemented yet')
                            end
                        end
                        obj.constraints.maneuverCases(iCase).csDefUB=zeros(1, 2*length(obj.maneuverCases(iCase).rollSurf));
                        obj.constraints.maneuverCases(iCase).csDefLB=zeros(1, 2*length(obj.maneuverCases(iCase).rollSurf));
                    end
                end
            end
           
            %% initialize Gust Cases design variables and bounds
            %per loadcase an Gla desVar is created for the filter coeffs.
            
            % if gustcases are present, add gust ctr desvar, constraints for buckling and
            % fail
            if ~isempty(obj.gustCases)
                 %initialize MLA design variables and bounds for gust cases 
                obj.desVar.gustCases(length(obj.gustCases)).mla=[];
                obj.lowerBounds.gustCases(length(obj.gustCases)).mla=[];
                obj.upperBounds.gustCases(length(obj.gustCases)).mla=[];
                for iCase=1:length(obj.gustCases)
                    if ~isempty(obj.gustCases(iCase).mlaSurf)
                        if isempty(obj.gustCases(iCase).mlaVal)
                            obj.desVar.gustCases(iCase).mla=zeros(1,length(obj.gustCases(iCase).mlaSurf));
                            obj.lowerBounds.gustCases(iCase).mla=obj.gustCases(iCase).mlaLB;
                            obj.upperBounds.gustCases(iCase).mla=obj.gustCases(iCase).mlaUB;
                        end
                    end
                end
          
                %initialize GLA design variables and bounds for gust cases 
                
                obj.desVar.gustCases(length(obj.gustCases)).gla=[];
                glaFlag=~cellfun(@isempty,{obj.gustCases.glaSurf});
                glaPresetFlag=~cellfun(@isempty,{obj.gustCases.glaVal});
                %check if there is a case where gla is a design variable
                if any(and(glaFlag,~glaPresetFlag))
                    obj.lowerBounds.gustCases(length(obj.gustCases)).gla=[];
                    obj.upperBounds.gustCases(length(obj.gustCases)).gla=[];
                    if ~obj.aeSSM.settings.nonlinActLimits
                        obj.constraints.gla(length(obj.gustCases)).defMax=[];
                        obj.constraints.gla(length(obj.gustCases)).defMin=[];
                        obj.constraints.gla(length(obj.gustCases)).rate=[];
                    end
                    for iCase=1:length(obj.gustCases)
                        if isempty(obj.gustCases(iCase).glaVal)
                            obj.desVar.gustCases(iCase).gla=zeros(1,size(obj.gustCases(iCase).glaSurf,2)*obj.settings.glaControllerOrder);
                            obj.lowerBounds.gustCases(iCase).gla=-Inf*ones(1,size(obj.gustCases(iCase).glaSurf,2)*obj.settings.glaControllerOrder);
                            obj.upperBounds.gustCases(iCase).gla=Inf*ones(1,size(obj.gustCases(iCase).glaSurf,2)*obj.settings.glaControllerOrder);     
                            
                            if ~obj.aeSSM.settings.nonlinActLimits
                                obj.constraints.gla(iCase).defMax=[zeros(1, size(obj.gustCases(iCase).glaSurf,2))];
                                obj.constraints.gla(iCase).defMin=[zeros(1, size(obj.gustCases(iCase).glaSurf,2))];
                                obj.constraints.gla(iCase).rate=[zeros(1, size(obj.gustCases(iCase).glaSurf,2))];
                            end
                        end
                    end
                end
                % constraints
                obj.constraints.gustCases.strength=[];
                for iBeam=1:length(obj.settings.beamsToOptimize) 
                    obj.constraints.gustCases.strength(iBeam).spFr=zeros(1,obj.str.beam(obj.settings.beamsToOptimize(iBeam)).nel);
                    obj.constraints.gustCases.strength(iBeam).spRe=zeros(1,obj.str.beam(obj.settings.beamsToOptimize(iBeam)).nel);
                    obj.constraints.gustCases.strength(iBeam).skLo=zeros(1,obj.str.beam(obj.settings.beamsToOptimize(iBeam)).nel);
                    obj.constraints.gustCases.strength(iBeam).skUp=zeros(1,obj.str.beam(obj.settings.beamsToOptimize(iBeam)).nel);
                    if ~isempty(obj.str.beam(obj.settings.beamsToOptimize(iBeam)).ribLocations)
                        obj.constraints.gustCases.buckl(iBeam).spFr=zeros(1,obj.str.beam(obj.settings.beamsToOptimize(iBeam)).nel*12);
                        obj.constraints.gustCases.buckl(iBeam).spRe=zeros(1,obj.str.beam(obj.settings.beamsToOptimize(iBeam)).nel*12);
                        obj.constraints.gustCases.buckl(iBeam).skLo=zeros(1,obj.str.beam(obj.settings.beamsToOptimize(iBeam)).nel*48);
                        obj.constraints.gustCases.buckl(iBeam).skUp=zeros(1,obj.str.beam(obj.settings.beamsToOptimize(iBeam)).nel*48);
                    end
                end
                obj.constraints.flutter=zeros(1,length(obj.gustCases));
            end
            % if divergence constraint is active; initialize for every case
            if obj.settings.divConstraint
                obj.constraints.div=[];
                if ~isempty(obj.maneuverCases)
                    obj.constraints.div.maneuverCases=zeros(1,length(obj.maneuverCases));
                end
                if ~isempty(obj.gustCases)
                    obj.constraints.div.gustCases=zeros(1,length(obj.gustCases));
                end
            end
                
               
            %% init design variable type vector
            level1Names=fieldnames(obj.desVar);
            types=cell(0);
            for iFieldLvl1=1:size(level1Names,1)
                level2Names=fieldnames(getfield(obj.desVar,level1Names{iFieldLvl1}));
                numLvl1=length(getfield(obj.desVar,level1Names{iFieldLvl1}));
                for iLvl1=1:numLvl1
                    for iFieldLvl2=1:size(level2Names)
                        nVar=length(getfield(getfield(obj.desVar,level1Names{iFieldLvl1},{iLvl1}),level2Names{iFieldLvl2}));
                        if nVar>0
                            varName=[level1Names{iFieldLvl1}, '_' num2str(iLvl1) '_' level2Names{iFieldLvl2}];
                            strings=cellstr([repmat([varName '_'],nVar,1)  num2str((1:nVar)','%04d')]);
                            types=[types; strings];
                        end
                    end
                end
            end
            obj.desVarTypes=types;
            
            obj.desVarScale=zeros(length(obj.desVarTypes),1);
            for iBeam=1:length(obj.settings.beamsToOptimize)
                obj.desVarScale(and(startsWith(obj.desVarTypes,['Str_' num2str(iBeam)]), contains(obj.desVarTypes,'tSkLo')))=obj.settings.thicknessScaleFactorSkin{iBeam};
                obj.desVarScale(and(startsWith(obj.desVarTypes,['Str_' num2str(iBeam)]), contains(obj.desVarTypes,'tSkUp')))=obj.settings.thicknessScaleFactorSkin{iBeam};
                obj.desVarScale(and(startsWith(obj.desVarTypes,['Str_' num2str(iBeam)]), contains(obj.desVarTypes,'tSpFr')))=obj.settings.thicknessScaleFactorSpar{iBeam};
                obj.desVarScale(and(startsWith(obj.desVarTypes,['Str_' num2str(iBeam)]), contains(obj.desVarTypes,'tSpRe')))=obj.settings.thicknessScaleFactorSpar{iBeam};
            end
            obj.desVarScale(and(startsWith(obj.desVarTypes,['Str_']), contains(obj.desVarTypes,'phi')))=obj.settings.stiffnessDirectionScaleFactor;
            obj.desVarScale(contains(obj.desVarTypes,'mla'))=obj.settings.mlaScaleFactor;
            obj.desVarScale(contains(obj.desVarTypes,'rollAlloc'))=obj.settings.rollAllocScaleFactor;
            obj.desVarScale(contains(obj.desVarTypes,'gla'))=obj.settings.glaScaleFactor;
            
            % initialize loadhulldata
            for iBeam=1:length(obj.settings.beamsToOptimize)
                    obj.loadHullData(iBeam).nodDisp=[];
                    obj.loadHullData(iBeam).info.gustCaseID=[];
                    obj.loadHullData(iBeam).info.manCaseID=[];
                    obj.loadHullData(iBeam).info.posNeg=[];
                    obj.loadHullData(iBeam).info.gustLength=[];
                    obj.loadHullData(iBeam).info.timeStep=[];
            end            
        end
        
        function obj=evaluate(obj, varargin)
            %% get flags
            if ~isempty(varargin)
                analysisType=varargin{1};
                caseId=varargin{2};
                fdFlag=varargin{3};
            else %no flags input -> evaluate completely
                analysisType=ones(4,1); %all four flags set to one
                caseId=0; %all cases
                fdFlag=0;
            end
            
            %% structure update and objective evaluation
            tic;
            if analysisType(1) % objective changes only when structure is recalculated
                obj=obj.update();
                obj=obj.evaluateObjective();
            end
            tstr=toc;
            
            %% maneuver analysis
            tic;
            if analysisType(2)
                obj=obj.evaluateManeuverCases(caseId);
            end
            tm=toc;
            tic;
            %% gust analysis
            if ~isempty(obj.gustCases)
                if analysisType(3)
                    obj=obj.evaluateGustCasesNominal(caseId);
                end
                if analysisType(4)
                    obj=obj.evaluateGustCasesDynamic(caseId, 0);
                end
                if or(analysisType(3), analysisType(4))
                    obj=obj.evaluateStructuralFailure(0);
                end                
            end
            tg=toc;
%             fprintf('Total Time: %.2f; STR Time: %.2f; MAN Time: %.2f; GUST Time: %.2f;\n',tstr+tm+tg,tstr,tm,tg)
        end
        
        function obj=update(obj)
            for iBeam=1:length(obj.settings.beamsToOptimize)
                obj.str.beam(iBeam).beam_element_failIdx = 1:obj.str.beam(iBeam).nel;
                beamID=obj.settings.beamsToOptimize(iBeam);
                %get control stations of this beam
                controlElementsThickness=obj.settings.thicknessControlStations{iBeam};
                %if symmetric gather all lengths of the right wing
                if obj.str.beam(beamID).is_sym
                    lengths=[obj.str.beam(beamID).beamelement(obj.str.beam(beamID).nel/2+1:end).le]';
                else
                    lengths=[obj.str.beam(beamID).beamelement(1:end).le]';
                end
                   
                lengthStations=[0; cumsum(lengths(1:end-1))]+lengths/2;
                
                
                
                if obj.str.beam(beamID).anisotropic
                    
                    topSkinThicknessVec=obj.desVar.Str(iBeam).tSkUp;
                    bottomSkinThicknessVec=obj.desVar.Str(iBeam).tSkLo;
                    frontSparThicknessVec=obj.desVar.Str(iBeam).tSpFr;
                    rearSparThicknessVec=obj.desVar.Str(iBeam).tSpRe;
                    topSkinThicknessesHalfWing=interp1(lengthStations(controlElementsThickness),topSkinThicknessVec,lengthStations,'linear','extrap');
                    bottomSkinThicknessesHalfWing=interp1(lengthStations(controlElementsThickness),bottomSkinThicknessVec,lengthStations,'linear','extrap');
                    frontSparThicknessesHalfWing=interp1(lengthStations(controlElementsThickness),frontSparThicknessVec,lengthStations,'linear','extrap');
                    rearSparThicknessesHalfWing=interp1(lengthStations(controlElementsThickness),rearSparThicknessVec,lengthStations,'linear','extrap');
                    if obj.str.beam(beamID).is_sym
                        allTopSkinThicknesses=[topSkinThicknessesHalfWing(end:-1:1); topSkinThicknessesHalfWing];
                        allBottomSkinThicknesses=[bottomSkinThicknessesHalfWing(end:-1:1); bottomSkinThicknessesHalfWing];
                        allFrontSparThicknesses=[frontSparThicknessesHalfWing(end:-1:1); frontSparThicknessesHalfWing];
                        allRearSparThicknesses=[rearSparThicknessesHalfWing(end:-1:1); rearSparThicknessesHalfWing];
                    else
                        allTopSkinThicknesses=[topSkinThicknessesHalfWing];
                        allBottomSkinThicknesses=[bottomSkinThicknessesHalfWing];
                        allFrontSparThicknesses=[frontSparThicknessesHalfWing];
                        allRearSparThicknesses=[rearSparThicknessesHalfWing];
                    end
                    
                    
                    %for anisotropic, the same has to be done for the angles
                    %init angle vectors with zeros as they are also needed
                    %when there is no angles in the design variables
                    allLowerSkinAngles=zeros(obj.str.beam(beamID).nel,1);
                    allUpperSkinAngles=zeros(obj.str.beam(beamID).nel,1);
                    allFrontSparAngles=zeros(obj.str.beam(beamID).nel,1);
                    allRearSparAngles=zeros(obj.str.beam(beamID).nel,1);
                    if  ~isempty(obj.settings.stiffnessDirectionControlStations)
                        if ~isempty(obj.settings.stiffnessDirectionControlStations{iBeam})
                            controlElementsStiffDir=obj.settings.stiffnessDirectionControlStations{iBeam};
                            lowerSkinAnglesVec=obj.desVar.Str(iBeam).phiSkLo;
                            upperSkinAnglesVec=obj.desVar.Str(iBeam).phiSkUp;
                            lowerSkinAnglesHalfWing=interp1(lengthStations(controlElementsStiffDir),lowerSkinAnglesVec,lengthStations,'next','extrap');
                            upperSkinAnglesHalfWing=interp1(lengthStations(controlElementsStiffDir),upperSkinAnglesVec,lengthStations,'next','extrap');
                            if obj.settings.sparTailoring
                                frontSparAnglesVec=obj.desVar.Str(iBeam).phiSpFr;
                                rearSparAnglesVec=obj.desVar.Str(iBeam).phiSpRe;
                                frontSparAnglesHalfWing=interp1(lengthStations(controlElementsStiffDir),frontSparAnglesVec,lengthStations,'next','extrap');
                                rearSparAnglesHalfWing=interp1(lengthStations(controlElementsStiffDir),rearSparAnglesVec,lengthStations,'next','extrap');
                            end
                            if obj.str.beam(beamID).is_sym
                                allLowerSkinAngles=[lowerSkinAnglesHalfWing(end:-1:1); -lowerSkinAnglesHalfWing];
                                allUpperSkinAngles=[upperSkinAnglesHalfWing(end:-1:1); -upperSkinAnglesHalfWing];
                                if obj.settings.sparTailoring
                                    allFrontSparAngles=[frontSparAnglesHalfWing(end:-1:1); -frontSparAnglesHalfWing];
                                    allRearSparAngles=[rearSparAnglesHalfWing(end:-1:1); -rearSparAnglesHalfWing];
                                end
                            else
                                allLowerSkinAngles=[lowerSkinAnglesHalfWing];
                                allUpperSkinAngles=[upperSkinAnglesHalfWing];
                                if obj.settings.sparTailoring
                                    allFrontSparAngles=[frontSparAnglesHalfWing];
                                    allRearSparAngles=[rearSparAnglesHalfWing];
                                end
                            end
                        end
                    end
                    if ~isempty(obj.settings.stiffnessDirectionPolynomialOrder)
                        if ~isempty(obj.settings.stiffnessDirectionPolynomialOrder{iBeam})
                            lowerSkinPolCoeff=obj.desVar.Str(iBeam).phiPolSkLo;
                            upperSkinPolCoeff=obj.desVar.Str(iBeam).phiPolSkUp;
                            lowerSkinAnglesHalfWing=zeros(length(lengthStations),1);
                            upperSkinAnglesHalfWing=zeros(length(lengthStations),1);
                            %map lengthstations on [-1;1]
                            statPol=interp1([lengthStations(1) lengthStations(end)],[-1 1],lengthStations);
                            for iPol=1:obj.settings.stiffnessDirectionPolynomialOrder{iBeam}
                                curPol=(((statPol+sqrt(statPol.^2-1)).^(iPol-1)+(statPol-sqrt(statPol.^2-1)).^(iPol-1))./2);
                                %add polynom of current order multiplied by
                                %coeff
                                lowerSkinAnglesHalfWing=curPol.*lowerSkinPolCoeff(iPol)+lowerSkinAnglesHalfWing;
                                upperSkinAnglesHalfWing=curPol.*upperSkinPolCoeff(iPol)+upperSkinAnglesHalfWing;
                            end
                            %apply bounds
                            lowerSkinAnglesHalfWing=max(min(lowerSkinAnglesHalfWing,obj.settings.stiffnessDirectionBounds(2)),obj.settings.stiffnessDirectionBounds(1));
                            upperSkinAnglesHalfWing=max(min(upperSkinAnglesHalfWing,obj.settings.stiffnessDirectionBounds(2)),obj.settings.stiffnessDirectionBounds(1));
                            
                            if obj.str.beam(beamID).is_sym
                                allLowerSkinAngles=[lowerSkinAnglesHalfWing(end:-1:1); -lowerSkinAnglesHalfWing];
                                allUpperSkinAngles=[upperSkinAnglesHalfWing(end:-1:1); -upperSkinAnglesHalfWing];
                            else
                                allLowerSkinAngles=[lowerSkinAnglesHalfWing];
                                allUpperSkinAngles=[upperSkinAnglesHalfWing];
                            end
                        end
                    end
                    if and(isempty(obj.settings.stiffnessDirectionPolynomialOrder),isempty(obj.settings.stiffnessDirectionControlStations))%no tailoring
                        for iEl = 1:obj.str.beam(beamID).nel
                            allLowerSkinAngles(iEl) = obj.str.beam(beamID).beamelement(iEl).crosssection.laminate_sk_lo.offset_angle;
                            allUpperSkinAngles(iEl) = obj.str.beam(beamID).beamelement(iEl).crosssection.laminate_sk_up.offset_angle;
                            if obj.settings.sparTailoring
                                allFrontSparAngles(iEl) = obj.str.beam(beamID).beamelement(iEl).crosssection.laminate_fs.offset_angle;
                                allRearSparAngles(iEl)  =obj.str.beam(beamID).beamelement(iEl).crosssection.laminate_rs.offset_angle;
                            end
                        end                        
                    end
                    %update beam elements for of anisotropic model
                    if obj.str.beam(beamID).is_sym
                        for iEl=1:obj.str.beam(beamID).nel/2
                            if obj.settings.sparTailoring
                                obj.str.beam(beamID).beamelement(iEl).crosssection=obj.str.beam(beamID).beamelement(iEl).crosssection.crosssection_update(allFrontSparThicknesses(iEl),allTopSkinThicknesses(iEl),allRearSparThicknesses(iEl),allBottomSkinThicknesses(iEl), allFrontSparAngles(iEl),allUpperSkinAngles(iEl),allRearSparAngles(iEl),allLowerSkinAngles(iEl));
                            else
                                obj.str.beam(beamID).beamelement(iEl).crosssection=obj.str.beam(beamID).beamelement(iEl).crosssection.crosssection_update(allFrontSparThicknesses(iEl),allTopSkinThicknesses(iEl),allRearSparThicknesses(iEl),allBottomSkinThicknesses(iEl), 0,allUpperSkinAngles(iEl),0,allLowerSkinAngles(iEl));
                            end
                            obj.str.beam(beamID).beamelement(iEl)=obj.str.beam(beamID).beamelement(iEl).f_calcCrossProp_anisotropic;                    
                        end
                        j=0;
                        for iEl=obj.str.beam(beamID).nel/2+1:obj.str.beam(beamID).nel
                            if obj.settings.sparTailoring
                                obj.str.beam(beamID).beamelement(iEl).crosssection=obj.str.beam(beamID).beamelement(iEl).crosssection.crosssection_update(allFrontSparThicknesses(iEl),allTopSkinThicknesses(iEl),allRearSparThicknesses(iEl),allBottomSkinThicknesses(iEl), allFrontSparAngles(iEl),allUpperSkinAngles(iEl),allRearSparAngles(iEl),allLowerSkinAngles(iEl));
                            else
                                obj.str.beam(beamID).beamelement(iEl).crosssection=obj.str.beam(beamID).beamelement(iEl).crosssection.crosssection_update(allFrontSparThicknesses(iEl),allTopSkinThicknesses(iEl),allRearSparThicknesses(iEl),allBottomSkinThicknesses(iEl), 0,allUpperSkinAngles(iEl),0,allLowerSkinAngles(iEl));
                            end
                            obj.str.beam(beamID).beamelement(iEl)=obj.str.beam(beamID).beamelement(iEl).f_calcCrossProp_anisotropic;                    
                            obj.str.beam(beamID).beamelement(iEl).crosssection.Se=abs(obj.str.beam(beamID).beamelement(iEl-(iEl-obj.str.beam(beamID).nel/2+j)).crosssection.Se).*sign(obj.str.beam(beamID).beamelement(iEl).crosssection.Se);
                            j=j+1;
                        end
                    else
                        for iEl=1:obj.str.beam(beamID).nel
                            if obj.settings.sparTailoring
                                obj.str.beam(beamID).beamelement(iEl).crosssection=obj.str.beam(beamID).beamelement(iEl).crosssection.crosssection_update(allFrontSparThicknesses(iEl),allTopSkinThicknesses(iEl),allRearSparThicknesses(iEl),allBottomSkinThicknesses(iEl), allFrontSparAngles(iEl),allUpperSkinAngles(iEl),allRearSparAngles(iEl),allLowerSkinAngles(iEl));
                            else
                                obj.str.beam(beamID).beamelement(iEl).crosssection=obj.str.beam(beamID).beamelement(iEl).crosssection.crosssection_update(allFrontSparThicknesses(iEl),allTopSkinThicknesses(iEl),allRearSparThicknesses(iEl),allBottomSkinThicknesses(iEl), 0,allUpperSkinAngles(iEl),0,allLowerSkinAngles(iEl));
                            end
                            obj.str.beam(beamID).beamelement(iEl).crosssection=obj.str.beam(beamID).beamelement(iEl).crosssection.crosssection_update(allFrontSparThicknesses(iEl),allTopSkinThicknesses(iEl),allRearSparThicknesses(iEl),allBottomSkinThicknesses(iEl), allFrontSparAngles(iEl),allUpperSkinAngles(iEl),allRearSparAngles(iEl),allLowerSkinAngles(iEl));
                            obj.str.beam(beamID).beamelement(iEl)=obj.str.beam(beamID).beamelement(iEl).f_calcCrossProp_anisotropic;                    
                        end
                    end
                else
                    skinThicknessVec=obj.desVar.Str(iBeam).tSk;
                    sparThicknessVec=obj.desVar.Str(iBeam).tSp;
                    skinThicknessesHalfWing=interp1(lengthStations(controlElementsThickness),skinThicknessVec,lengthStations,'linear','extrap');
                    sparThicknessesHalfWing=interp1(lengthStations(controlElementsThickness),sparThicknessVec,lengthStations,'linear','extrap');
                    if obj.str.beam(beamID).is_sym
                        allSkinThicknesses=[skinThicknessesHalfWing(end:-1:1); skinThicknessesHalfWing];
                        allSparThicknesses=[sparThicknessesHalfWing(end:-1:1); sparThicknessesHalfWing];
                    else
                        allSkinThicknesses=[skinThicknessesHalfWing];
                        allSparThicknesses=[sparThicknessesHalfWing];
                    end
                    %update beam elements for isotropic model
                    for iEl=1:obj.str.beam(beamID).nel
                        obj.str.beam(beamID).beamelement(iEl).crosssection.t_sk_lo=allSkinThicknesses(iEl);
                        obj.str.beam(beamID).beamelement(iEl).crosssection.t_sk_up=allSkinThicknesses(iEl);
                        obj.str.beam(beamID).beamelement(iEl).crosssection.t_sp_fr=allSparThicknesses(iEl);
                        obj.str.beam(beamID).beamelement(iEl).crosssection.t_sp_re=allSparThicknesses(iEl);
                        obj.str.beam(beamID).beamelement(iEl)=obj.str.beam(beamID).beamelement(iEl).f_calcCrossProp;
                    end
                end
                obj.str.beam(beamID).update_K=1;
                obj.str.beam(beamID).update_M=1;
            end
            obj.str=obj.str.f_calc_mass(obj.ac.weights) ;  
            
                if ~isempty(obj.str.modeshapes)
                    oldModeShapes=obj.str.modeshapes;
                end
            obj.str.modeshapes=[];            
            obj.str.modefrequencies=[];
            if obj.settings.clampedWing
                obj.str=obj.str.f_solve;
                obj.str=obj.str.f_solve_modes;
            else
                obj.str=obj.str.f_solve_free_modes(1,0);
                for iMode=1:size(obj.str.modeshapes,1)
                    signOld=sign(sum(oldModeShapes(:,iMode)));
                    signNew=sign(sum(obj.str.modeshapes(:,iMode)));
                    if signOld~=signNew
                        obj.str.modeshapes(:,iMode)=-obj.str.modeshapes(:,iMode);
                    end
                end
                
            end
            
            %set cg to all maneuver and gust cases
            %first set all deflections to zero
            obj.str.nodal_deflections=zeros(size(obj.str.node_coords,1)*6,1);
            cg= obj.str.f_compute_CG';
            
            
            %update AeSSM in case structure changed
            if ~isempty(obj.gustCases)
                obj.aeSSM=obj.aeSSM.updateStrQuick(obj.ac,obj.str);
            end
        end
              
        function obj=evaluateManeuverCases(obj, caseId)
            if caseId==0
                caseVec=1:length(obj.maneuverCases);
            else
                caseVec=caseId;
            end
            
            
            for iCase=caseVec
                %clear loadhulldata for this case if it is not empty
                for iBeam=1:length(obj.settings.beamsToOptimize)
                    idsToClear=obj.loadHullData(iBeam).info.manCaseID==iCase;
                    if ~isempty(idsToClear)
                        obj.loadHullData(iBeam).nodDisp(idsToClear,:)=[];
                        obj.loadHullData(iBeam).info.gustCaseID(idsToClear)=[];
                        obj.loadHullData(iBeam).info.manCaseID(idsToClear)=[];
                        obj.loadHullData(iBeam).info.posNeg(idsToClear)=[];
                        obj.loadHullData(iBeam).info.gustLength(idsToClear)=[];
                        obj.loadHullData(iBeam).info.timeStep(idsToClear)=[];
                    end
                end
                [obj, obj.maneuverCases(iCase), obj.maneuverVLM(iCase), dcp]=obj.aeroelasticTrim( obj.maneuverCases(iCase), obj.maneuverVLM(iCase), 'man', iCase);
                %add loadhull data (
                for iBeam=1:length(obj.settings.beamsToOptimize)
                    beamID=obj.settings.beamsToOptimize(iBeam);
                    obj.loadHullData(iBeam).nodDisp=[obj.loadHullData(iBeam).nodDisp; obj.str.nodal_deflections'];
                    obj.loadHullData(iBeam).info.gustCaseID=[obj.loadHullData(iBeam).info.gustCaseID; 0];
                    obj.loadHullData(iBeam).info.manCaseID=[obj.loadHullData(iBeam).info.manCaseID; iCase];
                    obj.loadHullData(iBeam).info.posNeg=[obj.loadHullData(iBeam).info.posNeg; 0];
                    obj.loadHullData(iBeam).info.gustLength=[obj.loadHullData(iBeam).info.gustLength; 0];
                    obj.loadHullData(iBeam).info.timeStep=[obj.loadHullData(iBeam).info.timeStep; 0];
                end
                %set divergence constraint
                obj.constraints.div.maneuverCases(iCase)=dcp;
                % set control surface constraints for roll cases
                if and(~isempty(obj.maneuverCases(iCase).rollSurf), isempty(obj.maneuverCases(iCase).rollAllocVal))
                    if ~isempty(obj.maneuverCases(iCase).mlaSurf)
                        if isempty(obj.maneuverCases(iCase).mlaVal)
                            totDefP=[obj.maneuverCases(iCase).rollDef+obj.desVar.maneuverCases(iCase).mla];
                            totDefN=[-obj.maneuverCases(iCase).rollDef+obj.desVar.maneuverCases(iCase).mla];
                        else
                            totDefP=[obj.maneuverCases(iCase).rollDef+obj.maneuverCases(iCase).mlaVal];
                            totDefN=[-obj.maneuverCases(iCase).rollDef+obj.maneuverCases(iCase).mlaVal];
                        end
                    else
                        totDefP=[obj.maneuverCases(iCase).rollDef];
                        totDefN=[-obj.maneuverCases(iCase).rollDef];
                    end
                    ub=obj.maneuverCases(iCase).mlaUB;
                    lb=obj.maneuverCases(iCase).mlaLB;
                    defRange=ub-lb;
                    obj.constraints.maneuverCases(iCase).csDefUB=[(totDefP-ub)./defRange (totDefN-ub)./defRange];
                    obj.constraints.maneuverCases(iCase).csDefLB=[-(totDefP-lb)./defRange -(totDefN-lb)./defRange];
                end
                %evaluate str constraints
                for iBeam=1:length(obj.settings.beamsToOptimize)
                    beamID=obj.settings.beamsToOptimize(iBeam);
                    obj.str.beam(beamID)=obj.str.beam(beamID).f_calc_stresses(0);
                    obj.str.beam(beamID)=obj.str.beam(beamID).f_solve_buckl;
                    if obj.str.beam(beamID).anisotropic
                        for iEl=1:obj.str.beam(beamID).nel
                            obj.constraints.maneuverCases(iCase).strength(iBeam).spFr(iEl)=obj.str.beam(beamID).beamelement(iEl).crosssection.laminate_fs.fail_idx_crit;
                            obj.constraints.maneuverCases(iCase).strength(iBeam).spRe(iEl)=obj.str.beam(beamID).beamelement(iEl).crosssection.laminate_rs.fail_idx_crit;
                            obj.constraints.maneuverCases(iCase).strength(iBeam).skLo(iEl)=obj.str.beam(beamID).beamelement(iEl).crosssection.laminate_sk_lo.fail_idx_crit;
                            obj.constraints.maneuverCases(iCase).strength(iBeam).skUp(iEl)=obj.str.beam(beamID).beamelement(iEl).crosssection.laminate_sk_up.fail_idx_crit;
                            fsIdx=obj.str.beam(iBeam).beamelement(iEl).crosssection.elements(:,4)==1;
                            skUpIdx=obj.str.beam(iBeam).beamelement(iEl).crosssection.elements(:,4)==2;
                            rsIdx=obj.str.beam(iBeam).beamelement(iEl).crosssection.elements(:,4)==3;
                            skLoIdx=obj.str.beam(iBeam).beamelement(iEl).crosssection.elements(:,4)==4;
                            if ~isempty(obj.str.beam(beamID).ribLocations)
                                obj.constraints.maneuverCases(iCase).buckl(iBeam).spFr(iEl)=max(max(max(obj.str.beam(beamID).beamelement(iEl).bucklingReserveFactors(:,fsIdx,:))));
                                obj.constraints.maneuverCases(iCase).buckl(iBeam).spRe(iEl)=max(max(max(obj.str.beam(beamID).beamelement(iEl).bucklingReserveFactors(:,rsIdx,:))));
                                obj.constraints.maneuverCases(iCase).buckl(iBeam).skLo(iEl)=max(max(max(obj.str.beam(beamID).beamelement(iEl).bucklingReserveFactors(:,skLoIdx,:))));
                                obj.constraints.maneuverCases(iCase).buckl(iBeam).skUp(iEl)=max(max(max(obj.str.beam(beamID).beamelement(iEl).bucklingReserveFactors(:,skUpIdx,:))));
                            end
                            
                        end
                        
                    else

                        for iEl=1:obj.str.beam(beamID).nel
                                obj.constraints.maneuverCases(iCase).strength(iBeam).spFr(iEl)=(obj.str.beam(beamID).beamelement(iEl).crosssection.sigma_sp_fr*1.5/obj.str.beam(beamID).beamelement(iEl).crosssection.tensile_yield_sp)-1;
                                obj.constraints.maneuverCases(iCase).strength(iBeam).spRe(iEl)=(obj.str.beam(beamID).beamelement(iEl).crosssection.sigma_sp_re*1.5/obj.str.beam(beamID).beamelement(iEl).crosssection.tensile_yield_sp)-1;
                                obj.constraints.maneuverCases(iCase).strength(iBeam).skLo(iEl)=(obj.str.beam(beamID).beamelement(iEl).crosssection.sigma_sk_lo*1.5/obj.str.beam(beamID).beamelement(iEl).crosssection.tensile_yield_sk_l)-1;
                                obj.constraints.maneuverCases(iCase).strength(iBeam).skUp(iEl)=(obj.str.beam(beamID).beamelement(iEl).crosssection.sigma_sk_up*1.5/obj.str.beam(beamID).beamelement(iEl).crosssection.tensile_yield_sk_u)-1;
                        end
                    end
                end
            end
        end
        
        function [obj, flightState, wingaero, dcp]=aeroelasticTrim(obj, flightState, wingaero, caseType, caseID)
            %% set MLA %%
             %set all cs to zero
            for iCs=1:length(obj.ac.control_surfaces_parents)
                obj.ac=obj.ac.f_set_control_surface(obj.ac.control_surfaces_parents(iCs).name,0);
            end
            %set mla
            if ~isempty(flightState.mlaSurf)
                for iSurf=1:length(flightState.mlaSurf)
                    %get symmetry and deflection type for mlaSurf
                    sym=0;
                    diffSymCond=0;
                    location=0;
                    parentID=0;
                    for iCs=1:length(obj.ac.control_surfaces_parents)
                        if strcmp(obj.ac.control_surfaces_parents(iCs).name,flightState.mlaSurf{iSurf})
                            sym=obj.ac.control_surfaces_parents(iCs).is_sym;
                            diffSymCond=~isequal(obj.ac.control_surfaces_parents(iCs).is_sym_defl,flightState.mlaSymDefl(iSurf));
                            location=obj.ac.control_surfaces_parents(iCs).pos;
                            parentID=iCs;
                        end
                    end

                    %check if deflection value for the mla were
                    %selected before the optimisation. i.e. mla are not
                    %design variables in that case

                    if ~isempty(flightState.mlaVal)
                        if location==2

                            obj.ac.control_surfaces_parents(parentID).delta=flightState.mlaVal(iSurf);
                        else
                            if and(sym,diffSymCond)
                                %then deflect left and right part seperately
                                obj.ac=obj.ac.f_set_control_surface([flightState.mlaSurf{iSurf} '_left'],-flightState.mlaVal(iSurf));
                                obj.ac=obj.ac.f_set_control_surface([flightState.mlaSurf{iSurf} '_right'],flightState.mlaVal(iSurf));
                                flightState.aircraft_state.control_deflections=obj.ac.control_deflections;
                            else
                                %otherwise just deflect normally
                                obj.ac=obj.ac.f_set_control_surface(flightState.mlaSurf{iSurf},flightState.mlaVal(iSurf));
                                flightState.aircraft_state.control_deflections=obj.ac.control_deflections;
                            end
                        end

                    else
                        %checking for mid control surfaces/spoilers (as
                        %they have to be deflected by normal vector
                        %rotation)
                        if location==2
                            % it is a handle class, so setting it to ac
                            % is enough; still effective in all vlms
                            if strcmp(caseType,'gust')
                                obj.ac.control_surfaces_parents(parentID).delta=obj.desVar.gustCases(caseID).mla(iSurf);
                            elseif strcmp(caseType,'man')
                                obj.ac.control_surfaces_parents(parentID).delta=obj.desVar.maneuverCases(caseID).mla(iSurf);
                            end
                        else
                            %check if both is true, symmetric cs and different
                            %symmetry condition
                            if strcmp(caseType,'gust')
                                surfDef=obj.desVar.gustCases(caseID).mla(iSurf);
                            elseif strcmp(caseType,'man')
                                surfDef=obj.desVar.maneuverCases(caseID).mla(iSurf);
                            end
                            
                            if and(sym,diffSymCond)
                                %then deflect left and right part seperately
                                obj.ac=obj.ac.f_set_control_surface([flightState.mlaSurf{iSurf} '_left'],-surfDef);
                                obj.ac=obj.ac.f_set_control_surface([flightState.mlaSurf{iSurf} '_right'],surfDef);
                                flightState.aircraft_state.control_deflections=obj.ac.control_deflections;
                            else
                                %otherwise just deflect normally
                                obj.ac=obj.ac.f_set_control_surface(flightState.mlaSurf{iSurf},surfDef);
                                flightState.aircraft_state.control_deflections=obj.ac.control_deflections;
                            end
                        end
                    end
                end
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% preprocessing aeroelastic loop %%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % gen start def vector
            defIn=obj.str.f_get_deflections;
            for i=1:length(obj.str.beam)
                defIn(i).def=defIn(i).def*0;
            end
%             prvDef=def;
% 
            % init parameters
            converged=0;
            j=0;
            error=0;
%             relaxationFactor=0.5;
            %reset alpha and trim surface deflection
            trimSurfId=find(strcmp(obj.ac.control_surfaces,'elevator'));
            flightState.aerodynamic_state.alpha=0;
            flightState.aircraft_state.control_deflections{trimSurfId}=0;
            %initial trim with startdef
            [obj, flightState,wingaero, defOut, diffDef]=f_aero_trim(obj,wingaero,flightState,defIn, caseID);
            %eta estimation
            [~, flightState2,~,~,~]=f_aero_trim(obj,wingaero,flightState,defOut, caseID);
            eta=min(0.7,max(0.3,abs(flightState.aerodynamic_state.alpha/flightState2.aerodynamic_state.alpha)));
            
            %history of angle of attack for debugging purposes
            alphaHist=flightState.aerodynamic_state.alpha;

            while or(~converged,j<2)
                j=j+1;
                
%                 eta=0.72;
                defInPrv=defIn;

                for i=1:length(obj.str.beam)
                    defIn(i).def=defInPrv(i).def+eta*(sqrt(5)-1)/2*diffDef(i).def;
                end

                [obj, flightState,wingaero, defOut, diffDef]=f_aero_trim(obj,wingaero,flightState,defIn, caseID);
          
                alphaHist=[alphaHist flightState.aerodynamic_state.alpha];
                error=sum( sqrt(( cell2mat(struct2cell(defOut'))-cell2mat(struct2cell(defIn'))).^2))/sum( sqrt(cell2mat(struct2cell(defOut')).^2) )*100;
                if  error<5
                    eta=max(eta,0.8);
                elseif error<10
                    eta=max(eta,0.55);
                elseif error<20
                    eta=max(eta,0.4);
                elseif error<30
                    eta=max(eta,0.3);
                end


                if((error<=obj.settings.aeroelasticSolverSettings.convergence_tol)&&(error>=0))
                    converged=1;
                elseif j>obj.settings.aeroelasticSolverSettings.max_it
                    disp('no convergence')
                    converged=1;
                end
            end
            change=abs(diff(abs(alphaHist)));
            %divergence constraint: does the deformation result in a lower
            %absolute Cl value
            dcp=-sign(wingaero.Cl)*(alphaHist(2)-alphaHist(1))/abs(alphaHist(1)); 
            flightState.aerodynamic_state=wingaero.state;
            flightState.def=defOut;
            %%%%%%%%%%%%%%%%%%%%%
            %% nested function %%
            %%%%%%%%%%%%%%%%%%%%%
            function [obj, flightState,wingaero,defOut, diffDef]=f_aero_trim(obj,wingaero,flightState,def,iCase)
                %get roll alloc from desvar or from case
                if ~isempty(obj.desVar.maneuverCases(iCase).rollAlloc)
                    rollAlloc=obj.desVar.maneuverCases(iCase).rollAlloc;
                else
                    rollAlloc=obj.maneuverCases(iCase).rollAllocVal;
                end
                %aero trim aircraft
                obj.ac=obj.ac.compute_deflected_grid(def);
                wingaero=wingaero.set_grid(obj.ac.grid_deflected, obj.ac.panels);
                if obj.settings.clampedWing %check if restrained
                    [obj.ac,flightState,wingaero]=trim_wing_fast(obj.ac,flightState,wingaero);
                else
                    [obj.ac,flightState,wingaero]=trim_aircraft_turbo(obj.ac,flightState,wingaero,def,rollAlloc);
%                     [obj.ac,flightState,wingaero]=trim_aircraft_fast(obj.ac,flightState,wingaero,def);
                end
                %compute loads on beam elements (distributed loads)
                obj.ac=obj.ac.compute_beam_forces(wingaero.F_body*1.0,wingaero.F_aero*0,obj.str);
                
                %set acceleration (needs to be changed for e.g. steady roll maneuver)
                if obj.settings.aeroelasticSolverSettings.CG_accelerations==1
                    acc=[0 0 -9.81 0 0 0]*flightState.load_factor;
                    obj.str=obj.str.f_set_acceleration(acc,wingaero.reference.p_ref);
                end
                
                % set distributed loads
                for iBeam=1:length(obj.str.beam)
                    if  isa(obj.str.beam(iBeam),'class_wing')
                        obj.str.beam(iBeam)=obj.str.beam(iBeam).f_set_aeroloads(obj.ac.wings(iBeam));
                    end
                end
                
                % solve structure
                if obj.settings.clampedWing
                    obj.str=obj.str.f_solve();
                else
                    obj.str=obj.str.f_solve_unrestrained();
                end

                defOut=obj.str.f_get_deflections;
                
                for iBeam=1:length(obj.str.beam)
                    diffDef(iBeam).def=defOut(iBeam).def-def(iBeam).def;
                end
            end


        end
        
        function obj=evaluateGustCasesNominal(obj,caseId)
            
            if caseId==0
                caseVec=1:length(obj.gustCases);
            else
                caseVec=caseId;
            end
            
            for iCase=caseVec
                %clear loadhulldata for this case if loadhulldata is not
                %empty
                for iBeam=1:length(obj.settings.beamsToOptimize)
                    idsToClear=and(obj.loadHullData(iBeam).info.gustCaseID==iCase, obj.loadHullData(iBeam).info.gustLength==0);
                    if ~isempty(idsToClear)
                        obj.loadHullData(iBeam).nodDisp(idsToClear,:)=[];
                        obj.loadHullData(iBeam).info.gustCaseID(idsToClear)=[];
                        obj.loadHullData(iBeam).info.manCaseID(idsToClear)=[];
                        obj.loadHullData(iBeam).info.posNeg(idsToClear)=[];
                        obj.loadHullData(iBeam).info.gustLength(idsToClear)=[];
                        obj.loadHullData(iBeam).info.timeStep(idsToClear)=[];
                    end
                end
                
                [obj, obj.gustCases(iCase), obj.gustVLM(iCase), dcp]=obj.aeroelasticTrim(obj.gustCases(iCase), obj.gustVLM(iCase), 'gust', iCase);
                %append nominal information to loadHullData and clear
                %previous nominal information
                for iBeam=1:length(obj.settings.beamsToOptimize)
                    beamID=obj.settings.beamsToOptimize(iBeam);
                    obj.loadHullData(iBeam).nodDisp=[obj.loadHullData(iBeam).nodDisp; obj.str.nodal_deflections'];
                    obj.loadHullData(iBeam).info.gustCaseID=[obj.loadHullData(iBeam).info.gustCaseID; iCase];
                    obj.loadHullData(iBeam).info.manCaseID=[obj.loadHullData(iBeam).info.manCaseID; 0];
                    obj.loadHullData(iBeam).info.posNeg=[obj.loadHullData(iBeam).info.posNeg; 0];
                    obj.loadHullData(iBeam).info.gustLength=[obj.loadHullData(iBeam).info.gustLength; 0];
                    obj.loadHullData(iBeam).info.timeStep=[obj.loadHullData(iBeam).info.timeStep; 0];
                end
					
                %set divergence constraint
                obj.constraints.div.gustCases(iCase)=dcp;
            end
        end
        
        function obj=evaluateGustCasesDynamic(obj,caseId, fdFlag)
            if caseId==0
                caseVec=1:length(obj.gustCases);
            else
                caseVec=caseId;
            end
            
            for iCase=caseVec
                %clear loadhulldata for this case if loadhulldata is not
                %empty
                
                if ~fdFlag
                    for iBeam=1:length(obj.settings.beamsToOptimize)
                        idsToClear=and(obj.loadHullData(iBeam).info.gustCaseID==iCase, obj.loadHullData(iBeam).info.gustLength~=0);
                        if ~isempty(idsToClear)
                            obj.loadHullData(iBeam).nodDisp(idsToClear,:)=[];
                            obj.loadHullData(iBeam).info.gustCaseID(idsToClear)=[];
                            obj.loadHullData(iBeam).info.manCaseID(idsToClear)=[];
                            obj.loadHullData(iBeam).info.gustLength(idsToClear)=[];
                            obj.loadHullData(iBeam).info.posNeg(idsToClear)=[];
                            obj.loadHullData(iBeam).info.timeStep(idsToClear)=[];
                        end
                    end
                end

                %determination of flightProfileAlleviationFactor
                Fgz=1-(obj.ac.Zmo/76200);
                R2=obj.ac.weights.MZFW/obj.ac.weights.MTOW;
                R1=obj.ac.weights.MLW/obj.ac.weights.MTOW;
                Fgm=sqrt(R2* tan(pi*R1/4));
                Fg0=0.5*(Fgz+Fgm);
                flightProfileAlleviationFactor=Fg0+(1-Fg0)/obj.ac.Zmo*obj.gustCases(iCase).h;

                idNom=and(obj.loadHullData(1).info.gustCaseID==iCase, obj.loadHullData(1).info.gustLength==0);
                for iBeam=1:length(obj.settings.beamsToOptimize)
                    %id of nominal case
                    beamID=obj.settings.beamsToOptimize(iBeam);
                end
                nomDisp=obj.loadHullData.nodDisp(idNom,:);
                
                %in case gla surf is specified, gla controller needs  to be
                %initialized
                %the filter coefficients depend on glaVal (preset) and if
                %empty, the coefficients are design variables
                if ~isempty(obj.gustCases(iCase).glaSurf)
                    outputName=strcat(obj.gustCases(iCase).glaSurf,'_sym')';
                    order=obj.settings.glaControllerOrder;
                    if ~isempty(obj.gustCases(iCase).glaVal)
                        %coeff are preset
                        coeff=obj.gustCases(iCase).glaVal;
                    else
                        %coeff are design variables
                        coeff=obj.desVar.gustCases(iCase).gla;
                    end
                    obj.aeSSM.ctr=class_FF_Controller(outputName,  order, coeff);
                    obj.aeSSM.ctr=obj.aeSSM.ctr.genSSM(obj.gustCases(iCase).glaTs);
                    %clear linSSM as controller changed when linear
                    %actuator limits are handeld by constraints
                    if ~obj.aeSSM.settings.nonlinActLimits
                    end
                end
                obj.aeSSM=obj.aeSSM.updateFdCgAndMass( obj.gustCases(iCase).aircraft_state.weight,obj.gustCases(iCase).aircraft_state.CG_ref,obj.str.f_compute_moment_of_inertia(obj.gustCases(iCase).aircraft_state.CG_ref));

                obj.aeSSM.linSSM=[];
                %unsteady aeroelastic simulation
                if obj.settings.closedLoopGust
                    obj.aeSSM=obj.aeSSM.runGustAnalysisClosedLoop(obj.gustCases(iCase),flightProfileAlleviationFactor,nomLoads);
                else
                    obj.aeSSM=obj.aeSSM.runGustAnalysis(obj.gustCases(iCase),flightProfileAlleviationFactor,obj.settings.nGustLengths);
                end
                % flutter check
                d=eig(obj.aeSSM.linSSM.a);
                
                obj.constraints.flutter(iCase)=constAgg(real(d(abs(d)>0.02)),50);
                
                
                %load dof counter required when multiple beams are optimized
                nLoadDofsPrv=0;
               
                %get hull of node_loadings_loc (Z-Mx)&(Mx-Mt) or stor all
                for iBeam=1:length(obj.settings.beamsToOptimize)
                    beamID=obj.settings.beamsToOptimize(iBeam);
                    %get loads out of gust simulation belonging to this beam
                    [nTimeSteps,nModes,nLengths]=size(obj.aeSSM.gustData.modes);
                    %reshape modal disp
                    modDisp=reshape(permute(obj.aeSSM.gustData.modes,[1 3 2]),nTimeSteps*nLengths,nModes);

                    if or(isempty(obj.aeSSM.ctr),~obj.aeSSM.settings.nonlinActLimits)
                        %add and subtract loads from nomloads
                        modDisp=[modDisp; -modDisp];
                    else
                        % if nonlinear actuator limits are used, the data
                        % already contains pos/neg gusts
                        nLengths=nLengths/2;
                    end
                    
                     % store all disp
                        nodDisp=(obj.aeSSM.strModel.modalBase*modDisp')';
                        [tVec,lVec,signVec]=ind2sub([nTimeSteps,nLengths,2],1:nTimeSteps*nLengths*2);

                        loadHullData{iBeam,iCase}={ [nomDisp+nodDisp],   tVec, lVec, signVec};
                    
                        
                end
    
                if ~obj.aeSSM.settings.nonlinActLimits
                    % calcualte deflection and rate constraints
                    if isempty(obj.gustCases(iCase).glaVal)
                        for iCs=1:length(obj.gustCases(iCase).glaSurf)
                            csId=(startsWith(obj.aeSSM.uvlm.csNames,obj.gustCases(iCase).glaSurf(iCs)));
                            csIdRightPart=(endsWith(obj.aeSSM.uvlm.csNames,'_1')); %first part is enough to look at
                            csId=find(and(csId,csIdRightPart));
                            csDef=squeeze(obj.aeSSM.gustData.csDef(:,csId,:));
                            rate=squeeze(obj.aeSSM.gustData.csRate(:,csId,:));
                            UB=obj.gustCases(iCase).glaUB(iCs)*pi/180;
                            LB=obj.gustCases(iCase).glaLB(iCs)*pi/180;
                            defRange=abs(UB-LB);
                            defPos=csDef;
                            defPos(defPos<0)=0;
                            defNeg=csDef;
                            defNeg(defNeg>0)=0;
                            defMax=(defPos-UB)/defRange;
                            defMin=(LB-defNeg)/defRange;
                            obj.constraints.gla(iCase).defMax(iCs)=constAgg(defMax(:),500);
                            obj.constraints.gla(iCase).defMin(iCs)=constAgg(defMin(:),500);
                            obj.constraints.gla(iCase).rate(iCs)=constAgg((abs(rate(:)))/(obj.aeSSM.act.rateLimit*pi/180)-1,500);

                        end
                    end
                end
            end
            if ~fdFlag
                %merge data of all cases for each beam 
                for iBeam=1:length(obj.settings.beamsToOptimize)
                    beamID=obj.settings.beamsToOptimize(iBeam);
                    if obj.str.beam(beamID).is_sym
                        nCuts=obj.str.beam(beamID).nel/2+1;
                    else
                        nCuts=obj.str.beam(beamID).nel+1;
                    end
                    %gather loads of all cases
                    allLengths=[];
                    allTimeSteps=[];
                    allNodDisp=[];
                    allPosNeg=[];
                    caseId=[];
                    for iCase=1:length(obj.gustCases)
                        allLengths=[allLengths;loadHullData{iBeam,iCase}{3}];
                        allTimeSteps=[allTimeSteps;loadHullData{iBeam,iCase}{2}];
                        allNodDisp=[allNodDisp;loadHullData{iBeam,iCase}{1}];
                        allPosNeg=[allPosNeg;loadHullData{iBeam,iCase}{4}];

                        caseId=[caseId; iCase*ones(size(loadHullData{iBeam,iCase}{1},1),1)];
                    end
                    
                    obj.loadHullData(iBeam).nodDisp=[obj.loadHullData(iBeam).nodDisp; allNodDisp];
                    obj.loadHullData(iBeam).info.gustCaseID=[obj.loadHullData(iBeam).info.gustCaseID; caseId];
                    obj.loadHullData(iBeam).info.manCaseID=[obj.loadHullData(iBeam).info.manCaseID; 0*caseId];
                    obj.loadHullData(iBeam).info.posNeg=[obj.loadHullData(iBeam).info.posNeg; allPosNeg'];
                    obj.loadHullData(iBeam).info.gustLength=[obj.loadHullData(iBeam).info.gustLength; allLengths'];
                    obj.loadHullData(iBeam).info.timeStep=[obj.loadHullData(iBeam).info.timeStep; allTimeSteps'];
                        
                    
                    
                end      
            end
            
                
            
        end
               
        function obj=evaluateStructuralFailure(obj, fdFlag)
           
            aggregationFactor=500;
            %evaluate failure
            for iBeamToOpt=1:length(obj.settings.beamsToOptimize)
                
                %subset of loadhulldata which is taken into account in this
                %evaluation of structural failure
                %only gust cases, and not their nominal loads
                subset=find(and(obj.loadHullData(iBeamToOpt).info.gustCaseID~=0, obj.loadHullData(iBeamToOpt).info.gustLength~=0));

                beamID=obj.settings.beamsToOptimize(iBeamToOpt);
                dofBeam=obj.str.sort_vec((obj.str.dof_node_beam(:,3)==beamID),1);
                nEl=obj.str.beam(beamID).nel;
                nP=size(subset,1);
                %create vector of nodal displacements with repeated nodes
                %2:end-1 to get nodal displacements stacked for each element
                elNodDispDof=reshape(sort([reshape(1:(nEl+1)*6,6,nEl+1) reshape(7:(nEl)*6,6,nEl-1)],2),nEl*2*6,1);
                elNodDisp=obj.loadHullData(iBeamToOpt).nodDisp(subset,dofBeam(elNodDispDof));
               

                for iEl=1:obj.str.beam(beamID).nel/2 %1:length(optimCase.str.beam(beamID).beamelement) 
                    %transform element nodal displacements to beam element frame
                    elDispDof=(iEl-1)*12+1:(iEl)*12;
                    elNodDispLoc=obj.str.beam(beamID).beamelement(iEl).T*elNodDisp(:,elDispDof)';
                    %calculate reaction forces
                    elNodRe=obj.str.beam(beamID).beamelement(iEl).elK*elNodDispLoc;
                    % transform to cross sectional modeller coordinate system
                    % exchange x and y
                    elNodRe=elNodRe([2 1 3 5 4 6 8 7 9 11 10 12],:);
                    % flip resulting y
                    elNodRe(2:3:12,:)=-elNodRe(2:3:12,:);

                    %element dofs in reaction force matrix
                    elReDof=(iEl-1)*12+1:(iEl)*12;
                    ribLoc1=obj.str.beam(beamID).beamelement(iEl).ribLocations(1);
                    ribLoc2=obj.str.beam(beamID).beamelement(iEl).ribLocations(end);
                    for iComp=1:4 %1:4
                        %shell element ids in this cross section component 
                        crossEl=find(obj.str.beam(beamID).beamelement(iEl).crosssection.elements(:,4)==iComp);
                        nCrossEl=length(crossEl);
                        %get gamma matrix rows of this component 
                        gammaRows=(min(crossEl)-1)*6+sort([[0:nCrossEl-1]*6+1 [0:nCrossEl-1]*6+2 [0:nCrossEl-1]*6+3]);

                        % calculate element strains due to reaction forces on both ends
                        elStrainEndA=(obj.str.beam(beamID).beamelement(iEl).crosssection.Gamma(gammaRows,:)*(-elNodRe(1:6,:)))';
                        elStrainEndB=(obj.str.beam(beamID).beamelement(iEl).crosssection.Gamma(gammaRows,:)*(elNodRe(7:12,:)))';
                        % element strains of this component
                        compElStrainEndA=[reshape(elStrainEndA()',3,nCrossEl*nP)];
                        compElStrainEndB=[reshape(elStrainEndB()',3,nCrossEl*nP)];
                        if or(iComp==2,iComp==4)
                            maxCompElStrainEndAx=max(compElStrainEndA(1,:));
                            maxCompElStrainEndAy=max(compElStrainEndA(2,:));
                            minCompElStrainEndAx=min(compElStrainEndA(1,:));
                            minCompElStrainEndAy=min(compElStrainEndA(2,:));
                            maxCompElStrainEndBx=max(compElStrainEndB(1,:));
                            maxCompElStrainEndBy=max(compElStrainEndB(2,:));
                            minCompElStrainEndBx=min(compElStrainEndB(1,:));
                            minCompElStrainEndBy=min(compElStrainEndB(2,:));

                            maxCompElStrainEndAxAgg=constAgg((compElStrainEndA(1,:)-minCompElStrainEndAx)/(maxCompElStrainEndAx-minCompElStrainEndAx),aggregationFactor)*(maxCompElStrainEndAx-minCompElStrainEndAx)+minCompElStrainEndAx;
                            maxCompElStrainEndAyAgg=constAgg((compElStrainEndA(2,:)-minCompElStrainEndAy)/(maxCompElStrainEndAy-minCompElStrainEndAy),aggregationFactor)*(maxCompElStrainEndAy-minCompElStrainEndAy)+minCompElStrainEndAy;
                            minCompElStrainEndAxAgg=-(constAgg((-compElStrainEndA(1,:)+maxCompElStrainEndAx)/(maxCompElStrainEndAx-minCompElStrainEndAx),aggregationFactor)*(maxCompElStrainEndAx-minCompElStrainEndAx)-maxCompElStrainEndAx);
                            minCompElStrainEndAyAgg=-(constAgg((-compElStrainEndA(2,:)+maxCompElStrainEndAy)/(maxCompElStrainEndAy-minCompElStrainEndAy),aggregationFactor)*(maxCompElStrainEndAy-minCompElStrainEndAy)-maxCompElStrainEndAy);

                            maxCompElStrainEndBxAgg=constAgg((compElStrainEndB(1,:)-minCompElStrainEndBx)/(maxCompElStrainEndBx-minCompElStrainEndBx),aggregationFactor)*(maxCompElStrainEndBx-minCompElStrainEndBx)+minCompElStrainEndBx;
                            maxCompElStrainEndByAgg=constAgg((compElStrainEndB(2,:)-minCompElStrainEndBy)/(maxCompElStrainEndBy-minCompElStrainEndBy),aggregationFactor)*(maxCompElStrainEndBy-minCompElStrainEndBy)+minCompElStrainEndBy;
                            minCompElStrainEndBxAgg=-(constAgg((-compElStrainEndB(1,:)+maxCompElStrainEndBx)/(maxCompElStrainEndBx-minCompElStrainEndBx),aggregationFactor)*(maxCompElStrainEndBx-minCompElStrainEndBx)-maxCompElStrainEndBx);
                            minCompElStrainEndByAgg=-(constAgg((-compElStrainEndB(2,:)+maxCompElStrainEndBy)/(maxCompElStrainEndBy-minCompElStrainEndBy),aggregationFactor)*(maxCompElStrainEndBy-minCompElStrainEndBy)-maxCompElStrainEndBy);
                            bucklStrainsEndA=[maxCompElStrainEndAxAgg minCompElStrainEndAxAgg; minCompElStrainEndAyAgg maxCompElStrainEndAyAgg;0 0];
                            bucklStrainsEndB=[maxCompElStrainEndBxAgg minCompElStrainEndBxAgg; minCompElStrainEndByAgg maxCompElStrainEndByAgg;0 0];
                            strainsAtRibs=[bucklStrainsEndA+ribLoc1.*(bucklStrainsEndB-bucklStrainsEndA) bucklStrainsEndA+ribLoc2.*(bucklStrainsEndB-bucklStrainsEndA)];
                        else
                            aggAbsShearA=constAgg(abs(compElStrainEndA(3,:))./max(abs(compElStrainEndA(3,:))),aggregationFactor)*max(abs(compElStrainEndA(3,:)));
                            aggAbsShearB=constAgg(abs(compElStrainEndB(3,:))./max(abs(compElStrainEndB(3,:))),aggregationFactor)*max(abs(compElStrainEndB(3,:)));
                            bucklStrainsEndA=[0; 0; aggAbsShearA];
                            bucklStrainsEndB=[0; 0; aggAbsShearB];
                            strainsAtRibs=bucklStrainsEndA+[ribLoc1.*(bucklStrainsEndB-bucklStrainsEndA) ribLoc2.*(bucklStrainsEndB-bucklStrainsEndA)];
                        end





                        % compute all principal Strains
                        compElStrainAvg=([compElStrainEndA(1,:) compElStrainEndB(1,:)]+[compElStrainEndA(2,:) compElStrainEndB(2,:)])/2;
                        compElStrainR=sqrt((([compElStrainEndA(1,:) compElStrainEndB(1,:)]-[compElStrainEndA(2,:) compElStrainEndB(2,:)])./2).^2+([compElStrainEndA(3,:) compElStrainEndB(3,:)]./2).^2);
                        compElStrainE1=compElStrainAvg+compElStrainR;
                        compElStrainE2=compElStrainAvg-compElStrainR;
                        % calc all failure indices
                        switch iComp
                            case 1
                                obj.str.beam(beamID).beamelement(iEl).crosssection.laminate_fs= ...
                                obj.str.beam(beamID).beamelement(iEl).crosssection.laminate_fs.calculate_safety_factor_from_princ(compElStrainE1, compElStrainE2, obj.str.beam(beamID).beamelement(iEl).crosssection.safety_factor, aggregationFactor);

                                r=obj.str.beam(beamID).beamelement(iEl).bucklElements(iComp).calcBucklingReserveFactors(strainsAtRibs, obj.str.beam(beamID).beamelement(iEl).crosssection.laminate_fs.ABD_stiff);
                                fail(iComp,iEl)=obj.str.beam(beamID).beamelement(iEl).crosssection.laminate_fs.fail_idx_crit;         
                                buckl(iComp,iEl)=constAgg(r(1,:),aggregationFactor)-1;
                                obj.str.beam(beamID).beamelement(iEl).crosssection.laminate_fs.fail_idx_arr=[];
                            case 2
                                obj.str.beam(beamID).beamelement(iEl).crosssection.laminate_sk_up= ...
                                obj.str.beam(beamID).beamelement(iEl).crosssection.laminate_sk_up.calculate_safety_factor_from_princ(compElStrainE1, compElStrainE2, obj.str.beam(beamID).beamelement(iEl).crosssection.safety_factor, aggregationFactor);

                                r=obj.str.beam(beamID).beamelement(iEl).bucklElements(iComp).calcBucklingReserveFactors(strainsAtRibs, obj.str.beam(beamID).beamelement(iEl).crosssection.laminate_sk_up.ABD_stiff);

                                fail(iComp,iEl)=obj.str.beam(beamID).beamelement(iEl).crosssection.laminate_sk_up.fail_idx_crit;        
                                buckl(iComp,iEl)=constAgg(r(1,:),aggregationFactor)-1;
                                obj.str.beam(beamID).beamelement(iEl).crosssection.laminate_sk_up.fail_idx_arr=[];

                            case 3
                                obj.str.beam(beamID).beamelement(iEl).crosssection.laminate_rs= ...
                                obj.str.beam(beamID).beamelement(iEl).crosssection.laminate_rs.calculate_safety_factor_from_princ(compElStrainE1, compElStrainE2, obj.str.beam(beamID).beamelement(iEl).crosssection.safety_factor, aggregationFactor);

                                r=obj.str.beam(beamID).beamelement(iEl).bucklElements(iComp).calcBucklingReserveFactors(strainsAtRibs, obj.str.beam(beamID).beamelement(iEl).crosssection.laminate_rs.ABD_stiff);

                                fail(iComp,iEl)=obj.str.beam(beamID).beamelement(iEl).crosssection.laminate_rs.fail_idx_crit;            
                                buckl(iComp,iEl)=constAgg(r(1,:),aggregationFactor)-1;
                                obj.str.beam(beamID).beamelement(iEl).crosssection.laminate_rs.fail_idx_arr=[];
                            case 4
                                obj.str.beam(beamID).beamelement(iEl).crosssection.laminate_sk_lo= ...
                                obj.str.beam(beamID).beamelement(iEl).crosssection.laminate_sk_lo.calculate_safety_factor_from_princ(compElStrainE1, compElStrainE2, obj.str.beam(beamID).beamelement(iEl).crosssection.safety_factor, aggregationFactor);
                                r=obj.str.beam(beamID).beamelement(iEl).bucklElements(iComp).calcBucklingReserveFactors(strainsAtRibs, obj.str.beam(beamID).beamelement(iEl).crosssection.laminate_sk_lo.ABD_stiff);

                                fail(iComp,iEl)=obj.str.beam(beamID).beamelement(iEl).crosssection.laminate_sk_lo.fail_idx_crit;        
                                buckl(iComp,iEl)=constAgg(r(1,:),aggregationFactor)-1;   
                                obj.str.beam(beamID).beamelement(iEl).crosssection.laminate_sk_lo.fail_idx_arr=[];
                        end
                    end
                end

                %store failure information in constraints
                obj.constraints.gustCases.strength(iBeamToOpt).spFr=[fail(1,1:end) fail(1,end:-1:1)];
                obj.constraints.gustCases.strength(iBeamToOpt).spRe=[fail(3,1:end) fail(3,end:-1:1)];
                obj.constraints.gustCases.strength(iBeamToOpt).skLo=[fail(4,1:end) fail(4,end:-1:1)];
                obj.constraints.gustCases.strength(iBeamToOpt).skUp=[fail(2,1:end) fail(2,end:-1:1)];
                if ~isempty(obj.str.beam(beamID).ribLocations)
                    obj.constraints.gustCases.buckl(iBeamToOpt).spFr=[buckl(1,1:end) buckl(1,end:-1:1)];
                    obj.constraints.gustCases.buckl(iBeamToOpt).spRe=[buckl(3,1:end) buckl(3,end:-1:1)];
                    obj.constraints.gustCases.buckl(iBeamToOpt).skLo=[buckl(4,1:end) buckl(4,end:-1:1)];
                    obj.constraints.gustCases.buckl(iBeamToOpt).skUp=[buckl(2,1:end) buckl(2,end:-1:1)];
                end
            end
            
        end
        
        function obj=evaluateObjective(obj)
            obj.objective=obj.str.beam(1).m_wingbox_total;
        end
        
        function obj=setScaledDesignVariables(obj, scaledX)
            %check input size
            if length(obj.getScaledDesignVariables)~=length(scaledX)
                disp('length of input is wrong')
            end

            % first all for first beam to optimize, then all for second
            % beam to optimize
            for iBeam=1:length(obj.settings.beamsToOptimize)
                nCtrThickness=length(obj.settings.thicknessControlStations{iBeam});
                nCtrDir=0;
                scaleSkin=obj.settings.thicknessScaleFactorSkin{iBeam};
                scaleSpar=obj.settings.thicknessScaleFactorSpar{iBeam};
                
                if obj.str.beam(obj.settings.beamsToOptimize(iBeam)).anisotropic
                    
                    % order of the design variables:
                    % tSkLo, tSkUp,tSpFr,tSpRe, phiSkLo,phiSkUp,phiSpFr, phiSpRe
                    %set thickness
                    obj.desVar.Str(iBeam).tSkLo=scaledX(1:nCtrThickness)./scaleSkin;
                    obj.desVar.Str(iBeam).tSkUp=scaledX(nCtrThickness+1:2*nCtrThickness)./scaleSkin;
                    obj.desVar.Str(iBeam).tSpFr=scaledX(2*nCtrThickness+1:3*nCtrThickness)./scaleSpar;
                    obj.desVar.Str(iBeam).tSpRe=scaledX(3*nCtrThickness+1:4*nCtrThickness)./scaleSpar;
                    %set angles if required
                    if obj.str.beam(obj.settings.beamsToOptimize(iBeam)).anisotropic
                        if ~isempty(obj.settings.stiffnessDirectionControlStations)
                            if ~isempty(obj.settings.stiffnessDirectionControlStations{iBeam})
                                nCtrDir=length(obj.settings.stiffnessDirectionControlStations{iBeam});
                                obj.desVar.Str(iBeam).phiSkLo=scaledX(4*nCtrThickness+1:4*nCtrThickness+nCtrDir)./obj.settings.stiffnessDirectionScaleFactor;
                                obj.desVar.Str(iBeam).phiSkUp=scaledX(4*nCtrThickness+nCtrDir+1:4*nCtrThickness+2*nCtrDir)./obj.settings.stiffnessDirectionScaleFactor;
                                if obj.settings.sparTailoring
                                    obj.desVar.Str(iBeam).phiSpFr=scaledX(4*nCtrThickness+2*nCtrDir+1:4*nCtrThickness+3*nCtrDir)./obj.settings.stiffnessDirectionScaleFactor;
                                    obj.desVar.Str(iBeam).phiSpRe=scaledX(4*nCtrThickness+3*nCtrDir+1:4*nCtrThickness+4*nCtrDir)./obj.settings.stiffnessDirectionScaleFactor;
                                end
                            end
                        elseif ~isempty(obj.settings.stiffnessDirectionPolynomialOrder)
                            if ~isempty(obj.settings.stiffnessDirectionPolynomialOrder{iBeam})
                                nCtrDir=obj.settings.stiffnessDirectionPolynomialOrder{iBeam};
                                obj.desVar.Str(iBeam).phiPolSkLo=scaledX(4*nCtrThickness+1:4*nCtrThickness+nCtrDir)./obj.settings.stiffnessDirectionScaleFactor;
                                obj.desVar.Str(iBeam).phiPolSkUp=scaledX(4*nCtrThickness+nCtrDir+1:4*nCtrThickness+2*+nCtrDir)./obj.settings.stiffnessDirectionScaleFactor;
                               
                            end
                        end
                    end
                    scaledX=scaledX(4*nCtrThickness+2*nCtrDir+2*obj.settings.sparTailoring*nCtrDir+1:end);
                else
                    %order Spar Skin
                    obj.desVar.Str(iBeam).tSp=scaledX(1:nCtrThickness)./scaleSpar;
                    obj.desVar.Str(iBeam).tSk=scaledX(nCtrThickness+1:2*nCtrThickness)./scaleSkin;
                    scaledX=scaledX(2*nCtrThickness+4*nCtrDir+1:end);
                end
            end
            %mla at man cases
            for iCase=1:length(obj.maneuverCases)
                if ~isempty(obj.maneuverCases(iCase).mlaSurf) && isempty(obj.maneuverCases(iCase).mlaVal)
                    n=length(obj.maneuverCases(iCase).mlaSurf);
                    obj.desVar.maneuverCases(iCase).mla=scaledX(1:n)/obj.settings.mlaScaleFactor;
                    scaledX=scaledX(n+1:end);
                    if ~isempty(obj.maneuverCases(iCase).rollSurf)
                        n=length(obj.maneuverCases(iCase).rollSurf);
                        obj.desVar.maneuverCases(iCase).rollAlloc=scaledX(1:n)/obj.settings.rollAllocScaleFactor;
                        scaledX=scaledX(n+1:end);
                    end
                end
                
            end
            %gust cases
            for iCase=1:length(obj.gustCases)
                %MLA
                if ~isempty(obj.gustCases(iCase).mlaSurf) && isempty(obj.gustCases(iCase).mlaVal)
                    n=length(obj.gustCases(iCase).mlaSurf);
                    obj.desVar.gustCases(iCase).mla=scaledX(1:n)/obj.settings.mlaScaleFactor;
                    scaledX=scaledX(n+1:end);
                end
                %GLA
                if ~isempty(obj.gustCases(iCase).glaSurf) && isempty(obj.gustCases(iCase).glaVal)
                    n=length(obj.gustCases(iCase).glaSurf)*obj.settings.glaControllerOrder;
                    obj.desVar.gustCases(iCase).gla=scaledX(1:n)/obj.settings.glaScaleFactor;
                    scaledX=scaledX(n+1:end);
                end
            end
        end
        
        function obj=setUnScaledDesignVariables(obj, unScaledX)
            % get scales
            scales=obj.desVarScale';
            
            % multiply unscaledX with scales
            scaledX=unScaledX.*scales;
            
            % call set scaled
            obj=obj.setScaledDesignVariables(scaledX);
            
        end
        
        function [x]=getScaledDesignVariables(obj)
            %total number of design variables is number of str control
            %stations times 2 (skin/spar) + the mla variables
            xStr=[];
            for iBeam=1:length(obj.settings.beamsToOptimize)
                
                scaleSkin=obj.settings.thicknessScaleFactorSkin{iBeam};
                scaleSpar=obj.settings.thicknessScaleFactorSpar{iBeam};
                if obj.str.beam(obj.settings.beamsToOptimize(iBeam)).anisotropic
                    % order of the design variables:
                    % tSkLo, tSkUp,tSpFr,tSpRe, phiSkLo,phiSkUp,phiSpFr, phiSpRe
                    
                    xStr=[xStr, obj.desVar.Str(iBeam).tSkLo.*scaleSkin, obj.desVar.Str(iBeam).tSkUp.*scaleSkin, obj.desVar.Str(iBeam).tSpFr.*scaleSpar,  obj.desVar.Str(iBeam).tSpRe.*scaleSpar];
                    if obj.str.beam(obj.settings.beamsToOptimize(iBeam)).anisotropic
                        if ~isempty(obj.settings.stiffnessDirectionControlStations)
                            if obj.settings.sparTailoring
                                xStr=[xStr, [obj.desVar.Str(iBeam).phiSkLo, obj.desVar.Str(iBeam).phiSkUp,obj.desVar.Str(iBeam).phiSpFr, obj.desVar.Str(iBeam).phiSpRe  ]*obj.settings.stiffnessDirectionScaleFactor];
                            else
                                xStr=[xStr, [obj.desVar.Str(iBeam).phiSkLo, obj.desVar.Str(iBeam).phiSkUp]*obj.settings.stiffnessDirectionScaleFactor];
                            end
                        end
                        if ~isempty(obj.settings.stiffnessDirectionPolynomialOrder)
                            if obj.settings.sparTailoring
                                xStr=[xStr, [obj.desVar.Str(iBeam).phiPolSkLo, obj.desVar.Str(iBeam).phiPolSkUp,obj.desVar.Str(iBeam).phiPolSpFr, obj.desVar.Str(iBeam).phiPolSpRe  ]*obj.settings.stiffnessDirectionScaleFactor];
                            else
                                xStr=[xStr, [obj.desVar.Str(iBeam).phiPolSkLo, obj.desVar.Str(iBeam).phiPolSkUp]*obj.settings.stiffnessDirectionScaleFactor];
                            end
                        end
                    end
                else
                    %order Spar Skin
                    xStr=[xStr, [obj.desVar.Str(iBeam).tSp.*scaleSpar, obj.desVar.Str(iBeam).tSk.*scaleSkin]];
                end
            end
            
            xMan=[];
            for iCase=1:length(obj.maneuverCases)
                if isempty(obj.maneuverCases(iCase).mlaVal)
                    xMan=[xMan, obj.desVar.maneuverCases(iCase).mla*obj.settings.mlaScaleFactor];
                end
                
                if ~isempty(obj.maneuverCases(iCase).rollSurf)
                    xMan=[xMan, obj.desVar.maneuverCases(iCase).rollAlloc*obj.settings.rollAllocScaleFactor];
                end
            end
            
            xGust=[];
            for iCase=1:length(obj.gustCases)
                if isempty(obj.gustCases(iCase).mlaVal)
                    xGust=[xGust, obj.desVar.gustCases(iCase).mla*obj.settings.mlaScaleFactor];
                end
                if isempty(obj.gustCases(iCase).glaVal)
                    xGust=[xGust, obj.desVar.gustCases(iCase).gla*obj.settings.glaScaleFactor];
                end
            end
            x=[xStr xMan xGust];
        end
                
        function [lb,ub]=getScaledBounds(obj)
            lbStr=[];
            ubStr=[];
            for iBeam=1:length(obj.settings.beamsToOptimize)
                scaleSkin=obj.settings.thicknessScaleFactorSkin{iBeam};
                scaleSpar=obj.settings.thicknessScaleFactorSpar{iBeam};
                if obj.str.beam(obj.settings.beamsToOptimize(iBeam)).anisotropic
                    % order of the design variables:
                    % tSkLo, tSkUp,tSpFr,tSpRe, phiSkLo,phiSkUp,phiSpFr, phiSpRe
                    
                    lbStr=[lbStr, [obj.lowerBounds.Str(iBeam).tSkLo.*scaleSkin, obj.lowerBounds.Str(iBeam).tSkUp.*scaleSkin, obj.lowerBounds.Str(iBeam).tSpFr.*scaleSpar, obj.lowerBounds.Str(iBeam).tSpRe.*scaleSpar]];
                    ubStr=[ubStr, [obj.upperBounds.Str(iBeam).tSkLo.*scaleSkin, obj.upperBounds.Str(iBeam).tSkUp.*scaleSkin, obj.upperBounds.Str(iBeam).tSpFr.*scaleSpar, obj.upperBounds.Str(iBeam).tSpRe.*scaleSpar]];
                    if obj.str.beam(obj.settings.beamsToOptimize(iBeam)).anisotropic
                        if ~isempty(obj.settings.stiffnessDirectionControlStations)
                            if obj.settings.sparTailoring
                                lbStr=[lbStr, [obj.lowerBounds.Str(iBeam).phiSkLo, obj.lowerBounds.Str(iBeam).phiSkUp,obj.lowerBounds.Str(iBeam).phiSpFr, obj.lowerBounds.Str(iBeam).phiSpRe  ]*obj.settings.stiffnessDirectionScaleFactor];
                                ubStr=[ubStr, [obj.upperBounds.Str(iBeam).phiSkLo, obj.upperBounds.Str(iBeam).phiSkUp,obj.upperBounds.Str(iBeam).phiSpFr, obj.upperBounds.Str(iBeam).phiSpRe  ]*obj.settings.stiffnessDirectionScaleFactor];
                            else
                                lbStr=[lbStr, [obj.lowerBounds.Str(iBeam).phiSkLo, obj.lowerBounds.Str(iBeam).phiSkUp]*obj.settings.stiffnessDirectionScaleFactor];
                                ubStr=[ubStr, [obj.upperBounds.Str(iBeam).phiSkLo, obj.upperBounds.Str(iBeam).phiSkUp]*obj.settings.stiffnessDirectionScaleFactor];
                            end
                        end
                        
                        if ~isempty(obj.settings.stiffnessDirectionPolynomialOrder)
                            if obj.settings.sparTailoring
                                lbStr=[lbStr, [obj.lowerBounds.Str(iBeam).phiPolSkLo, obj.lowerBounds.Str(iBeam).phiPolSkUp,obj.lowerBounds.Str(iBeam).phiPolSpFr, obj.lowerBounds.Str(iBeam).phiPolSpRe  ]*obj.settings.stiffnessDirectionScaleFactor];
                                ubStr=[ubStr, [obj.upperBounds.Str(iBeam).phiPolSkLo, obj.upperBounds.Str(iBeam).phiPolSkUp,obj.upperBounds.Str(iBeam).phiPolSpFr, obj.upperBounds.Str(iBeam).phiPolSpRe  ]*obj.settings.stiffnessDirectionScaleFactor];
                            else
                                lbStr=[lbStr, [obj.lowerBounds.Str(iBeam).phiPolSkLo, obj.lowerBounds.Str(iBeam).phiPolSkUp]*obj.settings.stiffnessDirectionScaleFactor];
                                ubStr=[ubStr, [obj.upperBounds.Str(iBeam).phiPolSkLo, obj.upperBounds.Str(iBeam).phiPolSkUp]*obj.settings.stiffnessDirectionScaleFactor];
                            end
                        end
                    end
                else
                    lbStr=[lbStr, [obj.lowerBounds.Str(iBeam).tSp.*scaleSpar, obj.lowerBounds.Str(iBeam).tSk.*scaleSkin]];
                    ubStr=[ubStr, [obj.upperBounds.Str(iBeam).tSp.*scaleSpar, obj.upperBounds.Str(iBeam).tSk.*scaleSkin]];
                end
            end
            
            lbMla=[];
            ubMla=[];
            for iCase=1:length(obj.maneuverCases)
                if isempty(obj.maneuverCases(iCase).mlaVal)
                    lbMla=[lbMla, obj.lowerBounds.maneuverCases(iCase).mla*obj.settings.mlaScaleFactor];
                    ubMla=[ubMla, obj.upperBounds.maneuverCases(iCase).mla*obj.settings.mlaScaleFactor];
                end
                if and(~isempty(obj.maneuverCases(iCase).rollSurf), isempty(obj.maneuverCases(iCase).rollAllocVal))
                    lbMla=[lbMla, obj.lowerBounds.maneuverCases(iCase).rollAlloc*obj.settings.rollAllocScaleFactor];
                    ubMla=[ubMla, obj.upperBounds.maneuverCases(iCase).rollAlloc*obj.settings.rollAllocScaleFactor];
                end
            end
            lbGust=[];
            ubGust=[];
            for iCase=1:length(obj.gustCases)
                if isempty(obj.gustCases(iCase).mlaVal)
                    lbGust=[lbGust, obj.lowerBounds.gustCases(iCase).mla*obj.settings.mlaScaleFactor];
                    ubGust=[ubGust, obj.upperBounds.gustCases(iCase).mla*obj.settings.mlaScaleFactor];
                end
                if ~isempty(obj.gustCases(iCase).glaSurf)
                    if isempty(obj.gustCases(iCase).glaVal)
                        lbGust=[lbGust, obj.lowerBounds.gustCases(iCase).gla*obj.settings.glaScaleFactor];
                        ubGust=[ubGust, obj.upperBounds.gustCases(iCase).gla*obj.settings.glaScaleFactor];
                    end
                end
            end
            lb=[lbStr, lbMla, lbGust];
            ub=[ubStr, ubMla, ubGust];
        end
        
        function [c]=getConstraints(obj)
            [cMan]=obj.getManConstraints();
            [cGust]=obj.getGustConstraints();
            c=[cMan cGust];
        end
        
        function [cMan]=getManConstraints(obj)
            %here it is possible to implement constraint aggregation
            %methods
            cMan=[];
            for iCase=1:length(obj.maneuverCases)
                for iBeam=1:length(obj.settings.beamsToOptimize)
                    cMan=[cMan ...
                        obj.constraints.maneuverCases(iCase).strength(iBeam).spFr...
                        obj.constraints.maneuverCases(iCase).strength(iBeam).spRe...
                        obj.constraints.maneuverCases(iCase).strength(iBeam).skLo...
                        obj.constraints.maneuverCases(iCase).strength(iBeam).skUp];
                    
                    cMan=[cMan ...
                        obj.constraints.maneuverCases(iCase).buckl(iBeam).spFr...
                        obj.constraints.maneuverCases(iCase).buckl(iBeam).spRe...
                        obj.constraints.maneuverCases(iCase).buckl(iBeam).skLo...
                        obj.constraints.maneuverCases(iCase).buckl(iBeam).skUp];
                end
                if obj.settings.divConstraint
                    cMan=[cMan obj.constraints.div.maneuverCases(iCase)+10^-4];
                end
                
                if and(~isempty(obj.maneuverCases(iCase).rollSurf), isempty(obj.maneuverCases(iCase).rollAllocVal))
                    cMan=[cMan obj.constraints.maneuverCases(iCase).csDefLB obj.constraints.maneuverCases(iCase).csDefUB];
                end
            end
        end
        
        function [cGust]=getGustConstraints(obj)
            cGust=[];
            if ~isempty(obj.gustCases)
                for iBeam=1:length(obj.settings.beamsToOptimize)
                     cGust=[cGust ...
                        obj.constraints.gustCases.strength(iBeam).spFr...
                        obj.constraints.gustCases.strength(iBeam).spRe...
                        obj.constraints.gustCases.strength(iBeam).skLo...
                        obj.constraints.gustCases.strength(iBeam).skUp];
                    
                    cGust=[cGust ...
                        obj.constraints.gustCases.buckl(iBeam).spFr...
                        obj.constraints.gustCases.buckl(iBeam).spRe...
                        obj.constraints.gustCases.buckl(iBeam).skLo...
                        obj.constraints.gustCases.buckl(iBeam).skUp];
                end
                if obj.settings.divConstraint
                    cGust=[cGust obj.constraints.div.gustCases+10^-4];
                end
                
                cGust=[cGust obj.constraints.flutter];
                if ~obj.aeSSM.settings.nonlinActLimits
                    glaFlag=~cellfun(@isempty,{obj.gustCases.glaSurf});
                    glaPresetFlag=~cellfun(@isempty,{obj.gustCases.glaVal});
                    if any(any(and(glaFlag,~glaPresetFlag)))
                        for iCase=1:length(obj.gustCases)
                            cGust=[cGust [obj.constraints.gla(iCase).defMax obj.constraints.gla(iCase).defMin obj.constraints.gla(iCase).rate]];
                        end
                    end
                end
            end                
        end
        
        function [o]=getScaledObjective(obj)
            o=obj.objective*obj.settings.objScaleFactor;
        end
        
        function []=plotConstraints(obj)
            %% init data struct
            nManeuverCases=length(obj.maneuverCases);
            nGustCases=length(obj.gustCases);
            nCases=nManeuverCases+(nGustCases>0);
            
            nBeams=length(obj.settings.beamsToOptimize);
            spFrData=cell(1,nBeams);
            spReData=cell(1,nBeams);
            skLoData=cell(1,nBeams);
            skUpData=cell(1,nBeams);
            spFrBucklData=cell(1,nBeams);
            spReBucklData=cell(1,nBeams);
            skLoBucklData=cell(1,nBeams);
            skUpBucklData=cell(1,nBeams);
            for iBeam=1:nBeams
                beamID=obj.settings.beamsToOptimize(iBeam);
                if nManeuverCases ~= 0 
                    names{iBeam}=cellstr([repmat('Maneuver Case Strength',nManeuverCases,1) num2str([1:nManeuverCases]','%i')]);
                    if ~isempty(obj.gustCases)
                        names{iBeam}=[names{iBeam}; cellstr(['Gust Cases Strength'])];
                    end
                    if ~isempty(obj.str.beam(beamID).ribLocations)
                        names{iBeam}=[names{iBeam}; cellstr([repmat('Maneuver Case Buckl',nManeuverCases,1) num2str([1:nManeuverCases]','%i')])];
                        if ~isempty(obj.gustCases)
                            names{iBeam}=[names{iBeam}; cellstr(['Gust Cases Buckl'])];
                        end
                    end
                else
                    if ~isempty(obj.gustCases)
                        names{iBeam}=[names{iBeam}; cellstr(['Gust Cases Strength'])];
                    end
                    if ~isempty(obj.str.beam(beamID).ribLocations)
                        names{iBeam}=[names{iBeam}; cellstr(['Gust Cases Buckl'])];
                    end
                end
                spFrData{iBeam}=zeros(nCases,obj.str.beam(iBeam).nel);
                spReData{iBeam}=zeros(nCases,obj.str.beam(iBeam).nel);
                skLoData{iBeam}=zeros(nCases,obj.str.beam(iBeam).nel);
                skUpData{iBeam}=zeros(nCases,obj.str.beam(iBeam).nel);
                    if ~isempty(obj.str.beam(beamID).ribLocations)
                        spFrBucklData{iBeam}=zeros(nCases,obj.str.beam(iBeam).nel);
                        spReBucklData{iBeam}=zeros(nCases,obj.str.beam(iBeam).nel);
                        skLoBucklData{iBeam}=zeros(nCases,obj.str.beam(iBeam).nel);
                        skUpBucklData{iBeam}=zeros(nCases,obj.str.beam(iBeam).nel);
                    end
            end

            %% data collection
            for iBeam=1:nBeams
                beamID=obj.settings.beamsToOptimize(iBeam);
                j=0;
                for i=1:length(obj.maneuverCases)
                    j=j+1;
                    spFrData{iBeam}(j,:)=obj.constraints.maneuverCases(i).strength(beamID).spFr;
                    spReData{iBeam}(j,:)=obj.constraints.maneuverCases(i).strength(beamID).spRe;
                    skLoData{iBeam}(j,:)=obj.constraints.maneuverCases(i).strength(beamID).skLo;
                    skUpData{iBeam}(j,:)=obj.constraints.maneuverCases(i).strength(beamID).skUp;
                    if ~isempty(obj.str.beam(beamID).ribLocations)
                        spFrBucklData{iBeam}(j,:)=obj.constraints.maneuverCases(i).buckl(beamID).spFr;
                        spReBucklData{iBeam}(j,:)=obj.constraints.maneuverCases(i).buckl(beamID).spRe;
                        skLoBucklData{iBeam}(j,:)=obj.constraints.maneuverCases(i).buckl(beamID).skLo;
                        skUpBucklData{iBeam}(j,:)=obj.constraints.maneuverCases(i).buckl(beamID).skUp;
                    end
                        
                end

                j=j+1;
                if ~isempty(obj.gustCases)
                    spFrData{iBeam}(j,:)=obj.constraints.gustCases.strength(beamID).spFr;
                    spReData{iBeam}(j,:)=obj.constraints.gustCases.strength(beamID).spRe;
                    skLoData{iBeam}(j,:)=obj.constraints.gustCases.strength(beamID).skLo;
                    skUpData{iBeam}(j,:)=obj.constraints.gustCases.strength(beamID).skUp;
                    if ~isempty(obj.str.beam(beamID).ribLocations)
                            spFrBucklData{iBeam}(j,:)=obj.constraints.gustCases.buckl(beamID).spFr;
                            spReBucklData{iBeam}(j,:)=obj.constraints.gustCases.buckl(beamID).spRe;
                            skLoBucklData{iBeam}(j,:)=obj.constraints.gustCases.buckl(beamID).skLo;
                            skUpBucklData{iBeam}(j,:)=obj.constraints.gustCases.buckl(beamID).skUp;
                    end
                end
                
            end
            %% data normalization

            %% data plotting
            for iBeam=1:nBeams
                beamID=obj.settings.beamsToOptimize(iBeam);
                nEl=obj.str.beam(beamID).nel;
                if obj.str.beam(beamID).anisotropic
                    sparIdxMinThick=[obj.settings.thicknessControlStations{iBeam}((round(obj.desVar.Str(iBeam).tSpFr,6)==obj.lowerBounds.Str(iBeam).tSpFr))];
                    skinIdxMinThick=[obj.settings.thicknessControlStations{iBeam}((round(obj.desVar.Str(iBeam).tSkUp,6)==obj.lowerBounds.Str(iBeam).tSkUp))];
                else
                    sparIdxMinThick=obj.settings.thicknessControlStations{iBeam}((round(obj.desVar.Str(iBeam).tSp,6)==obj.lowerBounds.Str(iBeam).tSp));
                    skinIdxMinThick=obj.settings.thicknessControlStations{iBeam}((round(obj.desVar.Str(iBeam).tSk,6)==obj.lowerBounds.Str(iBeam).tSk));
                end
                allSparIdxMinThick=[(nEl/2)+1-sparIdxMinThick nEl/2+sparIdxMinThick ];
                allSkinIdxMinThick=[(nEl/2)+1-skinIdxMinThick nEl/2+skinIdxMinThick ];

                figure;
                subplot(2,2,1)
                grid on; hold on; box on;
                title('Front Spar');
                plot(spFrData{iBeam}');
                plot(spFrBucklData{iBeam}','--');
                plot(allSparIdxMinThick, zeros(1,length(allSparIdxMinThick)),'ro');
                lim=axis;
                plot(repmat([(nEl/2+1)-obj.settings.thicknessControlStations{iBeam} nEl/2+obj.settings.thicknessControlStations{iBeam}],2,1), repmat([-1 1], length(obj.settings.thicknessControlStations{iBeam})*2,1)','--', 'color', [0.1 0.1 0.1])
                axis(lim);

                subplot(2,2,2)
                grid on; hold on; box on;
                title('Rear Spar');
                plot(spReData{iBeam}');
                plot(spReBucklData{iBeam}','--');
                plot(allSparIdxMinThick, zeros(1,length(allSparIdxMinThick)),'ro');
                lim=axis;
                plot(repmat([(nEl/2+1)-obj.settings.thicknessControlStations{iBeam} nEl/2+obj.settings.thicknessControlStations{iBeam}],2,1), repmat([-1 1], length(obj.settings.thicknessControlStations{iBeam})*2,1)','--', 'color', [0.1 0.1 0.1])
                axis(lim);

                subplot(2,2,3)
                grid on; hold on; box on;
                title('Upper Skin');
                plot(skUpData{iBeam}');
                plot(skUpBucklData{iBeam}','--');
                plot(allSkinIdxMinThick, zeros(1,length(allSkinIdxMinThick)),'ro');
                lim=axis;
                plot(repmat([(nEl/2+1)-obj.settings.thicknessControlStations{iBeam} nEl/2+obj.settings.thicknessControlStations{iBeam}],2,1), repmat([-1 1], length(obj.settings.thicknessControlStations{iBeam})*2,1)','--', 'color', [0.1 0.1 0.1])
                axis(lim);

                subplot(2,2,4)
                grid on; hold on; box on;
                title('Lower Skin');
                plot(skLoData{iBeam}');
                plot(skLoBucklData{iBeam}','--');
                plot(allSkinIdxMinThick, zeros(1,length(allSkinIdxMinThick)),'ro');
                lim=axis;
                plot(repmat([(nEl/2+1)-obj.settings.thicknessControlStations{iBeam} nEl/2+obj.settings.thicknessControlStations{iBeam}],2,1), repmat([-1 1], length(obj.settings.thicknessControlStations{iBeam})*2,1)','--', 'color', [0.1 0.1 0.1])
                axis(lim);
                legend([names{iBeam}; 'Minimum Thickness'; 'Str Control Stations']);
                
            end
        end
        
        function obj=computeGradients(obj,stepVec)
            
            %required constants
            nomConst=obj.getConstraints;
            nomObj=obj.getScaledObjective;
            nomDesVar=obj.getScaledDesignVariables;
            
            nConstr=length(nomConst);
            nDesVar=length(nomDesVar);
            
            %allocation
            constGradCalc=zeros(nDesVar,nConstr);
            objGradCalc=zeros(nDesVar,1);
            
            desVarMat=diag(stepVec)+repmat(nomDesVar,nDesVar,1);
            desVarMatBack=-diag(stepVec)+repmat(nomDesVar,nDesVar,1);
            %create a matrix with the flags for str man and gust nom/dyn
            %computation; create vector containing the man/gust case id
            %to be calculated (set zero when all cases have to be calculated)
            analysisType=zeros(nDesVar,4);
            caseId=zeros(nDesVar,1);
            for i=1:nDesVar
                splittedName=strsplit(obj.desVarTypes{i},'_');
                if strcmp(splittedName{1},'Str')
                    analysisType(i,:)=1; %all analyses have to be evaluated when the structure is modified
                elseif strcmp(splittedName{1},'maneuverCases')
                    caseId(i)=str2double(splittedName{2});
                    analysisType(i,:)=[0 1 0 0]; % only maneuvers have to be evaluated 
                elseif strcmp(splittedName{1},'gustCases')
                    caseId(i)=str2double(splittedName{2});
                    if strcmp(splittedName{3},'pta')
                        analysisType(i,:)=[0 0 1 0]; % only nominal case has to be evaluated
                    elseif strcmp(splittedName{3},'gla')
                        analysisType(i,:)=[0 0 0 1]; % only dynamic part has to be evaluated
                    end
                end
            end
            %parallel computation
            tic
            parfor i=1:nDesVar
                optimCasePerturbation=obj;
                if strcmp(obj.fdType,'backward')
                    optimCasePerturbation=optimCasePerturbation.setScaledDesignVariables(desVarMat(i,:));
                else
                    optimCasePerturbation=optimCasePerturbation.setScaledDesignVariables(desVarMat(i,:));
                end
                optimCasePerturbation=optimCasePerturbation.evaluate(analysisType(i,:),caseId(i),1);
                constPert=optimCasePerturbation.getConstraints;
                objPert=optimCasePerturbation.getScaledObjective;
                if strcmp(obj.fdType,'central')
                    optimCasePerturbation=obj;
                    optimCasePerturbation=optimCasePerturbation.setScaledDesignVariables(desVarMatBack(i,:));
                    optimCasePerturbation=optimCasePerturbation.evaluate(analysisType(i,:),caseId(i),1);
                    constPertBack=optimCasePerturbation.getConstraints;
                    objPertBack=optimCasePerturbation.getScaledObjective;
                    
                    constGradCalc(i,:)=(constPert-constPertBack)/stepVec(i)/2;
                    objGradCalc(i,1)=(objPert-objPertBack)/stepVec(i)/2;
                elseif strcmp(obj.fdType,'backward')
                    constGradCalc(i,:)=(constPert-nomConst)/-stepVec(i);
                    objGradCalc(i,1)=(objPert-nomObj)/-stepVec(i);
                else
                    constGradCalc(i,:)=(constPert-nomConst)/stepVec(i);
                    objGradCalc(i,1)=(objPert-nomObj)/stepVec(i);
                end
            end
            t2=toc;
            fprintf(['Total ' obj.fdType ' gradient computation time: %.2f'],t2);
            obj.constGrad=constGradCalc;
            obj.objGrad=objGradCalc;
            obj.gradStepVec=stepVec;
        end
        
        function obj=savePointInHistory(obj)
            if isempty(obj.history)
                % initialize history struct
                obj.history.constraints=[];
                obj.history.desVar=[];
                obj.history.objGrad=[];
                obj.history.constGrad=[];
                obj.history.objective=[];
                obj.history.gradStepVec=[];
            end
            obj.history.constraints=[obj.history.constraints; obj.getConstraints];
            obj.history.desVar=[obj.history.desVar; obj.getScaledDesignVariables];
            obj.history.objGrad=cat(3,obj.history.objGrad, obj.objGrad);
            obj.history.constGrad=cat(3,obj.history.constGrad, obj.constGrad);
            obj.history.objective=[obj.history.objective; obj.getScaledObjective];
            obj.history.gradStepVec=[obj.history.gradStepVec; obj.gradStepVec];
        end
            
        function obj=setTailVal(obj,controlElementsStiffDir,tailVal)
            for iBeam=1:length(obj.settings.beamsToOptimize)
                beamID=obj.settings.beamsToOptimize(iBeam);
                if obj.str.beam(beamID).is_sym
                    lengths=[obj.str.beam(beamID).beamelement(obj.str.beam(beamID).nel/2+1:end).le]';
                else
                    lengths=[obj.str.beam(beamID).beamelement(1:end).le]';
                end
                   
                lengthStations=[0; cumsum(lengths(1:end-1))]+lengths/2;
                
                

                nCtrStations=length(controlElementsStiffDir);
                lowerSkinAnglesVec=tailVal(1:nCtrStations);
                upperSkinAnglesVec=tailVal(nCtrStations+1:2*nCtrStations);
                frontSparAnglesVec=tailVal(nCtrStations*2+1:3*nCtrStations);
                rearSparAnglesVec=tailVal(nCtrStations*3+1:4*nCtrStations);

                lowerSkinAnglesHalfWing=interp1(lengthStations(controlElementsStiffDir),lowerSkinAnglesVec,lengthStations,'next','extrap');
                upperSkinAnglesHalfWing=interp1(lengthStations(controlElementsStiffDir),upperSkinAnglesVec,lengthStations,'next','extrap');
                frontSparAnglesHalfWing=interp1(lengthStations(controlElementsStiffDir),frontSparAnglesVec,lengthStations,'next','extrap');
                rearSparAnglesHalfWing=interp1(lengthStations(controlElementsStiffDir),rearSparAnglesVec,lengthStations,'next','extrap');
                if obj.str.beam(beamID).is_sym
                    allLowerSkinAngles=[lowerSkinAnglesHalfWing(end:-1:1); -lowerSkinAnglesHalfWing];
                    allUpperSkinAngles=[upperSkinAnglesHalfWing(end:-1:1); -upperSkinAnglesHalfWing];
                    allFrontSparAngles=[frontSparAnglesHalfWing(end:-1:1); -frontSparAnglesHalfWing];
                    allRearSparAngles=[rearSparAnglesHalfWing(end:-1:1); -rearSparAnglesHalfWing];
                else
                    allLowerSkinAngles=[lowerSkinAnglesHalfWing];
                    allUpperSkinAngles=[upperSkinAnglesHalfWing];
                    allFrontSparAngles=[frontSparAnglesHalfWing];
                    allRearSparAngles=[rearSparAnglesHalfWing];
                end
                for iEl = 1:obj.str.beam(beamID).nel
                    obj.str.beam(beamID).beamelement(iEl).crosssection.laminate_sk_lo.offset_angle=allLowerSkinAngles(iEl);
                    obj.str.beam(beamID).beamelement(iEl).crosssection.laminate_sk_up.offset_angle=allUpperSkinAngles(iEl);
                    obj.str.beam(beamID).beamelement(iEl).crosssection.laminate_fs.offset_angle=allFrontSparAngles(iEl);
                    obj.str.beam(beamID).beamelement(iEl).crosssection.laminate_rs.offset_angle=allRearSparAngles(iEl);
                end    
            end
        end
        
    end
    
end

