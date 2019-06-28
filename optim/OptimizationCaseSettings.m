classdef OptimizationCaseSettings
    %OPTIMIZATIONCASESETTINGS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        %array of id of beams which are changed during the optimization
        beamsToOptimize=[];
        %array of control stations (same size as beamsToOptimize) for
        %thickness
        thicknessControlStations={};
        %array of control stations (same size as beamsToOptimize) for
        %stiffness direction (only used in anisotropic models)
        stiffnessDirectionControlStations={};
        % order of highest chebyshev polynomial when tailoring is done by
        % shape functions
        stiffnessDirectionPolynomialOrder={};
        %aeroelastic solver settings
        aeroelasticSolverSettings@class_aeroelastic_solver_settings
        %scale factor for str design variables
        thicknessScaleFactor=15;
        %scale factor for skin and spar design variables
        thicknessScaleFactorSkin=[];
        thicknessScaleFactorSpar=[];
        %scale factor for mla design variables
        mlaScaleFactor=0.01;
        %scale factor for mla design variables
        rollAllocScaleFactor=1;
        %scale factor for gla design variables
        glaScaleFactor=0.1;
        %scale factor for stiffness direction offset angle
        stiffnessDirectionScaleFactor=0.01;
        %bounds for stiffness direction in degree
        stiffnessDirectionBounds=[-Inf Inf];
        %bounds for stiffness direction by polynomial shape functions
        stiffnessDirectionPolBounds=[-Inf Inf];
        %scale factor for objective 
        objScaleFactor=10^-5;
        %clamped wing (special case)
        clampedWing=0;
        % flag for closed loop gust analysis
        closedLoopGust=0;
        %divergence constraint active
        divConstraint=1;
        % number of gust lengths tested each gust case
        nGustLengths=15;
        % gla controller order per control surface
        glaControllerOrder=2^4;
        % flag for act/deactivating spar tailoring desVar
        sparTailoring=1;
    end
    
    methods
    end
    
end

