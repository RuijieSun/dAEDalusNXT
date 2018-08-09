%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%
classdef AeSSMSettings
    %AESSMSETTINGS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        %selected modes
        modes = 1:10
        %model Order Reduction method (currently only bPod)
        mor = 'bPod'
        %reduction order
        morOrder=50;
        %cheb order for cheb type bt
        morChebOrder=[20 5];
        %time settings for bod (tstep ttot)
        bpodT=[];
        %maximum error by output projection
        errPODProj=1e-3;
        %clamped or free AeSSM
        restrained = 0
        %points of interest (DOFs of free structure)
        poiDof=[1 2 3 4 5 6]
        %flag for gustInputs (gustZones specified in class_aircraft;)
        gustInputs =1
        %flag for rbmInputs
        rbmInputs =1
        %flag for controlInputs
        ctrInputs =1
        %flag for panel velocity inputs
        vPInputs =0
        % array of load stations for local load outputs
        strLoadStations=[1:6]
        %flag for deformed grid
        deformedGrid=1
        %rayleigh damping factor
        strRayleighDamping=0.00
        %proportional damping factor
        strPropDamping=0.00
        %exponential growth ratio for wake spacing (1 for linear)
        rW=1.15
        %additional wake length factor for last panel (last panel gets
        %stretched additionally to rW
        addLengthFactorLastPanel=20;
        %fixed wake panels (those rows are not stretched by rW)
        nFixedWakePanels=10;
        %wake length factor (drives number of states)
        wakeLengthFactor=1
        %wake panel density (drives number of states)
        wakePanelDensity=1
        %force type
        force='all';
        %consider induced component
        indComp=0;
        % aoa probe position relative to first gust zone
        aoaProbePos=-5; % in meters
        % contTurbSignal (uPCP and uSCP are stored once created in order to
        % enable optimization without randomness)
        contTurbSignal
    end
    
    methods
    end
    
end

