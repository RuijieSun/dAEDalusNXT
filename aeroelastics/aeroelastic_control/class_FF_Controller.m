%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%
classdef class_FF_Controller
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        % high pass characteristics of FF law (avoids steady aoa being damped)
        highPass = tf([1 0],[1 4]);
        % Delay of the feed forward controller
        delay = tf([1 -150 7500],[1 150 7500]);
        % order of the controller
        order=2^4;
        % ff coefficients
        coeff
        % outputnames
        outputName
        % ssm representation
        contSSM
        %timeStep
        Ts
    end
    
    methods
        function obj=class_FF_Controller(varargin)
            %check if order is equal to length(x)/length(inputName) and
            %length(x)/length(outputname)
            if length(varargin)==3
                outputName=varargin{1};
                order=varargin{2};
                coeff=varargin{3};
            elseif length(varargin)==2
                outputName=varargin{1};
                order=varargin{2};
                coeff=zeros(1,order*length(outputName));
            elseif length(varargin)==1
                outputName=varargin{1};
                order=obj.order;
                coeff=zeros(1,order*length(outputName));
            end
                
            if ~(length(coeff)/length(outputName)==order)
                disp('input error; length of inputName, n and x needs to be consistent')
            end
            obj.coeff=coeff;
            obj.order=order;
            obj.outputName=outputName;            
        end
        function obj=genSSM(obj,Ts)
            obj.Ts=Ts;
            FFController = tf(zeros(length(obj.outputName),1));
            denom=[1 zeros(1,obj.order-1)];
            nCs=length(obj.outputName);
            for iCS=1:length(obj.outputName)
                FFController(iCS)=tf(obj.coeff((iCS-1)*obj.order+1:iCS*obj.order),              denom,Ts);
            end
            
            opt = d2cOptions('Method','tustin','PrewarpFrequency',0.5);
            obj.contSSM = d2c(ss(FFController),opt)*obj.highPass*obj.delay;
            obj.contSSM.InputName = {'aoaProbeIn'};
            obj.contSSM.OutputName = obj.outputName;
            obj.contSSM.StateName = cellstr([repmat('ffCtr_',length(obj.contSSM.StateName),1) num2str([1:length(obj.contSSM.StateName)]','%04d')]);
        end
        function obj=bode(obj)
            N=obj.order;
            nCs=length(obj.outputName);
            figure, hold on;
            for iCs=1:nCs
                bode(tf([obj.coeff((iCs-1)+1:iCs*N)],[1 zeros(1,obj.order-1)],obj.Ts));
            end
            legend(obj.outputName);
            
        end
        function [varargout]=getCommands(obj,inputSignals,timeVec)
            for iSim=1:size(inputSignals,2)
                simRes(:,:,iSim)=lsim(obj.contSSM,inputSignals(:,iSim),timeVec);
            end
            if nargout==1
                varargout{1}=simRes;
            else
                figure
                for iCs=1:size(simRes,2)
                    subplot(size(simRes,2),2,(iCs-1)*2+1)
                    plot(timeVec,squeeze(simRes(:,iCs,:)))
                    ylabel(obj.outputName(iCs));
                    xlabel('time')
                    
                    subplot(size(simRes,2),2,(iCs-1)*2+2)
                    plot(timeVec(1:end-1),diff(squeeze(simRes(:,iCs,:)))/(timeVec(2)-timeVec(1)))
                    ylabel([obj.outputName(iCs) '_Dot']);
                    xlabel('time')
                    
                end
                suptitle('Commanded Def and Approx. Rate')
            end
                
        end
        
    end
    
end
