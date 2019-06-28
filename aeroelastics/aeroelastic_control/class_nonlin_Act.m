classdef class_nonlin_Act
    %CLASS_NONLIN_ACT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        inputName
        outputName
        maxRate=80; %degree per second
        ub=[] % upper deflection bound in degrees
        lb=[] % lower deflection bound in degrees
        DampCoeff=0.7;
        w0=10;
        asymFlag        
    end
    
    methods
        function obj=class_nonlin_Act(inputName, asymFlag)
            obj.inputName=inputName;
            obj.outputName=inputName;
            obj.asymFlag=asymFlag;
        end
        
        function data=simulate(obj,signal,tVec)
            data=simulateLimited2ndOrderActDiscrete(signal,tVec,obj.w0,obj.DampCoeff,obj.maxRate*pi/180,obj.lb*pi/180,obj.ub*pi/180);
            data=[data; data];
            data=data([1 4 2 5 3 6],:);
            if ~obj.asymFlag
                % invert signal for other side (deflected in same direction means opposite sign )
                data([2 4 6],:)=-data([2 4 6],:);
            end
        end
    end
    
end

