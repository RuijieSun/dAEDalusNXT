function [x,f,eflag,outpt, myOptimCase] = runOptimizationCase(myOptimCase,opts)

myOptimCase.fdType=opts.FiniteDifferenceType;
xLast = []; % Last x, evaluate was caled
xLastGrad =[]; %Last x, gradients were computed 
myObjective = []; % Use for objective at xLast
myConstraint = []; % Use for nonlinear inequality constraint
myObjGrad = []; % Use for objective at xLast
myConstGrad = []; % Use for nonlinear inequality constraint

initFdStepSize=opts.FiniteDifferenceStepSize;
minFdStepSize=opts.FiniteDifferenceStepSize;
opts.OutputFcn=@outfun;
opts.SpecifyObjectiveGradient=true;
opts.SpecifyConstraintGradient=true;
%the current state of the optimizer(function evaluation or finite
%difference evaluation) is passed to the objective and constraint function
%by evaluating the expression 'exist('funfcn')' within fmincon. Only if it
%is an function evaluation, the current state of myOptimCase is updated.
%The update is required to have faster trimming and store the current
%deformation case for all static cases
fun = @(x)objfun( x ); % the objective function, nested below
cfun = @(x)constr( x); % the constraint function, nested below

% Call fmincon
x0=myOptimCase.getScaledDesignVariables  ;
[lb, ub]=myOptimCase.getScaledBounds;
[x,f,eflag,outpt] = fmincon(fun,x0,[],[],[],[],lb,ub,cfun,opts);

    function [y,dy] = objfun( x)
        if ~isequal(x,xLast) % Check if computation is necessary
            myOptimCase=myOptimCase.setScaledDesignVariables(x);
            myOptimCase=myOptimCase.evaluate();
            myConstraint=myOptimCase.getConstraints;
            myObjective=myOptimCase.getScaledObjective;
            xLast = x;
        end
        % Now compute objective function
        y = myObjective;
        if nargout>1
            if ~isequal(x,xLastGrad)
                if ~isempty(xLastGrad)
                    stepVec=(x-xLastGrad);
                    stepVec(stepVec==0)=minFdStepSize; % repair zeros
                    idxToSmall=abs(stepVec)<minFdStepSize; %check step size
                    stepVec(idxToSmall)=sign(stepVec(idxToSmall))*minFdStepSize; % increase step size of too small steps
                    stepVec=ones(1,length(x))*minFdStepSize;
                    %calculate fd step direction based on change in last 5
                    %iterations; if change in last 5 is less than
                    %minFdStepSize, then go positive
                    nItMean=min(5, size(myOptimCase.history.desVar,1));
                    change=x-myOptimCase.history.desVar(nItMean,:);
                    change(abs(change)<minFdStepSize)=minFdStepSize;
                    stepDir=sign(change);
                    stepVec=stepDir.*ones(1,length(x))*minFdStepSize;
                else
                    stepVec=ones(1,length(x))*minFdStepSize;
                end
                % check if bounds are satisfied within finite differences
                % upper bound violated?
                idxUb=find((x+stepVec)>ub);
                for iUb=1:length(idxUb)
                % add current x
                    idxDesVar=idxUb(iUb);
                    if (ub(idxDesVar)-x(idxDesVar))>minFdStepSize
                        %if margin between x and upper bound is larger than 
                        %minimum step size, then set finite difference step 
                        %size so that we end up on the upper bound
                        stepVec(idxDesVar)=ub(idxDesVar)-x(idxDesVar);  
                    else
                        stepVec(idxDesVar)=-stepVec(idxDesVar);
                    end
                end
                % check if smaller than lower bound
                idxLb=find((x+stepVec)<lb);
                for iLb=1:length(idxLb)
                % add current x
                    idxDesVar=idxLb(iLb);
                    if (x(idxDesVar)-lb(idxDesVar))>minFdStepSize 
                        %if margin between x and lower bound is larger than 
                        %minimum step size, then set finite difference step 
                        %size so that we end up on the lower bound
                        stepVec(idxDesVar)=x(idxDesVar)-lb(idxDesVar); 
                    else %invert the step direction
                        stepVec(idxDesVar)=-stepVec(idxDesVar);
                    end
                end
                
                myOptimCase=myOptimCase.computeGradients(stepVec);
                myConstGrad=myOptimCase.constGrad;
                myObjGrad=myOptimCase.objGrad;
                xLastGrad=x;
            end
            dy=myObjGrad;
        end
    end

    function [c, ceq, dc, dceq] = constr( x)
        if ~isequal(x,xLast) % Check if computation is necessary
            myOptimCase=myOptimCase.setScaledDesignVariables(x);
            myOptimCase=myOptimCase.evaluate();
            myConstraint=myOptimCase.getConstraints;
            myObjective=myOptimCase.getScaledObjective;
            xLast = x;
        end
        % Now compute constraint functions
        c = myConstraint; % In this case, the computation is trivial
        ceq=[];
        if nargout>2
            if ~isequal(x,xLastGrad)
                if ~isempty(xLastGrad)
                    stepVec=(x-xLastGrad);
                    stepVec(stepVec==0)=minFdStepSize; % repair zeros
                    idxToSmall=abs(stepVec)<minFdStepSize; %check step size
                    stepVec(idxToSmall)=sign(stepVec(idxToSmall))*minFdStepSize; % increase step size of too small steps
                    stepVec=ones(1,length(x))*minFdStepSize;
                    %calculate fd step direction based on overall change of
                    %design variable; if change in last 5 is less than
                    %minFdStepSizea, then go forward
                    nItMean=min(5, size(myOptimCase.history.desVar,1));
                    change=x-myOptimCase.history.desVar(nItMean,:);
                    change(abs(change)<minFdStepSize)=minFdStepSize;
                    stepDir=sign(change);
                    stepVec=stepDir.*ones(1,length(x))*minFdStepSize;
                else
                    stepVec=ones(1,length(x))*initFdStepSize;
                end
                % check if bounds are satisfied within finite differences
                % upper bound violated?
                idxUb=find((x+stepVec)>ub);
                for iUb=1:length(idxUb)
                % add current x
                    idxDesVar=idxUb(iUb);
                    if (ub(idxDesVar)-x(idxDesVar))>minFdStepSize
                        %if margin between x and upper bound is larger than 
                        %minimum step size, then set finite difference step 
                        %size so that we end up on the upper bound
                        stepVec(idxDesVar)=ub(idxDesVar)-x(idxDesVar);  
                    else
                        stepVec(idxDesVar)=-stepVec(idxDesVar);
                    end
                end
                % check if smaller than lower bound
                idxLb=find((x+stepVec)<lb);
                for iLb=1:length(idxLb)
                % add current x
                    idxDesVar=idxLb(iLb);
                    if (x(idxDesVar)-lb(idxDesVar))>minFdStepSize 
                        %if margin between x and lower bound is larger than 
                        %minimum step size, then set finite difference step 
                        %size so that we end up on the lower bound
                        stepVec(idxDesVar)=x(idxDesVar)-lb(idxDesVar); 
                    else %invert the step direction
                        stepVec(idxDesVar)=-stepVec(idxDesVar);
                    end
                end
                



                myOptimCase=myOptimCase.computeGradients(stepVec);
                myConstGrad=myOptimCase.constGrad;
                myObjGrad=myOptimCase.objGrad;
                xLastGrad=x;
            end
            dc=myConstGrad;
            dceq=[];
        end
    end

    function stop = outfun(x,optimValues,state)
        stop=0;
        switch state
            case 'iter'
                  % Make updates to plot or guis as needed
                  figure(99)
                  subplot(2,2,1)
                  plot(optimValues.iteration,optimValues.fval,'bo');
                  subplot(2,2,2)
                  imagesc([optimValues.iteration optimValues.iteration],[1 size(x)],x')
                  subplot(2,2,3)
                  imagesc([optimValues.iteration optimValues.iteration],[1 size(myConstraint)],myConstraint')
                  subplot(2,2,4)
                  imagesc([optimValues.iteration optimValues.iteration],[1 size(optimValues.gradient)],optimValues.gradient)
                  pause(0.01);
                  myOptimCase=myOptimCase.savePointInHistory;
                  if isdeployed
                      diary off
                      diary on
                  end
                  %saving intermediate results
                  delete results_currentResults*
                  delete optimCase_currentResults*
                  saveResults(myOptimCase,['currentResults' num2str(optimValues.iteration)])
                  
            case 'interrupt'
                  % Probably no action here. Check conditions to see  
                  % whether optimization should quit.
            case 'init'
                  % Setup for plots or guis
                  figure(99)
                  close 99
                  figure(99)
                  subplot(2,2,1)
                  title('Objective');
                  hold on; box on;
                  subplot(2,2,2)
                  title('DesignVariables');
                  hold on; box on;
                  subplot(2,2,3)
                  title('Constraints');
                  hold on; box on;
                  subplot(2,2,4)
                  title('Gradient');
                  hold on; box on;
            case 'done'
                  % Cleanup of plots, guis, or final plot
        otherwise
        end
    end
end