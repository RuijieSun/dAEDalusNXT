%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%
function [ sysR ] = seriesPreserve( sys1, sys2,connectionGroups )
%SERIESPRESERVE connects 2 models and preserves inputs of second model which
%are not in the connectionGroups
conIds1=[];
conIds2=[];
for iG=1:length(connectionGroups)
    conIds1=[conIds1 getfield(sys1.OutputGroup,connectionGroups{iG})];
    conIds2=[conIds2 getfield(sys2.InputGroup,connectionGroups{iG})];
end
minIndex2=min(conIds2);
maxIndex2=max(conIds2);
% check if there is inputs of system 2 in between the connectiongroups:
if length(minIndex2:maxIndex2)>length(conIds2)
    % resort model 2 so that all connectioninputs are at the end
    sortIn=[setdiff(1:length(sys2.InputName),conIds2) conIds2];
    sys2=sys2(:,sortIn);
    %recompute conIds2
    conIds2=[];
    for iG=1:length(connectionGroups)
        conIds2=[conIds2 getfield(sys2.InputGroup,connectionGroups{iG})];
    end
    minIndex2=min(conIds2);
    maxIndex2=max(conIds2);
end
minIndex1=min(conIds1);
maxIndex1=max(conIds1);
% check if there is outputs of system 1 in between the connectiongroups:
if length(minIndex2:maxIndex2)>length(conIds2)
    disp('STOP HERE, the first model has outputs between the selected connection outputs, this does not yet work with seriesPreserve')
    return
end


sysPrv1=ss(eye(length(sys2.InputName)));
sysPrv1.OutputName=sys2.InputName;
sysPrv1.InputName=sys2.InputName;
sysPrv1.InputDelay=sys2.InputDelay;
sysPrv1.InputGroup=sys2.InputGroup;
sysPrv1.OutputGroup=sys2.InputGroup;
sysPost1=sysPrv1;

sysPrv1=sysPrv1(1:minIndex2-1,1:minIndex2-1);
sysPost1=sysPost1(maxIndex2+1:end,maxIndex2+1:end);

sys1b=append(sysPrv1,sys1,sysPost1);

sysPrv2=ss(eye(length(sys1.OutputName)));
sysPrv2.OutputName=sys1.OutputName;
sysPrv2.InputName=sys1.OutputName;
sysPrv2.InputGroup=sys1.OutputGroup;
sysPrv2.OutputGroup=sys1.OutputGroup;
sysPost2=sysPrv2;

sysPrv2=sysPrv2(1:minIndex1-1,1:minIndex1-1);
sysPost2=sysPost2(maxIndex1+1:end,maxIndex1+1:end);

sys2b=append(sysPrv2,sys2,sysPost2);

sysR=series(sys1b,sys2b,'name');
end

