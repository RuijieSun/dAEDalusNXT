%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%
function [ redSys,redBase ] = bpod( inSys, nOrder, bpodT, nPODModes)
%BPOD Summary of this function goes here
%   Detailed explanation goes here
stepsimtime=0:bpodT(1):bpodT(2); % only 6 Timesteps


%% empirical controllability
%generate snapshot data
fullSys=ss(inSys.a,inSys.b,ones(1,size(inSys.a,1)),0);
inputs=find(any(fullSys.b));  %<- only use inputs which have an influence, not sure if that works
[~, impT,impX]=impulse(fullSys(:,inputs),stepsimtime);  %<- improve by using sss toolbox?

% reorder data
impTStep=impT(3)-impT(2);
impXper=permute(impX,[2 1 3]);
impXper=impXper(1:order(fullSys),:,:);
xMat=reshape(impXper,order(fullSys),length(impT)*size(impX,3));

% multiply with quadrature coefficient
xMat=xMat.*sqrt(impTStep);

%estimate controllability 
wcemp=xMat*xMat';

%% get POD Modes
% only do when nPODModes is specified, only then the output projection is
% required
if ~isempty(nPODModes)
    xMatPodAllt=reshape(impXper,order(fullSys),length(impT)*size(impX,3));
    xMatPod=inSys.c*xMatPodAllt(:,1:end);
    [Vr,~,~]=svd(xMatPod);
    % reduction base
    nPODModes=size(inSys,1); %<- when using as many PODModes as inputs in the system, and the number of inputs is larger than the number of outputs, then there is no need for output projection!
    PodRed=Vr(:,1:nPODModes);      %add criteria for determination of the projection error by using the sum of eigenvalues (see paper)
    % PodRed=Vr*Vr';
end

%% empirical observability

% !! the old way of determining the x history with this.....
%   
%   fullSys_adj=ss(fullSys.a',fullSys.c'*PodRed,fullSys.b',0);
%   fullSys_adj=ss(fullSys.a',fullSys.c',fullSys.b',0);
%   inputs=find(any(fullSys_adj.b));
%   [~, impAdjT,impAdjX]=impulse(fullSys_adj(1,inputs),stepsimtime);
%
%.... was not good. Now add "state" outputs to adjoint system

if ~isempty(nPODModes)
    fullSys_adj=ss(inSys.a',inSys.c'*PodRed,eye(size(inSys.a)),0);
else
    fullSys_adj=ss(inSys.a',inSys.c',eye(size(inSys.a)),0);    
end

inputs=find(any(fullSys_adj.b));
[impAdjX, impAdjT]=impulse(fullSys_adj(:,inputs),stepsimtime);

% reorder data
impAdjTStep=impAdjT(3)-impAdjT(2);
impAdjXper=permute(impAdjX,[2 1 3]);
yMatAllt=reshape(impAdjXper,order(inSys),length(impAdjT)*size(impAdjX,3));
yMat=yMatAllt(:,1:end);
% multiply with quadrature coefficient
yMat=yMat.*sqrt(impAdjTStep);
% estimate observability
woemp=yMat*yMat';
%% balancing
% build transformation matrix
% option a
[redBaseF, redBaseVal]=eig(woemp*wcemp);
nOrderRec=rank(woemp*wcemp);
redBaseOnlyReal=redBaseF(:,find(imag(diag(redBaseVal))==0));
%how to find the best reduction order?
nOrderMax=size(redBaseOnlyReal,2);
if nOrder>nOrderMax
    fprintf('warning, specified order of %i too large, maximum order %i\n',nOrder, nOrderMax)
elseif isempty(nOrder)
    nOrder=nOrderMax;
    fprintf('warning, no order specified, using maximum order of %i, recommended Order %i\n',nOrder, nOrderRec)
end
redBase=real(redBaseOnlyReal(:,1:nOrder));     %<- workaround for complex eigenvalues -> also not sure if that works
% option b
% [a,b,c]=svd(yMat'*xMat,'econ');
% T=(xMat*c*b^-0.5);
% S=(b^-0.5*a'*yMat');
% redBase=T;


% actual reduction
redSys=ss(pinv(redBase)*inSys.a*redBase,  pinv(redBase)*inSys.b, inSys.c*redBase, inSys.d);
redSys.inputName=inSys.inputName;
redSys.outputName=inSys.outputName;
end

