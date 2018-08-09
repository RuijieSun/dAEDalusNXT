%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%
function [ redBase, redBaseInv,r, S] = bpod( inSys, nOrder, bpodT, errPODProj)
%BPOD Summary of this function goes here
%   Detailed explanation goes here
tic;
stepsimtime=0:bpodT(1):bpodT(2); % only 6 Timesteps

%% empirical controllability
%generate snapshot data
fullSys=ss(inSys.a,inSys.b,ones(1,size(inSys.a,1)),0);
inputs=find(any(fullSys.b));  %<- only use inputs which have an influence, not sure if that works
timeVecDivisions=15;
nStepPerDiv=10;
incrFact=0.5; % time step increases by 50 percent each division
idx=unique(round(logspace(0,log10(501),250)));
t0=bpodT(1);
sumStepSimTimeDiv=0;
t1=t0;
for iDivision=1:timeVecDivisions
    currTimestep=t0*(1+(iDivision-1)*incrFact);
    stepsimtimeDiv{iDivision}=t1:currTimestep:t1+(nStepPerDiv-1)*currTimestep;
    t1=t1+nStepPerDiv*currTimestep;
    sumStepSimTimeDiv=sumStepSimTimeDiv+length(stepsimtimeDiv{iDivision});
end
if size(inSys.a,1)*length(stepsimtime)*size(inSys.b,2)>5000*10^9
    impT=zeros(sumStepSimTimeDiv,1);
    impX=zeros(sumStepSimTimeDiv,size(inSys.a,1),size(inSys.b,2));
    nPast=0;
    for iDivision=1:timeVecDivisions
        nCurr=length(stepsimtimeDiv{iDivision});
        [~, impTcurr,impXcurr]=impulse(fullSys(:,inputs),stepsimtimeDiv{iDivision});
        impT(nPast+1:nPast+nCurr,:)=[ impTcurr];
        impX(nPast+1:nPast+nCurr,:,:)=[ impXcurr];
        nPast=nCurr+nPast;
    end
else
    [~, impT,impX]=impulse(fullSys(:,inputs),stepsimtime);  %<- improve by using sss toolbox?
    impT=impT(idx);
    impX=impX(idx,:,:);
end

%reorder data
impTStep=impT(3)-impT(2);
impXper=permute(impX,[2 1 3]);
impXper=impXper(1:order(fullSys),:,:);
xMat=reshape(impXper,order(fullSys),length(impT)*size(impX,3));

% multiply with quadrature coefficient
xMat=xMat.*sqrt(impTStep);



%estimate controllability 
% wcemp=xMat*xMat';
disp(['Computation of controllability (xMat) took ' num2str(toc) 's'])
tic;
%% get POD Modes
% only do when nPODModes is specified, only then the output projection is
% required
% if ~isempty(nPODModes)
if errPODProj>0
    [Vr,S,~]=svd(inSys.c*xMat(:,1:end), 'econ'); %< S contains singular values on diagonal, same as the sqrt of the eigenvalues of xMatPod*xMatPod';
    nPODModes=1;
    err=1;
    while err>errPODProj
        err=1-sum((diag(S(1:nPODModes,1:nPODModes))).^2)/sum((diag(S)).^2);
        nPODModes=nPODModes+1;
    end
    nPODModes=min(nPODModes,size(Vr,1));
    % reduction base
    PodRed=Vr(:,1:nPODModes);      
        % PodRed=Vr*Vr';
    % end
    disp(['Computation of ' num2str(nPODModes) ' POD modes (PodRed) took ' num2str(toc) 's'])
    tic;
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

if errPODProj>0
     fullSys_adj=ss(inSys.a',inSys.c'*PodRed,ones(1,size(inSys.a,1)),0);
else
     fullSys_adj=ss(inSys.a',inSys.c',ones(1,size(inSys.a,1)),0);    
end

inputs=find(any(fullSys_adj.b));
if size(fullSys_adj.a,1)*length(stepsimtime)*size(fullSys_adj.b,2)>5000*10^9
    impAdjT=zeros(sumStepSimTimeDiv,1);
    impAdjX=zeros(sumStepSimTimeDiv,size(fullSys_adj.a,1),size(fullSys_adj.b,2));
    nPast=0;
    for iDivision=1:timeVecDivisions
        nCurr=length(stepsimtimeDiv{iDivision});
        [~, impTcurr,impXcurr]=impulse(fullSys_adj(:,inputs),stepsimtimeDiv{iDivision});
        impAdjT(nPast+1:nPast+nCurr,:)=[ impTcurr];
        impAdjX(nPast+1:nPast+nCurr,:,:)=[ impXcurr];
        nPast=nCurr+nPast;
    end
else
    [~, impAdjT,impAdjX]=impulse(fullSys_adj(:,inputs),stepsimtime);
    impAdjT=impAdjT(idx);
    impAdjX=impAdjX(idx,:,:);
end
% reorder data
impAdjTStep=impAdjT(3)-impAdjT(2);
impAdjXper=permute(impAdjX,[2 1 3]);
yMatAllt=reshape(impAdjXper,order(inSys),length(impAdjT)*size(impAdjX,3));
yMat=yMatAllt(:,1:end);
% multiply with quadrature coefficient
yMat=yMat.*sqrt(impAdjTStep);
% estimate observability
% woemp=yMat*yMat';

disp(['Computation of observability (yMat) took ' num2str(toc) 's'])
tic;
%% balancing
% build transformation matrix
% option a
%% oldway:
% [redBaseF, redBaseVal]=eig(woemp*wcemp);
% nOrderRec=rank(woemp*wcemp);
% redBaseOnlyReal=redBaseF(:,find(imag(diag(redBaseVal))==0));
% %how to find the best reduction order?
% nOrderMax=size(redBaseOnlyReal,2);
% if nOrder>nOrderMax
%     fprintf('warning, specified order of %i too large, maximum order %i\n',nOrder, nOrderMax)
% elseif isempty(nOrder)
%     nOrder=nOrderMax;
%     fprintf('warning, no order specified,  maximum order of %i, using recommended Order %i\n',nOrder, nOrderRec)
% end
% redBase=real(redBaseOnlyReal(:,1:nOrderRec));     %<- workaround for complex eigenvalues -> also not sure if that works
%% new way:
% %selcting xMat entries
% timeIdx=unique(ceil(logspace(log10(stepsimtime(2)),log10(stepsimtime(end)),100)./stepsimtime(2)));
% allidxX=repmat(timeIdx,size(xMat,2)/length(stepsimtime),1)+repmat(((0:size(xMat,2)/length(stepsimtime)-1)*length(stepsimtime))',1,length(timeIdx));
% xMatSmall=xMat(:,allidxX(:));
% timeIdx=unique(ceil(logspace(log10(stepsimtime(2)),log10(stepsimtime(end)),300)./stepsimtime(2)));
% allidxY=repmat(timeIdx,size(yMat,2)/length(stepsimtime),1)+repmat(((0:size(yMat,2)/length(stepsimtime)-1)*length(stepsimtime))',1,length(timeIdx));
% yMatSmall=yMat(:,allidxY(:));
% 
%  xMatSmall=xMat*xMat'; %wcemp
%  yMatSmall=yMat*yMat'; %woemp
% 
%  r=rank(yMatSmall*xMatSmall);
% disp(['Computation of Rank (of product yMat xMat ' num2str(r) ') took ' num2str(toc) 's'])
% tic;
% [U2, S2, V2]=svd(yMatSmall*xMatSmall,'econ');
% S2=diag(S2);
% redBase=yMatSmall*xMatSmall*V2(:,:)*diag(S2.^-.5);
% redBaseInv=diag(S2.^-.5)*U2'*yMatSmall*xMatSmall;
% [V,D,W]= eig(yMatSmall*xMatSmall);
% redBase=real(V(:,1:r));
% redBaseInv=real(W(1:r,:));
%% real balanced truncation
% QUICK FIX get cholesky factors of wcemp and woemp by setting eigenvalues of wcemp
% and woemp to zero
% [Vo,Do]=eig(yMat*yMat');
% Do=diag(Do);
% Do(Do<0)=0;
% [Qo,Ro]=qr((Vo*diag(sqrt(Do)))');
% 
% yMatMod=(Qo*Ro)'; %ZQ
% 
% [Vc,Dc]=eig(xMat*xMat');
% Dc=diag(Dc);
% Dc(Dc<0)=0;
% [Qc,Rc]=qr((Vc*diag(sqrt(Dc)))');
% 
% xMatMod=(Qc*Rc)'; %ZP
% 
% % svd
% [U, S, V]=svd(yMatMod'*xMatMod,'econ');
% r=min(find(diag(S)<10e-10));
% r=size(S,1);
% %
% 
% redBase=xMatMod*V(:,1:r)*diag(diag(S(1:r,1:r)).^-.5);
% redBaseInv=(yMatMod*U(:,1:r)*diag(diag(S(1:r,1:r)).^-.5))';
%% real balanced truncation
% better way of getting "cholesky type" factors from wcemp and woemp:

% [Uo,So,~]=svd(yMat,'econ');
% [Qo,Ro]=qr((Uo*So)');
% yMatMod=(Qo*Ro)'; %ZQ
% 
% [Uc,Sc,~]=svd(xMat,'econ');
% [Qc,Rc]=qr((Uc*Sc)');
% xMatMod=(Qc*Rc)'; %ZP
% 
% [U, S, V]=svd(yMatMod'*xMatMod,'econ');
% r=size(S,1);
% redBase=xMatMod*V(:,1:r)*diag(diag(S(1:r,1:r)).^-.5);
% redBaseInv=(yMatMod*U(:,1:r)*diag(diag(S(1:r,1:r)).^-.5))';
% disp(['Computation of balancing transformation (redBase) took ' num2str(toc) 's'])
%% square root balance truncation
% better way of getting "cholesky type" factors from wcemp and woemp:

[Uo,So,~]=svd(yMat*yMat','econ');
yMatMod=Uo*So.^0.5; %ZQ

[Uc,Sc,~]=svd(xMat*xMat','econ');
xMatMod=Uc*Sc.^0.5; %ZP

[U, S, V]=svd(yMatMod'*xMatMod,'econ');

r=size(S,1);
 redBase=xMatMod*V(:,1:r)*diag(diag(S(1:r,1:r)).^-.5);
 redBaseInv=(yMatMod*U(:,1:r)*diag(diag(S(1:r,1:r)).^-.5))';

S=diag(S);
disp(['Computation of balancing transformation (redBase) took ' num2str(toc) 's'])

%% 
% option b
% [a,b,c]=svd(yMat'*xMat,'econ');
% T=(xMat*c*b^-0.5);
% S=(b^-0.5*a'*yMat');
% redBase=T;



end

