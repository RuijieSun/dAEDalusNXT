%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%
function [F] = GLA_Costfun(h, N, MassCoeffBend, StaticBending, TimeResponses_ContTurb)


%% simplified F computation

%% simplified F computation (without the section loop)
turbBending=[TimeResponses_ContTurb(:,4:6:end,end);zeros(N-1,length(StaticBending))]  ;

nCs=length(h)/N;
for iCs=1:nCs
    commandBendingAllSections(iCs,:,:)= conv2(h((iCs-1)*N+1:iCs*N),1,TimeResponses_ContTurb(:,4:6:end,iCs));
end
dynBendingAllSections=squeeze(sum(commandBendingAllSections,1))+turbBending;
bending=(abs(StaticBending)'+rms(dynBendingAllSections));
F=sum(MassCoeffBend.*bending);
end




