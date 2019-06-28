%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%
function [c , ceq]= GLA_Constraintsfun(h,N,Ts,GLAS_ActuatorModel,gust_inputs)

c=[];
for iLength=1:size(gust_inputs,2)
    gust_input1 = gust_inputs(:,iLength)';
    for iCs=1:length(GLAS_ActuatorModel.UB)
        command_gust = conv(h(  (iCs-1)*N+1:iCs*N), gust_input1);
        command_gust_dot=diff(command_gust)/Ts;


        c = [ c  command_gust-GLAS_ActuatorModel.UB(iCs)   ...       
                -command_gust+GLAS_ActuatorModel.LB(iCs)...             
                abs(command_gust_dot)-GLAS_ActuatorModel.maxRate   ];
    end
end
    ceq=[];
end

