%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%
function [c , ceq]= GLA_Constraintsfun(h,N,Ts,GLAS_ActuatorModel,gust_inputs)

gust_input1 = gust_inputs(:,end); % !!!Just looking at rate limitation for longest gust might not be sufficiently constraining!!! Simon: no its not.
gust_input1 =reshape(gust_inputs, (size(gust_inputs,1)*size(gust_inputs,2)) ,1);   %looking at all of them 
c=[];
for iCs=1:length(GLAS_ActuatorModel.UB)
    command_gust = conv(h(  (iCs-1)*N+1:iCs*N), gust_input1);
    command_gust_dot=diff(command_gust)/Ts;
    % output4_gust = conv(h(3*N+1:4*N), gust_input1);
    % output5_gust = conv(h(4*N+1:5*N), gust_input1);


    c = [ c;  command_gust-GLAS_ActuatorModel.UB(iCs);          
            -command_gust+GLAS_ActuatorModel.LB(iCs);             
            abs(command_gust_dot)-GLAS_ActuatorModel.maxRate;   ];
end
    ceq=[];
end

