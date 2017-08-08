%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%
function [RFA_MS_appx]=compute_appx_from_GAF(Q,k,aircraft,n_modes,name)
%% assuming control modes are at the end, gust modes are not implemented yet! careful!
Q=Q([1:n_modes end-length(aircraft.control_surfaces)+1:end],[1:n_modes end-length(aircraft.control_surfaces)+1:end],:);
k=k(1:end);
gamma=k;
[A0, A1, A2, Arest] = rogers_state_space_approximation(Q,k,gamma);
[A0ms,A1ms,A2ms,Ems,Dms,gamma_ms] = minimum_state_approximation_opt(Q,k,gamma,60);
RFA_MS_appx.A0=A0;
RFA_MS_appx.A1=A1;
RFA_MS_appx.A2=A2;
RFA_MS_appx.Arest=Arest;
RFA_MS_appx.A0ms=A0ms;
RFA_MS_appx.A1ms=A1ms;
RFA_MS_appx.A2ms=A2ms;
RFA_MS_appx.Ems=Ems;
RFA_MS_appx.Dms=Dms;
RFA_MS_appx.gamma_ms=gamma_ms;
save([aircraft.name '/' 'appx_' name],'RFA_MS_appx');

