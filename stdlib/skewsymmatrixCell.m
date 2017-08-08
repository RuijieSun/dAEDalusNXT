%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%
function [ out ] = skewsymmatrixCell( in )
%SKEWSYMMATRIX Summary of this function goes here
%   Detailed explanation goes here
for i=1:length(in)
  out{i}=[0 -in{i}(3) in{i}(2); in{i}(3) 0 -in{i}(1); -in{i}(2) in{i}(1) 0;];
end
end

