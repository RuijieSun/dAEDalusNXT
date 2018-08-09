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

  out=[0 -in(3) in(2); in(3) 0 -in(1); -in(2) in(1) 0;];

end

