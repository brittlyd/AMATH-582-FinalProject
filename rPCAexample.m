function [L, S] = rPCAexample(X, u, v)
% Inputs
% X: data matrix, columns are time and rows are space
% u,v: velocity

%requires supporting functions from github.com/ischerl/RPCA-PIV

%% Run RPCA
lambda = 1; % choose your sparsity constant value (1 is often a good starting point)
tol = 1e-7; % set your tolerance
maxIter = 1000; % set your maximum number of iterations
[L, S, ~] = inexact_alm_rpca(X, lambda, tol, maxIter); %L is lowrank, S is sparse

%% Reshape your u and v low rank and sparse data
% 
% Lu = reshape(L(1:end/2, : ), m(1), m(2), m(3)); Lv = reshape(L(1+end/2:end, : ), m(1), m(2), m(3));
% Su = reshape(S(1:end/2, : ), m(1), m(2), m(3)); Sv = reshape(S(1+end/2:end, : ), m(1), m(2), m(3));

%% Calculate turbulence Spectra
% w = zeroes(size(u));
% %for original data 
% [k, E] = turbspec(u, v, w); % if you do not have a third velocity dimension, use zeroes(size(u))
% %for low rank
% [k_L, E_L] = turbspec(Lu, Lv, Lw); % if you do not have a third velocity dimension, use zeroes(size(Lu))
