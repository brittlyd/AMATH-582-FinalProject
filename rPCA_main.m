function [L, N, mask_log] = rPCA_main(u,v,nx,ny, mask_ind, lambda, tol, maxIter)
% inputs
% u,v are velocity fields with dimensions nx, ny
% mask_ind is locations that should remain NaNs after processing
% v can be empty 



isv =  ~isempty(v); %check if there are two arrays
%% Remove mask, replace other NaNs with outliers

nt = size(u,3); % number of frames
if isv
    X = zeros((nx*ny-length(mask_ind))*2, nt); % data matrix
    for i = 1:nt
        [X(:,i),mask_log(:,:,i)] = stackPCA2(u(:,:,i),v(:,:,i),nx,ny,mask_ind);
    end
else
    X = zeros((nx*ny-length(mask_ind)), nt);
    for i = 1:nt
        [X(:,i),mask_log(:,:,i)] = stackPCA2(u(:,:,i),[],nx,ny,mask_ind);
    end
end

% avgX = mean(X,2)*ones(1,size(X,2));
% X = X - avgX;
%% rPCA
addpath(genpath('RPCA-PIV-master'));

[L, N, ~] = inexact_alm_rpca(X, lambda, tol, maxIter); %L is lowrank, N is noise


