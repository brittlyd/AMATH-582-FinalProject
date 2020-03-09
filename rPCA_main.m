clear all, close all

%% load data
filename = 'up1_1FULL.mat';
load(filename)
rotFields = fieldnames(up1_1FULL);
nRot = length(rotFields);
nx = up1_1FULL.deg_30.nx;
ny = up1_1FULL.deg_30.ny;
u = up1_1FULL.deg_30.u;
v = up1_1FULL.deg_30.v;
mask_ind = up1_1FULL.deg_30.mask_inds;
%% plots to check data input

figure
subplot(2,2,1)
pcolor(u(:,:,1))
shading flat
title ('first frame of u')

subplot(2,2,3)
pcolor(v(:,:,1))
shading flat
title ('first frame of v')

%check stack & restack commands
[Y,mask_log] = stackPCA(u(:,:,1),v(:,:,1),nx,ny,mask_ind);
[u_new,v_new] = unstackPCA(Y,nx,ny,mask_log, 1);

subplot(2,2,2)
pcolor(u_new)
shading flat
title ('first frame of restacked u')

subplot(2,2,4)
pcolor(v_new)
shading flat
title ('first frame of restacked v')
%% Remove mask, replace other NaNs with outliers

nt = size(u,3); % number of frames
X = zeros((nx*ny-length(mask_ind))*2, nt);
for i = 1:nt
    [X(:,i),mask_log(:,:,i)] = stackPCA(u(:,:,i),v(:,:,i),nx,ny,mask_ind);
end



%% rPCA
addpath RPCA-PIV-master;

lambda = 1; % choose your sparsity constant value (1 is often a good starting point)
tol = 1e-7; % set your tolerance
maxIter = 1000; % set your maximum number of iterations
[L, N, ~] = inexact_alm_rpca(X, lambda, tol, maxIter); %L is lowrank, N is noise

%% optional plots

% figure
% subplot(1,3,1)
% imagesc(X)
% title('Original data')
% subplot(1,3,2)
% pcolor(L)
% shading flat
% title('L')
% subplot(1,3,3)
% pcolor(N)
% shading flat
% title('N')


%% put mask back in
isv =  ~isempty(v); %check if there are two arrays
    
[uL,vL] = unstackPCA(L(:,1),nx,ny,mask_log(:,:,1),isv);
[uN,vN] = unstackPCA(N(:,1),nx,ny,mask_log(:,:,1),isv);
%%
%optional plots for checking lambda
figure
subplot(2,2,1)
pcolor(uL)
shading flat
title ('first column of uL')

subplot(2,2,2)
pcolor(uN)
shading flat
title ('first column of uN')

subplot(2,2,3)
pcolor(vL)
shading flat
title ('first column of vL')

subplot(2,2,4)
pcolor(vN)
shading flat
title ('first column of vN')
