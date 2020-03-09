clear; close all; clc;
%% Load data
load('up1_1Full.mat')

%%
%%% Pull data out of structure
u = up1_1FULL.deg_30.u(:,:,end);
v = up1_1FULL.deg_30.v(:,:,end);
nx = up1_1FULL.deg_30.nx;
ny = up1_1FULL.deg_30.ny;
mask_ind = up1_1FULL.deg_30.mask_inds;

%%% Stack the data, remove mask, and replace outliers
[X] = stackPCA(u,v,nx, ny, mask_ind);

%%% Put back mask NaNs, unstack the data
%[u_new] = unstackPCA(X,u,[],mask_ind);


