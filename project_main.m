
%% rPCA
clear all, close all
% load data, initialize
filename = 'up1_1FULL.mat';
load(filename)
rotFields = fieldnames(up1_1FULL);
nRot = length(rotFields);

% keep these variables outside loop if same for all rotations
nx = up1_1FULL.deg_30.nx; 
ny = up1_1FULL.deg_30.ny;
mask_ind = up1_1FULL.deg_30.mask_inds;

for iRot = 1:nRot
    u = up1_1FULL.(rotFields{iRot}).u;
    v = up1_1FULL.(rotFields{iRot}).v;
    lambda = 1; % choose your sparsity constant value (1 is often a good starting point)
    tol = 1e-7; % set your tolerance
    maxIter = 1000; % set your maximum number of iterations
    [up1_1FULL.(rotFields{iRot}).uL, up1_1FULL.(rotFields{iRot}).vL, up1_1FULL.(rotFields{iRot}).uN, up1_1FULL.(rotFields{iRot}).vN] ... 
        = rPCA_main(u,v,nx,ny, mask_ind, lambda, tol, maxIter);
end
%% SVD without rPCA

%load data
load('up1_1 Crop.mat')
rotFields = fieldnames(data);
nRot = length(rotFields);
nx = 52;
ny = 52;
Y = zeros(nx*ny, nRot);
for iRot = 1:nRot
    %Y(:,iRot) = reshape(data.(rotFields{iRot}).Interp.Vort_crop, [nx*ny 1]);
    %Y(:,iRot) = reshape(data.(rotFields{iRot}).InterpCommon.Vmag_crop, [nx*ny 1]);
    Y(:,iRot) = reshape(data.(rotFields{iRot}).Interp.Vmag_crop, [nx*ny 1]);
end

figure
subplot(1,2,1)
imagesc(Y)
title('Original data')

%% mean-subtract and fill in NaNs for SVD

Yavg = mean(Y,2,'omitnan'); % compute row mean for subtraction
Yms = Y-Yavg*ones(1,size(Y,2)); % Y mean-subtracted
Yfill = fillmissing(Yms, 'constant', 0); %replace NaN
subplot(1,2,2)
imagesc(Yfill)
title('Mean-subtracted, NaNs filled in')

[U,S,V] = svd(Yfill,'econ');
%% replace NaNs

Ynan = ismissing(Y);
nanRow = sum(Ynan,2); % number of NaNs in each row
th = 4; % threshold for number of NaNs allowed before row is masked
for iRow = 1:size(Y,1)
    if nanRow(iRow) > th
        U(iRow,:) = NaN;
    end
end
%% show modes

figure
set(gcf,'position',[331.8571  238.1429  822.8571  481.8571])
Xcrop=data.deg_48.InterpCommon.Xcrop;
Ycrop=data.deg_48.InterpCommon.Ycrop;
[ha, pos]= tight_subplot(2,2,[.05 .15],[.1 .01],[.1 .1]);
for k = 1:4
    axes(ha(k))
    pcolor(Xcrop,Ycrop,reshape(U(:,k), [nx ny]))
    hold on
    plot(data.deg_48.Interp.foil,'facecolor',[0 0 0],'facealpha',0.5...
        ,'edgecolor','none')
    title(strcat('mode ', num2str(k)))
    xlabel('x/c')
    ylabel('y/c')
    shading flat
    axis equal
    axis tight
    colorbar
    set(gca,'position',pos{k})
end

figure
subplot(1,5,1)
imagesc(reshape(Yavg,[nx ny]))
title('mean')
for k = 1:4
    subplot(1,5,k+1)
    imagesc(reshape(U(:,k), [nx ny]))
    title({'mode ',k})
end

figure
imagesc(V(:,1:4))
title('V')
