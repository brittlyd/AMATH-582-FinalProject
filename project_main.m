% clear all, close all

%load data
load('up1_1 Crop.mat')
rotFields = fieldnames(data);
nRot = length(rotFields);
nx = 52;
ny = 52;
Y = zeros(nx*ny, nRot);
for iRot = 1:nRot
    Y(:,iRot) = reshape(data.(rotFields{iRot}).Interp.Vort_crop, [nx*ny 1]);
end

Yavg = mean(Y,2,'omitnan'); % compute row mean for subtraction
figure
subplot(1,2,1)
imagesc(reshape(Yavg,[nx ny]))
title('average vorticity')
%colMean = mean(Y,1,'omitnan'); %compute column mean for NaN substitution
Yavg = fillmissing(Yavg, 'constant', 0); %replace NaN
subplot(1,2,2)
imagesc(reshape(Yavg,[nx ny]))
title('average vorticity, NaNs filled')

Yfill = fillmissing(Y, 'constant', 0); %replace NaN

[U,S,V] = svd(Yfill-Yavg*ones(1,size(Y,2)),'econ');

%% show modes

figure
set(gcf,'position',[331.8571  238.1429  822.8571  481.8571])
Xcrop=data.deg_48.InterpCommon.Xcrop;
Ycrop=data.deg_48.InterpCommon.Ycrop;
[ha pos]= tight_subplot(2,2,[.05 .15],[.1 .01],[.1 .1]);
for k = 1:4
    axes(ha(k))
    pcolor(Xcrop,Ycrop,reshape(U(:,k), [nx ny]))
    hold on
    plot(data.deg_48.InterpCommon.foil,'facecolor',[0 0 0],'facealpha',0.5...
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