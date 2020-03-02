clear all, close all

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