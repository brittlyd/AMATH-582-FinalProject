function [] = checkStack(u, v, nx ,ny, mask_ind)
%check stack & restack commands

figure
subplot(2,2,1)
pcolor(u(:,:,1))
shading flat
title ('first frame of u')

subplot(2,2,3)
pcolor(v(:,:,1))
shading flat
title ('first frame of v')

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