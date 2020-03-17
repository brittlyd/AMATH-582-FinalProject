function [] = lambdaCheck(uL,vL, uN, vN, lambda)

% plots for checking lambda
figure
subplot(2,2,1)
pcolor(uL)
shading flat
title ('first column of uL')
colorbar
%uLim = caxis;
caxis([-1.5 2])

subplot(2,2,2)
pcolor(uN)
shading flat
title ('first column of uN')
colorbar
%caxis(uLim)
caxis([-1.5 2])

subplot(2,2,3)
pcolor(vL)
shading flat
title ('first column of vL')
colorbar
%vLim = caxis;
caxis([-1.5 2])

subplot(2,2,4)
pcolor(vN)
shading flat
title ('first column of vN')
colorbar
%caxis(vLim)
caxis([-1.5 2])

sgtitle({'$\lambda = $' lambda }, 'Interpreter', 'latex')