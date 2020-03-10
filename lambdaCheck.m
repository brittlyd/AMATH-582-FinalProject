function [] = lambdaCheck(uL,vL, uN, vN)

% plots for checking lambda
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
