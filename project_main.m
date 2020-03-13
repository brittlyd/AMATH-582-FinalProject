

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
    lambda = 5; % choose your sparsity constant value (1 is often a good starting point)
    tol = 1e-7; % set your tolerance
    maxIter = 1000; % set your maximum number of iterations
    [up1_1FULL.(rotFields{iRot}).L, up1_1FULL.(rotFields{iRot}).N,mask_log] ... 
        = rPCA_main(u,v,nx,ny, mask_ind, lambda, tol, maxIter);
    % average all columns of L,N to get "rPCA cleaned" phase average
    up1_1FULL.(rotFields{iRot}).Lr = mean(up1_1FULL.(rotFields{iRot}).L,2);
    up1_1FULL.(rotFields{iRot}).Nr = mean(up1_1FULL.(rotFields{iRot}).N,2);
end
% checks for rPCA
% put mask back in
[uL,vL] = unstackPCA(up1_1FULL.deg_30.Lr,nx,ny,mask_log(:,:,1),1);
[uN,vN] = unstackPCA(up1_1FULL.deg_30.Nr,nx,ny,mask_log(:,:,1),1);
lambdaCheck(uL, vL,uN, vN, lambda)

figure
subplot(2,2,1)
pcolor(up1_1FULL.deg_30.u_avg)
title('u avg')
shading flat
caxis([-1.5 2])
colorbar

subplot(2,2,3)
pcolor(up1_1FULL.deg_30.v_avg)
title('v avg')
shading flat
caxis([-1.5 2])
colorbar
%% SVD without rPCA

run='Abby'; %what to append to all plot saving so things don't get overwritten between data sets
%load data
load("C:\Users\abber\Documents\School\Grad School\Winter 20\AMATH 582\Project\up1_1 Crop.mat")
rotFields = fieldnames(data);
nx = 52;
ny = 52;
uv = 1; %if 1 run for u and v, if 0 run for vmag
if uv
Y = zeros(nx*ny*2, nRot);
else
Y = zeros(nx*ny, nRot);
end
for iRot = 1:nRot
    %Y(:,iRot) = reshape(data.(rotFields{iRot}).Interp.Vort_crop, [nx*ny 1]);
    %Y(:,iRot) = reshape(data.(rotFields{iRot}).InterpCommon.Vmag_crop, [nx*ny 1]);
    if uv
        u = reshape(data.(rotFields{iRot}).interp.u_crop, [nx*ny 1]);
        v = reshape(data.(rotFields{iRot}).interp.v_crop, [nx*ny 1]);
        Y(:,iRot) = [u;v];
    else
        Y(:,iRot) = reshape(data.(rotFields{iRot}).interp.vmag_crop, [nx*ny 1]);
    end
end

figure
set(gcf,'position',1.0e+03 *[0.0016    0.2079    1.4600    0.5120])
xcrop=data.(rotFields{1}).interp.xcrop;
ycrop=data.(rotFields{1}).interp.ycrop;
[ha, pos]= tight_subplot(2,4,[0 0],[.01 .01],[.01 .01]);
p=1;
for n=1:ceil(nRot/8):nRot
    axes(ha(p))
    ax=gca;
    tmpU=data.(rotFields{n}).interp.u_crop;
    tmpV=data.(rotFields{n}).interp.v_crop;
    tmpVmag=data.(rotFields{n}).interp.vmag_crop;
    pcolor(data.(rotFields{n}).interp.xcrop,data.(rotFields{n}).interp.ycrop,tmpVmag)
    hold on
    quiver(data.(rotFields{n}).interp.xcrop,data.(rotFields{n}).interp.ycrop,tmpU,tmpV,2,'k')
    plot(data.(rotFields{n}).interp.foil,'facecolor',[0 0 0],'facealpha',0.5...
        ,'edgecolor','none')
    ax.XAxis.Visible='off';
    ax.YAxis.Visible='off';
    title(strcat('\theta= ',(rotFields{n}(strfind(rotFields{n},'_')+1:end))))
    axis equal
    axis tight
    shading interp
    caxis([0,3])
    set(gca,'position',pos{p})
    p=p+1;
end
print(gcf,strcat('velFields',run),'-dpng','-r600')
%% mean-subtract and fill in NaNs for SVD
Yavg = mean(Y,2,'omitnan'); %compute row mean for subtraction
Yms =Y-Yavg*ones(1,size(Y,2)); % Y mean-subtracted
Yms(isnan(Yms))=0; %replace NaN

%Remove values where the number of NaN in that row is greater than the
%threshold
Ynan = ismissing(Y);
nanRow = sum(Ynan,2); % number of NaNs in each row
th = 4; % threshold for number of NaNs allowed before row is masked
Yms(nanRow>th,:)=[];

% subplot(1,2,2)
% imagesc(Yfill)
% title('Mean-subtracted, NaNs filled in')

[U,S,V] = svd(Yms,'econ');
sig=diag(S);
energy=sig/sum(sig)*100;
for n=1:length(sig)
    energytotal(n)=sum(energy(1:n));
end
figure
pcolor(V)
if uv
    print(gcf,strcat('Vuv matrix',run),'-dpng','-r600')
else
    print(gcf,strcat('Vvmag matrix',run),'-dpng','-r600')
end

figure
set(gcf,'position',[227.3000  403.9000  720.5571  316.1000])
left_color = [0.9153    0.2816    0.2878];
right_color = [0 .5 .5];
set(gcf,'defaultAxesColorOrder',[left_color; right_color]);
yyaxis left
a=plot(energy,'o','markerfacecolor',left_color);
hold on
ind=find(energytotal>=90,1);
plot(ind,energy(ind),'o','markerfacecolor',[0 0 0],'markersize',12)
ylabel('% of energy')
axis tight
yyaxis right
plot(energytotal,'o','markerfacecolor',right_color)
hold on
plot(linspace(1,length(sig),100),90*ones(1,100),'--k','linewidth',1.75)
ylabel('% of Energy Captured')
xlabel('mode')
% set(gca, 'SortMethod', 'depth')
axis tight
grid on
print(gcf,strcat('Energy Cropped',run),'-dpng','-r600')

%% show modes
%put NaNs back in
Ureplaced=zeros(size(Y));
Ureplaced(nanRow>th,:)=NaN;
Ureplaced(~(nanRow>th),:)=U;

f=figure;
set(f,'position',[331.8571  238.1429  822.8571  481.8571])
[ha, pos]= tight_subplot(2,2,[.05 .15],[.1 .01],[.1 .1]);
f2=figure;
set(gcf,'position',[331.8571  238.1429  822.8571  481.8571])
[ha2, pos2]= tight_subplot(2,2,[.05 .15],[.1 .01],[.1 .1]);
for k = 1:4
    f
    axes(ha(k))
    %Put NaNs back in
    if uv
        mode=Ureplaced(1:end/2,k);
        mode=mode/max(mode,[],'all');
        mode2=Ureplaced(end/2+1:end,k);
        mode2=mode2/max(mode2,[],'all');
        titletext='u modes';
        titletext2='v modes';
    else
        mode=Ureplaced(:,k);
        mode=mode/max(mode,[],'all');
        titletext='vmag modes';
    end
    pcolor(xcrop,ycrop,reshape(mode, [nx ny]))
    hold on
    plot(data.deg_48.interp.foil,'facecolor',[0 0 0],'facealpha',0.5...
        ,'edgecolor','none')
    title(strcat('mode ', num2str(k)))
    xlabel('x/c')
    ylabel('y/c')
    shading interp
    axis equal
    axis tight
    colorbar
    set(gca,'position',pos{k})
    
    if uv
        f2
        axes(ha2(k))
        pcolor(xcrop,ycrop,reshape(mode2, [nx ny]))
        hold on
        plot(data.deg_48.interp.foil,'facecolor',[0 0 0],'facealpha',0.5...
            ,'edgecolor','none')
        title(strcat('mode ', num2str(k)))
        xlabel('x/c')
        ylabel('y/c')
        shading interp
        axis equal
        axis tight
        colorbar
        set(gca,'position',pos2{k})
    end
end
print(f,strcat(titletext,run),'-dpng','-r600')
if uv
   print(f2,strcat(titletext2,run),'-dpng','-r600')
end
%% Reconstruct Each Angle 
if ~uv
figure
set(gcf,'position',1.0e+03 *[0.0016    0.2079    1.4600    0.5120])
[ha, pos]= tight_subplot(2,4,[0 0],[.01 .01],[.01 .01]);
p=1;
indRecon=ceil(size(Ureplaced,2)/2); %Reconstructs with 50 perc. of the modes
% indRecon=ind; %Recontructs with 90 perc. of the energy
for n=1:ceil(nRot/8):nRot
    axes(ha(p))
    ax=gca;
    %Reconstruct and add mean back
    tmpVmag=Ureplaced(:,1:indRecon)*S(1:indRecon,1:indRecon)*V(n,1:indRecon)'+Yavg;
    pcolor(data.(rotFields{n}).interp.xcrop,data.(rotFields{n}).interp.ycrop...
        ,reshape(tmpVmag, [nx ny]))
    hold on
    plot(data.(rotFields{n}).interp.foil,'facecolor',[0 0 0],'facealpha',0.5...
        ,'edgecolor','none')
    ax.XAxis.Visible='off';
    ax.YAxis.Visible='off';
    title(strcat('\theta= ',rotFields{n}(strfind(rotFields{n},'_')+1:end)))
    axis equal
    axis tight
    shading interp
    caxis([0 3])
    set(gca,'position',pos{p})
    p=p+1;
end
print(gcf,strcat('velFieldsReconstruct',run),'-dpng','-r600')
end
