clear all, close all

%load data
load("C:\Users\abber\Documents\School\Grad School\Winter 20\AMATH 582\Project\up1_1 Crop.mat")
rotFields = fieldnames(data);
nRot = length(rotFields);
nx = 52;
ny = 52;
Y = zeros(nx*ny*2, nRot);
for iRot = 1:nRot
    %Y(:,iRot) = reshape(data.(rotFields{iRot}).Interp.Vort_crop, [nx*ny 1]);
    %Y(:,iRot) = reshape(data.(rotFields{iRot}).InterpCommon.Vmag_crop, [nx*ny 1]);
%     Y(:,iRot) = reshape(data.(rotFields{iRot}).Interp.Vmag_crop, [nx*ny 1]);
    U = reshape(data.(rotFields{iRot}).Interp.U_crop, [nx*ny 1]);
    V = reshape(data.(rotFields{iRot}).Interp.V_crop, [nx*ny 1]);
    Y(:,iRot) = [U;V];
end

figure
set(gcf,'position',1.0e+03 *[0.0016    0.2079    1.4600    0.5120])
Xcrop=data.(rotFields{1}).Interp.Xcrop;
Ycrop=data.(rotFields{1}).Interp.Ycrop;
[ha, pos]= tight_subplot(2,4,[0 0],[.01 .01],[.01 .01]);
p=1;
for n=1:ceil(nRot/8):nRot
    axes(ha(p))
    ax=gca;
    tmpU=data.(rotFields{n}).Interp.U_crop;
%     tmpU(mask)=NaN;
    tmpV=data.(rotFields{n}).Interp.V_crop;
%     tmpV(mask)=NaN;
    tmpVmag=data.(rotFields{n}).Interp.Vmag_crop;
%     tmpVmag(mask)=NaN;
    tmpVort=data.(rotFields{n}).Interp.Vort_crop;
%     tmpVmag(mask)=NaN;
    pcolor(data.(rotFields{n}).Interp.Xcrop,data.(rotFields{n}).Interp.Ycrop,tmpVort)
    hold on
    quiver(data.(rotFields{n}).Interp.Xcrop,data.(rotFields{n}).Interp.Ycrop,tmpU,tmpV,2,'k')
    plot(data.(rotFields{n}).Interp.foil,'facecolor',[0 0 0],'facealpha',0.5...
        ,'edgecolor','none')
    ax.XAxis.Visible='off';
    ax.YAxis.Visible='off';
    title(strcat('\theta= ',(rotFields{n}(strfind(rotFields{n},'_')+1:end))))
    axis equal
    axis tight
    shading interp
%     caxis([0,3])
    set(gca,'position',pos{p})
    p=p+1;
end
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
text(ind,energy(ind),strcat('90% of energy is captured, mode '...
    ,num2str(ind),'\rightarrow  '),'VerticalAlignment','bottom','Fontsize'...
    ,10,'HorizontalAlignment','right');
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
print(gcf,'Energy Cropped Abby','-dpng','-r600')

%% show modes
%put NaNs back in
Ureplaced=zeros(size(Y));
Ureplaced(nanRow>th,:)=NaN;
Ureplaced(~(nanRow>th),:)=U;

figure
set(gcf,'position',[331.8571  238.1429  822.8571  481.8571])
Xcrop=data.deg_48.InterpCommon.Xcrop;
Ycrop=data.deg_48.InterpCommon.Ycrop;
[ha, pos]= tight_subplot(2,2,[.05 .15],[.1 .01],[.1 .1]);
for k = 1:4
    axes(ha(k))
    %Put NaNs back in
    mode=Ureplaced(:,k);
%     mode=zeros(nx*ny,1);
%     mode(~mask)=Ureplaced(:,k);
%     mode(mask)=NaN;
    pcolor(Xcrop,Ycrop,reshape(Ureplaced(:,k), [nx ny]))
    hold on
    plot(data.deg_48.Interp.foil,'facecolor',[0 0 0],'facealpha',0.5...
        ,'edgecolor','none')
    title(strcat('mode ', num2str(k)))
    xlabel('x/c')
    ylabel('y/c')
    shading interp
    axis equal
    axis tight
    colorbar
    set(gca,'position',pos{k})
end

% figure
% subplot(1,5,1)
% imagesc(reshape(Yavg,[nx ny]))
% title('mean')
% for k = 1:4
%     subplot(1,5,k+1)
%     imagesc(reshape(U(:,k), [nx ny]))
%     title({'mode ',k})
% end
% 
% figure
% imagesc(V(:,1:4))
% title('V')
