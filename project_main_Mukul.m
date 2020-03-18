clear all, close all

run='Mukul';%what to append to all plot saving so things don't get 
%overwritten between data sets

%load data
load("C:\Users\abber\Documents\School\Grad School\Winter 20\AMATH 582\Project\Mukul flow_data Crop.mat")
nRot = length(data);
nRot = 32; %for just the upstream portion 
nx = 52;
ny = 52;
uv = 1; %if 1 run for u and v, if 0 run for vmag
if uv
Y = zeros(nx*ny*2, nRot);
else
Y = zeros(nx*ny, nRot);
end
for iRot = 1:nRot
    %Y(:,iRot) = reshape(data(iRot).Interp.Vort_crop, [nx*ny 1]);
    %Y(:,iRot) = reshape(data(iRot).InterpCommon.Vmag_crop, [nx*ny 1]);
    if uv
        u = reshape(data(iRot).Interp.U_crop, [nx*ny 1]);
        v = reshape(data(iRot).Interp.V_crop, [nx*ny 1]);
        Y(:,iRot) = [u;v];
    else
        Y(:,iRot) = reshape(data(iRot).Interp.Vmag_crop, [nx*ny 1]);
    end
end

figure
set(gcf,'position',1.0e+03 *[0.0016    0.0181    1.4600    0.7018])
Xcrop=data(1).Interp.Xcrop;
Ycrop=data(1).Interp.Ycrop;
[ha, pos]= tight_subplot(2,4,[0 0],[.01 .01],[.01 .01]);
p=1;
for n=2:4:nRot
    axes(ha(p))
    ax=gca;
    set(gca,'Visible','off')
    tmpU=data(n).Interp.U_crop;
    tmpV=data(n).Interp.V_crop;
    tmpVmag=data(n).Interp.Vmag_crop;
    tmpVort=data(n).Interp.Vort_crop;
    pcolor(data(n).Interp.Xcrop,data(n).Interp.Ycrop,tmpVmag)
    hold on
    quiver(data(n).Interp.Xcrop(1:2:end),data(n).Interp.Ycrop(1:2:end)...
        ,tmpU(1:2:end,1:2:end),tmpV(1:2:end,1:2:end),2,'k')
    plot(data(n).Interp.foil,'facecolor',[0 0 0],'facealpha',0.5...
        ,'edgecolor','none')
    ax.XAxis.Visible='off';
    ax.YAxis.Visible='off';
    title(strcat('\theta= ',num2str(data(n).theta)))
    axis equal
    axis tight
    shading interp
    caxis([0,3])
    set(gca,'position',pos{p})
    p=p+1;
end
c=colorbar('southoutside');
c.FontSize=12;
set(c,'position',[0.3870    0.4610    0.2201    0.0347])
set(get(c,'title'),'string','$(\frac{V_{mag}}{U_\infty})$','interpreter','latex');
print(gcf,strcat('velFields',run),'-dpng','-r600')

%% mean-subtract and remove NaNs for SVD

Yavg = mean(Y,2,'omitnan'); %compute row mean for subtraction
Yms =Y-Yavg*ones(1,size(Y,2)); % Y mean-subtracted
Yms(isnan(Yms))=0; %replace NaN

%Remove values where the number of NaN in that row is greater than the
%threshold
Ynan = ismissing(Y);
nanRow = sum(Ynan,2); % number of NaNs in each row
th = nRot/2; % threshold for number of NaNs allowed before row is masked
Yms(nanRow>th,:)=[];

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
set(f,'position',[331.8571   13.5714  822.8571  706.4286])
[ha, pos]= tight_subplot(3,2,[.05 .15],[.1 .01],[.1 .1]);
f2=figure;
set(gcf,'position',[331.8571   13.5714  822.8571  706.4286])
[ha2, pos2]= tight_subplot(3,2,[.05 .15],[.1 .01],[.1 .1]);
p=1;
pp=1;
for k = [1:5,ind]
    axes(ha(p))
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
    pcolor(Xcrop,Ycrop,reshape(mode, [nx ny]))
    hold on
    plot(data(1).Interp.foil,'facecolor',[0 0 0],'facealpha',0.5...
        ,'edgecolor','none')
    title(strcat('mode ', num2str(k)))
    xlabel('x/c')
    ylabel('y/c')
    shading interp
    axis equal
    axis tight
    colorbar
    set(gca,'position',pos{p})
    p=p+1;
    if uv
        f2
        axes(ha2(pp))
        pcolor(Xcrop,Ycrop,reshape(mode2, [nx ny]))
        hold on
        plot(data(1).Interp.foil,'facecolor',[0 0 0],'facealpha',0.5...
            ,'edgecolor','none')
        title(strcat('mode ', num2str(k)))
        xlabel('x/c')
        ylabel('y/c')
        shading interp
        axis equal
        axis tight
        colorbar
        set(gca,'position',pos2{pp})
        pp=pp+1;
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
for n=2:4:nRot
    axes(ha(p))
    ax=gca;
    %Reconstruct and add mean back
    tmpVmag=Ureplaced(:,1:indRecon)*S(1:indRecon,1:indRecon)*V(n,1:indRecon)'+Yavg;
    pcolor(data(n).Interp.Xcrop,data(n).Interp.Ycrop,reshape(tmpVmag, [nx ny]))
    hold on
    plot(data(n).Interp.foil,'facecolor',[0 0 0],'facealpha',0.5...
        ,'edgecolor','none')
    ax.XAxis.Visible='off';
    ax.YAxis.Visible='off';
    title(strcat('\theta= ',num2str(data(n).theta)))
    axis equal
    axis tight
    shading interp
    caxis([0 3])
    set(gca,'position',pos{p})
    p=p+1;
end
print(gcf,strcat('velFieldsReconstruct',run),'-dpng','-r600')
end




















