clear all, close all
%% Setup
%Choose if doing plain PCA or rPCA and with uv or with vmag
plain=1; %1 for plain PCA 0 for rPCA
uv=0; %1 for uv and 0 for vmag
workingfolder= 'C:\Users\abber\Documents\School\Grad School\Winter 20\AMATH 582\Project';
chord=4.06*10; %chord length in mm

if plain
    t=0; %truncation
    run='Run1'; %what to append to all plot saving so things don't get overwritten between data sets
    %load data
    load(fullfile(workingfolder,'up1_1 Crop.mat'))
    rotFields = fieldnames(data);
    nRot = length(rotFields);
    %prep for PCA
    nx = 52; %size of the cropped fields that go into the SVD
    ny = 52;
    if uv
        Y = zeros(nx*ny*2, nRot);
    else
        Y = zeros(nx*ny, nRot);
    end
    xcrop=data.(rotFields{1}).interp.xcrop;
    ycrop=data.(rotFields{1}).interp.ycrop;
end

%% rPCA
if ~plain
    t=2; %truncation
    run='Run1'; %what to append to all plot saving so things don't get overwritten between data sets
    
    % load data, initialize
    load(fullfile(workingfolder,'up1_1FULL.mat'))
    rotFields = fieldnames(up1_1FULL);
    nRot = length(rotFields);
    %prep for PCA
    nx = 52; %size of the cropped fields that go into the SVD
    ny = 52;
    if uv
        Y = zeros(nx*ny*2, nRot);
    else
        Y = zeros(nx*ny, nRot);
    end
    
    % keep these variables outside loop if same for all rotations
    nxFull = up1_1FULL.deg_30.nx;
    nyFull = up1_1FULL.deg_30.ny;
    
    for iRot = 1:nRot
        u = up1_1FULL.(rotFields{iRot}).u;
        v = up1_1FULL.(rotFields{iRot}).v;
        vmag = up1_1FULL.(rotFields{iRot}).vmag;
        mask_ind = up1_1FULL.(rotFields{iRot}).mask_inds;
        lambda = 5; % choose your sparsity constant value (1 is often a good starting point)
        tol = 1e-7; % set your tolerance
        maxIter = 1000; % set your maximum number of iterations
        if iRot>(7+t)
            fl=1; %because of data rotation we need to flip the mask left to right
        else
            fl=0;
        end
        
        % average all columns of L,N to get "rPCA cleaned" phase average
        if uv
            [up1_1FULL.(rotFields{iRot}).L, up1_1FULL.(rotFields{iRot}).N,mask_log(:,:,iRot)] ...
                = rPCA_main(u,v,nxFull,nyFull, mask_ind, lambda, tol, maxIter,fl);
            up1_1FULL.(rotFields{iRot}).Lr = mean(up1_1FULL.(rotFields{iRot}).L,2);
            up1_1FULL.(rotFields{iRot}).Nr = mean(up1_1FULL.(rotFields{iRot}).N,2);
            % put mask back in
            [up1_1FULL.(rotFields{iRot}).uL,up1_1FULL.(rotFields{iRot}).vL] = ...
                unstackPCA(up1_1FULL.(rotFields{iRot}).Lr,nxFull,nyFull,mask_log(:,:,iRot),1);
            [up1_1FULL.(rotFields{iRot}).uN,up1_1FULL.(rotFields{iRot}).vN] = ...
                unstackPCA(up1_1FULL.(rotFields{iRot}).Nr,nxFull,nyFull,mask_log(:,:,iRot),1);
        else
            [up1_1FULL.(rotFields{iRot}).L_vmag, up1_1FULL.(rotFields{iRot}).N_vmag,mask_log(:,:,iRot)] ...
                = rPCA_main(vmag,[],nxFull,nyFull, mask_ind, lambda, tol, maxIter,fl);
            up1_1FULL.(rotFields{iRot}).Lr_vmag = mean(up1_1FULL.(rotFields{iRot}).L_vmag,2);
            up1_1FULL.(rotFields{iRot}).Nr_vmag = mean(up1_1FULL.(rotFields{iRot}).N_vmag,2);
            % put mask back in
            [up1_1FULL.(rotFields{iRot}).vmagL] = ...
                unstackPCA(up1_1FULL.(rotFields{iRot}).Lr_vmag,nxFull,nyFull,mask_log(:,:,iRot),0);
            [up1_1FULL.(rotFields{iRot}).vmagN] = ...
                unstackPCA(up1_1FULL.(rotFields{iRot}).Nr_vmag,nxFull,nyFull,mask_log(:,:,iRot),0);
        end
    end
    
    % save post rPCA structure
    save(fullfile(workingfolder,'up1_1FULL rPCA'),'up1_1FULL')
    
    %Calculate Variance and plot
    if ~uv
        fvar=figure;
        set(gcf,'position',[633.5714  313.5714  806.8286  406.4286])
        
        subplot(1,3,1)
        pcolor(up1_1FULL.(rotFields{13}).x/chord,up1_1FULL.(rotFields{13}).y/chord...
            ,up1_1FULL.(rotFields{13}).vmag(:,:,1)')
        shading flat
        axis equal
        axis tight
        hold on
        plot(up1_1FULL.(rotFields{13}).foil,'facecolor',[0 0 0],'facealpha',0.5...
            ,'edgecolor','none')
        title('Single Frame')
        
        subplot(1,3,2)
        var_plain=nanmean((up1_1FULL.(rotFields{13}).vmag...
            -nanmean(up1_1FULL.(rotFields{13}).vmag,3)).^2,3);
        pcolor(up1_1FULL.(rotFields{13}).x/chord,up1_1FULL.(rotFields{13}).y/chord...
            ,(var_plain/max(var_plain,[],'all'))')
        shading flat
        axis equal
        axis tight
        hold on
        plot(up1_1FULL.(rotFields{13}).foil,'facecolor',[0 0 0],'facealpha',0.5...
            ,'edgecolor','none')
        title('\langle v''_{mag}^2 \rangle')
        
        var=mean((up1_1FULL.(rotFields{13}).L_vmag...
            -mean(up1_1FULL.(rotFields{13}).L_vmag,2)).^2,2);
        [var_vmag] = unstackPCA(var,nxFull,nyFull,mask_log(:,:,13),0);
        
        subplot(1,3,3)
        pcolor(up1_1FULL.(rotFields{13}).x/chord,up1_1FULL.(rotFields{13}).y/chord...
            ,(var_vmag/max(var_vmag,[],'all'))')
        shading flat
        axis equal
        axis tight
        hold on
        plot(up1_1FULL.(rotFields{13}).foil,'facecolor',[0 0 0],'facealpha',0.5...
            ,'edgecolor','none')
        title('\langle v''_{mag}^2 \rangle after rPCA')
        c=colorbar;
        set(c,'position',[0.9226    0.2804    0.0263    0.4742])
        sgtitle(strcat('\theta = '...
            , rotFields{13}(strfind(rotFields{13},'_')+1:end),char(176)))
        print(gcf,'rPCA variance','-dpng','-r600')
    end
    
    % Crops data for PCA
    DataCrop_rPCA(fullfile(workingfolder,'up1_1FULL rPCA'),t,0,uv)
    load(fullfile(workingfolder,'up1_1FULL rPCA Crop.mat'))
    rotFields = fieldnames(data);
    nRot = length(rotFields);
    xcrop=data.(rotFields{1}).interp.xcrop;
    ycrop=data.(rotFields{1}).interp.ycrop;

end
%% Form Y matrix

for iRot = 1:nRot
    if uv
        u = reshape(data.(rotFields{iRot}).interp.u_crop, [nx*ny 1]);
        v = reshape(data.(rotFields{iRot}).interp.v_crop, [nx*ny 1]);
        Y(:,iRot) = [u;v];
    else
        Y(:,iRot) = reshape(data.(rotFields{iRot}).interp.vmag_crop, [nx*ny 1]);
    end
end
%% Vel field plotting -rPCA
if ~plain
    figure('DefaultAxesFontsize', 16)
    ax=gca;
    pcolor(data.deg_137.interp.xcrop,data.deg_137.interp.ycrop,data.deg_137.interp.vmag_cropL)
    hold on
    plot(data.deg_137.interp.foil,'facecolor',[0 0 0],'facealpha',0.5...
        ,'edgecolor','none')
    axis equal
    axis tight
    shading interp
    caxis([-1,2.5])
    c=colorbar;
    c.FontSize=12;
    set(get(c,'title'),'string','$(\frac{V_{mag}}{U_\infty})$','interpreter','latex');
    xlabel('x/c')
    ylabel('y/c')
    title(['L   \lambda = ' num2str(lambda) ])
    
    figure('DefaultAxesFontsize', 16)
    ax=gca;
    pcolor(data.deg_137.interp.xcrop,data.deg_137.interp.ycrop,data.deg_137.interp.vmag_cropN)
    hold on
    plot(data.deg_137.interp.foil,'facecolor',[0 0 0],'facealpha',0.5...
        ,'edgecolor','none')
    axis equal
    axis tight
    shading interp
    caxis([-1,2.5])
    c=colorbar;
    c.FontSize=12;
    set(get(c,'title'),'string','$(\frac{V_{mag}}{U_\infty})$','interpreter','latex');
    xlabel('x/c')
    ylabel('y/c')
    title(['N   \lambda = ' num2str(lambda) ])
end
%% Vel field plotting- plain
if plain
figure
set(gcf,'position',1.0e+03 *[0.0016    0.2079    1.4600    0.5120])
[ha, pos]= tight_subplot(2,4,[0 0],[.01 .01],[.01 .01]);
p=1;
for n=1:ceil(nRot/8):nRot
    axes(ha(p))
    ax=gca;
    set(gca,'Visible','off')
    tmpU=data.(rotFields{n}).interp.u_crop;
    tmpV=data.(rotFields{n}).interp.v_crop;
    tmpVmag=data.(rotFields{n}).interp.vmag_crop;
    pcolor(data.(rotFields{n}).interp.xcrop,data.(rotFields{n}).interp.ycrop,tmpVmag)
    hold on
    quiver(data.(rotFields{n}).interp.xcrop(1:2:end)...
        ,data.(rotFields{n}).interp.ycrop(1:2:end),tmpU(1:2:end,1:2:end)...
        ,tmpV(1:2:end,1:2:end),2,'k')
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
c=colorbar;
c.FontSize=12;
set(c,'position',[0.7585    0.0513    0.0180    0.4074])
set(get(c,'title'),'string','$(\frac{V_{mag}}{U_\infty})$','interpreter','latex');
axes(ha(end))
set(gca,'Visible','off')
print(gcf,strcat('velFields',run),'-dpng','-r600')
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

[U,S,V] = svd(Yms,'econ');
sig=diag(S);
energy=sig/sum(sig)*100;
for n=1:length(sig)
    energytotal(n)=sum(energy(1:n));
end

%% Plot singular values & cumulative energy
figure
set(gcf,'position',[251.2857  403.9000  485.1429  316.1000])
left_color = [0.9153    0.2816    0.2878];
right_color = [0 .5 .5];
set(gcf,'defaultAxesColorOrder',[left_color; right_color]);
yyaxis left
set(gca,'fontsize',14)
a=plot(energy,'o','markerfacecolor',left_color);
hold on

ind=find(energytotal>=90,1);
plot(ind,energy(ind),'o','markerfacecolor',[0 0 0],'markersize',12)
ylabel('% of energy')
axis tight
ylim([0 35])
yyaxis right
set(gca,'fontsize',14)
plot(energytotal,'o','markerfacecolor',right_color)
hold on
plot(linspace(1,length(sig),100),90*ones(1,100),'--','color',right_color,'linewidth',1.75)
ylabel('Cumulative Energy %')
xlabel('mode')
axis tight
grid on
print(gcf,strcat('Energy Cropped',run),'-dpng','-r600')

%% show modes
%put NaNs back in
Ureplaced=zeros(size(Y));
Ureplaced(nanRow>th,:)=NaN;
Ureplaced(~(nanRow>th),:)=U;

f=figure;
set(f,'position',[ 331.8571   13.5714  760.0000  706.4286])
[ha, pos]= tight_subplot(3,2,[.05 .15],[.1 .01],[.1 .1]);
f2=figure;
set(f2,'position',[ 331.8571   13.5714  760.0000  706.4286])
[ha2, pos2]= tight_subplot(3,2,[.05 .15],[.1 .01],[.1 .1]);
j=1;
for k = [1:5,ind]
    f
    axes(ha(j))
    %Put NaNs back in
    if uv
        mode=Ureplaced(1:end/2,k);
        mode=mode/max(abs(mode),[],'all');
        mode2=Ureplaced(end/2+1:end,k);
        mode2=mode2/max(abs(mode2),[],'all');
        titletext='u modes';
        titletext2='v modes';
    else
        mode=Ureplaced(:,k);
        mode=mode/max(abs(mode),[],'all');
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
    caxis([-1 1])
    set(gca,'fontsize',12)
    set(gca,'position',pos{j})
    
    if uv
        f2
        axes(ha2(j))
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
        caxis([-1 1])
        set(gca,'fontsize',12)
        set(gca,'position',pos2{j})
    end
    j=j+1;
end
print(f,strcat(titletext,run),'-dpng','-r600')
if uv
    print(f2,strcat(titletext2,run),'-dpng','-r600')
end

%% plot how modes evolve with angle (right singular vectors)
angles = [48 57 66 75 84 93 98 101 110 119 137 146 155 164];

figure
p1 = plot(angles, V(:,1),'k-','Linewidth',[2]) ;
hold on
p2 = plot(angles, V(:,2),'k--','Linewidth',[2]) ;
plot(angles, V(:,3),'k:','Linewidth',[2])
legend('mode 1', 'mode 2', 'mode 3', 'Location', 'northwest')
ylabel('V')
xlabel('angle, degrees')
uistack(p2, 'top')
uistack(p1, 'top')
hold off
