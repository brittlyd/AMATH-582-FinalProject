function DataCrop(DataSet,truncateStart,truncateEnd)
load(DataSet);
chord=4.06*10; %mm
theta_p=6; %preset pitch angle in degrees
data=eval(DataSet);
names=fieldnames(data);

%truncate the data from the begining, end or both
if ~isempty(truncateStart)
    a=truncateStart;
    data=rmfield(data,names(1:a));
else
    a=1;
end

if ~isempty(truncateEnd)
    b=truncateEnd;
    data=rmfield(data,names(end-b+1:end));
else
    b=0;
end

xgrid={};
ygrid={};
findNaN=[];
for i=a+1:length(names)-b
    %Rotate Field
    angle=-str2double(names{i}(strfind(names{i},'_')+1:end))+theta_p;
    rot=[cosd(angle) -sind(angle); sind(angle) cosd(angle)];
    %move the centershaft to the origin to properly rotate
    [X,Y]=meshgrid(data.(names{i}).X-data.(names{i}).shaftX*chord...
        ,data.(names{i}).Y-data.(names{i}).shaftY*chord);
    X=reshape(X,[],1); 
    Y=reshape(Y,[],1);
    XYrot=rot*[X';Y'];
    Xrot=reshape(XYrot(1,:),[data.(names{i}).Ny data.(names{i}).Nx])...
        +data.(names{i}).shaftX*chord; %shift the centershaft back from the origin
    Yrot=reshape(XYrot(2,:),[data.(names{i}).Ny data.(names{i}).Nx])...
        +data.(names{i}).shaftY*chord;
    
    %Rotate Foil
    foilrot=rotate(data.(names{i}).foil,angle...
        ,[data.(names{i}).shaftX,data.(names{i}).shaftY]);
    
    %Crop Field
    [xbound,ybound]=boundingbox(foilrot);
    xbound(1)=xbound(1)+chord/15; %Grow the crop boundary and change Aspect Ratio
    xbound(2)=xbound(2)-chord/30;
    ybound(1)=ybound(1)+chord/60;
    ybound(2)=ybound(2)-chord/30;
    %Produce vericies to create polyshape
    [xbound,ybound]=meshgrid(xbound,ybound);
    xbound=reshape(xbound,[],1);
    xbound(3:end)=flipud(xbound(3:end));
    ybound=reshape(ybound,[],1);
    ybound(3:end)=flipud(ybound(3:end));
    cropbound=polyshape(xbound',ybound');
    xwidth=max(xbound)-min(xbound);
    ywidth=max(ybound)-min(ybound);
    %Determine Coordinates within the polyshape
    [in,out]=inpolygon(Xrot/chord,Yrot/chord,cropbound.Vertices(:,1)...
        ,cropbound.Vertices(:,2));
    ind=logical(in+out);
    U_trans=data.(names{i}).U_Avg';
    V_trans=data.(names{i}).V_Avg';
    Vmag_trans=data.(names{i}).Vmag_Avg';
    Vort_trans=data.(names{i}).Vort_Avg';
    Swirl_trans=data.(names{i}).Swirl_Avg';
    
    %Move all croped data so the bottom left corner of the bounding box is at (0,0) 
    Xcrop=Xrot(ind)-min(xbound)*chord;
    Ycrop=Yrot(ind)-min(ybound)*chord;
    U_crop=U_trans(ind);
    V_crop=V_trans(ind);
    Vmag_crop=Vmag_trans(ind);
    Vort_crop=Vort_trans(ind);
    Swirl_crop=Swirl_trans(ind);
    %Shift the foil vertices as well
    foilrot.Vertices(:,1)=foilrot.Vertices(:,1)-min(xbound);
    foilrot.Vertices(:,2)=foilrot.Vertices(:,2)-min(ybound);
    data.(names{i}).nonInterp.foil=foilrot;
    
    %Save the coordinates of the mask
    data.(names{i}).nonInterp.MaskCoordsX=Xcrop(isnan(Vort_crop));
    data.(names{i}).nonInterp.MaskCoordsY=Ycrop(isnan(Vort_crop));
    
    %Remove all NaN values
    data.(names{i}).nonInterp.Xcrop=Xcrop(~isnan(Vort_crop));
    data.(names{i}).nonInterp.Ycrop=Ycrop(~isnan(Vort_crop));
    data.(names{i}).nonInterp.U_crop=U_crop(~isnan(Vort_crop));
    data.(names{i}).nonInterp.V_crop=V_crop(~isnan(Vort_crop));
    data.(names{i}).nonInterp.Vmag_crop=Vmag_crop(~isnan(Vort_crop));
    data.(names{i}).nonInterp.Vort_crop=Vort_crop(~isnan(Vort_crop));
    data.(names{i}).nonInterp.Swirl_crop=Swirl_crop(~isnan(Vort_crop));
    
    %Interpolate to common grid
    N=52;
    xinterp=linspace(0,xwidth,N);
    yinterp=linspace(0,ywidth,N);
    data.(names{i}).Interp.foil=foilrot;
    [xgrid{i},ygrid{i}]=meshgrid(xinterp,yinterp);
    data.(names{i}).Interp.Xcrop=xgrid{i}(1,:);
    data.(names{i}).Interp.Ycrop=ygrid{i}(:,1);
    data.(names{i}).Interp.U_crop=griddata(Xcrop/chord,Ycrop/chord...
        ,U_crop,xgrid{i},ygrid{i});
    data.(names{i}).Interp.V_crop=griddata(Xcrop/chord,Ycrop/chord...
        ,V_crop,xgrid{i},ygrid{i});
    data.(names{i}).Interp.Vmag_crop=griddata(Xcrop/chord,Ycrop/chord...
        ,Vmag_crop,xgrid{i},ygrid{i});
    data.(names{i}).Interp.Vort_crop=griddata(Xcrop/chord,Ycrop/chord...
        ,Vort_crop,xgrid{i},ygrid{i});
    data.(names{i}).Interp.Swirl_crop=griddata(Xcrop/chord,Ycrop/chord...
        ,Swirl_crop,xgrid{i},ygrid{i});
    %find the NaN values in each field using vorticity because it has the
    %most NaNs
    findNaN(:,:,i)=isnan(data.(names{i}).Interp.Vort_crop); 
    %Save the coordinates of the mask
    data.(names{i}).Interp.MaskCoordsX=...
        reshape(xgrid{i}(logical(findNaN(:,:,i))),[],1);
    data.(names{i}).Interp.MaskCoordsY=...
        reshape(ygrid{i}(logical(findNaN(:,:,i))),[],1);
    %Save the foil
    data.(names{i}).Interp.foil=foilrot;
end

%Only keep the points that data exits at in every case
NaNsum=sum(findNaN,3);
for i=a+1:length(names)-b
    data.(names{i}).InterpCommon.Xcrop=data.(names{i}).Interp.Xcrop;
    data.(names{i}).InterpCommon.Ycrop=data.(names{i}).Interp.Ycrop;
    data.(names{i}).InterpCommon.U_crop=data.(names{i}).Interp.U_crop;
    data.(names{i}).InterpCommon.V_crop=data.(names{i}).Interp.V_crop;
    data.(names{i}).InterpCommon.Vmag_crop=data.(names{i}).Interp.Vmag_crop;
    data.(names{i}).InterpCommon.Vort_crop=data.(names{i}).Interp.Vort_crop;
    data.(names{i}).InterpCommon.Swirl_crop=data.(names{i}).Interp.Swirl_crop;
    
    figure(3); clf(3)
    subplot(1,2,1)
    pcolor(data.(names{i}).InterpCommon.Xcrop...
        ,data.(names{i}).InterpCommon.Ycrop,data.(names{i}).InterpCommon.Vort_crop)
    shading interp
    axis equal
    axis tight
    
    data.(names{i}).InterpCommon.U_crop(NaNsum~=0)=NaN;
    data.(names{i}).InterpCommon.V_crop(NaNsum~=0)=NaN;
    data.(names{i}).InterpCommon.Vmag_crop(NaNsum~=0)=NaN;
    data.(names{i}).InterpCommon.Vort_crop(NaNsum~=0)=NaN;
    data.(names{i}).InterpCommon.Swirl_crop(NaNsum~=0)=NaN;
    %Save the coordinates of the mask
    data.(names{i}).InterpCommon.MaskCoordsX=xgrid{i}(NaNsum~=0);
    data.(names{i}).InterpCommon.MaskCoordsY=ygrid{i}(NaNsum~=0);
    %Save the foil
    data.(names{i}).InterpCommon.foil=foilrot;
    
    subplot(1,2,2)
    pcolor(data.(names{i}).InterpCommon.Xcrop...
        ,data.(names{i}).InterpCommon.Ycrop,data.(names{i}).InterpCommon.Vort_crop)
    hold on
    plot(data.(names{i}).InterpCommon.foil)
    shading interp
    axis equal
    axis tight
    drawnow
    
    save(strcat(DataSet,' Crop'),'data')
end
end










