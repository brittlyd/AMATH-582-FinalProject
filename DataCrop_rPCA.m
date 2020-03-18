function DataCrop_rPCA(DataSet,truncateStart,truncateEnd,uv)
load(DataSet);
data=up1_1FULL;
chord=4.06*10; %mm
theta_p=6; %preset pitch angle in degrees
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
    [X,Y]=meshgrid(data.(names{i}).x-data.(names{i}).shaftx*chord...
        ,data.(names{i}).y-data.(names{i}).shafty*chord);
    X=reshape(X,[],1); 
    Y=reshape(Y,[],1);
    XYrot=rot*[X';Y'];
    Xrot=reshape(XYrot(1,:),[data.(names{i}).ny data.(names{i}).nx])...
        +data.(names{i}).shaftx*chord; %shift the centershaft back from the origin
    Yrot=reshape(XYrot(2,:),[data.(names{i}).ny data.(names{i}).nx])...
        +data.(names{i}).shafty*chord;
    
    %Rotate Foil
    foilrot=rotate(data.(names{i}).foil,angle...
        ,[data.(names{i}).shaftx,data.(names{i}).shafty]);
    
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
    if uv
    U_trans=data.(names{i}).uL';
    V_trans=data.(names{i}).vL';
    else
    Vmag_trans=data.(names{i}).vmagL';
    end
    
    %Move all croped data so the bottom left corner of the bounding box is at (0,0) 
    Xcrop=Xrot(ind)-min(xbound)*chord;
    Ycrop=Yrot(ind)-min(ybound)*chord;
    if uv
    U_crop=U_trans(ind);
    V_crop=V_trans(ind);
    else
    Vmag_crop=Vmag_trans(ind);
    end

    %Shift the foil vertices as well
    foilrot.Vertices(:,1)=foilrot.Vertices(:,1)-min(xbound);
    foilrot.Vertices(:,2)=foilrot.Vertices(:,2)-min(ybound);
    
    %Interpolate to common grid
    N=52;
    xinterp=linspace(0,xwidth,N);
    yinterp=linspace(0,ywidth,N);
    data.(names{i}).interp.foil=foilrot;
    [xgrid{i},ygrid{i}]=meshgrid(xinterp,yinterp);
    data.(names{i}).interp.xcrop=xgrid{i}(1,:);
    data.(names{i}).interp.ycrop=ygrid{i}(:,1);
    if uv
    data.(names{i}).interp.u_crop=griddata(Xcrop/chord,Ycrop/chord...
        ,U_crop,xgrid{i},ygrid{i});
    data.(names{i}).interp.v_crop=griddata(Xcrop/chord,Ycrop/chord...
        ,V_crop,xgrid{i},ygrid{i});
    else
    data.(names{i}).interp.vmag_crop=griddata(Xcrop/chord,Ycrop/chord...
        ,Vmag_crop,xgrid{i},ygrid{i});
    end
    data.(names{i}).interp.foil=foilrot;
end
save(strcat(DataSet,' Crop'),'data')
end










