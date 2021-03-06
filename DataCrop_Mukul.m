function DataCrop_Mukul(DataSet)
load("C:\Users\abber\Documents\School\Grad School\Winter 20\AMATH 582\Project\flow_data_confined_TSR_1_1.mat");
chord=4.06*10; %mm
theta_p=6; %preset pitch angle in degrees
data=eval(DataSet);

xgrid={};
ygrid={};
findNaN=[];
for i=1:length(data)
    %Rotate Field
    angle=-data(i).theta+theta_p;
    rot=[cosd(angle) -sind(angle); sind(angle) cosd(angle)];
    X=reshape(data(i).X,[],1); 
    Y=reshape(data(i).Y,[],1);
    XYrot=rot*[X';Y'];
    Xrot=reshape(XYrot(1,:),[data(i).Ny data(i).Nx]);
    Yrot=reshape(XYrot(2,:),[data(i).Ny data(i).Nx]);
    
    %Rotate Foil
    foilrot=rotate(data(i).foil,angle,[0,0]);
    
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
    [in,out]=inpolygon(Xrot,Yrot,cropbound.Vertices(:,1)...
        ,cropbound.Vertices(:,2));
    ind=logical(in+out);
    U_trans=data(i).u;
    V_trans=data(i).v;
    Vmag_trans=data(i).Umag;
    Vort_trans=data(i).vort;
    
    %Move all croped data so the bottom left corner of the bounding box is at (0,0) 
    Xcrop=Xrot(ind)-min(xbound);
    Ycrop=Yrot(ind)-min(ybound);
    U_crop=U_trans(ind);
    V_crop=V_trans(ind);
    Vmag_crop=Vmag_trans(ind);
    Vort_crop=Vort_trans(ind);
    %Shift the foil vertices as well
    foilrot.Vertices(:,1)=foilrot.Vertices(:,1)-min(xbound);
    foilrot.Vertices(:,2)=foilrot.Vertices(:,2)-min(ybound);
    data(i).nonInterp.foil=foilrot;
    
    %Save the coordinates inside the blade
    [inblade,outblade]=inpolygon(Xcrop,Ycrop,foilrot.Vertices(:,1)...
        ,foilrot.Vertices(:,2));
    indblade=logical(inblade+outblade);
    data(i).nonInterp.MaskCoordsX=Xcrop(indblade);
    data(i).nonInterp.MaskCoordsY=Ycrop(indblade);
    
    %Insert NaNs where the blade is
    U_crop(indblade)=NaN;
    V_crop(indblade)=NaN;
    Vmag_crop(indblade)=NaN;
    Vort_crop(indblade)=NaN;
    
    %Remove all NaN values
    data(i).nonInterp.Xcrop=Xcrop(indblade);
    data(i).nonInterp.Ycrop=Ycrop(indblade);
    data(i).nonInterp.U_crop=U_crop(indblade);
    data(i).nonInterp.V_crop=V_crop(indblade);
    data(i).nonInterp.Vmag_crop=Vmag_crop(indblade);
    data(i).nonInterp.Vort_crop=Vort_crop(indblade);
    
    %Interpolate to common grid
    N=52;
    xinterp=linspace(0,xwidth,N);
    yinterp=linspace(0,ywidth,N);
    data(i).Interp.foil=foilrot;
    [xgrid{i},ygrid{i}]=meshgrid(xinterp,yinterp);
    data(i).Interp.Xcrop=xgrid{i}(1,:);
    data(i).Interp.Ycrop=ygrid{i}(:,1);
    data(i).Interp.U_crop=griddata(Xcrop,Ycrop...
        ,U_crop,xgrid{i},ygrid{i});
    data(i).Interp.V_crop=griddata(Xcrop,Ycrop...
        ,V_crop,xgrid{i},ygrid{i});
    data(i).Interp.Vmag_crop=griddata(Xcrop,Ycrop...
        ,Vmag_crop,xgrid{i},ygrid{i});
    data(i).Interp.Vort_crop=griddata(Xcrop,Ycrop...
        ,Vort_crop,xgrid{i},ygrid{i});
    %find the NaN values in each field 
    findNaN(:,:,i)=isnan(data(i).Interp.Vort_crop); 
end 
for i=1:length(data)
    %Save the mask logicals
    findNaN=sum(findNaN,3);
    findNaN(findNaN~=0)=1;
    data(i).Interp.mask=...
        reshape(logical(findNaN),[],1);
    %Save the foil
    data(i).Interp.foil=foilrot;
    save(strcat(DataSet,' Crop'),'data')
end
end










