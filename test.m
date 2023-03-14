function test(app)
clc
appLock(app,'UNLOCK');
CalcPanel(app,'','','','off');


%step=266250;
%pathSc=MakePath(app,'SCF');
%scFnm=fullfile(pathSc,"LoopsSpaceCellsfile"+ step+".mat");
%sc=load(scFnm).sc;

%gr=grains('BASIC',step,'',app);
%sc = clusterVRangle(sc,gr);

%{
%step=app.N2EF.Value;
pathSc=MakePath(app,'SCF');
scFnm=fullfile(pathSc,"spaceCellsfile900000int5000.mat");
sc=load(scFnm).sc;
a=1;

%load gr
gr=grains('BASIC',step,'',app);
O=cat(1,sc.Loops.Order);
S=cat(1,sc.Loops.Size);

[av,od,~]= groupsummary(S, O, 'mean');

o=4+2*(1:400);
s=4+(1:400);
plot(o,s)
hold on
scatter(O,S)
plot(od,av)

figure
[av,od,~]= groupsummary(S./O, O, 'mean');
plot(o,s./o)
hold on
scatter(O,S./O)
plot(od,av)


lp=sc.Loops(find(O==36 & S./O<0.5,1));


d=cd;
d=fullfile(d,"loop");
mkdir(d)
%Identify all cells containing ALL these grains
grC=find(sum(ismember(sc.DelaunayT,lp.Grains),2)==4);
%get cells that DONOT make part of this cluster
grC=grC(~ismember(grC,lp.sCells)); 
cb=nchoosek(1:4,3);
for i=1:numel(grC)
    %get the grains
    tr=sc.DelaunayT(grC(i),:);
    fnm=fullfile(d,"Grains"+i+".vtk");
    vtkwrite(fnm,'unstructured_grid',...
            gr.Coord(tr,1),gr.Coord(tr,2),gr.Coord(tr,3),...
            'SCALARS','Radius',gr.Radius(tr),...
            'PRECISION',10);
    %transform then in cells
    tr=tr(cb);
    fnm=fullfile(d,"Cell"+i+".vtk");
    vtkwrite(fnm,'polydata',"TRIANGLE",...
        gr.Coord(:,1),gr.Coord(:,2),gr.Coord(:,3),tr);
    ismember(sort(tr,2),sort(lp.Vertices,2),'rows')
    
    
end

                
%{
v=lp.Vertices;
for i=1:length(v)
    vtc=lp.Vertices(i,:);
    fnm=fullfile(d,"Surface"+i+".vtk");
    vtkwrite(fnm,'polydata',"TRIANGLE",...
        gr.Coord(:,1),gr.Coord(:,2),gr.Coord(:,3),vtc);
    fnm=fullfile(d,"grs"+i+".vtk");
    vtkwrite(fnm,'unstructured_grid',...
            gr.Coord(vtc,1),gr.Coord(vtc,2),gr.Coord(vtc,3),...
            'SCALARS','Radius',gr.Radius(vtc),...
            'PRECISION',10);
end
%}
fnm=fullfile(d,"SurfaceTot.vtk");
vtkwrite(fnm,'polydata',"TRIANGLE",...
    gr.Coord(:,1),gr.Coord(:,2),gr.Coord(:,3),lp.Vertices);
v=unique(lp.Vertices);
fnm=fullfile(d,"grsTot.vtk");
vtkwrite(fnm,'unstructured_grid',...
            gr.Coord(v,1),gr.Coord(v,2),gr.Coord(v,3),...
            'SCALARS','Radius',gr.Radius(v),...
            'PRECISION',10);
%}

%plotClusterFC(app,'','')
% step=255000;
% gr=grains('basic',step,'',app);
% pathSc=[app.SavePath '/MatlabResultsFile/SpaceCellFiles/'];
% scFnm=[pathSc char("LoopsSpaceCellsfile"+step+".mat")];
% sc=load(scFnm).sc;
% lS = clusterStress(sc,step,app,gr);
% 
% save('sctestfile.mat','sc','-v7.3');
% save('grN.mat','gr','-v7.3');

% grNew=grains('basic',step,'',app);
% grOld=grains('basic',step-5000,'',app);
% sc = spaceCellSystem('TEST',step,grNew,grOld,app,'');
% save('sctestfile.mat','sc','-v7.3');
% save('grN.mat','grNew','-v7.3');
% save('grO.mat','grOld','-v7.3');

% load('spaceCellsfile250000int5000.mat')
% gr=grains('basic',250000,'',app)
% 
% vr=4/3*pi()*gr.Radius.^3./sc.GrainVolume;
% sum(vr>1)

% step=150000;
% gr=grains('STENSOR',step,'',app);
% sc = spaceCellSystem('TEST',step,gr,app,'');
% 
% save('sctestfile.mat','sc','-v7.3');
% save('grtestfile.mat','gr','-v7.3');
% vs=4/3*pi()*sc.Radius(1:NbG).^3;
% e=(sc.GrainVolume-vs)./vs;
% mnID=find(e=min(e));
% mxID=find(e=max(e));
% 
% [X,Y,Z]=sphere;
% %first grain
% f=figure;ax=axes(f);hold(ax,'on');
% r=sc.Radius(mnID);
% cel=cell2mat(vertexAttachments(sc.DelaunayT,mnID)); %ids of tetrahedron touching grain i
% cel=cel(ismember(cel,sc.GoodCells));
% cc=[sc.sCells(cel).Center]'; %get property
% cc=cc(cc(:,1)==0 | cc(:,1)==mnID | cc(:,2)==mnID, 3:end); %check correct positions not the ids
% [b,vol]=boundary(cc) %calculate volume
%                 
% surf(ax,X*r+gr.Coord(mnID,1),Y*r+gr.Coord(mnID,2),Z*r+gr.Coord(mnID,3))
% axis equal
% trisurf(ax,b,cc(:,1),cc(:,2),cc(:,3),'Facecolor','red','FaceAlpha',0.1)
% 
% %second grain
% g=figure;ax=axes(g);hold(ax,'on');
% r=sc.Radius(mxID);
% cel=cell2mat(vertexAttachments(sc.DelaunayT,mxID)); %ids of tetrahedron touching grain i
% cel=cel(ismember(cel,sc.GoodCells));
% cc=[sc.sCells(cel).Center]'; %get property
% cc=cc(cc(:,1)==0 | cc(:,1)==mxID | cc(:,2)==mxID, 3:end); %check correct positions not the ids
% [b,vol]=boundary(cc) %calculate volume
% surf(ax,X*r+gr.Coord(mxID,1),Y*r+gr.Coord(mxID,2),Z*r+gr.Coord(mxID,3))
% trisurf(ax,b,cc(:,1),cc(:,2),cc(:,3),'Facecolor','red','FaceAlpha',0.1)
% axis equal

% sum(sc.GrainVolume)
% N1=app.N1EF.Value;
% N2=app.N2EF.Value;
% interval=app.CalcInt.Value;
% nbFiles=ceil((N2-N1)/interval+1);
% V=zeros(nbFiles,5);
% for i=1:nbFiles
%     fprintf('Doing %d of %d\n',i,nbFiles)
%     step=min(N1+interval*(i-1),N2);
%     gr=grains('BASIC',step,'',app); %calculate loops
%     V(i,:) =[step avrgA(gr,step,app)];
% end
% 
% f=figure;ax=axes(f);hold(ax,'on');
% plot(ax,V(:,1),V(:,2))
% title(ax,"Max")
% f=figure;ax=axes(f);hold(ax,'on');
% plot(ax,V(:,1),V(:,3))
% title(ax,"Min")
% f=figure;ax=axes(f);hold(ax,'on');
% plot(ax,V(:,1),V(:,4))
% title(ax,"Mean")
% f=figure;ax=axes(f);hold(ax,'on');
% plot(ax,V(:,1),V(:,5))
% title(ax,"Md")

end

%{
%Loops order change

N1=app.N1EF.Value;
N2=app.N2EF.Value;
interval=app.CalcInt.Value;
nbFiles=ceil((N2-N1)/interval+1);

pathSc=MakePath(app,'SCF');
for j=1:nbFiles
    step=min(N1+interval*(j-1),N2);
    fprintf("Executing step "+step+"\n") 
    %load sc
    scFnm=fullfile(pathSc,"LoopsSpaceCellsfile"+step+".mat");
    sc=load(scFnm).sc;
    for i=1:numel(sc.Loops)
        %For each cluster
        lp=sc.Loops(i);
        %Identify all surfaces
        sf=sc.DelaunayT(lp.sCells,:);
        sf=cat(3, [sf(:,1) sf(:,2) sf(:,3)],...
            [sf(:,1) sf(:,2) sf(:,4)],...
            [sf(:,1) sf(:,3) sf(:,4)],...
            [sf(:,2) sf(:,3) sf(:,4)]);
        sf=sort(cat(1,sf(:,:,1),sf(:,:,2),sf(:,:,3),sf(:,:,4)),2);
        %get repeating surfaces
        [un,~,cnt]=unique(sf,'rows'); %get unique
        rpCE=accumarray(cnt,1);
        sf=un(rpCE==1,:); %keep not repeating values
        sc.Loops(i).Vertices=sf;
        sc.Loops(i).Order=size(sf,1);
    end
	save(scFnm,'sc','-v7.3');
end
%}