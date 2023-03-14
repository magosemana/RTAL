function plotClusterFC(app,step,clOrder)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if isempty(step);step=app.N2EF.Value;end
gr=grains('FORCECHAIN',step,'',app); %get grains
fnm=fullfile(MakePath(app,'SCF'),char("LoopsSpaceCellsfile"+step+".mat"));
sc=load(fnm).sc; %load Clusters

if isnumeric(clOrder)
    clID=find(cat(1,sc.Loops.Order)==clOrder,1);
    if isempty(clID)
       fprintf('No cluster of the asked order\n') 
    end
else
    [~,clID]=max(cat(1,sc.Loops.Order));
end
%clVRplot(sc,clID,gr)

clFCplot(app,sc,clID,gr)
end


function clVRplot(sc,clID,gr) %#ok<*DEFNU>
gC = goodCell(sc);
if isempty(sc.Loops(1).Volume)
    vl = perCellVolume(sc,gr,gC);
end
%Get a base sphere
bs=baseSphere();

%Turn base sphere vector into a 3D matrix, per grain
bsG=bs.*ones(1,1,sc.Loops(clID).Size);
%transforfm each page of the 3D matrix in a sphere
%represent each grain of the cluster
bsG=pagemtimes(bsG,permute(gr.Radius(sc.Loops(clID).Grains),[3,2,1]))+...
    permute(gr.Coord(sc.Loops(clID).Grains,:),[3,2,1]);
bsG=reshape(permute(bsG,[1,3,2]),[],3,1);

%Get the IDs of the cells that contain these grains and
%belong to the cluster l
ID = pointLocation(sc.DelaunayT,bsG);
lCl=ismember(ID,sc.Loops(clID).sCells);

%Need to add the point where the edges cross the grain boundary for more
%precision Vertices contain all exterior faces. Mix them to get get all
%exterior edges
ed=sc.Loops(clID).Vertices;
ed=unique(sort([ed(:,1) ed(:,2);ed(:,1) ed(:,3);ed(:,2) ed(:,3)],2),'rows');
edVec=gr.Coord(ed(:,2),:)-gr.Coord(ed(:,1),:);
edVec=edVec./vecnorm(edVec,2,2);
bsG2=[gr.Coord(ed(:,1),:)+edVec.*gr.Radius(ed(:,1),:);...
    gr.Coord(ed(:,2),:)-edVec.*gr.Radius(ed(:,2),:)];

%Calculate the volume of these points
bsG=[bsG(lCl,:);bsG2];
[bd,vv] = boundary(bsG);

%if total volume is empty calculate it
if isempty(sc.Loops(clID).Volume)
    sc.Loops(clID).Volume=sum(vl(ismember(gC,sc.Loops(clID).sCells)));
end
%calculate VR and save it
sc.Loops(clID).VoidRatio=vv/sc.Loops(clID).Volume;

figure;
trisurf(bd,bsG(:,1),bsG(:,2),bsG(:,3),'Facecolor','blue')
hold on
cel=sc.Loops(clID).sCells;

for i=1:numel(cel)
    g=gr.Coord(sc.DelaunayT(cel(i),:),:);
    [cbd,~] = boundary(g);
    trisurf(cbd,g(:,1),g(:,2),g(:,3),'Facecolor','red','FaceAlpha',0.1)
end


end

function clFCplot(app,sc,clID,gr)
clGr=sc.Loops(clID).Grains;
fc=gr.ForceChains;

% identify forcechains cointaing grs clGr
fcIds=zeros(numel(fc),1);
for i=1:numel(gr.ForceChains)
    chk=ismember(fc(i).IDs,clGr);
    if sum(chk)>0;fcIds(i)=1;end
end
fcIds=find(fcIds);
path=fullfile(MakePath(app),"Cluster"+clID);
if exist(path,'dir')==0;mkdir(path);end
%draw fc lines
vtkwrite(fullfile(path,"FcLines-All.vtk"),'polydata',...
    "LINESB",gr.Coord(:,1),gr.Coord(:,2),gr.Coord(:,3),cat(1,fc(fcIds).Lines));
%draw grains
id=unique(cat(1,fc(fcIds).IDs));
vtkwrite(fullfile(path,"FcGrains-All.vtk"),'unstructured_grid',...
    gr.Coord(id,1),gr.Coord(id,2),gr.Coord(id,3),...
    'SCALARS','Radius',gr.Radius(id),...
    'PRECISION',10);
for i=1:numel(fcIds)
    vtkwrite(fullfile(path,"FcLines"+i+".vtk"),'polydata',...
        "LINESB",gr.Coord(:,1),gr.Coord(:,2),gr.Coord(:,3),cat(1,fc(fcIds(i)).Lines));
    %draw grains
    id=unique(cat(1,fc(fcIds(i)).IDs));
    vtkwrite(fullfile(path,"FcGrains-All"+i+".vtk"),'unstructured_grid',...
        gr.Coord(id,1),gr.Coord(id,2),gr.Coord(id,3),...
        'SCALARS','Radius',gr.Radius(id),...
        'PRECISION',10);
end
%draw cluster
drawClst=cat(1,sc.Loops(clID).Vertices);
vtkwrite(fullfile(path,"Cluster.vtk"),'polydata',"TRIANGLE",...
    gr.Coord(:,1),gr.Coord(:,2),gr.Coord(:,3),drawClst);
%draw cluster's grains
id=sc.Loops(clID).Grains;
vtkwrite(fullfile(path,"ClusterGrains.vtk"),'unstructured_grid',...
    gr.Coord(id,1),gr.Coord(id,2),gr.Coord(id,3),...
    'SCALARS','Radius',gr.Radius(id),...
    'PRECISION',10);

end