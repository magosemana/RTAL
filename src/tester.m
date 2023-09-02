% clear
% clc
% load('grN.mat')
% load('sctestfile.mat')
% 

%{
NbG=gr.Nb;
vs=4/3*pi()*sc.Radius(1:NbG).^3;
e=(sc.GrainVolume-vs)./vs;
mn=min(e);
mx=max(e);
mnID=find(e==mn);
mxID=find(e==mx);

[X,Y,Z]=sphere;
%first grain
r=sc.Radius(mnID);
cel=cell2mat(vertexAttachments(sc.DelaunayT,mnID)); %ids of tetrahedron touching grain i
cel=cel(ismember(cel,sc.GoodCells));
cc=[sc.sCells(cel).Center]'; %get property
cc=cc(cc(:,1)==0 | cc(:,1)==mnID | cc(:,2)==mnID, 3:end); %check correct positions not the ids
[b,vol]=boundary(cc) %calculate volume

figure
surf(X*r+gr.Coord(mnID,1),Y*r+gr.Coord(mnID,2),Z*r+gr.Coord(mnID,3))
hold on
axis equal

trisurf(b,cc(:,1),cc(:,2),cc(:,3),'Facecolor','red','FaceAlpha',0.1)

%second grain
r=sc.Radius(mxID);
cel=cell2mat(vertexAttachments(sc.DelaunayT,mxID)); %ids of tetrahedron touching grain i
cel=cel(ismember(cel,sc.GoodCells));
cc=[sc.sCells(cel).Center]'; %get property
cc=cc(cc(:,1)==0 | cc(:,1)==mxID | cc(:,2)==mxID, 3:end); %check correct positions not the ids
[b,vol]=boundary(cc) %calculate volume
figure
surf(X*r+gr.Coord(mxID,1),Y*r+gr.Coord(mxID,2),Z*r+gr.Coord(mxID,3))
hold on
axis equal
trisurf(b,cc(:,1),cc(:,2),cc(:,3),'Facecolor','red','FaceAlpha',0.1)
%}