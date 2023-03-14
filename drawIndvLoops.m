function drawIndvLoops(step,interval,minlp,app)
%DRAWINDVLOOPS draw some loops of a certain step
% This function will load a spaceCellSystem object saved on a 3D loops
% calculation. Then it will get loops size distribution and plot some
% exemples of loops present on the given step.

%Ask for file
fnm="LoopsSpaceCellsfile"+step+".mat";
if nargin<4
    sc=load("../LoopsExemple/"+fnm).sc;
    gr=grains('STENSOR',step,'','');
    path="../LoopsExemple/Loop"+step;
else
    %Load file
    sc=load(app.SavePath+"/MatlabResultsFile/SpaceCellFiles/"+fnm).sc;
    gr=grains('STENSOR',step,'',app);
    path=app.SavePath+"/MatlabResultsFile/IndvLoop"+step;
end
if ~exist('interval','var')
    interval=0;
elseif isempty(interval)
    interval=0;
end
if ~exist('minlp','var')
    minlp=0;
elseif isempty(minlp)
    minlp=0;
end

%force analises
p=gr.PGStressTensor;
p=permute((p(1,1,:)+p(2,2,:)+p(3,3,:))/3,[3 1 2]);
%get loop sizes present into the file and choose which to plot
u=unique([sc.Loops.Order]);
u=u(u>=minlp); %remove smaller than minlp
vl=1:interval:size(u,2);
if numel(vl)~=1
    if vl(end)~=size(u,2);vl=[vl,size(u,2)];end
end
%Prepare file
if exist(path,'dir')==0
    mkdir(path)
end
%start saving loops
nmMods='A':'Z';
sz=size(nmMods,2);
for i=1:size(vl,2)
    fprintf("Saving Loop "+u(vl(i))+" \n")
    f=find([sc.Loops.Order]==u(vl(i)),1);
    lp=sc.Loops(f);
    cord=gr.Coord(lp.Grains,:);
    md=[nmMods(floor((i-1)/sz)+1) nmMods(i-sz*floor((i-1)/sz))];
    fnm=""+md+u(vl(i))+"Loop.vtk";
    vtkwrite(fnm,'polydata',path,"TRIANGLE",...
        gr.Coord(:,1),gr.Coord(:,2),gr.Coord(:,3),cat(1,lp.Vertices));
    fnm=""+md+u(vl(i))+"Grains.vtk";
    vtkwrite(fnm,'unstructured_grid',path,...
        cord(:,1),cord(:,2),cord(:,3),...
        'SCALARS','Radius',gr.Radius(lp.Grains),...
        'SCALARS','StressState',p(lp.Grains),...
        'PRECISION',10,'BINARY');
    
    wdth=max(cord(:,3))-min(cord(:,3));
    nb=10;
    for j=1:nb
        md=[nmMods(floor((j-1)/sz)+1) nmMods(j-sz*floor((j-1)/sz))];
        z=max(cord(:,3))-(j)*wdth/nb;
        fnm=""+md+"Grains.vtk";
        vtkwrite(fnm,'unstructured_grid',path,...
        cord(cord(:,3)<=z,1),cord(cord(:,3)<=z,2),cord(cord(:,3)<=z,3),...
        'SCALARS','Radius',gr.Radius(lp.Grains(cord(:,3)<=z)),...
        'PRECISION',10,'BINARY');
    end
end
%total grains plot
fnm="Grains_tot.vtk";
vtkwrite(fnm,'unstructured_grid',path,...
    gr.Coord(:,1),gr.Coord(:,2),gr.Coord(:,3),...
    'SCALARS','Radius',gr.Radius,...
    'SCALARS','StressState',p,...
    'PRECISION',10,'BINARY');
end