function drawIndvLoops(step,lpIntv,app)
%DRAWINDVLOOPS draw some loops of a certain step
% This function will load a spaceCellSystem object saved on a 3D loops
% calculation. Then it will get loops size distribution and plot some
% exemples of loops present on the given step.
% - step => timestep of the spaceCell to be loaded
% - lpIntv =>  interval of loops' orders to be saved
% - app => AnalyseExe app
%Ask for file
fnm="LoopsSpaceCellsfile"+step+".mat";
if nargin<3
    %If app is not called
    sc=load("../LoopsExemple/"+fnm).sc;
    gr=grains('STENSOR',step,'','');
    path="../LoopsExemple/Loop"+step;
else
    %Load file
    sc=load(app.SavePath+"/MatlabResultsFile/SpaceCellFiles/"+fnm).sc;
    gr=grains('STENSOR',step,'',app);
    path=app.SavePath+"/MatlabResultsFile/IndvLoop"+step;
end

%force analises
%p=gr.PGStressTensor;
%p=permute((p(1,1,:)+p(2,2,:)+p(3,3,:))/3,[3 1 2]);

%get loop sizes present into the file and choose which to plot
u=unique([sc.Loops.Order]);
chk=ismember(lpIntv,u);
if sum(chk)~=numel(lpIntv)
    f=find(~chk);
    for i=1:numel(f)
        [~,mn]=min(abs(u-lpIntv(f(i))));
        lpIntv(f)=u(mn);
    end
end

%u=u(u>=lpIntv); %remove smaller than minlp
%vl=1:interval:size(u,2);
%if numel(vl)~=1 && vl(end)~=size(u,2)
%    vl=[vl,size(u,2)];
%end
%Prepare file
[~,vl]=ismember(lpIntv,u);
if exist(path,'dir')==0
    mkdir(path)
end
%start saving loops
% nmMods='A':'Z';
% sz=size(nmMods,2);
for i=1:numel(vl)
    fprintf("Saving Loop "+u(vl(i))+" \n")
    chk=find([sc.Loops.Order]==u(vl(i)),1);
    lp=sc.Loops(chk);
    cord=gr.Coord(lp.Grains,:);
    fnm=fullfile(path,""+u(vl(i))+"Loop.vtk");
    vtkwrite(fnm,'polydata',"TRIANGLE",...
        gr.Coord(:,1),gr.Coord(:,2),gr.Coord(:,3),cat(1,lp.Vertices));
    fnm=fullfile(path,""+u(vl(i))+"Grains.vtk");
    vtkwrite(fnm,'unstructured_grid',...
        cord(:,1),cord(:,2),cord(:,3),...
        'SCALARS','Radius',gr.Radius(lp.Grains),...
        'PRECISION',10,'BINARY');
    
%     wdth=max(cord(:,3))-min(cord(:,3));
%     nb=10;
%     for j=1:nb
%         md=[nmMods(floor((j-1)/sz)+1) nmMods(j-sz*floor((j-1)/sz))];
%         z=max(cord(:,3))-(j)*wdth/nb;
%         fnm=fullfile(path,""+md+"Grains.vtk");
%         vtkwrite(fnm,'unstructured_grid',...
%             cord(cord(:,3)<=z,1),cord(cord(:,3)<=z,2),cord(cord(:,3)<=z,3),...
%             'SCALARS','Radius',gr.Radius(lp.Grains(cord(:,3)<=z)),...
%             'PRECISION',10,'BINARY');
%     end
end
%total grains plot
fnm=fullfile(path,"Grains_tot.vtk");
vtkwrite(fnm,'unstructured_grid',...
    gr.Coord(:,1),gr.Coord(:,2),gr.Coord(:,3),...
    'SCALARS','Radius',gr.Radius,...
    'PRECISION',10,'BINARY');
end