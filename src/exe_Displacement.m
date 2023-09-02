function exe_Displacement(app)
%MEANDISPLACEMENTTOVTK Read Displacement LIGGGHTS Txt files and export to Vtk
%   This function will, with the help of the class 'grainsdpl', read text
%   files, calculate the mean displacement of grains, substract from the
%   actual displacement and then write the data on VTK files for better
%   visualization.

%Load values
[N1,N2,interval,~,nbFiles] = createStepArray(app);
path=MakePath(app,'DPL');

if isequal(app.DplStep0Switch.Value,'On')
    step0=app.DplStep0EF.Value;
    if step0>N1
        warndlg("Default step0 is larger than initial"+...
            "calculation step. May lead to incorrect values.");return
    end
	gr0=grains('BASIC',step0,'',app);
end

%get value s
s=str2double(regexp(app.NeighborValuesEditField.Value,'\d*','match'));
%Turn on calculation pannel
app=CalcPanel(app,'',nbFiles,'Starting calculation','on');
for i=1:nbFiles
    if getappdata(app.CalculatingPanel,'Canceling')
        CalcPanel(app,'','','','off');
        warndlg('Calculation was canceled');return
    end
    step=min(N1+interval*(i-1),N2);
    app=CalcPanel(app,i,nbFiles,nbFiles);
    gr=grains('BASIC',step,'',app);
    if isempty(gr.Displacement);CalcPanel(app,'','','off');
        warndlg("Displacement is empty on step " +step);return
    end
    %take step0 into account
    if nargin>1;gr.Displacement=gr.Displacement-gr0.Displacement;end
    %check nb of grains and s defined values
    if sum(s>gr.Nb)>0;s=[s(s<gr.Nb) gr.Nb];end
    %calculate displacement in relation to the mean value
    dpl=gr.Displacement-mean(gr.Displacement,1);
    %calculate sLID using it's function
    slid = sLID(gr,s);
    sLIDtitle=sprintf(repmat('s-LID %d|',1,numel(s)),s);
    %vtk plot
    vtkwrite(path+"meandispl"+step+".vtk",'unstructured_grid',...
        gr.Coord(:,1), gr.Coord(:,2),gr.Coord(:,3),...
        'MATRIX',sLIDtitle,slid,... 
        'VECTORS','MeanDisplacements',dpl(:,1),dpl(:,2),dpl(:,3),...
        'PRECISION',5);
end
CalcPanel(app,i+1,nbFiles,'','off');
end
function slid = sLID(gr,s)
%SLID definition as ZHOU, TORDESILLAS 2021
% This function calculate the Local intrinsic dimensionality (LID) of each
% grain in the framework of the comulative displacement. This value is
% calculated for 's' closest neighbors of the grain.

slid=zeros(gr.Nb,numel(s));
for i=1:gr.Nb
    %find the closest neighbors (sort by smaller distance)
    dist=((gr.Coord(:,1)-gr.Coord(i,1)).^2+...
        (gr.Coord(:,2)-gr.Coord(i,3)).^2+...
        (gr.Coord(:,2)-gr.Coord(i,3)).^2).^(1/2);
    [~,I]=sort(dist);
    %get the displacement distance. Using I to put in the right order in
    %relation to the distance
    u=((gr.Displacement(I,1)-gr.Displacement(I,1)).^2+...
        (gr.Displacement(I,2)-gr.Displacement(I,3)).^2+...
        (gr.Displacement(I,2)-gr.Displacement(I,3)).^2).^(1/2);
    %for each s value sum the 's' closest neighbors
    for j=1:numel(s)
        slid(i,j)=-(sum(u(1:s(j))/u(s(j)))/s(j))^(-1);
    end    
end
end