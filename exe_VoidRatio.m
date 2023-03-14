function exe_VoidRatio(PD,app)
%VOIDRATIO Calculate the void ratio of the image
%   This function will calculate the void ratio for the test in many
%   diferent ways, depending on the input parameters.
%   mode will tell if we are calculation total or partial void ratio
%
%   INPUTS %%%%%
%   PartData will contain the partial data in the case of a partial
%   calculation, else it will be empty
%
%   CALCULATIONS %%%%%
%   Total calculation is made by calculating the total volume inside the
%   trial and comparing with the total volume of the grains.Easy and fast.
%
%   Partial 3D calculation is made by calculating the volume of the cut
%   rectangle and the volume of all particles that are inside it.
%
%   Partial 2D calculation is made by creating a image and calculating the
%   percentages of pixels of each color. First we create a plot containing
%   the area of the experiment in 'black' and the area of grains in 'grey'.
%   Then we will cut out the part that interest us and count how many
%   pixels are of each color. Givin us the voidratio.

%check Varargin files
if isequal(PD,'LOAD');voidLoad(app);return;end

%get needed variables
N1=app.N1EF.Value;
N2=app.N2EF.Value;
interval=app.CalcInt.Value;
if interval==app.IntervalEF.Value && app.SimType==3
    qst=app.QcstStep(2);
    if qst>N2
        stepArray=(N1:interval:N2)';
    elseif qst<=N1
        f=find(app.TrialData.Step==N1 | app.TrialData.Step==N2);
        stepArray=app.TrialData.Step(f(1):f(2));
    else
        stepArray=(N1:interval:qst)';
        f=find(app.TrialData.Step==qst | app.TrialData.Step==N2);
        stepArray=[stepArray(1:end-1);app.TrialData.Step(f(1):f(2))];
    end
else
    stepArray=(N1:interval:N2)';
    if stepArray(end)~=N2; stepArray=[stepArray;N2];end
end
nbFiles=numel(stepArray);
%file containing piston displacements
TD=app.TrialData;
if isempty(TD.Dz);return;end

%Important values
if isempty(PD);type="total";
else;type="partial";
end

%Check subdivision and do necessary calculations
%Start Vector containing wanted values
if app.SubdivisionButton.Value==1
    %Load important subdivision values
    type="subd";
    %read partial data variables
    lMax=PD.SubLin;	
    cMax=PD.SubCol;	    
    fnm=type+"VoidRatio"+N1+"to"+N2+"int"+interval+...
        "sub"+lMax+"x"+cMax+".mat";
    Results=zeros(nbFiles,2,lMax*cMax+1);
else
    lMax=1;cMax=1;
    fnm=type+"VoidRatio"+N1+"to"+N2+"int"+interval+".mat";
    Results=zeros(nbFiles,2,2);
end

app=CalcPanel(app,'',nbFiles,'Starting calculation','on');
%Begin Calculation

for i=1:nbFiles
    if getappdata(app.CalculatingPanel,'Canceling')
        CalcPanel(app,'','','','off');
        warndlg('Calculation was canceled');return
    end
    app=CalcPanel(app,i,nbFiles,stepArray(i));
    switch upper(type)
        case 'TOTAL'
            %Calculate grain and check if the files were well chosen
            gr=grains('BASIC',stepArray(i),'',app);
            if isempty(gr.Radius);CalcPanel(app,'','','','off');return;end
            
            %Calculate volumes - for total calculations (PartData is empty)
            %the function voroVolume returns scalar values
            [vS,vT]=voidVolume(type,gr,stepArray(i),app,'');
            vV=vT-vS;
            
            %Save Valuesmode
            Results(i,:,1)=[vV/vS vV/vT];
            
            %Second version of calculation, only taking into account
            %"goodCells" of spaceCellSystem object
            sc=spaceCellSystem("VR",stepArray(i),gr,app,'');
            [~,vs,vv]=totalVoidRatio(sc,gr);
            Results(i,:,2)=[vv/vs vv/(vs+vv)];
            
        case 'PARTIAL'
            %Partial 3D calcultions will be made in a similar way to the
            %total one. The volume of grains inide or partially inside the
            %rectangle will be taken into account and deducted from the
            %total. It will be done with the help of the "solidVolume"
            %function found in the later parts of this file.
            
            %Calculate grain and check if the files were well chosen
            [gr,PD]=grains('BASIC',stepArray(i),PD,app);
            if isempty(gr.Radius);CalcPanel(app,'','','','off');return;end
            
            %Calculate the volumes for partial calculations. The function
            %voroVolume will calculate the volume for the grains inside the
            %rectangle and add 0 to all the others.
            [vS,vT]=voidVolume(type,gr,stepArray(i),app,PD);
            vV=vT-vS;
            
            %Void calculation - save Step VoidRatio and Porosity
            Results(i,:,1)=[vV/vS vV/vT];
            
        case 'SUBD'
            %Partial 3D calcultions will be made in a similar way to the
            %total one. The volume of grains inide or partially inside the
            %rectangle will be taken into account and deducted from the
            %total. It will be done with the help of the "solidVolume"
            %function found in the later parts of this file.
            
            %calculation
            
            %Calculate grain and check if the files were well chosen
            [gr,PD]=grains('BASIC',stepArray(i),PD,app);
            if isempty(gr.Radius);CalcPanel(app,'','','','off');return;end
            
            %Calculate the volumes for partial calculations. The function
            %voroVolume will calculate the volume for the grains inside the
            %rectangle and add 0 to all the others.
            [vS,vT]=voidVolume('PARTIAL',gr,stepArray(i),app,PD);
            vV=vT-vS;
            
            %Void calculation - save Step VoidRatio and Porosity
            Results(i,:,1)=[vV/vS vV/vT];
            
            %same calculation for each subdivision if needed
            [vSsub,vTsub]=voidVolume(type,gr,stepArray(i),app,PD);
            vVsub=vTsub-vSsub;
            for j=1:lMax*cMax
                %Now we have all volumes we need to get the volume only of the grains that
                %were selected with PartData subdivision
                AG=PD.SubGrains{j};  %grains inside the subdv
                if isempty(AG) || vTsub(j)==0
                    Results(i,:,j+1)=[0 0 ];continue
                end
                
                %Void calculation - save Step VoidRatio and Porosity
                Results(i,:,j+1)=[vVsub(j)/vSsub(j) vVsub(j)/vTsub(j) ];
            end
    end
end
CalcPanel(app,i+1,nbFiles,'','off');

%Calculate the deformation
snExt = extStrains(app.TrialData,stepArray,N1,app,'dev');
ssExt = extStress(app.TrialData,stepArray,app);
consoStrain = extStrains(app.TrialData,app.ConsoStep,N1,app);

res.VoidRatio=Results;
res.Strain=snExt(:,end-2);
res.Pressure=ssExt(:,end);

pD=plotData("Normal",res,app,type,consoStrain(end-1));

fnm=fullfile(MakePath(app,'VOID'),fnm);
save(fnm,'pD','-v7.3');  

%Write down files
voidPlotter(app,'',pD);
end
function voidLoad(app)
%Load previously calculated file

[fnm,fPath]=MatLoader('VOID',app);
if fPath==0;return;end
if iscell(fnm)
    pD=plotData.empty(0,1);
    for i=1:numel(fnm)
        try pD(i)=load(fullfile(fPath,fnm{i})).pD;
        catch
            warndlg(['File ' fnm{i} ' not found or does not cointain '...
                '"pD" property, Try again.']);return
        end
        if isempty(pD(i).FileName)
            nm=fnm{i};
            pD(i).FileName=nm(1:find(nm=='.')-1);
        end
    end
else
    pD=load(fullfile(fPath,fnm)).pD;
end
%}
fprintf('File loaded. Starting calculations \n')


if isempty(pD(1).SubL)
    mode='NORMAL';
else
    mode='SUBD';
end

voidPlotter(app,mode,pD);
end
function voidPlotter(app,mode,pD)
%funciton to plot the values
if isequal(upper(mode),'SUBD')
    Nbplot=size(pD.Results,3);
    subV=[pD.SubL,pD.SubC];
    strg=strings(1,Nbplot);
    strg(1)="Full Selection";
    lin=1;col=1;
    for i=2:Nbplot
        if ((i-1)-lin*subV(2))>0
            lin=lin+1;
            col=1;
        end
        strg(i)="Subdivision l"+lin+" c"+col;
        col=col+1;
        
    end
    SubPlotter(app, 'VOID', pD.Results, strg,pD.ConsoStrain, subV);
    return;
end
if numel(pD)>1
    path=MakePath(app,'VOIDL');
else
    path=MakePath(app,'VOID');
end
png=".png";
%fig=".fig";
%Check title and legends option
if app.TitlesCB.Value;tit=1;else;tit=0;end
if app.LegendsCB.Value;leg=1;else;leg=0;end

%Create Figures and Axis
nb=3;
if size(pD,2)>1;nb=2*nb;end

f(nb)=figure;ax(nb)=axes(f(nb));hold(ax(nb),'on');
for i=1:(nb-1)
    f(i)=figure;ax(i)=axes(f(i));hold(ax(i),'on'); %#ok<LAXES>
end

C=app.PlotColors;
crPt=zeros(numel(pD),2); %create a line passing through last pts
for i=1:size(pD,2)
    res=pD(i).Results;
    if size(pD,2)>1
        optsA={'Color',C(i,:)};optsB={'Color',C(i,:)};
    else
        optsA={};optsB={};
    end
    for j=1:numel(pD(i).InfPts.q)
        optsA=[optsA,{'Pointx',pD(i).InfPts.ez(j)}]; %#ok<*AGROW>
        optsB=[optsB,{'Pointx',log10(pD(i).InfPts.p(j))}];
    end
    %Turn values in seconds
    plotMark(app,ax(1),res.Strain,res.VoidRatio(:,1,1),optsA{:}) % voidratio = f(Ez)
    plotMark(app,ax(2),res.Strain,res.VoidRatio(:,2,1),optsA{:}) % porosity = f(Ez)
    plotMark(app,ax(3),log10(res.Pressure),res.VoidRatio(:,1,1),optsB{:})   % voidratio = f(log(p))
    crPt(i,:)=[log10(res.Pressure(end)),res.VoidRatio(end,1,1)];
    if size(pD,2)>1;k=1;else;k=0;end
    plotMark(app,ax(3*k+1),res.Strain,res.VoidRatio(:,1,2),optsA{:}) % voidratio = f(Ez)
    plotMark(app,ax(3*k+2),res.Strain,res.VoidRatio(:,2,2),optsA{:}) % porosity = f(Ez)
    plotMark(app,ax(3*k+3),log10(res.Pressure),res.VoidRatio(:,1,2),optsB{:})   % voidratio = f(log(p))
    
end
if leg
    %If legend, do legemd following the number of files loaded
    if size(pD,2)>1
        for j=1:nb
            legend(ax(j),pD.FileName,'location','best')
        end
    else
        for j=1:nb
            legend(ax(j),'Whole specimen','Good Cells','location','best')
        end
    end
end
if size(pD,2)>2
    %if more than two files, create a fitted line through the end points -
    %tring to simulate the critical state
    ft=fit(crPt(:,1),crPt(:,2),'poly1');
    axes(ax(3))
    p=plot(ft,'--k',crPt(:,1),crPt(:,2));
    p(2).LineWidth=1;
    delete(ax(3).Children(2))
end

if size(pD,2)>1
    %if more then one file, make sure the two kinds of void ratio are in
    %the same scale
    lm=zeros(3,2);
    for i=1:3
        lm(i,:)=[min(ax(i).YLim(1),ax(i+3).YLim(1)),...
            max(ax(i).YLim(2),ax(i+3).YLim(2))] ;
    end
end

%plot1 Void Ratio = f(Ez)
i=1;
if tit;title(ax(i),'Evolution of Void Ratio');end
ylabel(ax(i),'Void raito')
xlabel(ax(i),'Strain z axis')
if size(pD,2)>1;ax(i).YLim=lm(i,:);end
% ax(1).YLim=[0.4,0.8];
fnm="VoidRatio";
saveas(f(i),fullfile(path,fnm+png));
i=i+1;

%plot2 Porosity = f(Ez)
if tit;title(ax(i),'Evolution of Porosity');end
ylabel(ax(i),'Porosity')
xlabel(ax(i),'Strain z axis')
if size(pD,2)>1;ax(i).YLim=lm(i,:);end
fnm="Porosity";
saveas(f(i),fullfile(path,fnm+png));
i=i+1;

%plot3 VR=f(log(p))
if tit;title(ax(i),'Void ratio in function of log(p)');end
ylabel(ax(i),'Void Ratio')
xlabel(ax(i),'log(p)')
if size(pD,2)>1;ax(i).YLim=lm(i,:);end
fnm="VoidRatio_LogP";
saveas(f(i),fullfile(path,fnm+png));
i=i+1;

if numel(pD)>1
    if tit;title(ax(i),'Evolution of Void Ratio');end
    ylabel(ax(i),'Void raito')
    xlabel(ax(i),'Strain z axis')
    if size(pD,2)>1;ax(i).YLim=lm(i-3,:);end
    % ax(1).YLim=[0.4,0.8];
    fnm="VoidRatio_GC";
    saveas(f(i),fullfile(path,fnm+png));
    i=i+1;
    
    %plot2 Porosity = f(Ez)
    if tit;title(ax(i),'Evolution of Porosity');end
    ylabel(ax(i),'Porosity')
    xlabel(ax(i),'Strain z axis')
    if size(pD,2)>1;ax(i).YLim=lm(i-3,:);end
    fnm="Porosity_GC";
    saveas(f(i),fullfile(path,fnm+png));
    i=i+1;
    
    %plot3 VR=f(log(p))
    if tit;title(ax(i),'Void ratio in function of log(p)');end
    ylabel(ax(i),'Void Ratio')
    xlabel(ax(i),'log(p)')
    if size(pD,2)>1;ax(i).YLim=lm(i-3,:);end
    fnm="VoidRatio_GC_LogP";
    saveas(f(i),fullfile(path,fnm+png));
    
end

%If legend is turned off, create a legend file that may be used as an
%outside legend
if ~leg && numel(pD)>1
    %Multi legend
    o=copyobj(f(1),0);
    l=legend(o.CurrentAxes,pD.FileName,'Orientation','horizontal');
    l.EdgeColor='none';
    set(o.CurrentAxes,'Visible','Off')
    % Set the figure Position using the normalized legend Position vector
    % as a multiplier to the figure's current position in pixels This sets
    % the figure to have the same size as the legend
    set(o,'Position',(get(l,'Position').*[0, 0, 1, 1].*get(o,'Position')));
    % The legend is still offset so set its normalized position vector to
    % fill the figure
    set(l,'Position',[0,0,1,1]);
    % Put the figure back in the middle screen area
    set(o, 'Position', get(o,'Position') + [500, 400, 0, 0]);
    saveas(o,fullfile(path,"Legend"+png));
    
    %Delete extra figures
    delete(o);
end

%Delete figures in the case of exe all
if isequal(app.LIGGGHTSAnalysisButtonGroup.SelectedObject,app.ExeAllButton)
    delete(f);
end
end

%Support functions
function [vS,vT]=voidVolume(mode,gr,step,app,PD)
%VOROVOLUME Calculates the total and the solid volume of the specimen
%   This function will create a delaunayTriangulation object using the
%   grains position and the 'fakegrains' placed in the walls position. Then
%   the voronoi of this triangulation will be calculated. Two diferent
%   calculations may now take place:
%   - for partial calculation each grain will have the its voronoi volume
%   calculated, the values returned are VECTORS.
%   - for total calculations the total volume will be accounted, the values
%   returned are SCALARS

%Geometric data
R=gr.Radius;

%Voronoi calculation
[vT,vTsub] = getVolume(app.TrialData,app,step,PD);

%Start volume calculations
switch upper(mode)
    case 'TOTAL'
        %Solid Volume
        if app.Bool3D
            vS=4/3*pi()*dot(R.*R,R);
        else
            vS=pi()*dot(R,R);
        end
    case 'PARTIAL'
        %Grains totally inside
        Ri=R(PD.GrainsInt);
        if app.Bool3D
            vS=4/3*pi()*dot(Ri.*Ri,Ri);
        else
            vS=pi()*dot(Ri,Ri);
        end
        %Grains partially inside
        idPi=setxor(PD.AllGrains,PD.GrainsInt);
        [vPi]=volumeInside(app,gr,idPi,PD.Vertices);
        
        %total solid volume
        vS=vS+sum(vPi);
    case 'SUBD'
        %cell volume calculation - for each voronoi cell we get its vertices IDs
        %and their position. With these values and the function 'boundary' the
        %volume is calculated.
        vS=zeros(size(PD.SubGrainsI,1),1); %will contain the volume of each grain
        for i=1:size(PD.SubGrainsI,1)
            %Grains totally inside
            Ri=R(PD.SubGrainsI{i});
            if app.Bool3D
                vSi=4/3*pi()*dot(Ri.*Ri,Ri);
            else
                vSi=pi()*dot(Ri,Ri);
            end
            %Grains partially inside
            idPi=setxor(PD.SubGrainsA{i},PD.SubGrainsI{i});
            vPi=volumeInside(app,gr,idPi,PD.SubVertices(:,:,i));
            vS(i)=vSi+sum(vPi);
        end
        vT=vTsub;
end

end
function vPi=volumeInside(app,gr,idPi,vert)
%VOLUMEINSIDE calculates the volume of each grain inside the rectangle
% This function will calculate points for the surface of each grain idPi
% that is partially inside the rectangle. These points will then be checked
% with the rectangle boundaries to find which ones are inside. The
% convexhull will then be calculated to get the volume.

%Check for dimension, prepare which colum to use in calculations
if app.Bool3D
    D=3;p=[2,3];
else
    D=2;p=[1,2];
end

dT=10;                      %delta angle of point creation
bs=baseSphere(dT,D);        %create base sphere
vPi=zeros(size(idPi,1),1);  %volume storing vector
for i=1:size(idPi,1)
    r=gr.Radius(idPi(i));   %radius
    cd=gr.Coord(idPi(i),:); %origin
    %multiply by radius and add origin
    grC=bs*r+cd;    
    %get grains inside
    ins=inpolygon(grC(:,p(1)),grC(:,p(2)),vert(:,1),vert(:,2));
    %get volume
    try [~,vPi(i)] = convhulln(grC(ins,:)) ;
    catch 
        %If there is not enough points inside to do the calculation,
        %consider vPi=0
        vPi(i)=0;
    end
end
end
