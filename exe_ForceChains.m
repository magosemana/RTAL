function exe_ForceChains(app,varargin)
%FORCECHAINS Will read LIGGGHTS files and export data on the Loops
%   This function will use the class 'grains' to read files and then
%   calculate the per grain internal stress tensor. We will then select
%   those that have values above the mean and similar directions to
%   determinate the force chains. Based on
%   Peters,Muthswamy,Wibowo,Tordesillas 2005

%% Preparation
%Check startup options
if nargin>1;FCLoad(app);return;end

%Load values
[N1,N2,interval,stepArray,nbFiles] = createStepArray(app);
fcPath=MakePath(app,'FC');

%Initiate variables for each type of calculation asked\
%base Calculation
if app.FCBaseCBox.Value
    %matrix to contain the data : Step - mean(Length) - max (Length) -
    %numel(fc) - numel(uniquefc) - numel(sigBranch) - nfc/HSTG  - num 3p
    bData=zeros(nbFiles,9);
    angData=cell(nbFiles,1);
    fcGrchck=zeros(app.NbGrains,1);
    %Prepare directories for VTK files
    pathVTKfc=fcPath+("ForceChains"+N1+"to"+N2+"int"+interval+"/");
    if exist(pathVTKfc,'dir')==0;mkdir(pathVTKfc);end
end
%Cluster
if app.FCClusterCBox.Value
    %Cell array containing Nx4 files on each cell. N representing each fc
    %calculated at each step and 4 being the Cluster size-mean/max/min
    %cluster order
    clData=cell(nbFiles,4);
    pathVTKcl(1)=fcPath+("ForceChainsClst"+N1+"to"+N2+"int"+interval+"/");
    if exist(pathVTKcl(1),'dir')==0;mkdir(pathVTKcl(1));end
end
%Grain VTK evolution
if app.FCGrVtkCBox.Value
    pathVTKcl(2)=fcPath+("ForceChainsGr"+N1+"to"+N2+"int"+interval+"/");
    if exist(pathVTKcl(2),'dir')==0;mkdir(pathVTKcl(2));end
end
%Bending and 
fcEleO=''; angO='';
bendRes=zeros(nbFiles,1);

%Bending+ClusterTransformation
%If cluster transformation also asked, first try loading file that is
%found in usual path. Else demand the user to choose the path
if app.FCBendCLTCBox.Value
    fcClTf=1;
    d=dir(MakePath(app,'LOOPCT','chk'));
    k='Loops-Transf';
    f=strncmp({d.name},k,numel(k));
    fnm=[MakePath(app,'LOOPCT','chk') d(f).name];
    try ctRes=load(fnm).ctRes;
    catch
        if isequal(app.LIGGGHTSAnalysisButtonGroup.SelectedObject,app.ExeAllButton)
            fprintf(['Cluster transformation file not found in usual...' ...
                'path, it will NOT be executed.\n']);fcClTf=0;
        else
            fprintf(['Cluster transformation file not found in usual...' ...
                'path, please identify the file locaiton.\n'])
            [fnm,fPath]=MatLoader('LOOPCT',app,'off');
            if fPath==0;warndlg('No files were chosen, Try again.');return;end
            ctRes=load([fnm,fPath]).ctRes;
        end
    end
    if fcClTf && numel(ctRes.Strain)~=(nbFiles)
        warndlg(['Loaded transformation file do not have the same '...
            'number of steps than this function. Check the chosen' ...
            ' file or calculation interval']);return
    end
else
    fcClTf=0;
end

%check dimesion
if app.Bool3D; D=3;else; D=2;end

%% Execution
%Turn on calculation pannel
app=CalcPanel(app,'',nbFiles,'Starting calculation','on');

%Start the Loop
for i=1:nbFiles
    if getappdata(app.CalculatingPanel,'Canceling')
        CalcPanel(app,'','','','off');
        warndlg('Calculation was canceled');return
    end
    step=stepArray(i);
    app=CalcPanel(app,i,nbFiles,step);
    %calculate force chian
    gr=grains('FORCECHAIN',step,'',app); 
    fc=gr.ForceChains;
    if isempty(fc)
        CalcPanel(app,'','','','off');
        warndlg("Force chain is empty on step " +step);return
    end
    %if any of the clusters is on, load spaceFiles of this step
    if app.FCClusterCBox.Value || fcClTf
        %load loops data
        scFnm=[MakePath(app,'SCF','check')...
            char("LoopsSpaceCellsfile"+step+".mat")];
        try sc=load(scFnm).sc;
        catch
            fprintf(['No Loops spacecell files, '...
                'please run loops calculation beforehand.\n'])
            CalcPanel(app,'','','','off');
            return
        end
    end
    
    %Base force chain data treattement
    if app.FCBaseCBox.Value
        %Time - mean(Length) - max (Length) - numel(fc) - ...
        %numel(uniquefc) - numel(sigBranch) - nfc/HSTG - nfc/Tot
        bData(i,:)=fcBase(fc,gr,step,pathVTKfc);
        angData{i}=cat(1,fc.Elevation);
        grIds=cat(1,fc.IDs);
        fcGrchck(grIds)=fcGrchck(grIds)+1;
    end
    %Force chain clusters
    if app.FCClusterCBox.Value
        [clData{i,1:4}]=fcClusters(app,fc,gr,sc,step,pathVTKcl);
        if isempty(clData{i,4});return;end
        %{FcCldata,nb fc grains,av(cl) fc grains,nb others grains, av(cl) others}
    end
    %Force chain bending and cluster transformation
    [bEv,fcEleO,angO,bID]=fcBending(app,fc,gr,fcEleO,angO,fcClTf);
    bendRes(i)=bEv; %save the number of bending events
    %Execute cluster transformation. Check in the old step (scO) if the
    %grains involved in bending event also changed in cluster order
    if fcClTf
        if i>1
            chg=fcBendingClTf(bID,scO,ctRes.ClTransf{i,1});
            ctRes.ClTransf{i,1}=chg;
        end
        scO=sc;
    end
    
end
CalcPanel(app,i+1,nbFiles,'','off');
%save(fnm,'grs','-v7.3');
%% Save and plot
%Prepare Strain data
consoStrain = extStrains(app.TrialData,app.ConsoStep,N1,app);
stn = extStrains(app.TrialData,stepArray,N1,app,'allCalc');
sts = extStress(app.TrialData,stepArray,app);
if app.FCBaseCBox.Value
    %Finish matrix
    res.Base=[stn(:,D) sts(:,end) bData(:,2:end) bendRes]; %Ez p vals
    res.Elevation=angData;
    fnm=MakePath(app,'FC')+"ForceChains-Base"+N1+"to"+N2+".mat";
    %Create plotData object and plot
    fcGrchck=fcGrchck/nbFiles;
    pD=plotData("Normal",res,app,'',consoStrain(D));
    save(fnm,'pD','-v7.3'); 
    fcPlotter(app,pD,fcGrchck);
    
    %create vtk state variables file
    fnm=pathVTKfc+"vtkVariables.txt";
    StepEzEvQP=[stepArray,stn(:,[D,D+1]),sts(:,end-1:end)];
    vtkLog(app,'FCbase',fnm,StepEzEvQP)
end
if app.FCClusterCBox.Value
    %Prepare object for plot
    res.Data=clData;
    res.Strain=stn(:,end-2); %strain Z
    res.Stress=sts(:,end-2); %stressZ
    res.Pressure=sts(:,end); %stress p
    fnm=MakePath(app,'FCCL')+"ForceChains-Cluster"+N1+"to"+N2+".mat";
    %Create plotData object and plot
    pD=plotData("Normal",res,app,'');
    save(fnm,'pD','-v7.3');%save
    fcClstPlotter(app,pD);
    
    %create vtk state variables file
    fnm=pathVTKcl(1)+"vtkVariables.txt";
    StepEzEvQP=[stepArray,stn(:,[D,D+1]),sts(:,end-1:end)];
    vtkLog(app,'FCbase',fnm,StepEzEvQP)
end
%ClusterTransformation
if fcClTf==1
    ctRes.NbClst=''; %remove non used values
    ctRes.NbCell='';
    pD=plotData("Normal",ctRes,app,''); 
    fnm=MakePath(app,'FCCLTF')+"ForceChains-CltTrf"+N1+"to"+N2+".mat";
    save(fnm,'ctRes','-v7.3'); %save
    fcTransfPlotter(plotData("FastCreation",pD,''),app)
end

end
%load functions
function FCLoad(app)
%FCLOAD load file for Force Chain functions
tf=1;
if app.FCBaseCBox.Value
    %Load file
    [fnm,fPath]=MatLoader('FC',app);
    if fPath==0;warndlg('No files were chosen, Try again.');return;end
    if iscell(fnm)
        nb=numel(fnm);
        pD(nb)=load(fullfile(fPath,fnm{end})).pD;
        for i=1:nb
            if i~=nb;pD(i)=load(fullfile(fPath,fnm{i})).pD;end
            fn=fnm{i};
            pD(i).FileName=fn(1:(find(ismember(fn,'.'),1,'last')-1));
        end
    else
        pD=load(fullfile(fPath,fnm)).pD;
    end
    if isempty(pD);return;end
    %Plot
    fcPlotter(app,pD);
    tf=0;
end
if app.FCClusterCBox.Value
    oldF=cd;
    [fnm,fPath]=MatLoader('FCCL',app);
    cd(oldF)
    if fPath==0;warndlg('No files were chosen, Try again.');return;end
    if iscell(fnm)
    	nb=numel(fnm);
        fprintf('Preparing images\n')
        for i=1:nb
            pD=load(fullfile(fPath,fnm{i})).pD;
            fn=fnm{i};
            pD.FileName=fn(1:(find(ismember(fn,'.'),1,'last')-1));
            fcClstPlotter(app,pD,pD.FileName);
            fprintf("File "+fnm{i}+ " done\n")
        end
        fprintf('Done\n')
    else
        pD=load(fullfile(fPath,fnm)).pD;
        fcClstPlotter(app,pD);
    end
    tf=0;
end
if tf
    warndlg('Bending function cannot be loaded or no other box marked.')
end
end
%plot functions
function fcPlotter(app,pD,varargin)
%FCPLOTTER Plot data for Force Chain calculations


% N1=app.N1EF.Value;
% N2=app.N2EF.Value;
% TD  = inflectionPoints(app.TrialData,app);
% pD.InfPts=TD.InfPts;
% fnm=MakePath(app,'FC')+"ForceChains-Base"+N1+"to"+N2+".mat";
% save(fnm,'pD','-v7.3');

if numel(pD)>1
    path=MakePath(app,'FCL');
else
    path=MakePath(app,'FC');
end
set(0,'defaultAxesFontSize',app.FontSizeEF.Value)
png=".png";
xlab='Axial strain';
x='e';
bend=1;
for i=1:numel(pD)
    if pD(i).SimType==3
       xlab='Mean pressure (kPA)';
       x='p';  
    end
    if size(pD(i).Results.Base,2)<11
        bend=0;
    end
end
%fig=".fig";
if numel(pD)<8
    C = app.PlotColors;
else
    C = graphClrCode(size(pD,2));%plot colorcode
end

%Check title and legends option
if app.TitlesCB.Value;tit=1;else;tit=0;end
if app.LegendsCB.Value;leg=1;else;leg=0;end

%Create Figures and Axis
nb=10+bend;
if numel(pD)==1;nb=nb+1;end
f(nb)=figure;ax(nb)=axes(f(nb));hold(ax(nb),'on');
for i=1:(nb-1)
    f(i)=figure;ax(i)=axes(f(i));hold(ax(i),'on'); %#ok<LAXES>
end

%Strain(1) - Stress(2) - Lngth MB(3) - Lngth SB(4) - nbFC(5) - ...
%nb branches(6) - nb SB(7) - nfc/HSTG(8) - nfc/Tot(9) - 3p(10) - Bending(11)
for i=1:size(pD,2)
    j=1;res=pD(i).Results.Base;
    v={};
    for a=1:numel(pD(i).InfPts.q)
        if strcmpi(x,'p')
            v=[v,{'Pointx',pD(i).InfPts.p(a)}]; %#ok<AGROW>
        else
            v=[v,{'Pointx',pD(i).InfPts.ez(a)}]; %#ok<AGROW>
        end
    end
    if strcmpi(x,'p');k=2;else;k=1;end %stress or strain X axis
    xAx=res(:,k);
    
    %first graph : number of force chains
    plotMark(app,ax(j),xAx,res(:,5),'Color',C(i,:),v{:});j=j+1;
    %second graph : mean multi-branch length
    plotMark(app,ax(j),xAx,res(:,3),'Color',C(i,:),v{:});j=j+1;
    %third graph : mean single-branch length
    plotMark(app,ax(j),xAx,res(:,4),'Color',C(i,:),v{:});j=j+1;
    %fourth graph : FC/HS
    plotMark(app,ax(j),xAx,res(:,8),'Color',C(i,:),v{:});j=j+1;
    %fifth graph : FC/Tot
    plotMark(app,ax(j),xAx,res(:,9),'Color',C(i,:),v{:});j=j+1;
    %sixth graph : nb signlebranch/total
    plotMark(app,ax(j),xAx,res(:,7)./...
        res(:,6),'Color',C(i,:),v{:});j=j+1;
    %seventh graph : 3p
    plotMark(app,ax(j),xAx,res(:,10),'Color',C(i,:),v{:});j=j+1;
    %eitght graph : total of branches
    plotMark(app,ax(j),xAx,res(:,6),'Color',C(i,:),v{:});j=j+1;
    %Bending
    if bend
%         if ~isempty(v) && pD(i).SimType==3
%             %fixing error in Qcst simuations where bending events should
%             %not be appearing but were because piston was advancing too
%             %fast
%             [~,f1]=min(abs(xAx-v{2}));
%             f2=find(res(:,end)<10,1,'last');
%             bl=zeros(size(res,1),1);
%             bl(f1:f2)=1;
%             fi=(res(:,end)>10 & logical(bl));
%             res(fi,end)=floor(rand(sum(fi),1)*10);
%         end
        plotMark(app,ax(j),xAx,res(:,end),'Color',C(i,:),v{:});j=j+1;
    end
    %Nb single branch and nb multible branch
    plotMark(app,ax(j),xAx,res(:,7),'Color',C(i,:),v{:});
    if numel(pD)==1
        plotMark(app,ax(j),xAx,(res(:,5)-res(:,7)),v{:});
        if leg
            legend(ax(j),'Single Branch','Multiple Branchs','Location','northeast')
        end
    else
        j=j+1;
        plotMark(app,ax(j),xAx,(res(:,5)-res(:,7)),'Color',C(i,:),v{:});
    end  
    
    
end

if leg && numel(pD)>1
    for i=1:nb
        legend(ax(i),pD.FileName,'Location','best')
    end
end

i=1;
%1 graph : number of force chains TREES
if tit;title(ax(i),'Evolution  of force chains');end
if strcmpi(x,'p');ax(i).XLim(1)=0;end
ylabel(ax(i),'Number')
xlabel(ax(i),xlab)
fnm="FC_Number";
saveas(f(i),fullfile(path,fnm+png));

i=i+1;
%2 graph : mean force chain length 
if tit;title(ax(i),'Evolution of the mean length of MB force chains');end
if strcmpi(x,'p');ax(i).XLim(1)=0;end
ylabel(ax(i),'Length of force chains')
xlabel(ax(i),xlab)
fnm="FC_Length_Mean_MB";
saveas(f(i),fullfile(path,fnm+png));

i=i+1;
%3 graph : mean force chain length 
if tit;title(ax(i),'Evolution of the mean length of SB force chains');end
if strcmpi(x,'p');ax(i).XLim(1)=0;end
ylabel(ax(i),'Length of force chains')
xlabel(ax(i),xlab)
fnm="FC_Length_Mean_SB";
saveas(f(i),fullfile(path,fnm+png));

i=i+1;
%4 graph :FC/HS
if tit;title(ax(i),'Ratio of FC/HS grains');end
if strcmpi(x,'p');ax(i).XLim(1)=0;end
ylabel(ax(i),'Ratio of chained grains')
xlabel(ax(i),xlab)
fnm="FC_Grain_HS_Ratio";
saveas(f(i),fullfile(path,fnm+png));

i=i+1;
%5 graph : FC/Tot
if tit;title(ax(i),'Ratio of FC/Tot grains');end
if strcmpi(x,'p');ax(i).XLim(1)=0;end
ylabel(ax(i),'Ratio of chained grains')
xlabel(ax(i),xlab)
fnm="FC_Grain_Tot_Ratio";
saveas(f(i),fullfile(path,fnm+png));

i=i+1;
%6 graph : Ratio of Single branch trees over total
if tit;title(ax(i),'Ratio of Single branch trees');end
if strcmpi(x,'p');ax(i).XLim(1)=0;end
ylabel(ax(i),'Ratio of force chains')
xlabel(ax(i),xlab)
fnm="FC_Branches_Ratio";
saveas(f(i),fullfile(path,fnm+png));

i=i+1;
%7 graph : Ratio of Single branch trees over total
if tit;title(ax(i),'Number of 3p');end
if strcmpi(x,'p');ax(i).XLim(1)=0;end
ylabel(ax(i),'Number of 3p')
xlabel(ax(i),xlab)
fnm="FC_3p_Number";
saveas(f(i),fullfile(path,fnm+png));
%saveas(f(nb),path+fnm+fig);

i=i+1;
%8graph : Ratio of Single branch trees over total
if tit;title(ax(i),'Total number of branches');end
if strcmpi(x,'p');ax(i).XLim(1)=0;end
ylabel(ax(i),'Number')
xlabel(ax(i),xlab)
fnm="FC_Branches";
saveas(f(i),fullfile(path,fnm+png));
%saveas(f(nb),path+fnm+fig);

i=i+1;
%9:
if bend
    %Bending results
    if tit;title(ax(i),'Evolution of Bending Event');end
    if strcmpi(x,'p');ax(i).XLim(1)=0;end
    ylabel(ax(i),'Number of bending events')
    xlabel(ax(i),xlab)
    fnm="FC_Bending";
    saveas(f(i),fullfile(path,fnm+png));
    i=i+1;
end

%10 to 12
if tit;title(ax(i),'Number of Force Chains');end
if strcmpi(x,'p');ax(i).XLim(1)=0;end
xlabel(ax(i),xlab)
if numel(pD)==1
    ylabel(ax(i),'Number')
    fnm="FC_Branch_Comp";
    saveas(f(i),fullfile(path,fnm+png));
    i=i+1;
    
    %tenth: Angle plot
    res=pD.Results.Elevation;
    dA=15;aCat=90/dA;
    ctc=cellfun(@length,res,'UniformOutput',false); %total nb of contacts
    ctc=cat(1,ctc{:});
    for a=1:aCat
        c=cellfun(@(x) sum(x>(a-1)*dA & x<=(a)*dA),res,'UniformOutput',false);
        c=cat(1,c{:});
        plotMark(app,ax(i),xAx,c,'Color',C(a,:),v{:});
        plotMark(app,ax(i+1),xAx,c./ctc,'Color',C(a,:),v{:});
    end
    
    lgStr(aCat)="";
    for a=1:aCat
        if a==aCat
            lgStr(a)=sprintf("[%05.2f , %05.2f]",(a-1)*dA,a*dA);
        else
            lgStr(a)=sprintf("[%05.2f , %05.2f[",(a-1)*dA,a*dA);
        end
    end

    if tit;title(ax(i),'Force chain elevation angle distribution');end
    if strcmpi(x,'p');ax(i).XLim(1)=0;end
    legend(ax(i),lgStr(:),'location','EastOutside')
    f(i).Position(3)=f(i).Position(3)*1.2;
    ylabel(ax(i),'Number of contacts')
    xlabel(ax(i),xlab)
    fnm="FC_Elevation_Nb";
    saveas(f(i),fullfile(path,fnm+png));
    
    if tit;title(ax(i+1),'Force chain elevation angle distribution');end
    if strcmpi(x,'p');ax(i+1).XLim(1)=0;end
    legend(ax(i+1),lgStr(:),'location','EastOutside')
    f(i+1).Position(3)=f(i+1).Position(3)*1.2;
    ylabel(ax(i+1),'Ratio of contacts')
    xlabel(ax(i+1),xlab)
    fnm="FC_Elevation_Rat";
    saveas(f(i+1),fullfile(path,fnm+png));
    
else
    ylabel(ax(i),'Number')
    fnm="FC_Sing_Branch";
    saveas(f(i),fullfile(path,fnm+png));
    i=i+1;
    %tenth :
    if tit;title(ax(i),'Number of Force Chains');end
    if strcmpi(x,'p');ax(i).XLim(1)=0;end
    ylabel(ax(i),'Number')
    xlabel(ax(i),xlab)
    fnm="FC_Multi_Branch";
    saveas(f(i),fullfile(path,fnm+png));
end

if nargin>2
    f(nb+1)=figure;ax(nb+1)=axes(f(nb+1));hold(ax(nb+1),'on');
    vals=varargin{1};
    scatter(ax(nb+1),1:numel(vals),vals,'x');
    if tit;title(ax(nb+1),'Ratio of time in force chain');end
    ylabel(ax(nb+1),'Ratio')
    xlabel(ax(nb+1),'Grain ID')
    fnm="Grain_Ration_FC";
    saveas(f(nb+1),fullfile(path,fnm+png));
end

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
function fcClstPlotter(app,pD,varargin)
%FCPLOTTER Plot data for Force Chain calculations

set(0,'defaultAxesFontSize',app.FontSizeEF.Value)
%Create Figures and Axis
nplg=9;
if size(pD.Results.Data,2)==4;nplg=nplg+9;end
f(nplg)=figure;ax(nplg)=axes(f(nplg));hold(ax(nplg),'on');
for i=1:(nplg-1)
    f(i)=figure;ax(i)=axes(f(i));hold(ax(i),'on'); %#ok<LAXES>
end

path=MakePath(app,'FCCL');
if nargin>2
    path=fullfile(path,varargin{1});
    if exist(path,'dir')==0;mkdir(path);end
end
png=".png";
xlab='Axial strain';
x='e';
for i=1:numel(pD)
    if pD(i).SimType==3
        xlab='Mean pressure (kPA)';
        x='p';
    end
end
%Check title and legends option
if app.TitlesCB.Value;tit=1;else;tit=0;end
if app.LegendsCB.Value;leg=1;else;leg=0;end

C=app.PlotColors;
if strcmpi(x,'p')
    xAx=pD(i).Results.Pressure;%strain
else
    xAx=pD(i).Results.Strain;%strain
end

%%%%%%%%%%%%%%%%%%%%%%%%% - PLOT DATA - %%%%%%%%%%%%%%%%%%%%%%%%%
lpData=pD.Results.Data(:,2);
lpData=cat(1,lpData{:});
totNb1=sum(lpData,3);
% clusters of each type
v={};
for k=1:numel(pD(1).InfPts.q)
    if strcmpi(x,'p')
        v=[v,{'Pointx',pD(1).InfPts.p(k)}]; %#ok<AGROW>
    else
        v=[v,{'Pointx',pD(1).InfPts.ez(k)}]; %#ok<AGROW>
    end
end
                    %%%% - NBPLOT - %%%%
j=1;
%plot FC cl category nb
plotMark(app,ax(j),xAx,lpData(:,1,1)./totNb1(:,1),'Color',C(1,:),v{:});
plotMark(app,ax(j),xAx,lpData(:,1,2)./totNb1(:,1),'Color',C(2,:),v{:});
plotMark(app,ax(j),xAx,lpData(:,1,3)./totNb1(:,1),'Color',C(3,:),v{:});
plotMark(app,ax(j),xAx,lpData(:,1,4)./totNb1(:,1),'Color',C(4,:),v{:});

j=j+1;
%plot FC cl category nb - NO 4
plotMark(app,ax(j),xAx,lpData(:,1,2)./totNb1(:,1),'Color',C(2,:),v{:});
plotMark(app,ax(j),xAx,lpData(:,1,3)./totNb1(:,1),'Color',C(3,:),v{:});
plotMark(app,ax(j),xAx,lpData(:,1,4)./totNb1(:,1),'Color',C(4,:),v{:});

j=j+1;
%plot NFC cl category nb
plotMark(app,ax(j),xAx,lpData(:,2,1)./totNb1(:,2),'Color',C(1,:),v{:});
plotMark(app,ax(j),xAx,lpData(:,2,2)./totNb1(:,2),'Color',C(2,:),v{:});
plotMark(app,ax(j),xAx,lpData(:,2,3)./totNb1(:,2),'Color',C(3,:),v{:});
plotMark(app,ax(j),xAx,lpData(:,2,4)./totNb1(:,2),'Color',C(4,:),v{:});

j=j+1;
%plot NFC cl category nb - NO 4
plotMark(app,ax(j),xAx,lpData(:,2,2)./totNb1(:,2),'Color',C(2,:),v{:});
plotMark(app,ax(j),xAx,lpData(:,2,3)./totNb1(:,2),'Color',C(3,:),v{:});
plotMark(app,ax(j),xAx,lpData(:,2,4)./totNb1(:,2),'Color',C(4,:),v{:});

                    %%%% - PERCELLPLOT - %%%%
lpData=pD.Results.Data(:,1);
lpData=cat(1,lpData{:});
totCl1=sum(lpData,3);
j=j+1;
%plot FC cl category PERCENTAGE
plotMark(app,ax(j),xAx,lpData(:,1,1)./totCl1(:,1),'Color',C(1,:),v{:});
plotMark(app,ax(j),xAx,lpData(:,1,2)./totCl1(:,1),'Color',C(2,:),v{:});
plotMark(app,ax(j),xAx,lpData(:,1,3)./totCl1(:,1),'Color',C(3,:),v{:});
plotMark(app,ax(j),xAx,lpData(:,1,4)./totCl1(:,1),'Color',C(4,:),v{:});

j=j+1;
%plot FC cl category PERCENTAGE - NO 4
plotMark(app,ax(j),xAx,lpData(:,1,2)./totCl1(:,1),'Color',C(2,:),v{:});
plotMark(app,ax(j),xAx,lpData(:,1,3)./totCl1(:,1),'Color',C(3,:),v{:});
plotMark(app,ax(j),xAx,lpData(:,1,4)./totCl1(:,1),'Color',C(4,:),v{:});

j=j+1;
%plot NFC cl category PERCENTAGE
plotMark(app,ax(j),xAx,lpData(:,2,1)./totCl1(:,2),'Color',C(1,:),v{:});
plotMark(app,ax(j),xAx,lpData(:,2,2)./totCl1(:,2),'Color',C(2,:),v{:});
plotMark(app,ax(j),xAx,lpData(:,2,3)./totCl1(:,2),'Color',C(3,:),v{:});
plotMark(app,ax(j),xAx,lpData(:,2,4)./totCl1(:,2),'Color',C(4,:),v{:});

j=j+1;
%plot NFC cl category PERCENTAGE - NO 4
plotMark(app,ax(j),xAx,lpData(:,2,2)./totCl1(:,2),'Color',C(2,:),v{:});
plotMark(app,ax(j),xAx,lpData(:,2,3)./totCl1(:,2),'Color',C(3,:),v{:});
plotMark(app,ax(j),xAx,lpData(:,2,4)./totCl1(:,2),'Color',C(4,:),v{:});


        %%%%%%%% - Calculation 2 - %%%%%%%%
if size(pD.Results.Data,2)==4
                        %%%% - NbPlot - %%%%
    lpData=pD.Results.Data(:,4);
    lpData=cat(1,lpData{:});
    totNb2=sum(lpData,3);
    j=j+1;
    %plot FC nb category PERCENTAGE
    plotMark(app,ax(j),xAx,lpData(:,1,1)./totNb2(:,1),'Color',C(1,:),v{:});
    plotMark(app,ax(j),xAx,lpData(:,1,2)./totNb2(:,1),'Color',C(2,:),v{:});
    plotMark(app,ax(j),xAx,lpData(:,1,3)./totNb2(:,1),'Color',C(3,:),v{:});
    plotMark(app,ax(j),xAx,lpData(:,1,4)./totNb2(:,1),'Color',C(4,:),v{:});
    
    j=j+1;
    %plot FC nb category PERCENTAGE - NO 4
    plotMark(app,ax(j),xAx,lpData(:,1,2)./totNb2(:,1),'Color',C(2,:),v{:});
    plotMark(app,ax(j),xAx,lpData(:,1,3)./totNb2(:,1),'Color',C(3,:),v{:});
    plotMark(app,ax(j),xAx,lpData(:,1,4)./totNb2(:,1),'Color',C(4,:),v{:});
    
    j=j+1;
    %plot NFC nb category PERCENTAGE
    plotMark(app,ax(j),xAx,lpData(:,2,1)./totNb2(:,2),'Color',C(1,:),v{:});
    plotMark(app,ax(j),xAx,lpData(:,2,2)./totNb2(:,2),'Color',C(2,:),v{:});
    plotMark(app,ax(j),xAx,lpData(:,2,3)./totNb2(:,2),'Color',C(3,:),v{:});
    plotMark(app,ax(j),xAx,lpData(:,2,4)./totNb2(:,2),'Color',C(4,:),v{:});
    
    j=j+1;
    %plot NFC nb category PERCENTAGE - NO 4
    plotMark(app,ax(j),xAx,lpData(:,2,2)./totNb2(:,2),'Color',C(2,:),v{:});
    plotMark(app,ax(j),xAx,lpData(:,2,3)./totNb1(:,2),'Color',C(3,:),v{:});
    plotMark(app,ax(j),xAx,lpData(:,2,4)./totNb2(:,2),'Color',C(4,:),v{:});

                        %%%% - PERCELLPLOT - %%%%
    lpData=pD.Results.Data(:,3);
    lpData=cat(1,lpData{:});
    totCl2=sum(lpData,3);
    j=j+1;
    %plot FC cl category PERCENTAGE
    plotMark(app,ax(j),xAx,lpData(:,1,1)./totCl2(:,1),'Color',C(1,:),v{:});
    plotMark(app,ax(j),xAx,lpData(:,1,2)./totCl2(:,1),'Color',C(2,:),v{:});
    plotMark(app,ax(j),xAx,lpData(:,1,3)./totCl2(:,1),'Color',C(3,:),v{:});
    plotMark(app,ax(j),xAx,lpData(:,1,4)./totCl2(:,1),'Color',C(4,:),v{:});
    
    j=j+1;
    %plot FC cl category PERCENTAGE - NO 4
    plotMark(app,ax(j),xAx,lpData(:,1,2)./totCl2(:,1),'Color',C(2,:),v{:});
    plotMark(app,ax(j),xAx,lpData(:,1,3)./totCl2(:,1),'Color',C(3,:),v{:});
    plotMark(app,ax(j),xAx,lpData(:,1,4)./totCl2(:,1),'Color',C(4,:),v{:});
    
    j=j+1;
    %plot NFC cl category PERCENTAGE
    plotMark(app,ax(j),xAx,lpData(:,2,1)./totCl2(:,2),'Color',C(1,:),v{:});
    plotMark(app,ax(j),xAx,lpData(:,2,2)./totCl2(:,2),'Color',C(2,:),v{:});
    plotMark(app,ax(j),xAx,lpData(:,2,3)./totCl2(:,2),'Color',C(3,:),v{:});
    plotMark(app,ax(j),xAx,lpData(:,2,4)./totCl2(:,2),'Color',C(4,:),v{:});
    
    j=j+1;
    %plot NFC cl category PERCENTAGE - NO 4
    plotMark(app,ax(j),xAx,lpData(:,2,2)./totCl2(:,2),'Color',C(2,:),v{:});
    plotMark(app,ax(j),xAx,lpData(:,2,3)./totCl2(:,2),'Color',C(3,:),v{:});
    plotMark(app,ax(j),xAx,lpData(:,2,4)./totCl2(:,2),'Color',C(4,:),v{:});

end

                    %%%% - RATIO OF FCC in relation to all - %%%%
        %%%%%%%% - Calculation 1 - %%%%%%%%
j=j+1;
plotMark(app,ax(j),xAx,totNb1(:,1)./sum(totNb1,2),'Color',C(1,:),v{:});
plotMark(app,ax(j),xAx,totCl1(:,1)./sum(totCl1,2),'Color',C(2,:),v{:});
        %%%%%%%% - Calculation 2 - %%%%%%%%
if size(pD.Results.Data,2)==4
    j=j+1;
    %plot the ratio of FCC in number and in cell in relation to total amount
    plotMark(app,ax(j),xAx,totNb2(:,1)./sum(totNb2,2),'Color',C(1,:),v{:});
    plotMark(app,ax(j),xAx,totCl2(:,1)./sum(totCl2,2),'Color',C(2,:),v{:});
end

%%%%%%%%%%%%%%%%%%%%%%%%% - LEGEND AND SAVE - %%%%%%%%%%%%%%%%%%%%%%%%%
if leg
    for j=1:min(nplg-1,16)
        if rem(j,2)==1
            legend(ax(j),"Order 4","Order 6","Order 8-20",...
                "Order 22+",'location','eastoutside')
        else
            legend(ax(j),"Order 6","Order 8-20","Order 22+",...
                'location','eastoutside')
        end
    end
end

                    %%%% - Nb plot 1 - %%%%
j=1;
%Prepare axis to get the same value between NFC and FC equal graphs
lim=ax(j).XLim;
%NUMBER - FCC
if tit;title(ax(j),'Distribution of FCC nb');end
ylabel(ax(j),'Ratio of FCC (number)')
ax(j).YLim =[0,1];
if strcmpi(x,'e')
    xticks(ax(j),lim(1):0.05:lim(2))
end
xlabel(ax(j),xlab)
fnm="FCC_Nb";
saveas(f(j),fullfile(path,fnm+png));
j=j+1;

%NUMBER - FCC - NO4
if tit;title(ax(j),'Distribution of FCC nb');end
ax(j).YLim =[0,max(ax(j).YLim(2),ax(j+2).YLim(2))];
if strcmpi(x,'e')
    xticks(ax(j),lim(1):0.05:lim(2))
end
ylabel(ax(j),'Ratio of FCC (number)')
xlabel(ax(j),xlab)
fnm="FCC_Nb_No4";
saveas(f(j),fullfile(path,fnm+png));
j=j+1;

%NUMBER - NFCC
if tit;title(ax(j),'Distribution of NFCC nb');end
ax(j).YLim =[0,1];
if strcmpi(x,'e')
    xticks(ax(j),lim(1):0.05:lim(2))
end
ylabel(ax(j),'Ratio of NFCC (number)')
xlabel(ax(j),xlab)
fnm="NFCC_Nb";
saveas(f(j),fullfile(path,fnm+png));
j=j+1;

%NUMBER - NFCC - NO4
if tit;title(ax(j),'Distribution of NFCC nb');end
ax(j).YLim =[0,max(ax(j).YLim(2),ax(j-2).YLim(2))];
if strcmpi(x,'e')
    xticks(ax(j),lim(1):0.05:lim(2))
end
ylabel(ax(j),'Ratio of NFCC (number)')
xlabel(ax(j),xlab)
fnm="NFCC_Nb_No4";
saveas(f(j),fullfile(path,fnm+png));
j=j+1;

                    %%%% - Cell plot 1 - %%%%
%CELL - FCC
if tit;title(ax(j),'Distribution of FCC cell');end
ax(j).YLim =[0,1];
if strcmpi(x,'e')
    xticks(ax(j),lim(1):0.05:lim(2))
end
ylabel(ax(j),'Ratio of FCC (cells)')
xlabel(ax(j),xlab)
fnm="FCC_Cell";
saveas(f(j),fullfile(path,fnm+png));
j=j+1;

%CELL - FCC - NO4
if tit;title(ax(j),'Distribution of FCC cell');end
ax(j).YLim =[0,max(ax(j).YLim(2),ax(j+2).YLim(2))];
if strcmpi(x,'e')
    xticks(ax(j),lim(1):0.05:lim(2))
end
ylabel(ax(j),'Ratio of FCC (cells)')
xlabel(ax(j),xlab)
fnm="FCC_Cell_No4";
saveas(f(j),fullfile(path,fnm+png));
j=j+1;

%CELL - NFCCE
if tit;title(ax(j),'Distribution of NFCC cell');end
ax(j).YLim =[0,1];
if strcmpi(x,'e')
    xticks(ax(j),lim(1):0.05:lim(2))
end
ylabel(ax(j),'Ratio of NFCC (cells)')
xlabel(ax(j),xlab)
fnm="NFCC_Cell";
saveas(f(j),fullfile(path,fnm+png));
j=j+1;

%CELL - NFCC - NO4
if tit;title(ax(j),'Distribution of NFCC cell');end
ax(j).YLim =[0,max(ax(j).YLim(2),ax(j-2).YLim(2))];
if strcmpi(x,'e')
    xticks(ax(j),lim(1):0.05:lim(2))
end
ylabel(ax(j),'Ratio of NFCC (cells)')
xlabel(ax(j),xlab)
fnm="NFCC_Cell_No4";
saveas(f(j),fullfile(path,fnm+png));
j=j+1;

if size(pD.Results.Data,2)==4
    %%%% - NB plot 2 - %%%%
    %Nb2 - FCC
    if tit;title(ax(j),'Distribution of FCC nb (2)');end
    ax(j).YLim =[0,1];
    if strcmpi(x,'e')
        xticks(ax(j),lim(1):0.05:lim(2))
    end
    ylabel(ax(j),'Ratio of FCC (number)')
    xlabel(ax(j),xlab)
    fnm="FCC_Nb2";
    saveas(f(j),fullfile(path,fnm+png));
    j=j+1;

    %Nb2 - FCC - NO4
    if tit;title(ax(j),'Distribution of FCC nb (2)');end
    ax(j).YLim =[0,max(ax(j).YLim(2),ax(j+2).YLim(2))];
    if strcmpi(x,'e')
        xticks(ax(j),lim(1):0.05:lim(2))
    end
    ylabel(ax(j),'Ratio of FCC (number)')
    xlabel(ax(j),xlab)
    fnm="FCC_Nb2_No4";
    saveas(f(j),fullfile(path,fnm+png));
    j=j+1;

    %Nb2 - NFCC
    if tit;title(ax(j),'Distribution of NFCC nb (2)');end
    ax(j).YLim =[0,1];
    if strcmpi(x,'e')
        xticks(ax(j),lim(1):0.05:lim(2))
    end
    ylabel(ax(j),'Ratio of NFCC (number)')
    xlabel(ax(j),xlab)
    fnm="NFCC_Nb2";
    saveas(f(j),fullfile(path,fnm+png));
    j=j+1;

    %Nb2 - NFCC - NO4
    if tit;title(ax(j),'Distribution of NFCC nb (2)');end
    ax(j).YLim =[0,max(ax(j).YLim(2),ax(j-2).YLim(2))];
    if strcmpi(x,'e')
        xticks(ax(j),lim(1):0.05:lim(2))
    end
    ylabel(ax(j),'Ratio of NFCC (number)')
    xlabel(ax(j),xlab)
    fnm="NFCC_Nb2_No4";
    saveas(f(j),fullfile(path,fnm+png));
    j=j+1;

    %%%% - Cell plot 2 - %%%%

    %CELL2 - FCC
    if tit;title(ax(j),'Distribution of FCC cells (2)');end
    ax(j).YLim =[0,1];
    if strcmpi(x,'e')
        xticks(ax(j),lim(1):0.05:lim(2))
    end
    ylabel(ax(j),'Ratio of FCC (cells)')
    xlabel(ax(j),xlab)
    fnm="FCC_Cell2";
    saveas(f(j),fullfile(path,fnm+png));
    j=j+1;

    %CELL2 - FCC - NO4
    if tit;title(ax(j),'Distribution of FCC cells (2)');end
    ax(j).YLim =[0,max(ax(j).YLim(2),ax(j+2).YLim(2))];
    if strcmpi(x,'e')
        xticks(ax(j),lim(1):0.05:lim(2))
    end
    ylabel(ax(j),'Ratio of FCC (cells)')
    xlabel(ax(j),xlab)
    fnm="FCC_Cell2_No4";
    saveas(f(j),fullfile(path,fnm+png));
    j=j+1;

    %CELL2 - NFCCE
    if tit;title(ax(j),'Distribution of NFCC cells (2)');end
    ax(j).YLim =[0,1];
    if strcmpi(x,'e')
        xticks(ax(j),lim(1):0.05:lim(2))
    end
    ylabel(ax(j),'Ratio of NFCC (cells)')
    xlabel(ax(j),xlab)
    fnm="NFCC_Cell2";
    saveas(f(j),fullfile(path,fnm+png));
    j=j+1;

    %CELL2 - NFCC - NO4
    if tit;title(ax(j),'Distribution of NFCC cells (2)');end
    ax(j).YLim =[0,max(ax(j).YLim(2),ax(j-2).YLim(2))];
    if strcmpi(x,'e')
        xticks(ax(j),lim(1):0.05:lim(2))
    end
    ylabel(ax(j),'Ratio of NFCC (cells)')
    xlabel(ax(j),xlab)
    fnm="NFCC_Cell2_No4";
    saveas(f(j),fullfile(path,fnm+png));

    %Ratio of FCC nb and cells
    j=j+2;
    ax(j).YLim =[0,max(ax(j).YLim(2),ax(j-1).YLim(2))];
    ax(j-1).YLim =[0,max(ax(j).YLim(2),ax(j-1).YLim(2))];
    if tit;title(ax(j),'Ratio of FCC Nb and Cells');end
    if strcmpi(x,'e')
        xticks(ax(j),lim(1):0.05:lim(2))
    end
    legend(ax(j),'Number','Cells','Location','southeast')
    ylabel(ax(j),'Ratio of FCC')
    xlabel(ax(j),xlab)
    fnm="Ratio_FCC2";
    saveas(f(j),fullfile(path,fnm+png));
    

end

j=j-1;
%Ratio of FCC nb and cells
if tit;title(ax(j),'Ratio of FCC Nb and Cells');end
if strcmpi(x,'e')
    xticks(ax(j),lim(1):0.05:lim(2))
end
legend(ax(j),'Number','Cells','Location','southeast')
ylabel(ax(j),'Ratio of FCC')
xlabel(ax(j),xlab)
fnm="Ratio_FCC";
saveas(f(j),fullfile(path,fnm+png));

%outside legend
if ~leg
    %Order 4+ Legend
    p=copyobj(f(1),0);
    l=legend(p.CurrentAxes,"Order 4","Order 6","Order 8-20",...
            "Order 22+",'Orientation','horizontal');
    l.EdgeColor='none';
    set(p.CurrentAxes,'Visible','Off')
    set(p,'Position',(get(l,'Position').*[0, 0, 1, 1].*get(p,'Position')));
    set(l,'Position',[0,0,1,1]);
    set(p, 'Position', get(p,'Position') + [500, 400, 0, 0]);
    saveas(p,fullfile(path,"Legend_Order"+png));
    
    %order 6+Legend
    q=copyobj(f(2),0);
    l=legend(q.CurrentAxes,"Order 6","Order 8-20",...
            "Order 22+",'Orientation','horizontal');
    l.EdgeColor='none';
    set(q.CurrentAxes,'Visible','Off')
    set(q,'Position',(get(l,'Position').*[0, 0, 1, 1].*get(q,'Position')));
    set(l,'Position',[0,0,1,1]);
    set(q, 'Position', get(q,'Position') + [500, 400, 0, 0]);
    saveas(q,fullfile(path,"Legend_Order_No4"+png));
    
    %Delete extra figures
    delete([q,p]);
end

%Delete figures in the case of exe all
if isequal(app.LIGGGHTSAnalysisButtonGroup.SelectedObject,...
        app.ExeAllButton) || nargin>2
    %delete(f);
    if numel(pD)==1;delete(f);end
end
end
function fcTransfPlotter(pD,app)
%CLSTRPLOTTER plot the evolution the Cl4 and Cl6 transformations

%linetype

%Check title and legends option
if app.TitlesCB.Value;tit=1;else;tit=0;end
if app.LegendsCB.Value;leg=1;else;leg=0;end

path=MakePath(app,'FCCLTF');
png='.png';
C=app.PlotColors;
set(0,'defaultAxesFontSize',app.FontSizeEF.Value)
%prepare variables for execution
ctRes=pD.Results;
cl=zeros(size(ctRes.ClTransf,1),3,4);%cluster data 
nbCll=NaN(size(ctRes.ClTransf,1),4); %nb of cells
for i=2:size(ctRes.ClTransf,1)
    %get the results and separate by which category they become
    r=ctRes.ClTransf{i,1};
    %r.Cl4 and r.Cl6 [NOrder,OldID,NewID]
    %while r.Cl8 [OOrder,NOrder,OldID, NewID]
    cl4=[(r.Cl4(:,1)==6),...
        (r.Cl4(:,1)>6 & r.Cl4(:,1)<22),...
    	(r.Cl4(:,1)>20)];
    cl6=[(r.Cl6(:,1)==4),...
    	(r.Cl6(:,1)>6 & r.Cl6(:,1)<22),...
    	(r.Cl6(:,1)>20)];
    cat8=(r.Cl8(:,1)>6 & r.Cl8(:,1)<22);
    cl8=[(r.Cl8(cat8,2)==4),...
    	(r.Cl8(cat8,2)==6),...
    	(r.Cl8(cat8,2)>20)];
    cl22=[(r.Cl8(~cat8,2)==4),...
    	(r.Cl8(~cat8,2)==6),...
    	(r.Cl8(~cat8,2)>6 & r.Cl8(~cat8,2)<22)];
    %all transformations
    cl(i,:,:)=cat(3,sum(cl4,1),sum(cl6,1),sum(cl8,1),sum(cl22,1));
    nbCll(i,:)=r.NbC;
end

plG1=2; %plot group 1 => number of groups of tiled figures
nb=4*(plG1);
f(4*(plG1))=figure;
f(4*(plG1)).Position(3)=2*f(4*(plG1)).Position(3);
tf(4*(plG1))=tiledlayout(f(4*(plG1)),1,2,'TileSpacing',...
    'tight','Padding','Compact');
for i=1:(4*(plG1)-1)
    f(i)=figure;
    f(i).Position(3)=2*f(i).Position(3);
    tf(i)=tiledlayout(f(i),1,2,'TileSpacing',...
        'tight','Padding','Compact');
end

i=1;
%% PLOTGROUP 1 transformatio ratio
%Cl4 => others  & others => Cl4 - - ratio
ax(i,1)=nexttile(tf(i));
plotMark(app,ax(i,1),ctRes.Strain,cl(:,1,1)./nbCll(:,1),'Color',C(2,:))
plotMark(app,ax(i,1),ctRes.Strain,cl(:,2,1)./nbCll(:,1),'Color',C(3,:))
plotMark(app,ax(i,1),ctRes.Strain,cl(:,3,1)./nbCll(:,1),'Color',C(4,:));
ax(i,2)=nexttile(tf(i)); 
plotMark(app,ax(i,2),ctRes.Strain,cl(:,1,2)./nbCll(:,1),'Color',C(2,:))
plotMark(app,ax(i,2),ctRes.Strain,cl(:,1,3)./nbCll(:,1),'Color',C(3,:))
plotMark(app,ax(i,2),ctRes.Strain,cl(:,1,4)./nbCll(:,1),'Color',C(4,:))
i=i+1;
%Cl6 => others & others => Cl6 - - ratio
ax(i,1)=nexttile(tf(i));
plotMark(app,ax(i,1),ctRes.Strain,cl(:,1,2)./nbCll(:,2),'Color',C(1,:))
plotMark(app,ax(i,1),ctRes.Strain,cl(:,2,2)./nbCll(:,2),'Color',C(3,:))
plotMark(app,ax(i,1),ctRes.Strain,cl(:,3,2)./nbCll(:,2),'Color',C(4,:));
ax(i,2)=nexttile(tf(i));
plotMark(app,ax(i,2),ctRes.Strain,cl(:,1,1)./nbCll(:,2),'Color',C(1,:))
plotMark(app,ax(i,2),ctRes.Strain,cl(:,2,3)./nbCll(:,2),'Color',C(3,:))
plotMark(app,ax(i,2),ctRes.Strain,cl(:,2,4)./nbCll(:,2),'Color',C(4,:));
i=i+1;
%Cl8 => others & others => Cl8 - - ratio
ax(i,1)=nexttile(tf(i));
plotMark(app,ax(i,1),ctRes.Strain,cl(:,1,3)./nbCll(:,3),'Color',C(1,:))
plotMark(app,ax(i,1),ctRes.Strain,cl(:,2,3)./nbCll(:,3),'Color',C(2,:))
plotMark(app,ax(i,1),ctRes.Strain,cl(:,3,3)./nbCll(:,3),'Color',C(4,:));
ax(i,2)=nexttile(tf(i));
plotMark(app,ax(i,2),ctRes.Strain,cl(:,2,1)./nbCll(:,3),'Color',C(1,:))
plotMark(app,ax(i,2),ctRes.Strain,cl(:,2,2)./nbCll(:,3),'Color',C(2,:))
plotMark(app,ax(i,2),ctRes.Strain,cl(:,3,4)./nbCll(:,3),'Color',C(4,:));
i=i+1;
%Cl22 => others & others => Cl22 - - ratio
ax(i,1)=nexttile(tf(i));
plotMark(app,ax(i,1),ctRes.Strain,cl(:,1,4)./nbCll(:,4),'Color',C(1,:))
plotMark(app,ax(i,1),ctRes.Strain,cl(:,2,4)./nbCll(:,4),'Color',C(2,:))
plotMark(app,ax(i,1),ctRes.Strain,cl(:,3,4)./nbCll(:,4),'Color',C(3,:));
ax(i,2)=nexttile(tf(i));
plotMark(app,ax(i,2),ctRes.Strain,cl(:,3,1)./nbCll(:,4),'Color',C(1,:))
plotMark(app,ax(i,2),ctRes.Strain,cl(:,3,2)./nbCll(:,4),'Color',C(2,:))
plotMark(app,ax(i,2),ctRes.Strain,cl(:,3,3)./nbCll(:,4),'Color',C(3,:));
i=i+1;

%% PLOTGROUP 2 transformation per nb
%Cl4 => others  & others => Cl4 - - nb
ax(i,1)=nexttile(tf(i));
plotMark(app,ax(i,1),ctRes.Strain,cl(:,1,1),'Color',C(2,:))
plotMark(app,ax(i,1),ctRes.Strain,cl(:,2,1),'Color',C(3,:))
plotMark(app,ax(i,1),ctRes.Strain,cl(:,3,1),'Color',C(4,:));
ax(i,2)=nexttile(tf(i)); 
plotMark(app,ax(i,2),ctRes.Strain,cl(:,1,2),'Color',C(2,:))
plotMark(app,ax(i,2),ctRes.Strain,cl(:,1,3),'Color',C(3,:))
plotMark(app,ax(i,2),ctRes.Strain,cl(:,1,4),'Color',C(4,:))
i=i+1;
%Cl6 => others & others => Cl6 - - nb
ax(i,1)=nexttile(tf(i));
plotMark(app,ax(i,1),ctRes.Strain,cl(:,1,2),'Color',C(1,:))
plotMark(app,ax(i,1),ctRes.Strain,cl(:,2,2),'Color',C(3,:))
plotMark(app,ax(i,1),ctRes.Strain,cl(:,3,2),'Color',C(4,:));
ax(i,2)=nexttile(tf(i));
plotMark(app,ax(i,2),ctRes.Strain,cl(:,1,1),'Color',C(1,:))
plotMark(app,ax(i,2),ctRes.Strain,cl(:,2,3),'Color',C(3,:))
plotMark(app,ax(i,2),ctRes.Strain,cl(:,2,4),'Color',C(4,:));
i=i+1;
%Cl8 => others & others => Cl8 - - nb
ax(i,1)=nexttile(tf(i));
plotMark(app,ax(i,1),ctRes.Strain,cl(:,1,3),'Color',C(1,:))
plotMark(app,ax(i,1),ctRes.Strain,cl(:,2,3),'Color',C(2,:))
plotMark(app,ax(i,1),ctRes.Strain,cl(:,3,3),'Color',C(4,:));
ax(i,2)=nexttile(tf(i));
plotMark(app,ax(i,2),ctRes.Strain,cl(:,2,1),'Color',C(1,:))
plotMark(app,ax(i,2),ctRes.Strain,cl(:,2,2),'Color',C(2,:))
plotMark(app,ax(i,2),ctRes.Strain,cl(:,3,4),'Color',C(4,:));
i=i+1;
%Cl22 => others & others => Cl22 - - nb
ax(i,1)=nexttile(tf(i));
plotMark(app,ax(i,1),ctRes.Strain,cl(:,1,4),'Color',C(1,:))
plotMark(app,ax(i,1),ctRes.Strain,cl(:,2,4),'Color',C(2,:))
plotMark(app,ax(i,1),ctRes.Strain,cl(:,3,4),'Color',C(3,:));
ax(i,2)=nexttile(tf(i));
plotMark(app,ax(i,2),ctRes.Strain,cl(:,3,1),'Color',C(1,:))
plotMark(app,ax(i,2),ctRes.Strain,cl(:,3,2),'Color',C(2,:))
plotMark(app,ax(i,2),ctRes.Strain,cl(:,3,3),'Color',C(3,:));

%% Titles and stuff
legs={["Cl6","Cl8-20","Cl22+"];
    ["Cl4","Cl8-20","Cl22+"];
    ["Cl4","Cl6","Cl22+"];
    ["Cl4","Cl6","Cl8-20"]};
ylab=["Ratio of Events";
    "Number of Events"];
titles=["Order 4 to Others","Others to 4";
    "Order 6 to Others","Others to 6";
    "Order 8-20 to Others","Others to 8-20";
    "Order 22+ to Others","Others to 22+";
    "Order 4 to Others","Others to 4";
    "Order 6 to Others","Others to 6";
    "Order 8-20 to Others","Others to 8-20";
    "Order 22+ to Others","Others to 22+"];
fnm=["Transf_O4";"Transf_Order6";
    "Transf_O8";"Transf_O22";
    "Transf_O4_Nb";"Transf_O6_Nb";
    "Transf_O8_Nb";"Transf_O22_Nb"];
for i=1:nb
    %Add legend in the for the tiledlayouts only
    if leg
        if i<=4*(plG1)
            if isprop(ax(i,2),'Xlim')
                lgd=legend(ax(i,1),legs{i-4*(ceil(i/4)-1)},'NumColumns',3);
                lgd.Layout.Tile = 'south';
            else
                legend(ax(i,1),legs{i-4*(ceil(i/4)-1)},'Location','best');
            end
        else
            legend(ax(i,1),legs{5},'Location','best');
        end
    end
    %add a horizontal line on zero
    yl=yline(ax(i,1),0,':','Color','#C1C1C1');
    set(get(get(yl,'Annotation'),'LegendInformation'),...
        'IconDisplayStyle','off');
    if isprop(ax(i,2),'Xlim')   %if tiled layout add on the other too
        yl=yline(ax(i,2),0,':','Color','#C1C1C1');
        set(get(get(yl,'Annotation'),'LegendInformation'),...
            'IconDisplayStyle','off');
    end
    if tit;title(ax(i,1),titles(i,1));end
    ylabel(ax(i,1),ylab(ceil(i/4)))
    xlabel(ax(i,1),'Axial Strain')
    if isprop(ax(i,2),'Xlim')
        if tit;title(ax(i,2),titles(i,2));end
        ylabel(ax(i,2),ylab(ceil(i/4)))
        xlabel(ax(i,2),'Axial Strain')
        %put them in the same vertical scale
        mx=max(ax(i,1).YLim(2),ax(i,2).YLim(2)); 
        ax(i,1).YLim(2)=mx;
        ax(i,2).YLim(2)=mx;
    end
    saveas(f(i),fullfile(path,fnm(i)+png));
end

if isequal(app.LIGGGHTSAnalysisButtonGroup.SelectedObject,...
        app.ExeAllButton)
    delete(f);
end
end
%execution function
function fcData=fcBase(fc,gr,step,pathVTK)
%FCNORMAL force chain normal calculation
% Uses previously calculated force chains to extract some data and plot vtk
% files containing information related to each grain stress and fc
% location.
id=cat(1,fc.IDs);
%write VTK files
sV=cat(1,gr.SingleGrain(id).StressVector);
if size(sV,2)<3
    sV=[zeros(numel(id),1) sV];
end
vtkwrite(pathVTK+"connectionsFC"+step+".vtk",'polydata',"LINESB",...
    gr.Coord(:,1),gr.Coord(:,2),gr.Coord(:,3),cat(1,fc.Lines));
vtkwrite(pathVTK+"grainsFC"+step+".vtk",'unstructured_grid',...
    gr.Coord(id,1),gr.Coord(id,2),gr.Coord(id,3),...
    'SCALARS','Radius',gr.Radius(id),...
    'SCALARS','Stress',cat(1,gr.SingleGrain(id).PrincipalStress),...
    'VECTORS','StressDir',sV(:,1),sV(:,2),sV(:,3),...
    'PRECISION',10);

sbFC=[fc.NbBranches]==1;
l=[fc.Length];
mnSB=mean(l(sbFC));
mnMB=mean(l(~sbFC));
%save important data
nfcg=numel(unique(cat(1,fc.IDs)));
fcData=[step mnMB mnSB numel(fc)...
    sum([fc.NbBranches]) sum(sbFC)...
    nfcg/gr.HighStrGrains nfcg/gr.Nb size(gr.ThreeP,1)];
end
function [fcClCell,fcClNb,fcClCell2,fcClNb2]=fcClusters(app,fc,gr,sc,step,pathVTK)
%FCCLUSTER join force chain an cluster calculation data
% Two calculations will be executed. The first considers clusters as part of 
% force chain if they share a closed contact. Second considers cluster any 
% cluster sharing a grain with the FC as belonging to it.
%
% This function will use a spaceCellSystem object that contains the data
% from a previous cluster calculation and match with a forcechain object to
% identify the clusters that make part of the latter.

%General values
fcGrains=sort(cat(1,fc.IDs));       %identify grains forming the force chain
nbCl=numel(sc.Loops);               %nb of clusters
nbCl4=sc.Clt4.nbCells;              %nb of clusters 4
Order=cat(1,sc.Loops.Order);        %cluster order
Cells=cat(1,sc.Loops.nbCells);      %ID of cluster cells
Size=cat(1,sc.Loops.Size);          %size of clusters

%%%%%%%%%%%%%%%%% Prepare first analysis %%%%%%%%%%%%%%%%%
%Create a list of all closed edges forming clusters, with the ID of the
%cluster in the third column and order as fourth
    %count the nb of closed edges per cluster
l=cell(nbCl,1);
[l{:}]=sc.Loops.ClEdges;
[nbCEdges,~]=cellfun(@size,l);
    %create a GrID1-GrID2-ClID-Order-nbCells vector for each closed contact
clInfo=[(1:nbCl)',Order,Cells]; %ClID-Order-nbCells
clCE=[cat(1,sc.Loops.ClEdges) repelem(clInfo,nbCEdges,1)]; %GrID1-GrID2-ClID-Order-nbCells
    %add cl4 to the analysis. All edges of Cl4 will be added and not only
    %the closed ones. But as it will be later compared to closed edges
    %belonging to force chains it will be automatically filtered
    %identify clusters4 edges
if isempty(sc.GoodCells);sc.GoodCells = goodCell(sc);end
cl4Id=sc.Clt4.sCells; %id Cl4 cell
cl4grs=sc.DelaunayT(cl4Id,:); %cl4 grains
cb=nchoosek(1:4,2); %combination
clCE4=[cl4grs(:,cb(1,:));...
    cl4grs(:,cb(2,:));...
    cl4grs(:,cb(3,:));...
    cl4grs(:,cb(4,:));...
    cl4grs(:,cb(5,:));...
    cl4grs(:,cb(6,:))];
    %creat a GrID1 - GrID2 - ClID - ClOrder vector for each cl4 edge
clGr4=[(nbCl+1:nbCl+nbCl4);ones(1,nbCl4)*4;ones(1,nbCl4)]';
clCE4=[clCE4 repmat(clGr4,6,1)];

    %merge both vectors
clCE=[clCE;clCE4];
    %only keep lines that contain grains from force chain
clCE=clCE(ismember(clCE(:,1),fcGrains) & ismember(clCE(:,2),fcGrains),:);

%%%%%%%%%%%%%%%%% Prepare second analysis %%%%%%%%%%%%%%%%%
%Create a list that resembles the previous ones, but for each individual
%grain insted of the contact %GrID1 - ClID - ClOrder 
clCEv2=[cat(1,sc.Loops.Grains) repelem(clInfo,Size,1)]; %GrID-ClID-Order-nbCells
    %add cl4s - first transform cl4gr matrix  in a vector of gr IDs
cl4grs=reshape(cl4grs',[],1);
    %for each gr repeat Cell ID Cell Order
clCE4=[cl4grs repmat(clGr4,4,1)];
    %join vectors
clCEv2=[clCEv2;clCE4];
    %only keep lines with FC grains
clCEv2=clCEv2(ismember(clCEv2(:,1),fcGrains),:);

%for each force chain, compare the force chain contacts with clCE matrix
%to identify the ID of the clusters forming it.
fcChk=ones(numel(fc),1);
fcLpInf=cell(numel(fc),1);    %stores ID,Order and NbCells of each loop per FC for calculation 2
for i=1:numel(fc)
    %%%% first calculation %%%%
    lp=ismember(clCE(:,1:2),fc(i).Lines,'rows');
    if isempty(lp)
        %if there are no clusters in the force chains, the fc must be
        %located around the wall and surrounded by 'fakegrains'. Thus it
        %should be removed from the calculation
        fcChk(i)=1;
    else
        lp=unique(clCE(lp,3:end),'rows');  %ClusterID - ClusterOrder - NbCells
        %max - min - mean order values
        fc(i).ClustersVals=[max(lp(:,2)),min(lp(:,2)),mean(lp(:,2))];
        fc(i).ClustersID=lp;
    end
    %%%% Second Calculation %%%%
    lp=ismember(clCE(:,1),fc(i).IDs); %compare grains with the list
    fcLpInf{i}=unique(clCEv2(lp,2:end),'rows');  %get unique lines
end

%remove excess data from the bad force chains
fcChk=logical(fcChk);       %turn into logical
fc=fc(fcChk);               %delete bad fc
idNew=sort(cat(1,fc.IDs));  %get the correct grain IDs

%vtk save
vtkwrite(pathVTK(1)+"connectionsFC"+step+".vtk",'polydata',"LINESB",...
    gr.Coord(:,1),gr.Coord(:,2),gr.Coord(:,3),cat(1,fc.Lines));
vtkwrite(pathVTK(1)+"grainsFC"+step+".vtk",'unstructured_grid',...
    gr.Coord(idNew,1),gr.Coord(idNew,2),gr.Coord(idNew,3),...
    'SCALARS','Radius',gr.Radius(idNew),...
    'PRECISION',10);

%FCC - count force chain clusters
    %get order of clusters beglonging to FC
fcClNb=unique(cat(1,fc.ClustersID),'rows');
results=zeros(1,4,4);

%small category FCC nb - NFCC nb - FCC cell - NFCC cell
chk=fcClNb(:,2)==4;
results(:,:,1)=[sum(chk),sc.Clt4.nbCells,sum(fcClNb(chk,3)),sc.Clt4.nbCells]; 
results(:,[2,4],1)=results(:,[2,4],1)-results(:,[1,3],1); 

%submedium category FCC nb - NFCC nb - FCC cell - NFCC cell
chk=fcClNb(:,2)==6;
chk2=Order==6;
results(:,:,2)=[sum(chk),sum(chk2),sum(fcClNb(chk,3)),sum(Cells(chk2))];
results(:,[2,4],2)=results(:,[2,4],2)-results(:,[1,3],2); 

%medium category FCC nb - NFCC nb - FCC cell - NFCC cell
chk=fcClNb(:,2)>6 & fcClNb(:,2)<21;
chk2=Order>6 & Order<21;
results(:,:,3)=[sum(chk),sum(chk2),sum(fcClNb(chk,3)),sum(Cells(chk2))];
results(:,[2,4],3)=results(:,[2,4],3)-results(:,[1,3],3); 

%large category FCC nb - NFCC nb - FCC cell - NFCC cell
chk=fcClNb(:,2)>21;
chk2=Order>21;
results(:,:,4)=[sum(chk),sum(chk2),sum(fcClNb(chk,3)),sum(Cells(chk2))];
results(:,[2,4],4)=results(:,[2,4],4)-results(:,[1,3],4); 

fcClNb=results(:,1:2,:);
fcClCell=results(:,3:4,:);


%%%%%%%%%%% Calculation 2 %%%%%%%%%%%%%%
%data analysis for the second calculation.
%first transform the cell vector in a matrix and get the unique lines
fcLpInf=unique(cat(1,fcLpInf{:}),'rows');

%FCC - count force chain clusters
    %get order of clusters beglonging to FC
results=zeros(1,4,4);

%small category FCC nb - NFCC nb - FCC cell - NFCC cell
chk=fcLpInf(:,2)==4;
results(:,:,1)=[sum(chk),sc.Clt4.nbCells,sum(fcLpInf(chk,3)),sc.Clt4.nbCells]; 
results(:,[2,4],1)=results(:,[2,4],1)-results(:,[1,3],1); 

%submedium category FCC nb - NFCC nb - FCC cell - NFCC cell
chk=fcLpInf(:,2)==6;
chk2=Order==6;
results(:,:,2)=[sum(chk),sum(chk2),sum(fcLpInf(chk,3)),sum(Cells(chk2))];
results(:,[2,4],2)=results(:,[2,4],2)-results(:,[1,3],2); 

%medium category FCC nb - NFCC nb - FCC cell - NFCC cell
chk=fcLpInf(:,2)>6 & fcLpInf(:,2)<21;
chk2=Order>6 & Order<21;
results(:,:,3)=[sum(chk),sum(chk2),sum(fcLpInf(chk,3)),sum(Cells(chk2))];
results(:,[2,4],3)=results(:,[2,4],3)-results(:,[1,3],3); 

%large category FCC nb - NFCC nb - FCC cell - NFCC cell
chk=fcLpInf(:,2)>21;
chk2=Order>21;
results(:,:,4)=[sum(chk),sum(chk2),sum(fcLpInf(chk,3)),sum(Cells(chk2))];
results(:,[2,4],4)=results(:,[2,4],4)-results(:,[1,3],4); 

fcClNb2=results(:,1:2,:);
fcClCell2=results(:,3:4,:);

%VTK Plot of individual grain evolution
if app.FCGrVtkCBox.Value
    fcGrainVTK(app,sc,cl4Id,nbCl,fc,gr,step,pathVTK(2))
end

end
function [bEv,fcEleO,angO,bID]=fcBending(app,fc,gr,fcEleO,angO,fcClTf)
%FCBENDING calculate the number of bending events
% A bending event takes place when the angle between 3 chained grains
% changes from one timestep to the other. To check the number of events
% first the forcechains will be separated in all possible elements of 3 and
% the angle they are forming will be measured. Next these angles will be
% checked with the previous timestep to see if bendling took place.

%Firts identify all possible 3 elements chains. Force chain branches are
%strings with the possible grain chains. So by analysing them all possible
%3 elements will be obtained
nbC=sum([fc.NbBranches]); %nb of branches
fcEleN=cell(nbC,1); %fc elements new - a cell array that will contain the data
k=1; %counter
%matrix to separate a chain into 3 elements chains. Get the max length of a
%force chain to set an upper limit and create a 3 column matrix containing
%IDs increasing by 1 that will allow the division of each branch into
%elements.
mtx=(1:(ceil(max([fc.Length])/10)*10))';
mtx=[mtx mtx+1 mtx+2]; 
for i=1:numel(fc)
    fci=fc(i);
    for j=1:numel(fci.Branches)
        %get the branch IDs
        bch=fci.Branches{j};
        if size(bch,1)>1; bch=bch';end %force column vector
        %prepare all possible combinations
        fcEleN{k}=bch(mtx(1:(numel(bch)-2),:));
        k=k+1;
    end
end
%Now get all the unique 3 element chains without changing the 3element
%order (ABC~=BCA for dist calculations)
fcEleN=cat(1,fcEleN{:});
[~,ia,~]=unique(sort(fcEleN,2),'rows');
fcEleN=fcEleN(ia,:); 
%calculate the angle through the acos of the cross product of the vectors
vect1=gr.Coord(fcEleN(:,2),:)-gr.Coord(fcEleN(:,1),:);
vect2=gr.Coord(fcEleN(:,2),:)-gr.Coord(fcEleN(:,3),:);
angN=acosd(dot(vect1,vect2,2)./(vecnorm(vect1').*vecnorm(vect2'))');

%Compare with previous values if fcEleO - fc elements old - vector was
%introduced in function call
if ~isempty(fcEleO)
    %identify trio of grains that survived inside a forcechain between
    %steps. location on Old - locO;location on New - locN;
    [locN,locO]=ismember(fcEleN,fcEleO,'rows');
    locO=locO(locO>0);
    %check if the angle between the surviving elements has changed of the
    %predetermined value app.FCBendAngleEF.Value. Sum the amount of values
    chck=abs(angN(locN)-angO(locO))>app.FCBendAngleEF.Value;
    bEv=sum(chck);
    %if cluster transformation was asked, return also the ID of grains
    %forming each 3-element chain that was bent in the "old" array 
    if fcClTf
        bID=fcEleO(locO,:);
        bID=bID(chck,:);
    else
        bID=[0 0 0];
    end
else
    bEv=0;bID=[0 0 0];
end
%save values from this T
fcEleO=fcEleN;
angO=angN;
end
function ctRes=fcBendingClTf(bID,sc,ctRes)
%FCBENDINGCLTF Check if bending events also correspond to cl transforation
% This function can only be executed if Cluster and Cluster transformation
% analysis were previously executed for this file. If the interval between 
% eachtimestep is changed this calculation may fail or return unreliable 
% results.
%  - Fist the cells of the loaded spaceCell (sc) file that contail the
% contacts that were bent (bID) will be identified. 
%  - Then a check will be performed on the
% results of a previous Cluster Transformation run (ctRes) to check if the
% identified cells had a change in cluster order.
%  - ctRes will be modified to

%get the contacts that were bent
uCtc=unique([bID(:,1:2);bID(:,2:3)],'rows');

%identify on sc the cells that share uCtc edges
gC=sum((sc.DelaunayT(:,:))>sc.NbG,2);
gC=find(gC==0);
spCell=sort(sc.DelaunayT(gC,:),2);
edges=cat(3,[spCell(:,1),spCell(:,2)],[spCell(:,1),spCell(:,3)],...
    [spCell(:,1),spCell(:,4)],[spCell(:,2),spCell(:,3)],...
    [spCell(:,2),spCell(:,4)],[spCell(:,3),spCell(:,4)]); %ngCx2x6
s=zeros(size(edges,1),1);
for i=1:6
    s=s+ismember(edges(:,:,i),uCtc,'rows');
end
cellID=find(s>0);

%For all Cl properties of ctRes, compare the OldID with identified cellID,
%only keeping the true values
ctRes.Cl4=ctRes.Cl4(ismember(ctRes.Cl4(:,2),cellID),:);
ctRes.Cl6=ctRes.Cl6(ismember(ctRes.Cl6(:,2),cellID),:);
ctRes.Cl8=ctRes.Cl8(ismember(ctRes.Cl8(:,3),cellID),:);

end
function fcGrainVTK(app,sc,cl4Id,L,fc,gr,step,pathVTK)
%FCGRAINVTK plot fc related to the specific grains asked
% This function will search for the grains with the IDs defined by the user
% and write in a vtk file only the force chains that contain it. For each
% step, for each grain, 4 files are created. First contain the force chain
% lines, second the fc grains, third clusters larger then 4 and last
% clusters 4.

grIds=str2double(regexp(app.FCGrVtkEF.Value,'\d*','match'));
%Create a Nx2 vectors with the first column the force chain ID and the
%second the grian belonging to the fc.
fcID=[repelem(1:numel(fc),[fc.NbGrains])',cat(1,fc.IDs)];
for i=1:numel(grIds)
    %Check if a force chain cointain the grain grIds(i) and get the Id of
    %this force chain. If it does not belong to a fc, continue
    idi=fcID(ismember(fcID(:,2),grIds(i)),1);
    if ~isempty(idi)
        %draw fc lines
        vtkwrite(fullfile(pathVTK,grIds(i)+"li"+step+".vtk"),'polydata',...
            "LINESB",gr.Coord(:,1),gr.Coord(:,2),gr.Coord(:,3),cat(1,fc(idi).Lines));
        %draw grains
        id=unique([fc(idi).IDs;grIds(i)]);
        vtkwrite(fullfile(pathVTK,grIds(i)+"gr"+step+".vtk"),'unstructured_grid',...
            gr.Coord(id,1),gr.Coord(id,2),gr.Coord(id,3),...
            'SCALARS','Radius',gr.Radius(id),...
            'PRECISION',10);
        %draw clusters : first get clusters belonging to the force chain.
        lpID=fc(idi).ClustersID;
        %start with cl>4
        cl=lpID(lpID(:,2)>4,1);
        if ~isempty(cl)
            drawClst=cat(1,sc.Loops(cl).Vertices);
            vtkwrite(fullfile(pathVTK,grIds(i)+"clL"+step+".vtk"),'polydata',"TRIANGLE",...
                gr.Coord(:,1),gr.Coord(:,2),gr.Coord(:,3),drawClst);
        else
            vtkwrite(fullfile(pathVTK,grIds(i)+"clL"+step+".vtk"),'polydata',...
                "NOLINES",gr.Coord(1,1),gr.Coord(1,2),gr.Coord(1,3));
        end
        %Cl4
        cl=lpID(lpID(:,2)==4,1);
        if ~isempty(cl)
            drawClst=sc.DelaunayT(cl4Id(cl-L),:); %Nx4 into 4Nx3
            drawClst=cat(1,[drawClst(:,1),drawClst(:,2),drawClst(:,3)],...
                [drawClst(:,1),drawClst(:,2),drawClst(:,4)],...
                [drawClst(:,1),drawClst(:,3),drawClst(:,4)],...
                [drawClst(:,2),drawClst(:,3),drawClst(:,4)]);
            vtkwrite(fullfile(pathVTK,grIds(i)+"clS"+step+".vtk"),'polydata',"TRIANGLE",...
                gr.Coord(:,1),gr.Coord(:,2),gr.Coord(:,3),drawClst);
        else
            vtkwrite(fullfile(pathVTK,grIds(i)+"clS"+step+".vtk"),'polydata',...
                "NOLINES",gr.Coord(1,1),gr.Coord(1,2),gr.Coord(1,3));
        end
        
    else
        vtkwrite(fullfile(pathVTK,grIds(i)+"li"+step+".vtk"),'polydata',...
            "NOLINES",gr.Coord(1,1),gr.Coord(1,2),gr.Coord(1,3));
        vtkwrite(fullfile(pathVTK,grIds(i)+"gr"+step+".vtk"),'unstructured_grid',...
            gr.Coord(grIds(i),1),gr.Coord(grIds(i),2),gr.Coord(grIds(i),3),...
            'SCALARS','Radius',gr.Radius(grIds(i)),...
            'PRECISION',10);
        vtkwrite(fullfile(pathVTK,grIds(i)+"clL"+step+".vtk"),'polydata',...
            "NOLINES",gr.Coord(1,1),gr.Coord(1,2),gr.Coord(1,3));
        vtkwrite(fullfile(pathVTK,grIds(i)+"clS"+step+".vtk"),'polydata',...
            "NOLINES",gr.Coord(1,1),gr.Coord(1,2),gr.Coord(1,3));
    end
end

end
%{
%%%%%%%%%%%%%%%%% Second calculation %%%%%%%%%%%%%%%%%%
% Count the the order of the cells located arround the force chains

    %get the cluster cells and add cluster order to it
clCe=[cat(2,sc.Loops.sCells)'...
    repelem(cat(1,sc.Loops.Order),cat(1,sc.Loops.nbCells),1)];
    %Fc contacts
fcCtc=cat(1,fc.Lines);
    %get gC from DT
gcDT=sc.DelaunayT(sc.GoodCells,:);
    %check matching edges between FC and GC thus obtaining the cells that
    %share a contact with the force chains
chk=zeros(numel(sc.GoodCells),1);
for i=1:6
    chk=chk+ismember(gcDT(:,cb(i,:)),fcCtc,'rows');
end
gcFC=sc.GoodCells(chk>1);
    %get the order of these cells
mb=ismember(clCe(:,1),gcFC);
fcO=clCe(mb,2);
    %count the number on each category knowing that the number of clusters
    %of order 4 is the values not matched between clCe(:,1) and gcFC.
%check cell per categories
clCl4=numel(gcFC)-sum(mb);     %small
clCl6=sum(fcO==6);             %submedium
clCl8=sum((fcO>6 & fcO<22));   %medium
clCl22=sum(fcO>20);            %large

%NFCC - non-force chain clusters are total minus FCC
clCl4(2)=sc.Clt4.nbCells-clCl4;
clCl6(2)=sum(C(O==6))-clCl6;            
clCl8(2)=sum(C(O>6 & O<21))-clCl8;	
clCl22(2)=sum(C(O>21))-clCl22;

fcClCell2=cat(3,clCl4,clCl6,clCl8,clCl22);
% Do the same but gettig all cells containing Fc grains
    %get good cells contianig at least one FC grains
gCgr=sc.GoodCells(sum(ismember(sc.DelaunayT(sc.GoodCells,:),cat(1,fc.IDs)),2)>0);
%get the order of these cells
mb=ismember(clCe(:,1),gCgr);
fcO=clCe(mb,2);
    %count the number on each category knowing that the number of clusters
    %of order 4 is the values not matched between clCe(:,1) and gcFC.
%check cell per categories
clCl4=numel(gCgr)-sum(mb);     %small
clCl6=sum(fcO==6);             %submedium
clCl8=sum((fcO>6 & fcO<22));   %medium
clCl22=sum(fcO>20);            %large

%NFCC - non-force chain clusters are total minus FCC
clCl4(2)=sc.Clt4.nbCells-clCl4;
clCl6(2)=sum(C(O==6))-clCl6;            
clCl8(2)=sum(C(O>6 & O<21))-clCl8;	
clCl22(2)=sum(C(O>21))-clCl22;

fcClCell3=cat(3,clCl4,clCl6,clCl8,clCl22);
%}
%}
