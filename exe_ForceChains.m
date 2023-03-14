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
    clData=cell(nbFiles,2);
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
        [clData{i,1:2}]=fcClusters(app,fc,gr,sc,step,pathVTKcl);
        if isempty(clData{i,2});return;end
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
C=app.PlotColors;

%Check title and legends option
if app.TitlesCB.Value;tit=1;else;tit=0;end
if app.LegendsCB.Value;leg=1;else;leg=0;end

%Create Figures and Axis
nb=9+bend;
if numel(pD)==1;nb=nb+1;end
f(nb)=figure;ax(nb)=axes(f(nb));hold(ax(nb),'on');
for i=1:(nb-1)
    f(i)=figure;ax(i)=axes(f(i));hold(ax(i),'on'); %#ok<LAXES>
end

%Time - Strain - mean(Length) - max (Length) - numel(fc) - ...
%numel(uniquefc) - numel(sigBranch) - nfc/HSTG - nfc/Tot
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
    if strcmpi(x,'p')
        xAx=res(:,2);%strain
    else
        xAx=res(:,1);%strain
    end
    %first graph : number of force chains TREES
    plotMark(app,ax(j),xAx,res(:,5),'Color',C(i,:),v{:});j=j+1;
    %third graph : mean force chain length
    plotMark(app,ax(j),xAx,res(:,3),'Color',C(i,:),v{:});j=j+1;
    %fourth graph : mean force chain length
    plotMark(app,ax(j),xAx,res(:,4),'Color',C(i,:),v{:});j=j+1;
    %fifth graph : FC/HS
    plotMark(app,ax(j),xAx,res(:,8),'Color',C(i,:),v{:});j=j+1;
    %sixth graph : FC/Tot
    plotMark(app,ax(j),xAx,res(:,9),'Color',C(i,:),v{:});j=j+1;
    %seventh graph : nb signlebranch/total
    plotMark(app,ax(j),xAx,100*res(:,7)./...
        res(:,5),'Color',C(i,:),v{:});j=j+1;
    %eight graph : 3p
    plotMark(app,ax(j),xAx,res(:,10),'Color',C(i,:),v{:});j=j+1;
    %ninith and possible tenth : Nb single branch and nb multible branch
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
    j=j+1;
    %Bending
    if bend
        plotMark(app,ax(j),xAx,res(:,end),'Color',C(i,:),v{:});
    end
end

if leg && numel(pD)>1
    for i=1:nb
        legend(ax(i),pD.FileName,'Location','best')
    end
end

i=1;
%first graph : number of force chains TREES
if tit;title(ax(i),'Evolution  of force chains');end
if strcmpi(x,'p');ax(i).XLim(1)=0;end
ylabel(ax(i),'Number')
xlabel(ax(i),xlab)
fnm="FC_Trees_Number";
saveas(f(i),fullfile(path,fnm+png));

i=i+1;
%third graph : mean force chain length 
if tit;title(ax(i),'Evolution of the mean length of MB force chains');end
if strcmpi(x,'p');ax(i).XLim(1)=0;end
ylabel(ax(i),'Number of grains')
xlabel(ax(i),xlab)
fnm="FC_Length_Mean_MB";
saveas(f(i),fullfile(path,fnm+png));

i=i+1;
%fourth graph : mean force chain length 
if tit;title(ax(i),'Evolution of the mean length of SB force chains');end
if strcmpi(x,'p');ax(i).XLim(1)=0;end
ylabel(ax(i),'Number of grains')
xlabel(ax(i),xlab)
fnm="FC_Length_Mean_SB";
saveas(f(i),fullfile(path,fnm+png));

i=i+1;
%fifth graph :FC/HS
if tit;title(ax(i),'Ratio of FC/HS grains');end
if strcmpi(x,'p');ax(i).XLim(1)=0;end
ylabel(ax(i),'Ratio of grains')
xlabel(ax(i),xlab)
fnm="FC_Grain_HS_Ratio";
saveas(f(i),fullfile(path,fnm+png));

i=i+1;
%sixth graph : FC/Tot
if tit;title(ax(i),'Ratio of FC/Tot grains');end
if strcmpi(x,'p');ax(i).XLim(1)=0;end
ylabel(ax(i),'Ratio of grains')
xlabel(ax(i),xlab)
fnm="FC_Grain_Tot_Ratio";
saveas(f(i),fullfile(path,fnm+png));

i=i+1;
%seventh graph : Ratio of Single branch trees over total
if tit;title(ax(i),'Ratio of Single branch trees');end
if strcmpi(x,'p');ax(i).XLim(1)=0;end
ylabel(ax(i),'Ratio of force chains')
xlabel(ax(i),xlab)
fnm="FC_Branches_Ratio";
saveas(f(i),fullfile(path,fnm+png));

i=i+1;
%eight graph : Ratio of Single branch trees over total
if tit;title(ax(i),'Number of 3p');end
if strcmpi(x,'p');ax(i).XLim(1)=0;end
ylabel(ax(i),'Number of 3p')
xlabel(ax(i),xlab)
fnm="FC_3p_Number";
saveas(f(i),fullfile(path,fnm+png));
%saveas(f(nb),path+fnm+fig);

i=i+1;
%ninith:
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

if bend
    i=i+1;
    %Bending results
    if tit;title(ax(i),'Evolution of Bending Event');end
    if strcmpi(x,'p');ax(i).XLim(1)=0;end
    ylabel(ax(i),'Number of bending events')
    xlabel(ax(i),xlab)
    fnm="FC_Bending";
    saveas(f(i),fullfile(path,fnm+png));
end

if nargin>2
    f(nb+1)=figure;ax(nb+1)=axes(f(nb+1));hold(ax(nb+1),'on');
    vals=varargin{1};
    scatter(ax(nb+1),1:numel(vals),vals,'x');
    if tit;title(ax(nb+1),'Ratio of time in force chain');end
    ylabel(ax(nb+1),'Ratio')
    xlabel(ax(nb+1),'Grain ID')
    fnm="FCRatioGrains";
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

%Create Figures and Axis
if size(pD,2)==1
    nplg=12;
    g(nplg)=figure;axG(nplg)=axes(g(nplg));hold(axG(nplg),'on');
    for i=1:(nplg-1)
        g(i)=figure;axG(i)=axes(g(i));hold(axG(i),'on'); %#ok<LAXES>
    end
end

% N1=app.N1EF.Value;
% N2=app.N2EF.Value;
% TD  = inflectionPoints(app.TrialData,app);
% pD.InfPts=TD.InfPts;
% fnm=MakePath(app,'FCCL')+"ForceChains-Cluster"+N1+"to"+N2+".mat";
% save(fnm,'pD','-v7.3');

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

if size(pD,2)==1
    C=app.PlotColors;
    if strcmpi(x,'p')
        xAx=pD(i).Results.Pressure;%strain
    else
        xAx=pD(i).Results.Strain;%strain
    end
    
    lpData=pD.Results.Data(:,2);
    lpData=cat(1,lpData{:});
    totCl=sum(lpData,3);
    % clusters of each type
    v={};
    for k=1:numel(pD(1).InfPts.q)
        if strcmpi(x,'p')
            v=[v,{'Pointx',pD(1).InfPts.p(k)}]; %#ok<AGROW>
        else
            v=[v,{'Pointx',pD(1).InfPts.ez(k)}]; %#ok<AGROW>
        end
    end
    j=1;
    %plot FC cl category
    plotMark(app,axG(j),xAx,lpData(:,1,1),'Color',C(1,:),v{:});
    plotMark(app,axG(j),xAx,lpData(:,1,2),'Color',C(2,:),v{:});
    plotMark(app,axG(j),xAx,lpData(:,1,3),'Color',C(3,:),v{:});
    plotMark(app,axG(j),xAx,lpData(:,1,4),'Color',C(4,:),v{:});

    j=j+1;
    %plot FC cl category - NO 4
    plotMark(app,axG(j),xAx,lpData(:,1,2),'Color',C(2,:),v{:});
    plotMark(app,axG(j),xAx,lpData(:,1,3),'Color',C(3,:),v{:});
    plotMark(app,axG(j),xAx,lpData(:,1,4),'Color',C(4,:),v{:});

    j=j+1;
    %plot NFC cl category
    plotMark(app,axG(j),xAx,lpData(:,2,1),'Color',C(1,:),v{:});
    plotMark(app,axG(j),xAx,lpData(:,2,2),'Color',C(2,:),v{:});
    plotMark(app,axG(j),xAx,lpData(:,2,3),'Color',C(3,:),v{:});
    plotMark(app,axG(j),xAx,lpData(:,2,4),'Color',C(4,:),v{:});

    j=j+1;
    %plot NFC cl category - NO 4
    plotMark(app,axG(j),xAx,lpData(:,2,2),'Color',C(2,:),v{:});
    plotMark(app,axG(j),xAx,lpData(:,2,3),'Color',C(3,:),v{:});
    plotMark(app,axG(j),xAx,lpData(:,2,4),'Color',C(4,:),v{:});

    j=j+1;
    %plot FC cl category PERCENTAGE
    plotMark(app,axG(j),xAx,lpData(:,1,1)./totCl(:,1),'Color',C(1,:),v{:});
    plotMark(app,axG(j),xAx,lpData(:,1,2)./totCl(:,1),'Color',C(2,:),v{:});
    plotMark(app,axG(j),xAx,lpData(:,1,3)./totCl(:,1),'Color',C(3,:),v{:});
    plotMark(app,axG(j),xAx,lpData(:,1,4)./totCl(:,1),'Color',C(4,:),v{:});

    j=j+1;
    %plot FC cl category PERCENTAGE - NO 4
    plotMark(app,axG(j),xAx,lpData(:,1,2)./totCl(:,1),'Color',C(2,:),v{:});
    plotMark(app,axG(j),xAx,lpData(:,1,3)./totCl(:,1),'Color',C(3,:),v{:});
    plotMark(app,axG(j),xAx,lpData(:,1,4)./totCl(:,1),'Color',C(4,:),v{:});

    j=j+1;
    %plot NFC cl category PERCENTAGE
    plotMark(app,axG(j),xAx,lpData(:,2,1)./totCl(:,2),'Color',C(1,:),v{:});
    plotMark(app,axG(j),xAx,lpData(:,2,2)./totCl(:,2),'Color',C(2,:),v{:});
    plotMark(app,axG(j),xAx,lpData(:,2,3)./totCl(:,2),'Color',C(3,:),v{:});
    plotMark(app,axG(j),xAx,lpData(:,2,4)./totCl(:,2),'Color',C(4,:),v{:});

    j=j+1;
    %plot NFC cl category PERCENTAGE - NO 4
    plotMark(app,axG(j),xAx,lpData(:,2,2)./totCl(:,2),'Color',C(2,:),v{:});
    plotMark(app,axG(j),xAx,lpData(:,2,3)./totCl(:,2),'Color',C(3,:),v{:});
    plotMark(app,axG(j),xAx,lpData(:,2,4)./totCl(:,2),'Color',C(4,:),v{:});
    
    %%%% - PERCELPLOT - %%%%
    lpData=pD.Results.Data(:,1);
    lpData=cat(1,lpData{:});
    totCl=sum(lpData,3);
    j=j+1;
    %plot FC cl category PERCENTAGE
    plotMark(app,axG(j),xAx,lpData(:,1,1)./totCl(:,1),'Color',C(1,:),v{:});
    plotMark(app,axG(j),xAx,lpData(:,1,2)./totCl(:,1),'Color',C(2,:),v{:});
    plotMark(app,axG(j),xAx,lpData(:,1,3)./totCl(:,1),'Color',C(3,:),v{:});
    plotMark(app,axG(j),xAx,lpData(:,1,4)./totCl(:,1),'Color',C(4,:),v{:});

    j=j+1;
    %plot FC cl category PERCENTAGE - NO 4
    plotMark(app,axG(j),xAx,lpData(:,1,2)./totCl(:,1),'Color',C(2,:),v{:});
    plotMark(app,axG(j),xAx,lpData(:,1,3)./totCl(:,1),'Color',C(3,:),v{:});
    plotMark(app,axG(j),xAx,lpData(:,1,4)./totCl(:,1),'Color',C(4,:),v{:});

    j=j+1;
    %plot NFC cl category PERCENTAGE
    plotMark(app,axG(j),xAx,lpData(:,2,1)./totCl(:,2),'Color',C(1,:),v{:});
    plotMark(app,axG(j),xAx,lpData(:,2,2)./totCl(:,2),'Color',C(2,:),v{:});
    plotMark(app,axG(j),xAx,lpData(:,2,3)./totCl(:,2),'Color',C(3,:),v{:});
    plotMark(app,axG(j),xAx,lpData(:,2,4)./totCl(:,2),'Color',C(4,:),v{:});

    j=j+1;
    %plot NFC cl category PERCENTAGE - NO 4
    plotMark(app,axG(j),xAx,lpData(:,2,2)./totCl(:,2),'Color',C(2,:),v{:});
    plotMark(app,axG(j),xAx,lpData(:,2,3)./totCl(:,2),'Color',C(3,:),v{:});
    plotMark(app,axG(j),xAx,lpData(:,2,4)./totCl(:,2),'Color',C(4,:),v{:});
    
    
    if leg
        for j=1:nplg
            if rem(j,2)==1
                legend(axG(j),"Order 4","Order 6","Order 8-20",...
                    "Order 22+",'location','eastoutside')
            else
                legend(axG(j),"Order 6","Order 8-20","Order 22+",...
                    'location','eastoutside')
            end
        end
    end
end


if size(pD,2)==1
    j=1;
    %Prepare axis to get the same value between NFC and FC equal graphs
    lim=axG(j).XLim;
    %NUMBER - FCC
    if tit;title(axG(j),'Distribution of FCC');end
    ylabel(axG(j),'Number of FCC')
    axG(j).YLim =[0,max(axG(j).YLim(2),axG(j+2).YLim(2))];
    if strcmpi(x,'e')
        xticks(axG(j),lim(1):0.05:lim(2))
    end
    xlabel(axG(j),xlab)
    fnm="FCC_Nb";
    saveas(g(j),fullfile(path,fnm+png));
    j=j+1;
    
    %NUMBER - FCC - NO4
    if tit;title(axG(j),'Distribution of FCC');end
    axG(j).YLim =[0,max(axG(j).YLim(2),axG(j+2).YLim(2))];
    if strcmpi(x,'e')
        xticks(axG(j),lim(1):0.05:lim(2))
    end
    ylabel(axG(j),'Number of FCC')
    xlabel(axG(j),xlab)
    fnm="FCC_Nb_No4";
    saveas(g(j),fullfile(path,fnm+png));
    j=j+1;
    
    %NUMBER - NFCC
    if tit;title(axG(j),'Distribution of NFCC');end
    axG(j).YLim =[0,max(axG(j).YLim(2),axG(j-2).YLim(2))];
    if strcmpi(x,'e')
        xticks(axG(j),lim(1):0.05:lim(2))
    end
    ylabel(axG(j),'Number of NFCC')
    xlabel(axG(j),xlab)
    fnm="NFCC_Nb";
    saveas(g(j),fullfile(path,fnm+png));
    j=j+1;
    
    %NUMBER - NFCC - NO4
    if tit;title(axG(j),'Distribution of NFCC');end
    axG(j).YLim =[0,max(axG(j).YLim(2),axG(j-2).YLim(2))];
    if strcmpi(x,'e')
        xticks(axG(j),lim(1):0.05:lim(2))
    end
    ylabel(axG(j),'Number of NFCC')
    xlabel(axG(j),xlab)
    fnm="NFCC_Nb_No4";
    saveas(g(j),fullfile(path,fnm+png));
    j=j+1;
    
    %PERCENTAGE - FCC
    if tit;title(axG(j),'Distribution of FCC');end
    axG(j).YLim =[0,1];
    if strcmpi(x,'e')
        xticks(axG(j),lim(1):0.05:lim(2))
    end
    ylabel(axG(j),'Ratio of FCC')
    xlabel(axG(j),xlab)
    fnm="FCC_Pct";
    saveas(g(j),fullfile(path,fnm+png));
    j=j+1;
    
    %PERCENTAGE - FCC - NO4
    if tit;title(axG(j),'Distribution of FCC');end
    axG(j).YLim =[0,max(axG(j).YLim(2),axG(j+2).YLim(2))];
    if strcmpi(x,'e')
        xticks(axG(j),lim(1):0.05:lim(2))
    end
    ylabel(axG(j),'Ratio of FCC')
    xlabel(axG(j),xlab)
    fnm="FCC_Pct_No4";
    saveas(g(j),fullfile(path,fnm+png));
    j=j+1;
    
    %PERCENTAGE - NFCC
    if tit;title(axG(j),'Distribution of NFCC');end
    axG(j).YLim =[0,1];
    if strcmpi(x,'e')
        xticks(axG(j),lim(1):0.05:lim(2))
    end
    ylabel(axG(j),'Ratio of NFCC')
    xlabel(axG(j),xlab)
    fnm="NFCC_Pct";
    saveas(g(j),fullfile(path,fnm+png));
    j=j+1;
    
    %PERCENTAGE - NFCC - NO4
    if tit;title(axG(j),'Distribution of NFCC');end
    axG(j).YLim =[0,max(axG(j).YLim(2),axG(j-2).YLim(2))];
    if strcmpi(x,'e')
        xticks(axG(j),lim(1):0.05:lim(2))
    end
    ylabel(axG(j),'Ratio of NFCC')
    xlabel(axG(j),xlab)
    fnm="NFCC_Pct_No4";
    saveas(g(j),fullfile(path,fnm+png));
    j=j+1;
    
    %CELL - FCC
    if tit;title(axG(j),'Distribution of FCC cells');end
    axG(j).YLim =[0,1];
    if strcmpi(x,'e')
        xticks(axG(j),lim(1):0.05:lim(2))
    end
    ylabel(axG(j),'Ratio of FCC cells')
    xlabel(axG(j),xlab)
    fnm="FCC_Cell";
    saveas(g(j),fullfile(path,fnm+png));
    j=j+1;
    
    %CELL - FCC - NO4
    if tit;title(axG(j),'Distribution of FCC cells');end
    axG(j).YLim =[0,max(axG(j).YLim(2),axG(j+2).YLim(2))];
    if strcmpi(x,'e')
        xticks(axG(j),lim(1):0.05:lim(2))
    end
    ylabel(axG(j),'Ratio of FCC')
    xlabel(axG(j),xlab)
    fnm="FCC_Cell_No4";
    saveas(g(j),fullfile(path,fnm+png));
    j=j+1;
    
    %CELL - NFCCE
    if tit;title(axG(j),'Distribution of NFCC cells');end
    axG(j).YLim =[0,1];
    if strcmpi(x,'e')
        xticks(axG(j),lim(1):0.05:lim(2))
    end
    ylabel(axG(j),'Ratio of NFCC cells')
    xlabel(axG(j),xlab)
    fnm="NFCC_Cell";
    saveas(g(j),fullfile(path,fnm+png));
    j=j+1;
    
    %CELL - NFCC - NO4
    if tit;title(axG(j),'Distribution of NFCC cells');end
    axG(j).YLim =[0,max(axG(j).YLim(2),axG(j-2).YLim(2))];
    if strcmpi(x,'e')
        xticks(axG(j),lim(1):0.05:lim(2))
    end
    ylabel(axG(j),'Ratio of NFCC cells')
    xlabel(axG(j),xlab)
    fnm="NFCC_Cell_No4_No4";
    saveas(g(j),fullfile(path,fnm+png));
end

%outside legend
if ~leg && numel(pD)==1
    %Order 4+ Legend
    p=copyobj(g(1),0);
    l=legend(p.CurrentAxes,"Order 4","Order 6","Order 8-20",...
            "Order 22+",'Orientation','horizontal');
    l.EdgeColor='none';
    set(p.CurrentAxes,'Visible','Off')
    set(p,'Position',(get(l,'Position').*[0, 0, 1, 1].*get(p,'Position')));
    set(l,'Position',[0,0,1,1]);
    set(p, 'Position', get(p,'Position') + [500, 400, 0, 0]);
    saveas(p,fullfile(path,"Legend_Order"+png));
    
    %order 6+Legend
    q=copyobj(g(2),0);
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
    if numel(pD)==1;delete(g);end
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
function [fcClCell,fcClNb]=fcClusters(app,fc,gr,sc,step,pathVTK)
%FCCLUSTER join force chain an cluster calculation data
% A cluster is considered as belonging in a force chain if at least one of
% its closed edges makes part of the force transmission.
%
% This function will use a spaceCellSystem object that contains the data
% from a previous cluster calculation and match with a forcechain object to
% identify the clusters that make part of the latter.

%General values
idGr=sort(cat(1,fc.IDs)); %identify grains forming the force chain
nbCl=numel(sc.Loops);nbCl4=sc.Clt4.nbCells;%nb of clusters
O=cat(1,sc.Loops.Order);C=cat(1,sc.Loops.nbCells);

%Create a list of all closed edges forming clusters, with the ID of the
%cluster in the third column and order as fourth
    %count the nb of closed edges per cluster
l=cell(nbCl,1);
[l{:}]=sc.Loops.ClEdges;
[nbCEdges,~]=cellfun(@size,l);
    %create a GrID1-GrID2-ClID-Order-nbCells vector for each closed contact
clGr=[(1:nbCl)',O,C]; %ClID-Order-nbCells
clCE=[cat(1,sc.Loops.ClEdges) repelem(clGr,nbCEdges,1)];
    %add cl4 to the analysis. All edges of Cl4 will be added and not only
    %the closed ones. But as it will be later compared to closed edges
    %belonging to force chains it will be automatically filtered
    %identify clusters4 edges
if isempty(sc.GoodCells);sc.GoodCells = goodCell(sc);end
cl4Id=sc.GoodCells(~ismember(sc.GoodCells,cat(2,sc.Loops.sCells)')); %id Cl4 cell location
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
clCE=clCE(ismember(clCE(:,1),idGr) & ismember(clCE(:,2),idGr),:);

%for each force chain, compare the force chain contacts with clCE matrix
%to identify the ID of the clusters forming it.
fcChk=ones(numel(fc),1);
for i=1:numel(fc)
    lp=ismember(clCE(:,1:2),fc(i).Lines,'rows');
    if isempty(lp)
        %if there are no clusters in the force chains, the fc must be
        %located around the wall and surrounded by 'fakegrains'. Thus it
        %should be removed from the calculation
        fcChk(i)=1;continue
    end
    lp=unique(clCE(lp,3:end),'rows');  %ClusterID - ClusterOrder - NbCells
    
    %max - min - mean order values
    fc(i).ClustersVals=[max(lp(:,2)),min(lp(:,2)),mean(lp(:,2))]; 
    fc(i).ClustersID=lp;
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
	%check nb per category
nbCl4=sum(fcClNb(:,2)==4);            %small
nbCl6=sum(fcClNb(:,2)==6);            %submedium
nbCl8=sum(fcClNb(:,2)>6 & fcClNb(:,2)<21);	%medium
nbCl22=sum(fcClNb(:,2)>21);           %large

%NFCC - non-force chain clusters are total minus FCC
nbCl4(2)=sc.Clt4.nbCells-nbCl4;
nbCl6(2)=sum(O==6)-nbCl6;            
nbCl8(2)=sum(O>6 & O<21)-nbCl8;	
nbCl22(2)=sum(O>21)-nbCl22;

%check cell per categories
clCl4=sum(fcClNb(fcClNb(:,2)==4,3));            %small
clCl6=sum(fcClNb(fcClNb(:,2)==6,3));            %submedium
clCl8=sum(fcClNb(fcClNb(:,2)>6 & fcClNb(:,2)<21,3));	%medium
clCl22=sum(fcClNb(fcClNb(:,2)>21,3));           %large

%NFCC - non-force chain clusters are total minus FCC
clCl4(2)=sc.Clt4.nbCells-clCl4;
clCl6(2)=sum(C(O==6))-clCl6;            
clCl8(2)=sum(C(O>6 & O<21))-clCl8;	
clCl22(2)=sum(C(O>21))-clCl22;

%concatenate them in third dimention so next step can go as next line
fcClNb=cat(3,nbCl4,nbCl6,nbCl8,nbCl22);
fcClCell=cat(3,clCl4,clCl6,clCl8,clCl22);

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
function [FcClData,clCat]=fcClustersOld(app,fc,gr,sc,step,pathVTK)
%FCCLUSTER join force chain an cluster calculation data
% A cluster is considered as belonging in a force chain if at least one of
% its closed edges makes part of the force transmission.
%
% This function will use a spaceCellSystem object that contains the data
% from a previous cluster calculation and match with a forcechain object to
% identify the clusters that make part of the latter.

%First part - find out for each grain the clusters it is attached to. 
idGr=sort(cat(1,fc.IDs));
sGr = findClusterGrain('ForceChain',sc,idGr,gr);


%for each force chain - get all the singleGrain objects that make part of
%the chain. Check then for Clsuters that have a common closed edge with the
%force chain
FcClData=zeros(numel(fc),5);
fcChk=ones(numel(fc),1);
for i=1:numel(fc)
    %get the position of the grains forming the fc(i) in the singleGrain
    %objects array created, as it only contains grains belonging to all FCs.
    grID=find(ismember(idGr,fc(i).IDs));
    clt=cat(1,sGr(grID).ClusterID);
    %check for repeating clusters IDs in the obtained grains
    [uClt,~,cnt]=unique(clt(:,1));	%unique value and their oder
    h=accumarray(cnt,ones(numel(cnt),1)); %get nb of repeated values
    uClt=uClt(h>1);  %get the ID of the clusters repeated more than once
    
    %update each grain's clusters only with the clusters inside the force
    %chain. As each grain can only belong to a single force chain, the
    %value save does not need to be done at the end.
    gr0=0; %check for gr with no clusters shared with contact
    for j=1:numel(grID)
        chk=ismember(sGr(grID(j)).ClusterID(:,1),uClt);
        if sum(chk)==0 || isempty(chk)
            %if chk is 0 means that the grain analysed does not share a
            %cluster with any of the other grains. The reason for this is
            %that the cells that would join them are in connection with
            %'fakegrains' and so removed from the cluster calculation. This
            %force chain then should be removed from the calculation
            gr0=1;
            break
        end
        sGr(grID(j)).ClusterID=sGr(grID(j)).ClusterID(chk,:);  %ClusterID - ClusterOrder
        sGr(grID(j)).ClusterVal=[max(sGr(grID(j)).ClusterID(:,2)),...
            min(sGr(grID(j)).ClusterID(:,2)),...
            mean(sGr(grID(j)).ClusterID(:,2))];
    end
    if gr0
        %get the locaiton of the force chain that must be removed
        fcChk(i)=0;
        continue
    end
    fc(i).ClustersVals=cat(1,sGr(grID).ClusterVal); %max - min - mean values
    fc(i).ClustersID=unique(cat(1,sGr(grID).ClusterID),'rows'); %ClusterID - ClusterOrder
    FcClData(i,:)=[fc(i).Length max(fc(i).ClustersVals(:,1)),...
        min(fc(i).ClustersVals(:,2)),mean(fc(i).ClustersVals(:,3)),...
        mean([sGr(grID).PrincipalStress])];
end
fcChk=logical(fcChk);%turn into logical

%remove data from the bad force chains
FcClData=FcClData(fcChk,:); %remove 0 lines
fc=fc(fcChk);               %delete bad fc
idNew=sort(cat(1,fc.IDs));  %get the correct grain IDs
clV=cat(1,sGr(ismember(idGr,idNew)).ClusterVal); %per grain cluster values

%vtk save
vtkwrite(pathVTK(1)+"connectionsFC"+step+".vtk",'polydata',"LINESB",...
    gr.Coord(:,1),gr.Coord(:,2),gr.Coord(:,3),cat(1,fc.Lines));
vtkwrite(pathVTK(1)+"grainsFC"+step+".vtk",'unstructured_grid',...
    gr.Coord(idNew,1),gr.Coord(idNew,2),gr.Coord(idNew,3),...
    'SCALARS','Radius',gr.Radius(idNew),...
    'SCALARS','MaxCl',clV(:,1),...
    'SCALARS','MinCl',clV(:,2),...
    'SCALARS','MeanCl',clV(:,3),...
    'PRECISION',10);

%Cluster categories
    %create an array containing ClID - ClOrder
nbCl=numel(sc.Loops);nbCl4=sc.Clt4.nbCells;
clGr=[(1:nbCl4+nbCl);
    cat(2,sc.Loops.Order),ones(1,nbCl4)*4]';
    %get unique Ids of Clusters in fc
fcClIds=cat(1,fc.ClustersID);
fcClIds=unique(fcClIds(:,1));
    %Transform these Ids into logical position
lgc(length(clGr))=false;
lgc(fcClIds)=true;
    %create a matrix where the first column is the cl in FC  and the other
    %NFC ones.
allCl=[clGr(:,2) clGr(:,2)];
allCl(~lgc,1)=0;
allCl(lgc,2)=0;

	%check categories
cl4=sum(allCl==4,1);            %small
cl6=sum(allCl==6,1);            %submedium
cl8=sum(allCl>6 & allCl<21,1);	%medium
cl22=sum(allCl>21,1);           %large
    %concatenate them in third dimention so next step can go as next line
clCat=cat(3,cl4,cl6,cl8,cl22);

if app.FCGrVtkCBox.Value
    fcGrainVTK(app,sc,cl4Id,nbCl,fc,gr,step,pathVTK(2))
end

end
%}
