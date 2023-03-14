function exe_StrainTensor(app,type,PD)
%STRAINTENSOR Calculates the evolution of the internal strain
%	This function has three modes of execution Global, PerCell or Both. All
%	of them using the spaceCellSystem class.
%   Global will return a graph with the evolution of the strain during the
%   expertiment.
%   PerCell will create a VTK file where each grain will have the value of
%   the strain on its neiborhood.
%   Both will do both calculations.

%Chose the time of calculation following TYPE
switch upper(type)
    case 'LOAD'
        strainLoad(app);return;
    case 'GLOBAL'
        calcType=1;
    case 'PERCELL'
        calcType=2;
    otherwise
        calcType=[1,2];
end

%Check prefix
if isempty(PD)
    prefix="total";mode='NORMAL';
elseif app.SubdivisionButton.Value
    prefix="subd";mode='SUBD';
    %read partial data variables
    lMax=PD.SubLin;
    cMax=PD.SubCol;
    nbSub=lMax*cMax;
else
    prefix="partial";mode='NORMAL';
end

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
    stepArray=(N1+interval:interval:N2)';
    if stepArray(end)~=N2; stepArray=[stepArray;N2];end
end
nbFiles=numel(stepArray)-1;

%Check if 3d or 2D
if app.Bool3D;D=3;else; D=2;end

%Prepare GLOBAL values
if ismember(1,calcType) 
    if ~app.SubdivisionButton.Value
        Results=zeros(nbFiles,D+3,3); %Step Ey Ez Ed Ev
        V_old='';
    else
        Results=zeros(nbFiles,D+3,nbSub); %Step Ex Ey Ez Ed Ev
    end
end
%Prepare PERCELL values
if ismember(2,calcType) 
    %creating a folder for the outputfiles (vtk)
    pathVTK=MakePath(app,'Strain')+"ShearBandsVTK"+N1+"to"+N2+"int"+...
        interval+"/";
    if exist(pathVTK,'dir')==0;mkdir(pathVTK);end
    % values calculated are incremental, thus will be cumulated in the
    % following vectors,
    pgEd=zeros(1,1,app.NbGrains);
    W2=zeros(1,1,app.NbGrains);
    %W2 and cluster calculation
    if app.StTW2ClCB.Value
        nbW2Gr=zeros(nbFiles,1);
        clW2Cat=zeros(nbFiles,4,2);
    end
    %Prepare per grain inertia calculation
    itvArray=stepArray(2:end)-stepArray(1:end-1);
    I=cell(nbFiles,1);
end

%Load first object containing grains data
grOld=grains('STENSOR',stepArray(1),'',app);
if isempty(grOld.Coord)
    warndlg(['gr.Coord is empty on step ' N1]);return;
end
stepArray=stepArray(2:end);

%create path
pathSc=MakePath(app,'SCF');

%Turn on calculation pannel
app=CalcPanel(app,'',nbFiles,'Starting calculation','on');
%Start Calculation
for i=1:(nbFiles)
    if getappdata(app.CalculatingPanel,'Canceling')
        CalcPanel(app,'','','','off');
        warndlg('Calculation was canceled');return
    end
    step=stepArray(i);
    
    if isa(app.CalculatingPanel,'double')
        app=CalcPanel(app,i,nbFiles,"NEWLINE");
    else
        app=CalcPanel(app,i,nbFiles,step);
    end
    %Load second object containing grains data
    [grNew,PD]=grains('STENSOR',step,PD,app);
    if isempty(grNew.Coord)
        CalcPanel(app,'','','','off');
        warndlg(['Gr.Coord is empty on step ' step]);return;
    end
    
    %check if spaceCell file exist from a previous execution. If not launch
    %the calculation
    fl=[pathSc char("spaceCellsfile"+step+"int"+interval+".mat")];
    if isfile(fl)
        %fprintf('Space cell %d file found and loaded \n',step)
        sc=load(fl).sc;sv=0;
    else
        sc = spaceCellSystem(type,step,grNew,grOld,app,PD);sv=1;
    end
    
    %GLOBAL CALCULATION
    if ismember(1,calcType)
        %check if the tensor exists or need to be calculated
        sc=gStrainTensor(sc,app,PD);
        %calculate total volume and save values into matrix
        if ~app.SubdivisionButton.Value
            %divide the strain value by hte volume of the previus step.
            %However for the first interval it is not calculated so divide
            %it by the volume of the actual step.
            if isempty(V_old);V_old=sum(sc.CellVol);end 
            stG(:,:,1)=(sc.GStTensor)/V_old;
            stG(:,:,2)=sum(sc.CellStn,3)/V_old;
            Results(i,1,3)=step;
            tri=2;
            V_old=sum(sc.CellVol);
        else
            [~,V_sub]=getVolume(app.TrialData,app,(step),PD); %v of previous step
            stG=(sc.GStTensor)./permute(V_sub,[3,2,1]);
            tri=size(Results,3);
        end
        %calculate incremental deviatoric and principal strain
        if D==3
            iEv=(stG(1,1,1:tri)+stG(2,2,1:tri)+stG(3,3,1:tri));
            iEd=sqrt(( (stG(2,2,1:tri)-stG(1,1,1:tri)).^2+...	%(sig1-sig2)^2 /2
                (stG(3,3,1:tri)-stG(1,1,1:tri)).^2+...       %(sig2-sig3)^2 /2
                (stG(1,1,1:tri)-stG(3,3,1:tri)).^2+...       %(sig3-sig1)^2 /2
                6*stG(1,2,1:tri).^2+...                    % 3*(sig12)^2
                6*stG(1,3,1:tri).^2+...                    % 3*(sig13)^2pcq
                6*stG(2,3,1:tri).^2 )/2);                  % 3*(sig23)^2
            Results(i,:,1:tri)=[ones(1,1,tri)*step...
                -[stG(1,1,1:tri) stG(2,2,1:tri) stG(3,3,1:tri) iEv ] iEd ]; %- for negative compression, GC standart
        else
            iEv=(stG(1,1,1:tri)+stG(2,2,1:tri));
            iEd=sqrt(( (stG(1,1,:)-stG(2,2,:)).^2+...     %(sig1-sig2)^2
                6*stG(1,2,:).*stG(1,2,:) )/2);          % 2*(sig12)^2
            Results(i,:,1:tri)=[ones(1,1,tri)*step...
                -[stG(1,1,1:tri) stG(2,2,1:tri) iEv ] iEd ]; %- for negative compression
        end
    end
    
    %PERCELL CALCULATION
    if ismember(2,calcType)
        %check if the tensor exists or need to be calculated
        if isempty(sc.PStTensor)
            sc=pgStrainTensor(sc);sv=1;
        end
        %Second order work. Per element ultiply the incremental strain and
        %stress. Then trainsform it from a 3x3xNbg to Nbgx1 by sum and
        %permtute
        incW2=(grNew.PGStressTensor-grOld.PGStressTensor).*sc.PStTensor;
        incW2=permute(sum(incW2,1:2),[3,1,2]);%./(grNew.Radius.^3*4/3*pi()); remove the volume division for now
        W2=W2+incW2;
        %calculate the cumulated strain
        iSTG=sc.PStTensor;  %incr grain strain (BAGI96)
        iSTC=sc.CellStn;   %incr cell strain (Interpolation)
        %Calculate the incremental deviatoric values
        if app.Bool3D % check existence of x direction piston
            iEd=sqrt(((iSTG(1,1,:)-iSTG(2,2,:)).^2+...
                (iSTG(2,2,:)-iSTG(3,3,:)).^2+ ...
                (iSTG(3,3,:)-iSTG(1,1,:)).^2+...
                6*iSTG(1,2,:).^2 + ...
                6*iSTG(1,3,:).^2 + ...
                6*iSTG(2,3,:).^2) /2);
            iEdC=sqrt(((iSTC(1,1,:)-iSTC(2,2,:)).^2+...
                (iSTC(2,2,:)-iSTC(3,3,:)).^2+ ...
                (iSTC(3,3,:)-iSTC(1,1,:)).^2+...
                6*iSTC(1,2,:).^2 + ...
                6*iSTC(1,3,:).^2 + ...
                6*iSTC(2,3,:).^2) /2);
        else
            iEd=sqrt(1/2*((iSTG(1,1,:)-iSTG(2,2,:)).^2+...
                6*iSTG(1,2,:).*iSTG(1,2,:) ));
            iEdC=sqrt(1/2*((iSTC(1,1,:)-iSTC(2,2,:)).^2+...
                6*iSTC(1,2,:).*iSTC(1,2,:) ));
        end
        %calculate the deviatoric value
        pgEd=pgEd+iEd;
        %Save the pergrain results in a matrix that will be writen in a
        %text file.
        
        %divide by cell volume
        iEdC=pagemtimes(iEdC,1./sc.CellVol); 

        cord=grNew.Coord;
        r=grNew.Radius;
        if isempty(PD)
            %ind=find(~isnan(sc.GrainVolume));
            ind=1:sc.NbG;
        else
            ind=PD.GrainsRectangle;
        end
        
        
        
        %Calculate the intertial number to verify if we still are in the
        %quasi-static regime
        P=permute(abs(grNew.PGStressTensor(1,1,:)+...
            grNew.PGStressTensor(2,2,:)+...
            grNew.PGStressTensor(3,3,:))/3,[3,2,1]);
        V=4/3*pi()*r.^3;
        I{i}=permute(iEd,[3,2,1])/(itvArray(i)*app.TimeStep)...
            .*sqrt(2600*V.^2./(P.*(2*r)));
        %save as vtk file
        if ~app.SubdivisionButton.Value
            %%%%% GRAIN VTK FILE %%%%
                %remove positive values and calculate a 90% negative
                %percentile of the values to be ploted. Thus a scale will
                %be created.
            pW2=W2;pW2(pW2>0)=0;
            k=sort(pW2);k(k==0)=[];
            k=k(max(floor(0.01*numel(k)),1));
            pW2=pW2/abs(k);
            pW2(pW2<-1)=-1;
            
            incpW2=incW2;incpW2(incpW2>0)=0;
            k=sort(incpW2);k(k==0)=[];
            k=k(max(floor(0.01*numel(k)),1));
            incpW2=incpW2/abs(k);
            incpW2(incpW2<-1)=-1;
            
            newFile=pathVTK+prefix+"ShearBand"+step+".vtk";
            vtkwrite(newFile,'unstructured_grid',...
                cord(ind,1),cord(ind,2),cord(ind,3),...
                'SCALARS','Radius',r(ind),...
                'SCALARS','Dev_Str',pgEd(ind),...
                'SCALARS','Dev_Str_Rat',pgEd(ind)/max(pgEd(ind)),...
                'SCALARS','Dev_Str_Inc',iEd(ind),...
                'SCALARS','W2',pW2(ind),...
                'SCALARS','W2_Inc',incpW2(ind),...
                'PRECISION',10,'BINARY');
            %'SCALARS','DevStrainRatio',pcEv(ind)/max(pcEv(ind)),...
            
            %%%%% CELL VTK FILE %%%%
            
            %Inc dev strain of cell with too low volume will be turned to 0
            cv=cat(1,sc.CellVol);
            cv=(cv<mean(cv)^2/max(cv));
            iEdC(cv)=0;
            k=sort(iEdC);
            k=k(floor(0.99*numel(k)));
            iEdC=iEdC/abs(k);
            iEdC(iEdC>1)=1;
                %calculate per cell VR
            [vr,~,~]=perCellVoidRatio(sc,grNew);
                %prepare to create vtk files
            newFile=pathVTK+prefix+"ShearBandCell"+step+".vtk";
            cGr=sc.CellGrn; %grain ditribution per cell
            cmb=nchoosek(1:(D+1),D);
            tri=cat(1,cGr(:,cmb(1,:),:),cGr(:,cmb(2,:),:),...
                cGr(:,cmb(3,:),:),cGr(:,cmb(4,:),:));
            tri=reshape(permute(tri,[2,1,3]),3,[],1)';
            iEdC=repelem(permute(iEdC,[3,2,1]),D+1);
            vr=repelem(vr,D+1);
            vtkwrite(newFile,'polydata',"TRIANGLE",...
                cord(:,1),cord(:,2),cord(:,3),tri,...
                'facescalar',2,"Dev_Str_Inc",iEdC,...
                "Void_Ratio",vr);
        else
            for l=1:lMax
                for tri=1:cMax
                    j=tri+cMax*(l-1);
                    if ismember(2,calcType)
                        newFile=pathVTK+prefix+"ShearBandl"+l+"c"+tri+"-"+step+".vtk";
                        vtkwrite(newFile,'unstructured_grid',...
                            cord(PD.SubGrains{j},1),cord(PD.SubGrains{j},2),cord(PD.SubGrains{j},3),...
                            'SCALARS','Radius',r(PD.SubGrains{j}),...
                            'SCALARS','Dev_Strain',Ed(PD.SubGrains{j}),...
                            'SCALARS','Dev_Strain_Inc',iEd(PD.SubGrains{j}),...
                            'SCALARS','Dev_Strain_Ratio',Ed(PD.SubGrains{j})/mean(Ed(PD.SubGrains{j})),...
                            'SCALARS','SecOrderWork',incW2(PD.SubGrains{j}),...
                            'PRECISION',10,'BINARY');
                    end
                end
            end
        end
        %Extra calculation : Cluster and second order work realtion
        if app.StTW2ClCB.Value
            [nbGr,clW2Gr]=w2Cluster(app,step,incW2);
            if isempty(nbGr);return;end
            nbW2Gr(i)=nbGr;
            clW2Cat(i,:,:)=clW2Gr;
        end
    end
    %Save SC in a file for a faster re execution if needed
    if sv==1 || sv==0
        sc=purge(sc,'Strain');
        save(fl,'sc','-v7.3');
    end
    %grNew becomes the next grOld
    grOld=grNew;
end
CalcPanel(app,i+1,nbFiles,'','off');

if ismember(2,calcType) && app.SubdivisionButton.Value
    %write down the subdivision selection as a vtk for easier analysis
    v=PD.SubVertices;
    pt=zeros(size(v,3)*2,3);    %contain the points
    cnct=zeros(size(v,3),2);  %contain the connection scheeme
    for i=1:size(v,3)
        if i==size(v,3)
            y=ismembertol(v(:,1,i-1),v(:,1,i));
            z=ismembertol(v(:,2,i-1),v(:,2,i));
            uV=unique(v(~(y & z),:,i),'rows');
        else
            y=ismembertol(v(:,1,i+1),v(:,1,i));
            z=ismembertol(v(:,2,i+1),v(:,2,i));
            uV=v((y & z),:,i);
        end
        pt((1+2*(i-1)):2*(i),:)=[+0.075,uV(1,:);
            +0.075,uV(2,:)];
        cnct(i,:) =[1 2]+(i-1)*2;
    end
    fnm=pathVTK+"subdivisionsLines.vtk";
    vtkwrite(fnm,'polydata',"LINESB",pt(:,1),pt(:,2),pt(:,3),cnct);
end

%Consolidation strain value
consoStrain = extStrains(app.TrialData,app.ConsoStep,N1,app);

%Write down Global results
if ismember(1,calcType)
    %Values calculated were incremental, need to sum them up
    Results(:,2:end,:)=cumsum(Results(:,2:end,:),1);
    
    %Calculate external strain
    if ~app.SubdivisionButton.Value
        Results(:,2:end,3)=extStrains(app.TrialData,Results(:,1,2),N1,app,'dev');
    end
    %Transform timestep into time
    Results(:,1,:)=Results(:,1,:)*app.TimeStep;
    
    %write down values
    
    res.Strain=Results;
    if exist('I','var')
        res.Inertia=I;
    else
        res.Inertia='';
    end 
    res.SimType=app.SimType;

    fnm=fullfile(MakePath(app,'STRAIN'),"Strain"+N1+"to"+N2+"int"...
        +interval+".mat");
    if ~app.SubdivisionButton.Value
        pD = plotData("Normal",res,app,prefix,consoStrain(end));
    else
        pD = plotData("Normal",res,app,prefix,consoStrain(end),lMax,cMax);
    end
    save(fnm,'pD','-v7.3');    %save
    %start plot
    strainPloter(app,mode,pD);
end
%create vtk state variables file
if ismember(2,calcType)
    fnm=pathVTK+"vtkVariables.txt";
    stress = extStress(app.TrialData,stepArray,app);
    strain = extStrains(app.TrialData,stepArray,N1,app,'allCalc');
    StepEzEvQP=[stepArray,strain(:,[D,D+1]),stress(:,end-1:end)];
    vtkLog(app,'FCbase',fnm,StepEzEvQP)
end
if app.StTW2ClCB.Value
    %Cluster and Second order work results
    w2Res.Nb=nbW2Gr/app.NbGrains;
    w2Res.Cluster=clW2Cat;
    strain= extStrains(app.TrialData,stepArray,N1,app);
    w2Res.Strain=strain(:,D);
    w2Res.InfPts=app.TrialData.InfPts;
    %create plotData object
    pD = plotData("Normal",w2Res,app,'',consoStrain(end));
    %save pD object
    fl=MakePath(app,'LOOPW2')+"Striain-ClW2"+N1+"to"+N2+"int"+interval+".mat";
    save(fl,'pD','-v7.3');
    %go to plot
    w2ClustPlotter(app,pD);
end
end
%execution function
function [nbW2Gr,clW2Gr]=w2Cluster(app,step,W2)
%W2CLUSTER Calculates the clusters near negative W2 values

clW2Gr=zeros(1,4,2);
%Load Cluster data
scFnm=[app.SavePath '/MatlabResultsFile/SpaceCellFiles/'...
    char("LoopsSpaceCellsfile"+step+".mat")];
try sc=load(scFnm).sc;
catch
    fprintf(['No Loops spacecell files, '...
        'please run loops calculation beforehand.'])
    nbW2Gr='';clW2Gr='';return
end

%create a 3 column 2D matrix containing the ID of the cluster, the order
%and each grain. So for every grain of each Cluster the first two elements
%of the line is the same.  %ClusterID - ClusterOrder - grID
Cl=numel(sc.Loops);
clGr=repelem([(1:Cl)',cat(1,sc.Loops.Order)],cat(1,sc.Loops.Size),1);
clGr=[clGr cat(1,sc.Loops.Grains)];

%do the same for the Cluster4 category %ClusterID - grID
cl4Id=sc.GoodCells(~ismember(sc.GoodCells,cat(2,sc.Loops.sCells)')); %id Cl4 cell location
clGr4=repelem((1+Cl):(numel(cl4Id)+Cl),4)'; %create ClusterID starting from L+1
clGr4=[clGr4 reshape(sc.DelaunayT(cl4Id,:)',[numel(cl4Id)*4,1]) ];

%Get both matrixes in one vector
clGr = [clGr;clGr4(:,1) ones(length(clGr4),1)*4 clGr4(:,2)];  %ClusterID - ClusterOrder - grID

%Get ID of negative W2 grains
id=find(W2<0);
nbW2Gr=numel(id); %nb of negative W2 grains
%get cl IDs of grains in contact to negative W2 grains
pos=ismember(clGr(:,3),id);
clW2id=unique(clGr(pos,1:2),'rows'); %ClusterID - ClusterOrder
cl=clW2id(:,2);
    %count categories
clW2Gr(1,:,1)=sum([cl==4,cl==6,(cl>6 & cl<21),cl>21],1);
%nb of non W2 clusters = total - w2clusters
cl=cat(1,sc.Loops.Order);
clW2Gr(1,:,2)=[sc.Loop4,sum([cl==6,cl>6 & cl<21,cl>21],1)]...
    -clW2Gr(1,:,1);
end
%load functions
function strainLoad(app)
%STRAINLOAD Loader for the strain function

if app.StTW2ClCB.Value
    bsPth=MakePath(app,'LOOPW2','check');
    if exist(bsPth,'dir')==0
        fprintf('Cluster W2 usual path not found \n')
        bsPth='';
    end    
    [fnm,Path,~] = uigetfile( ...
        {'*.mat','Matlab files(*.mat)';'*.*',  'All Files (*.*)'},...
        'Select file for loading','Multiselect','on',bsPth);
    if Path==0
        warndlg('No files were chosen, Try again.')
    else
        pD=load([Path fnm]).pD;
        w2ClustPlotter(app,pD)
    end
    return
end

%Load file
[fnm,fPath]=MatLoader('STRAIN',app);
if fPath==0;return;end
if iscell(fnm)
    pD=plotData.empty(0,1);
    mode='MULTI';
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
    if size(pD.Results,3)<4
        mode='NORMAL';
    else
        mode='SUBD';
    end
end

%Plot
if numel(pD)==1
    strainPloter(app,mode,pD)
else
    fprintf('Executing many\n')
    for i=1:numel(pD)
        mode='NORMAL';
        path=fullfile(MakePath(app,'STRAIN'),pD(i).FileName);
        if exist(path,'dir')==0;mkdir(path);end
        strainPloter(app,mode,pD(i),path)
    end
    fprintf('Done\n')
end

end
%plot functions
function strainPloter(app,mode,pD,path)
%STRAINPLOTTER Plot the strain values

%Subdivision
if isequal(upper(mode),'SUBD')
    subV=[pD.SubL,pD.SubC];
    labels=strings(1,subV(1)*subV(2));
    for l=1:subV(1)
        for c=1:subV(2)
            labels((l-1)*subV(2)+c)="Sub l"+l+" c"+c;
        end
    end
    SubPlotter(app,'STRAIN',pD.Results.Strain,labels,pD.ConsoTime,subV);
    return;
end

%Check title and legends option
if app.TitlesCB.Value;tit=1;else;tit=0;end
if app.LegendsCB.Value;leg=1;else;leg=0;end

if ~exist('path','var');path=MakePath(app,'STRAIN');end
%testplotshit(pD)
fi=pD.Results.Strain;
%CHECK DIMENSION
if pD(1).Bool3D;D=3;else; D=2;end
try lw=app.PlotWidthEF.Value;
catch
    lw=1.5;
end
%Create Figures and Axis
if D==3
    nb=5;
else
    nb=4;
end
if isequal(upper(mode),'MULTI');nb=nb+1;end

f(nb)=figure;ax(nb)=axes(f(nb));hold(ax(nb),'on');
for i=1:(nb-1)
    f(i)=figure;ax(i)=axes(f(i));hold(ax(i),'on'); %#ok<LAXES>
end

switch upper(mode)
    case 'NORMAL'
        %LEGENDS
        str=["BAGI96","LinInterp","ExternalCalculation"];
        lty=["-^","-s","-o"];
        for i=1:size(fi,3)
            if D==3
                plot(ax(5),fi(:,1,i),fi(:,2,i),lty(i),'LineWidth',lw)  %Ex=f(t)
            end
            plot(ax(1),fi(:,1,i),fi(:,D,i),lty(i),'LineWidth',lw)      %Ey=f(t)
            plot(ax(2),fi(:,1,i),fi(:,D+1,i),lty(i),'LineWidth',lw)    %Ez=f(t)
            plot(ax(3),fi(:,D+1,i),fi(:,D+2,i),lty(i),'LineWidth',lw)    %Ev=f(t)
            plot(ax(4),fi(:,D+1,i),fi(:,D+3,i),lty(i),'LineWidth',lw)    %Ed=f(t)      
        end
    case 'MULTI'
        %PLOTS
        C=app.PlotColors;
        for i=1:size(pD,2)
            fi=pD(i).Results.Strain;
            if D==3
                plot(ax(5),fi(:,1,1),fi(:,2,1),'Color',C(i,:))           %Ex=f(t)
            end
            plotMark(app,ax(1),fi(:,1,1),fi(:,D,1),'Color',C(i,:))       %Ey=f(t)
            plotMark(app,ax(2),fi(:,1,1),fi(:,D+1,1),'Color',C(i,:))     %Ez=f(t)
            plotMark(app,ax(3),fi(:,D+1,1),-fi(:,D+2,1),'Color',C(i,:))	 %Ev=f(Ez)
            plotMark(app,ax(4),fi(:,D+1,1),fi(:,D+3,1),'Color',C(i,:))	 %Ed=f(Ez)
            plotMark(app,ax(end),fi(:,D+1,3),-fi(:,D+2,1),'Color',C(i,:))  %Ev=f(EzEXT)
        end
        %LEGENDS
        str=cat(1,pD.FileName);
end

%LEGENDS
if leg
    for i=1:numel(ax)
        legend(ax(i),str,'location','best')
    end
end

%Prepare path
png=".png";
%PLOT Ex=f(t)
%Check if 3d or 2D for X calculation
if D==3
    if tit;title(ax(5),'Evolution of the Strain X');end
    ylabel(ax(5),'Strain X')
    xlabel(ax(5),'Time(s)')
    fnm=pD(1).Prefix+"StrainX";
    saveas(f(5),fullfile(path,fnm+png));
end

if tit
    title(ax(1),'Evolution of the Strain Y')
    title(ax(2),'Evolution of the Strain Z')
    title(ax(3),'Evolution of the Volumetric Strain ')
    title(ax(4),'Evolution of the Deviatoric Strain ')
end

%%PLOT Ey=f(t)
ylabel(ax(1),'Strain Y')
xlabel(ax(1),'Time(s)')
fnm=pD(1).Prefix+"Strain_Y";
saveas(f(1),fullfile(path,fnm+png));

%%PLOT Ez=f(t)
ylabel(ax(2),'Axial Strain')
xlabel(ax(2),'Time(s)')
fnm=pD(1).Prefix+"Strain_Z";
saveas(f(2),fullfile(path,fnm+png));

%%PLOT Ev=f(Ez)
ylabel(ax(3),'Volumetric Strain')
ax(3).YDir='reverse';
xlabel(ax(3),'Axial Strain')
fnm=pD(1).Prefix+"Strain_Vol";
saveas(f(3),fullfile(path,fnm+png));

%%PLOT Ed=f(Ez)
ylabel(ax(4),'Deviatoric Strain')
xlabel(ax(4),'Axial Strain')
fnm=pD(1).Prefix+"Strain_Dev";
saveas(f(4),fullfile(path,fnm+png));

if isequal(upper(mode),'MULTI')
    %%PLOT Ev=f(EzEXT)
    ax(end).YDir='reverse';
    if tit;title(ax(end),'Evolution of the Volumetric Strain ');end
    ylabel(ax(end),'Volumetric Strain')
    xlabel(ax(end),'Axial Strain')
    if leg;legend(ax(end),str,'location','best');end
    fnm=pD(1).Prefix+"Strain_Vol_Ext";
    saveas(f(end),fullfile(path,fnm+png));
end

%PLot inertal term evolution if possible
if ~isempty(pD.Results.Inertia) && numel(pD)==1
    I=cat(2,pD.Results.Inertia{:});
    %calculate mean value
    mnI=mean(I,1);
    %calculate standart deviation in log scale
    eT=log10(std(I,1));
    %prepare scatter
    nbg=numel(pD.Results.Inertia{1});
    str=pD.Results.Strain(:,D+1,end);
    sctI=[cat(1,pD.Results.Inertia{:})  repelem(str,nbg)  ];
    %prepare inf points
    vez={};
    for i=1:numel(pD.InfPts.q)
        vez=[vez,{'Pointx',pD.InfPts.ez(i)}]; %#ok<AGROW>
    end
    %plot 1 scatter
    k=figure;axk=axes(k);hold(axk,'on');
    scatter(axk,sctI(:,2),log10(sctI(:,1)))
    plotMark(app,axk,str,log10(mnI),vez{:})
    yl=yline(axk,-3,':','Color','k','LineWidth',1.5);
    set(get(get(yl,'Annotation'),'LegendInformation'),...
        'IconDisplayStyle','off');
    yl=yline(axk,-2,':','Color','k','LineWidth',1.5);
    set(get(get(yl,'Annotation'),'LegendInformation'),...
        'IconDisplayStyle','off');
    ylabel(axk,'log(Inertial Number)')
    xlabel(axk,'Axial strain')
    fnm=pD(1).Prefix+"Inertia_Scat";
    saveas(k,fullfile(path,fnm+png));
    
    %plot 2 standart dev
    k(2)=figure;axk(2)=axes(k(2));hold(axk(2),'on');
    errorbar(axk(2),str,log10(mnI),-eT,eT)
    plotMark(app,axk(2),str,log10(mnI),vez{:})
    yl=yline(axk(2),-3,':','Color','k','LineWidth',1);
    set(get(get(yl,'Annotation'),'LegendInformation'),...
        'IconDisplayStyle','off');
    yl=yline(axk(2),-2,':','Color','k','LineWidth',1);
    set(get(get(yl,'Annotation'),'LegendInformation'),...
        'IconDisplayStyle','off');
    ylabel(axk(2),'log(Inertial Number)')
    xlabel(axk(2),'Axial Strain')
    %axk(2).YLim=[-5,-.5];
    fnm=pD(1).Prefix+"Inertia_Std";
    saveas(k(2),fullfile(path,fnm+png));
    
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
if isequal(app.LIGGGHTSAnalysisButtonGroup.SelectedObject,...
        app.ExeAllButton) || nargin>3
    delete(f);
    try delete(k);
    catch 
    end
end
end
function w2ClustPlotter(app,pD)
%W2CLUSTPLOTTER Plots the distribution of clusters in negative W2 grains

%create figures
nb=9;
f(nb)=figure;ax(nb)=axes(f(nb));hold(ax(nb),'on');
for i=1:(nb-1)
    f(i)=figure;ax(i)=axes(f(i));hold(ax(i),'on'); %#ok<LAXES>
end
%Load variables
cl=pD.Results.Cluster;
str=pD.Results.Strain;
rNb=pD.Results.Nb;
%path
path=MakePath(app,'LOOPW2');
%Check title and legends option
if app.TitlesCB.Value;tit=1;else;tit=0;end
if app.LegendsCB.Value;leg=1;else;leg=0;end

%Clusters 
C=[0 0.4470 0.7410; %matlab default colors
    0.8500 0.3250 0.0980;
    0.9290 0.6940 0.1250;
    0.4940 0.1840 0.5560];
for i=1:2
    for j=1:2
        if j==1
            plotMark(app,ax(4*(i-1)+j),str,cl(:,1,i),'Color',C(1,:))%Order 4=f(Ez)
        end
        plotMark(app,ax(4*(i-1)+j),str,cl(:,2,i),'Color',C(2,:)) %Order 6=f(Ez)
        plotMark(app,ax(4*(i-1)+j),str,cl(:,3,i),'Color',C(3,:)) %Order 8-20=f(Ez)
        plotMark(app,ax(4*(i-1)+j),str,cl(:,4,i),'Color',C(4,:)) %Order 22=f(Ez)
    end
    totCl=sum(cl(:,:,i),2);
    %graph ratio
    for j=3:4
        if j==3
            plotMark(app,ax(4*(i-1)+j),str,cl(:,1,i)./totCl,'Color',C(1,:))%Order 4=f(Ez)
        end
        plotMark(app,ax(4*(i-1)+j),str,cl(:,2,i)./totCl,'Color',C(2,:)) %Order 6=f(Ez)
        plotMark(app,ax(4*(i-1)+j),str,cl(:,3,i)./totCl,'Color',C(3,:)) %Order 8-20=f(Ez)
        plotMark(app,ax(4*(i-1)+j),str,cl(:,4,i)./totCl,'Color',C(4,:)) %Order 22=f(Ez)
    end
    if leg
        legend(ax(4*(i-1)+1),"Order 4","Order 6","Order 8-20",...
            "Order 22+",'location','best')
        legend(ax(4*(i-1)+2),"Order 6","Order 8-20",...
            "Order 22+",'location','best')
        legend(ax(4*(i-1)+3),"Order 4","Order 6","Order 8-20",...
            "Order 22+",'location','best')
        legend(ax(4*(i-1)+4),"Order 6","Order 8-20",...
            "Order 22+",'location','best')
    end
end
%nb of W2grains
plotMark(app,ax(end),str,rNb);

infP=pD.InfPts;
for k=1:numel(infP) %nb of pts
    if infP(k)==0 || infP(k)>=str(end) || infP(k)<=str(1)
        continue;
    end
    for j=1:nb %nb of plots
        xl=xline(ax(j),infP(k),'--','Color','#C1C1C1');
        set(get(get(xl,'Annotation'),'LegendInformation'),...
            'IconDisplayStyle','off');
        x2=xline(ax(j),infP(k),'--','Color','#C1C1C1');
        set(get(get(x2,'Annotation'),'LegendInformation'),...
            'IconDisplayStyle','off');
    end
end
    
%titles
if tit
    strg=["Cluster for negative W2 grains";
        "Cluster for positive W2 grains";
        "Ratio of negative W2 grains"];
    for i=1:nb
        title(ax(i),strg(ceil(i/4)))
    end
end
%get correct yscales between fc and nfc grains
for i=1:((nb-1)/2)
    mx=max(cat(2,ax(i+[0,4]).YLim));
    ax(i).YLim=[0 mx];
    ax(i+4).YLim=[0 mx];
end
%Labels and save
fnm=["W2Cl";"W2ClN4";"W2ClPct";"W2ClN4Pct";...
    "NW2Cl";"NW2ClN4";"NW2ClPct";"NW2ClN4Pct"];
for i=1:(nb-1)
    if ismember(i,[1,2,5,6])
        ylabel(ax(i),'Number of Cluster')
    else
        ylabel(ax(i),'Ratio of Cluster')
    end
    xlabel(ax(i),'Axial Strain')
    f(i).Position(3:4)=app.Figsize;
    saveas(f(i),fullfile(path,fnm(i)+".png"));
end

%last figure
ylabel(ax(end),'Ratio of grains')
xlabel(ax(end),'Axial Strain')
f(end).Position(3:4)=app.Figsize;
saveas(ax(end),fullfile(path,"NbW2Gr.png"));
%Delete figures in the case of exe all
if isequal(app.LIGGGHTSAnalysisButtonGroup.SelectedObject,app.ExeAllButton)
    delete(f);
end
end
%{
function testplotshit(pD)
file=pD.Results;
%CHECK DIMENSION
if pD.Bool3D;D=3;else; D=2;end
nb=3;

f(nb)=figure;ax(nb)=axes(f(nb));hold(ax(nb),'on');
for i=1:(nb-1)
    f(i)=figure;ax(i)=axes(f(i));hold(ax(i),'on'); %#ok<LAXES>
end
x1=file(2:end,2,1)-file(1:end-1,2,1);
x2=file(2:end,2,end)-file(1:end-1,2,end);
%plot(ax(1),file(1:end-1,1,1),x1,'-')
plot(ax(1),file(1:end-1,1,2),x2./x1,'-x')
title(ax(1),'ratio dEX Bagi/Ext')
%Ey=f(t)
y1=file(2:end,D,1)-file(1:end-1,D,1);
y2=file(2:end,D,end)-file(1:end-1,D,end);
%plot(ax(2),file(1:end-1,1,1),y1,'-')
plot(ax(2),file(1:end-1,1,2),y2./y1,'-x')
title(ax(2),'ratio dEY Bagi/Ext')
%Ez=f(t)
z1=file(2:end,D+1,1)-file(1:end-1,D+1,1);
z2=file(2:end,D+1,end)-file(1:end-1,D+1,end);
%plot(ax(3),file(1:end-1,1,1),z1,'-')
plot(ax(3),file(1:end-1,1,2),z2./z1,'-x')
title(ax(3),'ratio dEZ Bagi/Ext')

end
%}