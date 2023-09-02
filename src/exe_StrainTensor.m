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
        strainLoad(app,PD);return;
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
[N1,N2,interval,stepArray,nbFiles] = createStepArray(app);
%strain is the relative deformation between two steps, thus the nb of
%calculations is given by :
nbFiles=nbFiles-1;

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
    %creating a folder for the outputfiles (vtk) if vtk files were asked
    if app.StnVTKCB.Value
        pathVTK=MakePath(app,'Strain')+"ShearBandsVTK"+N1+"to"+N2+"int"+...
            interval+"/";
        if exist(pathVTK,'dir')==0;mkdir(pathVTK);end
    end
    % values calculated are incremental, thus will be cumulated in the
    % following vectors,
    pgEd=zeros(1,1,app.NbGrains);
    W2=zeros(1,1,app.NbGrains);
    %Cluster and strain calculation
    if app.StnClusterCB.Value
        resClStr=zeros(2000,nbFiles,2);
    end
    %Prepare per grain inertia calculation
    itvArray=stepArray(2:end)-stepArray(1:end-1);
    I=cell(nbFiles,1);

    %Try loading per grain average cluster value
    aveCl=0;
    fnm=fullfile(MakePath(app,'LOOPAC'),"Average_Cluste_pGrain_"+N1+"to"+N2+"int"+interval+".txt");
    if isfile(fnm)
        %Prepare a matrix that will contain the aveCluster and W2 values
        %for each grain
        pgRes=zeros(app.NbGrains,nbFiles,2);
        aveCl=1;
        opts = detectImportOptions(fnm,'NumHeaderLines',4);
        opts.Delimiter='|';
        opts.VariableNamesLine = 5;
        Table = readtable(fnm,opts);
        pgRes(:,:,2)=Table{2:end,3:end};
    end
end

%Load first object containing grains data
[grT,PDT]=grains('STENSOR',stepArray(1),PD,app);
if isempty(grT.Coord)
    warndlg(['gr.Coord is empty on step ' N1]);return;
end

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
    stepT=stepArray(i);
    stepTdT=stepArray(i+1);
    
    if isa(app.CalculatingPanel,'double')
        app=CalcPanel(app,i,nbFiles,"NEWLINE");
    else
        app=CalcPanel(app,i,nbFiles,stepTdT);
    end
    %Load second object containing grains data
    [grTdT,PDTdT]=grains('STENSOR',stepTdT,PD,app);
    if isempty(grTdT.Coord)
        CalcPanel(app,'','','','off');
        warndlg(['Gr.Coord is empty on step ' stepTdT]);return;
    end
    
    %check if spaceCell file exist from a previous execution. If not launch
    %the calculation
    fl=[pathSc char("spaceCellsfile"+stepTdT+"int"+interval+".mat")];
    if isfile(fl)
        %fprintf('Space cell %d file found and loaded \n',step)
        sc=load(fl).sc;sv=0;
    else
        %Strain calculation should be made using the displacement between
        %T-TdT of the grain located in T
        sc = spaceCellSystem(type,stepArray(i),grTdT,grT,app,PDT);sv=1;
    end
    
    %GLOBAL CALCULATION
    if ismember(1,calcType)
        %check if the tensor exists or need to be calculated
        sc=gStrainTensor(sc,app,PDT);
        %calculate total volume and save values into matrix
        if ~app.SubdivisionButton.Value
            %divide the strain value by hte volume of the previus step.
            %However for the first interval it is not calculated so divide
            %it by the volume of the actual step.
            if isempty(V_old);V_old=sum(sc.CellVol);end 
            stG(:,:,1)=(sc.GStTensor)/V_old;
            stG(:,:,2)=sum(sc.CellStn,3)/V_old;
            Results(i,1,3)=stepTdT;
            tri=2;
            V_old=sum(sc.CellVol);
        else
            [~,V_sub]=getVolume(app.TrialData,app,(stepTdT),PDT); %v of previous step
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
            Results(i,:,1:tri)=[ones(1,1,tri)*stepTdT...
                -[stG(1,1,1:tri) stG(2,2,1:tri) stG(3,3,1:tri) iEv ] iEd ]; %- for negative compression, GC standart
        else
            iEv=(stG(1,1,1:tri)+stG(2,2,1:tri));
            iEd=sqrt(( (stG(1,1,:)-stG(2,2,:)).^2+...     %(sig1-sig2)^2
                6*stG(1,2,:).*stG(1,2,:) )/2);          % 2*(sig12)^2
            Results(i,:,1:tri)=[ones(1,1,tri)*stepTdT...
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
        %permtute. PGStressTensor should be divided by vol of the sphere to
        %return the correct unity for W2
        incW2=(grTdT.PGStressTensor-grT.PGStressTensor).*sc.PStTensor;
        incW2=permute(sum(incW2,1:2),[3,1,2]);%./(grT.Radius.^3*4/3*pi()); remove the volume division for now
        W2=W2+incW2;
        %calculate the cumulated strain
        iSTG=sc.PStTensor;  %incr grain strain (BAGI96)
        iSTC=sc.CellStn;    %incr cell strain (Interpolation)
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

        %Extra calculation : Cluster and second order work realtion
        if app.StnClusterCB.Value
            rClStr=cellStrainCluster(app,stepT, sc, iEdC);
            if isempty(rClStr);return;end
            resClStr(1:size(rClStr,1),i,:)=rClStr;
        end
        
        cord=grT.Coord;
        r=grT.Radius;
        if aveCl
            %IF average cluste file was found, save per grain W2 results
            pgRes(:,i,1)=incW2;
        end
                
        %Calculate the intertial number to verify if we still are in the
        %quasi-static regime
        P=permute(abs(grT.PGStressTensor(1,1,:)+...
            grT.PGStressTensor(2,2,:)+...
            grT.PGStressTensor(3,3,:))/3,[3,2,1]);
        V=4/3*pi()*r.^3;
        I{i}=permute(iEd,[3,2,1])/(itvArray(i)*app.TimeStep)...
            .*sqrt(2600*V.^2./(P.*(2*r)));
        %save as vtk file -- skip plot if asked
        if ~app.SubdivisionButton.Value && app.StnVTKCB.Value
            %%%%% GRAIN VTK FILE %%%%
                %remove positive values and calculate a 90% negative
                %percentile of the values to be ploted. Thus a scale will
                %be created.

            if isempty(PDT)
                %ind=find(~isnan(sc.GrainVolume));
                ind=1:sc.NbG;
            else
                ind=PDT.GrainsRectangle;
            end
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
            
            %prepare Scalar data
            svData={'SCALARS','Radius',r(ind),...
                'SCALARS','Dev_Str',pgEd(ind),...
                'SCALARS','Dev_Str_Rat',pgEd(ind)/max(pgEd(ind)),...
                'SCALARS','Dev_Str_Inc',iEd(ind),...
                'SCALARS','W2',pW2(ind),...
                'SCALARS','W2_Inc',incpW2(ind)};

            if aveCl
                %Add ave cluster
                svData={svData,'SCALARS','W2_Inc',pgRes(ind,i,2)}; %#ok<AGROW>
            end

            newFile=fullfile(pathVTK,prefix+"ShearBand"+stepTdT+".vtk");
            vtkwrite(newFile,'unstructured_grid',...
                cord(ind,1),cord(ind,2),cord(ind,3),...
                svData{:},...
                'PRECISION',10,'BINARY');
            %'SCALARS','DevStrainRatio',pcEv(ind)/max(pcEv(ind)),...
            
            %%%%% CELL VTK FILE %%%%
            
            %Inc dev strain of cell with too low volume will be turned to 0
            cv=cat(1,sc.CellVol);
            cv=(cv<mean(cv)^2/max(cv));
            iEdC(cv)=0;
            %get value that represents  99% of cells and turn into the max
            %value
            k=sort(iEdC);
            k=k(floor(0.99*numel(k)));
            iEdC=iEdC/abs(k);
            iEdC(iEdC>1)=1;
                %calculate per cell VR
            [vr,~,~]=perCellVoidRatio(sc,grT);
                %prepare to create vtk files
            newFile=fullfile(pathVTK,prefix+"ShearBandCell"+stepTdT+".vtk");
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
        elseif app.StnVTKCB.Value
            for l=1:lMax
                for tri=1:cMax
                    j=tri+cMax*(l-1);
                    %prepare Scalar data
                    svData={'SCALARS','Radius',r(PDT.SubGrains{j}),...
                        'SCALARS','Dev_Strain',Ed(PDT.SubGrains{j}),...
                        'SCALARS','Dev_Strain_Inc',iEd(PDT.SubGrains{j}),...
                        'SCALARS','Dev_Strain_Ratio',Ed(PDT.SubGrains{j})/mean(Ed(PDT.SubGrains{j})),...
                        'SCALARS','SecOrderWork',incW2(PDT.SubGrains{j})};

                    if aveCl
                        %Add ave cluster
                        svData={svData,'SCALARS','W2_Inc',pgRes((PDT.SubGrains{j}),i,2)}; %#ok<AGROW>
                    end
                    newFile=fullfile(pathVTK,prefix+"ShearBandl"+l+"c"+tri+"-"+stepTdT+".vtk");
                    vtkwrite(newFile,'unstructured_grid',...
                        cord(PDT.SubGrains{j},1),cord(PDT.SubGrains{j},2),cord(PDT.SubGrains{j},3),...
                        svData{:},...
                        'PRECISION',10,'BINARY');
                end
            end
        end
    end
    %Save SC in a file for a faster re execution if needed
    if sv==1 || sv==0
        sc=purge(sc,'Strain');
        save(fl,'sc','-v7.3');
    end
    %grTdT becomes the next grT
    grT=grTdT;
    PDT=PDTdT;
end
CalcPanel(app,i+1,nbFiles,'','off');

if ismember(2,calcType) && app.SubdivisionButton.Value
    %Create VTK lines separating the subdvision for better analysis
    v=PDT.SubVertices;
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
    fnm=fullfile(pathVTK,"subdivisionsLines.vtk");
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
    fnm=fullfile(pathVTK,"vtkVariables.txt");
    stepArray=stepArray(2:end);
    stress = extStress(app.TrialData,stepArray,app);
    strain = extStrains(app.TrialData,stepArray,N1,app,'allCalc');
    StepEzEvQP=[stepArray,strain(:,[D,D+1]),stress(:,end-1:end)];
    vtkLog(app,'FCbase',fnm,StepEzEvQP)

    if aveCl
        acRes.Strain=strain(:,D);
        acRes.Stress=stress(:,end);
        acRes.W2=pgRes(:,:,1);
        acRes.AveCluster=pgRes(:,:,2);
        pD = plotData("Normal",acRes,app,'',consoStrain(end));
        fl=fullfile(MakePath(app,'STNW2AC'),"Strain-ClW2"+N1+"to"+N2+"int"+interval+".mat");
        save(fl,'pD','-v7.3');
        %go to plot
        aveClW2Plotter(app,pD);
    end
    if app.StnClusterCB.Value
        %Cluster and Second order work results
        stnCl.SC = resClStr;
        stnCl.Strain = strain(:,D);
        stnCl.Stress = stress(:,end);
        %create plotData object
        pD = plotData("Normal",stnCl,app,'',consoStrain(end));
        %save pD object
        fl=fullfile(MakePath(app,'STNCL'),"Strain-Cl"+N1+"to"+N2+"int"+interval+".mat");
        save(fl,'pD','-v7.3');
        %go to plot
        strainClustPlotter(app,pD);
    end
end
end
%execution function
function resClStr=cellStrainCluster(app,step, scStn, iEdC)
%Load Cluster data
scFnm=[app.SavePath '/MatlabResultsFile/SpaceCellFiles/'...
    char("LoopsSpaceCellsfile"+step+".mat")];
try scClu=load(scFnm).sc;
catch
    fprintf(['No Loops spacecell files, '...
        'please run loops calculation beforehand.'])
    resClStr='';return
end

%%% Prepare strain Data
%Inc dev strain of cell with too low volume will be removed from
%calculation - ratio (min vol)/(mean vol) should be the same as
%(mean)/(max)
cv=cat(1,scStn.CellVol);
cv=(cv<mean(cv)^2/max(cv));
iEdC(cv)=0;
%get value that represents  99% of cells and turn into the max
%value
k=sort(iEdC);
k=k(floor(0.99*numel(k)));
iEdC=iEdC/abs(k);
%transform from a 1x1xN to Nx1x1
iEdC=permute(iEdC,[3,2,1]);

%%% Prepare Cluster data
gC=scClu.GoodCells;
if isempty(gC)
    gC = goodCell(scClu);
end
%Get clusters cells
cC=cat(2,scClu.Loops.sCells)';
%Transform them into goodcell Ids
[~,cC]=ismember(cC,gC);
%add cluster order to it
clCe=[cC repelem(cat(1,scClu.Loops.Order),cat(1,scClu.Loops.nbCells),1)];
%add cluster 4s
clCe=sortrows([clCe;[cat(1,scClu.Clt4.sCells) 4*ones(scClu.Clt4.nbCells,1)]]);

maxO=max(cat(1,scClu.Loops.Order));

%prepare results vector : count how many cells of each order are above or
%below the dev stn limit
stnL=.10*abs(k);
resClStr=zeros(maxO/2,1,2);
for i=1:2
    %get cells above or below dev strain limit
    if i==1
        chk=iEdC>=stnL;
    else
        chk=iEdC<stnL;
    end
    %get clusters of these cells
    cls=clCe(chk,2);
    %count repetitions
    cls=accumarray(cls,ones(size(cls)));
    %save into results
    resClStr(1:size(cls,1)/2,:,i)=cls(2:2:size(cls,1));
end


end
%load functions
function strainLoad(app,type)
%STRAINLOAD Loader for the strain function
switch upper(type)
    case 'STRAIN'
        %Load file
        [fnm,fPath]=MatLoader(type,app);
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
            fprintf('Executing several plots\n')
            for i=1:numel(pD)
                mode='NORMAL';
                path=fullfile(MakePath(app,'STRAIN'),pD(i).FileName);
                if exist(path,'dir')==0;mkdir(path);end
                strainPloter(app,mode,pD(i),path)
            end
            fprintf('Done\n')
        end
    case 'STNCL'
        bsPth=MakePath(app,type,'check');
        if exist(bsPth,'dir')==0
            fprintf('Strain-Cluster usual path not found \n')
            bsPth='';
        end
        [fnm,Path,~] = uigetfile( ...
            {'*.mat','Matlab files(*.mat)';'*.*',  'All Files (*.*)'},...
            'Select file for loading','Multiselect','on',bsPth);
        if Path==0
            warndlg('No files were chosen, Try again.')
        else
            pD=load(fullfile(Path,fnm)).pD;
            strainClustPlotter(app,pD)
        end
    case 'STNW2AC'
        bsPth=MakePath(app,type,'check');
        if exist(bsPth,'dir')==0
            fprintf('Average cluster W2 usual path not found \n')
            bsPth='';
        end
        [fnm,Path,~] = uigetfile( ...
            {'*.mat','Matlab files(*.mat)';'*.*',  'All Files (*.*)'},...
            'Select file for loading','Multiselect','on',bsPth);
        if Path==0
            warndlg('No files were chosen, Try again.')
        else
            pD=load(fullfile(Path,fnm)).pD;
            aveClW2Plotter(app,pD)
        end

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
set(0,'defaultAxesFontSize',app.FontSizeEF.Value)

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
        if numel(pD)<8
            C = app.PlotColors;
        else
            C = graphClrCode(size(pD,2));%plot colorcode
        end
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

%Plot inertial term evolution if possible
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
function aveClW2Plotter(app,pD)
%STRAINCLUSTPLOTTER Plots the distribution of clusters cells and high dev
%strain
%Check title and legends option

if app.TitlesCB.Value;tit=1;else;tit=0;end

%If any of the loaded files is a Qcst simultaion, the plot will be done
%using pressure as x-axis

%prepare paths
path=MakePath(app,'STNW2AC');
png=".png";

%Plot histogram for last step, inflection points and a pt
%betweem inflection and last step
k=numel(pD.InfPts.q)+2;
pts=[zeros(1,k-1) size(pD.Results.AveCluster,2)];
for i=1:numel(pD.InfPts.q)
    [~,pos]=min(abs(pD.Results.Strain-pD.InfPts.ez(i)));
    pts(i)=pos;
end
pts(end-1)=ceil((pts(end)+pts(end-2))/2);
fntSz=app.FontSizeEF.Value;


%Prepare Figures
nbPlt=4;
nb=nbPlt*k+1;
f(nb)=figure;ax(nb)=axes(f(nb));hold(ax(nb),'on');ax(nb).FontSize = fntSz;
for i=1:(nb-1)
    f(i)=figure;ax(i)=axes(f(i));hold(ax(i),'on'); %#ok<LAXES>
    ax(i).FontSize = fntSz;
end

maxCl=30;
%Plot histograms
for i=1:k
    aveCl=pD.Results.AveCluster(:,pts(i));  %pg ave cl value
    w2=pD.Results.W2(:,pts(i));             %pg w2 value
    %limit aveCl to 30 to reduce nb of bins
    aveCl(aveCl>maxCl)=maxCl;
    hstgEdges=(3.5:1:(ceil(max(aveCl))+0.5));
    %plot histogram of all grains
    histogram(ax((i-1)*nbPlt+1),aveCl,hstgEdges,'Normalization','probability') %'DisplayStyle','stairs'
    %plot W2<0 histrograms
    histogram(ax((i-1)*nbPlt+2),aveCl(w2<0),hstgEdges,'Normalization','probability')
    %plot W2>=0 histograms
    histogram(ax((i-1)*nbPlt+3),aveCl(w2>=0),hstgEdges,'Normalization','probability')
    %Plot both w2 lines
    h1=histcounts(aveCl(w2<0),'BinEdges',hstgEdges,'Normalization','probability');
    h2=histcounts(aveCl(w2>=0),'BinEdges',hstgEdges,'Normalization','probability');
    x=4:1:hstgEdges(end)-0.5;
    plot(ax((i-1)*nbPlt+4),x,h1,x,h2,"LineWidth", app.PlotWidthEF.Value)

    if tit; title(ax((i-1)*nbPlt+1),"Average cluster at strain " + pD.Results.Strain(pts(i)));end
    xlabel(ax((i-1)*nbPlt+1),'Average Cluster Order')
    ylabel(ax((i-1)*nbPlt+1),'Numeber of grains')
    ax((i-1)*nbPlt+1).XTick=4:max(2,ceil((ax((i-1)*nbPlt+1).XLim(2)-4)/20)*2):2000;
    fnm="St_"+pD.Results.Strain(pts(i))+"_Total_Histogram";
    saveas(f((i-1)*nbPlt+1),fullfile(path,fnm+png));

    if tit; title(ax((i-1)*nbPlt+2),"Average cluster at strain " + pD.Results.Strain(pts(i))+" (W2<0)");end
    xlabel(ax((i-1)*nbPlt+2),'Average Cluster Order')
    ylabel(ax((i-1)*nbPlt+2),'Numeber of grains')
    ax((i-1)*nbPlt+2).XTick=4:max(2,ceil((ax((i-1)*nbPlt+2).XLim(2)-4)/20)*2):2000;
    fnm="St_"+pD.Results.Strain(pts(i))+"_W2_Neg_Histogram";
    saveas(f((i-1)*nbPlt+2),fullfile(path,fnm+png));

    if tit; title(ax((i-1)*nbPlt+3),"Average cluster at strain " + pD.Results.Strain(pts(i))+" (W2>=0)");end
    xlabel(ax((i-1)*nbPlt+3),'Average Cluster Order')
    ylabel(ax((i-1)*nbPlt+3),'Numeber of grains')
    ax((i-1)*nbPlt+3).XTick=4:max(2,ceil((ax((i-1)*nbPlt+3).XLim(2)-4)/20)*2):2000;
    fnm="St_"+pD.Results.Strain(pts(i))+"_W2_Pos_Histogram";
    saveas(f((i-1)*nbPlt+3),fullfile(path,fnm+png));

    if tit; title(ax((i-1)*nbPlt+4),"Average cluster at strain " + pD.Results.Strain(pts(i))+" (W2>=0)");end
    xlabel(ax((i-1)*nbPlt+4),'Average Cluster Order')
    ylabel(ax((i-1)*nbPlt+4),'Numeber of grains')
    legend(ax((i-1)*nbPlt+4),"W2<0","W2>=0")
    ax((i-1)*nbPlt+4).XTick=4:max(2,ceil((ax((i-1)*nbPlt+4).XLim(2)-4)/20)*2):2000;
    ax((i-1)*nbPlt+4).XLim(1)=4;
    fnm="St_"+pD.Results.Strain(pts(i))+"_Curve";
    saveas(f((i-1)*nbPlt+4),fullfile(path,fnm+png));

end

%plot evolution of nb of w2<0 grains
if pD.SimType==3
    xlab='Mean pressure (kPA)';
    x=pD.Results.Stress;
else
    xlab='Axial strain';
    x=pD.Results.Strain;
end
plot(ax(end),x,sum(pD.Results.W2<0,1)/size(pD.Results.W2,1),"LineWidth",app.PlotWidthEF.Value)
if tit; title(ax(end),"Ratio of W2<0 grains");end
ylabel(ax(end),'Ratio of Grains')
xlabel(ax(end),xlab)
fnm="Evol_W2_Nb";
saveas(f(end),fullfile(path,fnm+png));

%Delete figures in the case of exe all
if isequal(app.LIGGGHTSAnalysisButtonGroup.SelectedObject,app.ExeAllButton)
    delete(f);
end
end
function strainClustPlotter(app,pD)
%AVECLW2PLOTTER Function to plot aveCl and W2 data
%Check title and legends option

if app.TitlesCB.Value;tit=1;else;tit=0;end

%If any of the loaded files is a Qcst simultaion, the plot will be done
%using pressure as x-axis

%prepare paths
path=MakePath(app,'STNCL');
png=".png";

%Plot histogram for last step, inflection points and a pt
%betweem inflection and last step
k=numel(pD.InfPts.q)+2;
pts=[zeros(1,k-1) size(pD.Results.SC,2)];
for i=1:numel(pD.InfPts.q)
    [~,pos]=min(abs(pD.Results.Strain-pD.InfPts.ez(i)));
    pts(i)=pos;
end
pts(end-1)=ceil((pts(end)+pts(end-2))/2);
fntSz=app.FontSizeEF.Value;

%Prepare Figures
nbPlt=1;
nb=nbPlt*k;
f(nb)=figure;ax(nb)=axes(f(nb));hold(ax(nb),'on');ax(nb).FontSize = fntSz;

for i=1:(nb-1)
    f(i)=figure;ax(i)=axes(f(i));hold(ax(i),'on'); %#ok<LAXES>
    ax(i).FontSize = fntSz;
end

%maxCl=30;
%Plot histograms
for i=1:k
    clStn=pD.Results.SC(:,pts(i),:);  %pg ave cl value
    % find last 0
    clStn=clStn(1:find(clStn(:,:,1)~=0 | clStn(:,:,2)~=0,1,"last")+1,:,:);
    %limit aveCl to 30 to reduce nb of bins
    %clStn(clStn>maxCl)=maxCl;
    O=(2:1:size(clStn,1))*2;
    plot(ax((i-1)*nbPlt+1),O,clStn(2:end,:,1)/sum(clStn(2:end,:,1)),...
        O,clStn(2:end,:,2)/sum(clStn(2:end,:,2)),...
        "LineWidth", app.PlotWidthEF.Value)

    if tit; title(ax((i-1)*nbPlt+1),"Cluster cell distribution at strain: " + pD.Results.Strain(pts(i)));end
    xlabel(ax((i-1)*nbPlt+1),'Cluster Order')
    ylabel(ax((i-1)*nbPlt+1),'Ratio of grains')
    legend(ax((i-1)*nbPlt+1),"Above lim","Below lim")
    ax((i-1)*nbPlt+1).XTick=4:max(2,ceil((ax((i-1)*nbPlt+1).XLim(2)-4)/20)*2):2000;
    ax((i-1)*nbPlt+1).XLim(1)=4;
    fnm="St_"+pD.Results.Strain(pts(i))+"_Cell_Cluster";
    saveas(f((i-1)*nbPlt+1),fullfile(path,fnm+png));

end

end
