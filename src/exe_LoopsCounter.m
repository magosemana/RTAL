function exe_LoopsCounter(app,varargin)
%LoopsCounter Will read LIGGGHTS files and export data on the Loops
%   This function will use the class 'grains' to read files and then
%   calculate the Loops cointained in the trial. We will then give the
%   possibility of writing these loops in a VTK file for better visibility.

%% Part 1 - Preparation
if nargin>1
    if isequal(upper(varargin{2}),'LOAD');LoopsLoad(app,varargin{1});return;end %check Load
end

PartData='';vtkOF=0;
if nargin>1;PartData=varargin{1};end
if nargin>2;vtkOF=varargin{2};end
%Check startup options
if isempty(PartData);prefix="total";
else; prefix="partial";
end

%check dimesion
D= app.Bool3D+2;

%Load values
[N1,N2,interval,stepArray,nbFiles] = createStepArray(app);

%prepare SpaceCellFiles directory
pathSc=MakePath(app,'SCF');

%start variables
aniOF=0;defOF=0;vrOF=0;strOF=0;trOF=0;stnOF=0;aveOF=0;doVol=1;
%Check VTK
if isequal(upper(vtkOF),'ON')
    vtkOF=1;
    if D==2
        fileName=["Loop3_","Loop4_","Loop5_","Loop+6_"];
        dirName=prefix+"LoopsVTK"+N1+"to"+N2+"int"+interval;
    else
        dirName=prefix+"ClusterVTK"+N1+"to"+N2+"int"+interval;
    end
    pathVTK=fullfile(app.SavePath,"MatlabResultsFile",dirName,"/");
    if exist(pathVTK,'dir')==0
        mkdir(pathVTK)
    end
else
    vtkOF=0;
end
%prepare base files
if D==2
    baseRes=zeros(nbFiles,7); %steps (1) Loops(4) MeanVal(2)
else
    maxCl=10000;
    baseRes=zeros(maxCl/2,nbFiles,4); %on 3rd dim : order/size/cells
    %prepare anisotropy
    if isequal(app.LPAniSwitch.Value,'On')
        aniOF=1;
        aniRes=cell(nbFiles,1);
    end
    %prepare Deformability
    if isequal(app.LPDefSwitch.Value,'On')
        defOF=1;
        defRes=cell(nbFiles,2);
    end
    %prepare VoidRatio
    if isequal(app.LPVRSwitch.Value,'On')
        vrOF=1;
        vrData=zeros(maxCl,nbFiles,2);
        vrTot=zeros(nbFiles,1);
        resVRZ=cell(nbFiles,2);
    end
    %prepare Stress
    if isequal(app.LPStrSwitch.Value,'On')
        strOF=1;
        strRes=cell(nbFiles,1);
    end
    %prepare ClusterTransformation
    if isequal(app.LPCTrSwitch.Value,'On')
        trOF=1;
        trRes=cell(nbFiles,2);
    end
    %Prepare Strain
    if isequal(app.LPStnSwitch.Value,'On')
        stnOF=1;
    end
    %prepare Average cluster order
    if isequal(app.LPACSwitch.Value,'On')
        aveOF=1;
        aveCl=zeros(app.NbGrains,nbFiles);
        pathAve=MakePath(app,'LOOPAC');
    end
end

%% Part 2 - Execution

%Turn on calculation pannel
app=CalcPanel(app,'',nbFiles,'Starting calculation','on');
%Start the Loop
for i=1:nbFiles
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
    gr=grains('BASIC',step,PartData,app); %calculate loops
    if isempty(PartData)
        PD='';
    else
        PD = grainListing(PartData,gr);
    end
    
    %Check if the spaceCellSystem exists, else calculate the sCS for the
    %desired type
    scFnm=fullfile(pathSc,"LoopsSpaceCellsfile"+step+".mat");
    if isfile(scFnm)
        sc=load(scFnm).sc;sv=0;
        if stnOF==1 && ~isequal(sc.Interval,interval)
            %if the interval of this calculation is different from hte
            %interval of a previos one, the calculation must be done again
            %and sc will be saved with the new valeus
            calcStn=1;  %calculate strain
            sv=1;       %prepar to save file
            sc.Interval=interval; %update interval value
        else
            calcStn=0; %else just read it
        end
    else
        sc = spaceCellSystem("Loops",step,gr,app,PD);sv=1;calcStn=1;
        sc.Interval=interval;
    end
    
    %check if loops exist
    if isempty(sc.Loops)
        CalcPanel(app,'','','','off');
        warndlg("ERROR:Loops Matrix empty on step :"+step);return;
    end
    
    %get size in a vector
    lOd=cat(1,sc.Loops.Order);   %loops order
    if D==2
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 2D %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Save step
        baseRes(i,1,1)=step;
        %count number of loops of same sizes, starting at loop 3 in 2D
        baseRes(i,6:7)=[mean(lOd),(2*gr.Z)/(gr.Z-2)];%calculate mean Loop
        %check partial data
        cord=gr.Coord;
        if ~isempty(PD)
            cord=gr.Coord(PD.GrainsRectangle,:);
        end
        %Loop to separate grains per size
        for j=1:4
            %find ID of lines with 3:+6 grains
            if j==4
                Ci=(lOd>=(j+D));
            else
                Ci=(lOd==(j+D));
            end
            baseRes(i,j+1,1)=sum(Ci);
            
            %VtkFiles
            if ~vtkOF; continue; end %if vtkOff we skip the vtk writing
            %copy only the correct lines
            drawLoops=cat(1,sc.Loops(Ci).Vertices);
            fnm=fullfile(pathVTK,fileName(j)+step+".vtk");
            if sum(Ci)==0
                %sometimes we will not have loops of a certain size, but we still
                %need to write downthe vtk file for it to show well on paraview.
                vtkwrite(fnm,'polydata','NoElements',...
                    zeros(size(cord,1),1),cord(:,2),cord(:,3),[0 0 0])
            else
                vtkwrite(fnm,'polydata',"TRIANGLE",...
                    zeros(size(cord,1),1),cord(:,2),cord(:,3),drawLoops);
            end
        end
    else 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 3D %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%% STRAIN %%%%%%
        %Strain - first calculation because it add Volume properties to
        %Clusters that can be used later
        if stnOF
            if i==1
                if calcStn || isempty(sc.Loops(1).Volume)
                    sc = clusterVolume(sc,gr);sv=1;
                end
            else
                if calcStn || isempty(sc.Loops(1).iStrain)
                    sc.Radius=gr.Radius;
                    sc.Displacements = gr.Coord-grOld.Coord;
                    sc=clusterStrain(sc,gr);sv=1;
                end
                stn=cat(3,sc.Loops.iStrain);
                %calc deviatoric strain for each cluster
                stn=sqrt(( (stn(2,2,:)-stn(1,1,:)).^2+...	%(sig1-sig2)^2 /2
                    (stn(3,3,:)-stn(1,1,:)).^2+...          %(sig2-sig3)^2 /2
                    (stn(1,1,:)-stn(3,3,:)).^2+...          %(sig3-sig1)^2 /2
                    6*stn(1,2,:).^2+...                     % 3*(sig12)^2
                    6*stn(1,3,:).^2+...                     % 3*(sig13)^2pcq
                    6*stn(2,3,:).^2 )/2);                   % 3*(sig23)^2
                stn=permute(stn,[3,2,1]); %turn 1x1xNl into Nlx1x1
            end
        end
        %%%%%% BASE VALS %%%%%%
        %Also get size in a vector for the Cluster
        lSz=cat(1,sc.Loops.Size);   %Cluster size
        nbC=cat(1,sc.Loops.nbCells);%nb of cells
        clV=cat(1,sc.Loops.Volume);%volume
        if isempty(clV)
            %in the case strain was not called in this istance or never for
            %these previous files Volume property needs to be recalculated
            sc = clusterVolume(sc,gr);
            clV=cat(1,sc.Loops.Volume);
        end
        %get maxsize and order
        baseRes(1,i,1) = max(lOd); %max Cluster order
        baseRes(1,i,2) = max(lSz); %max Cluster size
        %check partial data
        cord=gr.Coord;
        if ~isempty(PD)
            cord=gr.Coord(PD.GrainsRectangle,:);
        end
        %Count unique Clusters order
        %u are the unique valus and ic their position at lOd. We count
        %the diferent ics and save the amount of them in the positions
        %given by u/2 (since u is always even it works perfectly).
        lOd(lOd>maxCl)=maxCl;
        [u,~,ic]=unique(lOd);
        baseRes(u/2,i,1) =accumarray(ic,1);
        baseRes(2,i,1) = sc.Clt4.nbCells; %Cluster4 is no longer being counted in lSz
        %use the same count to take the nbCells into account
        baseRes(u/2,i,3) =accumarray(ic,nbC);
        baseRes(2,i,3) = sc.Clt4.nbCells; %nbC4=nbLp4
        %and Volume
        if doVol
            baseRes(u/2,i,4) =accumarray(ic,clV);
            baseRes(2,i,4)=sc.TotalVolume-sum(clV);
        end
        %Count unique clusters size - almost the same but no u/2
        lSz(lSz>maxCl)=maxCl;
        [u,~,ic]=unique(lSz);
        baseRes(u,i,2) =accumarray(ic,1);
        baseRes(2,i,2) = sc.Clt4.nbCells; %Cluster4 is no longer being counted in lSz
        
        %%%%%% ANISOTROPY %%%%%%
        if aniOF
            an=clusterAnisotropy(sc,gr);
            aniRes{i} =an;
        end
        %%%%%% VOID RATIO %%%%%%
        if vrOF
            %analyse VR
            %mean value per step
            if isempty(sc.Loops(end).VoidRatio)
                sc = clusterVR(sc,gr);sv=1;
            end
            clVR=[cat(1,sc.Loops.Order),cat(1,sc.Loops.VoidRatio)];
            [avrgVR,lpOd]= groupsummary(clVR(:,2), clVR(:,1), @mean);
            vrData(lpOd/2,i,1)=avrgVR;
            vrData(2,i,1)=mean(sc.Clt4.VoidRatio);
            %analyse Coordination nbr
            [avrgZ,lpOd]= groupsummary(cat(1,sc.Loops.Z), lOd, @mean);
            vrData(lpOd/2,i,2)=avrgZ;
            vrData(2,i,2)=mean(sc.Clt4.Z);
            %Save into cell array
            resVRZ{i,1}=[clVR,...
                cat(1,sc.Loops.Z), 1-cat(1,sc.Loops.Deformability)];
            resVRZ{i,2}=[ones(sc.Clt4.nbCells,1)*4, sc.Clt4.VoidRatio,...
                sc.Clt4.Z, 1-sc.Clt4.Z/3];
            %save total vr
            if isempty(sc.TotalVR)
                [~,vS,vV]=totalVoidRatio(sc,gr);
                vrTot(i)=vV/vS;
            else
                vrTot(i)=sc.TotalVR;
            end
            
        end
        %%%%%% DEFORMABILITY %%%%%%
        if defOF
            %deformability in singleLoops is calculated by Nclosed/Ntotal,
            %but Nopen/Ntotal has a better meaning as a Perfect loop 4 will
            %have a deformability of 0 (non deformable). Thus the def is
            %recalculated as (lp.Deformability).
            % def contain : [Order, Size, Def, Nedges, Nopen]
            def=cat(1,sc.Loops.Deformability);          %deformability
            nbCe=cellfun(@height, {sc.Loops.ClEdges})';  %nb closed edges
            nTot=nbCe./(1-def);                         %total edges
            defRes{i,1}=[lOd,lSz,def,nTot,nTot-nbCe];   %clusters O>4
            [nbD4,d4]=groupcounts(1-cat(1,sc.Clt4.Z)/3);
            defRes{i,2}=[d4,nbD4];                      %clusters O===4
        end
        %%%%%% STRESS %%%%%%
        if strOF
            if isempty(sc.Loops(1).Stress)
                sc=clusterStress(sc,step,app,gr);sv=1;
            end
            %Order - q - p - vol
            lS=[cat(1,sc.Loops.Order);ones(sc.Clt4.nbCells,1)*4];
            strRes{i}=[lS [cat(1,sc.Loops.Stress) cat(1,sc.Loops.Volume);
                sc.Clt4.Stress,sc.Clt4.Volume]];
        end
        %%%%%% TRANSFORMATION %%%%%%
        if trOF
            %identify the cells from the sc that are clusters cat 4 or 6
            %and save in a new 'cl' object
            cl=clIdentifier(sc);
            
            %if strain was not called as a result, calculate the per
            %cell strain to analyse with transformation
            if isempty(sc.CellStn)
                if i>1;sc = pcStrain(sc,gr,grOld);end
            elseif isempty(sc.CellVol) && i>1
                sc.CellVol = perCellVolume(sc,gr);
            end
            %Analysis can only be made between two steps, so pass the first
            if i>1
                trRes{i,1}=clusterTransformation(sc,cl,gr,clOld);
                trRes{i,2}=pagemtimes(sc.CellStn,sc.CellVol);
            end
            %Save values for next step
            clOld=cl;%grOld=gr;
        end
        %%%%%% AVE GRAIN CATEGORY %%%%%%
        if aveOF
            gr = aveCluster(gr,sc);
            fnm=fullfile(pathAve,"Clusters_"+step+".vtk");
            vtkwrite(fnm,'unstructured_grid',...
                gr.Coord(:,1),gr.Coord(:,2),gr.Coord(:,3),...
                'SCALARS','Radius',gr.Radius,...
                'SCALARS','Ave_Cat',gr.AveCluster,...
                'PRECISION',10);
            aveCl(:,i)=gr.AveCluster;
        end
        %%%%%% BASE VTK WRITING %%%%%%
        if vtkOF
            fnm=fullfile(pathVTK,"Clusters_"+step+".vtk");
            drawClst=cat(1,sc.Loops.Vertices);
            ord=cat(1,sc.Loops.Order);
            sz=cat(1,sc.Loops.Size);
            nbscVal=(3+aniOF+vrOF+strOF*2+stnOF);
            scalVal=cell(nbscVal*2,1);   
            
            scalVal(1:6)={"Order",repelem(ord,ord),...
                "OrderRat",repelem(ord/max(ord),ord),...
                "Size",repelem(sz,ord)};
            scalp=7;
            if aniOF
                %get only the elevation
                anV=abs(an.Angles(:,2,1)); %negative and pos angles are the sme
                scalVal(scalp:scalp+1)={"Elevation",repelem(anV(:,1),ord)};
                scalp=scalp+2;
            end
            if vrOF
                scalVal(scalp:scalp+1)={"Void_Ratio",...
                    repelem(cat(1,sc.Loops.VoidRatio),ord)};
                scalp=scalp+2;
            end
            if strOF
                st=cat(1,sc.Loops.Stress);
                scalVal(scalp:scalp+1)={"Stress_p",repelem(st(:,1),ord),...
                    "Stress_q",repelem(st(:,2),ord)};
                scalp=scalp+4;
            end
            if stnOF && i>1
                scalVal(scalp:scalp+1)={"Strain_Dev",repelem(stn,ord)};
            elseif stnOF
                scalVal(scalp:scalp+1)={"Strain_Dev",...
                    repelem(zeros(size(ord)),ord)};
            end
            vtkwrite(fnm,'polydata',"TRIANGLE",...
                cord(:,1),cord(:,2),cord(:,3),drawClst,...
                'facescalar',nbscVal,scalVal{:});
        end
        
        if trOF || strOF;grOld=gr;end
    end%end if==2
    
    %save file if it was created here
    if sv
        sc=purge(sc,'LOOP');
        save(scFnm,'sc','-v7.3');
    end
    
end
CalcPanel(app,i+1,nbFiles,'','off');

%% Part 3 - Save and plot
%Calculate stress and strain values
strain = extStrains(app.TrialData,stepArray,N1,app,'allCalc');
stress = extStress(app.TrialData,stepArray,app);
consoStrain = extStrains(app.TrialData,app.ConsoStep,N1,app);

%In most cases the curbs x axis will be the axial strain. However for the
%q-cst tests it must be replaced for the mean pressure value.
xAxis=strain(:,D); 
if app.SimType==3;xAxis=stress(:,end);end %if Qcst axStr=p;

%Create the final vector adding the strain column just after the steps one
if D==2
    resBs=[baseRes(:,1)*app.TimeStep.*ones(1,1,2),...
        xAxis.*ones(1,1,2),baseRes(:,2:end,:)];
else
    %remove excess values - Order
    f1=find(sum(baseRes(:,:,1),2),1,'last');
    resO=[2*(1:f1)' baseRes(1:f1,:,1)];
    resC=[2*(1:f1)' baseRes(1:f1,:,3)]; %do the same for cells
    if doVol;resV=[2*(1:f1)' baseRes(1:f1,:,4)];end %and volume
    %remove excess values - Size
    f2=find(sum(baseRes(:,:,2),2),1,'last');
    resS=[(3:f2)' baseRes([1:2,5:f2],:,2)]; %lines 3 and 4 are empty and need to be removed
    %Create obj containing values
    resBs.Order=resO;
    resBs.Size=resS;
    resBs.Cells=resC;
    resBs.Volume=resV;
    resBs.Strain=strain(:,D);
    resBs.Pressure=stress(:,end); 
    resBs.SimType=app.SimType;
end

%Base loops plot
pD=plotData("Normal",resBs,app,prefix,consoStrain(D));
fnm=fullfile(MakePath(app,'LOOPS'),...
    "Loops-Base"+N1+"to"+N2+"int"+interval+".mat");
save(fnm,'pD','-v7.3');  
lPlotter(pD,app);

if vtkOF
    fnm=fullfile(pathVTK,"vtkVariables.txt");
    StepEzEvQP=[stepArray,strain(:,[D,D+1]),stress(:,end-1:end)];
    vtkLog(app,'Loops',fnm,StepEzEvQP)
end

%plot anisotropy data
if aniOF
    aRes.Anisotropy=aniRes;
    aRes.Strain=strain(:,D);
    aRes.Pressure=stress(:,end); 
    aRes.QEv=[stress(:,end-1),strain(:,end-1)];
    pD=plotData("Normal",aRes,app,prefix);
    fnm=fullfile(MakePath(app,'LOOPA'),...
        "Loops-Anisotropy"+N1+"to"+N2+"int"+interval+".mat");
    save(fnm,'pD','-v7.3');  
    clAniPlotter(pD,app)
end
%plot Deformability data
if defOF
    dRes.SOD=defRes;          %Order - Size - Deformability
    dRes.Strain=strain(:,D);
    dRes.Pressure=stress(:,end); 
    pD=plotData("Normal",dRes,app,prefix);
    fnm=fullfile(MakePath(app,'LOOPD'),...
        "Loops-Dfblty"+N1+"to"+N2+"int"+interval+".mat");
    save(fnm,'pD','-v7.3');   %save
    clDefPlotter(pD,app)
end
%plot Void Ratio data
if vrOF
    %find last value of loopsize
    f=find(sum(sum(vrData,2),3),1,'last');
    vrData=vrData(1:f,:,:);
    %add cluster order indicator to first line after rotating it
    vrData=cat(2,(2:2:2*(size(vrData,1)))'.*ones(1,1,2),vrData);

    %Save resVRZ file
    %Create an object
    resBs.Mean=vrData;
    resBs.ClusterVRZ=resVRZ;
    resBs.TotalVR=vrTot;
    pD=plotData("Normal",resBs,app,prefix);
    fnm=fullfile(MakePath(app,'LOOPVR'),...
        "Loops-VR"+N1+"to"+N2+"int"+interval+".mat");
    save(fnm,'pD','-v7.3');
    %plot
    clVRPlotter(pD,app)
end
%stress plot
if strOF
    stRes.Stress=strRes;
    stRes.Strain=xAxis;
    pD=plotData("FastCreation",stRes,app,prefix);
    fnm=fullfile(MakePath(app,'LOOPST'),...
        "Loops-Stress"+N1+"to"+N2+"int"+interval+".mat");
    save(fnm,'pD','-v7.3');  %save
    clStrPlotter(pD,app)
end
%ClusterTransformation
if trOF
    ctRes.ClTransf=trRes;       %cl transf results
    ctRes.Stress=strRes;        %axial stress per step
    ctRes.Strain=xAxis;         %axial strain per step
    ctRes.NbClst=resO;          %nb cluster per order
    ctRes.NbCell=resC;          %nb cells per order
    ctRes.SimType=app.SimType;  %sim type
    pD=plotData("Normal",ctRes,app,prefix);
    fnm=fullfile(MakePath(app,'LOOPCT'),...
        "Loops-Transf"+N1+"to"+N2+"int"+interval+".mat");
    save(fnm,'pD','-v7.3');    %save
    clTransfPlotter(pD,app)

end

if aveOF
    %Add first colum as the nb of grain and first line as step array
    aveCL=[ [0 stepArray'] ; [(1:app.NbGrains)'  aveCl]];
    fnm=fullfile(pathAve,"Average_Cluste_pGrain_"+N1+"to"+N2+"int"+interval+".txt");
    fid = fopen(fnm, 'w');
    fprintf(fid, '## Average Cluster per Grain ##\n');
    fprintf(fid, 'Step_Init|Step_Final|Interval\n');
    fprintf(fid, '%d|%d|%d\n',N1,N2,interval);
    fprintf(fid, '\n');
    fprintf(fid,['V%d' repmat('|V%d',1,nbFiles) '\n'] ,0:nbFiles);
    fprintf(fid,['%d' repmat('|%d',1,nbFiles) '\n'] ,aveCL');
    fclose(fid);
end
%Copy forcechains.vtk files, interesting only on 2D cases
if ~vtkOF || D==3; return;end
%copy grains VTK files inside the directory for better visualization
OldF=cd;
cd(app.VTKpathEF.Value);
for i=1:nbFiles
    step=min(N1+interval*(i-1),N2);
    fnm="partcles_"+step+".vtk";
    fnm2="forcechain"+step+".vtk";
    copyfile(fnm,pathVTK);
    copyfile(fnm2,pathVTK);
end
cd(OldF);
end

%load functions
function LoopsLoad(app,mode)
%LOOPSLOAD load already calculated Loops
%	Read a file previously calculated by this function and prepare to plot
%	it using the LoopsPlotter function
%
%   Variable PartData in this case contain a value to separate the three
%   diferent kinds of load. Nothing to do with the actual partialData
%   Objects (it was just easier to do)

%Load file
switch upper(mode)
    case 'NORMAL'
        %Normal LOAD
        %{
        %OLD LOAD PROCEDURE FROM TXT FILES
        pD=FileLoader(app,'LOOPS');
        
        if isempty(pD);return;end
        if size(pD,2)>1
            pD(1).Prefix='multi';
        end
        for i=1:size(pD,2)
            file=pD(i).Results;
            f=find(sum(isnan(file),2));
            if numel(f)<2
                warndlg(['Wrong file or not up to date. Please run'...
                    ' the calculation again or choose another file']);
                return;
            end
            res.Order=file(1:f(1)-1,:);
            res.Size=file(f(1)+1:f(2)-1,:);
            if numel(f)==2
                res.Cells=file(f(2)+1:end,:);
                res.Volume='';
            else
                res.Cells=file(f(2)+1:f(3)-1,:);
                res.Volume=file(f(3)+1:end,:);
            end
            pD(i).Results=res;
        end
        %}
        [fnm,fPath]=MatLoader('LOOPS',app);
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
        lPlotter(pD,app)
        fprintf('Calculations finished \n')
    case 'VOIDRATIO'
        %VoidRatio Plot
        [fnm,fPath]=MatLoader('LOOPVR',app);
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
        fprintf('File loaded. Starting calculations. \n')
        fprintf('It may take a few minutes. \n')
        pD(1).Prefix="Load";
        clVRPlotter(pD,app)
        fprintf('Calculations finished \n')
    case 'DEFORMABILITY'
        [fnm,fPath]=MatLoader('LOOPD',app);
        if fPath==0;warndlg('No files were chosen, Try again.');return;end
        %create plot data objects
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
            pD=load(fPath+""+fnm).pD;
        end
        fprintf('File loaded. Starting calculations \n')
        fprintf('It may take a few minutes \n')
        clDefPlotter(pD,app)
        fprintf('Calculations finished \n')
    case 'STRESS'
        [fnm,fPath]=MatLoader('LOOPST',app,'off');
        if fPath==0;warndlg('No files were chosen, Try again.');return;end
        stRes=load(fullfile(fPath,fnm)).pD;
        clStrPlotter(plotData("FastCreation",stRes,''),app)
    case 'CLTRF'
        [fnm,fPath]=MatLoader('LOOPCT',app,'off');
        if fPath==0;warndlg('No files were chosen, Try again.');return;end
        try pD=load(fullfile(fPath,fnm)).pD;
        catch
            ctRes=load(fullfile(fPath,fnm)).ctRes;
            pD=plotData("Normal",ctRes,app,'Load');
            save(fullfile(fPath,fnm),'pD','-v7.3');    %save
        end
        clTransfPlotter(pD,app)
    case 'ANISOTROPY'
        [fnm,fPath]=MatLoader('LOOPA',app,'on');
        if fPath==0;warndlg('No files were chosen, Try again.');return;end
        if iscell(fnm) && numel(fnm)>1
            for i=1:numel(fnm)
                n=fnm{i};
                pD=load(fullfile(fPath,n)).pD;
                pD.Prefix=n(1:strfind(n,'.mat')-1);
                clAniPlotter(pD,app)
            end
        else
            pD=load(fullfile(fPath,fnm)).pD;
            pD.Prefix="Load";
            clAniPlotter(pD,app)
        end
end %end Switch
end %end load

%plot functions
function lPlotter(pD,app)
%LOOPSPLOTTER plot the evolution of loops using SubPlotter applicatioon
%   2D - Prepare the values calculated by the main function plot the
%   evolution of the Loops in the SubPlotter aplication
%   3D - plot directly the important Curbes
%check if pD.Results is a structure or not
if ~isa(pD(1).Results,'struct')
    LGraph=pD(1).Results;
    %prepare to plot
    grpLabel=["Loop3","Loop4","Loop5","Loop+6"];
    SubPlotter(app, 'LOOPS', LGraph, grpLabel,pD.ConsoStrain, '');
    return
end

if numel(pD)>1
    pltIndv=0;
    %Ask user if, for multi-load, plot results of individual cluster files
    answer = questdlg('Plot cluster results for each file separately?', ...
    	'Extra figures', ...
    	'Yes','No');
    % Handle response
    switch answer
        case 'Yes'
            pltIndv=1;
        case 'No'
            pltIndv=0;
    end
    
end
% N1=app.N1EF.Value;
% N2=app.N2EF.Value;
% interval=app.CalcInt.Value;
% TD  = inflectionPoints(app.TrialData,app);
% pD.InfPts=TD.InfPts;
% fnm=MakePath(app,'LOOPS')+"Loops-Base"+N1+"to"+N2+"int"+interval+".mat";
% save(fnm,'pD','-v7.3');

%matlab default colors to mantain every class in the same clor
if numel(pD)<8
    C = app.PlotColors;
else
    C = graphClrCode(size(pD,2));%plot colorcode
end

%Check title and legends option
if app.TitlesCB.Value;tit=1;else;tit=0;end
if app.LegendsCB.Value;leg=1;else;leg=0;end
set(0,'defaultAxesFontSize',app.FontSizeEF.Value)

%If any of the loaded files is a Qcst simultaion, the plot will be done
%using pressure as x-axis
xlab='Axial strain';
x='e';
for i=1:numel(pD)
    if pD(i).SimType==3
        xlab='Mean pressure (kPA)';
        x='p';
    end
end
%prepare paths
pathA=MakePath(app,'LOOPL');
pathB=MakePath(app,'LOOPS');
png=".png";
%fig=".fig";

%%%%%%%%%%%%%%%% Multifiles plot %%%%%%%%%%%%%%%%

if numel(pD)>1
    nb=18;
    fA(nb)=figure;axA(nb)=axes(fA(nb));hold(axA(nb),'on');
    for i=1:(nb-1)
        fA(i)=figure;axA(i)=axes(fA(i));hold(axA(i),'on');
    end
    data=zeros(numel(pD),4,4);
    for i=1:numel(pD)
        if strcmpi(x,'p')
            xaxis= pD(i).Results.Pressure;
            [~,pkPos]=min(abs(xaxis-pD(i).InfPts.p(end)));
        else
            xaxis= pD(i).Results.Strain;
            [~,pkPos]=min(abs(xaxis-pD(i).InfPts.ez(end)));
        end
        opts={'Color',C(i,:)};
        for j=1:numel(pD(i).InfPts.q)
            if strcmpi(x,'p')
                opts=[opts,{'Pointx',pD(i).InfPts.p(j)}]; %#ok<AGROW>
            else
                opts=[opts,{'Pointx',pD(i).InfPts.ez(j)}]; %#ok<AGROW>
            end
        end
        rO=(pD(i).Results.Order(1:end,2:end));totO=sum(rO,1); %Orders
        rC=(pD(i).Results.Cells(1:end,2:end));totC=sum(rC,1); %Cells
        rC=rC./(ones(size(rC,1),1)*totC);
        %Order
        j=1;
        %Numbers of Orders
        plotMark(app,axA(j),xaxis,rO(2,:),opts{:});j=j+1;               %Order 4=f(Ez)
        plotMark(app,axA(j),xaxis,rO(3,:),opts{:});j=j+1;               %Order 6=f(Ez)
        plotMark(app,axA(j),xaxis,sum(rO(4:10,:),1),opts{:});j=j+1;     %Order 8-20=f(Ez)
        plotMark(app,axA(j),xaxis,sum(rO(11:end,:),1),opts{:});j=j+1;   %Order 22=f(Ez)
        %ratio of Orders
        plotMark(app,axA(j),xaxis,rO(2,:)./totO,opts{:});j=j+1;             %Order 4=f(Ez)
        plotMark(app,axA(j),xaxis,rO(3,:)./totO,opts{:});j=j+1;             %Order 6=f(Ez)
        plotMark(app,axA(j),xaxis,sum(rO(4:10,:)./totO,1),opts{:});j=j+1;   %Order 8-20=f(Ez)
        plotMark(app,axA(j),xaxis,sum(rO(11:end,:)./totO,1),opts{:});j=j+1; %Order 22=f(Ez)
        %MaxOrder=f(Ez)
        plotMark(app,axA(j),xaxis,rO(1,:),opts{:});j=j+1; %MaxOrder=f(Ez)
        %ratio of Cells
        plotMark(app,axA(j),xaxis,rC(2,:),opts{:});j=j+1;             %Ratio Order 4=f(Ez)
        plotMark(app,axA(j),xaxis,rC(3,:),opts{:});j=j+1;             %Ratio Order 6=f(Ez)
        plotMark(app,axA(j),xaxis,sum(rC(4:19,:),1),opts{:});j=j+1;   %Ratio Order 8-20=f(Ez)
        plotMark(app,axA(j),xaxis,sum(rC(11:end,:),1),opts{:});       %Ratio Order 22=f(Ez)
        %Size
        rS=(pD(i).Results.Size(1:end,2:end));
        j=j+1;
        plotMark(app,axA(j),xaxis,rS(1,:),opts{:});

        %Create a mark a max q/p and final state for each specimen for each
        %cluster category
        j=j+1;
        %Order 4
        data(i,:,1)=[xaxis(pkPos),rC(2,pkPos+1),xaxis(end),rC(2,end)]; 
        scatter(axA(j),data(i,[1,3],1),data(i,[2,4],1) );j=j+1;
        %Order 6
        data(i,:,2)=[xaxis(pkPos),rC(3,pkPos+1),xaxis(end),rC(3,end)];
        scatter(axA(j),data(i,[1,3],2),data(i,[2,4],2) );j=j+1;
        %Order 8
        data(i,:,3)=[xaxis(pkPos),sum(rC(4:10,pkPos+1),1),xaxis(end),sum(rC(4:10,end),1)];
        scatter(axA(j),data(i,[1,3],3),data(i,[2,4],3) );j=j+1;
        %Order 22
        data(i,:,4)=[xaxis(pkPos),sum(rC(11:end,pkPos+1),1),xaxis(end),sum(rC(11:end,end),1)];
        scatter(axA(j),data(i,[1,3],4),data(i,[2,4],4) )

    end

    %create log law with data matrix
    for ii=1:4
        if ii==1
            mp=[0.5,0,0,0,0];            %Cl4
            side=["left","left"];
        elseif ii==2
            mp=[0.5,0,-0.0005,0,0];           %Cl6
            side=["left","left"];
        elseif ii==3
            mp=[0.6,0,0,0,0.005];  %Cl8
            side=["left","left"];
        else
            mp=[0.3,10,0.005,10,0];%Cl22
            side=["left","left"];
        end
        %Failure state
        [ft,gof]=fit(data(:,1,ii),data(:,2,ii),fittype('a*log10(x)+b')); 
        x=[min(data(:,1,ii)),max(data(:,1,ii))];y=[min(data(:,2,ii)),max(data(:,2,ii))];
        axes(axA(nb-4+ii))
        p=plot(ft,'--k',x,y);
        p(2).LineWidth=1;
        delete(axA(nb-4+ii).Children(2))
        xx=(x(1)+x(2))*mp(1);
        text(axA(nb-4+ii),xx+mp(2),(ft.a*log10(xx)+ft.b)+mp(3),...
            {'',"     Failure state, R2 = "+ sprintf("%02.f ",gof.rsquare*100)+ "%"},...
            'HorizontalAlignment',side(1),'Color','black','FontSize',12)

        %Final state
        [ft,gof]=fit(data(:,3,ii),data(:,4,ii),fittype('a*log10(x)+b'));
        x=[min(data(:,3,ii)),max(data(:,3,ii))];y=[min(data(:,4,ii)),max(data(:,4,ii))];
        axes(axA(nb-4+ii))
        p=plot(ft,':k',x,y);
        p(2).LineWidth=1;
        delete(axA(nb-4+ii).Children(2))
        xx=(x(1)+x(2))*mp(1);
        text(axA(nb-4+ii),xx+mp(4),(ft.a*log10(xx)+ft.b)+mp(5),...
            {'',"     Final state, R2 = " + sprintf("%02.f",gof.rsquare*100) + "%"},...
            'HorizontalAlignment',side(2),'Color','black','FontSize',12)

    end


    leglocation=["northeast";"northeast";"southeast";
        "southeast";"northeast";"northeast";"southeast";
        "southeast";"northwest";
        "northeast";"northeast";"southeast";"southeast";
        "northwest";
        "northeast";"northeast";"southeast";"southeast"];
    titles=["Number of Clusters Order 4";
        "Number of Clusters Order 6";
        "Number of Clusters Order 8-20";
        "Number of Clusters Order 22+";
        "Ratio of Clusters Order 4";
        "Ratio of Clusters Order 6";
        "Ratio of Clusters Order 8-20";
        "Ratio of Clusters Order 22+";
        "Maximal Cluster Order";
        "Ratio of Cells Order 4";
        "Ratio of Cells Order 6";
        "Ratio of Cells Order 8-20";
        "Ratio of Cells Order 22+";
        "Maximal Cluster Size";
        "Ratio of Cells Order 4";
        "Ratio of Cells Order 6";
        "Ratio of Cells Order 8-20";
        "Ratio of Cells Order 22+"];
    ylab=["Number of Clusters";"Number of Clusters";
        "Number of Clusters";"Number of Clusters";
        "Ratio of Clusters";"Ratio of Clusters";
        "Ratio of Clusters";"Ratio of Clusters";
        "Max Order";
        "Ratio of Cells";"Ratio of Cells";
        "Ratio of Cells";"Ratio of Cells";
        "Max Size";        
        "Ratio of Cells";"Ratio of Cells";
        "Ratio of Cells";"Ratio of Cells"];
    fnm=["Cluster_Order_Nb_4";"Cluster_Order_Nb_6";
        "Cluster_Order_Nb_8";"Cluster_Order_Nb_22";
        "Cluster_Order_Pct_4";"Cluster_Order_Pct_6";
        "Cluster_Order_Pct_8";"Cluster_Order_Pct_22";
        "Cluster_Order_Max";
        "Cluster_Cell_4";"Cluster_Cell_6";
        "Cluster_Cell_8";"Cluster_Cell_22";
        "Cluster_Size_Max";
        "Cluster_Law_4";"Cluster_Law_6";
        "Cluster_Law_8";"Cluster_Law_22"];
    for i=1:nb
        if leg && i<nb-4;legend(axA(i),pD.FileName,'location',leglocation(i));end
        if tit;title(axA(i),titles(i));end
        ylabel(axA(i),ylab(i))
        xlabel(axA(i),xlab)
        if strcmpi(x,'p');axA(i).XLim(1)=0;end
        %         if i==7
        %             axA(i).YLim=[0.78,0.94];
        %         elseif i==8
        %             axA(i).YLim=[0.05,0.12];
        %         elseif i==9
        %             axA(i).YLim=[0.01,0.11];
        %         elseif i==10
        %             axA(i).YLim=[0,0.014];
        %         end
        saveas(fA(i),fullfile(pathA,pD(1).Prefix+fnm(i)+png));
        if i==1
            if (fA(i).Position(3:4)~=app.Figsize)
                app.Figsize=fA(i).Position(3:4);
            end
        end
    end
end

%%%%%%%%%%%%%%%% Single files plot (also run on multifiles, for each of them) %%%%%%%%%%%%%%%%
if numel(pD)==1 || pltIndv
    %FIGURE CREATION
    %figures for Order and size
    nOrd=3; %number of figures
    nb=nOrd*size(pD,2);
    fB(nb)=figure;axB(nb)=axes(fB(nb));hold(axB(nb),'on'); %order
    %fC(nb)=figure;axC(nb)=axes(fC(nb));hold(axC(nb),'on');
    for i=1:(nb-1)
        fB(i)=figure;axB(i)=axes(fB(i));hold(axB(i),'on'); %#ok<*LAXES>
        %fC(i)=figure;axC(i)=axes(fC(i));hold(axC(i),'on');
    end
    %figure Percentage of Cl Order, cells and volume
    %%check if volume exists in all pD objects
    doVol=1;nPct=6;
    for j=1:size(pD,2)
        if isempty(pD(j).Results.Volume);doVol=0;nPct=4;end
    end
    nb=nPct*size(pD,2);
    fD(nb)=figure;axD(nb)=axes(fD(nb));hold(axD(nb),'on'); %order
    for i=1:(nb-1)
        fD(i)=figure;axD(i)=axes(fD(i));hold(axD(i),'on'); %#ok<*LAXES>
    end
    %figures comparaision Size-Order, nb of Cells
    nCmp=2;
    nb=nCmp*size(pD,2);
    fE(nb)=figure;axE(nb)=axes(fE(nb));hold(axE(nb),'on');
    for i=1:(nb-1)
        fE(i)=figure;axE(i)=axes(fE(i));hold(axE(i),'on');
    end

    %PLOTTER
    for j=1:size(pD,2)
        %Either 1 file loaded or multifile
        if size(pD,2)==1
            pfx=pD.Prefix;
            path=pathB;
        else
            pfx=pD(j).FileName;
            path=pathA;
        end
        if strcmpi(x,'p')
            xaxis= pD(j).Results.Pressure;
        else
            xaxis= pD(j).Results.Strain;
        end
        opts={'Color',[0 0 0]}; %color will be changed later
        for i=1:numel(pD(j).InfPts.q)
            if strcmpi(x,'p')
                opts=[opts,{'Pointx',pD(j).InfPts.p(i)}]; %#ok<AGROW>
            else
                opts=[opts,{'Pointx',pD(j).InfPts.ez(i)}]; %#ok<AGROW>
            end
        end
        %get all the different results : order, size, nb of cells and strain.
        rO=(pD(j).Results.Order);      %order (vector)
        totO=sum(rO(2:end,2:end),1);    %number of clusters (scalar)
        rS=(pD(j).Results.Size);       %size (vector)
        rC=(pD(j).Results.Cells);      %cells (vector)
        totC=sum(rC(2:end,2:end),1);    %number of cells (scalar)
        if doVol
            rV=(pD(j).Results.Volume);     %volume (vector)
            totV=sum(rV(2:end,2:end),1);    %total volume (scalar)
        end
        for i=1:4
            switch i
                case {1,2}
                    pO=rO(i+1,2:end);
                    pS=rS(i+1,2:end);
                    pC=rC(i+1,2:end);
                    if doVol;pV=rV(i+1,2:end);end
                case 3
                    pO=sum(rO(4:10,2:end),1);
                    pS=sum(rS(4:10,2:end),1);
                    pC=sum(rC(4:10,2:end),1);
                    if doVol;pV=sum(rV(4:10,2:end),1);end
                case 4
                    pO=sum(rO(11:end,2:end),1);
                    pS=sum(rS(11:end,2:end),1);
                    pC=sum(rC(11:end,2:end),1);
                    if doVol;pV=sum(rV(11:end,2:end),1);end
            end
            opts{2}=C(i,:);
            %CLUSTER NUMBER
            plotMark(app,axB((j-1)*nOrd+1),xaxis,pO,opts{:}); %Order =f(Ez)
            plotMark(app,axD((j-1)*nPct+1),xaxis,pO./totO,opts{:}); %Order/totO =f(Ez)
            plotMark(app,axD((j-1)*nPct+3),xaxis,pC./totC,opts{:}); %Cells/totC =f(Ez)
            if doVol %Vol/totVol =f(Ez)
                plotMark(app,axD((j-1)*nPct+5),xaxis,pV./totV,opts{:});
            end
            %CLUSTER NUMBER NO4
            if i~=1
                plotMark(app,axB((j-1)*nOrd+2),xaxis,pO,opts{:});
                plotMark(app,axD((j-1)*nPct+2),xaxis,pO./totO,opts{:});
                plotMark(app,axD((j-1)*nPct+4),xaxis,pC./totC,opts{:});
                if doVol
                    plotMark(app,axD((j-1)*nPct+6),xaxis,pV./totV,opts{:});
                end
            end
            %CLUSTER ORDER Vs SIZE
            plotMark(app,axE((j-1)*nCmp+1),xaxis,pS-pO,opts{:}); %Cluster size-Order =f(Ez)
        end
        opts{2}=C(1,:);
        plotMark(app,axB((j-1)*nOrd+3),xaxis,rO(1,2:end),opts{:}); %maximal order
        %     if size(pD,2)==1
        %         pct=cell(1,2);
        %         pct{1,1}=percentile(rO,p);
        %         %pct{1,2}=percentile(rS,p);
        %     end
        %     plotMark(app,axB((j-1)*nOrd+4),xaxis,pct{j,1},'Color',C(i,:)); %percentile order
        plotMark(app,axE((j-1)*nCmp+2),xaxis,totC,opts{:}); %nb of cells

        if leg
            %Order
            legend(axB((j-1)*nOrd+1),"Small","Submedium","Medium",...
                "Large",'location','northeast')
            legend(axB((j-1)*nOrd+2),"Submedium","Medium",...
                "Large",'location','northeast')
            %Size
            %         legend(axC((j-1)*npl+1),"Size 4","Size 5","Size 6-12",...
            %             "Size 13+",'location','eastoutside')
            %         legend(axC((j-1)*npl+2),"Size 5","Size 6-12",...
            %             "Size 13+",'location','eastoutside')
            %Order ratio
            legend(axD((j-1)*nPct+1),"Small","Submedium","Medium",...
                "Large",'location','east')
            legend(axD((j-1)*nPct+2),"Submedium","Medium",...
                "Large",'location','east')
            %Cell ratio
            legend(axD((j-1)*nPct+3),"Small","Submedium","Medium",...
                "Large",'location','northeast')
            legend(axD((j-1)*nPct+4),"Submedium","Medium",...
                "Large",'location','southeast')
            if doVol
                %Volume ratio
                legend(axD((j-1)*nPct+5),"Small","Submedium","Medium",...
                    "Large",'location','northeast')
                legend(axD((j-1)*nPct+6),"Submedium","Medium",...
                    "Large",'location','southeast')
            end
            %size/order
            legend(axE((j-1)*nCmp+1),"S4/O4","S5/O6","S6-12/O8-20",...
                "S13+/O22+",'location','northeast')
        end

        %titles
        if tit
            title(axB((j-1)*nOrd+1),'Number per Order Category')
            title(axB((j-1)*nOrd+2),'Number per Order Category')
            title(axB((j-1)*nOrd+3),'Max Cluster Order')
            %       title(axB((j-1)*nOrd+4),['Cluster Order ' num2str(p) ' Percentile'])
            %       title(axC((j-1)*npl+1),'Number per Size Category')
            %       title(axC((j-1)*npl+2),'Number per Size Category')
            %       title(axC((j-1)*npl+3),'Max Cluster Size')
            %       title(axC((j-1)*npl+4),['Cluster Size ' num2str(p) ' Percentile'])
            title(axD((j-1)*nPct+1),'Number Ratio per Order Category')
            title(axD((j-1)*nPct+2),'Number Ratio per Order Category')
            title(axD((j-1)*nPct+3),'Cells Ratio per Order Category')
            title(axD((j-1)*nPct+4),'Cells Ratio per Order Category')
            if doVol
                title(axD((j-1)*nPct+5),'Volume Ratio per Order Category')
                title(axD((j-1)*nPct+6),'Volume Ratio per Order Category')
            end
            title(axE((j-1)*nCmp+1),'Evolutio of Size-Order')
            title(axE((j-1)*nCmp+2),'Evolutio of nb Cells')
        end
        %plot1 nb ORDER=f(Ez)
        ylabel(axB((j-1)*nOrd+1),'Number of Cluster')
        xlabel(axB((j-1)*nOrd+1),xlab)
        if size(pD,2)>1
            axB((j-1)*nOrd+1).YLim(2)=axA(1).YLim(2);
            yticks(axB((j-1)*nOrd+1),'auto');
            fB((j-1)*nOrd+1).Position(3:4)=app.Figsize;
        end
        if strcmpi(x,'p');axB((j-1)*nOrd+1).XLim(1)=0;end
        fnm=pfx+"Cl_Order";
        saveas(fB((j-1)*nOrd+1),fullfile(path,fnm+png));

        %plot2 nb ORDER=f(Ez) no Clusters 4
        ylabel(axB((j-1)*nOrd+2),'Number of Clusters')
        xlabel(axB((j-1)*nOrd+2),xlab)
        if size(pD,2)>1
            axB((j-1)*nOrd+2).YLim(2)=max(axA(2).YLim(2),axA(3).YLim(2));
            yticks(axB((j-1)*nOrd+2),'auto');
            fB((j-1)*nOrd+2).Position(3:4)=app.Figsize;
        end
        if strcmpi(x,'p');axB((j-1)*nOrd+2).XLim(1)=0;end
        fnm=pfx+"Cl_Order_N4";
        saveas(fB((j-1)*nOrd+2),fullfile(path,fnm+png));

        %plot3 ORDER Max
        ylabel(axB((j-1)*nOrd+3),'Max Order')
        xlabel(axB((j-1)*nOrd+3),xlab)
        if size(pD,2)>1
            axB((j-1)*nOrd+3).YLim(2)=axA(9).YLim(2);
            yticks(axB((j-1)*nOrd+3),'auto');
            fB((j-1)*nOrd+3).Position(3:4)=app.Figsize;
        end
        if strcmpi(x,'p');axB((j-1)*nOrd+3).XLim(1)=0;end
        fnm=pfx+"Cl_Order_Max";
        saveas(fB((j-1)*nOrd+3),fullfile(path,fnm+png));

        %     %plot4 ORDER percentile
        %     ylabel(axB((j-1)*nOrd+4),'Order 99 percentile')
        %     xlabel(axB((j-1)*nOrd+4),xlab)
        %     if size(pD,2)>1
        %         axB((j-1)*nOrd+4).YLim(2)=axA(6).YLim(2);
        %         yticks(axB((j-1)*nOrd+4),'auto');
        %         fB((j-1)*nOrd+4).Position(3:4)=app.Figsize;
        %     end
        %     if strcmpi(x,'p');axB((j-1)*nOrd+1).XLim(1)=0;end
        %     fnm=pfx+"Cl_Order_Pctile";
        %     saveas(fB((j-1)*nOrd+4),fullfile(path,fnm+png));
        %     %saveas(fB((j-1)*npl+4),path+fnm+fig);

        %     %plot1 SIZE Clusters=f(Ez)
        %     ylabel(axC((j-1)*npl+1),'Number of Clusters')
        %     xlabel(axC((j-1)*npl+1),xlab)
        %     if size(pD,2)>1
        %         axC((j-1)*npl+1).YLim(2)=axA(6).YLim(2);
        %         yticks(axC((j-1)*npl+1),'auto');
        %         fC((j-1)*npl+1).Position(3:4)=app.Figsize;
        %     end
        %     fnm=pfx+"ClSize";
        %     saveas(fC((j-1)*npl+1),fullfile(path,fnm+png));
        %     %saveas(fC((j-1)*npl+1),path+fnm+fig);
        %
        %     %plot2 SIZE Clusters=f(Ez) no Clusters 4
        %     ylabel(axC((j-1)*npl+2),'Number of Clusters')
        %     xlabel(axC((j-1)*npl+2),xlab)
        %     if size(pD,2)>1
        %         axC((j-1)*npl+2).YLim(2)=max(axA(7).YLim(2),axA(8).YLim(2));
        %         yticks(axC((j-1)*npl+2),'auto');
        %         fC((j-1)*npl+2).Position(3:4)=app.Figsize;
        %     end
        %     fnm=pfx+"ClSizeN4";
        %     saveas(fC((j-1)*npl+2),fullfile(path,fnm+png));
        %     %saveas(fC((j-1)*npl+2),path+fnm+fig);
        %
        %     %plot3 SIZE Max Clustersize
        %     ylabel(axC((j-1)*npl+3),'Max Size')
        %     xlabel(axC((j-1)*npl+3),xlab)
        %     if size(pD,2)>1
        %         axC((j-1)*npl+3).YLim(2)=axA(10).YLim(2);
        %         yticks(axC((j-1)*npl+3),'auto');
        %         fC((j-1)*npl+3).Position(3:4)=app.Figsize;
        %     end
        %     fnm=pfx+"ClSizeMax";
        %     saveas(fC((j-1)*npl+3),fullfile(path,fnm+png));
        %     %saveas(fC((j-1)*npl+3),path+fnm+fig);
        %
        %     %plot4 SIZE percentile
        %     ylabel(axC((j-1)*npl+4),'Size')
        %     xlabel(axC((j-1)*npl+4),xlab)
        %     if size(pD,2)>1
        %         axC((j-1)*npl+4).YLim(2)=axA(12).YLim(2);
        %         yticks(axC((j-1)*npl+4),'auto');
        %         fC((j-1)*npl+4).Position(3:4)=app.Figsize;
        %     end
        %     fnm=pfx+"ClSizePctile";
        %     saveas(fC((j-1)*npl+4),fullfile(path,fnm+png));
        %     %saveas(fC((j-1)*npl+4),path+fnm+fig);

        %Ratio RATIO Nb=f(Ez)
        ylabel(axD((j-1)*nPct+1),'Cluster Ratio')
        xlabel(axD((j-1)*nPct+1),xlab)
        if size(pD,2)>1
            axB((j-1)*nOrd+1).YLim(2)=axA(10).YLim(2);
            yticks(axB((j-1)*nOrd+1),'auto');
            fD((j-1)*nPct+1).Position(3:4)=app.Figsize;
        end
        if strcmpi(x,'p');axD((j-1)*nPct+1).XLim(1)=0;end
        fnm=pfx+"Cl_Order_Pct";
        saveas(fD((j-1)*nPct+1),fullfile(path,fnm+png));

        %Ratio RATIO Nb=f(Ez) no Clusters 4
        ylabel(axD((j-1)*nPct+2),'Clusters Ratio')
        xlabel(axD((j-1)*nPct+2),xlab)
        if size(pD,2)>1
            axB((j-1)*nOrd+2).YLim(2)=max(axA(11).YLim(2),axA(12).YLim(2));
            yticks(axB((j-1)*nOrd+2),'auto');
            fD((j-1)*nPct+2).Position(3:4)=app.Figsize;
        end
        if strcmpi(x,'p');axD((j-1)*nPct+2).XLim(1)=0;end
        fnm=pfx+"Cl_Order_Pct_N4";
        saveas(fD((j-1)*nPct+2),fullfile(path,fnm+png));

        %Ratio RATIO CELLS=f(Ez)
        ylabel(axD((j-1)*nPct+3),'Cells Ratio')
        xlabel(axD((j-1)*nPct+3),xlab)
        if size(pD,2)>1
            fD((j-1)*nPct+3).Position(3:4)=app.Figsize;
        end
        if strcmpi(x,'p');axD((j-1)*nPct+3).XLim(1)=0;end
        fnm=pfx+"Cl_Cells_Pct";
        saveas(fD((j-1)*nPct+3),fullfile(path,fnm+png));

        %Ratio RATIO CELLS=f(Ez) no Clusters 4
        ylabel(axD((j-1)*nPct+4),'Cells Ratio')
        xlabel(axD((j-1)*nPct+4),xlab)
        if size(pD,2)>1
            fD((j-1)*nPct+4).Position(3:4)=app.Figsize;
        end
        if strcmpi(x,'p');axD((j-1)*nPct+4).XLim(1)=0;end
        fnm=pfx+"Cl_Cells_Pct_N4";
        saveas(fD((j-1)*nPct+4),fullfile(path,fnm+png));

        if doVol
            %Ratio VOLUME =f(Ez)
            ylabel(axD((j-1)*nPct+5),'Volume Ratio')
            xlabel(axD((j-1)*nPct+5),xlab)
            if size(pD,2)>1
                fD((j-1)*nPct+5).Position(3:4)=app.Figsize;
            end
            if strcmpi(x,'p');axD((j-1)*nPct+5).XLim(1)=0;end
            fnm=pfx+"Cl_Volume_Pct";
            saveas(fD((j-1)*nPct+5),fullfile(path,fnm+png));

            %Ratio VOLUME =f(Ez) no Clusters 4
            ylabel(axD((j-1)*nPct+6),'Volume Ratio')
            xlabel(axD((j-1)*nPct+6),xlab)
            if size(pD,2)>1
                fD((j-1)*nPct+6).Position(3:4)=app.Figsize;
            end
            if strcmpi(x,'p');axD((j-1)*nPct+6).XLim(1)=0;end
            fnm=pfx+"Cl_Volume_Pct_N4";
            saveas(fD((j-1)*nPct+6),fullfile(path,fnm+png));
        end

        %SIZE/ORDER
        ylabel(axE((j-1)*nCmp+1),'Size-Order')
        xlabel(axE((j-1)*nCmp+1),xlab)
        if strcmpi(x,'p');axE((j-1)*nCmp+1).XLim(1)=0;end
        fnm=pfx+"Cl_Size_Ord";
        fE((j-1)*nCmp+1).Position(3:4)=app.Figsize;
        saveas(fE((j-1)*nCmp+1),fullfile(path,fnm+png));
        %saveas(fE((j-1)*nCmp+1),path+fnm+fig);

        %Total Nb of Cells
        ylabel(axE((j-1)*nCmp+2),'Nb of Cells')
        xlabel(axE((j-1)*nCmp+2),xlab)
        if strcmpi(x,'p');axE((j-1)*nCmp+2).XLim(1)=0;end
        fnm=pfx+"Cl_Nb_of_Cells";
        fE((j-1)*nCmp+2).Position(3:4)=app.Figsize;
        saveas(fE((j-1)*nCmp+2),fullfile(path,fnm+png));
        %saveas(fE((j-1)*nCmp+2),path+fnm+fig);

    end

end

%If legend is turned off, create a legend file that may be used as an
%outside legend
if ~leg
    %Multi legend
    if numel(pD)>1
        ii=copyobj(fA(1),0);
        l=legend(ii.CurrentAxes,pD.FileName,'NumColumns',6,...
            'Orientation','horizontal');
        l.EdgeColor='none';
        set(ii.CurrentAxes,'Visible','Off')
        % Set the figure Position using the normalized legend Position vector
        % as a multiplier to the figure's current position in pixels This sets
        % the figure to have the same size as the legend
        set(ii,'Position',(get(l,'Position').*[0, 0, 1, 1].*get(ii,'Position')));
        % The legend is still offset so set its normalized position vector to
        % fill the figure
        set(l,'Position',[0,0,1,1]);
        % Put the figure back in the middle screen area
        set(ii, 'Position', get(ii,'Position') + [500, 400, 0, 0]);
        saveas(ii,pathA+"Legend_Multi"+png);
        delete(ii)
    end
    if numel(pD)==1 || pltIndv
        %Order 4+ Legend
        p=copyobj(fB(1),0);
        l=legend(p.CurrentAxes,"Small","Submedium","Medium",...
            "Large",'Orientation','horizontal');
        l.EdgeColor='none';
        set(p.CurrentAxes,'Visible','Off')
        set(p,'Position',(get(l,'Position').*[0, 0, 1, 1].*get(p,'Position')));
        set(l,'Position',[0,0,1,1]);
        set(p, 'Position', get(p,'Position') + [500, 400, 0, 0]);
        saveas(p,pathA+"Legend_Order"+png);

        %order 6+Legend
        q=copyobj(fB(2),0);
        l=legend(q.CurrentAxes,"Submedium","Medium",...
            "Large",'Orientation','horizontal');
        l.EdgeColor='none';
        set(q.CurrentAxes,'Visible','Off')
        set(q,'Position',(get(l,'Position').*[0, 0, 1, 1].*get(q,'Position')));
        set(l,'Position',[0,0,1,1]);
        set(q, 'Position', get(q,'Position') + [500, 400, 0, 0]);
        saveas(q,pathA+"Legend_Order_No4"+png);

        %Delete extra figures
        delete([p,q])
    end
end

if isequal(app.LIGGGHTSAnalysisButtonGroup.SelectedObject,...
        app.ExeAllButton) || size(pD,2)>1 
    try delete(fB);catch;end
    try delete(fD);catch;end
    try delete(fE);catch;end
end
end
function clAniPlotter(pD,app)
%CLANIPLOTTER plot the evolution of anisotropy of each loop
% Prepare the values calculated by the main function to plot into different
% graphs. 
% Elevation, analysis is made from 0 to 90 degrees for the Azimuth from and
% -90 to 90.
% Surface and Gravity are ploted side by side, while fabric is seperatadly.
% They represent two different properties of the clusters.


xlab='Axial strain';
x='e';
xVals=pD.Results.Strain;
if pD.SimType==3
    xlab='Mean pressure (kPA)';
    x='p';
    xVals=pD.Results.Pressure;
end
lw={'LineWidth',1.5};
set(0,'defaultAxesFontSize',app.FontSizeEF.Value)

% Initial variables
res=pD.Results.Anisotropy{1};%Nstep x {Nloops x [Order SurfVal GrvtVal] x Category}
aniType= ["Surface","Gravity","Fabric"];     %legend for anisotropy type
clCat = ["Cl6","Cl8","Cl22"];%legend for cluster categories
aD = ["Elevation","Azimuth"];       %legend for angle type
nCag=6;                     %Number of angle categories
nAt=size(res.Angles,2);     %nb of anisotropy types 
nAa=size(res.Angles,3);     %nb of anisotropy angles
dA=90/nCag;                 %nb of angle divisions
nclCat=size(clCat,2);       %nb of cluster categories
clrRGB = graphClrCode(nCag);%plot colorcode
clrRGBL= graphClrCode(2*nCag);%plot colorcode
lgStr=string.empty(0,1);    %plot legends
infP = pD.InfPts;
%plot legends
for aA=-nCag+1:nCag
    if aA==nCag
        lgStr(aA+nCag)=sprintf("[%05.2f , %05.2f]",(aA-1)*dA,aA*dA);
    else
        lgStr(aA+nCag)=sprintf("[%05.2f , %05.2f[",(aA-1)*dA,aA*dA);
    end
end

%Check title and legends option
if app.TitlesCB.Value;tit=1;else;tit=0;end
if app.LegendsCB.Value;leg=1;else;leg=0;end

%Separate angles by category. Using a 5D matrix :
%Step X ClCat X AnglCat X AniType X AngleDirection
angPerOrdNbr=zeros(numel(pD.Results),nclCat,nCag,nAt,nAa);
angPerOrdRat=zeros(numel(pD.Results),nclCat,nCag,nAt,nAa);
devRes=zeros(nclCat,nAt,2);
for s=1:numel(pD.Results.Anisotropy) %step
    res=pD.Results.Anisotropy{s};
    ang=res.Angles; %get only the angle values
    %divide them by 'dA' to make easy to categorise
    ang=floor(ang/dA);
    clOd=res.Order;
    for c=1:nclCat
        %for each Loop category
        if c<2
            Ci=find(clOd==6); %6
        elseif c==2
            Ci=find(clOd>=10 & clOd<=20);
        else
            Ci=find(clOd>20);
            %if there are too few values do not plot
            if numel(Ci)<20
                angPerOrdNbr(s,c,:,:,:)=NaN;
                angPerOrdRat(s,c,:,:,:)=NaN;
                %devRes{cl,2}=[devRes{cl,2}; NaN(1,nCag)];
                continue
            end
        end
        angS=ang(Ci,:,:);
        for aA=1:2*nCag
            %for each angle category calculate the nb and ratio of clusters
            %on the angle direction
            if aA==2*nCag
                chk=(angS==aA-7 | angS==aA-6); %catch the case angle=90 too
            else
                chk=(angS==aA-7);
            end
            angPerOrdNbr(s,c,aA,:,:)=sum(chk,1);
            angPerOrdRat(s,c,aA,:,:)=sum(chk,1)/size(Ci,1);
            %for the deviatoric value we will separate it following the
            %elevation categories
%             for aT=1:nAty
%                 elev=chk(:,aT,1);%elevation only
%                 devRes{cl,1,aA,aT}=[devRes{cl,1,aA,aT};...
%                     [ones(sum(elev),1)*strain(s) dev(Ci(elev),:)]];
%                 devRes{cl,2,aA,aT}=[devRes{cl,2,aA,aT};...
%                     strain(s) mean(dev(Ci(elev),:),1)];
%             end
        end
        %prepare deviatoric values
        dev=res.Values(:,:,1);
        devRes(s,:,:)=cat(3,mean(dev,1),var(dev,0,1));
    end
end
%PLOTGROUP 1 - 
nAtn=nAt-1; %take into account the dif between the two properties
if strcmp(pD.Prefix,"Load")
    path=MakePath(app,'LOOPAL');del=0;
elseif strcmp(pD.Prefix,"total")
    path=MakePath(app,'LOOPA');del=0;
else
    path=fullfile(MakePath(app,'LOOPA'),pD.Prefix);del=1;
    if exist(path,'dir')==0;mkdir(path);end
end
png=".png";
lTyp=["-";"-|";"-x"];
    %initiate variables
f=gobjects(nclCat,8);    %figures
axf=gobjects(nAt,4);    %axes
li=min(2,ceil(nAtn/2)); %nb of lines per subplot
cl=ceil(nAtn/li);       %nb of columns per subplot
aA1=1;aA2=2;            %separate Azimuth from elevaiton
tits=[aD(aA1),aD(aA2),aD(aA1),aD(aA2)]; %spec for title/name
nm=["Ratio";"Number"];                 %spec for title/name
    %two cells one with half the legends (pos values) and one with all
lSt={1:numel(lgStr),... 
    numel(lgStr):-1:nCag+1};
for c=1:nclCat
    %for each loopsize we create a figure containing a tiledlayout
    f(c,1)=figure;f(c,2)=figure;f(c,3)=figure;f(c,4)=figure;
    %f(i).Position([3,4])=f(i).Position([3,4])*2;
    tf(4)=tiledlayout(f(c,4),li,cl,'TileSpacing','tight','Padding','Compact');
    tf(1)=tiledlayout(f(c,1),li,cl,'TileSpacing','tight','Padding','Compact');
    tf(2)=tiledlayout(f(c,2),li,cl,'TileSpacing','tight','Padding','Compact');
    tf(3)=tiledlayout(f(c,3),li,cl,'TileSpacing','tight','Padding','Compact');
    plRat1=permute(angPerOrdRat(:,c,:,:,aA1),[1,3,4,2,5]); %elevation
    plRat2=permute(angPerOrdRat(:,c,:,:,aA2),[1,3,4,2,5]); %azimuth
    plNbr1=permute(angPerOrdNbr(:,c,:,:,aA1),[1,3,4,2,5]); %elevation
    plNbr2=permute(angPerOrdNbr(:,c,:,:,aA2),[1,3,4,2,5]); %azimuth
    for aT=1:nAt %for each anitype
        if aT==nAt
            f(c,5)=figure;axf(aT,1)=axes(f(c,5));hold(axf(aT,1),'on');
            f(c,6)=figure;axf(aT,2)=axes(f(c,6));hold(axf(aT,2),'on');
            f(c,7)=figure;axf(aT,3)=axes(f(c,7));hold(axf(aT,3),'on');
            f(c,8)=figure;axf(aT,4)=axes(f(c,8));hold(axf(aT,4),'on');
        else
            axf(aT,1)=nexttile(tf(1));hold(axf(aT,1),'on'); 
            axf(aT,2)=nexttile(tf(2));hold(axf(aT,2),'on');
            axf(aT,3)=nexttile(tf(3));hold(axf(aT,3),'on');
            axf(aT,4)=nexttile(tf(4));hold(axf(aT,4),'on');
        end

        for p=1:2*nCag
            if p<=nCag
                %Elevation ratio 
                plot(axf(aT,1),xVals,sum(plRat1(:,[p,2*nCag+1-p],aT),2),...
                    lTyp(rem(p-1,numel(lTyp))+1),'Color',clrRGB(p,:),lw{:})
                %Elevation number 
                plot(axf(aT,3),xVals,sum(plNbr1(:,[p,2*nCag+1-p],aT),2),...
                    lTyp(rem(p-1,numel(lTyp))+1),'Color',clrRGB(p,:),lw{:})
            end
            %Azimuth ratio 
            plot(axf(aT,2),xVals,plRat2(:,p,aT),...
                lTyp(rem(p-1,numel(lTyp))+1),'Color',clrRGBL(p,:),lw{:})
            %Azimuth number 
            plot(axf(aT,4),xVals,plNbr2(:,p,aT),...
                lTyp(rem(p-1,numel(lTyp))+1),'Color',clrRGBL(p,:),lw{:})
        end
        if tit
            title(axf(aT,1),aniType(aT)+"-"+clCat(c))
            title(axf(aT,2),aniType(aT)+"-"+clCat(c))
            title(axf(aT,3),aniType(aT)+"-"+clCat(c))
            title(axf(aT,4),aniType(aT)+"-"+clCat(c))    
        end
    end
    
    %Save plots
    for i=1:4
        %put all the graphs in the sime scale and add inflection points lines
        yl=cat(1,axf(1:nAtn,i).YLim);
        yl=[min(yl(:,1)) max(yl(:,2))];
        [axf(1:nAtn,i).YLim]=deal(yl);
        xl=cat(1,axf(1:nAtn,i).XLim);
        xl=[min(xl(:,1)) max(xl(:,2))];
        [axf(1:nAtn,i).XLim]=deal(xl);
        for aT=1:nAt
            if strcmp(x,'p')
                pts=infP.p;
            else
                pts=infP.ez;
            end
            for k=1:numel(pts) %nb of pts
                if pts(k)>=max(xVals) || pts(k)<=min(xVals);continue;end
                xl=xline(axf(aT,i),pts(k),'--','Color',...
                    '#C1C1C1',lw{:});
                set(get(get(xl,'Annotation'),'LegendInformation'),...
                    'IconDisplayStyle','off');
            end
        end
        
        %SAVE SURFACE AND GRAVITY anisotropy
        ylabel(tf(i),nm(ceil(i/2)) + " of Clusters")
        xlabel(tf(i),xlab)
        if tit;title(tf(i), tits(i)+" anisotropy for the "+clCat(c));end
        f(c,i).Position(3:4) = f(c,i).Position(3:4).*[min(2,cl),min(2,li)] ;
        if leg
            lgd=legend(axf(1,i),lgStr(lSt{rem(i,2)+1}),'NumColumns',1);
            lgd.Layout.Tile = 'east';
        end
        fnm=tits(i)+"_"+clCat(c)+"_"+nm(ceil(i/2)); 
        saveas(f(c,i),fullfile(path,fnm+png));

        %SAVE FABRIC anisotropy
        ylabel(axf(end,i),'Ratio of Clusters')
        xlabel(axf(end,i),xlab)
        if tit
            title(axf(end,i), aniType(end)+" "+tits(i)+...
                " anisotropy for the "+clCat(c))
        end
        if leg
            legend(axf(end,i),lgStr(lSt{rem(i,2)+1}),'NumColumns',1,...
                'Location','eastoutside');
        end
        f(c,i+4).Position(3) = f(c,i+4).Position(3)+100;
        saveas(f(c,i+4),path+fnm+"_"+aniType(end)+png);
    end
    
end

%PLOTGROUP 2 - Deviatoric
%prepare figures
% h(nclCat)=figure;
% th(nclCat)=tiledlayout(h(nclCat),1,nAty,'TileSpacing',...
%     'tight','Padding','Compact');
% axh=gobjects(1,nAty);
% for c=1:nclCat %each cluster size category
%     if c~=nclCat
%         h(c)=figure;
%         th(c)=tiledlayout(h(c),1,nAty,'TileSpacing',...
%            'tight','Padding','Compact');
%     end
%     for aT=1:nAty
%         axh(aT)=nexttile(th(c));hold(axh(aT),'on');
%         %plot
%         plot(axh(aT),strain,devRes(:,aT,1))
%         ylabel(th(c),'Mean value')
%         yyaxis(axh(aT),'right') 
%         plot(axh(aT),strain,devRes(:,aT,2))
%         ylabel(th(c),'Variance')
%         title(aniType(aT))
%     end
%     %Fix axis
%     yl=cat(1,axh(:).YLim);
%     yl=[min(yl(:,1)) max(yl(:,2))];
%     for aT=1:nAty
%         axh(aT).YLim=yl;
%     end
%     title(th(c), "Deviatoric anisotropy for the "+clCat(c))
%     xlabel(th(c),xlab)
%     h(c).Position(3) = 2*h(c).Position(3) + 250;
%     fnm="deviatoric of "+clCat(c);
%     saveas(h(c),fullfile(path,fnm+png));
% end

if ~leg
    %Legend of Elevation
    p=copyobj(f(1,1),0);
    l=legend(p.CurrentAxes,lgStr(lSt{2}),'NumColumns',1);
    l.EdgeColor='none';
    set(p.CurrentAxes,'Visible','Off')
    set(p,'Position',(get(l,'Position').*[0, 0, 1, 1].*get(p,'Position')));
    set(l,'Position',[0,0,1,1]);
    set(p, 'Position', get(p,'Position') + [500, 400, 0, 0]);
    saveas(p,path+"Legend_Elev"+png);
    delete(p);
    %Legend of Azimuth
    p=copyobj(f(1,2),0);
    l=legend(p.CurrentAxes,lgStr(lSt{1}),'NumColumns',1);
    l.EdgeColor='none';
    set(p.CurrentAxes,'Visible','Off')
    set(p,'Position',(get(l,'Position').*[0, 0, 1, 1].*get(p,'Position')));
    set(l,'Position',[0,0,1,1]);
    set(p, 'Position', get(p,'Position') + [500, 400, 0, 0]);
    saveas(p,path+"Legend_Azi"+png);
    delete(p);
end
if isequal(app.LIGGGHTSAnalysisButtonGroup.SelectedObject,app.ExeAllButton)...
        || del==1
    try delete([f,h]);
    catch ME
        if strcmp(ME.identifier,'MATLAB:UndefinedFunction')
            delete(f);
        else
            rethrow(ME)
        end
    end
else
    %go for Rosete plotter
    ClusterAniPlotter(pD,aniType)
end
end
function clDefPlotter(pD,app)
%CLDEFPLOTTER plot Loop variation, Size/Order and Deformability
%    Prepare the values calculated by the main function plot the evolution
%    of the derivative of loops, deformability and Size vs Order values

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

try lw=app.PlotWidthEF.Value;
catch
    lw=1.5;
end
set(0,'defaultAxesFontSize',app.FontSizeEF.Value)
boolNbEdges=1;

%Prepare variables
nF=numel(pD);%nb of files loaded
nb=5;   %nb of variables to plot
%prepare files
if nF>1 
    %if more than one file was loaded
    Res=cell(nF,1);
    defRes=double.empty(0,3);           %def results for O>4
    defRes4=[(0:4)'/6 , zeros(5,1)];    %def results for O==4
    resNbE=double.empty(0,3);           %account nb of edges if this is calculated
    for j=1:nF
        r=cat(1,pD(j).Results.SOD{:,1});
        [avrgS,axO]= groupsummary(r(:,2), r(:,1), @mean);
        [avrgSoO,~]= groupsummary(r(:,2)./r(:,1), r(:,1), @mean);
        [avrgDfO,~]= groupsummary(r(:,3), r(:,1), @mean);
        Res{j}=[axO,avrgS,avrgSoO,avrgDfO]; %Order - AveSize - AveRat - AveDef
        %defRes=cat(1,defRes,unique(r,'rows'));
        if size(pD(j).Results.SOD,2)==1
            %save results for O==4
            cl4=r(:,1)==4;
            [nbD4,d4]=groupcounts(r(cl4,3));
            [chk,ia]=ismember(defRes4(:,1),d4);
            defRes4(chk,2)=defRes4(chk,2)+nbD4(ia(ia>0));
            %save other results
            defRes=cat(1,defRes,r(~cl4,1:3));
            %if one of the files do not have edges info 
            boolNbEdges=0;
        else
            %save other results
            defRes=cat(1,defRes,r(:,1:3));
            %save results 4
            r4=cat(1,pD(j).Results.SOD{:,2});
            [nbD4,d4]= groupsummary(r4(:,2), r4(:,1), @sum);
            [chk,ia]=ismember(defRes4(:,1),d4);
            defRes4(chk,2)=defRes4(chk,2)+nbD4(ia(ia>0));
            %nb of edges
            resNbE=cat(1,resNbE,r(:,[1,4,5]));
        end
    end
    nb=nb+3;   %plotting only distribution for the total file
else
    r=cat(1,pD.Results.SOD{:,1});
    if size(pD.Results.SOD,2)==1
        %save results for O==4
        cl4=r(:,1)==4;
        [nbD4,d4]=groupcounts(r(cl4,3));
        defRes4=[d4,nbD4];
        %save other results
        defRes=r(~cl4,1:3);
        boolNbEdges=0;
    else
        %save other results
        defRes=r(:,1:3);
        %save results 4
        r4=cat(1,pD.Results.SOD{:,2});
        [nbD4,d4]= groupsummary(r4(:,2), r4(:,1), @sum);
        defRes4=[d4,nbD4];
        %nb of edges
        resNbE=r(:,[1,4,5]);
    end
end
nb=nb+2*boolNbEdges;

%Start figures
f(nb)=figure;ax(nb)=axes(f(nb));hold(ax(nb),'on');
for i=1:(nb-1)
    f(i)=figure;ax(i)=axes(f(i));hold(ax(i),'on');
end

%Start Plots
path=MakePath(app,'LOOPDL');
png=".png";

%%%%%%%%%%%%%%%%%%%%%%%%%% SCATTER SIZE=F(ORDER) %%%%%%%%%%%%%%%%%%%%%%%%%%
SnO=[4,4;defRes(defRes(:,1)>0,1:2)];
i=1;
[avrgS,axS]= groupsummary(SnO(:,2), SnO(:,1), @mean);
%if max(axS)>300;f(i).Position(3)=2*f(i).Position(3);end
sct=unique(SnO,'rows');
scatter(ax(i),sct(:,1),sct(:,2),'x','MarkerEdgeColor',[0.6 0.6 0.6])
plot(ax(i),axS,avrgS,'k','LineWidth',lw)
ly=ax(i).YLim(2)*.5;
w=1:(-1/(numel(axS)-1)):0;
ft=fit(axS,avrgS,'poly1','Weights',w);
axes(ax(i))
p=plot(ft,'--k',axS,avrgS);
p(2).LineWidth=lw;
text(ax(i),ly,ly*ft.p1+ft.p2,...
    {'',[' \leftarrow S=' num2str(ft.p1,'%.3f') '*O+' num2str(ft.p2,'%.3f')]},...
    'HorizontalAlignment','left','Color','black')
delete(ax(i).Children(3))
hLeg=findobj(f(i),'type','legend');
set(hLeg,'visible','off')

%%%%%%%%%%%%%%%%%%%%%% SCATTER (SIZE/ORDER)=F(ORDER) %%%%%%%%%%%%%%%%%%%%%%
so=SnO(:,2)./SnO(:,1);
i=i+1;
[avrgS,axO]= groupsummary(so, SnO(:,1), @mean);
%if max(axSoO)>300;f(i).Position(3)=2*f(i).Position(3);end
sct=unique([SnO(:,1),so],'rows');
scatter(ax(i),sct(:,1),sct(:,2),'x','MarkerEdgeColor',[0.6 0.6 0.6])
plot(ax(i),axO,avrgS,'k','LineWidth',lw)
estO=4:2:max(axO);
plot(ax(i),estO,(2+estO/2)./estO,'r','LineWidth',lw)

%%%%%%%%%%%%%%%%%%%%%%%%%%% SCATTER DEF=F(ORDER) %%%%%%%%%%%%%%%%%%%%%%%%%%%
i=i+1;
[avrgDfO,axDfO]= groupsummary(defRes(:,3), defRes(:,1), @mean);
%add average 4 def
avrgDfO=[dot(defRes4(:,1),defRes4(:,2))/sum(defRes4(:,2));avrgDfO];
axDfO=[4;axDfO];
    %get unique values and plot
udef=[ones(5,1)*4,(0:4)'/6; unique([defRes(:,1),defRes(:,3)],'rows')];
scatter(ax(i),udef(:,1),udef(:,2),'x','MarkerEdgeColor',[0.6 0.6 0.6])
    %get min and max values for each order
maxD= groupsummary(defRes(:,3), defRes(:,1), @max); maxD=[4/6;maxD];
minD= groupsummary(defRes(:,3), defRes(:,1), @min); minD=[0;minD];
    %check for places where there is no dif between vals and average
chk=find((avrgDfO-minD)./avrgDfO<=0.01 | (maxD-avrgDfO)./avrgDfO<=0.01);
m=minD;mx=maxD;
for l=1:numel(chk)
    if chk(l)<=3
        m(1)=mean(mink(minD(1:7),3));
        mx(1)=mean(maxk(maxD(1:5),3));
    elseif chk(l)>=(numel(minD)-2)
        m(end)=mean(mink(minD(end-6:end),3));
        mx(end)=mean(maxk(maxD(end-6:end),3));
    else
        m(chk(l))=mean(mink(minD(chk(l)-3:chk(l)+3),3));
        mx(chk(l))=mean(maxk(maxD(chk(l)-3:chk(l)+3),3));
    end
end
minD=m;maxD=mx;
w=ones(numel(axDfO),1);
w(axDfO>220)=500;
    %modify them, fit values
ftmax=fit(axDfO,maxD,fittype({'x','1'}),'Weight',w);
ftmin=fit(axDfO,minD,fittype({'x','log(x)','1'}));
axes(ax(i))
p=plot(ftmax,'-b',axDfO,maxD);
p(2).LineWidth=lw;
p=plot(ftmin,'-b',axDfO,minD);
p(2).LineWidth=lw;
    %remove the autogenerated legends
hLeg=findobj(f(i),'type','legend');
set(hLeg,'visible','off')
plot(ax(i),axDfO,avrgDfO,'k','LineWidth',lw)
    %remove the generated pts
delete(ax(i).Children([3,5]))

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Count Def values %%%%%%%%%%%%%%%%%%%%%%%%%%%%
i=i+1;
%count the nb of unique deformability values per order
cntDef=unique(defRes(:,[1,3]),'rows');
nbUDefO = groupcounts(cntDef(:,1));
%count the number of def values per order
[nbDefO,axO] = groupcounts(defRes(:,1));
plot(ax(i),axO,nbUDefO,'Color','#000000','LineWidth',lw)
%plot the ratio between nb of unique values over total nb o values. 0 or
%close indicates many repetitions; 1 or close means no repetition.
yyaxis(ax(i),'right')
plot(ax(i),axO,nbUDefO./nbDefO,'Color','#a1a1a1','LineWidth',lw)
ax(i).YAxis(1).Color='k';
ax(i).YAxis(2).Color='#a1a1a1';

%%%%%%%%%%%%%%%%%%%%%%%%% Granulometry like curve %%%%%%%%%%%%%%%%%%%%%%%%%
i=i+1;
%nb values in each category has already been calculated as nbDefO in 
%previous section, so only the curve need to be created.
pct=0.99;
cs6=cumsum(nbDefO);cs6=cs6/cs6(end);
[~,mn]=min(abs(cs6-pct));O6=axO(mn);
p=plot(ax(i),axO,cs6,'LineWidth',lw);
xl=xline(ax(i),O6,"--",'Color',p.Color);
set(get(get(xl,'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','off');

cs8=cumsum(nbDefO(axO>6));cs8=cs8/cs8(end);
[~,mn]=min(abs(cs8-pct));O8=axO(mn+sum(~(axO>6)));
p=plot(ax(i),axO(axO>6),cs8,'LineWidth',lw);
xl=xline(ax(i),O8,"--",'Color',p.Color);
set(get(get(xl,'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','off')

cs20=cumsum(nbDefO(axO>20));cs20=cs20/cs20(end);
[~,mn]=min(abs(cs20-pct)); O20=axO(mn+sum(~(axO>20)));
p=plot(ax(i),axO(axO>20),cs20,'LineWidth',lw);
xl=xline(ax(i),O20,"--",'Color',p.Color);
set(get(get(xl,'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','off')

legend(ax(i),"O>4","O>6","O>20","location","southeast")
ax(i).XLim(1)=6;
yl=yline(ax(i),pct,"--","Color","#A1A1A1");
set(get(get(yl,'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','off');

text(ax(i),ax(i).XLim(2)*.3,ax(i).YLim(2)*.5,...
        {'',"99% of clusters of O>4  is Order "+num2str(O6),...
        "99% of clusters of O>6  is Order "+num2str(O8),...
        "99% of clusters of O>20 is Order "+num2str(O20)},...
        'HorizontalAlignment','left','Color','black')
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Number of Edges %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if boolNbEdges
    %%%%%%%%%%%%%%%%% Total Edges %%%%%%%%%%%%%%%%%
    i=i+1;
    [avrgN,axNfO]= groupsummary(resNbE(:,2), resNbE(:,1), @mean);
    %values
    scatter(ax(i),resNbE(:,1),resNbE(:,2),'x','MarkerEdgeColor',[0.6 0.6 0.6])
    %get min and max values for each order
    maxD= groupsummary(resNbE(:,2), resNbE(:,1), @max);
    minD= groupsummary(resNbE(:,2), resNbE(:,1), @min);
    %check for places where there is no dif between vals and average
    chk=find((avrgN-minD)./avrgN<=0.01 | (maxD-avrgN)./avrgN<=0.01);
    m=minD;mx=maxD;
    for l=1:numel(chk)
        if chk(l)<=3
            m(1)=mean(mink(minD(1:7),3));
            mx(1)=mean(maxk(maxD(1:5),3));
        elseif chk(l)>=(numel(minD)-2)
            m(end)=mean(mink(minD(end-6:end),3));
            mx(end)=mean(maxk(maxD(end-6:end),3));
        else
            m(chk(l))=mean(mink(minD(chk(l)-3:chk(l)+3),3));
            mx(chk(l))=mean(maxk(maxD(chk(l)-3:chk(l)+3),3));
        end
    end
    minD=m;maxD=mx;
    axDfO=axDfO(2:end);
    plot(ax(i),axNfO,avrgN,'k','LineWidth',lw) %mean line
%     ftMax=fit(axDfO,maxD,fittype({'x','log(x)','1'}));
%     ftMin=fit(axDfO,minD,fittype({'x','log(x)','1'}));
    ftMax=fit(axNfO,maxD,'poly1');
    ftMin=fit(axNfO,minD,'poly1');
    axes(ax(i))
    p1=plot(ftMax,'b',axNfO,maxD);
    p1(2).LineWidth=lw;
    p2=plot(ftMin,'b',axNfO,minD);
    p2(2).LineWidth=lw;
    text(ax(i),ax(i).XLim(2)*.05,ax(i).YLim(2)*.95,...
        {'',"N(max)="+num2str(ftMax.p1,'%.3f')+"*O+"+num2str(ftMax.p2,'%.3f'),...
        "N(min)="+num2str(ftMin.p1,'%.3f')+"*O+"+num2str(ftMin.p2,'%.3f')},...
        'HorizontalAlignment','left','Color','blue')
    hLeg=findobj(f(i),'type','legend');
    set(hLeg,'visible','off')
    delete(ax(i).Children([3,5]))

    %%%%%%%%%%%%%%%%% Open Edges %%%%%%%%%%%%%%%%%
    i=i+1;
    [avrgN,axNfO]= groupsummary(resNbE(:,3), resNbE(:,1), @mean);
    %values
    scatter(ax(i),resNbE(:,1),resNbE(:,3),'x','MarkerEdgeColor',[0.6 0.6 0.6])
    %get min and max values for each order
    maxD= groupsummary(resNbE(:,3), resNbE(:,1), @max);
    minD= groupsummary(resNbE(:,3), resNbE(:,1), @min);
    %check for places where there is no dif between vals and average
    chk=find((avrgN-minD)./avrgN<=0.01 | (maxD-avrgN)./avrgN<=0.01);
    m=minD;mx=maxD;
    for l=1:numel(chk)
        if chk(l)<=3
            m(1)=mean(mink(minD(1:7),3));
            mx(1)=mean(maxk(maxD(1:5),3));
        elseif chk(l)>=(numel(minD)-2)
            m(end)=mean(mink(minD(end-6:end),3));
            mx(end)=mean(maxk(maxD(end-6:end),3));
        else
            m(chk(l))=mean(mink(minD(chk(l)-3:chk(l)+3),3));
            mx(chk(l))=mean(maxk(maxD(chk(l)-3:chk(l)+3),3));
        end
    end
    minD=m;maxD=mx;
    plot(ax(i),axNfO,avrgN,'k','LineWidth',lw) %mean line
    ftMax2=fit(axNfO,maxD,'poly1');
    ftMin2=fit(axNfO,minD,'poly1');
%     ftMax2=fit(axDfO,maxD,fittype({'x','log(x)','1'}));
%     ftMin2=fit(axDfO,minD,fittype({'x','log(x)','1'}));
    axes(ax(i))
    p1=plot(ftMax2,'b',axNfO,maxD);
    p1(2).LineWidth=lw;
    p2=plot(ftMin2,'b',axNfO,minD);
    p2(2).LineWidth=lw;
    ax(i).YLim=ax(i-1).YLim;
    text(ax(i),ax(i).XLim(2)*.05,ax(i).YLim(2)*.95,...
        {'',"N(max)="+num2str(ftMax2.p1,'%.3f')+"*O+"+num2str(ftMax2.p2,'%.3f'),...
        "N(min)="+num2str(ftMin2.p1,'%.3f')+"*O+"+num2str(ftMin2.p2,'%.3f')},...
        'HorizontalAlignment','left','Color','blue')
    hLeg=findobj(f(i),'type','legend');
    set(hLeg,'visible','off')
    delete(ax(i).Children([3,5]))
    
    o=copyobj(f(3),0);axo=o.CurrentAxes;
    vs=4:2:max(axNfO);
    cb1=(ftMax2.p1*vs+ftMax2.p2)./(ftMin.p1*vs+ftMin.p2);
    cb2=(ftMax2.p1*vs+ftMax2.p2)./(ftMax.p1*vs+ftMax.p2);
    cb3=(ftMin2.p1*vs+ftMin2.p2)./(ftMax.p1*vs+ftMax.p2);
    cb4=(ftMin2.p1*vs+ftMin2.p2)./(ftMin.p1*vs+ftMin.p2);
    
    plot(axo,vs(cb1<1 & cb1>0.3),cb1(cb1<1 & cb1>0.3),...
        vs(cb2<1 & cb2>0.3),cb2(cb2<1 & cb2>0.3),...
        vs(cb3<1 & cb3>0.3),cb3(cb3<1 & cb3>0.3),...
        vs(cb4<1 & cb4>0.3),cb4(cb4<1 & cb4>0.3),'LineWidth',lw);    
    xlabel(axo,'Order')
    ylabel(axo,'Deformability')
    saveas(o,fullfile(path,'Def_Edges_Law.fig'))
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Save files %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Title
if tit
    i=1;title(ax(i),'Evolution of Size and Order')
    i=i+1;title(ax(i),'Evolution of Size over Order')
    i=i+1;title(ax(i),'Evolution of Deformability')
    i=i+1;title(ax(i),'Evolution of Nb of Unique Def values')
    i=i+1;title(ax(i),'Cumulative ratio of orders')
    if boolNbEdges
        i=i+1;title(ax(i),'Evolution of total number of Edges')
        i=i+1;title(ax(i),'Evolution of number of open edges')        
    end
end
%Scatter size v order
i=1;
xlabel(ax(i),'Order')
ylabel(ax(i),'Size')
fnm="Size_f_Order_Scat";i=1;
saveas(f(i),fullfile(path,fnm+png));

%Scatter size/order v order
i=i+1;
xlabel(ax(i),'Order')
ylabel(ax(i),'Size/Order')
if max(axDfO)>400
    L=ax(i).XLim;
    ax(i).XLim=[0,200];
    fnm="Size_o_Order_Scat_LowO";
    saveas(f(i),fullfile(path,fnm+png));
    ax(i).XLim=L;
    %f(i).Position(3)=2*f(2).Position(3);
end
fnm="Size_o_Order_Scat";
saveas(f(i),fullfile(path,fnm+png));

%Scatter deformability
i=i+1;
xlabel(ax(i),'Order')
ylabel(ax(i),'Deformability')
if max(axDfO)>400
    L=ax(i).XLim;
    ax(i).XLim=[0,200];
    fnm="Def_Order_Scat_LowO";
    saveas(f(i),fullfile(path,fnm+png));
    ax(i).XLim=L;
    %f(i).Position(3)=2*f(2).Position(3);
end
fnm="Def_Order_Scat";
saveas(f(i),fullfile(path,fnm+png));

%Unique def values
i=i+1;
xlabel(ax(i),'Order')
yyaxis(ax(i), 'left')
ylabel(ax(i),'Nb of Unique Deformabilities')
yyaxis(ax(i), 'right')
ylabel(ax(i),'Ratio of Unique Deformability')
fnm="Def_Order_Unique";
saveas(f(i),fullfile(path,fnm+png));

%Cumulative ratio
i=i+1;
xlabel(ax(i),'Order')
ylabel(ax(i),'Cumulative Ratio')
fnm="Cum_Rat_Order";
saveas(f(i),fullfile(path,fnm+png));

if boolNbEdges
    %Total edges
    i=i+1;
    xlabel(ax(i),'Order')
    ylabel(ax(i),'Number of Edges')
    fnm="Edges_Total";
    saveas(f(i),fullfile(path,fnm+png));
    
    %Open Edges
    i=i+1;
    xlabel(ax(i),'Order')
    ylabel(ax(i),'Number of Open Edges')
    fnm="Edges_Open";
    saveas(f(i),fullfile(path,fnm+png));
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% If multifiles %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Save mean values of each category
if nF>1
    % plot mean in the 'mean comparaison graphs'
    lt=["-",":","+"];
    for j=1:nF
        r=Res{j};
        plot(ax(i+1),r(:,1),r(:,2),lt(ceil(j/7)),'LineWidth',lw)
        plot(ax(i+2),r(:,1),r(:,3),lt(ceil(j/7)),'LineWidth',lw)
        plot(ax(i+3),r(:,1),r(:,4),lt(ceil(j/7)),'LineWidth',lw)
    end
    %Title
    if tit
        title(ax(i+1),'Evolution of Size and Order')
        title(ax(i+2),'Evolution of Size over Order')
        title(ax(i+3),'Evolution of Deformability')
    end
    if leg
        legend(ax(i+1),pD.FileName,'location','best')
        legend(ax(i+2),pD.FileName,'location','best')
        legend(ax(i+3),pD.FileName,'location','best')
    end
    %Size-order
    ax(i+1).YLim=ax(1).YLim;
    ax(i+1).XLim=ax(1).XLim;
    xlabel(ax(i+1),'Order')
    ylabel(ax(i+1),'Size')
    %Size-order
    ax(i+2).YLim=ax(2).YLim;
    ax(i+2).XLim=ax(2).XLim;
    xlabel(ax(i+2),'Order')
    ylabel(ax(i+2),'Size/Order')
    %Def-order
    ax(i+3).YLim=ax(3).YLim;
    ax(i+3).XLim=ax(3).XLim;
    xlabel(ax(i+3),'Order')
    ylabel(ax(i+3),'Deformability')

    %Plot Saves
    fnm="Size_f_Order_Multi";
    saveas(f(i+1),fullfile(path,fnm+png));
    fnm="Size_o_Order_Multi";
    saveas(f(i+2),fullfile(path,fnm+png));
    fnm="Def_Order_Multi";
    saveas(f(i+3),fullfile(path,fnm+png));
end

%%%%%%%%%%%%%%%%%%%% Def evolution per order over time %%%%%%%%%%%%%%%%%%%%
lp='';%6:2:20;
if ~isempty(lp)
    %calculate the ave deformability for each step of each file
    aveDef=cell(numel(pD),1);
    for i=1:numel(pD)
        ave=zeros(numel(lp),numel(pD(i).Results.SOD));
        for j=1:numel(pD(i).Results.SOD)
            r=pD(i).Results.SOD{j};
            r=r(ismember(r(:,1),lp),:);
            [avrgDfO,avO]= groupsummary(r(:,3), r(:,1), @mean);
            ave(ismember(avO,lp),j)=avrgDfO;
        end
        aveDef{i}=ave;
    end
    if numel(pD)>6
        C=app.PlotColors;
    else
        C = graphClrCode(numel(pD));%plot colorcode
    end
    nb=numel(lp);
    path=MakePath(app,'LOOPDCL');
    g(nb)=figure;axG(nb)=axes(g(nb));hold(axG(nb),'on');
    for i=1:(nb-1)
        g(i)=figure;axG(i)=axes(g(i));hold(axG(i),'on');
    end
    
    for i=1:size(pD,2)
        r=aveDef{i};
        opts={'Color',C(i,:)};
        for k=1:numel(pD(i).InfPts.q)
            if strcmpi(x,'p')
                opts=[opts,{'Pointx',pD(i).InfPts.p(k)}]; %#ok<AGROW>
                xAx=pD(i).Results.Pressure;
            else
                opts=[opts,{'Pointx',pD(i).InfPts.ez(k)}]; %#ok<AGROW>
                xAx=pD(i).Results.Strain;
            end
        end
        for j=1:numel(lp)
            %find no zero values
            chk=r(j,:)~=0;
            if sum(chk)==0,continue;end
            %plot it
            plotMark(app,axG(j),xAx(chk),...
                r(j,chk),opts{:})
        end
    end
    
    for j=1:numel(lp)
        if size(pD,2)>1
            if leg
                legend(axG(j),pD(:).FileName,'location','best')
            end
        end
        %Loop e = f(Ez)
        if tit;title(axG(j),"Cluster "+ lp(j) +" Deformability");end
        if strcmpi(x,'p');axG(j).XLim(1)=0;end
        ylabel(axG(j),'Deformability')
        xlabel(axG(j),xlab)
        fnm="VR_Cluster"+lp(j)+".png";
        saveas(g(j),fullfile(path,fnm));
    end
end


if isequal(app.LIGGGHTSAnalysisButtonGroup.SelectedObject,app.ExeAllButton)
    delete(f);
    if ~isempty(lp);delete(g);end
    if exist(o,'var');delete(o);end
end
end
function clVRPlotter(pD,app)
%CLVRPLOTTER plot the evolution the VR of loops each timestep
% This function will take the results obtained for the Cluster void ratio
% and plot many different interesting curbs. Firstly, the pD.Results has 2
% types of values saved in different properties. 

% Data is contained in the following properties of pD
% - pD.Results.Mean - N x N x 2 matrix. Contains the mean void ratio and  
% coordination number per cluster order for each calculated step.
% - pD.Results.ClusterVRZ - cell matrix Nx2. Inside each cell there is a
% matrix containing the Order x VR x Z for each cluster identified. Each
% line of the cell matrix corresponds to a calculated step, while the first
% column contains all clusters larger 4 and the second all clusters 4.


%boolean controlling background plot of total vr
boolbg=1;

%linetype
if isequal(app.CourbePointsSwitch.Value,'On')
    lType="-+";
else
    lType="-";
end
try lw=app.PlotWidthEF.Value;
catch
    lw=1.5;
end
set(0,'defaultAxesFontSize',app.FontSizeEF.Value)

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
if pD(1).Prefix=="Load"
    path=MakePath(app,'LOOPVRL');
else
    path=MakePath(app,'LOOPVR');
end


%%%%%%%%%%%%%%%%%%%%% PLOT SCATTER (Z) - DENSITY (VR) %%%%%%%%%%%%%%%%%%%%%
nbT=2; nbPl=4;
nb=nbT*nbPl*size(pD,2);
f(nb)=figure;ax(nb)=axes(f(nb));hold(ax(nb),'on');
for i=1:(nb-1)
    f(i)=figure;ax(i)=axes(f(i));hold(ax(i),'on');
end

for i=1:size(pD,2)
    %Choose steps to plot : begining, infpts and end
    k=size(pD(i).Results.ClusterVRZ,1);
    str=[1,0,0,k];
    for j=1:numel(pD(i).InfPts.ez)
        [~,str(j+1)]=min(abs(pD(i).Results.Strain-pD(i).InfPts.ez(j)));
        if numel(pD(i).InfPts.ez)==1
            str(3)=floor((str(2)+str(4))/4)*2;
        end
    end
    for j=1:nbPl
        n=(j-1)*nbT+(i-1)*nbPl*nbT;
        vrz=[pD(i).Results.ClusterVRZ{str(j),1};...
            pD(i).Results.ClusterVRZ{str(j),2}];
        maxO=max(vrz(:,1));
        vrzMean=pD(i).Results.Mean(2:maxO/2,[1,str(j)+1],:);
        vrzMean=vrzMean(vrzMean(:,2,1)~=0,:,:);

        %plot VR density points
        [X,Y,map]=pointDensity(vrz(:,1:2),2,0.05);
        surface(ax(n+1),X,Y,map,'EdgeCOlor',[.65 .65 .65])
        plot3(ax(n+1),vrzMean(:,1,1),vrzMean(:,2,1),...
            ones(numel(vrzMean(:,2,1)),1),"k"+lType,'LineWidth',lw)   %plot mean line
        colormap(ax(n+1),flip(autumn))
        if leg;colorbar(ax(n+1));end

        %plot Z
        scatter(ax(n+2),vrz(:,1),vrz(:,3),...
            'x','MarkerEdgeColor',[0.6 0.6 0.6]) %scatter all pts
        plot(ax(n+2),vrzMean(:,1,1),vrzMean(:,2,2),...
            "k"+lType,'LineWidth',lw)   %plot mean line

        %save files
        strain=num2str(pD(i).Results.Strain(str(j)),'%.3f');
        if tit;title(ax(n+1),"Void Ratio per Order ("+strain+")");end
        ax(n+1).YLim=[0,min([4,max(cat(2,ax((i-1)*6+(1:4)).YLim))])];
        ax(n+1).XTick=4:max(8,ceil((ax(n+1).XLim(2)-4)/20)*2):2000;
        ylabel(ax(n+1),'Void Ratio')
        xlabel(ax(n+1),'Cluster Order')

        if tit;title(ax(n+2),"Coordination per Order ("+strain+")");end
        ylabel(ax(n+2),'Coordination')
        xlabel(ax(n+2),'Cluster Order')
        ax(n+2).XTick=4:max(8,ceil((ax(n+2).XLim(2)-4)/20)*2):2000;

        prefix="";
        if numel(pD)>1
            prefix=pD(i).FileName; 
        end
        fnm=prefix+"Cluster_VR_"+strain+".png";
        saveas(f(n+1),fullfile(path,fnm));
        fnm=prefix+"Cluster_Z_"+strain+".png";
        saveas(f(n+2),fullfile(path,fnm));
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOT SURFACE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if numel(pD)==1
    Z=pD(1).Results.Mean(3:end,2:end,1); %VR
    Y=pD(1).Results.Mean(3:end,1,1);     %Order
    Y=Y*ones(1,size(Z,2));
    if strcmpi(x,'p')
        X=pD(1).Results.Pressure; %Evolution axis data
    else
        X=pD(1).Results.Strain;
    end
    X=ones(size(Z,1),1)*X';
    Z(Z==0)=NaN;
    l=figure;
    surf(X,Y,Z);
    xlabel(l.CurrentAxes,xlab)
    ylabel(l.CurrentAxes,'Cluster Order')
    zlabel(l.CurrentAxes,'Void Ratio')
    fnm="Surfaceplot.png";
    l.CurrentAxes.YLim(2)=100;
    saveas(l,fullfile(path,fnm));
end


%%%%%%%%%%%%%%%%%%%%%%%%%% PLOT AVE VR EVOLUTION %%%%%%%%%%%%%%%%%%%%%%%%%%
if numel(pD)<8
    C=app.PlotColors;
else
    C = graphClrCode(size(pD,2));%plot colorcode
end
if numel(pD)>1;nb=8;else;nb=2;end

h(nb)=figure;axH(nb)=axes(h(nb));hold(axH(nb),'on');
for i=1:(nb-1)
    h(i)=figure;axH(i)=axes(h(i));hold(axH(i),'on');
end
%plot cluster VR evolution in f of Ez and p
for i=1:size(pD,2)
    if numel(pD)>1
        optsP={'Color',C(i,:)};
        optsE=optsP;
    else
        optsP={};optsE={};
    end
    for k=1:numel(pD(i).InfPts.q)
        optsP=[optsP,{'Pointx',pD(i).InfPts.p(k)}]; %#ok<AGROW>
        optsE=[optsE,{'Pointx',pD(i).InfPts.ez(k)}]; %#ok<AGROW>
    end
    xAxp=pD(i).Results.Pressure;
    xAxE=pD(i).Results.Strain;

    %Get mean values per step - For small and submedium the mean per order
    %can be taken. However, for medium and large, several orders are in
    %play so the mean must be recalculated/
    VRmn=pD(i).Results.Mean(:,:,1);
    VRZ=pD(i).Results.ClusterVRZ(:,1);
    try VRtot=pD(i).Results.TotalVR;
    catch
        boolbg=0;
    end
    %plot values removing 0s (no cluster in that pt)
        %small category
    ct=VRmn(VRmn(:,1)==4,2:end,1);
    z=~(ct==0);
    vrP=1;vrE=2;
    plotMark(app,axH(vrP),xAxp(z),ct(z),optsP{:});
    plotMark(app,axH(vrE),xAxE(z),ct(z),optsE{:});
    if boolbg && numel(pD)>1
        plotMark(app,axH(vrP),xAxp(z),VRtot(z),'Color',C(i,:),'Background','Nolegend');
        plotMark(app,axH(vrE),xAxE(z),VRtot(z),'Color',C(i,:),'Background','Nolegend');
    end
        %submedium ctegory
    ct=VRmn(VRmn(:,1)==6,2:end,1);
    z=~(ct==0);
    if numel(pD)>1;vrP=vrP+2;vrE=vrE+2;end
    plotMark(app,axH(vrP),xAxp(z),ct(z),optsP{:});
    plotMark(app,axH(vrE),xAxE(z),ct(z),optsE{:});
    if boolbg && numel(pD)>1
        plotMark(app,axH(vrP),xAxp(z),VRtot(z),'Color',C(i,:),'Background','Nolegend');
        plotMark(app,axH(vrE),xAxE(z),VRtot(z),'Color',C(i,:),'Background','Nolegend');
    end
        %medium ctegory
    VR8=cellfun(@(X) mean(X(X(:,1)>6 & X(:,1)<21,2)),VRZ,'UniformOutput',true);
    z=~(VR8==0);
    if numel(pD)>1;vrP=vrP+2;vrE=vrE+2;end
    plotMark(app,axH(vrP),xAxp(z),VR8(z),optsP{:});
    plotMark(app,axH(vrE),xAxE(z),VR8(z),optsE{:});
    if boolbg && numel(pD)>1
        plotMark(app,axH(vrP),xAxp(z),VRtot(z),'Color',C(i,:),'Background','Nolegend');
        plotMark(app,axH(vrE),xAxE(z),VRtot(z),'Color',C(i,:),'Background','Nolegend');
    end
        %large ctegory
    VR22=cellfun(@(X) mean(X(X(:,1)>20,2)),VRZ,'UniformOutput',true);
    z=~(VR22==0);
    if numel(pD)>1;vrP=vrP+2;vrE=vrE+2;end
    plotMark(app,axH(vrP),xAxp(z),VR22(z),optsP{:});
    plotMark(app,axH(vrE),xAxE(z),VR22(z),optsE{:});
    if boolbg && numel(pD)>1
        plotMark(app,axH(vrP),xAxp(z),VRtot(z),'Color',C(i,:),'Background','Nolegend');
        plotMark(app,axH(vrE),xAxE(z),VRtot(z),'Color',C(i,:),'Background','Nolegend');
    elseif boolbg
        plotMark(app,axH(1),xAxp(z),VRtot(z),'Background');
        plotMark(app,axH(2),xAxE(z),VRtot(z),'Background');
    end

end

%Add legends or create legend file
if size(pD,2)>1
    if leg
        for l=1:nb
            legend(axH(l),pD(:).FileName,'location','best')
        end
    else
        o=copyobj(h(1),0);
        l=legend(o.CurrentAxes,pD.FileName);
        l.EdgeColor='none';
        set(o.CurrentAxes,'Visible','Off')
        % Set the figure Position using the normalized legend Position vector
        % as a multiplier to the figure's current position in pixels This sets
        % the figure to have the same size as the legend
        set(o,'Position',(get(l,'Position').*[0, 0, 1, 1].*get(o,'Position')));
        % The legend is still offset so set its normalized position vector to
        % fill the figure
        set(l,'Position',[0,0,.9999,.9999]);
        % Put the figure back in the middle screen area
        set(o, 'Position', get(o,'Position') + [500, 400, 0, 0]);
        saveas(o,fullfile(path,"Legend.png"));
        
        %Delete extra figures
        delete(o);
    end
    yl=cat(1,axH.YLim);
    yl=[min(yl(:,1)) ,max(yl(:,2))];
    
    j=1;
    %Small VR = f(p)
    if tit;title(axH(j),"Average small cluster Void Ratio");end
    axH(j).XLim(1)=0;
    axH(j).YLim=yl;
    ylabel(axH(j),'Void Ratio')
    xlabel(axH(j),'Mean pressure (kPA)')
    fnm="VR_Small_P.png";
    saveas(h(j),fullfile(path,fnm));
    j=j+1;
    
    %Small VR = f(E)
    if tit;title(axH(j),"Average small cluster Void Ratio");end
    axH(j).YLim=yl;
    ylabel(axH(j),'Void Ratio')
    xlabel(axH(j),'Axial Strain')
    fnm="VR_Small_E.png";
    saveas(h(j),fullfile(path,fnm));
    j=j+1;
    
    %Submed VR = f(p)
    if tit;title(axH(j),"Average submedium cluster Void Ratio");end
    axH(j).XLim(1)=0;
    axH(j).YLim=yl;
    ylabel(axH(j),'Void Ratio')
    xlabel(axH(j),'Mean pressure (kPA)')
    fnm="VR_Submedium_P.png";
    saveas(h(j),fullfile(path,fnm));
    j=j+1;
    
    %Submed VR = f(E)
    if tit;title(axH(j),"Average submedium cluster Void Ratio");end
    axH(j).YLim=yl;
    ylabel(axH(j),'Void Ratio')
    xlabel(axH(j),'Axial Strain')
    fnm="VR_Submedium_E.png";
    saveas(h(j),fullfile(path,fnm));
    j=j+1;
    
    %Med VR = f(p)
    if tit;title(axH(j),"Average medium cluster Void Ratio");end
    axH(j).XLim(1)=0;
    axH(j).YLim=yl;
    ylabel(axH(j),'Void Ratio')
    xlabel(axH(j),'Mean pressure (kPA)')
    fnm="VR_Medium_P.png";
    saveas(h(j),fullfile(path,fnm));
    j=j+1;
    
    %Med VR = f(E)
    if tit;title(axH(j),"Average medium cluster Void Ratio");end
    axH(j).YLim=yl;
    ylabel(axH(j),'Void Ratio')
    xlabel(axH(j),'Axial Strain')
    fnm="VR_Medium_E.png";
    saveas(h(j),fullfile(path,fnm));
    j=j+1;
    
    %Larg VR = f(p)
    if tit;title(axH(j),"Large cluster Void Ratio");end
    axH(j).XLim(1)=0;
    axH(j).YLim=yl;
    ylabel(axH(j),'Void Ratio')
    xlabel(axH(j),'Mean pressure (kPA)')
    fnm="VR_Large_P.png";
    saveas(h(j),fullfile(path,fnm));
    j=j+1;
    
    %Larg VR = f(E)
    if tit;title(axH(j),"Large cluster Void Ratio");end
    axH(j).YLim=yl;
    ylabel(axH(j),'Void Ratio')
    xlabel(axH(j),'Axial Strain')
    fnm="VR_Large_E.png";
    saveas(h(j),fullfile(path,fnm));
    
else
     if leg
        legend(axH(1),"Small","Submedium","Medium","Large","Specimen",'location','best')
        legend(axH(2),"Small","Submedium","Medium","Large","Specimen",'location','best')
     end
    j=1;
    %cat VR = f(p)
    if tit;title(axH(j),"Average cluster Void Ratio");end
    ylabel(axH(j),'Void Ratio')
    xlabel(axH(j),'Mean pressure (kPA)')
    fnm="VR_Small_P.png";
    saveas(h(j),fullfile(path,fnm));
    j=j+1;
    
    %cat VR = f(E)
    if tit;title(axH(j),"Average cluster Void Ratio");end
    ylabel(axH(j),'Void Ratio')
    xlabel(axH(j),'Axial Strain')
    fnm="VR_Small_E.png";
    saveas(h(j),fullfile(path,fnm));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% VR per Cl order %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plot VR of the demanded Cluster's orders
lp='';%4:2:20;
if ~isempty(lp)
    nb=numel(lp);
    if nb<8
        C = app.PlotColors;
    else
        C = graphClrCode(nb);%plot colorcode
    end
    path=MakePath(app,'LOOPVRCL');
    g(nb)=figure;axG(nb)=axes(g(nb));hold(axG(nb),'on');
    for i=1:(nb-1)
        g(i)=figure;axG(i)=axes(g(i));hold(axG(i),'on');
    end
    %plot cluster VR evolution in f of Ez and Ezpic
        %pD(i).Results.Mean first colum~order, first line ~strain. 3 dim :
        %1 VR, 2 Z
    for j=1:numel(lp)
        for i=1:size(pD,2)
            line=find(pD(i).Results.Mean(:,1,1)==lp(j),1);
            if isempty(line); continue;end
            opts={'Color',C(i,:)};
            for k=1:numel(pD(i).InfPts.q)
                if strcmpi(x,'p')
                    opts=[opts,{'Pointx',pD(i).InfPts.p(k)}]; %#ok<AGROW>
                    xAx=pD(i).Results.Pressure;
                else
                    opts=[opts,{'Pointx',pD(i).InfPts.ez(k)}]; %#ok<AGROW>
                    xAx=pD(i).Results.Strain;
                end
            end
            %find no zero values
            no0ln=find(pD(i).Results.Mean(line,2:end,1))+1;
            %plot it
            plotMark(app,axG(j),xAx,...
                pD(i).Results.Mean(line,no0ln,1),opts{:})
%             %find Ez for min (e~=0)
%             mine=min(pD(i).Results.Mean(line,...
%                 (pD(i).Results.Mean(line,:,1)>0),1));
%             fd=find(pD(i).Results.Mean(line,:,1)==mine,1,'first');
%             %Obtain the first point where the strain is positive to divide
%             %all and try putting all different consolidation pressures in
%             %the same 'axis'
%             fd=find(pD(i).Results.Mean(line,2:end,1)>0,1,'first')+1;
%             plotMark(app,axG(2*j),pD(i).Results.Mean(1,2:end,1),...
%                 pD(i).Results.Mean(line,2:end,1)/...
%                 pD(i).Results.Mean(line,fd,1),opts{:})
%             %correct the limits
%             xlim(axG(2*j-1),[pD(i).Results.Mean(1,1,1)...
%                 pD(i).Results.Mean(end,1,1)*1.05])
        end
        
        if size(pD,2)>1
            if leg
                legend(axG(j),pD(:).FileName,'location','best')
%                 legend(axG(2*j),pD(:).FileName,'location','best')
            end
        end
        
        %Loop e = f(Ez)
        if tit;title(axG(j),"Cluster "+ lp(j) +" Void Ratio as f(Ez)");end
        if strcmpi(x,'p');axG(j).XLim(1)=0;end
        ylabel(axG(j),'Void Ratio')
        xlabel(axG(j),xlab)
        fnm="VR_Cluster"+lp(j)+".png";
        saveas(g(j),fullfile(path,fnm));
%         
%         %Loop e/e0 = f(Ez/Ez0)
%         if tit;title(axG(2*j),"Cluster "+ lp(j) +" Void Ratio as f(Ez)");end
%         ylabel(axG(2*j),'Void Ratio')
%         xlabel(axG(2*j),xlab)
%         fnm="VR_Cluster"+lp(j)+"_Ratio";
%         saveas(g(2*j),fullfile(path,fnm+png));
    end
end

if isequal(app.LIGGGHTSAnalysisButtonGroup.SelectedObject,app.ExeAllButton)
    delete([f,h]);
    if~isempty(lp);delete(g);end
    if numel(pD)==1;delete(l);end
elseif numel(pD)>2
    delete(f);
% else
%     ax(3).XLim(2)=l1;
%     ax(4).XLim(2)=l2;
end
end
function clStrPlotter(pD,app)
%CLSTRPLOTTER plot the evolution the loop stress during the test

if app.SimType==3
   xlab='Mean pressure (kPA)';
else
   xlab='Axial strain';
end

if isequal(app.CourbePointsSwitch.Value,'On')
    lType="-+";
else
    lType="-";
end
%Check title and legends option
if app.TitlesCB.Value;tit=1;else;tit=0;end
if app.LegendsCB.Value;leg=1;else;leg=0;end

C=app.PlotColors;
set(0,'defaultAxesFontSize',app.FontSizeEF.Value)
nb=3;
f(nb)=figure;ax(nb)=axes(f(nb));hold(ax(nb),'on');
for i=1:(nb-1)
    f(i)=figure;ax(i)=axes(f(i));hold(ax(i),'on');
end
%Plot scatter of data Stress x Cl order
strs=cat(1,pD.Results.Stress{:}); %ClusterOrder - p - q
scatter(ax(1),strs(:,1),strs(:,3))

%plot mean stress per category
str=zeros(numel(pD.Results.Stress),4);
for i=1:numel(pD.Results.Stress)
    t=pD.Results.Stress{i};
    str(i,:)=[mean(t(t(:,1)==4,3)), mean(t(t(:,1)==6,3)),...
        mean(t(t(:,1)>6 & t(:,1)<21,3)), mean(t(t(:,1)>21,3))];
end
strain=pD.Results.Strain;
for i=1:2
    if i==1
        plot(ax(i+1),strain,str(:,1),lType,'Color',C(1,:))
    end
    plot(ax(i+1),strain,str(:,2),lType,'Color',C(2,:))
    plot(ax(i+1),strain,str(:,3),lType,'Color',C(3,:))
    plot(ax(i+1),strain,str(:,4),lType,'Color',C(4,:))
end
if leg
    legend(ax(2),"Small","Submedium","Medium","Large",...
        'location','best')
    legend(ax(3),"Submedium","Medium","Large",...
        'location','best')
end
%titles
if tit
    strg=["Per Cluster stress data";
        "Mean Stress per category";
        "Mean Stress per category"];
    for i=1:nb
        title(ax(i),strg(ceil(i)))
    end
end

path=MakePath(app,'LOOPST');
%Labels and save
fnm=["W2Cl";"W2ClN4";"W2ClPct";];
for i=1:nb
    if i==1
        ylabel(ax(i),'Stress (kPa)')
        xlabel(ax(i),'Order')
    else
        ylabel(ax(i),'Stress (kPa)')
        xlabel(ax(i),xlab)
    end
    f(i).Position(3:4)=app.Figsize;
    saveas(f(i),path+fnm(i)+".png");
end
if isequal(app.LIGGGHTSAnalysisButtonGroup.SelectedObject,app.ExeAllButton)
    delete(f);
end
end
function clTransfPlotter(pD,app)
%CLSTRPLOTTER plot the evolution the Cl4 and Cl6 transformations
% start it
if app.SimType==3
   xlab='Mean pressure (kPA)';
else
   xlab='Axial strain';
end

%matlab default colors to mantain every class in the same clor
C = app.PlotColors;
set(0,'defaultAxesFontSize',app.FontSizeEF.Value)

%Check title and legends option
if app.TitlesCB.Value;tit=1;else;tit=0;end
if app.LegendsCB.Value;leg=1;else;leg=0;end

path=MakePath(app,'LOOPCT');

png='.png';
%prepare variables for execution
ctRes=pD.Results;            %save
cl=zeros(size(ctRes.ClTransf,1),3,4);%cluster data 
clqStn=zeros(size(ctRes.ClTransf,1),3,4);%cluster data, high strain
clpStn=zeros(size(ctRes.ClTransf,1),3,4);%cluster data, high strain
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
    
    %Only high (than ave) dev strain transformations
    %calculate the deviatoric strain of cells
    stn=ctRes.ClTransf{i,2};
    p=(stn(1,1,:)+stn(2,2,:)+stn(3,3,:))/3;
    p=find(p>mean(p));
    q=sqrt(( (stn(2,2,:)-stn(1,1,:)).^2+...	%(sig1-sig2)^2 /2
                    (stn(3,3,:)-stn(1,1,:)).^2+...          %(sig2-sig3)^2 /2
                    (stn(1,1,:)-stn(3,3,:)).^2+...          %(sig3-sig1)^2 /2
                    6*stn(1,2,:).^2+...                     % 3*(sig12)^2
                    6*stn(1,3,:).^2+...                     % 3*(sig13)^2pcq
                    6*stn(2,3,:).^2 )/2);                   % 3*(sig23)^2
    q=find(q>mean(q)); %get id of cells with high dev strain
    %get only results from these
    clpStn(i,:,:)=cat(3,sum(cl4(ismember(r.Cl4(:,end),p),:),1),...
        sum(cl6(ismember(r.Cl6(:,end),p),:),1),...
        sum(cl8(ismember(r.Cl8(cat8,end),p),:),1),...
        sum(cl22(ismember(r.Cl8(~cat8,end),p),:),1));
    %get only results from these
    clqStn(i,:,:)=cat(3,sum(cl4(ismember(r.Cl4(:,end),q),:),1),...
        sum(cl6(ismember(r.Cl6(:,end),q),:),1),...
        sum(cl8(ismember(r.Cl8(cat8,end),q),:),1),...
        sum(cl22(ismember(r.Cl8(~cat8,end),q),:),1));
    
end
plG1=4; %plot group 1 => number of groups of tiled figures
plG2=2; %plot group 2 => number of groups of simple figures
nb=4*(plG1+plG2);
%Create figures
f(nb)=figure;ax(nb,1)=axes(f(nb));hold(ax(nb),'on');
for i=4*(plG1)+1:(nb-1)
    f(i)=figure;ax(i,1)=axes(f(i));hold(ax(i),'on');
end

%Create tiled layouts
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

opts={};
for k=1:numel(pD.InfPts.q)
    if app.SimType==3
        opts=[opts,{'Pointx',pD.InfPts.p(k)}]; %#ok<AGROW>
    else
        opts=[opts,{'Pointx',pD.InfPts.ez(k)}]; %#ok<AGROW>
    end
end

i=1;
%%%%%%%%%%%%%%%% PLOTGROUP 1 transformatio ratio %%%%%%%%%%%%%%%%
%Cl4 => others  & others => Cl4 - - ratio
ax(i,1)=nexttile(tf(i));hold(ax(i,1),'on');
plotMark(app,ax(i,1),ctRes.Strain,cl(:,1,1)./nbCll(:,1),'Color',C(2,:),opts{:});
plotMark(app,ax(i,1),ctRes.Strain,cl(:,2,1)./nbCll(:,1),'Color',C(3,:),opts{:});
plotMark(app,ax(i,1),ctRes.Strain,cl(:,3,1)./nbCll(:,1),'Color',C(4,:),opts{:});
ax(i,2)=nexttile(tf(i));hold(ax(i,2),'on');
plotMark(app,ax(i,2),ctRes.Strain,cl(:,1,2)./nbCll(:,1),'Color',C(2,:),opts{:});
plotMark(app,ax(i,2),ctRes.Strain,cl(:,1,3)./nbCll(:,1),'Color',C(3,:),opts{:});
plotMark(app,ax(i,2),ctRes.Strain,cl(:,1,4)./nbCll(:,1),'Color',C(4,:),opts{:});
i=i+1;
%Cl6 => others & others => Cl6 - - ratio
ax(i,1)=nexttile(tf(i));hold(ax(i,1),'on');
plotMark(app,ax(i,1),ctRes.Strain,cl(:,1,2)./nbCll(:,2),'Color',C(1,:),opts{:});
plotMark(app,ax(i,1),ctRes.Strain,cl(:,2,2)./nbCll(:,2),'Color',C(3,:),opts{:});
plotMark(app,ax(i,1),ctRes.Strain,cl(:,3,2)./nbCll(:,2),'Color',C(4,:),opts{:});
ax(i,2)=nexttile(tf(i));hold(ax(i,2),'on');
plotMark(app,ax(i,2),ctRes.Strain,cl(:,1,1)./nbCll(:,2),'Color',C(1,:),opts{:});
plotMark(app,ax(i,2),ctRes.Strain,cl(:,2,3)./nbCll(:,2),'Color',C(3,:),opts{:});
plotMark(app,ax(i,2),ctRes.Strain,cl(:,2,4)./nbCll(:,2),'Color',C(4,:),opts{:});

i=i+1;
%Cl8 => others & others => Cl8 - - ratio
ax(i,1)=nexttile(tf(i));hold(ax(i,1),'on');
plotMark(app,ax(i,1),ctRes.Strain,cl(:,1,3)./nbCll(:,3),'Color',C(1,:),opts{:});
plotMark(app,ax(i,1),ctRes.Strain,cl(:,2,3)./nbCll(:,3),'Color',C(2,:),opts{:});
plotMark(app,ax(i,1),ctRes.Strain,cl(:,3,3)./nbCll(:,3),'Color',C(4,:),opts{:});
ax(i,2)=nexttile(tf(i));hold(ax(i,2),'on');
plotMark(app,ax(i,2),ctRes.Strain,cl(:,2,1)./nbCll(:,3),'Color',C(1,:),opts{:});
plotMark(app,ax(i,2),ctRes.Strain,cl(:,2,2)./nbCll(:,3),'Color',C(2,:),opts{:});
plotMark(app,ax(i,2),ctRes.Strain,cl(:,3,4)./nbCll(:,3),'Color',C(4,:),opts{:});
i=i+1;
%Cl22 => others & others => Cl22 - - ratio
ax(i,1)=nexttile(tf(i));hold(ax(i,1),'on');
plotMark(app,ax(i,1),ctRes.Strain,cl(:,1,4)./nbCll(:,4),'Color',C(1,:),opts{:});
plotMark(app,ax(i,1),ctRes.Strain,cl(:,2,4)./nbCll(:,4),'Color',C(2,:),opts{:});
plotMark(app,ax(i,1),ctRes.Strain,cl(:,3,4)./nbCll(:,4),'Color',C(3,:),opts{:});
ax(i,2)=nexttile(tf(i));hold(ax(i,2),'on');
plotMark(app,ax(i,2),ctRes.Strain,cl(:,3,1)./nbCll(:,4),'Color',C(1,:),opts{:});
plotMark(app,ax(i,2),ctRes.Strain,cl(:,3,2)./nbCll(:,4),'Color',C(2,:),opts{:});
plotMark(app,ax(i,2),ctRes.Strain,cl(:,3,3)./nbCll(:,4),'Color',C(3,:),opts{:});
i=i+1;

%%%%%%%%%%%%%%%% PLOTGROUP 2 transformation per nb %%%%%%%%%%%%%%%%
%Cl4 => others  & others => Cl4 - - nb
ax(i,1)=nexttile(tf(i));hold(ax(i,1),'on');
plotMark(app,ax(i,1),ctRes.Strain,cl(:,1,1),'Color',C(2,:),opts{:});
plotMark(app,ax(i,1),ctRes.Strain,cl(:,2,1),'Color',C(3,:),opts{:});
plotMark(app,ax(i,1),ctRes.Strain,cl(:,3,1),'Color',C(4,:),opts{:});
ax(i,2)=nexttile(tf(i));hold(ax(i,2),'on');
plotMark(app,ax(i,2),ctRes.Strain,cl(:,1,2),'Color',C(2,:),opts{:});
plotMark(app,ax(i,2),ctRes.Strain,cl(:,1,3),'Color',C(3,:),opts{:});
plotMark(app,ax(i,2),ctRes.Strain,cl(:,1,4),'Color',C(4,:),opts{:});
i=i+1;


%Cl6 => others & others => Cl6 - - nb
ax(i,1)=nexttile(tf(i));hold(ax(i,1),'on');
plotMark(app,ax(i,1),ctRes.Strain,cl(:,1,2),'Color',C(1,:),opts{:});
plotMark(app,ax(i,1),ctRes.Strain,cl(:,2,2),'Color',C(3,:),opts{:});
plotMark(app,ax(i,1),ctRes.Strain,cl(:,3,2),'Color',C(4,:),opts{:});
ax(i,2)=nexttile(tf(i));hold(ax(i,2),'on');
plotMark(app,ax(i,2),ctRes.Strain,cl(:,1,1),'Color',C(1,:),opts{:});
plotMark(app,ax(i,2),ctRes.Strain,cl(:,2,3),'Color',C(3,:),opts{:});
plotMark(app,ax(i,2),ctRes.Strain,cl(:,2,4),'Color',C(4,:),opts{:});
i=i+1;

%Cl8 => others & others => Cl8 - - nb
ax(i,1)=nexttile(tf(i));hold(ax(i,1),'on');
plotMark(app,ax(i,1),ctRes.Strain,cl(:,1,3),'Color',C(1,:),opts{:});
plotMark(app,ax(i,1),ctRes.Strain,cl(:,2,3),'Color',C(2,:),opts{:});
plotMark(app,ax(i,1),ctRes.Strain,cl(:,3,3),'Color',C(4,:),opts{:});
ax(i,2)=nexttile(tf(i));hold(ax(i,2),'on');
plotMark(app,ax(i,2),ctRes.Strain,cl(:,2,1),'Color',C(1,:),opts{:});
plotMark(app,ax(i,2),ctRes.Strain,cl(:,2,2),'Color',C(2,:),opts{:});
plotMark(app,ax(i,2),ctRes.Strain,cl(:,3,4),'Color',C(4,:),opts{:});
i=i+1;

%Cl22 => others & others => Cl22 - - nb
ax(i,1)=nexttile(tf(i));hold(ax(i,1),'on');
plotMark(app,ax(i,1),ctRes.Strain,cl(:,1,4),'Color',C(1,:),opts{:});
plotMark(app,ax(i,1),ctRes.Strain,cl(:,2,4),'Color',C(2,:),opts{:});
plotMark(app,ax(i,1),ctRes.Strain,cl(:,3,4),'Color',C(3,:),opts{:});
ax(i,2)=nexttile(tf(i));hold(ax(i,2),'on');
plotMark(app,ax(i,2),ctRes.Strain,cl(:,3,1),'Color',C(1,:),opts{:});
plotMark(app,ax(i,2),ctRes.Strain,cl(:,3,2),'Color',C(2,:),opts{:});
plotMark(app,ax(i,2),ctRes.Strain,cl(:,3,3),'Color',C(3,:),opts{:});
i=i+1;

%%%%%%%%%%%%%%%% PLOTGROUP 3 - High p strain vaules %%%%%%%%%%%%%%%%
%Cl4 => others  & others => Cl4 - - nb
ax(i,1)=nexttile(tf(i));hold(ax(i,1),'on');
plotMark(app,ax(i,1),ctRes.Strain,clpStn(:,1,1),'Color',C(2,:),opts{:});
plotMark(app,ax(i,1),ctRes.Strain,clpStn(:,2,1),'Color',C(3,:),opts{:});
plotMark(app,ax(i,1),ctRes.Strain,clpStn(:,3,1),'Color',C(4,:),opts{:});
ax(i,2)=nexttile(tf(i)); hold(ax(i,2),'on');
plotMark(app,ax(i,2),ctRes.Strain,clpStn(:,1,2),'Color',C(2,:),opts{:});
plotMark(app,ax(i,2),ctRes.Strain,clpStn(:,1,3),'Color',C(3,:),opts{:});
plotMark(app,ax(i,2),ctRes.Strain,clpStn(:,1,4),'Color',C(4,:),opts{:});
i=i+1;
%cl6 => others & others => cl6 - - nb
ax(i,1)=nexttile(tf(i));hold(ax(i,1),'on');
plotMark(app,ax(i,1),ctRes.Strain,clpStn(:,1,2),'Color',C(1,:),opts{:});
plotMark(app,ax(i,1),ctRes.Strain,clpStn(:,2,2),'Color',C(3,:),opts{:});
plotMark(app,ax(i,1),ctRes.Strain,clpStn(:,3,2),'Color',C(4,:),opts{:});
ax(i,2)=nexttile(tf(i));hold(ax(i,2),'on');
plotMark(app,ax(i,2),ctRes.Strain,clpStn(:,1,1),'Color',C(1,:),opts{:});
plotMark(app,ax(i,2),ctRes.Strain,clpStn(:,2,3),'Color',C(3,:),opts{:});
plotMark(app,ax(i,2),ctRes.Strain,clpStn(:,2,4),'Color',C(4,:),opts{:});
i=i+1;
%cl8 => others & others => cl8 - - nb
ax(i,1)=nexttile(tf(i));hold(ax(i,1),'on');
plotMark(app,ax(i,1),ctRes.Strain,clpStn(:,1,3),'Color',C(1,:),opts{:});
plotMark(app,ax(i,1),ctRes.Strain,clpStn(:,2,3),'Color',C(2,:),opts{:});
plotMark(app,ax(i,1),ctRes.Strain,clpStn(:,3,3),'Color',C(4,:),opts{:});
ax(i,2)=nexttile(tf(i));hold(ax(i,2),'on');
plotMark(app,ax(i,2),ctRes.Strain,clpStn(:,2,1),'Color',C(1,:),opts{:});
plotMark(app,ax(i,2),ctRes.Strain,clpStn(:,2,2),'Color',C(2,:),opts{:});
plotMark(app,ax(i,2),ctRes.Strain,clpStn(:,3,4),'Color',C(4,:),opts{:});
i=i+1;

%cl22 => others & others => cl22 - - nb
ax(i,1)=nexttile(tf(i));hold(ax(i,1),'on');
plotMark(app,ax(i,1),ctRes.Strain,clpStn(:,1,4),'Color',C(1,:),opts{:});
plotMark(app,ax(i,1),ctRes.Strain,clpStn(:,2,4),'Color',C(2,:),opts{:});
plotMark(app,ax(i,1),ctRes.Strain,clpStn(:,3,4),'Color',C(3,:),opts{:});
ax(i,2)=nexttile(tf(i));hold(ax(i,2),'on');
plotMark(app,ax(i,2),ctRes.Strain,clpStn(:,3,1),'Color',C(1,:),opts{:});
plotMark(app,ax(i,2),ctRes.Strain,clpStn(:,3,2),'Color',C(2,:),opts{:});
plotMark(app,ax(i,2),ctRes.Strain,clpStn(:,3,3),'Color',C(3,:),opts{:});
i=i+1;

%%%%%%%%%%%%%%%% PLOTGROUP 4 - High q strain vaules %%%%%%%%%%%%%%%%
%Cl4 => others  & others => Cl4 - - nb
ax(i,1)=nexttile(tf(i));hold(ax(i,1),'on');
plotMark(app,ax(i,1),ctRes.Strain,clqStn(:,1,1),'Color',C(2,:),opts{:});
plotMark(app,ax(i,1),ctRes.Strain,clqStn(:,2,1),'Color',C(3,:),opts{:});
plotMark(app,ax(i,1),ctRes.Strain,clqStn(:,3,1),'Color',C(4,:),opts{:});
ax(i,2)=nexttile(tf(i)); hold(ax(i,2),'on');
plotMark(app,ax(i,2),ctRes.Strain,clqStn(:,1,2),'Color',C(2,:),opts{:});
plotMark(app,ax(i,2),ctRes.Strain,clqStn(:,1,3),'Color',C(3,:),opts{:});
plotMark(app,ax(i,2),ctRes.Strain,clqStn(:,1,4),'Color',C(4,:),opts{:});
i=i+1;
%cl6 => others & others => cl6 - - nb
ax(i,1)=nexttile(tf(i));hold(ax(i,1),'on');
plotMark(app,ax(i,1),ctRes.Strain,clqStn(:,1,2),'Color',C(1,:),opts{:});
plotMark(app,ax(i,1),ctRes.Strain,clqStn(:,2,2),'Color',C(3,:),opts{:});
plotMark(app,ax(i,1),ctRes.Strain,clqStn(:,3,2),'Color',C(4,:),opts{:});
ax(i,2)=nexttile(tf(i));hold(ax(i,2),'on');
plotMark(app,ax(i,2),ctRes.Strain,clqStn(:,1,1),'Color',C(1,:),opts{:});
plotMark(app,ax(i,2),ctRes.Strain,clqStn(:,2,3),'Color',C(3,:),opts{:});
plotMark(app,ax(i,2),ctRes.Strain,clqStn(:,2,4),'Color',C(4,:),opts{:});
i=i+1;
%cl8 => others & others => cl8 - - nb
ax(i,1)=nexttile(tf(i));hold(ax(i,1),'on');
plotMark(app,ax(i,1),ctRes.Strain,clqStn(:,1,3),'Color',C(1,:),opts{:});
plotMark(app,ax(i,1),ctRes.Strain,clqStn(:,2,3),'Color',C(2,:),opts{:});
plotMark(app,ax(i,1),ctRes.Strain,clqStn(:,3,3),'Color',C(4,:),opts{:});
ax(i,2)=nexttile(tf(i));hold(ax(i,2),'on');
plotMark(app,ax(i,2),ctRes.Strain,clqStn(:,2,1),'Color',C(1,:),opts{:});
plotMark(app,ax(i,2),ctRes.Strain,clqStn(:,2,2),'Color',C(2,:),opts{:});
plotMark(app,ax(i,2),ctRes.Strain,clqStn(:,3,4),'Color',C(4,:),opts{:});
i=i+1;

%cl22 => others & others => cl22 - - nb
ax(i,1)=nexttile(tf(i));hold(ax(i,1),'on');
plotMark(app,ax(i,1),ctRes.Strain,clqStn(:,1,4),'Color',C(1,:),opts{:});
plotMark(app,ax(i,1),ctRes.Strain,clqStn(:,2,4),'Color',C(2,:),opts{:});
plotMark(app,ax(i,1),ctRes.Strain,clqStn(:,3,4),'Color',C(3,:),opts{:});
ax(i,2)=nexttile(tf(i));hold(ax(i,2),'on');
plotMark(app,ax(i,2),ctRes.Strain,clqStn(:,3,1),'Color',C(1,:),opts{:});
plotMark(app,ax(i,2),ctRes.Strain,clqStn(:,3,2),'Color',C(2,:),opts{:});
plotMark(app,ax(i,2),ctRes.Strain,clqStn(:,3,3),'Color',C(3,:),opts{:});
i=i+1;

%%%%%%%%%%%%%%%% PLOTGROUP 5 - Delta values %%%%%%%%%%%%%%%%
ord=ctRes.NbCell(3:end,2:end);
ord=ord(:,2:end)-ord(:,1:end-1);
%Delta Cl4
plotMark(app,ax(i,1),ctRes.Strain,(sum(cl(:,1,2:4),3)...
    -sum(cl(:,1:3,1),2)),'Color',C(1,:),opts{:});
plotMark(app,ax(i,1),ctRes.Strain,[0 ord(1,:)],'Color',C(2,:),opts{:});
i=i+1;
%Delta Cl6
plotMark(app,ax(i,1),ctRes.Strain,(cl(:,1,1)+sum(cl(:,2,3:4),3)...
    -sum(cl(:,1:3,2),2)),'Color',C(1,:),opts{:});
plotMark(app,ax(i,1),ctRes.Strain,[0 ord(2,:)],'Color',C(2,:),opts{:});
i=i+1;
%Delta Cl8
plotMark(app,ax(i,1),ctRes.Strain,(cl(:,3,4)+sum(cl(:,2,1:2),3)...
    -sum(cl(:,1:3,3),2)),'Color',C(1,:),opts{:});
plotMark(app,ax(i,1),ctRes.Strain,[0 sum(ord(3:9,:),1)],'Color',C(2,:),opts{:});
i=i+1;
%Delta Cl22
plotMark(app,ax(i,1),ctRes.Strain,(sum(cl(:,3,1:3),3)...
    -sum(cl(:,1:3,4),2)),'Color',C(1,:),opts{:});
plotMark(app,ax(i,1),ctRes.Strain,[0 sum(ord(10:end,:),1)],'Color',C(2,:),opts{:});
i=i+1;
%%%%%%%%%%%%%%%% PLOTGROUP 6  %%%%%%%%%%%%%%%%
% (others => Cl4) -(Cl4 => others)
plotMark(app,ax(i,1),ctRes.Strain,cl(:,1,2)-cl(:,1,1),'Color',C(2,:),opts{:});
plotMark(app,ax(i,1),ctRes.Strain,cl(:,1,3)-cl(:,2,1),'Color',C(3,:),opts{:});
plotMark(app,ax(i,1),ctRes.Strain,cl(:,1,4)-cl(:,3,1),'Color',C(4,:),opts{:});
i=i+1;
% (others => Cl6) - (Cl6 => others)
plotMark(app,ax(i,1),ctRes.Strain,cl(:,1,1)-cl(:,1,2),'Color',C(1,:),opts{:});
plotMark(app,ax(i,1),ctRes.Strain,cl(:,2,3)-cl(:,2,2),'Color',C(3,:),opts{:});
plotMark(app,ax(i,1),ctRes.Strain,cl(:,2,4)-cl(:,3,2),'Color',C(4,:),opts{:});
i=i+1;
% (others => Cl8) - (Cl8 => others)
plotMark(app,ax(i,1),ctRes.Strain,cl(:,2,1)-cl(:,1,3),'Color',C(1,:),opts{:});
plotMark(app,ax(i,1),ctRes.Strain,cl(:,2,2)-cl(:,2,3),'Color',C(2,:),opts{:});
plotMark(app,ax(i,1),ctRes.Strain,cl(:,3,4)-cl(:,3,3),'Color',C(4,:),opts{:});
i=i+1;
% (others => Cl22) - (Cl22 => others)
plotMark(app,ax(i,1),ctRes.Strain,cl(:,3,1)-cl(:,1,4),'Color',C(1,:),opts{:});
plotMark(app,ax(i,1),ctRes.Strain,cl(:,3,2)-cl(:,2,4),'Color',C(2,:),opts{:});
plotMark(app,ax(i,1),ctRes.Strain,cl(:,3,3)-cl(:,3,4),'Color',C(3,:),opts{:});

%%%%%%%%%%%%%%%% Titles and stuff %%%%%%%%%%%%%%%%
legs={["Cl6","Cl8-20","Cl22+"];
    ["Cl4","Cl8-20","Cl22+"];
    ["Cl4","Cl6","Cl22+"];
    ["Cl4","Cl6","Cl8-20"];
    ["Transformation","Delta Numbers"]};
ylab=["Ratio of Events";
    "Number of Events";
    "Number of High Vol Strain Events";
    "Number of High Dev Strain Events";
    "Number of Events";
    "Number of Events"];
titles=["Order 4 to Others","Others to 4";
    "Order 6 to Others","Others to 6";
    "Order 8-20 to Others","Others to 8-20";
    "Order 22+ to Others","Others to 22+";
    "Order 4 to Others","Others to 4";
    "Order 6 to Others","Others to 6";
    "Order 8-20 to Others","Others to 8-20";
    "Order 22+ to Others","Others to 22+";
    "Order 4 to Others","Others to 4";
    "Order 6 to Others","Others to 6";
    "Order 8-20 to Others","Others to 8-20";
    "Order 22+ to Others","Others to 22+";
    "Order 4 to Others","Others to 4";
    "Order 6 to Others","Others to 6";
    "Order 8-20 to Others","Others to 8-20";
    "Order 22+ to Others","Others to 22+";
    "Delta Curb Order 4","";
    "Delta Curb Order 6","";
    "Delta Curb Order 8-20","";
    "Delta Curb Order 22+","";
    "Small Category Transf","";
    "Submedium Category Transf","";
    "Medium Category Transf","";
    "Large Category Transf",""];
fnm=["Transf_Base_O4";"Transf_Base_O6";
    "Transf_Base_O8";"Transf_Base_O22";
    "Transf_Nb_O4";"Transf_Nb_O6";
    "Transf_Nb_O8";"Transf_Nb_O22";
    "Transf_p_O4";"Transf_p_O6";
    "Transf_p_O8";"Transf_p_O22";
    "Transf_q_O4";"Transf_q_O6";
    "Transf_q_O8";"Transf_q_O22";
    "Transf_Curb_O4";"Transf_Curb_O6";
    "Transf_Curb_O8";"Transf_Curb_O22";
    "Transf_Tot_O4";"Transf_Tot_O6";
    "Transf_Tot_O8";"Transf_Tot_O22"];

%%%%%%%%%%%%%%%% Save files %%%%%%%%%%%%%%%%
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
        elseif i>4*(plG1) && i<=5*(plG1)
            legend(ax(i,1),legs{5},'Location','best');
        else
            legend(ax(i,1),legs{i-5*(plG1)},'Location','best');
        end
    end
    %add a horizontal line on zero
    yl=yline(ax(i,1),0,':k');
    set(get(get(yl,'Annotation'),'LegendInformation'),...
        'IconDisplayStyle','off');
    if isprop(ax(i,2),'Xlim')   %if tiled layout add on the other too
        yl=yline(ax(i,2),0,':k');
        set(get(get(yl,'Annotation'),'LegendInformation'),...
            'IconDisplayStyle','off');
    end
    if tit;title(ax(i,1),titles(i,1));end
    ylabel(ax(i,1),ylab(ceil(i/4)))
    xlabel(ax(i,1),xlab)
    if isprop(ax(i,2),'Xlim')
        if tit;title(ax(i,2),titles(i,2));end
        ylabel(ax(i,2),ylab(ceil(i/4)))
        xlabel(ax(i,2),xlab)
        %put them in the same vertical scale
        mx=max(ax(i,1).YLim(2),ax(i,2).YLim(2)); 
        ax(i,1).YLim(2)=mx;
        ax(i,2).YLim(2)=mx;
    end
    saveas(f(i),string(path)+fnm(i)+png);
end

%Put paired graphs in the same vertical scale
for i=1:4
   
end
if isequal(app.LIGGGHTSAnalysisButtonGroup.SelectedObject,...
        app.ExeAllButton)
    delete(f);
end
end
%support functions
function cl=clIdentifier(sc)
%CLIDENTIFIER identify ot which cluster each cell belongs to
% Function used on cluster transformaiton calculation. For Cluster identified,
% the cell it is formed from will be identified. Thus an object containing
% 3 lists of cell IDs will be returned with which cluster each cell belongs
% to. This object will have 3 properties with diferent cluster categories
% (4,6 and 8+) because of the different procedures used to identify them.
DT=sort(sc.DelaunayT(:,:),2);
ord=cat(1,sc.Loops.Order);
%cl4=sc.GoodCells(~ismember(sc.GoodCells,cat(2,sc.Loops.sCells)')); %cl4 cells id
cl4=sc.Clt4.sCells;
cl.cl4Id=DT(cl4,:);                 %matrix of cl4 grains
cl6=find(ord==6); %ID of clusters 6
cl.cl6Id=sort(cat(2,sc.Loops(cl6).Grains),1)';%matrix of cl6 grainsIds
c6c=cell(numel(cl6),1);             %cl6 cells
[c6c{:}]=deal(sc.Loops(cl6).sCells);
cl.cl6cells=c6c;                    %cell array containing cl6 cells

cl.cl6=cl6; %cl6 Clusters Ids

%do the same for the others
cl8=find(ord~=6);
nbCl8=numel(cl8); %nb of clusters 8+
mxSz=max(cat(1,sc.Loops.Size));
cl8Id=NaN(nbCl8,mxSz); %create base for saving grains ID
c8c=cell(nbCl8,1); %create base for saving cells IDs
for i=1:nbCl8
    cl8Id(i,1:sc.Loops(cl8(i)).Size)=sc.Loops(cl8(i)).Grains;
    c8c{i}=sc.Loops(cl8(i)).sCells;
end
%save properties
cl.cl8Id=sort(cl8Id,2); %cl8grIds
cl.cl8cells=c8c;        %clcells
cl.cl8Order=ord(cl8);   %cl8 order
cl.cl8=cl8;             %cl8 cluster Ids
cl.srtDT=DT;            %Delaunay triangulation

end
function [X,Y,map]=pointDensity(XY,dX,dY)
%POINTDENSITY calculate a density map centered arround points XY
% - XY is a Nx2 matrix containign the data
% - dX is the interval between each point of the first column of XY
% - dY is interval between each point of the second column of XY
% Returns the necessary values to plot the point density surface using the
% following commands :
%   surface(X,Y,map)
%   colormap(flip(autumn))

%Find out minimum and maximum values of XY and create two vectors following
%the demanded interval
mx=max(XY,[],1);
mn=min(XY,[],1);
X=mn(1):dX:mx(1)+dX;
Y=mn(2):dY:mx(2)+dY;

% Create the map by checking the number of points inside each rectangular
% of size dX dY centered at X,Y
map=NaN(numel(Y),numel(X));
for i=1:(numel(X)-1)
    %check X
    chk=XY(:,1)>(X(i)-dX/2) & XY(:,1)<=(X(i)+dX/2);
    for j=1:(numel(Y)-1)
        %check Y
        v=sum(XY(chk,2)>(Y(j)-dY/2) & XY(chk,2)<=(Y(j)+dY/2));
        if v>0
            %add a 0 to surrounding positions to validate 
            map(j:j+1,i:i+1)=0; 
            map(j,i)=v;
        end 
        %if there is no other point above the check, stop checking
        if sum(XY(chk,2)>(Y(j)-dY/2))==0;break;end
    end
end

%divide the map by the maximum of each column to obtain a 0 to 1 scale per
%X value. Check for NaN values before and after to be sure they are not
%changed. The 0s are necessary to plot the graph correclty
chk=isnan(map); %save previous NaN positions
map=map./(ones(numel(Y),1)*max(map,[],1));
map(logical(isnan(map)-chk))=0; %if any other NaN is created turn into 0

%remove a half interval from each to center each square around X and Y
%values.
X=X-dX/2;
Y=Y-dY/2;
end
%{
function [grid,x,y]=distGraph(data,C1,C2,C3,C4,type)
%DISTGRAPH creates a point density graph to be plot
% This function is then followed by contourf, colormap and colorbar to
% create a density graph. It takes the imput of 4 variables to configurate
% the map. A grid is created of rectangular pieces. Then the number of
% points inside each piece will be accounted. To finish a measure of
% influence of neighbor squares will be accounted.
%  - C1 controls the width X of the subdivisions
%  - C2 controls the height Y of the subdivisions
%  - C3 controls the influence of the rect value in the sum
%  - C4 controls the influence of the neighbor rect value in the sum

% Size vs Order constants : 0.5,1,1.6,0.05
% Def vs Order : 0.5,100,1.6,0.05
% Def vs Size : 1,100,1.6,0.05
%Count stances of unique data
[u,~,cnt]=unique(data,'rows');  %get unique rows and pos
data=[u,accumarray(cnt,1)];     %add the count of how many as third column

%Prepare values
minVx=min(data(:,1));
maxVx=max(data(:,1));                   
divVx=ceil(maxVx-minVx)*C1;
intVx=(maxVx-minVx)/(divVx);

minVy=min(data(:,2));           %minimal
maxVy=max(data(:,2));           %maximal
divVy=ceil(maxVy-minVy)*C2;     %nb of divisions 
intVy=(maxVy-minVy)/(divVy);	%interval between divisions

%sum the influence of each point that is inside each of the small squares
rect=zeros(divVx,divVy);
for i=1:divVx
    %get only the values of the correct order
    for j=1:divVy
        %create the representative rectangle
        polyV=[minVx+intVx*(i-1),minVy+intVy*(j-1);
            minVx+intVx*(i-1),minVy+intVy*(j);
            minVx+intVx*(i),minVy+intVy*(j);
            minVx+intVx*(i),minVy+intVy*(j-1);
            minVx+intVx*(i-1),minVy+intVy*(j-1)];
        test=inpolygon(data(:,1),data(:,2),polyV(:,1),polyV(:,2));
        rect(i,j)=sum(data(test,3));
    end 
end
%calculate the 'Real' value taking a bit of influence of neighbors squares
grid=zeros(divVx,divVy);
for i=1:divVx
    k=0;
    for j=1:divVy
        v=C3*rect(i,j);
        %if v==0;continue;end
        if i~=1;v=v+C4*rect(i-1,j);k=k+1;end
        if i~=divVx;v=v+C4*rect(i+1,j);k=k+1;end
        if j~=1;v=v+C4*rect(i,j-1);k=k+1;end
        if j~=divVy;v=v+C4*rect(i,j+1);k=k+1;end
        if i~=1 && j~=1;v=v+C4*rect(i-1,j-1);k=k+1;end
        if i~=1 && j~=divVy;v=v+C4*rect(i-1,j+1);k=k+1;end
        if i~=divVx && j~=1;v=v+C4*rect(i+1,j-1);k=k+1;end
        if i~=divVx && j~=divVy;v=v+C4*rect(i+1,j+1);k=k+1;end
        v=v/(C3+k*C4);
        if v>0. && v<1; v=1;end %for points that are isolated make sure they are still 1
        grid(i,j)=v;
    end
end
y=minVy:intVy:maxVy;
x=minVx:intVx:maxVx;

if isequal(type,'pcolor')
    x=[x(1)-3*intVx/2, x-intVx/2];
    y=[y,y(end)+intVy];
    grid(grid<1)=NaN;
    grid=[ NaN(1,size(grid,2)+2);...
        NaN(size(grid,1),1) grid NaN(size(grid,1),1);...
        NaN(1,size(grid,2)+2)];
else
    y=y(1:(end-1));
    x=x(1:(end-1))-intVx/2;
    grid(grid<1)=0;
end
grid=log10(grid)';
end
%}
%{
function prc=percentile(file,p)
%PERCENTILE calculates cl order representing p percentile of the Loops data
% This function will return for each column of the "file" entered through
% the loopsPlotter function a value representing the closest order in which
% the order is larger than p percet of the total orders found.

O=file(2:end,1); %order array
%make a cumulative sum of all the loop quantites and then transpos it to
%make it easier to work
f=cumsum(file(2:end,2:end),2)';
%transform p in decimal if it is in percentage
if p>=1;p=p/100;end 
%get a percentage matrix
pct=(f./(ones(size(f,1),1)*f(end,:)));
%find the line of the first occurence a value >=p.
[~,li]=max(pct>=p,[],1);
%create a vector with the position of the first value larger than p and the
%previous one to check witch one is closest to p
oPos=li+(0:size(f,2)-1)*size(f,1); %transform line in overall position
oPos=[pct(oPos-1);pct(oPos)];      % values of pos-1 and pos
[~,oPos]=min(abs(oPos-p),[],1); %check if the minimum is in the line 1 or 2

%return order located in pos or pos-1 depending where the value is closest.
prc=O(li-(oPos-2));
end
%}
%{
OLD VR SCATTER ALL DATA
for j=1:nb
    % FIRST SCATTER j=1 => VR
    % SECOND SCATTER j=1 => Z
    VRZ=double.empty(0,2);   %[Ord VR Z]
    VRZ4=double.empty(0,2);  %[Ord VR Z]
    order=NaN(5000,size(pD,2));   %Ord x numel(pD)
    minVRZ=NaN(5000,size(pD,2));   %min x numel(pD)
    maxVRZ=NaN(5000,size(pD,2));   %max x numel(pD)
    aveVRZ=zeros(5000,size(pD,2)); %mean x numel(pD)
    nbrVRZ=zeros(5000,size(pD,2)); %nb x numel(pD)
    % all points plot
    for i=1:size(pD,2)
        vrz=cat(1,pD(i).Results.ClusterVRZ{:,1});
        if j==2;vrz4=cat(1,pD(i).Results.ClusterVRZ{:,2});end
        mx= groupsummary(vrz(:,j+1), vrz(:,1), 'max');
        mn= groupsummary(vrz(:,j+1), vrz(:,1), 'min');
        [av,od,nb]= groupsummary(vrz(:,j+1), vrz(:,1), 'mean');
        order(od/2,i)=od;
        minVRZ(od/2,i)=mn;
        maxVRZ(od/2,i)=mx;
        aveVRZ(od/2,i)=av;
        nbrVRZ(od/2,i)=nb;
        VRZ=unique([VRZ;vrz(:,[1,j+1])],'rows');
        if j==2;VRZ4=unique([VRZ4;vrz4(:,[1,j+1])],'rows');end
    end
    if j==2
        %do not add the VR data to 
        VRZ=[VRZ;VRZ4]; %#ok<AGROW>
    end
    
    %Transform mmm into a column vector removing exces values (second dim)
    if numel(pD)>1
        order=max(order,[],2);
        minVRZ=min(minVRZ,[],2);
        maxVRZ=max(maxVRZ,[],2);
        aveVRZ=sum(aveVRZ.*nbrVRZ,2)./sum(nbrVRZ,2); 
    end
    
    %remove exces lines from mmm (full 0 lines)
    order=order(order>0);
    minVRZ=minVRZ(order/2,:);maxVRZ=maxVRZ(order/2,:);
    aveVRZ=aveVRZ(order/2,:);
    %plot
    scatter(ax(j),VRZ(:,1),VRZ(:,2),'x','MarkerEdgeColor',[0.6 0.6 0.6]) %scatter all pts
    plot(ax(j),order(aveVRZ>0),aveVRZ(aveVRZ>0),"k"+lType,'LineWidth',lw)   %plot mean line

    %check for places where there is no dif between vals and average for a mean
    %value to be calculated with values around it
    chk=find((aveVRZ-minVRZ)./aveVRZ<=0.01 | (maxVRZ-aveVRZ)./aveVRZ<=0.01);
    m=minVRZ;mx=maxVRZ;
    for l=1:numel(chk)
        if chk(l)<=3
            m(chk(l))=min(minVRZ(1:7));
            mx(chk(l))=max(maxVRZ(1:7));
        elseif chk(l)>=(numel(minVRZ)-20)
            m(chk(l))=min(minVRZ(end-20:end));
            mx(chk(l))=max(maxVRZ(end-20:end));
        else
            m(chk(l))=min(minVRZ(chk(l)-10:chk(l)+10));
            mx(chk(l))=max(maxVRZ(chk(l)-10:chk(l)+10));
        end
    end
    minVRZ=m;maxVRZ=mx;
    %modify them, fit values
    ftmax=fit(order,maxVRZ,fittype({'log(x)','x','1'}));
    ftmin=fit(order,minVRZ,fittype({'log(x)','x','1'}));
    axes(ax(j))
    p=plot(ftmax,'-b',order,maxVRZ); %plot max line
    p(2).LineWidth=lw;
    p=plot(ftmin,'-b',order,minVRZ); %plot min line
    p(2).LineWidth=lw;
    %remove the generated pts
    delete(ax(j).Children([2,4])) %remove points
    %remove the autogenerated legends
    hLeg=findobj(f(j),'type','legend');
    set(hLeg,'visible','off')

end
%plot save
if tit;title(ax(1),"Void Ratio per Order");end
ylabel(ax(1),'Void Ratio')
xlabel(ax(1),'Cluster Order')
if max(order)>200
    L=ax(1).XLim;
    ax(1).XLim=[0,200];
    fnm="Cluster_VR_LowO.png";
    saveas(f(1),fullfile(path,fnm));
    ax(1).XLim=L;
    %f(1).Position(3)=2*f(1).Position(3);
end
fnm="Cluster_VR.png";
saveas(f(1),fullfile(path,fnm));

if tit;title(ax(2),"Coordination per Order");end
ylabel(ax(2),'Coordination')
xlabel(ax(2),'Cluster Order')
if max(order)>200
    L=ax(2).XLim;
    ax(2).XLim=[0,200];
    fnm="Cluster_Z_LowO.png";
    saveas(f(2),fullfile(path,fnm));
    ax(2).XLim=L;
    %f(2).Position(3)=2*f(2).Position(3);
end
fnm="Cluster_Z.png";
saveas(f(2),fullfile(path,fnm));
%}

