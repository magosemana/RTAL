function exe_Anisotropy(PD,app)
%ANISOTROPY Calculates the variation in the fabric anisotropie and the mean 
%coordination number 
%	This function does nothing particularly impressive. It will read the
%	data created by the class 'grains', calculate the eignenvalues of the
%	fabric tensor and plot the important values.

%Check load
if isequal(PD,'LOAD');AniLoad(app);return;end

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

%Check prefix
if isempty(PD);prefix="total";
else;prefix="partial";
end
%Check subdivision and do necessary calculations
if app.SubdivisionButton.Value
    prefix="subd";mode='SUB';
    %read partial data variables
    lMax=PD.SubLin;	
    cMax=PD.SubCol;	    
    nbSub=lMax*cMax;        %total number of extra calculations per file
    %Creating matrix variable storing results
    % contain Step Z Ag Am NbCtctTot NbCtctVert NbCtctHorz
    Results=zeros(nbFiles,7,nbSub); 
else
    nbSub=1;mode='NORMAL';
    %Creating matrix variable storing results
    Results=zeros(nbFiles,9);
end

%Turn on calculation pannel
app=CalcPanel(app,'',nbFiles,'Starting calculation','on');
for i=1:nbFiles
    if getappdata(app.CalculatingPanel,'Canceling')
        CalcPanel(app,'','','','off');
        warndlg('Calculation was canceled');return
    end
    step=stepArray(i);
    app=CalcPanel(app,i,nbFiles,step);
    [gr,PD]=grains('CONTACT',step,PD,app);
    if isempty(gr.FabricTensor) && isempty(gr.PGFabricTensor)
        CalcPanel(app,'','','','off');
        warndlg("Fabric tensor is empty on step "+ step);return;
    end
    %no subdvision calculation
    if app.SubdivisionButton.Value==0 
        Results(i,1)=step;
        FT=gr.FabricTensor; %geometric, force
        e = [eig(FT(:,:,1)),eig(FT(:,:,2)),eig(FT(:,:,3)),eig(FT(:,:,4))];
        e=sort(e,1);
        %coord number
        Results(i,2)=mean(gr.ContactNubr);
        %fabric anisotropy
        if app.Bool3D
            Results(i,3:6)=sqrt(((e(1,:)-e(2,:)).^2+...
                (e(2,:)-e(3,:)).^2+(e(3,:)-e(1,:)).^2)/2)';
        else
            Results(i,3:6)=(e(2,:)-e(1,:))';
        end
        %Contact nb and directions
        Results(i,7:9)=gr.ContactsDir;
    else 
        %subdvision calculation
        pgFT=gr.PGFabricTensor;
        for j=1:nbSub
            subCtcts=gr.ContactNubr(PD.SubGrains{j});
            FTsub=sum(pgFT(:,:,PD.SubGrains{j},:),3);
            FTsub=permute(FTsub,[4,1,2,3])/sum(subCtcts);
            Results(i,2,j)=mean(subCtcts); % Z
            %fabric anisotropy
            if app.Bool3D
                Results(i,3:4,j)=sqrt((FTsub(:,1,1)-FTsub(:,2,2)).^2+...  %(sig1-sig2)^2
                    (FTsub(:,2,2)-FTsub(:,3,3)).^2+...  %(sig2-sig3)^2
                    (FTsub(:,3,3)-FTsub(:,1,1)).^2+...  %(sig3-sig1)^2
                    6*FTsub(:,1,2).^2 +...              % 2*(sig12)^2
                    6*FTsub(:,1,3).^2 +...              % 2*(sig13)^2
                    6*FTsub(:,2,3).^2 )';               % 2*(sig23)^2
            else
                Results(i,3:4,j)=sqrt((...
                    (FTsub(:,1,1)-FTsub(:,2,2)).^2+...  %(sig1-sig2)^2
                    6*FTsub(:,1,2).*FTsub(:,1,2) ))';	% 2*(sig12)^2';
            end
            %Contact nb and directions
            Results(i,5:7,j)=sum(gr.ContactsDir(PD.SubGrains{j},:),1);
        end
        Results(i,1,:)=step*ones(1,1,nbSub);
    end
end
CalcPanel(app,i+1,nbFiles,'','off');

%Calculate the deformation - the last value will be the consolidation step
%as strain
infP=[0,0,0,0];
%take mean pressure if Qcst else strain
if app.SimType==3
    stress = extStress(app.TrialData,Results(:,1,1),app);
    xAxis=stress(:,end);
else
    strain = extStrains(app.TrialData,Results(:,1,1),N1,app);
    xAxis=strain(:,end); 
end 
consoStrain = extStrains(app.TrialData,app.ConsoStep,N1,app);

%Transform Timesteps into Time
Results(:,1,:)=Results(:,1,:)*app.TimeStep;

%Start writing file - some informations about the simulation will also be
%added to the file
path=MakePath(app,'ANI');
fnm="Anisotropy"+N1+"to"+N2+"int"+interval+string(app.FormatEF.Value);
fid = fopen(path+fnm, 'w');
fprintf(fid,'## Evolution of Anisotropy and Coordination Number##\n');
fprintf(fid,'Simulation | Important | Values\n');
fprintf(fid,'%d|%d|%d\n',N1,N2,interval);
fprintf(fid,'%d|%d|%d\n',app.Bool3D,app.checkPiston,app.TimeStep);
fprintf(fid,'%d|%d|%d\n',consoStrain(end),app.ConsoStep*app.TimeStep,infP(1));
fprintf(fid,'%d|%d|%d\n',infP(2:4));
    
if app.SubdivisionButton.Value==0
    fprintf(fid, '-1|-1|-1\n');
    fprintf(fid,['From step %d to %d with interval %d - '...
        char(prefix) '\n'],N1,N2,interval);
    fprintf(fid, ['Step | Coordination | AnisotropyG | AnisotropyM |'...
       'AnisotropyFn | AnisotropyFt | NbContacts | NbcHorizontal | '...
       'NbcVertical | StrainZ\n']);
    Results=[Results xAxis];
    fprintf(fid,'%d|%d|%d|%d|%d|%d|%d|%d|%d|%d\n',Results');
    %Create plotData object
    pD = plotData("Normal",Results,app,prefix,consoStrain(end));
else
    fprintf(fid, '-1|%d|%d\n',lMax,cMax);
    fprintf(fid, ['From step %d to %d, with %d subdivisions - %d Lines'...
        ' %d Columns\n'],N1, N2,nbSub,lMax,cMax);
    fprintf(fid, ['Step | Coordination | AnisotropyG | AnisotropyF |'...
       ' NbContacts | NbcHorizontal | NbcVertical | StrainZ\n']);
    Results=[Results,xAxis.*ones(nbFiles,1,nbSub)];
    %for the subdivision we will separate the results by a line of NaN
    for i=1:nbSub
        fprintf(fid,'%d|%d|%d|%d|%d|%d|%d|%d\n',Results(:,:,i)');
        fprintf(fid,'NaN|NaN|NaN|NaN|NaN|NaN|NaN|NaN\n');
    end
    %Create plotData object
    pD = plotData("Normal",Results,app,prefix,consoStrain(end),lMax,cMax);
end
fclose(fid);
AniPlotter(app,mode,pD);
end
function AniLoad(app)
%Load file
pD=FileLoader(app,'ANI');
if isempty(pD);return;end
%read prefix
if size(pD,2)>1
    mode='MULTI';
else
    if pD.SubL<0
        mode='NORMAL';
    else
        mode='SUB';
        %Reformat file, transform the file matrix, from a big 2D matrix
        %separated by NaN values into a 3D matrix
        file=pD.Results;
        f=find(sum(isnan(file),2));
        fl=file(1:f(1)-1,:);
        for i=2:size(f,1)
            fl=cat(3,fl,file(f(i-1)+1:f(i)-1,:));
        end
        pD.Results=fl;
    end
end
%Plot
AniPlotter(app,mode,pD);
end
function AniPlotter(app,mode,pD)
%If subdivision, start SubPlotter
switch upper(mode)
    case 'SUB'
        subV=[pD.SubL,pD.SubC];
        labels=strings(1,subV(1)*subV(2));
        for l=1:subV(1)
            for c=1:subV(2)
                labels((l-1)*subV(2)+c)="Sub l"+l+" c"+c;
            end
        end
        SubPlotter(app,'ANI',pD.Results,labels,pD.ConsoStrain,subV);
        return;
end

if app.SimType==3
   xlab='Mean pressure (kPA)';
else
   xlab='Axial strain';
end


%Check title and legends option
if app.TitlesCB.Value;tit=1;else;tit=0;end
if app.LegendsCB.Value;leg=1;else;leg=0;end

%Create Figures and Axis
if isequal(upper(mode),'MULTI');nb=8;else;nb=6;end
f(nb)=figure;ax(nb)=axes(f(nb));hold(ax(nb),'on');
for i=1:(nb-1)
    f(i)=figure;ax(i)=axes(f(i));hold(ax(i),'on'); %#ok<LAXES>
end

path=MakePath(app,'ANI');
png=".png";
C=app.PlotColors;
switch upper(mode)
    case 'MULTI'   
        for i=1:size(pD,2)
            fl=pD(i).Results;j=1;
            %coordination
            plotMark(app,ax(j),fl(:,end),fl(:,2),'Color',C(i,:));j=j+1;
            %Anisotropy
            plotMark(app,ax(j),fl(:,end),fl(:,3),'Color',C(i,:));j=j+1;%geo
            plotMark(app,ax(j),fl(:,end),fl(:,4),'Color',C(i,:));j=j+1;%mec
            plotMark(app,ax(j),fl(:,end),fl(:,5),'Color',C(i,:));j=j+1;%normal F
            plotMark(app,ax(j),fl(:,end),fl(:,6),'Color',C(i,:));j=j+1;%tang F
            plotMark(app,ax(j),fl(:,end),fl(:,7),'Color',C(i,:));j=j+1;%tot ct
            plotMark(app,ax(j),fl(:,end),fl(:,8),'Color',C(i,:));j=j+1;%ver ct
            plotMark(app,ax(j),fl(:,end),fl(:,9),'Color',C(i,:))       %hor ct
        end
        if leg
            for j=1:nb
                legend(ax(j),pD.FileName,'location','eastoutside')
            end
        end
        prefix=pD.Prefix;j=1;
        %saveplots
        %plot1 Z=f(Ez)
        if tit;title(ax(j),'Variation of the Coodination Number');end
        ylabel(ax(j),'Coordination Number (Z)')
        xlabel(ax(j),xlab)
        fnm=prefix+"CoordNumber";
        saveas(f(j),fullfile(path,fnm+png));
        j=j+1;
        
        %plot2 Ag=f(Ez)
        if tit;title(ax(j),'Evolution of Geometric Anisotropy');end
        ylabel(ax(j),'Anisotropy')
        xlabel(ax(j),xlab)
        fnm=prefix+"GeoAnisotropy";
        saveas(f(j),fullfile(path,fnm+png));
        j=j+1;
        
        %plot3 Am=f(t)
        if tit;title(ax(j),'Evolution of Mechanical Anisotropy');end
        ylabel(ax(j),'Anisotropy')
        xlabel(ax(j),xlab)
        fnm=prefix+"MechAnisotropy";
        saveas(f(j),fullfile(path,fnm+png));
        j=j+1;
        
        %plot3 An et At=f(t)
        if tit;title(ax(j),'Evolution of Normal Anisotropy');end
        ylabel(ax(j),'Anisotropy')
        xlabel(ax(j),xlab)
        fnm=prefix+"NormAnisotropy";
        saveas(f(j),fullfile(path,fnm+png));
        j=j+1;
        
        %plot3 An et At=f(t)
        if tit;title(ax(j),'Evolution of Tangencial Anisotropy');end
        ylabel(ax(j),'Anisotropy')
        xlabel(ax(j),xlab)
        fnm=prefix+"TanAnisotropy";
        saveas(f(j),fullfile(path,fnm+png));
        j=j+1;
        
        %plot4 TotCtcts=f(t)
        if tit;title(ax(j),'Variation of the total nb of contacts');end
        ylabel(ax(j),'Contacts')
        xlabel(ax(j),xlab)
        fnm=prefix+"nbContacts";
        saveas(f(j),fullfile(path,fnm+png));
        j=j+1;
        
        %plot5 VCtcts=f(t)
        if tit;title(ax(j),'Variation of the number of vertical contacts');end
        ylabel(ax(j),'Contacts percentage')
        xlabel(ax(j),xlab)
        fnm=prefix+"nbVContacts";
        saveas(f(j),fullfile(path,fnm+png));
        j=j+1;
        
        %plot6 HCtcts=f(t)
        if tit;title(ax(j),'Variation of the number of horizontal contacts');end
        ylabel(ax(j),'Contacts percentage')
        xlabel(ax(j),xlab)
        fnm=prefix+"nbHContacts";
        saveas(f(j),fullfile(path,fnm+png));
        
    otherwise
%         if pD.ConsoTime>pD.N1*pD.TimeStep && pD.ConsoTime<pD.N2*pD.TimeStep
%             tf=1;
%         else
%             tf=0;
%         end
        j=1;
        %Z=f(Ez)
        plotMark(app,ax(j),pD.Results(:,end),pD.Results(:,2));
        j=j+1;
        %TotCtcts=f(Ez)
        plotMark(app,ax(j),pD.Results(:,end),pD.Results(:,7));
        j=j+1;
        %Ver & Hor Ctcts=f(Ez)
        plotMark(app,ax(j),pD.Results(:,end),pD.Results(:,8)./pD.Results(:,7))
        plotMark(app,ax(j),pD.Results(:,end),pD.Results(:,9)./pD.Results(:,7));
        j=j+1;
        %A=f(Ez)
        plotMark(app,ax(j),pD.Results(:,end),pD.Results(:,3))
        plotMark(app,ax(j),pD.Results(:,end),pD.Results(:,4));
        j=j+1;
        %An & At=f(Ez)
        plotMark(app,ax(j),pD.Results(:,end),pD.Results(:,5))
        plotMark(app,ax(j),pD.Results(:,end),pD.Results(:,6));
        j=j+1;
        %An & At=f(Ez)
        plotMark(app,ax(j),pD.Results(:,end),sum(pD.Results(:,[3,5,6]),2))
        
%         if tf
%             for j=1:nb
%                 xline(ax(1),pD.ConsoStrain,'Color','k')
%             end
%         end
        
        if leg
            j=3;
            legend(ax(j),'Vertical','Horizontal','location','best');j=j+1;
            legend(ax(j),'Geometric','Mechanical','location','best');j=j+1;
            legend(ax(j),'Normal','Tangencial','location','best');
        end
        %create suffix for filenames
        suffix=pD.N1+"to"+pD.N2+"int"+pD.Interval;
        
        j=1;
        %saveplots
        %Z=f(Ez)
        if tit;title(ax(j),'Evolution of Coodination Number');end
        ylabel(ax(j),'Coordination Number')
        xlabel(ax(j),xlab)
        fnm=pD.Prefix(1)+"CoordNumber"+suffix;
        saveas(f(j),fullfile(path,fnm+png));
        j=j+1;
        %TotCtcts=f(Ez)
        if tit;title(ax(j),'Evolution of the Nb of contacts');end
        ylabel(ax(j),'Number of Contacts')
        xlabel(ax(j),xlab)
        fnm=pD.Prefix(1)+"TotCtcts"+suffix;
        saveas(f(j),fullfile(path,fnm+png));
        j=j+1;
        %Ver & Hor Ctcts=f(Ez)
        if tit;title(ax(j),'Evolution of the Nb of contacts');end
        ylabel(ax(j),'Ratio')
        xlabel(ax(j),xlab)
        fnm=pD.Prefix(1)+"HorVerCtcts"+suffix;
        saveas(f(j),fullfile(path,fnm+png));
        j=j+1;
        %A=f(Ez)
        if tit;title(ax(j),'Evolution of Anisotropy');end
        ylabel(ax(j),'Anisotropy')
        xlabel(ax(j),xlab)
        fnm=pD.Prefix(1)+"FabricAnisotropy"+suffix;
        saveas(f(j),fullfile(path,fnm+png));
        j=j+1;
        %An & At=f(Ez)
        if tit;title(ax(j),'Evolution of Normal & Tangecial Anisotropy');end
        ylabel(ax(j),'Anisotropy')
        xlabel(ax(j),xlab)
        fnm=pD.Prefix(1)+"NormTanAnisotropy"+suffix;
        saveas(f(j),fullfile(path,fnm+png));
        j=j+1;
        %0.5(Ag+An+At)=f(Ez)
        if tit;title(ax(j),'Evolution of 0.5(Ag+An+At)');end
        ylabel(ax(j),'Anisotropy')
        xlabel(ax(j),xlab)
        fnm=pD.Prefix(1)+"SumAnisotropy"+suffix;
        saveas(f(j),fullfile(path,fnm+png));
end

%outside legend
if ~leg && numel(pD)>1
    %Multi legend
    o=copyobj(f(1),0);
    l=legend(o.CurrentAxes,pD.FileName,'Orientation','horizontal');
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
    saveas(o,path+"Legend"+png);
    
    %Delete extra figures
    delete(o);
end
if isequal(app.LIGGGHTSAnalysisButtonGroup.SelectedObject,app.ExeAllButton)
    delete(f);
end
end
