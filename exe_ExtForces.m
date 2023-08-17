function exe_ExtForces(mode,app)
%EXTERNALFORCES  calculate the Stress-Strain graphs for the external forces
%   In this function we will only need the files that contain the efforts
%   and displacements of the pistons. We will use the class 'Trial data' to
%   read the files and then we will analyse the data.
%   We are interested in the graphs that belong to the area between the
%   begining of the compression and when the moment when the vertical
%   piston have compressed up to 25% of its original height. We are not
%   analysing the consolidation part.

%check for Load
if isequal(upper(mode),'LOAD');extLoad(app);return;end

%start calculation pannel
CalcPanel(app,'','','Starting calculation','on');

%Getting variables
TD=app.TrialData;
if isempty(TD.Dz);return;end

%Load values
N1=app.N1EF.Value;
N2=app.N2EF.Value;
interval=app.CalcInt.Value;
if interval==app.IntervalEF.Value
    f=find(TD.Step==N1 | TD.Step==N2);
    steps=TD.Step(f(1):f(2));
else
    steps=(N1:interval:N2)';
    if steps(end)~=N2; steps=[steps;N2];end
end
snExt = extStrains(TD,steps,N1,app,'dev');
ssExt = extStress(TD,steps,app,'separatePist');
consoStrain = extStrains(TD,app.ConsoStep,N1,app);
%calculate inertial term
gr=grains('BASIC',N1,'',app);
itv=steps(2:end)-steps(1:end-1);
devStrRat=(snExt(2:end,end)-[0;snExt(2:end-1,end)])./(itv*app.TimeStep);
mass=2600*(4/3*pi()*sum(gr.Radius.^3));
d=(TD.boxZ+TD.Dz(ismember(TD.Step,steps(2:end))));
I=(devStrRat.*sqrt(mass./(ssExt(2:end,end)*10^3.*d)));
%save values
Res.Inertia=I;
Res.Stress=ssExt;
Res.Strain=snExt;
Res.Step=steps;

fnm=fullfile(MakePath(app,'FORCEXT','check'),'Drained_M_Values.txt');
if isfile(fnm) && app.SimType~=1
    opts = detectImportOptions(fnm);
    tab = readtable(fnm,opts);
    tab=tab{:,:};
    Res.SlopeM=tab;
else
    Res.SlopeM='';
end
%Create plotData object

pD = plotData("Normal",Res,app,'Normal',consoStrain(end));
fnm=MakePath(app,'FORCEXT')+"ExternalForces"+N1+"to"+N2+"int"+interval+".mat";
save(fnm,'pD','-v7.3');

CalcPanel(app,'','','','off');
extPlot(app,pD);
%remove calculation pannel
end
function extLoad(app)
%Load file
%Ask for file
[fnm,fPath]=MatLoader('FORCEXT',app);
if fPath==0;return;end
%read prefix

if ~isa(fnm,'cell')
    pD=load([fPath fnm]).pD;
else
    nb=numel(fnm);
    pD(nb)=load(fullfile(fPath,fnm{end})).pD;
    for i=1:nb
        if i~=nb;pD(i)=load(fullfile(fPath,fnm{i})).pD;end
        fn=fnm{i};
        pD(i).FileName=fn(1:(find(ismember(fn,'.'),1,'last')-1));
    end
end
fnm=fullfile(fPath,'Drained_M_Values.txt');
if isfile(fnm) && (max([pD.SimType])~=1)
    opts = detectImportOptions(fnm);
    opts.Delimiter='|';
    opts.VariableNamesLine = 3;
    opts.DataLines=[4 Inf];
    tab = readtable(fnm,opts);
    tab=tab{:,1:2};
    %fuck matlab file reading funcitons
    if iscell(tab)
        if ischar(tab{1,1})
            tab=cellfun(@str2num, tab);
        end
    end
    pD(1).Results.SlopeM=tab;
end
pD(1).Prefix='Load';
%Plot
extPlot(app,pD);
end
function extPlot(app,pD)
%EXTPLOT Plot data related to extforcefiles

%Greek characters for laendbels
sigma=(char(963));
epsilon=(char(949));
%delta=(char(916));

%check all pD have the same amount of columns (9 3D, 6 2D)
nC=size(pD(1).Results.Stress,2);
for i=2:numel(pD)
    if size(pD(i).Results.Stress,2)~=nC
        warndlg("Loaded files do not have the same amount of columns,"+...
            "cannot compare 2D to 3D files.");
        return
    end
end
set(0,'defaultAxesFontSize',app.FontSizeEF.Value)

%variable to know the correct column of data
if nC==9
    clmnP=5;
else
    clmnP=4;
end

%Check title and legends option
if app.TitlesCB.Value;tit=1;else;tit=0;end
if app.LegendsCB.Value;leg=1;else;leg=0;end

%Path
if strcmpi(pD(1).Prefix,'LOAD')
    path=MakePath(app,'FORCEXTL');
else
    path=MakePath(app,'FORCEXT');
end
pathB=MakePath(app,'FORCEXT');
png=".png";

%Base number of plots
nb=15;
%If we are in 2D the unit of the stress should be kPa/m not kPa
if size(pD(1).Results.Stress,2)>6
    unit='[kPa]';
else
    unit='[kPa/m]';nb=nb-1;
end

%Create Figures and Axis
if numel(pD)>1;nb=nb-5;end
%Add interstitial pressure calculation if all are undrained
iPre=1;
for i=1:numel(pD);if pD(i).SimType~=2;iPre=0;end;end
nb=nb+iPre;
f(nb)=figure;ax(nb)=axes(f(nb));hold(ax(nb),'on');
for i=1:(nb-1)
    f(i)=figure;ax(i)=axes(f(i));hold(ax(i),'on'); %#ok<*LAXES> 
end

%check if Mfile was loaded. If it wasnt, do the necessary calculations.
if ~isempty(pD(1).Results.SlopeM)
    Mc=pD(1).Results.SlopeM(1);
    Mp=pD(1).Results.SlopeM(2);
    doM=0;
else
    doM=1;Mc=0;Mp=0;
end

%matlab base graph colors
if numel(pD)<8
    C = app.PlotColors;
else
    C = graphClrCode(size(pD,2));%plot colorcode
end
%create vectors to store value for Mslope calculations
Mvals=zeros(numel(pD),2,2); 
for j=1:numel(pD)
    vez={};vp={};
    for i=1:numel(pD(j).InfPts.q)
        vp=[vp,{'Pointx',pD(j).InfPts.p(i)}]; %#ok<AGROW>
        vez=[vez,{'Pointx',pD(j).InfPts.ez(i)}]; %#ok<AGROW>
    end
    %Stress : [sigX,sigY,sigZ,q,p]  Strain : [Ex,Ey,Ez,Ev,Ed]
    Res=pD(j).Results;i=1;
    %q=f(Ez)
    plotMark(app,ax(i),Res.Strain(:,end-2),Res.Stress(:,clmnP-1),...
        'Color',C(j,:),vez{:});i=i+1;
    %graph p=f(Ez)
    plotMark(app,ax(i),Res.Strain(:,end-2),Res.Stress(:,clmnP),...
        'Color',C(j,:),vez{:});i=i+1;
    %graph q=f(p) + Slope M
    plotMark(app,ax(i),Res.Stress(:,clmnP),Res.Stress(:,clmnP-1),...
        'Color',C(j,:),vp{:});
        %M_CSL - Simtype 1 and 3
    %add slope M that connects the origin to the steady state marked by
    %the last point on the graph
    Mvals(j,:,1)=[Res.Stress(end,clmnP-1),Res.Stress(end,clmnP)];%[q,p]
        %M_PEAK - All
    %add slope M that connects the origin to the failure point marked by
    %the max q/p
    [~,mxPos]=max(Res.Stress(:,clmnP-1)./Res.Stress(:,clmnP));
        
    Mvals(j,:,2)=[Res.Stress(mxPos,clmnP-1),Res.Stress(mxPos,clmnP)];%[q,p]
    if j==numel(pD)
        ly=ax(i).YLim(2);lx=ax(i).XLim(2);
        %CRITICAL STATE line
        if  (~doM) || (pD(j).SimType~=3 && doM)
            srtM=sortrows(Mvals(:,:,1));
            if Mc==0;Mc=mean(srtM(:,1)./srtM(:,2));end
            phi=asind(3*Mc/(6+Mc));
            if ly*(1/Mc)>lx
                p2=plot(ax(i),[0 lx],[0 lx*Mc],'--k');
            else
                p2=plot(ax(i),[0 ly*(1/Mc)],[0 ly],'--k');
            end
            text(ax(i),max(srtM(:,2))*2/3,max(srtM(:,2))*Mc*2/3,...
                {'',[' \leftarrow M_{CSL}=' num2str(Mc,'%.3f')],...
                ['      \Phi_{CSL}=' num2str(phi,'%.2f') '\circ']},...
                'HorizontalAlignment','left','Color','black')
            set(get(get(p2,'Annotation'),'LegendInformation'),...
                'IconDisplayStyle','off');
        end
        %PEAK line
        srtM=sortrows(Mvals(:,:,2));
        if Mp==0;Mp=mean(srtM(:,1)./srtM(:,2));end
        phi=asind(3*Mp/(6+Mp));
        if ly*(1/Mp)>lx
            p2=plot(ax(i),[0 lx],[0 lx*Mp],'-.k');
        else
            p2=plot(ax(i),[0 ly*(1/Mp)],[0 ly],'-.k');
        end
        qp=(max(srtM(:,2))*2/3);
        if (qp/ax(i).XLim(2))<0.3
            qp=max(srtM(:,2));
        end
        text(ax(i),qp,qp*Mp,...
            {'',['M_{FL}=' num2str(Mp,'%.3f') ' \rightarrow  '],...
            ['\Phi_{FL}=' num2str(phi,'%.2f') '\circ    ']},...
            'HorizontalAlignment','right','Color','black')
        set(get(get(p2,'Annotation'),'LegendInformation'),...
            'IconDisplayStyle','off');
    end
    i=i+1;
    
    %graph q/p=f(Ez)
    plotMark(app,ax(i),Res.Strain(:,end-2),...
        Res.Stress(:,clmnP-1)./Res.Stress(:,clmnP),'Color',C(j,:),vez{:});i=i+1;
    %graph Ev=f(Ez)
    plotMark(app,ax(i),Res.Strain(:,end-2),Res.Strain(:,end-1),...
        'Color',C(j,:),vez{:});i=i+1;
    %graph Ev=f(p)
    plotMark(app,ax(i),Res.Stress(:,clmnP),Res.Strain(:,end-1),...
        'Color',C(j,:),vp{:});i=i+1;
    %graph Ev=f(Ed)
    plotMark(app,ax(i),Res.Strain(:,end),Res.Strain(:,end-1),...
        'Color',C(j,:),vp{:});i=i+1;
    %graph Ed=f(Ez)
    plotMark(app,ax(i),Res.Strain(:,end-2),Res.Strain(:,end),...
        'Color',C(j,:),vez{:});i=i+1;

    %Second order work calculation
    stn=Res.Strain(2:end,1:end-2)-Res.Strain(1:end-1,1:end-2);
    sts=Res.Stress(2:end,1:clmnP-2)-Res.Stress(1:end-1,1:clmnP-2);
    dW2=[0;sum(stn.*sts,2)];
    %normalized W2 if wanted
    %dW2=[0;sum(stn.*sts,2)./(sqrt(sum(sts.^2,2)).*sqrt(sum(stn.^2,2)))];
    plotMark(app,ax(i),Res.Strain(:,end-2),dW2,...
        'Color',C(j,:),vez{:});i=i+1;
    %stop here for multiload
    if numel(pD)==1
        %graph sigZ =f(Ez)
        plotMark(app,ax(i),Res.Strain(:,end-2),Res.Stress(:,clmnP-2),...
            'Color',C(j,:),vez{:});i=i+1;
        %graph sigY1 & sigY2 =f(EZ)
        plot(ax(i),Res.Strain(:,end-2),Res.Stress(:,end),'LineWidth',1.5);
        plot(ax(i),Res.Strain(:,end-2),Res.Stress(:,end-1),'+','LineWidth',1.5);
        i=i+1;
        if size(Res.Stress,2)>6
            %graph sigX1 & sigX2 =f(EZ)
            plot(ax(i),Res.Strain(:,end-2),Res.Stress(:,end-2),'LineWidth',1.5);
            plot(ax(i),Res.Strain(:,end-2),Res.Stress(:,end-3),'+','LineWidth',1.5);
            i=i+1;
            %graph(Ez and -Ex-Ey)
            plot(ax(i),Res.Strain(:,end-2),Res.Strain(:,end-2),'LineWidth',1.5);
            plot(ax(i),Res.Strain(:,end-2),-sum(Res.Strain(:,end-4:end-3),2),'LineWidth',1.5);
            i=i+1;
            %plot sigX sigY and sigZ in a graph
            plot(ax(i),Res.Strain(:,end-2),Res.Stress(:,1),...
                Res.Strain(:,end-2),Res.Stress(:,2),...
                Res.Strain(:,end-2),Res.Stress(:,3),'LineWidth',1.5);
        else
            %plot sigX sigY and sigZ in a graph
            plot(ax(end),Res.Strain(:,end-2),Res.Stress(:,1),...
                Res.Strain(:,end-2),Res.Stress(:,2),'LineWidth',1.5);
        end
        
        i=i+1;
    end
    
    %interts pressure
    if iPre
        row=(Res.Step==round(pD(j).ConsoTime/pD(j).TimeStep));
        consP=Res.Stress(row,end);consQ= Res.Stress(row,end-1);
        newQ=consQ+3*(Res.Stress(:,end)-consP);
        plotMark(app,ax(i),Res.Strain(:,end-2),newQ,'Color',C(j,:))
        i=i+1;
    end
    fi=find(Res.Inertia==0);
    for k=1:numel(fi)
        Res.Inertia(fi(k))=Res.Inertia(min(fi(k)-1,1));
    end

    %graph Inertia
    plotMark(app,ax(i),Res.Strain(2:end,end-2),log10(Res.Inertia),...
        'Color',C(j,:),vez{:});
end

if max([pD.SimType])==1 && doM
    %save M values if drained to be used in other cases
    fnm=fullfile(pathB,"Drained_M_Values.txt");
    fid = fopen(fnm, 'w');
    fprintf(fid, '## M Drained Values ##\n');
    fprintf(fid, ['Used on other types of simulations to plot peak '...
        ' and crit state lines\n']);
    fprintf(fid,'Mc|Mp\n');
    fprintf(fid, '%d|%d\n',Mc,Mp);
    fclose(fid);
end

if leg && numel(pD)>1
    switch pD(j).SimType
        case 1
            legloc=["northeast","northeast","southeast","southeast","southeast",...
                "northeast","northeast","southeast","northeast","southeast"];
        case 2
            legloc=["southeast","southeast","southeast","southeast","northwest",...
                "northwest","northwest","southeast","southeast","southeast","southeast"];
        case 3
            legloc=["southeast","northeast","southeast","southeast","southeast",...
                "northeast","northeast","southeast","northeast","southeast"];
    end
    
    for i=1:nb
        legend(ax(i),pD.FileName,'location',legloc(i))
    end
end
%Titles,labels and save graphs
i=1;
%q=f(Ez)
if tit; title(ax(i),['q=f(' epsilon 'z)']);end
ylabel(ax(i),['Deviatoric stress (q) ' unit])
xlabel(ax(i),'Axial Strain');
filename="Stress_q";
saveas(f(i),fullfile(path,filename+png));
i=i+1;

%p=f(Ez)
if tit; title(ax(i),['p=f(' epsilon 'z)']);end
ylabel(ax(i),['Mean stress (p) ' unit])
xlabel(ax(i),'Axial Strain');
ax(i).YLim(1)=0;
filename="Stress_p";
saveas(f(i),fullfile(path,filename+png));
i=i+1;


%q=f(p) + Pente M rupture
if tit; title(ax(i),'q=f(p)');end
ylabel(ax(i),['Deviatoric stress (q) ' unit])
xlabel(ax(i),['Mean stress (p) ' unit])
ax(i).XLim(1)=0;
filename="Stress_q=f(p)";
saveas(f(i),fullfile(path,filename+png));
i=i+1;

%q/p=f(Ez) 
if tit; title(ax(i),'q/p=f(Ez)');end
ylabel(ax(i),'Relative strength (q/p)')
xlabel(ax(i),'Axial Strain')
filename="Stress_qoverp";
saveas(f(i),fullfile(path,filename+png));
i=i+1;

%Ev=f(Ez)
yl=yline(ax(i),0,':','Color','k','LineWidth',1);
set(get(get(yl,'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','off');
if tit; title(ax(i),[epsilon 'v=f(' epsilon 'z)']);end
ylabel(ax(i),'Volumetric Strain')
xlabel(ax(i),'Axial Strain');
ax(i).YDir='reverse';
filename="Strain_Vol_Ez";
saveas(f(i),fullfile(path,filename+png));
i=i+1;

%Ev=f(p)
yl=yline(ax(i),0,':','Color','k','LineWidth',1);
set(get(get(yl,'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','off');
if tit; title(ax(i),[epsilon 'v=f(p)' ]);end
ylabel(ax(i),'Volumetric Strain')
xlabel(ax(i),['Mean stress (p) ' unit]);
ax(i).YDir='reverse';
ax(i).XLim(1)=0;
filename="Strain_Volum_p";
saveas(f(i),fullfile(path,filename+png));
i=i+1;

%Ev=f(Ed)
yl=yline(ax(i),0,':','Color','k','LineWidth',1);
set(get(get(yl,'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','off');
if tit; title(ax(i),[epsilon 'v=f(' epsilon 'd)' ]);end
ylabel(ax(i),'Volumetric Strain')
xlabel(ax(i),'Deviatoric Strain');
ax(i).YDir='reverse';
ax(i).XLim(1)=0;
filename="Strain_Volum_Ed";
saveas(f(i),fullfile(path,filename+png));
i=i+1;

%Ed=f(Ez)
if tit; title(ax(i),[epsilon 'd=f(' epsilon 'z)']);end
ylabel(ax(i),'Deviatoric Strain')
xlabel(ax(i),'Axial Strain');
filename="Strain_Dev";
saveas(f(i),fullfile(path,filename+png));
i=i+1;

%Ed=f(Ez)
yl=yline(ax(i),0,':','Color','k','LineWidth',1);
set(get(get(yl,'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','off');
if tit; title(ax(i),'Second order work');end
ylabel(ax(i),'Second order work')
xlabel(ax(i),'Axial Strain');
filename="W2";
saveas(f(i),fullfile(path,filename+png));
i=i+1;

if numel(pD)==1
    %sigZ=f(Ez)
    if tit; title(ax(i),[sigma 'Z=f(' epsilon '1)']);end
    ylabel(ax(i),['Stress Z ' unit])
    xlabel(ax(i),'Axial Strain');
    filename="StressZ=f(EZ)";
    saveas(f(i),fullfile(path,filename+png));
    i=i+1;

    %sigY1 & sigY2 =f(Ez)
    if tit; title(ax(i),[sigma 'Y=f(' epsilon 'z)']);end
    ylabel(ax(i),[sigma 'Y ' unit]);
    xlabel(ax(i),'Axial Strain');
    legend(ax(i),'Direction Y+','Direction Y-');
    filename="StressY=f(E1)";
    saveas(f(i),fullfile(path,filename+png));
    i=i+1;

    %sigX1 & sigX2 =f(Ez);
    if size(pD(1).Results.Stress,2)>6
        if tit; title(ax(i),[sigma 'X=f(' epsilon 'z)']);end
        ylabel(ax(i),[sigma 'X ' unit]);
        xlabel(ax(i),'Axial Strain');
        legend(ax(i),'Direction X+','Direction X-');
        filename="StressX=f(EZ)";
        saveas(f(i),fullfile(path,filename+png));
        i=i+1;
    end

    %Ez-(Ex+Ey)=f(Ez)
    if tit; title(ax(i),[epsilon 'Z vs -' epsilon '(X and Y)']);end
    ylabel(ax(i),epsilon);
    xlabel(ax(i),'Axial Strain');
    legend(ax(i),'Z','-(Y+X)');
    filename="DeformationCompare";
    saveas(f(i),fullfile(path,filename+png));
    i=i+1;

    %sigX,sigY,sizZ=f(Ez)
    if tit; title(ax(i),'Principal stress comparison');end
    if size(pD(1).Results.Stress,2)>6
        legend(ax(i),"Sigma X","Sigma Y","Sigma Z")
    else
        legend(ax(i),"Sigma Y","Sigma Z")
    end
    ylabel(ax(i),'Stress [MPa]')
    xlabel(ax(i),'Axial Strain');
    filename="Stress_Comp";
    saveas(f(i),fullfile(path,filename+png));
    i=i+1;

end



%u'= 3/1qtop -sig   - intert pressure
if iPre
    if tit; title(ax(i),[sigma char("'") ' =f(Ez)']);end
    ylabel(ax(i),'Interstitial pressure [kPa]');
    xlabel(ax(i),'Axial Strain');
    filename="StressInterst";
    saveas(f(i),fullfile(path,filename+png));
    i=i+1;
end

%Inertia
yl=yline(ax(i),-3,':','Color','k','LineWidth',1);
set(get(get(yl,'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','off');
yl=yline(ax(i),-2,':','Color','k','LineWidth',1);
set(get(get(yl,'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','off');
if tit; title(ax(i),'Spercimen inertia evolution');end
%ax(i).YLim=[-7,-1]
ylabel(ax(i),'log(Inertial Number)')
xlabel(ax(i),'Axial Strain');
filename="Inertia";
saveas(f(i),fullfile(path,filename+png));


%Qcst comparaision plot
if pD(j).SimType==3 && numel(pD)==1
    f(end+1)=figure;
    %title(f(end),"Qcst Stresses Comparasion");
    tf=tiledlayout(f(end),2,3,'TileSpacing','tight','Padding','Compact');
    
    %q=f(Ez)
    axf=nexttile(tf);
    plotMark(app,axf,Res.Strain(:,end-2),Res.Stress(:,clmnP-1),vez{:});
    ylabel(axf,'Deviatoric Stress');
    xlabel(axf,'Axial Strain');
    %graph q=f(p)
    axf=nexttile(tf);
    plotMark(app,axf,Res.Stress(:,clmnP),Res.Stress(:,clmnP-1),vp{:});
    ylabel(axf,'Deviatoric Stress');
    xlabel(axf,'Mean Stress');
    %jump a tile
    t=nexttile(tf);
    set(t,'visible','off')
    %graph Ev=f(Ez)
    axf=nexttile(tf);
    plotMark(app,axf,Res.Strain(:,end-2),Res.Strain(:,end-1),vez{:});
    axf.YDir='reverse';
    ylabel(axf,'Volumetric Strain');
    xlabel(axf,'Axial Strain');
    %graph Ev=f(p)
    axf=nexttile(tf);
    plotMark(app,axf,Res.Stress(:,clmnP),Res.Strain(:,end-1),vp{:});
    axf.YDir='reverse';
    ylabel(axf,'Volumetric Strain');
    xlabel(axf,'Mean Stress');
    %graph Ev=f(sig3)
    axf=nexttile(tf);
    plotMark(app,axf,Res.Stress(:,1),Res.Strain(:,end-1));
    axf.YDir='reverse';
    ylabel(axf,'Volumetric Strain');
    xlabel(axf,'Stress X');
    
    f(end).Position(3:4)=[3*f(end).Position(3),2*f(end).Position(4)];
    %Save file
    saveas(f(end),fullfile(path,"Qcst-Stresses"+png));
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
    set(l,'Position',[0,0,.9999,.9999]);
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
