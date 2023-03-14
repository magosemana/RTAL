function exe_EdgesRatio(app,mode)
%EDGESRATIO Summary of this function goes here
%   Detailed explanation goes here
%Turn on calculation pannel

if nargin>1
    if isequal(upper(mode),'LOAD');loadEdges(app);return;end %check Load
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
    stepArray=(N1:interval:N2)';
    if stepArray(end)~=N2; stepArray=[stepArray;N2];end
end
nbFiles=numel(stepArray);
%check dimesion
if app.Bool3D; D=3;else; D=2;end

res=zeros(nbFiles,3);
app=CalcPanel(app,'',nbFiles,'Starting calculation','on');
for i=1:nbFiles
    if getappdata(app.CalculatingPanel,'Canceling')
        CalcPanel(app,'','','','off');
        warndlg('Calculation was canceled');return
    end
    step=stepArray(i);
    app=CalcPanel(app,i,nbFiles,step);
    gr=grains('BASIC',step,app.PartialData,app); %calculate loops
    if isempty(app.PartialData)
        PD='';
    else
        PD = grainListing(app.PartialData,gr);
    end
    sc = spaceCellSystem("Edges",step,gr,app,PD);
    if isempty(sc.OpenEdges)
        CalcPanel(app,'','','','off');
        warndlg("ERROR:No closed/open edges found :"+step);return;
    end
    res(i,:)=[step size(sc.ClosedEdges,1) size(sc.OpenEdges,1) ];
end
CalcPanel(app,i+1,nbFiles,'','off');

%Calculate strain values and add to results
if app.SimType==3
    stress = extStress(app.TrialData,res(:,1),app);
    xAxis=stress(:,end);
else
    strain = extStrains(app.TrialData,res(:,1,1),N1,app);
    xAxis=strain(:,end); 
end 
consoStrain = extStrains(app.TrialData,app.ConsoStep,N1,app);
res=[xAxis res(:,2:3)];
%Save calculated values
%Inflection Points
infP=[0,0,0,0];
%Start file to store values
fnm=MakePath(app,'EDG')+"Edges"+N1+"to"+N2+"int"+interval+".txt";
fid = fopen(fnm, 'w');
fprintf(fid, '## Evolution of nb of Edges ##\n');
fprintf(fid, 'Simulation|Important|Values\n');
fprintf(fid, '%d|%d|%d\n',N1,N2,interval);
fprintf(fid, '%d|%d|%d\n',app.Bool3D,app.checkPiston,app.TimeStep);
fprintf(fid, '%d|%d|%d\n',consoStrain(end),app.ConsoStep*app.TimeStep,infP(1));
fprintf(fid, '%d|%d|%d\n',infP(2:4));
fprintf(fid, '-1|-1|-1\n');
fprintf(fid,'From step %d to %d with interval %d - \n',N1,N2,interval);
fprintf(fid, 'Strain|ClosedEdges|OpenEdges \n');
fprintf(fid,'%d|%d|%d \n',res');
fclose(fid);

%go to plot
plotEdges(app,plotData("Normal",res,app,'',consoStrain(D)))
end
function loadEdges(app)
pD=FileLoader(app,'EDG');
if isempty(pD);return;end
plotEdges(app,pD);
end
function plotEdges(app,pD)
%PLOT PART
%Check title and legends option
if app.TitlesCB.Value;tit=1;else;tit=0;end
if app.LegendsCB.Value;leg=1;else;leg=0;end

path=MakePath(app,'EDG');
png=".png";
C=app.PlotColors;

if numel(pD)==1
    nb=3;
    l=["Closed Edes","Open Edges"];
    mp=0;
else
    nb=5;
    l=[pD.FileName];
    mp=1;
end
f(nb)=figure;ax(nb)=axes(f(nb));hold(ax(nb),'on');
for i=1:(nb-1)
    f(i)=figure;ax(i)=axes(f(i));hold(ax(i),'on'); %#ok<LAXES>
end

%plot pcentage 
for j=1:numel(pD)
    tot=sum(pD(j).Results(:,2:3),2);i=1;c=j;
    plotMark(app,ax(i),pD(j).Results(:,1),pD(j).Results(:,2),'Color',C(c,:)) %Clsd
    plotMark(app,ax(i+1),pD(j).Results(:,1),pD(j).Results(:,2)./tot,'Color',C(c,:)) %Clsd ratio
    plotMark(app,ax(i+2),pD(j).Results(:,1),tot/tot(1),'Color',C(c,:)) %Total

    if numel(pD)==1;c=2;end
    plotMark(app,ax(i+3*mp),pD(j).Results(:,1),pD(j).Results(:,3),'Color',C(c,:)) %Open
    plotMark(app,ax(i+1+3*mp),pD(j).Results(:,1),pD(j).Results(:,3)./tot,'Color',C(c,:)) %Open ratio
    
    
end
if leg
    for i=1:(nb-1+mp)
        legend(ax(i),l,"Location","Southoutside","Orientation","Horizontal")
    end
end

i=1;
if tit;title(ax(i),'Evolution of number of Edges');end
ylabel(ax(i),'Number of Edges')
xlabel(ax(i),'Axial Strain')
saveas(f(i),fullfile(path,"Edges_Nb"+png));
i=i+1;
if tit;title(ax(i),'Evolution of ratio of Edges');end
ylabel(ax(i),'Ratio of Edges')
xlabel(ax(i),'Axial Strain')
saveas(f(i),fullfile(path,"Edges_Ratio"+png));
i=i+1;
if tit;title(ax(i),'Ratio of total edges compared to beggining');end
ylabel(ax(i),'Ratio of Edges')
xlabel(ax(i),'Axial Strain')
saveas(f(i),fullfile(path,"Edges_Total"+png));
i=i+1;
if nb==5
    if tit;title(ax(i),'Evolution of number of Edges');end
    ylabel(ax(i),'Number of Edges')
    xlabel(ax(i),'Axial Strain')
    saveas(f(i),fullfile(path,"Edges_Nb_Open"+png));
    i=i+1;
    if tit;title(ax(i),'Evolution of ratio of Edges');end
    ylabel(ax(i),'Ratio of Edges')
    xlabel(ax(i),'Axial Strain')
    saveas(f(i),fullfile(path,"Edges_Rat_Open"+png));
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

