function exe_BasicInfo(app)
%DDISTRIBUTION Distribution of diameters
%	This function does nothing particularly impressive. It will read the
%	data created by the class 'grains', calculate the eignenvalues of the
%	fabric tensor and plot the important values.

gr=grains('BASIC',app.N1EF.Value,'',app);
TD=app.TrialData;

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
path=MakePath(app,'BSCINF');
png=".png";

%Check title and legends option
if app.TitlesCB.Value;tit=1;else;tit=0;end
if app.LegendsCB.Value;leg=1;else;leg=0;end

%Prepare histogram plot
D=gr.Radius*2000; %calculate diamiter in mm
[uD,~,iuD]=unique(D);

nb=5;
if size(uD,1)>1;nb=nb+3;end
if app.Bool3D;nb=nb+1;end
f(nb)=figure;ax(nb)=axes(f(nb));hold(ax(nb),'on');
for i=1:(nb-1)
    f(i)=figure;ax(i)=axes(f(i));hold(ax(i),'on'); %#ok<LAXES>
end
i=1;
if size(uD,1)>1 %jump on monodisp case
    dDif=(max(uD)-min(uD))/((size(uD,1)-1)*2);
    cnt = accumarray(iuD,1); %count the nb of each Diameter

    %calculate total mass and divide by total mass
    mass=1/6*pi()*(uD.*uD.*uD)*2600*10^-6; %calculate mass of each D in grams
    cntM=cnt.*mass;
    cntM=cntM./sum(cntM);
    
    %Mass Histogram
    bar(ax(i),uD,cnt,dDif) 
    if tit;title(ax(i),'Diameter distribution');end
    ylabel(ax(i),'Number')
    xlabel(ax(i),'Diameter (mm)')
    fnm="HistogramDistributionNB";
    saveas(f(i),fullfile(path,fnm+png));
    i=i+1;

    %Mass Histogram
    bar(ax(i),uD,cntM,dDif) 
    if tit;title(ax(i),'Diameter distribution');end
    ylabel(ax(i),'Mass percentage')
    xlabel(ax(i),'Diameter (mm)')
    fnm="HistogramDistributionMASS";
    saveas(f(i),fullfile(path,fnm+png));
    i=i+1;
    %Prepare passing percentage
    for j=2:size(cntM,1)
        cntM(j)=cntM(j)+cntM(j-1);
    end
    %Mass Passing percentage
    plotMark(app,ax(i),[2*uD(1)-uD(2);uD],[0;cntM])
    if tit;title(ax(i),'Diameter distribution : Pass percentage');end
    ylabel(ax(i),'Mass percentage')
    xlabel(ax(i),'Diameter (mm)')
    fnm="PassTDistribution";
    saveas(f(i),fullfile(path,fnm+png));
    i=i+1;
end
if app.Bool3D;dim=3;else;dim=2;end
%Third plot - volume evolution 
V = getVolume(TD,app,TD.Step,'');
plotMark(app,ax(i),TD.Step*app.TimeStep,V)
xline(ax(i),app.ConsoStep*app.TimeStep,'Color','k');
if tit;title(ax(i),'Volume evolution');end
ylabel(ax(i),'Volume (m3)');
xlabel(ax(i),'Time(s)');
fnm="VolumeEvolution";
saveas(f(i),fullfile(path,fnm+png));
i=i+1;

%dEv comp part
l=find(TD.Step==app.ConsoStep);
V=V(l:end);
dV=(V(1:end-1)-V(2:end))./V(1:end-1);
dV=[0;-cumsum(dV,1)];
snExt = extStrains(TD,TD.Step(l:end),app.ConsoStep,app,'dev');
plot(ax(i),snExt(:,dim),dV,lType,'LineWidth',lw);
dV=(V(1:end-1)-V(2:end))./V(1);
dV=[0;-cumsum(dV,1)];
plot(ax(i),snExt(:,dim),dV,'o','LineWidth',lw)
%dV2ord=((snExt(:,1)+1).*(snExt(:,2)+1).*(snExt(:,3)+1)-1);
%plot(ax(i),snExt(:,dim),dV2ord,lType)
plot(ax(i),snExt(:,dim),-snExt(:,end-1),lType)
if tit;title(ax(i),'dVolume evolution');end
ylabel(ax(i),'dVolume');
xlabel(ax(i),'Time(s)');
if leg;legend(ax(i),'dV/V','dV/Vo','sum dExyz');end
fnm="dVolumeEvolution";
saveas(f(i),fullfile(path,fnm+png));
i=i+1;

%Third plot - Displacement piston Z
plotMark(app,ax(i),TD.Step*app.TimeStep,TD.Dz)
xline(ax(i),app.ConsoStep*app.TimeStep,'Color','k');
if tit;title(ax(i),'Displacement of the Piston Z');end
ylabel(ax(i),'Vertical displacement (m)');
xlabel(ax(i),'Time(s)');
fnm="PzDisplacement";
saveas(f(i),fullfile(path,fnm+png));
i=i+1;

if app.Bool3D
    %Fourth plot - Displacement lateral piston
    plot(ax(i),TD.Step*app.TimeStep,TD.Dy1,'LineWidth',lw)
    plot(ax(i),TD.Step*app.TimeStep,-TD.Dy2,'o','LineWidth',lw)
    plot(ax(i),TD.Step*app.TimeStep,TD.Dx1,'LineWidth',lw)
    plot(ax(i),TD.Step*app.TimeStep,-TD.Dx2,'^','LineWidth',lw)
    xline(ax(i),app.ConsoStep*app.TimeStep,'Color','k');
    if tit;title(ax(i),'Displacement of lateral Pistons');end
    ylabel(ax(i),'Lateral displacement (m)');
    xlabel(ax(i),'Time(s)');
    if leg
        legend(ax(i),'Piston Y1','(-1)*Piston Y2','Piston X1',...
            '(-1)*Piston X2','location','eastoutside');
    end
    fnm="LateralPistonDpl";
    saveas(f(i),fullfile(path,fnm+png));
    
    i=i+1;
    %Fifth plot - normalize the fourth plot to match take in the accoun
    %their diferent widths
    pX=app.boxXEF.Value;
    pY=app.boxYEF.Value;
    nX=pY/(pX+pY);
    nY=pX/(pX+pY);
    plot(ax(i),TD.Step*app.TimeStep,TD.Dy1*nY,'LineWidth',lw)
    plot(ax(i),TD.Step*app.TimeStep,-TD.Dy2*nY,'o','LineWidth',lw)
    plot(ax(i),TD.Step*app.TimeStep,TD.Dx1*nX,'LineWidth',lw)
    plot(ax(i),TD.Step*app.TimeStep,-TD.Dx2*nX,'^','LineWidth',lw)
    xline(ax(i),app.ConsoStep*app.TimeStep,'Color','k');
    if tit;title(ax(i),'Normalized Displacement of lateral Pistons');end
    ylabel(ax(i),'Lateral displacement (m2)');
    xlabel(ax(i),'Time(s)');
    if leg
        legend(ax(i),'Piston Y1','(-1)*Piston Y2','Piston X1',...
            '(-1)*Piston X2','location','eastoutside');
    end
    fnm="LateralPistonDplNorm";
    saveas(f(i),fullfile(path,fnm+png));
    i=i+1;
else
    %Fourth plot - Displacement lateral piston
    plot(ax(i),TD.Step*app.TimeStep,TD.Dy1,'LineWidth',lw)
    plot(ax(i),TD.Step*app.TimeStep,-TD.Dy2,'o','LineWidth',lw)
    xline(ax(i),app.ConsoStep*app.TimeStep,'Color','k');
    if tit;title(ax(i),'Displacement of the Piston Y');end
    ylabel(ax(i),'Lateral displacement (m)');
    xlabel(ax(i),'Time(s)');
    if leg
        legend(ax(i),'Piston Y1','(-1)*Piston Y2','location','eastoutside');
    end
    fnm="LateralPistonDisplacement";
    saveas(f(i),fullfile(path,fnm+png));
    i=i+1;
end

%Evolution of surfaces
l=TD.boxY-TD.Dy1+TD.Dy2;
h=TD.boxZ+TD.Dz;
if app.checkPiston
    w=TD.boxX-TD.Dx1+TD.Dx2;
    plot(ax(i),TD.Step*app.TimeStep,h.*l/(h(1)*l(1)),lType,'LineWidth',lw)
else
    w=TD.boxX;
end
plot(ax(i),TD.Step*app.TimeStep,h.*w/(h(1)*w(1)),'LineWidth',lw)
plot(ax(i),TD.Step*app.TimeStep,w.*l/(w(1)*l(1)),'o','LineWidth',lw)
if tit;title(ax(i),'Surface Size Evolution');end
ylabel(ax(i),'Ratio from intial value');
xlabel(ax(i),'Time(s)');
if leg
    if app.checkPiston
        legend(ax(i),'SurfaceX','SurfaceY','SurfaceZ','location','eastoutside');
    else
        legend(ax(i),'SurfaceY','SurfaceZ','location','eastoutside');
    end
end
fnm="LateralPistonDisplacement";
saveas(f(i),fullfile(path,fnm+png));


%Delete figures in the case of exe all
if isequal(app.LIGGGHTSAnalysisButtonGroup.SelectedObject,app.ExeAllButton)
    delete(f);
end
end
