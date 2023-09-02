function exe_GetHSG(app)
%FORCECHAINS Will read LIGGGHTS files and export data on the Loops
%   This function will use the class 'grains' to read files and then
%   calculate the per grain internal stress tensor. We will then select
%   those that have values above the mean and similar directions to
%   determinate the force chains. Based on
%   Peters,Muthswamy,Wibowo,Tordesillas 2005

%Load values
[N1,N2,interval,~,nbFiles] = createStepArray(app);

%Turn on calculation pannel
app=CalcPanel(app,'',nbFiles,'Starting calculation','on');

%Start the Loop
for i=1:nbFiles
    if getappdata(app.CalculatingPanel,'Canceling')
        CalcPanel(app,'','','','off');
        warndlg('Calculation was canceled');return
    end
    step=min(N1+interval*(i-1),N2);
    app=CalcPanel(app,i,nbFiles,step);
    gr=grains('STENSOR',step,'',app); %calculate tensor
    if isempty(gr.PGStressTensor)
        CalcPanel(app,'','','','off');
        warndlg("Stress tensor is empty on step " +step);return
    end
    pgST=gr.PGStressTensor./permute((4*pi()/3)*gr.Radius.^3,[3,2,1]);
    if i==1
        %initialize variables
        nG=size(pgST,3); %nb of grains
        hsgID=zeros(nG,1);
    end
    pStress=ones(nG,1);
    for j=1:nG
        eVal=eig(-pgST(:,:,j));  %get eig vals
        pStress(j)=max(eVal);
    end
    hsgID=hsgID+(pStress(:,1)>=mean(pStress(:,1)));
end
CalcPanel(app,i+1,nbFiles,'','off');
set(0,'defaultAxesFontSize',app.FontSizeEF.Value)

%Prepare Strain data
f=figure;ax=axes(f);hold(ax,'on');
scatter(ax,1:nG,hsgID/nbFiles,'x');
title(ax(1),'Ratio of time each grain is HS')
ylabel(ax(1),'Ratio')
xlabel(ax(1),'Grain ID')
fnm="FCNumberTrees";
saveas(f(1),path+fnm+png);
saveas(f(1),path+fnm+fig);
end

