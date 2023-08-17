function plotMark(app,ax,xData,yData,varargin)
%MARKPOINTS Mark points in the curb and axis given as argument
%   This function will be responsable to plot most of the graphs in the
%   app. It will check if app properties CourbePointsSwitch and MarkSwitch
%   are on and behave acordingly
%   If CourbePointsSwitch is 'on' it will mark each point of calculation in
%   the curb with a '+'
%   If MarkSwitch is 'on' it will mark with a triangle and square the
%   beginign and end of the curb. Further modifiers can be added :
%       - if 'Peak' is added as a varargin a circle will mark the highest Y
%       value.
%       - If 'Point' is added, followed by [x,y], a mark  will be added at
%       the point [x,y]. 
%       - If 'Pointx' is added, followed by x value, a mark will be added
%       at the correspinding y of the curb (closest one )
%       - For the latter only 3 marks have been added in the code
%       ["*",'o','+']
%   If 'Color' modified is added, curve color will be set as the next
%   argument
%   If 'Background' modifier is added, curve will be plot as grey in the
%   background and will be removed from the legend

%Make sure axis is on hold
hold(ax,'on')

%keyword checks
chkLT=strcmp(varargin,'LineType');
chkCl=strcmp(varargin,'Color');
chkBG=strcmp(varargin,'Background');
chkLeg=strcmp(varargin,'Nolegend');

%Choose linetype
if sum(chkLT)==1
    lType=varargin{find(chkLT)+1};
else
    if isequal(app.CourbePointsSwitch.Value,'On')
        lType="-+";
    else
        lType="-";
    end
end

%Get plot width
try lw=app.PlotWidthEF.Value;
catch
    lw=1.5;
end

%Plot
if sum(chkBG)==1
    if sum(chkCl)==1
        cv=plot(ax,xData,yData,'--','Color',varargin{find(chkCl)+1},...
            'LineWidth',max(0.5,lw-.5));
    else
        cv=plot(ax,xData,yData,'--','Color','#A1A1A1',...
            'LineWidth',max(0.5,lw-.5));
    end
    if sum(chkLeg)==1
        set(get(get(cv,'Annotation'),'LegendInformation'),...
            'IconDisplayStyle','off');
    end
    return
else
    if sum(chkCl)==1
        cv=plot(ax,xData,yData,lType,'Color',varargin{find(chkCl)+1},...
            'LineWidth',lw);
    else
        cv=plot(ax,xData,yData,lType,'LineWidth',lw);
    end
end
if sum(chkLeg)==1
    set(get(get(cv,'Annotation'),'LegendInformation'),...
                'IconDisplayStyle','off');
end

%Add Marks
if isequal(app.MarkSwitch.Value,'On')
    %If 'Color' was called, check if there is 4 arguments. If yes, the
    %latter refers to the alpha
    alp=1;
    if any(chkCl)
        C=varargin{find(chkCl)+1};
        if numel(C)>3;alp=C(4);end
    end
    opts={'MarkerEdgeColor',cv.Color,'LineWidth',lw,'MarkerEdgeAlpha',alp};
    %Mark first point as square
    sct=scatter(ax,cv.XData(1),cv.YData(1),'^',opts{:});
    set(get(get(sct,'Annotation'),'LegendInformation'),...
        'IconDisplayStyle','off');
    %mark l;ast point as Triangle
    sct=scatter(ax,cv.XData(end),cv.YData(end),'s',opts{:});
    set(get(get(sct,'Annotation'),'LegendInformation'),...
        'IconDisplayStyle','off');
    %find max point, if it is different
    if any(strcmpi(varargin,'PEAK'))
        %peak may be negative or positive
        [~,mxp]=max(cv.YData);
        [~,mnp]=min(cv.YData);
        dmxp1=abs(cv.XData(mxp)-cv.XData(1));
        dmxp2=abs(cv.XData(mxp)-cv.XData(end));
        dmnp1=abs(cv.XData(mnp)-cv.XData(1));
        dmnp2=abs(cv.XData(mnp)-cv.XData(end));
        %if value is ~= from first and last point mark it
        if (mxp~=1 && mxp~=numel(cv.YData)) && ...
                (min(dmxp1,dmxp2)>min(dmnp1,dmnp2))
            sct=scatter(ax,cv.XData(mxp),cv.YData(mxp),'o',opts{:});
            set(get(get(sct,'Annotation'),'LegendInformation'),...
                'IconDisplayStyle','off');
        elseif mnp~=1 && mnp~=numel(cv.YData)
            sct=scatter(ax,cv.XData(mnp),cv.YData(mnp),'o',opts{:});
            set(get(get(sct,'Annotation'),'LegendInformation'),...
                'IconDisplayStyle','off');
        end
    end
    mk=["diamond","o","pentagram","hexagram","x"];k=1;
    if any(strcmpi(varargin,'POINT'))
        %Point must be followe by two arguments, x and y of the chosen
        %point to be plot
        f=find(strcmpi(varargin,'POINT'));
        if numel(varargin)<max(f)+3
            fprintf("plotMark function was called with Point argument "+...
                " but without the required argumetns for the execution\n")
        else
            for i=1:numel(f)
                sct=scatter(ax,varargin{f(i)+1},varargin{f(i)+2},mk(k),opts{:});
                set(get(get(sct,'Annotation'),'LegendInformation'),...
                    'IconDisplayStyle','off');
                k=k+1;if k>numel(mk);k=1;end
            end
        end
    end
    if any(strcmpi(varargin,'POINTX'))
        %Pointx must be followed by one argument, the x position of the
        %point to the plot. The actual value to be added will be the 
        %closest point from xData found.
        
        f=find(strcmpi(varargin,'POINTx'));
        if numel(varargin)<max(f)+1
            fprintf("plotMark function was called with Point argument "+...
                " but without the required argumetns for the execution\n")
        else
            for i=1:numel(f)
                [~,mnpos]=min(abs(xData-varargin{f(i)+1})); %find the
                sct=scatter(ax,xData(mnpos),yData(mnpos),mk(k),opts{:});
                set(get(get(sct,'Annotation'),'LegendInformation'),...
                    'IconDisplayStyle','off');
                k=k+1;if k>numel(mk);k=1;end
            end
        end
    end
end
end

