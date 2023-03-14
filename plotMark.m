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

%Make sure axis is on hold
hold(ax,'on')

%Check what kind of kine to plot
if isequal(app.CourbePointsSwitch.Value,'On')
    lType="-+";
else
    lType="-";
end

try lw=app.PlotWidthEF.Value;
catch
    lw=1.5;
end

%Plot
chk=strcmp(varargin,'Color');
if sum(chk)==1
    curb=plot(ax,xData,yData,lType,'Color',varargin{find(chk)+1},...
        'LineWidth',lw);
else
    curb=plot(ax,xData,yData,lType,'LineWidth',lw);
end

%Add Initial/Final/Peak values
if isequal(app.MarkSwitch.Value,'On')
    opts={'MarkerEdgeColor',curb.Color,'LineWidth',lw};
    %Mark first point as square
    p=scatter(ax,curb.XData(1),curb.YData(1),'^',opts{:});
    set(get(get(p,'Annotation'),'LegendInformation'),...
        'IconDisplayStyle','off');
    %mark l;ast point as Triangle
    p=scatter(ax,curb.XData(end),curb.YData(end),'s',opts{:});
    set(get(get(p,'Annotation'),'LegendInformation'),...
        'IconDisplayStyle','off');
    %find max point, if it is different
    if any(strcmpi(varargin,'PEAK'))
        %peak may be negative or positive
        [~,mxp]=max(curb.YData);
        [~,mnp]=min(curb.YData);
        dmxp1=abs(curb.XData(mxp)-curb.XData(1));
        dmxp2=abs(curb.XData(mxp)-curb.XData(end));
        dmnp1=abs(curb.XData(mnp)-curb.XData(1));
        dmnp2=abs(curb.XData(mnp)-curb.XData(end));
        %if value is ~= from first and last point mark it
        if (mxp~=1 && mxp~=numel(curb.YData)) && ...
                (min(dmxp1,dmxp2)>min(dmnp1,dmnp2))
            p=scatter(ax,curb.XData(mxp),curb.YData(mxp),'o',opts{:});
            set(get(get(p,'Annotation'),'LegendInformation'),...
                'IconDisplayStyle','off');
        elseif mnp~=1 && mnp~=numel(curb.YData)
            p=scatter(ax,curb.XData(mnp),curb.YData(mnp),opts{:});
            set(get(get(p,'Annotation'),'LegendInformation'),...
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
                p=scatter(ax,varargin{f(i)+1},varargin{f(i)+2},mk(k),opts{:});
                set(get(get(p,'Annotation'),'LegendInformation'),...
                    'IconDisplayStyle','off');
                k=k+1;
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
                p=scatter(ax,xData(mnpos),yData(mnpos),mk(k),opts{:});
                set(get(get(p,'Annotation'),'LegendInformation'),...
                    'IconDisplayStyle','off');
                k=k+1;
            end
        end
    end
end
end

