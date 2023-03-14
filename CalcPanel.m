function app=CalcPanel(app,i,nbStep,stepstr,varargin)
%CALCULATINGPANEL Control the and update the calculating pannel
if isempty(app.TimeTracker)
    tic
    app.TimeTracker=zeros(2,1);
    h=0;m=0;s=0;
else
    t=toc;tic;
    app.TimeTracker(1)=(app.TimeTracker(1)*app.TimeTracker(2)+t)/...
        (app.TimeTracker(2)+1);
    app.TimeTracker(2)=app.TimeTracker(2)+1;
    remT=seconds(ceil((nbStep+1-i)*app.TimeTracker(1)));
    remT.Format='hh:mm:ss';
    [h,m,s]=hms(remT);
end

if isa(app.CalculatingPanel,'double') 
    %NOAPPEXECUTION
    if ~isempty(i)
        pct=floor((i-1)/nbStep*100);
        if i==1 || strcmpi(stepstr,'NEWLINE')
            fprintf('Executing - %02d%% - %02d:%02d:%02d\n',pct,h,m,s);
        else
            fprintf([repmat('\b', 1, 15) '%02d%% - %02d:%02d:%02d\n'],pct,h,m,s);
        end
    end
    return
else
    if strcmpi(varargin,'on')
        setappdata(app.CalculatingPanel,'Canceling',0);
        app.CalculatingPanel.set('Visible','on');
    end
    if strcmpi(varargin,'off')
        app.CalculatingPanel.set('Visible','off');
        app.TimeTracker='';
        return
    end
    if isempty(i)
        app.CalculatingPanel.set('Title',"Calculating");
        app.TimeRem.set('Text','Estimating remaining time...');
    else
        pct=floor((i-1)/nbStep*100);
        app.CalculatingPanel.set('Title',"Calculating - "...
            +pct+"%");
        app.TimeRem.set('Text',sprintf('Time Remaining %02d:%02d:%02d',h,m,s));
    end
    app.StepTracker.set('Text',"Step "+stepstr);
    drawnow
end
end

