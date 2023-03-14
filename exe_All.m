function exe_All(app,varargin)
%EXE_ALL executes all the selected calculations
%this function will execute all the selected calculations in one go. This
%way the user does not need to do one calculation at a time.

%Partial Data
switch app.PartialSwitch.Value
    case 'Total'
        PD='';
    case 'Partial'
        PD=app.PartialData;
end

%Start calculations
%Anisotropy
if app.Ani_EA_CB.Value
    fprintf('Executing Anisotropy\n');
    exe_Anisotropy(PD,app);
end
if getappdata(app.CalculatingPanel,'Canceling');return;end
%Loops
if app.Loop_EA_CB.Value
    fprintf('Executing Loops\n');
    exe_LoopsCounter(app,PD,'on');
end
if getappdata(app.CalculatingPanel,'Canceling');return;end

%Void Ratio
if app.VR_EA_CB.Value
    fprintf('Executing Void Ratio\n');
    exe_VoidRatio(PD,app);
end
if getappdata(app.CalculatingPanel,'Canceling');return;end

%Internal Stress
if app.SInt_EA_CB.Value
    fprintf('Executing Internal Stress\n');
    exe_IntForces(PD,app);
end
if getappdata(app.CalculatingPanel,'Canceling');return;end

%External Stress
if app.SOut_EA_CB.Value
    fprintf('Executing External Stress\n');
    exe_ExtForces(string(app.ExtFSwitch.Value),app);
end
if getappdata(app.CalculatingPanel,'Canceling');return;end

%Strain Tensor
if app.StrT_EA_CB.Value
    fprintf('Executing Internal Strain\n');
    exe_StrainTensor(app,app.StrT_EA_DropDown.Value,PD);
end
if getappdata(app.CalculatingPanel,'Canceling');return;end

%Basic information
if app.BI_EA_CB.Value
    fprintf('Executing Basic Information\n');
    exe_BasicInfo(app)
end

end

