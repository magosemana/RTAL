function exe_NoApp(path)
% NOAPPEXE will execute the calculations without the use of the app, this
% way it can be donne one after another without pause.

%% Variable definitions
%Start Step & end step - values must be well chosen to match the values of
%the simulation.
%   If N1==0 & N2==0 the calculation will be done in the entire simulation
%   If N1=='cs', calculation will stop with the end of consolidation
%   If N1=='cp', calculation will start after the end of consolidation
%   If N1=='dr' and Qcst simulation, run only after drained triaxial part.
%   works as 'cp' for other cases
%   If interval==0, the minimum interval will be chosen (defined by the
%   data of the simulation).
%   If interval== 1-9, the minimum interval will be multiplied by this
%   value.
N1='cp';
N2=0;
interval=0;

%Basic Calculations
    %Value of 0 means will NOT be executed while a value of 1.
Anisotropy=0;
BasicInfo=0;
ExtForces=1;
IntForces=0;
VoidRatio=1;
EdgeRatio=0;

%Calculations with multiple options
ForceChain={'BASE','CL'};
    %ForceChain accepts 0 (no calculation) or a cell array. Array values
    %must be chosen from the following (not necessaraly on order) :
    %'BASE' for the base fc study,
    %'CL' for cluster and force chain,
    %'BNDANGLE' for fc bending limit angle,
    %'CLTF' for fc bending+cluster transformation
    %'GRN' for following a certaing grain, followed by chosen gr ids in a
    %vector form [id1,id2..idn].
    %'CNST' to change parameters, followed by a 2value vector in the form
    %[Angle, multFactor].If no values determined default [45,1] will
    %be used.
    % EX: {'BASE','GRN',[1,383,1547],'CLTF','CL','CNST',[30,1.1]}
StrainTensor=3;
    %1 for Global, 2 PerCell, 3 Both.  W2Cluster can be done adding a 4 ,
    %however PerCell need to be executed in the same run as [2,4] or [3,4].
Loops={'VTK','ANI','DEF','VR'};         
    %Loops accepts 0 (no calculation) or a cell array. Array values
    %must be chosen from the following (not necessaraly on order) :
    %'VTK' for vtk file save,
    %'DEF' for cluster deformability,
    %'TRF' for cluster transformation,
    %'ANI' for cluster anisotropy,
    %'STN' for cluster strain,
    %'STS' for cluster stress.
    %'VR'  for cluster Voidratio.
    % EX: {'VTK','DEF','TRF','ANI','STN','STS','VR'}

%Partial Calculation
%   Value of 1 means it will be executed 
%   Value of 0 means will NOT be executed
%   If partial is to be executed, a rectangle need to be chosen. To LOAD a
%   rectangle it's name have to be given and the rectangle need to be
%   present on the MatlabResultsFile folder. To choose the values rectname
%   has to be 'Values' and the rectInfo needs to be completed. 
Partial=0;
rectname='Value';	%'filename.txt' to load or 'Value' to use numeric values   
rectInfo=[0.01,...	%originY
        0.055,...   %originZ
        0.35,...    %length (y)
        0.2,...     %height (z)
        0,...       %width  (x)
        55];        %rotation angle

%Subdivision
%   Value of 1 means the partial calculation will be subdvided
%   Value of 0 means the partial calculation will NOT be subdvided
%   Lines and Columns controls the number of subdvided squares that will be
%   in the calculation.
Sub=1;          %#ok<*NASGU>
Lines=10;        
Columns=1;      

%Plot options
fSz=[560,420];  %figure size
tit=0;          %Plot titles. 1 for yes, 0 for no
leg=0;          %Plot legends. 1 for yes, 0 for no
mrk=1;
pts=0;

%% Preparation --- modify with caution ---
%In this part an object will be created with all the necessary information
%to simulate the app and make all funciton work properly. Modifications in
%this part may break this function


%Start the 'app' object that will substitute the application and contain
%the values.
if ~isa(path,'char');path=convertStringsToChars(path);end
app.CalculatingPanel=0;
app.FormatEF.Value='.txt';
app.GrainsEF.Value='grains';
app.PistEF.Value='servoForce';
app.GrainsForceEF.Value='contForce';
app.Figsize=fSz;
app.SimType='';
app.TitlesCB.Value=tit;
app.LegendsCB.Value=leg;
app.TimeTracker='';
if mrk;of='On';else;of='Off';end
app.MarkSwitch.Value=of;
if pts;of='On';else;of='Off';end
app.CourbePointsSwitch.Value=of;
if nargin>0
    if path(end) ~= '/'; path = [path '/'];end
else
    path='..';
end
[app,stp1,stp12,stp2,intV]=logReader(app,path);
if isempty(app.SavePath);return;end
fprintf(['Path : ' path ' \n']);
%Step information
app.FirstStepEF.Value=stp1;
app.ConsoStep=stp12;
app.LastStepEF.Value=stp2;
app.IntervalEF.Value=intV;
%create the TrialData and check 3D or 2D
app=FirstValues(app);   

%Necessary values
app.ExeAllButton=1;
app.LIGGGHTSAnalysisButtonGroup.SelectedObject=app.ExeAllButton;
if isa(N1,'double')
    if N1==0 %#ok<*BDSCI>
        N1=app.FirstStepEF.Value;
        fprintf("N1 defined as "+N1+"\n");
    end
    if N2==0
        N2=app.LastStepEF.Value;
        fprintf("N2 defined as "+N2+"\n");
    end
    app.ExtFSwitch.Value='Compression';
else
    switch N1
        case 'cs'
            fprintf('Executing consolidation part only\n');
            N1=app.FirstStepEF.Value;
            N2=app.ConsoStep;
        case 'cp'
            fprintf('Executing compression part only\n');
            N1=app.ConsoStep;
            N2=app.LastStepEF.Value;
        case 'dr'
            if app.SimType==3
                fprintf('Executing after drained part only\n');
                N1=app.QcstStep(1);
                N2=app.LastStepEF.Value;
            else
                fprintf('Executing compression part only\n');
                N1=app.ConsoStep;
                N2=app.LastStepEF.Value;
            end
        otherwise
            fprintf('N1 was badly defined, please \n');
            return
    end
end
if N1==0 || N2==0
    fprintf('Problem in step definition, N1=%d, N2=%d\n',N1,N2);return
end
app.N1EF.Value=N1;
app.N2EF.Value=N2;
if interval==0
    interval=app.IntervalEF.Value;
elseif interval<20
    interval=interval*app.IntervalEF.Value;
else
    a=ceil(interval/app.IntervalEF.Value)*app.IntervalEF.Value;
end
app.CalcInt.Value=interval;
fprintf('Calculation interval equal %d \n',interval);

%Partial and Subdivision
if Partial
    switch rectname
        case 'Triaxial'
            PD=partialData('DATA','',app,[0,0,app.boxYEF.Value,...
                app.boxZEF.Value],app.boxXEF.Value,0);
        case 'Value'
            PD=partialData('DATA','',app,[rectInfo(1),rectInfo(2),...
                rectInfo(3),rectInfo(4)],rectInfo(5),rectInfo(6));
        otherwise
            
    end
    app.PartialSwitch.Value='Partial'; %#ok<*UNRCH>
    app.SubdivisionButton.Value=Sub;
    if Sub
        app.ColEF=Columns;
        PD.SubCol=Columns;
        app.LinesEF=Lines;
        PD.SubLin=Lines;
        PD = subPoints(PD);
    end
else
    app.PartialSwitch.Value='Total';
    app.SubdivisionButton.Value=0;
    PD='';
end
app.PartialData=PD;

%% Execution --- modify with caution ---
%Using the created object the functions will be lauched following the
%parameters defined by the user in the first part.Modifications in this
%part may break this function

%START CALCULATIONS
fprintf('Start Calculations\n');


%Loops - start with loops as other calc may depend on this results
    %1 for VTK files, 2 deformability, 3 for Cl Transformation, 4 for
    % Stress, 5 for Anisotropy and 6 Void ratio.
if isa(Loops,'cell') %{"VTK","DEF","TRF","ANI","STN","STS","VR"}
    app.LPAniSwitch.Value='Off';
    app.LPDefSwitch.Value='Off';
    app.LPCTrSwitch.Value='Off';
    app.LPStrSwitch.Value='Off';
    app.LPVRSwitch.Value='Off';
    app.LPStnSwitch.Value='Off';
    vtk='Off';
    opt='';
    if sum(strcmpi(Loops,'VTK'));vtk='On';opt=' Vtk';end
    if sum(strcmpi(Loops,'DEF'));app.LPDefSwitch.Value='On';opt=[opt ' Defbt'];end
    if sum(strcmpi(Loops,'TRF'));app.LPCTrSwitch.Value='On';opt=[opt ' ClstTr'];end
    if sum(strcmpi(Loops,'ANI'));app.LPAniSwitch.Value='On';opt=[opt ' Ani'];end
    if sum(strcmpi(Loops,'STN'));app.LPStnSwitch.Value='On';opt=[opt ' Stn'];end
    if sum(strcmpi(Loops,'STS'));app.LPStrSwitch.Value='On';opt=[opt ' Str'];end
    if sum(strcmpi(Loops,'VR' ));app.LPVRSwitch.Value='On';opt=[opt ' VR'];end
    app.LPPctEF.Value=99;
    fprintf('Executing Loops\n');
    if ~isempty(opt)
        fprintf(['Extra Options :' opt '\n']);
    end
    exe_LoopsCounter(app,PD,vtk);
end

%Internal Stress
if IntForces
    fprintf('Executing Internal Stress\n');
    exe_StressTensor(PD,app);
end

%Strain Tensor - second strain as other calc may depend on this results
if StrainTensor~=0
    type='';
    if ismember(1,StrainTensor)
        type='GLOBAL';
    elseif ismember(2,StrainTensor)
        type='PERCELL';
    elseif ismember(3,StrainTensor)
        type='BOTH';
    end
    if ismember(StrainTensor,4)
        app.StTW2ClCB.Value=1;
    else
        app.StTW2ClCB.Value=0;
    end
    if ~isempty(type)
        fprintf('Executing Strain\n');
        fprintf(['Options : ' type '\n']);
        exe_StrainTensor(app,type,PD);
    end
end

%Anisotropy
if Anisotropy
    fprintf('Executing Anisotropy\n');
    exe_Anisotropy(PD,app);
end

%External Stress
if ExtForces
    fprintf('Executing External Stress\n');
    exe_ExtForces('',app);
end

%EdgeRatio
if EdgeRatio
    fprintf('Executing Edge Ratio\n');
    exe_EdgesRatio(app);
end

% Void Ratio
if VoidRatio
    fprintf('Executing Void Ratio\n');
    exe_VoidRatio(PD,app);
end

if isa(ForceChain,'cell')
    %check number of elements of ForceChain
    %{'BASE','BEND','CLTF','CLST','GRN',[3,10,15],'CNST',[30,1.1]}
    exFC=0;
    opt='';
    if any(strcmpi(ForceChain,'BASE'))
        of=1;exFC=1;opt=[opt ' BASE'];
    else
        of=0;
    end
    app.FCBaseCBox.Value=of;
    if any(strcmpi(ForceChain,'CL'))
        of=1;exFC=1;opt=[opt ' CL'];
    else
        of=0;
    end
    app.FCClusterCBox.Value=of;
    if any(strcmp(ForceChain,'CLTF'))
        opt=[opt ' CLTF'];of=1;
    else
        of=0;
    end
    app.FCBendCLTCBox.Value=of;
    chk=strcmpi(ForceChain,'GRN');
    if any(chk)
        exFC=1;opt=[opt ' GRN'];
        v=ForceChain{find(chk)+1};
        app.FCGrVtkCBox.Value=1;
        app.FCGrVtkEF.Value=sprintf(['[' repmat('%d,',1,numel(v)-1) '%d]'],v);
    else
        app.FCGrVtkCBox.Value=0;
    end
    chk=strcmpi(ForceChain,'CNST');
    if any(chk)
        prm=ForceChain{find(chk)+1};
        if numel(prm)~=2
            fprintf('Wrong def of Force Chain parameters\n');exFC=0;
        else
            app.FCAngEF.Value=prm(1);
            app.FCForceCoefEF.Value=prm(2);
        end
        exFC=1;opt=[opt ' CNST ang=' num2str(prm(1))...
            ' Coef=' num2str(prm(2)) ' '];
    else
        app.FCAngEF.Value=45;
        app.FCForceCoefEF.Value=1;
    end
    chk=strcmpi(ForceChain,'BNDANGLE');
    if any(chk)
        app.FCBendAngleEF.Value=ForceChain{find(chk)+1};
    else
        app.FCBendAngleEF.Value=0.5;
    end
    if exFC
        fprintf('Executing Force Chain\n');
        if ~isempty(opt)
            fprintf(['Extra Options :' opt '\n']);
        end
        exe_ForceChains(app)
    end
end

%Basic Information
if BasicInfo
    fprintf('Executing Basic Information\n');
    exe_BasicInfo(app)
end

end
%% Sup  port Functions --- modify with caution ---
%Following funcitons are used to simulate the existence of the matlab
%aplication. It will create an 'app' object with the necessary properties
%for the well execution of the functions.
function app=FirstValues(app)
%FIRSTVALUES read pistons information, check the dimension
%   This function will look for the file containing the
%   information about the displacement of the pistons. It will
%   then create a 'trialData' object and save it as a property.
%   Also it will run the 'grains' class once to find out if the
%   simulation was made into 2D or 3D. This calculation is made
%   inside this class because we can analyse grains positions.

TD=trialData(app);
if isempty(TD.Step)
    warndlg(['Files containing servo data could not be found. ' ...
        'Check Settings values to make sure the app will work']);return;
end
app.checkPiston=TD.checkPiston;
app.TrialData= inflectionPoints(TD,app);

%Check first and last step - take into account the case where
%consolidation data is missing
if app.FirstStepEF.Value<TD.Step(1)
    app.FirstStepEF.Value=TD.Step(1);
end
if app.LastStepEF.Value>TD.Step(end)
    app.LastStepEF.Value=TD.Step(end);
end

%check if 3D or 2D. The calculation is made inside de grain class,
%So we just call it once.
gr=grains('BASIC',app.FirstStepEF.Value,'',app);
if isempty(gr.Coord);return;end
%Check dimesions - check if 3D or 2D - diference between X
%coordinates must be smaller then the max radius for 2D
if abs(max(gr.Coord(:,1))-min(gr.Coord(:,1)))>max(gr.Radius)
    app.Bool3D=1; %3D
else
    app.Bool3D=0; %2D
end
app.NbGrains=gr.Nb;
app.Figsize=[560,420];
end
