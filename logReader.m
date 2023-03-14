function [app,stp1,stp12,stp2,intV]=logReader(app,fPath,fName)
%LOGREADER read log file and save important variables
%   This function will read the 'LiggghtstoMatlab.txt' file that is
%   created after the simulation. The data that interest us is
%   situated at the beggining at the file. The data will be
%   identified and saved on the correct position.

if nargin==2
    filenm=[fPath 'LiggghtstoMatlab.txt'];
else
    filenm=[fPath fName];
end
try opts = detectImportOptions(filenm,'FileType','text','NumHeaderLines',2,'Delimiter','=');
catch
    warndlg("cound not find file "+filenm);
    app.SavePath='';
    return;
end
%Save Path
app.SavePath=fPath;
%load file
opts.VariableTypes={'char', 'char'};
M = readtable(filenm,opts); %1-11
sz=min(30,size(M,1));
dV=zeros(1,5); %position of piston before compression
dVO=zeros(1,5); %position of piston before consolidaiton
for i=1:sz
    switch cell2mat(M{i,1})
            %Geometric values
        case "Width";app.boxXEF.Value=str2double(cell2mat(M{i,2}));
        case "Length";app.boxYEF.Value=str2double(cell2mat(M{i,2}));
        case "Height";app.boxZEF.Value=str2double(cell2mat(M{i,2}));
        case "ExpWidth";app.boxXEF.Value=str2double(cell2mat(M{i,2}));
        case "ExpLength";app.boxYEF.Value=str2double(cell2mat(M{i,2}));
        case "ExpHeight";app.boxZEF.Value=str2double(cell2mat(M{i,2}));
            %displacement values
        case "dw1";dV(1)=str2double(cell2mat(M{i,2}));
        case "dw2";dV(2)=str2double(cell2mat(M{i,2}));
        case "dl1";dV(3)=str2double(cell2mat(M{i,2}));
        case "dl2";dV(4)=str2double(cell2mat(M{i,2}));
        case "dh";dV(5)=str2double(cell2mat(M{i,2}));
        case "dw1O";dVO(1)=str2double(cell2mat(M{i,2}));
        case "dw2O";dVO(2)=str2double(cell2mat(M{i,2}));
        case "dl1O";dVO(3)=str2double(cell2mat(M{i,2}));
        case "dl2O";dVO(4)=str2double(cell2mat(M{i,2}));
        case "dhO";dVO(5)=str2double(cell2mat(M{i,2}));
            %Step and time
        case "Timestep"; app.TimeStep=str2double(cell2mat(M{i,2}));
        case "StartCons";stp1=str2double(cell2mat(M{i,2}));
        case "EndCons";stp12=str2double(cell2mat(M{i,2}));
        case "EndComp";stp2=str2double(cell2mat(M{i,2}));
        case "Intv";intV=str2double(cell2mat(M{i,2}));
        case "QcstInit";qcst(1)=str2double(cell2mat(M{i,2}));
        case "QcstRupt";qcst(2)=str2double(cell2mat(M{i,2}));
            %folder Details
        case "vtkFolder";app.VtkFolder{1}=[fPath cell2mat(M{i,2})];     %VTK files folder
        case "genFolder";app.DataFolder{1}=[fPath cell2mat(M{i,2})];     %general files folder
        case "vtkCnsF";app.VtkFolder{1}=[fPath cell2mat(M{i,2})];
        case "vtkCmpF";app.VtkFolder{2}=[fPath cell2mat(M{i,2})];
        case "resCnsF";app.DataFolder{1}=[fPath cell2mat(M{i,2})];
        case "resCmpF";app.DataFolder{2}=[fPath cell2mat(M{i,2})];
            %execution details
        case "ExecutionType";app.ExeType=str2double(cell2mat(M{i,2}));
        case "SimType";st=cell2mat(M{i,2});
    end
end
if numel(app.VtkFolder)==1
    %update settings if there is only 1 element
    app.VTKpathEF.Value=app.VtkFolder{1};
    app.GVpathEF.Value=app.DataFolder{1};
else
    stp2=stp2+stp12;
end
app.ConsoDpl=[dVO;dV];
app.QcstStep='';
if ~exist('st','var')
    app.SimType=1;
else
    switch upper(st)
        case {'ISOCHORIC','UNDRAINED'}
            app.SimType=2;
        case 'QCST'
            app.SimType=3;
            app.QcstStep=qcst;
        otherwise
            app.SimType=1;
    end
end

%matlab default colors to mantain every class in the same clor
app.PlotColors=[0 0.4470 0.7410;
    0.8500 0.3250 0.0980;
    0.9290 0.6940 0.1250;
    0.4940 0.1840 0.5560;
    0.4460 0.6740 0.1880;
    0.3010 0.7450 0.9330;
    0.6350 0.0780 0.1840];

end