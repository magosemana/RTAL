function [N1,N2,interval,stepArray,nbFiles] = createStepArray(app)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
N1=app.N1EF.Value;
N2=app.N2EF.Value;
interval=app.CalcInt.Value;
% if interval==app.IntervalEF.Value && app.SimType==3
%     qst=app.QcstStep(2);
%     if qst>N2
%         stepArray=(N1:interval:N2)';
%     elseif qst<=N1
%         f=find(app.TrialData.Step==N1 | app.TrialData.Step==N2);
%         stepArray=app.TrialData.Step(f(1):f(2));
%     else
%         stepArray=(N1:interval:qst)';
%         f=find(app.TrialData.Step==qst | app.TrialData.Step==N2);
%         stepArray=[stepArray(1:end-1);app.TrialData.Step(f(1):f(2))];
%     end
% else
    stepArray=(N1:interval:N2)';
    if stepArray(end)~=N2; stepArray=[stepArray;N2];end
% end
nbFiles=numel(stepArray);
end