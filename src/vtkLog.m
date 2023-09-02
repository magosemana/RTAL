function vtkLog(app,type,fnm,StepEzEvQP)
%VTKLOG create a log containing data related to the vtk files created
% Upon checking the images in paraview, it is important to know at wich
% state the specimen is in. Thus values lique Ez, Ev, p and q will be
% exported in a table to be placed next to the vtk files, thus making it
% easier to identify these points.

fid = fopen(fnm, 'w');
fprintf(fid, '## Table containing the VTK files specimen information ##\n');
if app.SimType==1
    tp='Drained';
elseif app.SimType==2
    tp='Drained';
else
    tp='Qcst';
end
fprintf(fid, ['Simulation type :' tp '\n']);
fprintf(fid, ['VtK type :' char(type) '\n']);
fprintf(fid,'NoImage|Step|Ez|Ev|q|p\n');
fprintf(fid,['%d' repmat('|%d',1,5) '\n'],[1:size(StepEzEvQP,1);StepEzEvQP']);
fclose(fid);

end

