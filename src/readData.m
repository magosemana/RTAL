function rData = readData(type,app,step)
%READDATA Contain scripts to read LIGGGHTS data
%   Detailed explanation goes here

if numel(app.DataFolder)==2 && step>app.ConsoStep
    pth=app.DataFolder{2};
    step=step-app.ConsoStep;
else
    pth=app.DataFolder{1};
end
       
rData='';
switch upper(type)
    case 'GRAINS'
        fnm=fullfile(pth,app.GrainsEF.Value+string(step)+app.FormatEF.Value);
        %start reading file
        try opts = detectImportOptions(fnm);
        catch
            warndlg("ERROR IN GRAINS :File "+fnm+" could "+...
                "not be found. Try changing inputs or the name of "+...
                "the folder the files are found");
            return;
        end
        opts.DataLines = 10;
        opts.VariableNames={'ID','x','y','z','Radius','Z','ux','uy','uz',...
            'sx','sy','sz','sxy','sxz','syz','keng'};
        tab = readtable(fnm,opts);
        try rData=tab{:,:};
        catch ME
            if (strcmp(ME.identifier,'MATLAB:table:ExtractDataIncompatibleTypeError'))
                opts = detectImportOptions(fnm);
                opts.DataLines = 10;
                opts.VariableNames={'ID','x','y','z','Radius','Z','ux','uy',...
                    'uz','sx','sy','sz','sxy','sxz','syz','omgx','omgy','omgz'};
                tab = readtable(fnm,opts);
                rData=tab{:,:};
            else
               rethrow(ME) 
            end
        end
        rData=sortrows(rData);
    case 'GRAINSTEST'
        fnm=step;
        %start reading file
        try opts = detectImportOptions(fnm);
        catch
            warndlg("ERROR IN GRAINS :File "+fnm+" could "+...
                "not be found. Try changing inputs or the name of "+...
                "the folder the files are found");
            return;
        end
        opts.DataLines = 10;
        opts.VariableNames={'ID','x','y','z','Radius','Z','ux','uy','uz',...
            'sx','sy','sz','sxy','sxz','syz','keng'};
        tab = readtable(fnm,opts);
        rData=tab{:,:};
        rData=sortrows(rData);
    case 'CONTACT'
        fnm=fullfile(pth,app.GrainsForceEF.Value+...
            string(step)+app.FormatEF.Value);
        try opts = detectImportOptions(fnm);
        catch
            warndlg("ERROR IN CONTACT :File "+fnm+" could not "+...
                "be found. Try changing inputs or the name of the "+...
                "folder the files are found"); 
            return;
        end
        opts.DataLines = 10;
        opts.VariableNames={'ID1','ID2','Fx','Fy','Fz','Ctpx','Ctpy',...
            'Ctpz','Overlap'};
        tab= readtable(fnm,opts);
        rData=tab{:,:};
        
    case 'WCONTACT'
        fnm=fullfile(pth,app.GrainsForceEF.Value+"W"+...
            string(step)+app.FormatEF.Value);
        try opts = detectImportOptions(fnm);
        catch
            warndlg("ERROR IN CONTACT :File "+fnm+" could not "+...
                "be found. Try changing inputs or the name of the "+...
                "folder the files are found"); 
            return;
        end
        opts.DataLines = 10;
        opts.VariableNames={'ID1','ID2','Fx','Fy','Fz','Ctpx','Ctpy',...
            'Ctpz'};
        tab= readtable(fnm,opts);
        rData=tab{:,:};
        
    case 'PISTONS'
        fnm=fullfile(pth,app.PistEF.Value +""+ app.FormatEF.Value);
        %load file
        try opts = detectImportOptions(fnm);
        catch
            warndlg(['File could not be found. Try changing inputs'...
                ' or the name of the folder they files are found']);
            return;
        end
        if size(opts.VariableTypes,2)==7
            opts.VariableNames={'Step','Fz','Dz','Fy1','Dy1','Fy2','Dy2'};
        else
            opts.VariableNames={'Step','Fz','Dz','Fy1','Dy1','Fy2','Dy2',...
                'Fx1','Dx1','Fx2','Dx2'};
        end
        opts.DataLines = 3;
        tab = readtable(fnm,opts);
        rData=tab{:,:};
        rData=sortrows(rData);
        %chech the existence a second piston file (Qcst after rupture)
        fnm=fullfile(pth,app.PistEF.Value +"2"+ app.FormatEF.Value);
        try opts = detectImportOptions(fnm);
        catch
            return;
        end
        if size(opts.VariableTypes,2)==7
            opts.VariableNames={'Step','Fz','Dz','Fy1','Dy1','Fy2','Dy2'};
        else
            opts.VariableNames={'Step','Fz','Dz','Fy1','Dy1','Fy2','Dy2',...
                'Fx1','Dx1','Fx2','Dx2'};
        end
        opts.DataLines = 3;
        tab = readtable(fnm,opts);
        rData2=tab{:,:};
        rData2=sortrows(rData2);
        %Check if there is a repeated step
        if rData2(1,1)==rData(end,1)
            rData=[rData(1:end-1,:);rData2];
        else
            rData=[rData;rData2];
        end
    case 'PISTONSYADE'
        fnm=fullfile(pth,app.PistEF.Value +""+ app.FormatEF.Value);
        %load file
        try opts = detectImportOptions(fnm);
        catch
            warndlg(['File could not be found. Try changing inputs'...
                ' or the name of the folder they files are found']);
            return;
        end
        opts.VariableNames={'Step','Sigx','Sigy','Sigz','Uz','Uy1',...
            'Uy2','Ux1','Ux2'};
        opts.DataLines = 3;
        tab = readtable(fnm,opts);
        rData=tab{:,:};
        rData=sortrows(rData);

    otherwise
        warndlg('Wrong readData Type');return;
end
end

