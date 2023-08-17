function pD=FileLoader(app,type)
%FILELOADER loads file for other funtions. It will behave diferently
%depending on the "type" being called.

%Define values that are function dependant: file location, number of
%columns and multi-load on/off.
switch upper(type)
    case "ANI";multi='on';nbCol=[9,10];
    case "CONTACT";multi='off';nbCol=4;
    case "INTERNALFORCES";type='Forces';multi='on';nbCol=[8,13];
    case "EDG";multi='on';nbCol=0;
    case "LOOPS";multi='on';nbCol=0;
    case "LOOPVR";multi='on';nbCol=0;
    case "VOID";multi='on';nbCol=5;
    case "STRAIN";multi='on';nbCol=[5,6];
    case "FC";multi='on';nbCol=9;
end

bsPth=MakePath(app,type,'check');
if exist(bsPth,'dir')==0
    fprintf('Could not find the usual path \n');
    bsPth=MakePath(app);
end
%Ask for file
[fnm,fPath,~] = uigetfile( ...
    {'*.txt','Text files(*.txt)';'*.*','All Files (*.*)'},...
    'Select file for loading','Multiselect',multi,bsPth);
%Check if file is not empty
if fPath==0;pD='';return;end

switch class(fnm)
    case 'char' %one file selection
        %Load Results
        opts = detectImportOptions(fullfile(fPath,fnm),'NumHeaderLines',8);
        opts.Delimiter='|';
        opts.VariableNamesLine = 9;
        Table = readtable(fullfile(fPath,fnm),opts);
        
        for i=1:size(Table,2)%count how many columns of data we have
            if isnan(Table{1,i}) || ~isequal(class(Table{1,i}),'double')
                i=i-1;break; 
            end
        end
        if ~ismember(i,nbCol) && nbCol>0
            k=warndlg("Could not correclty determinate the"+...
                " columns of the file");
            pD='';waitfor(k);return;
        end
        file=Table{1:end,1:i};

        %Load other information
        opts2 = detectImportOptions(fullfile(fPath,fnm),'NumHeaderLines',1);
        opts2.Delimiter='|';
        opts2.DataLines=[3 7];
        opts2.VariableNames = {'A','B','C'};
        opts2 = setvartype(opts2,{'double'});
        Table2 = readtable(fullfile(fPath,fnm),opts2,'ReadVariableNames',false);
        file2=Table2{1:5,1:3};
        
        %remove format from the file name
        fnm=fnm(1:(find(ismember(fnm,'.'),1,'last')-1));
            
        %create PLOTDATA object
        pD = plotData("Load",file,file2(1,1),file2(1,2),file2(1,3),...
            file2(2,1),file2(2,2),file2(2,3),...
            file2(3,1),file2(3,2),...                   %consopt information
            [file2(3,3),file2(4,1),file2(4,2),file2(4,3)],... %inflection pts
            file2(5,2),file2(5,3),fnm,fPath);     %subdivision information
          
    case 'cell' %many file selection
        NbFiles=size(fnm,2);
        %Check all file names to see if they match, there cant be any 'sub'
        %files in the middle of the others
        for i=1:NbFiles
            if sum(strfind(upper(fnm{i}),'subd'))~=0
                k=warndlg('Subdivision files dont work on multi selection');
                pD='';waitfor(k);return;
            end
        end
        
        %If there is only 1 number in each of the names of files, change 
        %their order for the crescent one.
        nnm=numel(fnm);
        if nnm>1
            fail=0; %checker to see if it is needed or not
            nb=zeros(nnm,1);
            for i=1:nnm
                f=fnm{i};
                num = regexp(f,'(-)?\d+(\.\d+)?(e(-|+)\d+)?','match');
                if numel(num)==1
                    nb(i)=str2double(num{1});
                else
                    fail=1;
                end
            end
            if fail==0
                [~,pos]=sort(nb);
                fnm=fnm(pos);
            end
        end

        %Start creating the cells containing the matrix with data
        pD=plotData.empty(NbFiles,0);
        for i=1:NbFiles
            %Load file
            opts = detectImportOptions(fullfile(fPath,fnm{i}),'NumHeaderLines',8);
            opts.Delimiter='|';
            opts.VariableNamesLine = 9;
            opts.VariableNamingRule='preserve';
            Table = readtable(fullfile(fPath,fnm{i}),opts);
            
            for j=1:size(Table,2)%count how many columns of data we have
                if ~isequal(class(Table{1,j}),'double')
                    j=j-1;break;  %#ok<*FXSET>
                end
                if isnan(Table{1,j})
                    j=j-1;break;  %#ok<*FXSET>
                end
            end
            if ~ismember(j,nbCol) &&  nbCol>0
                k=warndlg("Could not correclty determinate the"+...
                    " columns of the file");
            	pD='';waitfor(k);return;
            end
            file=Table{1:end,1:j};
            
            %Load other information
            opts2 = detectImportOptions(fullfile(fPath,fnm{i}),'NumHeaderLines',1);
            opts2.Delimiter='|';
            opts2.DataLines=[3 7];
            opts2.VariableNames = {'A','B','C'};
            opts2 = setvartype(opts2,{'double'});
            Table2 = readtable(fullfile(fPath,fnm{i}),opts2);
            file2=Table2{1:5,1:3};
            
            %remove format from the file name
            fn=fnm{i};
            fn=fn(1:(find(ismember(fn,'.'),1,'last')-1));
            
            %create PLOTDATA object
            pD(i) = plotData("Load",file,file2(1,1),file2(1,2),file2(1,3),...
                file2(2,1),file2(2,2),file2(2,3),...
                file2(3,1),file2(3,2),...                   %consopt information
                [file2(3,3),file2(4,1),file2(4,2),file2(4,3)],... %inflection pts
                file2(5,2),file2(5,3),""+fn,fPath);     %subdivision information
        end
        
end

end

