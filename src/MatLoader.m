function [fnm,fPath]=MatLoader(type,app,multi)
%Matloader loads file for other funtions. It will behave diferently
%depending on the "type" being called.

%Define values that are function dependant: file location, number of
%columns and multi-load on/off.
bsPth=MakePath(app,type,'check');
if exist(bsPth,'dir')==0
    fprintf('Usual path not found \n')
    bsPth='';   
end
%check if multi variable was called, else chose 'on' as default value
if nargin<3;multi='on';end
[fnm,fPath,~] = uigetfile( ...
    {'*.mat','Matlab files(*.mat)';'*.*','All Files (*.*)'},...
    'Select file for loading','Multiselect',multi,bsPth);

%If there is only 1 number in each of the names of files, change their
%order for the crescent one.
nnm=numel(fnm);
if iscell(fnm) && nnm>1 
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

end

