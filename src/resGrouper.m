function resGrouper(tgPath,orPath,mode)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

app.SavePath=orPath;
if~isfolder(tgPath);mkdir(tgPath);end
apptgt.SavePath=tgPath;
nFnm=split(orPath,'/');
nFnm=nFnm(end);
fprintf("Executing "+nFnm+"\n")
for j=1:numel(mode)
    %get directory and check if it exists
    p=MakePath(app,mode(j),'');
    if ~isfolder(p)
        fprintf("Folder of "+mode(j)+" was not found for "+nFnm+"\n")
        continue;
    end
    f=dir(p(1)); %get files inside of it
    f([f.isdir]) = [];  %remove folders types
    for k=1:numel(f)
        if strcmp(f(k).name(end-3:end),'.mat') ||...
                strcmp(f(k).name(end-3:end),'.txt')
            pathDs=MakePath(apptgt,mode(j));
            fl=fullfile(p,f(k).name);
            copyfile(fl,pathDs)
            movefile(fullfile(pathDs,f(k).name),...
                fullfile(pathDs, nFnm + f(k).name(end-3:end)) );
        end
       
    end
end
fprintf("Finished "+nFnm+"\n")

end

