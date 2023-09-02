classdef forceChain
    %FORCECHAIN object containing information of a grain loop
    %   This object is formed by calculating the ForceChains inside
    %   exe_ForceChains file.
    
    properties
        %Base properties
        AveStress       %Average major principal stress of grains
        Branches        %All single chains inside the force chain
        Elevation       %Elavation angle of each contact in the chain
        IDs             %Grains Id's
        Length          %Force chain length
        Lines           %All contact (easy vtk drawing)
        NbGrains        %Nb of grains
        NbBranches      %Nb of unique chains
        %Cluster - Force chain properties
        ClustersID      %ID of Clusters belonging to fc
        ClustersVals    %Per grain max/min/mean cluster values

    end
    
    methods
        function fc = forceChain(sGr)
            %FORCECHAIN Construct an instance of this class
            %   Force chain object containing all the SingleGrains objects
            %   that form it and its length. The length is calculated
            %   through the method fcLength.
            if ~nargin; fc.IDs =0;return;end
            fc.IDs=cat(1,sGr.ID);
            fc.NbGrains=numel(sGr);
            [fc.Lines,pos]=unique(cat(1,sGr.ContactLine),'rows');
            ang=cat(1,sGr.Elevation);
            fc.Elevation=ang(pos,:);
            fc.AveStress=mean([sGr.PrincipalStress]);
            %calculate the length
            [lngth,branch] = forcechainLength(fc,sGr);
            fc.Length=lngth;
            fc.NbBranches=numel(branch);
            fc.Branches=branch;
        end
        function [lngth,uBranch] = forcechainLength(fc,sGr)
            %FCLENGTH calculates the length of the given force chain
            %   Goes trough the singleGrains objects foming the force chain
            %   searching for the extremities of the force chain (grains
            %   with one contact). If a force chain has exactly two
            %   extremities the length is the number of grains, if it has
            %   more the possible trajects must be found.
            
            %get the number of contacts of each grain. 1 is extremity, 2
            %midle and 3 or more bifurcations.
            if sum(cat(1,sGr.ContactNb)>2)==0
                %no bifurcation, length is the nb of grains
                lngth=numel(sGr);
                uBranch={[fc.IDs]};
                return;
            end
            
            %find the probable starting points of the fc. As contacts IDs
            %are put in order of fc contact direction, grains present in
            %the first column of fc.Lines but not the second begin the
            %chain.
            ids=ismember(unique(fc.Lines(:,1)),unique(fc.Lines(:,2)));
            ids=fc.Lines(~ids,1);
            %Start checking all different paths of the forcechain
            uBranch=cell.empty(0,1); %will contain all unique force chains
            fcLine=1;               %line counter
            lngth=0;
            for i=1:numel(ids)
                %start the line with ID of one starting point
                newChain=ids(i);    
                check=0;
                while check~=1
                    %get next grains that the last element of newChain is
                    %connected to. Look int othe first column of fc.Lines
                    %for the last grain added to newChain.
                    nextV=fc.Lines((fc.Lines(:,1)==newChain(end)),2);
                    %new member cannot be already used in the chain (block
                    %loops)
                    nextV=nextV(~ismember(nextV,newChain));
                    if numel(nextV)~=0
                        if numel(nextV)>1
                            %if bifurcation create a new line in a cell to
                            %be analysed later
                            for k=2:numel(nextV)
                                uBranch{end+1}=[newChain,nextV(k)]; %#ok<*AGROW>
                            end
                        end
                        %add new grain to fc branch and redo loop
                        newChain=[newChain,nextV(1)];
                    else
                        %If there are no more grains to add to the fc
                        %branch. Then save it and check if any bifurcation
                        %points were found.
                        uBranch{fcLine}=newChain;
                        lngth=max(lngth,numel(newChain));
                        if fcLine==numel(uBranch)
                            %if no, leave loop
                            check=1;
                        else
                            %If yes redo the loop until they are all done.
                            newChain=uBranch{fcLine+1};
                        end
                        fcLine=fcLine+1;
                    end 
                end
            end
            %return fc branch and Legnth as results
        end
    end
end

