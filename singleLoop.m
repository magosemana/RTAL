classdef singleLoop
    %SINGLELOOP object containing information of Loops and Clusters
    %   This object is formed by calculating the Loops inside a spaceCell
    %   system.
    
    properties
        %All loop properties
        ClEdges         %Closed edges belonging to the Loop - Nx2 matrix
        Vertices        %IDs of grains forming the limits (closed edg or surf) of the loop or cluster, 
                        %separed in trios for surface plot - Nx3 matrix
        Volume          %total Loop volume/area - double
        sCells          %ID of spaceCells forming the loop - Nx1 vector
        nbCells         %Nb of spaceCells forming the loop - scalar
        Size            %Number of grains - double
        Order           %Number of exterior surfaces/edges - double
        Deformability   %Nb O edges / N tot edges
        
        
        %3D only properties
        Grains          %Ids of grains forming the loops (only 3D) - Nx1 vector
        NbGrains        %Number of grains forming the loop - double
        Z               %Loop coordination nbr
        
        %Extra properties
        Stress          %Loop Stress - double
        iStrain         %Loop incremental Strain - double
        Transformation  %Loop Void Ratio - Ncellsx2 
        VoidRatio       %Loop Void Ratio - double
    end
    
    methods
        function sl = singleLoop(D,sc,clID,sCarray,varargin)
            %SINGLELOOP Construct an instance of this class
            %   This function saves the caracteristics of the loops
            %   obtained through the calculation. 
            %
            %   -Vertices will contain each triangle, wich is formed by a
            %   trio of vertices that form the delaunayTesselation. It was
            %   chosen this way because it is easy for drawing.
            %   -Cells contain the ID of the spaceCells forming the loops
            %   -Order contains the loopsize in 2D, and the nb of virtual
            %   surfaces in 3D
            if ~nargin %create a empty object
                sl.Order=0;
                sl.sCells =0;
                return
            elseif D==4
                %Shortcut - the input names mean nothing, profiting of the
                %existence of this class to create an object containing all
                %cluster 4 data
                sl.Order=D;
                sl.nbCells=sc;
                sl.sCells=clID;
                sl.Z=sCarray;
                return
            end
            %Vertices Ids,Cell Ids,Loopsize,AnisotropyEdges,Closed Surfaces
            %only on 3D. 
            if D==3
                %Two pathologies may appear on cluster search. "pseudo
                %closed frontiers" and "phantom cluster 4 "
                
                %Phatology 1 : it may happen that two cells that share
                %a closed surface are inside a Cluster. This is then a
                %"pseudo closed surface" that must be removed from the
                %calculation. 
                cS=sort(cat(1,sCarray.ClosedS),2);
                if length(cS)>200
                   a=1; 
                end
                [un,~,cnt]=unique(cS,'rows'); %get unique
                rpCE=accumarray(cnt,1);
                cS=un(rpCE==1,:); %keep not repeating values
                
               
                %Pathology 2 : sometimes, Cl4 of very low volume appear
                %inside other clusters. These respect the 'closed
                %frontiers' definition but when analysed mechanicaly should
                %not exist. Thus they have to be removed and incorporated
                %into the cluster. 
                    %First, get from the DT cells containg ALL grains
                grs=unique(sc.DelaunayT(clID,:));
                if size(grs,2)>1;grs=grs';end %force it to be a column vector
                 %{
                    %Identify all cells containing ALL these grains
                grC=find(sum(ismember(sc.DelaunayT,grs),2)==4);
                    %check if any of those is a Cl4
                is4=find(ismember(grC,sc.Clt4.sCells));
                
                for i=1:numel(is4)
                    %check if all surfaces from the cl4 correspond to
                    %closed surfaces from this cluster
                    grCl4=sc.DelaunayT(grC(is4(i)),:);
                    sf=sort(grCl4(nchoosek(1:4,3)),2);
                    [mb,mb2]=ismember(sf,cS,'rows');
                    if sum(mb==4)
                        %if all surfaces match closed surfaces, add the 
                        %ID of the Cl4 to the cluster
                        clID=[clID,grC(is4(i))]; %#ok<AGROW>
                        cS=cS(~mb2); %turn into "pseudo-closed surfaces"
                    end
                end
                %}
                
                %Cluster order calculation : (D+1)*NbCell - 2*NbOFront or
                %nb of closed surfaces
                %clOd=(D+1)*numel(cltI)...
                %    -sum([sCarray.NbOCon],'all')-2*sum(reapCE>1);
                clOd=size(cS,1);
                
                %Closed edges : for Z calculation
                cE=unique(cat(1,sCarray.CEdges),'rows');
                nbCE=size(cE,1);
                
                %Total unique edges : get all edges from all cells,
                %calculate the uniques to get the amount. Value used for
                %deformability indice
                edges=varargin{1};
                nbTot=reshape(permute(edges,[3 1 2]),[],2); % transform 3D matrix int 2D keeping the correct edges
                nbTot=size(unique(nbTot,'rows'),1);
                
                %Save values needed on calculations
                sl.ClEdges=cE; %unique closed edges
                sl.Grains=grs;      %unique grains
                sl.Order=clOd;      %cluster Order
                sl.sCells = clID;   %ID of cells inside the cluster
                sl.nbCells=numel(clID); %nb of cells
                sl.Size=numel(grs);   %Loop size (nb of grains)
                sl.Vertices=cS;     %vertices (unique closed surafeces used for drawing)
                sl.Z=2*nbCE/sl.Size;%coord nb
                sl.Deformability =1 - nbCE/nbTot; %Deformability =OE/Tot
            else
                %2D Loops
                sl.Vertices=varargin{1};            %vertices for drawing
                sl.Size=(D+1)*numel(loopi)...
                    -sum([sCarray.NbOCon],'all');   %size calculation
                sl.Order=sl.Size;                   %order 
                sl.ClEdges=unique(cat(1,sCarray.CEdges),'rows'); %unique closed edges

            end
            
        end
        function sl = loopStrain(sl,cellStn,cellVol)
            %takes the strain from each cell belonging to the Cluster and
            %calcuates the value of the cluster strain.
            sl.iStrain=sum(cellStn,3)/sum(cellVol);
            sl.Volume=sum(cellVol);
        end
    end
end

