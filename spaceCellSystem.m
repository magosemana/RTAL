classdef spaceCellSystem
    %SPACECELLSYSTEM this class will be used to calculate the strain tensor
    %in diferent points of the experiment and the loops formed by the grains
    %
    %   This function uses the formulation proposed by BAGI 1996. She
    %   proposes the definition of two system to analyse the granular
    %   material morphology. The first one is the material cell (better
    %   represented by the Voronoi tesselation) from which the stress
    %   tensor definition (Love-Webber) can be obtained. The second one is
    %   the space cell system (based on the Delaunay Tesselation). The
    %   latter one has closed contacts and open contacts that are used to
    %   obtain the displacement gradient tensor and later the strain
    %   tensor.
    %   These definitions are complementary to each other. The first one
    %   takes the grain as the center of the cell, and all points that
    %   belong to the cell are the points that are closer to it's center
    %   than to any other cell. The second one takes the 'void' as the
    %   center of the cell, and connects all the grains that are close to
    %   it than any other cell. If both are drawn together the vertices of
    %   one is the center of the other.
    %
    %%STRAINTENSOR
    %   The calculation will be made for 2D and 3D simulations. Depending
    %   on the dimension the calculation must be made differently and the
    %   meaning of each calculate differs by the order of one dimension.
    %     -For 2D, the grains will be diveded into triangles formed by
    %     grouping 3 vertices. For each triangle, each vertice will be
    %     atributed a vector wich's size is the length of the distance
    %     between the two other vertices of the triangle. This vector is is
    %     then used to calculate the 'Complementary area vector', wich is
    %     used to calculate the strain.
    %
    %     -For 3D, the grains will be diveded into 3D tetrahedrons formed
    %     by grouping 4 vertices. For each tetrahedrons, each vertice will
    %     be atributed a vector wich's size is the area of the surface
    %     formed by the three other vertices of the tetrahedrons. This
    %     vector is is then used to calculate the 'Complementary area
    %     vector', wich is used to calculate the strain.
    %
    %%LOOPS
    %   2D
    %   As said previously each cell of the Space Cell system is centered
    %   in a void and it's vertices are the closest grains to it than to
    %   any other. Using the properties of the Space Cell, it can be shown
    %   that all edges are either closed contacts (that actually physically
    %   exists) or open contacts, connecting grains that are the closest
    %   to a cell's center but not really in touch. Thus, two cells that
    %   share a open contact are actually part of the same 'void' space
    %   in the simulation. If we join all the cells that share open
    %   contacts between them we will have for all individual 'void's in
    %   the simulation all the grains that form the boundary between then
    %   and the others.
    %
    %   This function will then search for all the cells that have matching
    %   open contacts and join them into a big cell, obtaining this way
    %   the ID of all the grains that form it's boundary.
    %
    %   3D
    %   The same principle of 2D calculation apply, but this time the
    %   element that makes frontier between two cells is a surface.
    %   Therefore the definition of open and closed surface needs to be
    %   created. A open surface is a surface formed only by open edges, and
    %   a closed surface must have at least one closed edge.
    
    properties
        %BASE PROPERTIES
        D                   %Dimension the calculation      - integer (2 or 3)
        DelaunayT           %DelaunayTriangulation 'object' using the grain centers
        GoodCells           %Indice of 'spaceCell' used     - Nx1
        NbC                 %Nb of spaceCells               - integer
        NbG                 %Nb of grains                   - integer
        NbGc                %number of good cells           - integer
        Radius              %Grains radius                  - NbGx1 array
        sCells              %'spaceCells' objects           - NbCx1
        Step                %Calculation Step               - integer
        TotalVolume         %Total DT volume/area           - double
        Interval            %Interval between relative calc - integer
        
        %STRAIN TENSOR PROPERTIES
        Displacements       %Displacements of grains        - NbGxD matrix
        EdgeID              %grain Ids forming edges        - NbEx2 matix
        NbE                 %TotalNb of edges between grains- integer
        GrainVolume         %Voronoi volume of each grain   - NbGx1 Vector
            %BAGI96
        CompAVector         %Complementary area vector      - NxD matrix
        GStTensor           %Global Strain tensor           - DxD matrix
        PStTensor           %Per Grain Strain tensor        - NbGxD matrix
            %Linear Interpolation
        CellStn             %Per Cell Strain value          - 3x3xNbGc matrix
        CellVol             %Per Cell Volume                - NbGcx1
        CellGrn             %Grains Forming each cell       - NbGcx4 matrix
        
        %LOOPS PROPERTIES
        Loops               %SingleLoops objects            - object vector
        Clt4                %SingleLoops obj                - object
        BadIds              %Ids of "fakepoints"            - Nx1 vector
        ClosedEdges         %Closed Edges                   - Nx2 matrix
        OpenEdges           %Open Edges                     - Nx2 matrix
        TotalVR             %Void ratio of the specime      - double
    end
    
    methods
        %Calling function
        function sc = spaceCellSystem(mode,step,gr,varargin)
            %SPACELL Construct an instance of this class
            %   Construct the fist part of the space cell by creating the
            %   delaunay tesselation. Also it atributes properties that
            %   will later be used. These are the number of grains (NbG),
            %   nuber of spaceCells (NbC), calculation's dimension (D) and
            %   the each grains displacement between both steps.
            
            %Check the arguments given upon calling the function
            if ~nargin
                sc.D=0;
                sc.NbC=0;
                sc.NbG=0;
                sc.Step=-1;
                return;
            end
            switch upper(mode)
                case {"GLOBAL","PERCELL","BOTH"}
                    %strain calculation
                    n=2;
                case {"LOOPS","EDGES","VR"}
                    %Cluster/loops calculation
                    n=1;
                otherwise
                    if nargin==6;n=2;s="Strain";else;n=1;s="Loops";end
                    fprintf("No spaceCellSystem mode chosen, advancing"...
                        +" as "+s+" format. May cause errors\n")
            end
            %check app
            app=varargin{n};
            if app.Bool3D;D=3;else; D=2;end %Check if 3d or 2D
            vert=gr.Coord(:,(4-D):3);
            PD=varargin{n+1};
            if n==2
                %grain object at next interval used to calculate the
                %displacement of each particle.
                grTdT=varargin{1};
            end
            %To better identify and make the connection with actual grains
            %better, 'fakepoints' that follow the geometry of the walls
            %will be added to the delaunay triangulation. These points will
            %later have to be removed from all calculations.
            pts=wallFakepoints(app.TrialData,app,step,max(gr.Radius),PD);
            
            %Delaunay tesselation
            DT=delaunayTriangulation([vert;pts]);
            
            %atribute properties to the object
            sc.D=D;
            sc.DelaunayT=DT;
            sc.NbC=size(DT,1);      %number of cells
            sc.NbG=size(vert,1);    %number of grains
            sc.Step=step;
            %calculate total volume and grains at the exterior of the triangulation
            [~,sc.TotalVolume]=convhulln(vert);
            
            %extra calculations depending on 'mode' calling variable
            switch upper(mode)
                case {"GLOBAL","PERCELL","BOTH"}
                    fprintf('Strain calculation at step %d\n',sc.Step)
                    sc.Radius=[gr.Radius;zeros(size(pts,1),1)];
                    sc.Displacements = [grTdT.Coord(:,(4-D):3)-vert;zeros(size(pts))];
                    sc = pcStrain(sc,gr);
                    sc=strainCells(sc);                 %create spaceCells objects
                    sc=compAreaVector(sc);              %calculate complementary area vectors
                case {"LOOPS"}
                    sc=findOpenEdges(sc,step,app,PD);	%find all open contacts
                    if D==2;sc=findLoops(sc);           %join all cells that form Loops on 2D
                    else;sc=findCluster(sc);            %join all cells that form Loops on 3D
                    end
                case "EDGES"
                    sc=findOpenEdges(sc,step,app,PD);
                case "TEST"
                    fprintf('Test calculation at step %d\n',sc.Step)
                    sc.Radius=[gr.Radius;zeros(size(pts,1),1)];
                    sc.Displacements = [grTdT.Coord(:,(4-D):3)-vert;zeros(size(pts))];
                    %sc=grainVolume(sc);
            end
        end
        %Strain
        function sc = strainCells(sc)
            %strainCells create all single cell objects
            %   When called will create all single cell objects wich are
            %   needed for the strain calculation. The only value that
            %   iterest for this calculation is the surface vector so this
            %   will be the only value assigned to cells
            
            %create a surface combination matrix
            cmb=nchoosek(1:sc.D+1,sc.D);
            extN=zeros(sc.NbC,4,4);
            for i=1:size(cmb,1)
                sf=sc.DelaunayT(:,cmb(i,:));
                %Calculate the normal of the edge/surface between the
                %grains
                n=0.5*cross(sc.DelaunayT.Points(sf(:,2),:)-...
                    sc.DelaunayT.Points(sf(:,1),:),...
                    sc.DelaunayT.Points(sf(:,3),:)-...
                    sc.DelaunayT.Points(sf(:,1),:)); %NbGcx3 matrix
                %Make sure the normal points towards the exterior by
                %analysing the angle between the normal and the a vector
                %from one any grain from the surf/edge to the opposing
                %grain. If the angleis inferior to 90 degrees then n=-n;
                opId=sc.DelaunayT(:,setxor(1:sc.D+1,cmb(i,:))); %opposing ID to the surf/edge
                v=sc.DelaunayT.Points(opId,:)...
                    -sc.DelaunayT.Points(sf(:,1),:); %NbGcx3 matrix
                alpha=acosd(dot(n,v,2)./(vecnorm(n,2,2).*vecnorm(v,2,2)));
                n(alpha<90,:)=-n(alpha<90,:);  %NbGcx3 matrix
                %Save the Opposing grain ID and the normal to be used in
                %the next calculations
                extN(:,:,i)=[opId,n];
            end
            %Create cell object for all goodcells and assign properties
            sCarray(sc.NbC)=spaceCell();  %array to store the 'spaceCell' objects
            [sCarray.SurfaceV]=disperse(permute(extN,[3,2,1]),[1,2]);
            %save only non zero values
            sc.sCells=sCarray;
        end
        function sc = compAreaVector(sc)
            %compAreaVector calculates the complementary area vector
            %   This function will use a previously created 'spaceCell'
            %   objects to calculate the complementary area vector for each
            %   grain contact (closed or open) as defined by BAGI 96. This
            %   vector is used to take into account the position of the
            %   surround grains in relation to each edge 'ij'.
            
            %Good cell vector, find the ID of the tetrahedrons that :
            gC=sum((sc.DelaunayT(:,:))>sc.NbG,2);
            gC=find(gC<4); %contain at least 1 real grain
            %gC=find(gC==0); %contain only real grain cells
            %gC=find(gC<3); %contain at least 1 real edge
            sc.GoodCells=gC;
            sc.NbGc=size(gC,1);
            
            %load edgesID and get only the ones that does contains at least
            %one closed grain
            edg=edges(sc.DelaunayT);
            edg=edg(edg(:,1)<=sc.NbG & edg(:,2)<=sc.NbG,:);
            nbE=size(edg,1);  %nb of edges
            d=zeros(nbE,sc.D);  %complementary area vector for each edge
            for i=1:nbE %for each edge
                %get all the 'spaceCell' objects that share the same edge
                ThIDs=cell2mat(edgeAttachments(sc.DelaunayT,edg(i,1),edg(i,2)));
                ThIDs=find(ismember(sc.GoodCells,ThIDs)); %find position in goodcells vector
                if isempty(ThIDs);continue;end
                %sv cointans the  surface vectors for a 'spaceCell'.
                %The first column being the ID and the others the
                %coordinates. Now the complementary area vector d will
                %be calculated .The definiton of d, for a contact
                %between grains m an n, is d^mn=sum(b^n-b^m).
                sv=cat(3,sc.sCells(sc.GoodCells(ThIDs)).SurfaceV);
                e1=cat(3,sv(:,1,:)==edg(i,1));
                e2=cat(3,sv(:,1,:)==edg(i,2));
                s=reshape(permute(sv(:,2:end,:).*e1,[2,1,3]),sc.D,[],1)'-...
                    reshape(permute(sv(:,2:end,:).*e2,[2,1,3]),sc.D,[],1)';
                d(i,:)=sum(s,1);
            end
            %save found values in the properties
            sc.CompAVector=d/(sc.D*(sc.D+1));
            sc.NbE=nbE;
            sc.EdgeID=edg;
            fprintf('CAV done\n')
        end
        function sc = gStrainTensor(sc,app,PD)
            %GSTRAINTENSOR calculate the global incremental strain tensor
            %   Calculate the symetrical part of the 'Displacement Gradient
            %   tensor' to get the strain tensor. Then the principal values
            %   is calculated and returned. The values returned need to be
            %   divided by the total volume of the simulation!
            %
            %   Strain tensor of a edge pq is the multiplication of two
            %   terms. The first one is a COLUMN vector (3x1) formed by the
            %   diference between the displacement of the grain q and p.
            %   The second one is a LINE vector (1x3) formed by the
            %   complement area vector of the same edge. DxN*NxD => DxD
            %   tensor
            
            %In the case of partial calculation we must get only the values
            %of the grains that are inside the rectangle
            if isempty(PD)
                %total calculation
                t=(sc.Displacements(sc.EdgeID(:,2),:)...
                    -sc.Displacements(sc.EdgeID(:,1),:))'*sc.CompAVector;
            elseif app.SubdivisionButton.Value
                %subdivision
                t=zeros(sc.D,sc.D,PD.SubCol*PD.SubLin);
                for i=1:(PD.SubCol*PD.SubLin)
                    chck=ismember(sc.EdgeID,PD.SubGrains{i});
                    lin=chck(:,1) & chck(:,2); %should external border contacts also be included?
                    t(:,:,i)=(sc.Displacements(sc.EdgeID(lin,2),:)...
                        -sc.Displacements(sc.EdgeID(lin,1),:))'*sc.CompAVector(lin,:);
                end
            else
                %partial calculation -get specific grains
                chck=ismember(sc.EdgeID,PD.GrainsRectangle);
                lin=chck(:,1) & chck(:,2); %should external border contacts also be included?
                t=(sc.Displacements(sc.EdgeID(lin,2),:)...
                    -sc.Displacements(sc.EdgeID(lin,1),:))'*sc.CompAVector(lin,:);
            end
            
            %get symetrical part. To get strain value still need to be
            %divided by the volume of the requested area.
            sc.GStTensor=(t+permute(t,[2,1,3]))/2; %t+t' for any dim matrix
            if ~isempty(sc.DelaunayT)
                fprintf('Global Tensor done\n')
            end
        end
        function sc = pgStrainTensor(sc)
            %PCSTRAINTENSOR calculate per grain incremental strain tensor
            %   This function will use the previously calculated
            %   complementary area vector to get the strain tensor of the
            %   region around each grain. To do so we will calculate the
            %   product of the incremental displacement x the CAV of each
            %   edge connected to the grain. However this value has to be
            %   pondered by the volume of both grains sharing the edge,
            %   thus resulkting i the following formulation.
            %
            %   e^g= sum( d^mg * du^mg /(V^g +V^m)) for all grains m
            %   connected to the grain g.
            
            %start vectors that will contain calculated data
            pgST=zeros(sc.D,sc.D,sc.NbG);	%per cell strain tensor
            %for each Real Grain :
            for i=1:sc.NbG
                %New Method :
                %Get the IDs of the edges containing the
                %grain i
                edgID=find(sum(sc.EdgeID==i,2));
                %To construct a 'per grain' value, each edge
                %d*u value will be pondered by the mass of the grain over
                %the mass of both grains in the edge.
                tjk=(sc.Displacements(sc.EdgeID(edgID,2),:)...
                    -sc.Displacements(sc.EdgeID(edgID,1),:))'...
                    *( sc.CompAVector(edgID,:)./(4/3*pi()*...
                    (sc.Radius(sc.EdgeID(edgID,2)).^3+...
                    sc.Radius(sc.EdgeID(edgID,1)).^3)) );
                %get symetrical part and divide by the volume
                pgST(:,:,i)=(tjk+tjk')/2;
                
                %{
                %Old Method : Get the IDs of the edges containing the
                %grain i
                edgID=find(sum(sc.EdgeID==i,2));
                %Strain tensor of a edge pq is the multiplication of two
                %terms. The first one is a COLUMN vector (3x1) formed by
                %the diference between the displacement of the grain q and
                %p. The second one is a LINE vector (1x3) formed by the
                %complement area vector of the same edge. DxN*NxD => DxD
                %tensor
                tjk=(sc.Displacements(sc.EdgeID(edgID,2),:)...
                    -sc.Displacements(sc.EdgeID(edgID,1),:))'...
                    *sc.CompAVector(edgID,:);
                vol=4/3*pi()*sc.Radius(i)^3;
                %get symetrical part and divide by the volume
                pgST(:,:,i)=(tjk+tjk')/(2*vol);
                %}
            end
            sc.PStTensor=pgST;
            if ~isempty(sc.DelaunayT)
                fprintf('Per Cell Tensor done\n')
            end
        end
        function [sc,gC] = pcStrain(sc,gr,grOld)
            %PCSTRAIN calculates the strain per cell
            % This function will not create singleCells but apply the same
            % principles of that function. Volume is calculated thanks to
            % the determinant of 3 edges, surface vector is calculated from
            % the cross product of 2 vector (direciton later checked), and
            % surface displacement as average dpl of the three grains
            % forming it.
            
            %Get goodcells again
            gC=sum((sc.DelaunayT(:,:))>sc.NbG,2);
            gC=find(gC==0); %all cells that contain no fakepts
            %Calculate displacements if it was not done before
            if nargin==3;sc.Displacements=gr.Coord - grOld.Coord;end
            %Calculate Volume
            sc.CellVol=permute(perCellVolume(sc,gr,gC),[3,2,1]); %per cell volume 1x1xNbGc
            grID=sc.DelaunayT(gC,:); %gc grains
            %create a surface combination matrix
            cmb=nchoosek(1:sc.D+1,sc.D);
            pcStn=zeros(3,3,size(gC,1));
            %for each surface combination do:
            for i=1:size(cmb,1)
                sf=grID(:,cmb(i,:));
                %calculate the average displacement of the grains forming
                %each surface
                u=(sc.Displacements(sf(:,1),:)+...
                    sc.Displacements(sf(:,2),:)+...
                    sc.Displacements(sf(:,3),:))/3; %NbGcx3 matrix
                %calculate the normal of the edge/surface between the
                %grains
                n=0.5*cross(gr.Coord(sf(:,2),:)-gr.Coord(sf(:,1),:),...
                    gr.Coord(sf(:,3),:)-gr.Coord(sf(:,1),:)); %NbGcx3 matrix
                %Make sure the normal points towards the exterior.
                %Compare n with a vector between the grain not in the
                %surface/edge and one in the edge. Check if the angle with
                %the normal is inferior to 90 degrees. If Yes n=-n;
                v=gr.Coord(grID(:,setxor(1:sc.D+1,cmb(i,:))),:)...
                    -gr.Coord(sf(:,1),:); %NbGcx3 matrix
                alpha=acosd(dot(n,v,2)./(vecnorm(n,2,2).*vecnorm(v,2,2)));
                n(alpha<90,:)=-n(alpha<90,:);  %NbGcx3 matrix
                %calculate the contribution of this surface
                pcStn=pcStn+pagemtimes(permute(n,[2,3,1]),permute(u,[3,2,1]));
            end
            sc.CellStn=(pcStn+permute(pcStn,[2,1,3]))/2; %symetric part
            sc.CellGrn=permute(grID,[3,2,1]);   %1x4xNbGc
        end
        %Loops/Cluster creation
        function sc = findOpenEdges(sc,step,app,PD)
            %FINDVCONTACTS find open contacts in the sclist
            %    This function will load the information of the contacts
            %    obtained through LIGGGHTS, it will then compare with the
            %    edges obtained through the Delaunay Triangulation to check
            %    wich of these edges are Open.
            
            %First part : get closed contacts from LIGGGHTS
            ctData = readData('CONTACT',app,step);
            if isempty(ctData);return;end
            ctData=ctData(:,1:2);
            ctData=sortrows(sort(ctData,2));
            if ~isempty(PD)
                %in the partial calculation case, as all the grains are not
                %being used, their ID for the closed contacts obtained through
                %LIGGGHTS is diferent of their ID obtained through the DT. So
                %a 'translation' is needed. To do so, the GrainsRect property
                %will be use. By doing 'find' we get a vector containing the
                %'real ID' of the grain placed in the position that represent
                %the ID on the triangulation
                
                ctData=ctData(ismember(ctData(:,1),PD.GrainsRectangle) &...
                    ismember(ctData(:,2),PD.GrainsRectangle),:);
                %the edges containing only the right grains were found, they
                %must be transformed into the new ones
                M2=zeros(size(ctData));
                for i=1:size(PD.GrainsRectangle,1)
                    f=find(ctData==PD.GrainsRectangle(i));
                    M2(f)=i*ones(size(f));
                end
                ctData=M2;
            end
            edgID=edges(sc.DelaunayT);
            chk=ismember(ctData,edgID,'rows');
            if sum(chk)~=size(ctData,1)
                %Problem if all contacts are not shown in the DT. Add a
                %warning here to be printed in the command pannel
                fprintf('#### PROBLEM - NOT ALL CONTACTS FOUND ON DT ####\n')
                f=find(~chk);
                nf=numel(f);
                fprintf(['Ids %d : ' repmat('%d ',1,nf) '\n' ],[nf;f])
                save( [MakePath(app,'LOOPER') 'step'...
                    num2str(step) '.mat'] ,'sc','-v7.3');
            end
            %Second part : get virtual-real contacts
            %The grains that are in contact with the wall may have 'open'
            %contacts with each others, but they actually should form loops.
            %To be able to take this into account we will check in the DT
            %connectivity list for cells with only 2 (2 or 3 in 3D) real
            %grains inside. We will get then all these edges and add them to
            %the 'ClosedEdge' matrix to be used in the loops calculation.
            cnl=sum((sc.DelaunayT(:,:))<=sc.NbG,2);  %nb of real grains per tetrahedron
            rvc=sort(sc.DelaunayT(cnl==2,:),2); %get the lines with at least 2 reals
            rvc=rvc(:,1:2); %get from these only the id of the 'real'
            if sc.D>2
                rvc2=sort(sc.DelaunayT(cnl==3,:),2); %get the lines with only 3 real and 1 fake
                rvc2=rvc2(:,1:3);
                for i=1:size(rvc2,1) %recombine the real to get the contacts between them
                    rvc=[rvc;nchoosek(sort(rvc2(i,:)),2)];
                end
            end
            rvc=unique(sortrows(sort([rvc;ctData],2)),'rows');
            
            %get all edges containing ONLY real values
            edgID=edgID(edgID(:,1)<=sc.NbG & edgID(:,2)<=sc.NbG,:);
            
            %get a logical vector containing the 0 for closed contacts and 1 for open ones
            vE=~ismember(edgID,rvc,'rows');
            
            %remove grains that are near
            
            %Good cell vector : find the ID of the tetrahedrons that do not
            %contain 'fakepoints' or grains in contact with it. Start by
            %getting the cells that contain "fakepoints" and all the grains
            %that make cells with these
            %bId=sum((sc.DelaunayT(:,:))>sc.NbG,2);
            %bId=unique(sc.DelaunayT(bId>0,:));
            %sc.BadIds=bId;
            %Approach 2 : get all the grains that are fake points
            bId=(sc.NbG+1:numel(unique(sc.DelaunayT(:,:))));
            sc.BadIds=bId;
            %get good cells
            gC=sum((ismember(sc.DelaunayT(:,:),bId)),2); %bad cells
            gC=find(gC==0); %all cells that do not contain those grains
            sc.GoodCells=gC;
            sc.NbGc=size(gC,1);
            
            %save the closed and open contacts.
            sc.ClosedEdges=rvc;
            sc.OpenEdges=edgID(vE,:);
        end
        function [oScon,nbOs,clsS]=findOpenSurf(sc,cID,oE,cE)
            %OPENSURF return the ID of the cells connected to the cID by Open surfaces
            %   This function will use the knwoled of open and closed edges to
            %   determinate the IDs of grains forming open cells for a given
            %   tetrahedron, if there are any.
            %
            %   vScon - ID of the cells connected to cID by open surfaces
            %   nbVs - number of open surfaces (or nb of vScon, same thing)
            %   aniEdg - edges that may be used on anisotropy calculation
            %   clsS - retunr the closed surfaces used on paraview drawing and
            %   anisotropy
            
            %Nb open edges
            nbvE=size(oE,1);
            nbOs=0;
            oScon=double.empty(0,1);
            aniEdg=double.empty(0,2);
            clsS=double.empty(0,sc.D);
            switch nbvE
                case 3
                    %one open surface if 3 vE formed by only 3 grains
                    if (max(size(unique(oE)))~=3);return;end
                    atC=edgeAttachments(sc.DelaunayT,oE);
                    oScon=setxor(intersect(atC{1},atC{2}),cID);
                    aniEdg=cE;
                case 4
                    %one open surface if one of the grains is repeated at least 3x
                    [uVert,~,cnt]=unique(oE);   %unique value and their oder
                    h=histcounts(cnt);
                    if max(h)<3;return;end
                    %Get the neighboring cell ID
                    h=ismember(oE,uVert(h==2)); %get the position of the ones repeated 2 times
                    chk=(h(:,1) | h(:,2));      %get ID of vE containing the ones repeated 2 times
                    atC=edgeAttachments(sc.DelaunayT,oE(chk,:));	%get the cells attached to all these edges
                    oScon=setxor(intersect(atC{1},atC{2}),cID);     %get the intersecting cell
                    aniEdg=[cE;oE(~chk,:)];     %closed edges + the one not used in the surface
                case 5
                    %Two open surfaces and 2 grains are going to repeated 3x
                    [uVert,~,cnt]=unique(oE);   %unique value and their oder
                    h=histcounts(cnt);          %get nb of repeated values
                    h=find(h==3);               %get the postion of the ones repeated 3x
                    %get the cells that are attached to this edge and intersect with
                    %those that are neighboring cId.
                    atC=edgeAttachments(sc.DelaunayT,[uVert(h(1)) uVert(h(2))]); %
                    oScon=intersect(neighbors(sc.DelaunayT,cID),atC{:});
                    clsS=[[cE uVert(h(1))];[cE uVert(h(2))]];
                case 6
                    %All four surfaces are open and so all the cells neighboring the
                    %cID are connected to it and none closed surfaces are present.
                    oScon = neighbors(sc.DelaunayT,cID);
            end
            oScon=oScon(ismember(oScon,sc.GoodCells)); %compare with good cells
            nbOs=size(oScon,2);
            
            if size(aniEdg,1)==3
                %create the external closed surfaces for case 3 and 4 above
                vert=setxor(aniEdg(2,:),aniEdg(3,:));
                clsS=sort([[aniEdg(1,:) vert(1)];
                    [aniEdg(1,:) vert(2)];
                    union(aniEdg(2,:),aniEdg(3,:))],2);
            end
            
        end
        function sc = findLoops(sc)
            %FINDLOOPS This funciton will calculate all the loops formed
            %   This function will use the properties of the space cell
            %   system to find the loops formed by the grains. The
            %   loops can be defined as a void that is surrounded by
            %   grains. On the space cell system two types of contact are
            %   defined, closed and open. To find the loops, all cells
            %   that have a shared open edge will be joined.
            
            %Good cell vector : find the ID of the tetrahedrons that do not
            %contain 'fakepoints'
            gC=sc.GoodCells;
            
            %Create loop storing variables
            lArray(size(gC,1))=singleLoop();lp=1;
            
            %Get all contacts belonging to each of the good cells and check
            %how many of them are Closed or Open
            spCell=sort(sc.DelaunayT(gC,:),2);
            edges=cat(3,[spCell(:,1),spCell(:,2)],[spCell(:,1),spCell(:,3)],...
                [spCell(:,2),spCell(:,3)]);
            chkV=ones(size(gC,1),3);
            for i=1:3
                chkV(:,i)=ismember(edges(:,:,i),sc.OpenEdges,'rows');
            end
            edges=permute(edges,[3,2,1]); %permute for easier acces later
            
            %Separate the cells into 2 types.
            %The first one are direct loop 3, and can already be written.
            %The other types need to check for the open surfaces
            g1=find(sum(chkV,2)==0);     %id of the gC vector
            g2=find(sum(chkV,2)>0);      %id of the gC vector
            
            %save all g1 loop 3
            for i=1:size(g1,1)
                rcID=gC(g1(i)); %actual cell ID
                lps=singleLoop(sc.DelaunayT(rcID,:),rcID,3,...
                    edges(:,:,g1(i)));
                lArray(lp)=lps;         %save LoopObject
                lp=lp+1;
            end
            
            %Prepare the spaceCells for the loopJoining
            sCarray(size(g2,1))=spaceCell();
            for i=1:size(g2,1)
                rcID=gC(g2(i)); %actual cell ID
                %get open edges connected to g2(i)
                vE=edges((chkV(g2(i),:)==1),:,g2(i));
                aniE=edges((chkV(g2(i),:)==0),:,g2(i));
                %check cells connected to them
                atC=edgeAttachments(sc.DelaunayT,vE);
                neighID=setxor([atC{:}],rcID);
                sCarray(i)= spaceCell('LOOPS',neighID,rcID,aniE);
            end
            
            fprintf('Joining Loops\n')
            %Create variables
            usd=zeros(sc.NbC,1);u=1; %contain sigleCells that were already used
            for i=1:size(sCarray,2) %for each spaceCell created
                loopI=sCarray(i).ID;
                %check if cell was already joined to other Loop
                if sum(usd==loopI)==1 ;continue;end
                %get the cells connected to cID through open edges
                vECon=sCarray(i).OpenCon;
                nbVc=size(vECon,2); %nb of grains connected
                %Start loop joining
                if nbVc>0
                    nbVc2=0;
                    %To be sure all cells were added, the number of unique
                    %virutal edges between two interations must stay
                    %constant.
                    while nbVc2~=nbVc
                        nbVc2=nbVc;
                        %get the cells that share the open edges,
                        %there are always two cells sharing a edge
                        ea=[loopI,vECon];
                        %get non repeated values
                        loopI=unique(ea);
                        loopI=loopI(ismember(loopI,gC));
                        %make sure all vertices are loopi on a line : 1xN vectors
                        if size(loopI,1)>size(loopI,2);loopI=loopI';end
                        %calculate open edges for all the 'spaceCell's
                        %composing the loop
                        scLpI=ismember([sCarray.ID],loopI);
                        vECon=unique([sCarray(scLpI).OpenCon]);
                        nbVc=size(vECon,2);
                    end
                end
                %Save found loop into the cell array by creating a
                %singleLoop object.
                lps=singleLoop(sc.D,sc.DelaunayT(loopI,:),loopI,sCarray(scLpI));
                lArray(lp)=lps;         %save LoopObject
                lp=lp+1;
                %Save spaceCells that were used into the 'usd' array, so
                %loop is not repeated
                s=sum(loopI>0);
                usd(u:(u-1+s))=loopI(1:s);
                u=u+s;
            end
            fprintf('Loops finished\n')
            %remove the excess prealocated values of lArray
            lArray=lArray(1:(lp-1));
            sc.Loops=lArray;
        end
        function sc = findCluster(sc)
            %FINDCLUSTER This funciton will calculate all the loops formed
            %   This function will use the properties of the space cell
            %   system to find the clusters formed by the grains. The
            %   clusters can be defined as a void that is surrounded by
            %   grains. On the space cell system two types of contact are
            %   defined, closed and open. These contact define open and
            %   closed surfaces. To find the clusters, all cells that have
            %   a shared open surface will be joined.
            %
            %   Clusters of size 4 are no longer being as individual
            %   singleLoop objects to save in total space. A single object
            %   will be created containign all necessary information and
            %   will be stored in sc.Cl4 property.
            fprintf('Calculating open surfaces at %d\n',sc.Step)
            
            %Good cell vector : find the ID of the tetrahedrons that do not
            %contain 'fakepoints'
            gC=sc.GoodCells;
            
            %Create loop storing variables
            clArray(size(gC,1))=singleLoop();lp=1;
            
            %As there are too many Cluster4, to save speed and memory they
            %will not be stored into loops, they will only be counted
            %inside the following variables :
            clt4=0;                     %nb of Clt4
            clt4ID=zeros(numel(gC),1);  %ID of cell forming Clt4
            clt4z=zeros(numel(gC),1);   %Coordination nb of each Clt4
            %Get all contacts belonging to each of the good cells and check
            %how many of them are Closed or Open
            spCell=sort(sc.DelaunayT(gC,:),2);
            edges=cat(3,[spCell(:,1),spCell(:,2)],[spCell(:,1),spCell(:,3)],...
                [spCell(:,1),spCell(:,4)],[spCell(:,2),spCell(:,3)],...
                [spCell(:,2),spCell(:,4)],[spCell(:,3),spCell(:,4)]); %ngCx2x6
            %check Open edges
            chkO=ones(size(gC,1),6);
            for i=1:6
                chkO(:,i)=ismember(edges(:,:,i),sc.OpenEdges,'rows');
            end
            edges=permute(edges,[3,2,1]); %permute for easier acces later 6x2xngC
            
            %Separate the cells into 2 types. The first one contain cells
            %with less than 2 open edges, thus no open surface. The second
            %group may have open surfaces thus needing additional
            %treatment.
            g1=find(sum(chkO,2)<3);     %id of the gC vector
            g2=find(sum(chkO,2)>=3);    %id of the gC vector
            
            %all g1 are cluster4
            clt4=clt4+size(g1,1);
            clt4ID(1:numel(g1))=gC(g1);
            clt4z(1:numel(g1))=2*(6-sum(chkO(g1,:),2))/4;
            
            %Prepare the spaceCells objects
            scArray(size(g2,1))=spaceCell(); c=1;
            sCI=zeros(size(g2,1),1);
            for i=1:size(g2,1)
                cID=gC(g2(i)); %actual cell ID
                %get open and closed edges connected to g2(i)
                oE=edges((chkO(g2(i),:)==1),:,g2(i));
                cE=edges((chkO(g2(i),:)==0),:,g2(i));
                %check the open surfaces
                [oSCon,nbOs,clsS]=findOpenSurf(sc,cID,oE,cE);
                if nbOs>0
                    %if at least one open surface create a spaceCell object
                    %to prepare for cluster merging
                    scArray(c)= spaceCell("LOOPS",oSCon,cID,cE,clsS);
                    sCI(c)=g2(i); %save the indice of the sCarray corresponding to 'edges'
                    c=c+1;
                else
                    %clt 4 no longer being saved as a singleLoop to save
                    %time and space
                    clt4=clt4+1;
                    clt4ID(clt4)=cID;
                    clt4z(clt4)=2*size(cE,1)/4;
                end
            end
            %remove excess values from cl4 data and save on sc so I can
            %acces inside of singleLoop call
            clt4ID=clt4ID(1:clt4);
            clt4z=clt4z(1:clt4);
            sc.Clt4=singleLoop(4,clt4,clt4ID,clt4z);
            
            %remove the excess prealocated values of sCarray
            scArray=scArray(1:(c-1));
            edges=edges(:,:,sCI(1:(c-1))); %remove the edges of cells belonging to Cl4
            fprintf('Mergin cells into Clusters\n')
            %Create variables
            usd=zeros(sc.NbC,1);u=1; %contain sigleCells that were already used
            for i=1:size(scArray,2) %for each spaceCell created
                cltI=scArray(i).ID;
                %check if cell was already joined to other Cluster
                if sum(usd==cltI)==1 ;continue;end
                %get the cells connected to cID through open surfaces
                oSCon=scArray(i).OpenCon;
                nbOs=size(oSCon,2); %nb of grains connected
                %Start Cluster joining
                if nbOs>0
                    nbOc2=0;
                    %To be sure all cells were added, the number of unique
                    %virutal surfaces between two interations must stay
                    %constant.
                    while nbOc2~=nbOs
                        nbOc2=nbOs;
                        %get non repeated values of the cells that share
                        %the open surfaces, there are always two cells
                        %sharing a surface
                        cltI=unique([cltI,oSCon]);
                        %check goodcells to be sure
                        cltI=cltI(ismember(cltI,gC));
                        %make sure all vertices are clustI on a colum : Nx1 vectors
                        if size(cltI,1)>1;cltI=cltI';end
                        %get all the cells connected to the ones in clustI
                        scOI=ismember([scArray.ID],cltI);
                        oSCon=unique([scArray(scOI).OpenCon]);
                        %update value and redo
                        nbOs=size(oSCon,2);
                    end
                end
                %Save found cluster into the cell array by creating a
                %singleLoop object.
                lps=singleLoop(sc.D,sc,cltI,scArray(scOI),edges(:,:,scOI));
                clArray(lp)=lps;         %save LoopObject
                lp=lp+1;
                
                %Save spaceCells that were used into the 'usd' array, so
                %cluster is not repeated
                s=sum(cltI>0);
                usd(u:(u-1+s))=cltI(1:s);
                u=u+s;
            end
            fprintf('Cluster identification finished\n')
            %remove the excess prealocated values and save properties.
            %Cluster4 values will be stored in a single loop containig all
            %information on number, cell IDs, coordination number. Creating
            %a singleLoop object for each cell would take a very large
            %space.
            clArray=clArray(1:(lp-1));
            sc.Loops=clArray;
            %take into account the cl4 modifications
            %mb=ismember(clt4ID,cat(2,clArray.sCells));
            %sc.Clt4.nbCells=sc.Clt4.nbCells-sum(mb);
            %sc.Clt4.sCells=sc.Clt4.sCells(~mb);
            %sc.Clt4.Z=sc.Clt4.Z(~mb);
        end
        %Clusters analysis
        function lA = clusterAnisotropy(sc,gr)
            %LOOPANISOTROPY calculate the anisotropy of loops
            % This function will calculate, for each singleLoop object
            % obtained through the spaceCellSystem class, a tensor that
            % will characterize the anisotropy of the loop.
            lps=find(cat(1,sc.Loops.Order)>4);
            aT=3; %anisotropy types
            vC=2; %values calculated per category
            lA.Angles=zeros(numel(lps),aT,vC);
            lA.Values=zeros(numel(lps),aT,vC);
            lA.Order=cat(1,sc.Loops(lps).Order);
            r3=gr.Radius.^3;
            for i=1:numel(lps)
                sl=sc.Loops(lps(i));
                %Surface Anisotropy : find normal to each surf through
                %cross product and the value obtained is a vector with the
                %area of the surface and normal to it. The sense of the
                %vector (+ or -) does not influence the value because we
                %calculate n^2 later. Divide one of n by the its norm so
                %the tensor we obtain is calculated as area*normal*normal.
                %Later it is divided by its trace to get a unitaire tensor.
                
                surf=sl.Vertices; %id of each surface Nx3 vector
                n=cross(gr.Coord(surf(:,1),:)-gr.Coord(surf(:,2),:),...
                    gr.Coord(surf(:,1),:)-gr.Coord(surf(:,3),:));
                t=n'*(n/norm(n))/size(n,1); %area*normal*normal
                [Ve,Di]=eig(t/trace(t));    %eig of unitaire tensor
                Di=diag(Di);
                pos=find(abs(Di)==min(abs(Di),[],'all')); %find column smallest eig val
                ele(1)=atand(Ve(3,pos)/sqrt(Ve(1,pos)^2+Ve(2,pos)^2)); %w/sqrt(x^2+y^2)
                azi(1)=atand(Ve(2,pos)/Ve(1,pos));       %y/x
                dev(1)=sqrt((Di(1)-Di(2)).^2+... 	%(x-y)^2
                    (Di(2)-Di(3)).^2+...            %(y-z)^2
                    (Di(3)-Di(1)).^2);              %(z-x)^2
                eiv(1)=Di(pos);
                
                %Center of Gravity : get the center of mass G of the loop.
                %Since all grains have the same density it is in function
                %of r^3.
                r3g=r3(sl.Grains);
                clGr=gr.Coord(sl.Grains,:);
                G=sum(clGr.*r3g,1)/sum(r3g);
                n=clGr-G; %vector Gr center to cluster center
                t=n'*n; %tensor calculation
                [Ve,Di]=eig(t/trace(t)); %transforming into a unitaire tensor
                Di=diag(Di);
                pos=find(abs(Di)==max(abs(Di),[],'all')); %find column largest eig val
                ele(2)=atand(Ve(3,pos)/sqrt(Ve(1,pos)^2+Ve(2,pos)^2)); %w/sqrt(x^2+y^2)
                azi(2)=atand(Ve(2,pos)/Ve(1,pos));       %y/x
                dev(2)=sqrt((Di(1)-Di(2)).^2+... 	%(x-y)^2
                    (Di(2)-Di(3)).^2+...            %(y-z)^2
                    (Di(3)-Di(1)).^2);              %(z-x)^2
                eiv(2)=Di(pos);
                
                %Fabric Tensor : get the closed edges in each cluster and
                %calculate the fabric tensor.
                cd=gr.Coord(sl.ClEdges(:,1),:)-gr.Coord(sl.ClEdges(:,2),:);
                cd=cd./vecnorm(cd,2,2); %unitaire transf
                t=cd*cd'/size(cd,1);
                [Ve,Di]=eig(t/trace(t));    %eig of unitaire tensor
                Di=diag(Di);
                pos=find(abs(Di)==max(abs(Di),[],'all')); %find column largest eig val
                ele(3)=atand(Ve(3,pos)/sqrt(Ve(1,pos)^2+Ve(2,pos)^2)); %w/sqrt(x^2+y^2)
                azi(3)=atand(Ve(2,pos)/Ve(1,pos));       %y/x
                dev(3)=sqrt((Di(1)-Di(2)).^2+... 	%(x-y)^2
                    (Di(2)-Di(3)).^2+...            %(y-z)^2
                    (Di(3)-Di(1)).^2);              %(z-x)^2
                eiv(3)=Di(pos);
                
                %angles to principal values %Nloops,ManiTypes,Kvalues
                lA.Angles(i,:,1)=ele;     %elevation angle
                lA.Angles(i,:,2)=azi;     %azimuth angle
                lA.Values(i,:,1)=dev;	  %deviatoric
                lA.Values(i,:,2)=eiv;     %eig
            end
            
        end
        function sc = clusterStress(sc,step,app,gr)
            %LOOPSTRESS Calculates the stress for each loop
            %First part : get contacts information from LIGGGHTS
            fprintf('Begin Cluster Stress \n')
            
            rData = readData('CONTACT',app,step);
            if isempty(rData);return;end
            rData=rData(:,1:5);%Nx5x1 matrix
            
            %get only the contac Ids for later
            ctcID=rData(:,1:2);
            %Get the branch vectors of these contacts
            n=gr.Coord(ctcID(:,2),:)-gr.Coord(ctcID(:,1),:); %Nx3x1 matrix
            n=n./vecnorm(n,2,2);    %calculate the normal
            n=permute(n,[3,2,1]); %1x3xN matrix
            %do the same permute with forces
            f=permute(rData(:,3:5),[3,2,1]); %1x3xN matrix
            %lxf calculation
            lf=(pagemtimes(permute(n,[2,1,3]),f)+...
                pagemtimes(permute(f,[2,1,3]),n))/2; %3x3xN matrix
            %2x calculation to get the symetrical tensor
            
            %Get goodcells
            %Get goodcells again
            gC=goodCell(sc);
            
            %calculate the volume of all cells by the determinant of a
            %matrix (check spaceCell object for more information). To do
            %all in one go the calculation will be vectorised.
            volMat=perCellVolume(sc,gr,gC);
            
            %Now, for each cluster, get the lf value of the contacts with
            %exterior grains and divide it by 2x the cluster volume. The
            %same is needed for Cl4 category
            %lp4ID=gC(~ismember(gC,cat(2,sc.Loops.sCells)')); %id Cl4 cell location
            lp4ID=sc.Clt4.sCells;
            lp=numel(sc.Loops);
            lpNb=lp+numel(lp4ID);
            vol=zeros(lpNb,1);
            strT=zeros(3,3,lpNb);
            for i=1:lpNb
                if i>lp
                    clCl=lp4ID(i-lp);
                    clGr=sc.DelaunayT(lp4ID(i-lp),:);
                else
                    %Larger Clusters
                    clCl=sc.Loops(i).sCells;
                    clGr=unique(sc.Loops(i).Vertices);
                end
                %need to transform DelaunayT id into Goodcell ID for next
                %line
                V=sum(volMat(ismember(gC,clCl)));
                %check contacts done with the exterior. So Only one of the
                %grains in the loop can be in the contact list.
                chk=(ismember(ctcID(:,1),clGr)...
                    +ismember(ctcID(:,2),clGr))==1;
                cts=ctcID(chk,:);
                [~,pos]=ismember(cts,clGr);
                r=permute(gr.Radius(sum(pos,2)),[3,2,1]);
                %calculate the tensor
                strT(:,:,i)=sum(lf(:,:,chk).*r,3)/V;
                vol(i)=V;
            end
            
            %Deviatoric and Principal calculation
            q=sqrt( ((strT(1,1,:)-strT(2,2,:)).^2+...	%(sig1-sig2)^2 /2
                (strT(2,2,:)-strT(3,3,:)).^2+...        %(sig2-sig3)^2 /2
                (strT(3,3,:)-strT(1,1,:)).^2+...        %(sig3-sig1)^2 /2
                6*strT(1,2,:).^2 +...                   %3*(sig12)^2
                6*strT(1,3,:).^2 +...                   %3*(sig13)^2
                6*strT(2,3,:).^2) /2);                  %3*(sig23)^2
            p=(strT(1,1,:)+strT(2,2,:)+strT(3,3,:))/3;
            str=permute([p,q],[3,2,1]);
            %Create a matrix to return the results in the following form:
            [sc.Loops.Stress]=disperse(str(1:lp,:),2);
            [sc.Loops.Volume]=disperse(vol(1:lp,:),1);
            sc.Clt4.Stress=str(lp+1:end,:);
            sc.Clt4.Volume=vol(lp+1:end,:);
            fprintf('Cluster Stress finished \n')
        end
        function sc = clusterStrain(sc,gr)
            %LOOPSTRAIN calculates the strain of each cluster
            %Get per cell strain
            fprintf('Begin Cluster Strain \n')
            [sc,gC] = pcStrain(sc,gr);
            %Add values to Loops
            for l=1:numel(sc.Loops)
                [~,p]=ismember(sc.Loops(l).sCells,gC);
                sc.Loops(l) = loopStrain(sc.Loops(l),...
                    sc.CellStn(:,:,p),...
                    sc.CellVol(:,:,p));
            end
            fprintf('Cluster Strain finished \n')
        end
        function sc = clusterVolume(sc,gr)
            %LOOPSTRAIN calculates the strain of each cluster
            %Get per cell strain
            fprintf('Begin Cluster Volume \n')
            if isempty(sc.GoodCells)
                sc.GoodCells=goodCell(sc);
            end
            vl = perCellVolume(sc,gr,sc.GoodCells);
            %Add values to Loops
            for l=1:numel(sc.Loops)
                [~,p]=ismember(sc.Loops(l).sCells,sc.GoodCells);
                sc.Loops(l).Volume=sum(vl(p));
            end
            fprintf('Cluster Volume finished \n')
        end
        function cT = clusterTransformation(sc,clN,gr,clO)
            %CLUSTERTRANSFORMATION Count the transformation of cluster
            %orders
            % This fucntion will take as imputs two diferent
            % spaceCellSystems objects from two different timesteps and
            % analyse what the clusters of the previous step became in the
            % new one.
            %
            % PS : By calling this function with Old and New reversed the
            % analysis is then made in reverse. We obtain wich clusters in
            % the new step were created by the one in the old.
            
            fprintf('Begin Cluster Transf \n')
            %Compare clusters scOld=>scNew to see which ones were
            %destroyed. Get the members of the old that are not present on
            %the new to check what happened to the grains in the New
            difCl4=find(~ismember(clO.cl4Id,clN.cl4Id,'rows'));
            difCl6=find(~ismember(clO.cl6Id,clN.cl6Id,'rows'));
            %in the case of cl8+ the matrixes may have different sizes. Nan
            %values will be added to the larger one to allow a direct
            %comparison of clusters.
            if size(clO.cl8Id,2)~=size(clN.cl8Id,2)
                mx=max(size(clO.cl8Id,2),size(clN.cl8Id,2));
                clOld8Id=[clO.cl8Id,...
                    NaN(size(clO.cl8Id,1),mx-size(clO.cl8Id,2))];
                cl8Id=[clN.cl8Id,...
                    NaN(size(clN.cl8Id,1),mx-size(clN.cl8Id,2))];
                difCl8=find(~ismember(clOld8Id,cl8Id,'rows'));
            else
                difCl8=find(~ismember(clO.cl8Id,clN.cl8Id,'rows'));
            end
            %Find out, for each grain, which cluster it is part of in the
            %scNew. First identify for each grain wich cluster it makes
            %part of in the scNew (through findClusterGrain function) then
            %look to the old cluster on the scOld and check what the grains
            %have become
            [~,clCell]= findClusterGrain('ClustTransf',sc,1:gr.Nb);
            
            %Prepare variables
            resCl4=zeros(2*numel(difCl4),3);k4=1;
            resCl6=zeros(2*numel(difCl6),3);k6=1;n6=numel(difCl6);
            resCl8=zeros(5*numel(difCl8),4);k8=1;n8=numel(difCl8);
            
            %Start checking for Cl4
            for i=1:numel(difCl4)
                %get the grains forming the Old cell
                gri=clO.cl4Id(difCl4(i),:);
                %check if the a cell containing all of the 4 still exists
                %on the scNew
                [chk,ncID]=ismember(gri,clN.srtDT,'rows');
                if chk
                    %if it was not broken between timesteps check to which
                    %cluster it belongs to. [NOrder,OldID, NewID]
                    resCl4(k4,:)=[clCell(clCell(:,3)==ncID,2),difCl4(i),ncID];
                    k4=k4+1;
                else
                    %If the cell was broken try to find to which clusters
                    %the grains now belong to.
                    avePt=sum(gr.Coord(gri,:).*gr.Radius(gri).^3,1)/...
                        sum(gr.Radius(gri).^3);
                    ncID=pointLocation(sc.DelaunayT,avePt);
                    clt=clCell(clCell(:,3)==ncID,1:2);
                    %if it is empty it means the cell it belongs to
                    %does not make part of the sc.Loops cells, so it is
                    %a different cell of order 4
                    if isempty(clt);continue;end
                    resCl4(k4,:)=[clt(:,2),difCl4(i),ncID];
                    k4=k4+1;
                end
            end
            
            %Start checking for Cl6 and Cl8+. The procedure is the same but
            %the data pool changes between both. So two 'if' will separate
            %wheter we analyse one or the other.
            for i=1:(n6+n8)
                %get the cells forming the cl
                if i>n6
                    cel=clO.cl8cells{difCl8(i-n6)};
                else
                    cel=clO.cl6cells{difCl6(i)};
                end
                if isrow(cel);cel = cel';end %force column
                %Check for surviving cells between both steps by comparing
                %grains with the newDT. chk is a logical, ncID is the
                %position of the cell in the newDT (new cell ID)
                [chk,ncID]=ismember(clO.srtDT(cel,:),clN.srtDT,'rows');
                
                %Part 1 : for those that still exists get the clusters
                %it belong to. Check on the clCell vector to which cluster
                %the cells ncID belong to.
                unbCl=double.empty(0,3);
                if sum(chk)>0 %at least 1 unbroken cell
                    ncID=ncID(chk);
                    c=cel(chk);
                    [ckClNew,posClN]=ismember(ncID,clCell(:,3));
                    unbCl=[clCell(posClN(posClN>0),2),...
                        c(ckClNew),...
                        ncID(ckClNew)]; %[NOrder,OldID, NewID] vector
                    
                    %if any of the surviving cells did not find a match in
                    %clCell, it means it is a cluster 4. Thus it need to be
                    %added to the unbCl vector
                    nb=numel(ncID)-sum(ckClNew);
                    if sum(~ckClNew)>0
                        unbCl=[unbCl;
                            [ones(nb,1)*4,c(~ckClNew),ncID(~ckClNew)]]; %#ok<*AGROW>;
                    end
                end
                
                %Part 2 : cells that did not survive
                if sum(~chk)~=0
                    %If any of the cells was broken try to find to which
                    %clusters the grains now belong to.
                    f=find(~chk);
                    brkCl=zeros(numel(f),3);
                    for j=1:numel(f)
                        %get the grains that belong to the cell that was
                        %broken
                        gri=clO.srtDT(cel(f(j)),:);
                        %get coordinates of the mean point between the
                        %grains forming the broken cell
                        avePt=sum(gr.Coord(gri,:).*gr.Radius(gri).^3,1)/...
                            sum(gr.Radius(gri).^3);
                        idC=pointLocation(sc.DelaunayT,avePt);
                        %Identify the order of the cluster the cell makes
                        %part of
                        clt=clCell(clCell(:,3)==idC,2); %[Order]
                        %if no value was found it is a cl4
                        if isempty(clt);clt=4;end
                        %save [NOrd,OId,NId] and go next
                        brkCl(j,:)=[clt,cel(f(j)),idC];
                    end
                    %Remove the 'unique' thus transforming the calculation
                    %in a event value, per cell transformation. Thus 2 Cl4
                    %becoming a Cl6 is the same value of a Cl6 becoming 2
                    %Cl4. %unbCl=unique([unbCl;brkCl]);%
                    unbCl=[unbCl;brkCl];
                end
                n=size(unbCl,1);
                %save into calculation values
                if i>n6
                    %for Cl8+ add the order of the cluster it is coming
                    %from. Transform into [OOrder,NOrder,OldID, NewID]
                    resCl8(k8:(k8+n-1),:)=...
                        [ones(n,1)*clO.cl8Order(difCl8(i-n6)),unbCl];
                    k8=k8+n;
                else
                    resCl6(k6:(k6+n-1),:)=unbCl; %[NOrder,OldID, NewID]
                    k6=k6+n;
                end
                %Save the transformation value into the cluster object
                
            end
            
            %Calculate the number of cells in each category
            nbC=[size(clO.cl4Id,1),...
                numel(cat(2,clO.cl6cells{:})),...
                numel(cat(2,clO.cl8cells{clO.cl8Order<22})),...
                numel(cat(2,clO.cl8cells{clO.cl8Order>20}))];
            %Save results in the properties of the cltrf object that will
            %be returned.
            cT.Cl4=resCl4(1:(k4-1),:);
            cT.Cl6=resCl6(1:(k6-1),:);
            cT.Cl8=resCl8(1:(k8-1),:);
            cT.NbC=nbC;
            fprintf('Cluster Transf finished \n')
        end
        function sc = clusterVR(sc,gr)
            %CLUSTERVOIDRATIO calculate each cluster Void Ratio
            % This function go through all clusters inside the simulation
            % and determinate the volume of the grains inside each cell.
            % Then for each loop the Void Ratio will be calculated using
            % the total volume of the cells and the solid volume previously
            % calculated.
            %
            % Departing from the base sphere element, the points located
            % inside the analyzed cluster will be identified. As points
            % created by basesphere are well distributed in the sphere
            % surface, each point will represent a share of the volume of
            % the grain. Thus the solid volume is calculated by accounting
            % the number of points inside.
            fprintf('Begin Cluster VR \n')
            %check total volume calculation
            gC = goodCell(sc);
            vt = perCellVolume(sc,gr,gC);
            %For each good cell, calculate the solid volume inside of it
            [~,vs,~]=perCellVoidRatio(sc,gr,gC,vt);
            
            %Calculate the void ratio of each cluster
            for cl=1:numel(sc.Loops)
                lgC=ismember(gC,sc.Loops(cl).sCells);
                lvs=sum(vs(lgC));
                %if total volume is empty calculate it
                if isempty(sc.Loops(cl).Volume)
                    sc.Loops(cl).Volume=sum(vt(lgC));
                end
                %calculate VR and save it
                sc.Loops(cl).VoidRatio=(sc.Loops(cl).Volume-lvs)/lvs;
            end
            %Identify cluster4 cells
            lgC=~ismember(gC,cat(2,sc.Loops.sCells));
            sc.Clt4.VoidRatio=(vt(lgC)-vs(lgC))./vs(lgC);
            %calculate total vr
            sc.TotalVR=(sum(vt)-sum(vs))/sum(vs);
            fprintf('Cluster VR finished\n')
        end
        %Support functions
        function vl = perCellVolume(sc,gr,gC)
            %PERCELLVOLUME calculates the volume of each gC cell
            % Calculate the volume of cells by the determinant of a matrix
            % formed by 3 vectors departin from the same point. To do all
            % in one go the calculation will be vectorised.
            
            %If gc was not given as argument, take all cells that contain
            %only real grains (no fakepts)
            if nargin==2 || isempty(gC)
                gC = goodCell(sc);
            end
            gcGrID=sc.DelaunayT(gC,:);
            cord=permute(gr.Coord,[3,2,1]); %1x3xN
            mat=[cord(1,:,gcGrID(:,1))-cord(1,:,gcGrID(:,4));...
                cord(1,:,gcGrID(:,2))-cord(1,:,gcGrID(:,4));...
                cord(1,:,gcGrID(:,3))-cord(1,:,gcGrID(:,4));];
            mat=reshape(mat,[],size(mat,3)); %vectorize the matrix
            vl=abs(1/6*(mat(1,:).*mat(5,:).*mat(9,:)...
                +mat(2,:).*mat(6,:).*mat(7,:)...
                +mat(3,:).*mat(4,:).*mat(8,:)...
                -mat(3,:).*mat(5,:).*mat(7,:)...
                -mat(6,:).*mat(8,:).*mat(1,:)...
                -mat(2,:).*mat(4,:).*mat(9,:))'); %Nx1 vector containin det
        end
        function [vr,vs,vv]=totalVoidRatio(sc,gr)
            %TOTALVOIDRATIO calc the void ratio of the goodCell domain
            % This function will aproximate the void ratio inside the the
            % totality of goodCells, thus reducing the influence of the
            % external walls.
            
            %The total volume can be obtained as the sum of goodCell
            %volumes
            gC = goodCell(sc);
            vl = perCellVolume(sc,gr,gC);
            
            %The solid volume will be obtained as the sum of the volume of
            %the spheres partially inside + those completely inside. The
            %first group are those that share a cell with any fakeGrains,
            %the sencd is the opposite
            bC= find(~ismember(1:sc.NbC,gC));      %bad cells
            g=unique(sc.DelaunayT(bC,:));
            g=g(ismember(g,1:sc.NbG));      %grain IDs in the periphery
            
            %Grains completely inside
            vsIN=sum(4/3*pi()*gr.Radius(~ismember(1:sc.NbG,g)).^3);
            %Grains partially inside volume is divede by 2
            vsOUT=sum(4/3*pi()*gr.Radius(ismember(1:sc.NbG,g)).^3)/2;
            
            %Final calculations
            vs=vsIN+vsOUT;      %solid volume
            vv=sum(vl)-vs;      %void volume
            vr=vv/vs;           %void ratio
        end
        function [vr,vs,vv]=perCellVoidRatio(sc,gr,gC,vl)
            %PERCELLVOIDRATIO calculate the void ratio of each goodcell
            % This function will use the
            if nargin==2
                gC = goodCell(sc);
                vl = perCellVolume(sc,gr,gC);
            end
            %For each good cell, calculate the solid volume inside of it
            gcgrID=sc.DelaunayT(gC,:);
            vs=zeros(numel(gC),1);
            for i=1:4
                ids=setdiff(1:4,i);
                r=gr.Radius(gcgrID(:,i));
                R1=gr.Coord(gcgrID(:,ids(1)),:)-gr.Coord(gcgrID(:,i),:);
                R1=R1./vecnorm(R1,2,2).*r;
                R2=gr.Coord(gcgrID(:,ids(2)),:)-gr.Coord(gcgrID(:,i),:);
                R2=R2./vecnorm(R2,2,2).*r;
                R3=gr.Coord(gcgrID(:,ids(3)),:)-gr.Coord(gcgrID(:,i),:);
                R3=R3./vecnorm(R3,2,2).*r;
                
                multi= abs(dot(R1,cross(R2,R3,2),2))./...
                    (r.*r.*r+dot(R1,R2,2).*r...
                    +dot(R1,R3,2).*r+dot(R2,R3,2).*r);
                sf=2*atan(multi);
                vs=vs+(gr.Radius(gcgrID(:,i)).^3).*sf/3;
            end
            %Calculate void volume and void ratio
            vv=vl-vs;
            vr=vv./vs;
        end
        function sc = grainVolume(sc)
            %GRAINVOLUME calculate each grains voronoi volume
            % Normally this funciton is already calculated inside
            % pgStrainTensor. But in the specific case that only
            % gStrainTensor was called, to save disk space by deleting the
            % sCells and DelaunayT from this class we will execute the
            % volume separately. This way if in the future we will be able
            % to calculate pgStrainTensor if asked.
            gV=zeros(sc.NbG,1);
            for i=1:sc.NbG
                %Volume : calculate the boundary of the 'Center' of the
                %tetrahedrons to which the grain belongs. To do so the
                %property 'center' of the spaceCell will be used. After
                %getting the property the correct lines must be chosen.
                %They are the ones where Ids are 0 0 or one of them equal
                %to i.
                cel=cell2mat(vertexAttachments(sc.DelaunayT,i)); %ids of tetrahedron touching grain i
                cel=cel(ismember(cel,sc.GoodCells));
                if isempty(cel)
                    fprintf("Gr %d has no cells formed by only real grains\n",i);
                    gV(i)=NaN;
                    continue;
                end
                cc=[sc.sCells(cel).Center]'; %get property
                cc=cc(cc(:,1)==0 | cc(:,1)==i | cc(:,2)==i, 3:end); %check correct positions not the ids
                [~,vol]=boundary(cc); %calculate volume
                gV(i)=vol;
            end
            sc.GrainVolume=gV;
        end
        function gC = goodCell(sc)
            %GOODCELL will create the goodcell vector that is needed in many funcitons
            % This function was created because in many condition the goodcell vector
            % is purged from the spaceCell file to save space upon save. Thus many
            % functions need to recall this quantity. This function will then normalize
            % the way the goodCell is created.
            gC=sum((sc.DelaunayT(:,:))>sc.NbG,2);
            gC=find(gC==0);
        end
        %purge
        function sc = purge(sc,mode)
            %PURGE remove exces properties to save space and decrease
            %executing time
            switch upper(mode)
                case 'STRAIN'
                    %if isempty(sc.GrainVolume);sc=grainVolume(sc);end
                    sc.ClosedEdges= [];
                    sc.GoodCells=[];
                    sc.Loops= [];
                    sc.OpenEdges= [];
                    sc.sCells= [];
                case {'LOOP','LOOPS'}
                    sc.CellGrn=[];
                    sc.CellStn=[];
                    sc.CellVol=[];
                    sc.ClosedEdges= [];
                    sc.Displacements=[];
                    sc.EdgeID= [];
                    sc.GoodCells=[];
                    sc.GStTensor= [];
                    sc.GrainVolume=[];
                    sc.OpenEdges= [];
                    sc.PStTensor= [];
                    sc.Radius= [];
                    sc.sCells= [];
            end
        end
    end%meth end
end%class end
