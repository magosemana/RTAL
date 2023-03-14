classdef spaceCell
    %SPACECELL single object used to form a space cell system
    %   This class will create a tetrahedron or triangle object containing
    %   the information necessary to do the strain tensor calculation
    %   following the proposition of BAGI 1996.
    
    properties
        %Global properties
        Volume      %Volume (3D) or area (2D) of the cell - double
        
        %Strain related properties
        Center      %Approximation of the center of the cell - 1xD vector.
        SurfaceV    %Surface vector - DxD (vertId x y z) or (vertId y z)
        Strain      %Strain through linear interpolation, not divided by V
        vIDs        %Vertices Ids
        
        %Loop Related properties
        CEdges      %Cell's closed edges - Nx2 vector
        ClosedS     %ID of grains forming the closed surfaces 
        ID          %Line of the delaunayTriangulation this cell represents - integer
        NbOCon      %Number of virtual surface - scalar
        OpenCon     %ID of the cells connected by a open surface to this one
        SolidVol    %Volume of solid inside the cell - double

    end
    
    methods
        %Calling Function
        function sC = spaceCell(mode,varargin)
            %SPACECELL Construct an instance of this class
            %   Creates the representant of a cell (triangle or
            %   tetrahedron) with 4 caracteristcs calculated with the hel
            %   of support functions.
            if ~nargin %empty cell creation
                sC.Center=zeros(1,2);
                sC.ID=0;
                sC.OpenCon=double.empty(0,3);
                sC.SolidVol=0;
                return
            end
            switch upper(mode)
                case "STRAIN" %Args : sc.D,vID,vPos,vRad,Dplc
                    %create surface vectors
                    sC.SurfaceV=surfV(varargin{1},varargin{2},varargin{3});
                    %Following part only is needed on 'PerCell' calculation
                    sC.Volume = cellVolume(varargin{1},varargin{3});
                    sC.Center = cellCenter(varargin{1},varargin{2},...
                        varargin{3},varargin{4});
                    if nargin>5
                        sC=cellStrain(sC,varargin{1},varargin{2},...
                            varargin{5});
                        sC.vIDs=varargin{2};
                    end
                case "LOOPS" %Args : OpCon,IDs,CEdges,CSurf
                    sC.OpenCon=varargin{1};
                    sC.NbOCon=numel(varargin{1});
                    sC.ID=varargin{2};
                    sC.CEdges=varargin{3};
                    if nargin> 4
                        sC.ClosedS=varargin{4};
                    end
                case "CELLSTRAIN"
                    sC.SurfaceV=surfV(varargin{1},varargin{2},varargin{3});
                    sC.Volume = cellVolume(varargin{1},varargin{3});
                    sC=cellStrain(sC,varargin{1},varargin{2},...
                        varargin{5});
                    sC.vIDs=varargin{2};
                case "VOLUME" %Args : sc,vPos,sVol
                    sC.Volume = cellVolume(varargin{1},varargin{2});
                    sC.SolidVol=varargin{3};
            end
        end
    end 
end

%Suport functions
function sfV = surfV(D,vID,vPos)
%CELLAN calculate the surface vector of the cell
%   Surface vector is the outward normal vector to each surface. Its norm
%   is the area of the surface, its ID is the grain point that does not
%   belong to it.
%
%   - 3D Tetrahedrons
%   This function uses the properties of the cross product to calculate the
%   area of the triangle formed by the points. With three point we can
%   obtain two vectors. The area of the triangle formed by these points is
%   then half of the norm of vector obtained through the cross product of
%   the two formed vectors. 
%
%   - 2D Triangles
%   For this version the surface vector coresponds to the outward normal of
%   the EDGE of the triangle and its norm is the LENGTH.
sfV=zeros(D);
if D==3
    %3D Tetrahedrons
    for i=1:4
        f=setxor(1:4,i);
        n=0.5*cross(vPos(f(2),:)-vPos(f(1),:),...
            vPos(f(3),:)-vPos(f(1),:));
        %check if the found normal is pointing outside
        %the volume, if it isn't multiply for -1. First
        %we need to get the fourth grain, the one not
        %used on the normal calculation.
        v=vPos(i,:)-vPos(f(1),:);
        alfa=acosd(dot(n,v)/(norm(n)*norm(v)));
        if alfa<90;n=-n;end
        sfV(i,:)=n;
    end
else
    %2D Triangles
    for i=1:3
        f=setxor(1:3,i);
        %get normal vector bi = (Yj-Yk;Xk-Yj)
        n=[(vPos(f(1),2)-vPos(f(2),2)),...
            (vPos(f(2),1)-vPos(f(1),1))];
        %get the right norm, it must be ||bi||=length(jk)
        n=n*norm(vPos(f(1),:)-vPos(f(2),:))/norm(n);
        %check if the found normal is pointing outside the
        %surface, if it isn't multiply for -1.
        v=vPos(i,:)-vPos(f(1),:);
        alfa=acosd(dot(n,v)/(norm(n)*norm(v)));
        if alfa<90;n=-n;end
        sfV(i,:)=n;
    end
end
sfV=[vID' sfV];
end
function v = cellVolume(D,vPos)
%CELLVOLUME calculates the volume of the tetrahedron (or triangle if 2D)
%   - 3D
%   It is a sixth of the determinant of the vectors formed by the
%   substraction of 3 vertices by the last one.
%   - 2D 
%   On 2d we will use the same calculation made for the surface vector in
%   3D

if D==3
    v = 1/6*abs(det([vPos(1,:)-vPos(4,:);vPos(2,:)-vPos(4,:);...
        vPos(3,:)-vPos(4,:)]));
else
    v=1/2*norm(cross([0 (vPos(2,:)-vPos(1,:))],[0 (vPos(3,:)-vPos(1,:))]));
end
end
function c = cellCenter(D,vID,vPos,vRad)
%CELLCENTER calculates the center of the cell and returns it's value and
%the midpoints of each edge
%   This function will calculate the position of the 'center' of the cell.
%   It will take into account the radius of the grains to determine the
%   mid-point,wheighted mean value of each edge. These mid-points will then
%   be averaged to get the center of the cell. All these values are used in
%   the approximation of the voronoi cell of each grain. 

%   The result is a (2+D)xN vector. The lines contain Id1 Id2 (Px) Py Pz.
%   For the values of 'centerpoints' between the edges. The Ids  of the
%   first column will be 0 0 as it represent the middle point of the entire
%   tetrahedron.

%   The 'fakepoints' represent the walls of the experiment, so when a
%   'fakepoint' is present on the cell, it will be taken into accont to the
%   calculation of the voronoi volume.
nb0=(vRad==0); %check the number of 'fakegrains'
c=double.empty(2+D,0);
%if there is at least 2 real grains in 
if sum(~nb0)>1
    vP2=vPos(~nb0,:);   %get the real grains position
    vRad=vRad(~nb0,:);  %get the real grains radius
    e=nchoosek(1:size(vP2,1),2); %recombine their IDs
    %calculate the midpoints of the edges belonging to the realgrains
    p=(vRad(e(:,1)).*vP2(e(:,1),:)+vRad(e(:,2)).*vP2(e(:,2),:))./...
        (vRad(e(:,1))+vRad(e(:,2)));
    %calculte the center of the cell
    c=[zeros(1,2) mean(p,1)]';
    %add the 'midpoints' p to the returning vector
    c=[c [vID(e) p]'];
    
end
if sum(nb0)>0 %add fakepoints as an acceptable limit
	c=[c [zeros(sum(nb0),2) (vPos(nb0,:))]'];
end

end
function sC=cellStrain(sC,D,vID,grDpl)
%CELLSTRAN Calculates strain of the cell through gauss+linear interpolation
%Strain is calculated as E_ij=1/V*sum(ave(dlp)_i*SurfV_j); the sum is over
%all surfaces, ave(dpl) is the average of the displacement of all three
%grains and SurfV is the surface vector of the surface.

%Combination of Ids
c = nchoosek(1:(D+1),D); %combination
p = perms(1:(D+1));
[~,pos]=ismember(c,p(:,1:end-1),'rows');
nc=vID(p(pos,end)); %opposite id of the combination

%Average displacements of surface
ave=grDpl(permute(c,[3,2,1]),:);    %get the values
ave=permute(reshape(ave',3,3,4),[2,1,3]); %turn into a 3D matrix
ave=mean(ave,1);%average into first dimesion obtaining a 1xDxD+1 matrix

%get the surface vector in the good order
[~,pos]=ismember(sC.SurfaceV(:,1),nc);
sV=permute(sC.SurfaceV(pos,2:end),[2,3,1]); %permute into a Dx1xD+1 matrix

%multiply both and sum into the third dimension to obtain the tensor
t=sum(pagemtimes(sV,ave),3);

%Symetric part
t=(t+t')/2;

%divide by the volume and sum
sC.Strain=t;
end

