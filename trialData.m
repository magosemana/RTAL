classdef trialData
    %TRIALEF object that contains the information obout the test
    %   The object created through this class will contain all the
    %   information needed about the test : the geometry and the pistons
    %   forces and displacements.
    
    properties
        %vectorial values
        Step
        Dz      %vertical displacement Z, sense -z
        Fz      %vertical force Z, sense -z
        Dy1     %horizontal displacement Y, sense +y
        Fy1     %horizontal force Y1, sense +
        Dy2     %horizontal displacement Y, sense -y
        Fy2     %horizontal force Y2, sense -y
        Dx1     %horizontal displacement X, sense +x
        Fx1     %horizontal force X1, sense +x
        Dx2     %horizontal displacement X, sense -x
        Fx2     %horizontal force X2, sense -x
        dE      %incremental strain in lagrangian framework
        InfPts  %Inflection points of Ev and q
        
        %scalar valors
        boxX    %box width  - x
        boxY 	%box length - y
        boxZ 	%box height - z - also pistonY height
        
        %boolean
        checkPiston
    end
    
    methods
        %Calling Functions
        function TD = trialData(app,trialData)
            %check witch folder to load and filename
            if numel(app.DataFolder)==1
                rData = readData('PISTONS',app,app.ConsoStep);
            else
                rData = readData('PISTONS',app,app.ConsoStep);
                rData2 = readData('PISTONS',app,app.ConsoStep+1);
                rData2(:,1)=rData2(:,1)+rData(end,1); %correct comp steps
                rData=[rData;rData2(2:end,:)];
            end
            if size(rData,2)>7;TD.checkPiston=1;else;TD.checkPiston=0;end
            
            %store data into properties
            if exist('trialData','var')
                TD=trialData;
                %Check if last TD element is the same of the new one that
                %will be added. If yes remove copy from the selection
                if TD.Step(end) == rData(1,1);k=2;
                else;k=1;
                end
                TD.Step=[TD.Step;rData(k:end,1)];
                TD.Fz=[TD.Fz;rData(k:end,2)];
                TD.Dz=[TD.Dz;rData(k:end,3)];
                TD.Fy1=[TD.Fy1;rData(k:end,4)];
                TD.Dy1=[TD.Dy1;rData(k:end,5)];
                TD.Fy2=[TD.Fy2;rData(k:end,6)];
                TD.Dy2=[TD.Dy2;rData(k:end,7)];
                if TD.checkPiston
                    TD.Fx1=[TD.Fx1;rData(k:end,8)];
                    TD.Dx1=[TD.Dx1;rData(k:end,9)];
                    TD.Fx2=[TD.Fx2;rData(k:end,10)];
                    TD.Dx2=[TD.Dx2;rData(k:end,11)];
                end
            else
                TD.Step=rData(:,1);
                TD.Fz=rData(:,2);
                TD.Dz=rData(:,3);
                TD.Fy1=rData(:,4);
                TD.Dy1=rData(:,5);
                TD.Fy2=rData(:,6);
                TD.Dy2=rData(:,7);
                if TD.checkPiston
                    TD.Fx1=rData(:,8);
                    TD.Dx1=rData(:,9);
                    TD.Fx2=rData(:,10);
                    TD.Dx2=rData(:,11);
                end
                %save geometric data
                TD.boxX = app.boxXEF.Value;     %box width - x
                TD.boxY = app.boxYEF.Value;     %box length - y
                TD.boxZ = app.boxZEF.Value;     %box height - z - also pistonY height
            end
            if ~isempty(app.SimType)
                %On LIGGGHTS, after consolidation, some pistons change
                %control from stress to displacement. On this process the
                %position is reset. This part recalculates the position
                %correctly.
                
                %correct pos before consolidation
                TD.Dz=TD.Dz+app.ConsoDpl(1,5);
                TD.Dx1=TD.Dx1+app.ConsoDpl(1,1);
                TD.Dx2=TD.Dx2+app.ConsoDpl(1,2);
                TD.Dy1=TD.Dy1+app.ConsoDpl(1,3);
                TD.Dy2=TD.Dy2+app.ConsoDpl(1,4);
                
                %correct pos after consolidation. First find the postion of
                %the step transition
                pos=find((TD.Step)>app.ConsoStep,1);
                if app.ExeType==2;pos=pos-1;end
                TD.Dz(pos:end)=TD.Dz(pos:end)+app.ConsoDpl(2,5);
                TD.Dx1(pos:end)=TD.Dx1(pos:end)+app.ConsoDpl(2,1);
                TD.Dx2(pos:end)=TD.Dx2(pos:end)+app.ConsoDpl(2,2);
                TD.Dy1(pos:end)=TD.Dy1(pos:end)+app.ConsoDpl(2,3);
                TD.Dy2(pos:end)=TD.Dy2(pos:end)+app.ConsoDpl(2,4);
            end
        end
        %Calculation Functions
        function snExt = extStrains(TD,steps,step0,app,varargin)
            %EXTSTRAINS calculates external strains
            %   This function use the data contained in the 'trialData'
            %   object to calculate the strain for all the calling 'steps'.
            %   The strain will be calculated starting at the step0. This
            %   function returns a NxD matrix, where N is the number of
            %   steps and D the dimension of calculation.
            
            if step0==0;step0=min(TD.Step);end
            
            %check the full lengths arrays
            if app.checkPiston
                w=TD.boxX-TD.Dx1+TD.Dx2;
            end
            l=TD.boxY-TD.Dy1+TD.Dy2;
            h=TD.boxZ+TD.Dz;
            %incremental strains in a lagrangian point of reference :
            %calculated in function of last step not in function of
            %initial state
            if app.checkPiston
                incE=[w,l,h];
            else
                incE=[l,h];
            end
            incE=(incE(2:end,:)-incE(1:end-1,:))./incE(1:end-1,:);
            incE=[zeros(1,size(incE,2));incE];
            
            %Deviatoric and principal if asked
            if nargin>4
                if app.checkPiston
                    Ed=sqrt(((incE(:,1)-incE(:,2)).^2 +...
                        (incE(:,2)-incE(:,3)).^2 +...
                        (incE(:,3)-incE(:,1)).^2 )/2);
                    Ev=((incE(:,1)+1).*(incE(:,2)+1).*(incE(:,3)+1)-1);
                    %Ev=sum(incE,2);
                else
                    Ed=(incE(:,2)-incE(:,1))/2;
                    Ev=((incE(:,1)+1).*(incE(:,2)+1)-1);
                    %Ev=sum(incE,2);
                end
                incE=[incE Ev -Ed]; %ed tont take int the GC convetion -
            end
            
            %calculate strain - negative for GC convention
            snExt=-cumsum(incE,1);
            
            %remove step0 Value and get only asked steps
            snExt=snExt((ismember(TD.Step,steps)),:)...
                -snExt((TD.Step==step0),:);
        end
        function ssExt = extStress(TD,steps,app)
            %EXTSTRESS calculates external stress
            %   This function use the data contained in the 'trialData'
            %   object to calculate the stress for all the calling 'steps'.
            %   It will also calculate the principal and deviatoric stress.
            %   This function returns a Nx(D+2) matrix, where N is the
            %   number of steps and D the dimension of calculation.
            
            %logical vector containing the position of the steps we want to use
            stId=ismember(TD.Step,steps);
            
            %%%%Calculate the surface of the pistons
            boxZarray=TD.boxZ+TD.Dz;
            boxYarray=TD.boxY-TD.Dy1+TD.Dy2;
            
            %Calculate the surface of pistons, *1000 to work on kPa
            if app.checkPiston %take pistonX area into account
                boxXarray=TD.boxX-TD.Dx1+TD.Dx2;
                SurfPX=boxZarray.*boxYarray*1000;
                SurfPY=boxXarray.*boxZarray*1000;
                SurfPZ=boxXarray.*boxYarray*1000;
            else
                if app.Bool3D
                    SurfPY=TD.boxX*boxZarray*1000;
                    SurfPZ=TD.boxX*boxYarray*1000;
                else
                    SurfPY=boxZarray*1000;
                    SurfPZ=boxYarray*1000;
                end
            end
            
            %Calculate the stress, Forces/Surfaces
            sigZ=TD.Fz(stId)./SurfPZ(stId);
            sigY1=-TD.Fy1(stId)./SurfPY(stId);
            sigY2=TD.Fy2(stId)./SurfPY(stId);
            sigY=(sigY1+sigY2)/2;
            if app.checkPiston %check existence of x direction pistons
                sigX1=-TD.Fx1(stId)./SurfPX(stId);    %per element division
                sigX2=TD.Fx2(stId)./SurfPX(stId);     %per element division
                sigX=(sigX1+sigX2)/2;                       %mean value
            end
            
            %Calculate the stress principal and deviatoric stresses
            if app.checkPiston
                q=sqrt(1/2*((sigX-sigY).^2+(sigY-sigZ).^2+...
                    (sigZ-sigX).^2 ));
                p=(sigX+sigY+sigZ)/3;
            else
                q=(sigZ-sigY)/2;
                p=(sigY+sigZ)/2;
            end
            
            %Returning Value
            if app.checkPiston
                ssExt=[sigX,sigY,sigZ,q,p];
            else
                ssExt=[sigY,sigZ,q,p];
            end
        end
        function TD  = inflectionPoints(TD,app)
            %VALUABLEPOINTS Calculates Ez of important points
            % This function will load all the external strain and stress
            % values. Then it will look for points of inflection of the
            % courbs of Ev (volumetric strain) and q (deviatoric stress).
            % These are important points that have a deep value on the
            % evolution of the simulation.
            
            csStp=find(TD.Step==app.ConsoStep);
            stp=TD.Step(csStp:end);
            snExt = extStrains(TD,stp,app.ConsoStep,app,1);
            ssExt = extStress(TD,stp,app);
            TD.InfPts.q='';
            TD.InfPts.p='';
            TD.InfPts.ez='';
            TD.InfPts.ev='';
            switch app.SimType
                case 1
                    %{
                    %OLD EXECUTION THROUGH MEAN VALUES. NOW I WILL JUST GET
                    MIN AND MAX EV AND Q
                    %get Ez, Ev and q
                    ez=snExt(:,app.Bool3D+2);
                    ev=snExt(:,end-1);
                    q=ssExt(:,end-1);
                    p=ssExt(:,end-1);
                    
                    %smooth the courbs with the values around them
                    ev=smoothdata(ev,'movmean',[15,15]);
                    q=smoothdata(q,'movmean',[15,15]);
                    
                    %Calculate the diference of values (q(t+1)-q(t)) and smooth the
                    %data with a gaussian courb
                    dEz=diff(ez);
                    ez=(ez(1:end-1)+ez(2:end))/2;
                    p=(p(1:end-1)+p(2:end))/2;
                    qbs=(q(1:end-1)+q(2:end))/2;
                    dEv = smoothdata(diff(ev)./dEz,'gaussian');
                    dQ = smoothdata(diff(q)./dEz,'gaussian');
                    
                    %get the first two minimal values of dEv. Find first negative
                    %value, then find first positive value after it.
                    v1=find(dEv<0,1);
                    v2=find(dEv(v1:end)>0,1)+v1-1;
                    if isempty(v2)
                        InfE=v1;
                    else
                        InfE=[v1;v2];
                    end
                    
                    v3=find(dQ<0,1);
                    v4=find(dQ(v3:end)>0,1)+v3-1;
                    if isempty(v4)
                        InfQ=v3;
                    else
                        InfQ=[v3;v4];
                    end
                    TD.InfPts.e=ez([InfE;InfQ]);
                    TD.InfPts.p=p([InfE;InfQ]);
                    TD.InfPts.q=qbs([InfE;InfQ]);
                    %}
                    [~,pmq]=max(ssExt(:,end-1));    %max q location
                    [~,pmev]=max(snExt(:,end-1));   %max ev location
                    
                    TD.InfPts.ez=snExt([pmev,pmq],end-2);
                    TD.InfPts.ev=snExt([pmev,pmq],end-1);
                    TD.InfPts.p=ssExt([pmev,pmq],end);
                    TD.InfPts.q=ssExt([pmev,pmq],end-1);
                    
                case 2
                    %In the case of Undrained files, inflection point will
                    %be considered as the point where qp is at its max
                    qp=ssExt(:,end-1)./ssExt(:,end);
                    [~,f]=max(qp);
                    %cmp=(abs(mx-qp)/mx)<0.005;
                    %f(1)=find(cmp,1,'last');
                    %f(2)=find(cmp,1,'first');
                    
                    TD.InfPts.ez=snExt(f,end-2);
                    TD.InfPts.ev=snExt(f,end-1);
                    TD.InfPts.p=ssExt(f,end);
                    TD.InfPts.q=ssExt(f,end-1);
                case 3
                    %on Qcst simulations : step where drained triaxial
                    %phase ends AND point where q values drops of 1%
                    
                    fA=find(stp==app.QcstStep(1));
                    
                    [mxQ,mxPos]=max(ssExt(: ,end-1));
                    fB=find((ssExt(mxPos:end,end-1)/mxQ)>0.99,1,'last')+mxPos;
                    
%                     % changed to q/p max too
%                     qp=ssExt(:,end-1)./ssExt(:,end);
%                     [~,f]=max(qp);
                    TD.InfPts.ez=snExt([fA,fB],end-2);
                    TD.InfPts.ev=snExt([fA,fB],end-1);
                    TD.InfPts.p=ssExt([fA,fB],end);
                    TD.InfPts.q=ssExt([fA,fB],end-1);
            end   
            
        end
        function pts = wallFakepoints(TD,app,steps,R,PartData,varargin)
            %WALLFAKEPOINTS points placed on wall location
            %   This function will calculate the position of the walls and
            %   pistons. It will then create points spaced of R. These
            %   points will then be used to increase the accuracy of the
            %   Delaunay Triangulation and the Voronoi calculations by
            %   imposing a exterior limit to which the grains will be
            %   connected.
            
            %Get the rows position
            row=(TD.Step==steps); %row in wich STEP N1 is found
            
            %The value Rm add a distance of Rm towards the exterior of the
            %experiment to all added FakePoints. So the points will be
            %placed a Rm distance away of the actual placement of the walls.
            if nargin>5
                Rm=varargin{1};
            else
                Rm=0;
            end
            
            %Take into account the changing piston position.
            height=TD.boxZ+TD.Dz(row); %actual height of the piston
            
            if isempty(PartData)
                %Vertices of the simulation box
                base=(TD.boxY-TD.Dy1(row)+TD.Dy2(row)); %actual size of the base
                Y=[TD.Dy1(row)-Rm,TD.Dy1(row)-Rm, base+TD.Dy1(row)+Rm,...
                    base+TD.Dy1(row)+Rm,TD.Dy1(row)-Rm];
                Z=[-Rm,height+Rm,height+Rm,-Rm,-Rm];
            else
                %Partial data vertices increased to take into account
                %grains that are partially inside the rectangle. Then turn
                %into a 1xN vector
                pos=PartData.RectPos;
                if PartData.AngleRotation==0
                    Y=([pos(1),pos(1),pos(1)+pos(3),pos(1)+pos(3),pos(1)] + 2*R*[-1,-1,1,1,-1]);
                    Z=([pos(2),pos(2)+pos(4),pos(2)+pos(4),pos(2),pos(2)] + 2*R*[-1,1,1,-1,-1]);
                    %check if vertice of partial data is depassing the pZ
                    ch=find(Z>(height+Rm));
                    if size(ch,2)>0
                        Z(ch)=ones(size(ch))*(height+Rm);
                    end
                else
                    %When the rectangle is rotated we first calculate the Y
                    %and Z values in relation to the rectangle's center.
                    %Then we rotate back and translate to the correct
                    %origin.
                    Y=((2*R+pos(3)/2)*[-1,-1,1,1,-1]);
                    Z=((2*R+pos(4)/2)*[-1,1,1,-1,-1]);
                end
            end
            
            %populate these edges with vertices
            pts=zeros(2);p=1;
            for i=1:(size(Z,2)-1)
                %odd values does vertical lines while even does horizontal.
                %'rem' returns the rest. If it is odd rest is 1 else 0
                if rem(i,2)
                    %odd - vertical lines
                    Z2=(min(Z(i:i+1)):(R):max(Z(i:i+1)))';
                    if Z(i)>Z(i+1);Z2=flip(Z2);end %flip if necessary
                    if Z2(end)~=Z(i+1);Z2=[Z2;Z(i+1)];end %#ok<AGROW> %add value if necessary
                    nb=size(Z2,1); %nb of points
                    pts(p:(p+nb-1),:)=[ones(nb,1)*Y(i) Z2];p=p+nb;
                else
                    %even - horizontal lines
                    Y2=(min(Y(i:i+1)):(R):max(Y(i:i+1)))';
                    if Y(i)>Y(i+1);Y2=flip(Y2);end %flip if necessary
                    if Y2(end)~=Y(i+1);Y2=[Y2;Y(i+1)];end %#ok<AGROW> %add end value if necessary
                    nb=size(Y2,1); %nb of points
                    pts(p:(p+nb-1),:)=[Y2 ones(nb,1)*Z(i)];p=p+nb;
                end
            end
            
            %Check third dimension
            if app.Bool3D
                nbPts=size(pts,1);p=1; %nb of points already created
                if app.checkPiston
                    X=[-(TD.boxX/2-TD.Dx1(row)+Rm) TD.boxX/2+TD.Dx2(row)+Rm];
                else
                    X=[-(TD.boxX/2+Rm) TD.boxX/2+Rm];
                end
                if ~isempty(PartData)
                    if PartData.Xlimit~=0
                        X=[-1 1]*(PartData.Xlimit/2+Rm);
                    end
                end
                X2=X(1):R:X(2);
                if X2(end)~=X(2);X2=[X2 X(2)];end %add end value if necessary
                nb=size(X2,2); %nb of points X
                %Mid points
                %For each X we have to repeat all the Y Z calculated
                %before, we will save them in a nez vector pts2.
                pts2=zeros(1,3);
                for i=2:(nb-1)
                    pts2(p:(p+nbPts-1),:)=[ones(nbPts,1)*X2(i) pts];p=p+nbPts;
                end
                pts=pts2;
                
                %Now we need to add the PistonX equivalent points (walls)
                y=(min(Y):(R):max(Y))';
                if y(end)~=max(Y);y=[y;max(Y)];end %add end value if necessary
                z=(min(Z):(R):max(Z))';
                if z(end)~=max(Z);z=[z;max(Z)];end %add end value if necessary
                
                pts2=[size(y,1)*size(z,1),2];k=1;
                for i=1:size(y,1)
                    pts2(k:(k+size(z,1)-1),:)=[y(i)*ones(size(z,1),1) z];
                    k=k+size(z,1);
                end
                sz=size(pts2,1);
                pts=[pts; [ones(sz,1)*X2(1) pts2];[ones(sz,1)*X2(end) pts2]];
                D=3;
            else
                D=2;
            end
            
            %Base correction for the angled rectangle case
            if ~isempty(PartData)
                if PartData.AngleRotation~=0
                    %rotation
                    pts(:,(D-1):D)=pts(:,(D-1):D)*(PartData.RotationMatrix)';
                    %translation
                    pts(:,D-1)=pts(:,D-1)+(pos(1)+pos(3)/2);
                    pts(:,D)=pts(:,D)+(pos(2)+pos(4)/2);
                    ch=find(pts(:,D)>(height+Rm));
                    if size(ch,1)>0
                        pts(ch,D)=ones(size(ch))*(height+Rm);
                    end
                end
            end
            
            pts=unique(pts,'rows');
        end
        function [V,subV] = getVolume(TD,app,steps,PD)
            %GETVOLUME calculates the volume of the simulation
            %   This function use the data contained in the 'trialData'
            %   object to calculate the volume or area of the simulation
            %   for all 'steps'. This calcualtion takes into account the
            %   existence of a 'partialData' object calculating only the
            %   volume or area selected by the user.
            
            %find the column containing piston that for the step we are
            %calculating.
            row=ismember(TD.Step,steps);
            subV='';  
            %Volume calculations
            if isempty(PD) %check for partial data
                V=(TD.boxY-TD.Dy1(row)+TD.Dy2(row)).*(TD.boxZ+TD.Dz(row));
                if app.Bool3D
                    if app.checkPiston 
                        V=V.*(TD.boxX-TD.Dx1(row)+TD.Dx2(row));
                    else
                        V=V.*TD.boxX;
                    end
                end 
            else
                %Calculate the Y Z area by comparing two polyshapes. This
                %way we dont need to check if the rectangle depass the box
                %dimension in any point.
                rect=polyshape([PD.Vertices]);
                sz=max(size(steps)); %nb of steps
                if sz==1 %only 1 Step asked
                    box=polyshape([TD.Dy1(row),0;...
                        TD.Dy1(row),TD.boxZ+TD.Dz(row);...
                        TD.boxY+TD.Dy2(row),TD.boxZ+TD.Dz(row);...
                        TD.boxY+TD.Dy2(row),0;
                        TD.Dy1(row),0]);
                    V = area(intersect(box,rect));
                else %many steps asked
                    boxDim=[permute(TD.Dy1(row),[3,2,1]),zeros(1,1,sz);...
                        zeros(1,1,sz),permute(TD.boxZ+TD.Dz(row),[3,2,1]);...
                        permute(TD.boxY+TD.Dy2(row),[3,2,1]),permute(TD.boxZ+TD.Dz(row),[3,2,1]);...
                        permute(TD.boxY+TD.Dy2(row),[3,2,1]),zeros(1,1,sz)];
                    V=zeros(sz,1);
                    for i=1:sz
                        box=polyshape(boxDim(:,:,i));
                        V(i) = area(intersect(box,rect));
                    end
                end
                
                %X dimension
                if app.Bool3D
                    if app.checkPiston %add X piston volume
                        width=(TD.boxX-TD.Dx1(row)+TD.Dx2(row));
                    else
                        width=TD.boxX;
                    end
                    if PD.Xlimit>0
                        width=min(width,PD.Xlimit);
                    end
                    V=V.*width;
                end
                
                %Check the existence of subdivisions
                if app.SubdivisionButton.Value && sz==1
                    subV=zeros(PD.SubLin*PD.SubCol,1);
                    for i=1:(PD.SubLin*PD.SubCol)
                        rect=polyshape(PD.SubVertices(:,:,i));
                        subV(i) = area(intersect(box,rect));
                    end
                    if app.Bool3D;subV=subV*width;end
                end
            end
        end

        %Getters Functions
        function w = getWidth(TD,steps) %x
            %GETWIDTH return the width of the box at the specified steps
            if isempty(TD.boxPX)
                w=TD.boxX;
            else
                w=TD.boxX+TD.Dx1-TD.Dx2;
                w=w(ismember(TD.Step,steps));
                
            end
        end
        function l = getLength(TD,steps) %y
            %GETLENGTH return the length of the box at the specified steps
            l=TD.boxY+TD.Dy1-TD.Dy2;
            l=l(ismember(TD.Step,steps));
        end
        function h = getHeight(TD,steps) %z
            %GETHEIGHT return the height of the box at the specified steps
            h=TD.boxZ+TD.Dz;
            h=h(ismember(TD.Step,steps));
        end
        
    end
end