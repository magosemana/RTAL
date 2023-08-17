classdef grains
    %GRAINS contains the variables for each grain at at instant T
    %   This is the most important class of the simulation.
    
    properties
        %Properties read from files
        ContactNubr     %grain's number of contacts	- Nx1 array
        Coord           %grain's X,Y,Z coordinates	- Nx3 matrix
        Displacement    %grain's displacement       - Nx3 matrix
        ID              %grain's ID                 - Nx1 array
        KEnergy         %grain's ID Knectic energy  - Nx1 array
        Nb              %Number of grains           - scalar
        Radius          %grain's radius             - Nx1 array
        PGStressTensor  %Per grain stress tensor    - 3x3xN matrix
        
        %Properties calculated through this funciton
        AveCluster      %grain's average cluster category from [1:4] - Nx2 matrix
        ContactsDir     %grain's contacts number, total vert and horiz. Nbgrainsx3 matrix
        ContactFNormals %grain's contact normals calculated with the force vectors
        FabricTensor    %model's fabric tensor, geometric and force - 3x3x2 matrix
        ForceChains     %forceChain object array - Nx1
        HighStrGrains   %nb of highly stressed grains
        PGFabricTensor  %Per grain fabric tensor, geometric and force - 3x3XNgrainsx2 matrix
        SingleGrain     %singleGrain objects
        StressedGrains  %number of grains highly stressed - scalar
        ThreeP          %minimal force chain element - Nx3 matrix
        Z               %coordination number (average contact) - scalar
        
    end
    methods
        %Calling function
        function [gr,PartialData] = grains(mode,step,PartialData,app)
            if ~nargin
                gr.ID=0;
                return
            end
            %check which folder to load and filename
            if nargin>2
                rData = readData('GRAINS',app,step);
            else
                rData = readData('GRAINSTEST',app,step);
                PartialData='';
            end
            if isempty(rData);return;end
            
%             %the following code turn a 3D model into a 2D model taking only
%             %the grains in the surface located in positive x.
%             %NOT BEING USED FOR NOW
%             if rData==0
%                 maxR=max(rData(:,5));       %maximum radius
%                 boxWidth=app.boxXEF.Value;
%                 %find the grains that are close to the wall
%                 xlim=boxWidth/2-maxR*1.5;
%                 rData=rData(rData(:,2)>xlim,:);
%             end
            %sort following the Ids
            rData=sortrows(rData);
            %Save object properties
            gr.ID=rData(:,1);
            gr.Nb=numel(gr.ID);
            if max(gr.ID)~=gr.Nb
                %Functions that take into account the contacts work as if
                %the ID of the grain is the position it is on the 
                fprintf("Number of grains %d differs from the max ID %d\n",...
                    gr.Nb,max(gr.ID))
            end
            gr.Coord=rData(:,2:4);
            gr.Radius=rData(:,5);
            ContNbr=rData(:,6);
            gr.Displacement=rData(:,7:9);
            %gr.KEnergy=M(:,16);
            
            %Check for partial selection.
            if ~isempty(PartialData)
                %identify grains inside the chosen rectangle.
                PartialData = grainListing(PartialData,gr);
                if app.SubdivisionButton.Value==0
                    %we will limit the contact number to only the grains with
                    %at least one part inside the rectangle.
                    ContNbr=ContNbr(PartialData.GrainsRectangle);
                end
            end
            %Save contact number properties
            gr.ContactNubr=ContNbr;
            gr.Z=mean(ContNbr);
            %Check for extra calculations. varargin{2} must contain the
            %name of the second file to be read, containing contact forces.
            switch upper(mode)
                case 'BASIC'; return;
                case {'STENSOR','FORCECHAIN'}
                    %transform pgST from a Ngrainsx6 matrix in a 3x3xNgrains
                    pgST=permute(rData(:,10:15),[2 3 1]); %6x1xNgrains
                    pgST=[pgST(1,:,:) pgST(4,:,:) pgST(5,:,:);...
                        pgST(4,:,:) pgST(2,:,:) pgST(6,:,:);...
                        pgST(5,:,:) pgST(6,:,:) pgST(3,:,:)];
                    gr.PGStressTensor=pgST;
                    %gr = grStress(gr,step,PartialData,app);
                    if isequal(upper(mode),'FORCECHAIN')
                        gr = forcechains(gr,step,app);
                    end
                case 'CONTACT'
                    gr = contact(gr,step,PartialData,app);
                case 'CONTACTNBR'      %contact
                    gr = contactNbr(gr,PartialData,app);
            end
        end
        %Object functions
        function gr = contact(gr,step,PartialData,app)
            %CONTACT this function will read the file containing the ID of
            %the grains in contact and the force between them obtained from
            %LIGGGHTS. It will then calculate the following things :
            %
            %angleForce : is a vector containing all the angles of the
            %grains's forces.
            %
            %nbContacts : matrix containing total, vertical and horizontal
            %contacts. NbGrainsx3 matrix
            %
            %The following tensors are obtained through multiplication of
            %column x line vectors:
            %
            %Fabric Tensor : this tensor will contain information about the
            %contact distribution. Obtained through the multiplication of
            %the contact normal. Symetric tensor.
            %
            %Per Grain Fabric Tensor : just like the Fabric tensor but
            %divided for each diferent grain.
            
            %get filename and read files related to contact forces
            rData = readData('CONTACT',app,step);
            if isempty(rData);return;end
            rData=rData(:,1:5);
            
            %Check partial data to keep only contacts that at least one of
            %the grains is inside the rectangle
            if ~isempty(PartialData)
                grList=PartialData.GrainsRectangle;
                mbm=ismember(rData(:,[1,2]),grList);
                rData=rData(mbm(:,1)|mbm(:,2),:);
            else
                grList=unique(rData(:,1:2));
            end
            
            %Prepare base values
            nbCtcs=size(rData,1);       %nb of contacts
            nbGr=size(gr.Radius,1); %nb of grains
            if app.Bool3D;D=3;else;D=2;end
            %Matrixes and variables that will contain the data.
            pcFT=zeros(D,D,nbCtcs,4);	%per Contact fabric tensor
            %(geometric, mechanical,normal force and tangencial force anisotropy)
            
            %%% PER CONTACT VALUES %%%
            %FORCE rDatanormal
            F=rData(:,6-D:5);        %all forces
            normF=vecnorm(F,2,2);   %force magnitude
            f=F./normF;             %unitary force vector
            %BRANCH VECTOR normal
            l=(gr.Coord(rData(:,2),4-D:3)-gr.Coord(rData(:,1),4-D:3)); %all branch
            l=l./vecnorm(l,2,2); %unitary branch vector
            
            %ANG between F and l
            cosang=sum(f.*l,2); %cosinus of ang between F and l (dot product)
            fn=normF.*cosang.*l;%force projected into the branch vector
            ft=F-fn;	%rest of the force
            
            %HOR and VER forces
            if app.Bool3D %atan(Fz/sqrt(Fx^2+Fy^2))
                angleB=abs(atan(l(:,3)./sqrt(l(:,1).^2+l(:,2).^2)))>=(pi()/4);
            else %atan(Fz/Fy)
                angleB=abs(atan(l(:,2)./l(:,1)))>=(pi()/4);
            end
            
            %TENSORS per contact - 3x3xNctctsx2
            l=permute(l,[3,2,1]);	%contact ani
            pcFT(:,:,:,1)=pagemtimes(permute(l,[2,1,3]),l);
            nF=permute(f,[3,2,1]);
            pcFT(:,:,:,2)=pagemtimes(permute(nF,[2,1,3]),nF);
            %pcFT(:,:,:,2)=pagemtimes(permute(F,[2,3,1]),nF);
            pcFT(:,:,:,3)=pagemtimes(permute(fn,[2,3,1]),l);
            pcFT(:,:,:,4)=pagemtimes(permute(ft,[2,3,1]),l);
            
            if ~app.SubdivisionButton.Value
                cDir=[nbCtcs, sum(angleB) sum(~angleB)];
                FT=permute(sum(pcFT,3),[1 2 4 3]);
            else
                %%% PER GRAIN VALUES %%%
                pcDir=zeros(nbCtcs,3);  %nb of Total, Vertical and Horizontal contacts per contact
                pcDir(angleB,[1,2])=1;	%Vertical contacts
                pcDir(~angleB,[1,3])=1;	%Horizontal contacts
                cDir=zeros(nbGr,3);     %nb of Total, Vertical and Horizontal contacts per grain
                pgFT=zeros(D,D,nbGr,2);	%per grain fabric tensor
                for i=1:size(grList,1)
                    %find to wich contacts the grain i belongs to.
                    ctcI=ismember(rData(:,[1,2]),grList(i));
                    ctcI=ctcI(:,1)|ctcI(:,2);
                    %sum contact directions and fabric tensor that belongs
                    %to the certain grain
                    cDir(grList(i),:)=sum(pcDir(ctcI,:),1);
                    pgFT(:,:,grList(i),:)=sum(pcFT(:,:,ctcI,:),3);
                end
            end
            %Save values onto the object
            gr.ContactFNormals=[rData(:,1:2) f];
            gr.ContactsDir=cDir;
            if ~app.SubdivisionButton.Value
                FT=FT/cDir(1);
                FT(:,:,2)=FT(:,:,2)/norm(FT(:,:,2));
                FT(:,:,3)=FT(:,:,3)/norm(FT(:,:,3));
                FT(:,:,4)=FT(:,:,4)/norm(FT(:,:,4));
                gr.FabricTensor=FT;
            else
                gr.PGFabricTensor=pgFT;
            end
        end
        function gr = forcechains(gr,step,app)
            %SINGLEGRAIN this function will create singleGrain objects
            %that will be used to facilitate the forcechain calculations
            
            %First step : find the grains that have higher than average
            %stress value. pgST needs to be divided by the volume of each
            %grain actually be the stress
            pgST=gr.PGStressTensor./permute((4*pi()/3)*gr.Radius.^3,[3,2,1]);
            sGr(size(gr.Radius,1))=singleGrain(); %prepare singleGrain object array
            nG=size(pgST,3); %nb of grains
            for i=1:nG
                sGr(i)=basicInfo(sGr(i),i,-pgST(:,:,i),gr.Coord(i,:));
            end
            pStress=cat(1,sGr.PrincipalStress);
            k=app.FCForceCoefEF.Value;%user chosen value to define HSG
            hsgID=find(pStress(:,1)>=k*mean(pStress(:,1))); %ID of Highly Stressed Grains
            gr.HighStrGrains=numel(hsgID); %nb of HSG
            
            %Load contact forces
            rData = readData('CONTACT',app,step);
            if isempty(rData);return;end
            ctc=rData(:,1:2);
            
            %Keep only the contacts between HSG grains
            %Keep only the contacts between HSG grains
            ctc=ctc(sum(ismember(ctc,hsgID),2)>1,:);
            
            %For a contact to be in a force chain, the most compressive
            %principal stress of each grain in the chain must point in the
            %direction of the next grain in the chain in both directions,
            %forward(A=>B) and backwards(B=>A).
            chcker=cosd(app.FCAngEF.Value); %max acceptable value
            %Get branch vectors
            l=cat(1,sGr(ctc(:,2)).Position)-cat(1,sGr(ctc(:,1)).Position);
            %Check forward
            stF=cat(1,sGr(ctc(:,1)).StressVector);
            checkF = dot(l,stF,2)./(vecnorm(l,2,2).*vecnorm(stF,2,2));
            %Check Backwards
            stB=cat(1,sGr(ctc(:,2)).StressVector);
            checkB = dot(l,stB,2)./(vecnorm(l,2,2).*vecnorm(stB,2,2));
            %Get contacts that validate both passes
            pass=(abs(checkF)>=chcker & abs(checkB)>=chcker);
            ctc=ctc(pass,:);
            checkF=checkF(pass,:);
            l=l(pass,:);
            %Identify if the contact direction A=>B is the same of the
            %stress direction. If the angle of checkF is [90,270]
            %(negative cosinus) then it is not and the contact need to
            %be fliped to B=>A
            f=(checkF<0);
            ctc(f,:)=[ctc(f,2) ctc(f,1)];
            %Calculate elevation of the contact
            ang=abs(atand(l(:,3)./(sqrt(l(:,1).^2+l(:,2)).^2)));
            
            %for the contact is validated as a part of a force chain the
            %contact information of both singleGrains object is updated.
            for i=1:length(hsgID)
                grID=hsgID(i);
                f=(ctc(:,1)==grID | ctc(:,2)==grID);
                if ~any(f);continue;end %no contacts identified
                sGr(grID).Contact=setdiff(unique(ctc(f,:)),grID);
                sGr(grID).ContactLine=ctc(f,:);
                sGr(grID).ContactNb=sum(f);
                sGr(grID).Elevation=ang(f);
            end
            
            %Obtain the three element force chains
            gr = tpfcElements(gr,ctc);
            
            %check all grains that have at least 1 contact with HSG
            cIds=cat(1,sGr.ContactNb);
            cIds=find(cIds>0);
            
            %start checking for chains
            nG=numel(cIds); %nb of HSG with at least 1 ctc
            fc(ceil(nG/3))=forceChain();k=1;    %force chain object and counter
            usd=zeros(nG,1);u=1;                %used array and counter
            for i=1:nG
                %Check if grain was already used by this loop or if the
                %grain has no HSG contacts
                if isempty(sGr(cIds(i)).Contact);continue;end
                if ismember(sGr(cIds(i)).ID,usd);continue;end
                nbGrCon=1;    %numer of HSG in contact
                grCon=[sGr(cIds(i)).ID; sGr(cIds(i)).Contact]; %add grs contacts ID
                nbCnt2=numel(grCon); %added sGrain contacts to the count
                while nbGrCon~=nbCnt2
                    nbGrCon=nbCnt2;
                    grCon=unique(cat(1,sGr(grCon).Contact));
                    %make sure all vertices are clustI on a colum : Nx1 vectors
                    if size(grCon,1)>1;grCon=grCon';end
                    nbCnt2=numel(grCon);
                end
                
                %force chains are formed bt the conjuction of at least
                %three grains
                if nbGrCon >2
                    fc(k) = forceChain(sGr(grCon)); k=k+1;
                end
                %aded used grains to the used array so they do not repeat
                s=numel(grCon);
                usd(u:(u-1+s))=grCon;
                u=u+s;
            end
            gr.SingleGrain=sGr;
            gr.ForceChains=fc(1:(k-1));
        end
        function gr = tpfcElements(gr,ctc)
            %TPFCELEMENTS calculates all elementary fc partitions
            % Uses the fc validated contact list and grains ID to return a
            % list containing the ID of grains forming threep (3 particle
            % force chains). To do so, for each line of the contact list,
            % try to find other contact conecting backwards or forward to
            % it. Save found threep and go next contact. However do not
            % check any previously analysed contact to save time.
            % Function will be based on the fact that the contact list is
            % ordered in the direction of the force propagation. Thus the
            % line i contain the ID of two grains with the force going from
            % column 1 to column 2. A three element containing line i must
            % have another connection either starting with element of
            % column 2 or ending with element of column 1.

            threep=zeros(5*size(ctc,1),3);k=1;
            for i=1:(size(ctc,1)-1)
                %Check forward, if there is any contact that starts with
                %the grain ending the contact i
                f=find(ctc(i+1:end,1)==ctc(i,2));
                if ~isempty(f)
                    t=numel(f);
                    threep(k:k+t-1,:)=[ones(t,1)*ctc(i,:) ctc(f,2)];
                    k=k+t;
                end
                %Check backwards, if there is any contact that ends with
                %the grain starting the contact i
                f=find(ctc(i+1:end,2)==ctc(i,1));
                if ~isempty(f)
                    t=numel(f);
                    threep(k:k+t-1,:)=[ctc(f,1) ones(numel(f),1)*ctc(i,:)];
                end
            end
            threep=unique(threep,'rows');
            gr.ThreeP=threep;
        end
        function gr = contactNbr(gr,PartialData,app)
            %CONTACTNBR     calculates an aproximation of the real contact
            %values by removing all the contacts with the walls. This will
            %be made by identifying the grains that are in the exterior of
            %the simulation.
            
            %get the position (line) of grains in the exterior
            if app.Bool3D
                k=convhull(gr.Coord);
            else
                k=convhull(gr.Coord(:,2:3));
            end
            %Vector containing the IDs of the grains that are on the
            %exterior
            grID=gr.ID(unique(k));
            
            %Load grains variables
            GrCnbr=gr.ContactNubr; % grains contact number
            if isempty(PartialData)
                IDs=gr.ID;
            else
                IDs=PartialData.AllGrains;
            end
            
            [~,ia]=intersect(IDs,grID);%get
            GrCnbr(ia)=[];
            %Save files
            gr.Z=mean(GrCnbr);
            gr.ContactNubr=GrCnbr;
        end
        function gr = aveCluster(gr,sc)
            %AVECLUSTER return the ave cluster per grain
            gC=sc.GoodCells;
            if isempty(gC)
                gC = goodCell(sc);
            end
            %Get clusters cells
            cC=cat(2,sc.Loops.sCells)';
            %Transform them into goodcell Ids
            [~,cC]=ismember(cC,gC);
            %add cluster order to it
            clCe=[cC repelem(cat(1,sc.Loops.Order),cat(1,sc.Loops.nbCells),1)];
            %get a list grain forming good cells only from DT
            gcDT=sc.DelaunayT(gC,:);
            aveCl=zeros(numel(gr.ID),1);
            for i=1:numel(gr.ID) 
                %get good cells belonging to grain i
                grCl=find(sum(gcDT==i,2)==1);
                %mathc it with clCe vector
                [mnb,pos]=ismember(grCl,clCe(:,1));
                %calculate average
                if sum(mnb)>1
                    aveCl(i)=( sum(~mnb)*4 + sum(clCe(pos(mnb),2)) )/numel(mnb);
                else
                    aveCl(i)=4;
                end
            end
            gr.AveCluster=aveCl;
        end
        function gr = purge(gr)
            %PURGE remove excess values from grain object for saving
            gr.Displacement='';
            gr.ContactNubr='';
            gr.Z='';
        end
        %Testing
        function [V] = avrgA(gr,step,app)
            %get filename and read files related to contact forces
            rData = readData('CONTACT',app,step);
            if isempty(rData);return;end
            rData=rData(:,1:5);
            D=3;
            %%% PER CONTACT VALUES %%%
            %FORCE rDatanormal
            F=rData(:,6-D:5);        %all forces
            normF=vecnorm(F,2,2);   %force magnitude
            %BRANCH VECTOR normal
            d=(gr.Coord(rData(:,2),4-D:3)-gr.Coord(rData(:,1),4-D:3));
            r1=gr.Radius(rData(:,1));r2=gr.Radius(rData(:,2));
            d=vecnorm(d,2,2)-r1-r2;
            md=mean(2*d./(r1+r2));
            R=r1.*r2./(r1+r2);
            A=2*pi()*d.*(2*R);
            V=abs(normF./A);
            V=V(V>0);
            V=[max(V) min(V) mean(V) md];
        end
    end %end method
end %end class

%{
        function gr = grStress(gr,step,PartialData,app)
            %get filename and read files related to contact forces
            rData = readData('CONTACT',app,step);
            if isempty(rData);return;end
            rData=rData(:,1:5);
            
            if app.Bool3D;D=3;else;D=2;end
            %Check partial data to keep only contacts that at least one of
            %the grains is inside the rectangle
            if ~isempty(PartialData)
                grList=PartialData.GrainsRectangle;
                mbm=ismember(rData(:,[1,2]),grList);
                rData=rData(mbm(:,1)|mbm(:,2),:);
            else
               % grList=unique(rData(:,1:2));
            end
            
            %branch vectors - N v 3 vector
            l=(gr.Coord(rData(:,2),4-D:3)-gr.Coord(rData(:,1),4-D:3));
            %force vector - N v 3 vector
            f=rData(:,6-D:5);
            strGr= f'*l;
            
            %wall - grain contact
            rData = readData('WCONTACT',app,step);
            if isempty(rData);return;end
            %branch vectors - N v 3 vector
            l=(gr.Coord(rData(:,2),4-D:3)-rData(:,9-D:8));
            %force vector - N v 3 vector
            f=rData(:,6-D:5);
            strW=f'*l;
            
            %save to property
            gr.StressTensor=strGr-strW;
        end
%}