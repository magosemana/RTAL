function [PA,PR,PS,PG,PGT]=AppGraph(app,varargin)
%APPGRAPH function will control app's Axes drawing on it accordingly
                
% Variables

%File related
N1=app.N1EF.Value;

%Geometry related
boxY=app.boxYEF.Value;
boxZ=app.boxZEF.Value;

%Program related
PD=app.PartialData;
TD=app.TrialData;
if isempty(TD.Dz);return;end

%plot related
ax=app.UIAxes;hold(ax,'on');
set(ax, 'XLimSpec', 'Tight');
set(ax, 'YLimSpec', 'Tight');
PR=app.PlotRect;
PS=app.PlotSubdv;
PG=app.PlotGrains;
PGT=app.PlotGrainsTot;
PA=app.PlotArea;

%Base Drawing
if any(strcmpi(varargin, 'BASE'))
    %base Plot
    row=(TD.Step==N1); %row in wich STEP N1 is found
    %Calculate Experiment walls for drawing
    height=boxZ+TD.Dz(row); %actual height of the piston
    Y=[TD.Dy1(row),TD.Dy1(row),(boxY+TD.Dy2(row)),(boxY+TD.Dy2(row)),TD.Dy1(row)];
    Z=[0,height,height,0,0];
    %first we will delete the existing plot/property then we will redraw it
    %so we can have the right handles in the app properties
    delete(PA);PA=plot(ax,Y,Z,'k','LineWidth',3);
end

%Rectangle plot
    %check if all conditions are met to start the rectangle partial draw
if any(strcmpi(varargin, 'RECTANGLE'))
    if ~isempty(PD) && isequal(upper(app.PartialSwitch.Value),'PARTIAL')
        if ~isempty(PD.RectPos)
            Vert=PD.Vertices;
            delete(PR);
            PR=plot(ax,[Vert(:,1);Vert(1,1)],[Vert(:,2);Vert(1,2)],'b');
        end
    end

    %Subdivision for rectangular Plot
    if app.SubdivisionButton.Value && exist('Vert','var')
        pos=PD.RectPos;
        %Subdivision Parts
        L=app.LinesEditField.Value; %Number of lines of area subdivision
        C=app.ColEditField.Value; %Number of columns of area subdivision
        Width=pos(3)/C; %width of area
        Height=pos(4)/L;%Height of area

        %Take rotation into account
        if PD.AngleRotation ~=0
            %rectangle position calculation
            RM=PD.RotationMatrix;
            %Calculate the rotation center of the rectangle on the
            %common base
            R=pos(1:2)+pos(3:4)/2;
            %Calculate the rotation center of the rectangle on a
            %rotated base
            R=(R*RM);
            %Calculate one vertice in the rotated base, we will work
            %with one vertice and Width and height onwards
            pos(1:2)=R-pos(3:4)/2;
        end 
        delete(PS);PS=hggroup(ax);
        for i=1:(C-1) %draw lines
            Y3=[pos(1),pos(1)]+i*Width;
            Z3=[pos(2),pos(2)+pos(4)];
            if PD.AngleRotation~=0 
                M=[Y3',Z3']*RM';
                Y3=M(:,1)';Z3=M(:,2)';
            end
            plot(Y3,Z3,'b','Parent',PS);
        end
        for j=1:(L-1) %draw columns
            Y3=[pos(1),pos(1)+pos(3)];
            Z3=[pos(2),pos(2)]+j*Height;
            if PD.AngleRotation~=0 
                M=[Y3',Z3']*RM';
                Y3=M(:,1)';Z3=M(:,2)';
            end
            plot(ax,Y3,Z3,'b','Parent',PS);
        end
    end
end

%Grains plot
if any(strcmpi(varargin, 'GRAINS'))
    if app.GrainsButton.Value
        %Get grains data
        gr=grains('BASIC',N1,'',app);
        Y=gr.Coord(:,2);Z=gr.Coord(:,3);R=gr.Radius;
        %check Partial
        if ~isempty(PD)
            PD = grainListing(PD,gr);
            Grl=PD.GrainsRectangle;
        end
        hold(ax,'on')
        NbGr=size(Y,1);
        CalcPanel(app,'','','on');
        delete(PG);PG=hggroup(ax); %partial grains
        delete(PGT);PGT=hggroup(ax); %allgrains
        for i=1:NbGr 
            y=Y(i);z=Z(i);r=R(i);%Calculates all grains
            rectangle('Position',[(y-r) (z-r) 2*(r) 2*(r)],'Curvature',[1,1],'Parent',PGT)
            if ~isempty(PD)
                if ~isempty(intersect(i,Grl))%Calculates grains inside partial
                    rectangle('Position',[(y-r) (z-r) 2*(r) 2*(r)],'Curvature',[1,1],'Parent',PG)
                end
            end
        end
        CalcPanel(app,'','','off');
    end
end
%Freehand plot
end

