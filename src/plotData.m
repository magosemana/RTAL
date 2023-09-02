classdef plotData
    %PLOTDATA will contain information that is necessary for plotting
    %   This class will create an object that will be used to
    
    properties
        Bool3D      %Boolean check 3D
        BoolPx      %Boolean piston X
        ConsoStrain %Strain where consolidation ended
        ConsoTime   %Time where consolidation ended
        Interval    %Calculation interval
        InfPts      %Inflection points of Ev and q
        FileName    %File name in case of "Load"
        N1          %Calculation starting step
        N2          %Calculation endig step
        Results     %Calculation results
        Path        %Path for the file
        Prefix      %Prefix for file saving
        SubL        %Subdivision lines (if exists)
        SubC        %Subdivision columns (if exists)
        TimeStep    %Timestep value in seconds
        SimType     %simulation type
        SlopeM      %used only on load ExtF
    end
    
    methods
        function pD = plotData(type,results,varargin)
            %PLOTDATA Construct an instance of this class
            %   Atributes the given information to the object. This object
            %   will then be used to plot diferent types of calculations.
            if ~nargin
                pD.Bool3D=0;
                return
            end
            pD.Results=results;
            switch type
                case "Normal"
                    app=varargin{1};
                    pD.Prefix=varargin{2};
                    %save app atributes into the pD
                    pD.Bool3D=app.Bool3D;
                    pD.BoolPx=app.checkPiston;
                    pD.ConsoTime=app.ConsoStep*app.TimeStep;
                    pD.Interval=app.CalcInt.Value;
                    pD.N1=app.N1EF.Value;
                    pD.N2=app.N2EF.Value;
                    pD.InfPts=app.TrialData.InfPts;
                    pD.TimeStep=app.TimeStep;
                    pD.SimType=app.SimType;
                    if nargin>4;pD.ConsoStrain=varargin{3};end
                    if nargin>5  %called on a subdivision calculation
                        pD.SubL=varargin{4};
                        pD.SubC=varargin{5};
                    end
                case "Load"
                    pD.N1=varargin{1};
                    pD.N2=varargin{2};
                    pD.Interval=varargin{3};
                    pD.Bool3D=varargin{4};
                    pD.BoolPx=varargin{5};
                    pD.TimeStep=varargin{6};
                    pD.ConsoStrain=varargin{7};
                    pD.ConsoTime=varargin{8};
                    pD.InfPts=varargin{9};
                    pD.SubL=varargin{10};
                    pD.SubC=varargin{11};
                    pD.FileName=varargin{12};
                    pD.Prefix="Load";
                    pD.Path=varargin{13};
                case "FastCreation"
                    pD.Prefix=varargin{1};
            end
            
        end
    end
end
