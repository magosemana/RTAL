classdef singleGrain
    %SINGLEGRAIN object used to create force chains
    %   This object is created to facilitate the search for force chains.    
    properties
        %Force chain Values
        Contact         %ID of Highly stressed grains in contact with this grain
        ContactLine     %Lines formed by grain`s contact (for easy VTK draw)
        ContactNb       %Number of contacts (for easy fc length calc)
        ID              %ID of the grain
        Elevation       %Elevation angle of contacts
        Position        %Grain's position
        PrincipalStress %Max Eignen value of the stress tensor
        StressVector    %Major Stress vector
        %Cluster+Force chain values
        ClusterID       %used in loops and force chain analysis
        ClusterGrp      %used in loops and force chain analysis
        ClusterVal      %used in loops and force chain analysis
    end
    
    methods
        function sGr = singleGrain(sgID,str,pos)
            %SINGLEGRAIN Construct an instance of this class
            %   Object used in the force chain calculation. GrID, stress
            %   tensor and position needed for the creation. Contact
            %   information are added later in the force chain calculation
            %   algorythm.
            sGr.Contact=double.empty(0,1);
            sGr.ContactLine=double.empty(0,2);
            sGr.ContactNb=0;
            if ~nargin
                sGr.ID=-1;
            else
                sGr=basicInfo(sGr,sgID,str,pos);
            end
        end
        function sGr=basicInfo(sGr,sgID,str,pos)
            %BASICINFO Fill the previously constructed object with the
            %grain information
            sGr.ID=sgID;
            [dir,eVal]=eig(str);  %get eig vals
            [mx,I]=max(diag(eVal));
            sGr.PrincipalStress=mx;
            sGr.StressVector=dir(:,I)';
            sGr.Position=pos((4-size(str,1)):3);
        end
    end
end

