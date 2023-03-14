function [sGr,clCell] = findClusterGrain(type,sc,idGr,varargin)
%FINDCLUSTERGRAIN return cluster ID belonging to each grain in idGr
% Return a singleGrain object containing the ID n orer of the clusters that
% it is attached to.

Cl=numel(sc.Loops);
switch type
    case 'ForceChain'
        sGr=varargin{1}.SingleGrain(idGr); %profit the existence in "grains" and use as support
    case 'ClustTransf'
        %sGr(numel(idGr))=singleGrain();
        %do the same as the first part but for the cells now %clID - clOrder - cellID
        clCell=repelem([(1:Cl)',cat(1,sc.Loops.Order)],cat(1,sc.Loops.nbCells),1);
        clCell=[clCell cat(2,sc.Loops.sCells)'];
        %removed the single grain analysis for now.
        sGr='';
        return;
end
%Make sure Goodcell vector exists
if isempty(sc.GoodCells);sc.GoodCells = goodCell(sc);end

%First part - find out for each grain the clusters it is attached to.
%To do so first create a 3 column 2D matrix containing the ID of the
%cluster, the order and each grain. So for every grain of a Cluster the
%first two elements of the line is the same.  %ClusterID - ClusterOrder -
%grID
clGr=repelem([(1:Cl)',cat(1,sc.Loops.Order)],cat(1,sc.Loops.Size),1);
clGr=[clGr cat(1,sc.Loops.Grains)];
%do the same for the Cluster4 category %ClusterID - grID

cl4Id=sc.GoodCells(~ismember(sc.GoodCells,cat(2,sc.Loops.sCells)')); %id Cl4 cell location
clGr4=repelem((1+Cl):(numel(cl4Id)+Cl),4)'; %create ClusterID starting from L+1
clGr4=[clGr4 reshape(sc.DelaunayT(cl4Id,:)',[numel(cl4Id)*4,1]) ];
%Get both matrixes in one
clGr = [clGr;clGr4(:,1) ones(length(clGr4),1)*4 clGr4(:,2)];  %ClusterID - ClusterOrder - grID
%get the IDs of all the grains belonging to force chains

%Second part - for each grain create a singleGrain object that will contain
%the ClusterID and Order of clusters that contain this grain
for i=1:numel(idGr)
    sGr(i).ID=idGr(i);
    sGr(i).ClusterID=clGr(clGr(:,3)==idGr(i),1:2); %ClusterID - ClusterOrder
end

end

