function clrRGB = graphClrCode(nbPts)
%GRAPHCLRCODE return RGB triplets for 'nbPts' nice plot colors 
%   'jet' colorcode from matlab uses the following intervals to choose
%   colors : [0,0,1] => [0,1,1] => [1,1,0] => [1,0,0]. Use these values as
%   bases to return spreadout colors for the graph. Also sum(RGB)<=2.
clrRGB=[0,0,1;
    0,1,1;
    1,0,0];
d=size(clrRGB,1);
if nbPts<=d;return;end
%nb of extra points in to get between each intv
pts=zeros(d-1,1);
pts=pts+floor((nbPts-d)/(d-1)); %exact divisions
r=rem(nbPts-d,d-1);           %extra from divisions
pts(1:r)=pts(1:r)+1;
%distribute values and get the triplets
rgbs=cell(d-1,1);
for i=1:d-1
    int=1/(pts(i)+1);
    dif=clrRGB(i+1,:)-clrRGB(i,:);
    if i==d-1
        rgbs{i}=ones(pts(i)+2,1)*clrRGB(i,:)+((0:int:1)')*dif;
    else
        rgbs{i}=ones(pts(i)+1,1)*clrRGB(i,:)+((0:int:(1-10^-9))')*dif;
    end
end
clrRGB=cat(1,rgbs{:});
end