function varargout = disperse(x,dim)
% DISPERSE was created so that you can stop doing things like this:
%
%   x1 = array(1); % ...repetitive assignments from an array
%   x2 = array(2);
%   x3 = array(3);
%   x4 = array(4);
%
% and start doing things like this:
%
%   [x1 x2 x3 x4] = disperse(array);
%
% DISPERSE generalizes to arbitrary dimensions, and is extended to follow
% analogous behavior on cell arrays and structure arrays. See the html
% documentation for more details and examples.
%
% Example:
%   Grab the color channels from an RGB image:
%   [r g b] = disperse(im);
% Sam Hallman
% shallman@uci.edu
% May 26, 2010
% num2cell on column vectors is problematic
if ndims(x)==2 && size(x,2)==1 %#ok<ISMAT>
    x = x';
end
if isnumeric(x) || ischar(x) || islogical(x) || isstruct(x)
    varargout = num2cell(x,dim);
elseif iscell(x)
    if size(x,1) == 1
        varargout = x;
    else
        varargout = num2cell(x,dim);
    end
else
    error('unknown data type');
end
% Sam Hallman (2022). disperse 
%(https://www.mathworks.com/matlabcentral/fileexchange/33866-disperse),
%MATLAB Central File Exchange. Retrieved June 9, 2022.
%later modified by Joao Chueire
end