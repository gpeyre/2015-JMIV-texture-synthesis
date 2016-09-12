function h = imshow(im, prop, value)
% IMSHOW        - display an image.
%
%   Basic replacement function for the image processing toolbox function.
%
%   Usage:
%       imshow(im);
%       imshow(im, 'Parent', ax);
%       h = imshow(...);
%
%   Description:
%       Display an image.
%       Other arguments (like 'Parent') can be provided.

    % Optional arguments
    if (exist('prop', 'var') && ~strcmpi(prop, 'parent'))
        error('invalid second argument');
    elseif (exist('value', 'var'))
        ax = value;
    else
        ax = gca;
    end;
    
    % Call the function
    h = image(im, 'Parent', ax);
    axis(ax, 'equal', 'off');
    
end