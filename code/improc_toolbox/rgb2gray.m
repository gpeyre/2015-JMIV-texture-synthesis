function Y = rgb2gray(im)
% RGB 2 GRAY    - convert a color image to grayscale (using Luma).
%
%   Usage:
%       Y = rgb2gray(im);
%
%   Note:
%       the result is the luma of the inial image.

    coefs = [.299 .587 .114];
    h = reshape(coefs, 1, 1, 3);
    Y = convn(im, h, 'valid');
end