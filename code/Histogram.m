classdef Histogram
% HISTOGRAM     - provide functions to manipulate image histograms.
%
%   See also:
%       Histogram.demo
%       Histogram.transfer
%
%   Contact:
%       Guillaume.Tartavel @ telecom-paristech.fr
%       Downloadable on perso.enst.fr/~tartavel/


methods (Static)

    function [imF, E] = transfer(im, imRef)
    % TRANFER   - transfer the histogram of an image to another.
    %
    %   Usage:
    %       imF = Histogram.transfer(im, imRef);
    %       [imF, E] = ...
    %
    %   Arguments:
    %       im      - an image (gray or color)
    %       imRef   - another image
    %
    %   Output:
    %       imF     - image like 'im' with the hsitogram of 'imRef'.
    %       E       - squared distance between the histograms.
    %
    %   Notes:
    %       - for gray images, the transfer is exact.
    %       - for color images, it is a fast approximation using a PCA.
    %       - for gray images with the same number of pixels:
    %           imF = arg min || im - imBis ||
    %         among all imBis obtained by permutation of imRef's pixels.

        % Deal with color
        if (size(im, 3) == 1 && size(imRef, 3) > 1)
            im = repmat(im, size(imRef, 3));
        elseif (size(im, 3) > 1 && size(imRef, 3) == 1)
            im = rgb2gray(im);
        end;

        % Get the values as a vector
        X = reshape(im, [], size(im, 3))';
        Xref = reshape(imRef, [], size(imRef, 3))';

        % Transfer values and reshape back
        if (size(im, 3) == 1)
            [Xf, E] = Histogram.p_transferPoints1d(X, Xref);
        else
            [Xf, E] = Histogram.p_transferPointsPca(X, Xref);
        end;

        % Reshape back to an image
        imF = reshape(Xf', size(im));
    end;

    function d2 = distance2(im, imRef)
    % DISTANCE 2    - compute the distance between the histogram.
    %
    %   Usage:
    %       d2 = Histogram.distance2(im, imRef);
    %
    %   Argument:
    %       im, imRef   - images of the same size
    %
    %   Return:
    %       the (squared) distance between the histogram, approximated
    %
    %   Note:
    %       this function works only for images of the same size
    %
    %   See also:
    %       Histogram.transfer

        % Check size
        if (any(size(im) ~= size(imRef)))
            error('Histogram distance only implemented for images of same size');
        end;

        % Get pixel values
        Nc = size(im, 3);
        Xa = reshape(im, [], Nc)';
        Xb = reshape(imRef, [], Nc)';
        X = [Xa, Xb];

        % PCA basis
        [base, o] = eig(X * X');
        Ya = sort(base' * Xa, 2);
        Yb = sort(base' * Xb, 2);

        % Distance
        d2 = sum((Ya(:) - Yb(:)) .^ 2);
    end;

    function demo()
    % DEMO      - launch a demo of the histogram transfer.
    %
    %   Usage:
    %       Histogram.demo

        % Name of the images to be loaded
        imNameA = 'mandril512.png';
        imNameB = 'peppers.png';

        % Load images an gray versions
        imcA = double(imread(imNameA)) / 255;
        imcB = double(imread(imNameB)) / 255;
        imA = rgb2gray(imcA);
        imB = rgb2gray(imcB);

        % Histogram transfers
        imF = Histogram.transfer(imA, imB);
        imcF = Histogram.transfer(imcA, imcB);

        % Display
        clf;
        subplot 331; imshow(imcA);  title('Image');
        subplot 332; imshow(imcB);  title('Target histogram');
        subplot 333; imshow(imcF);  title('Tranfered (approximated)');
        subplot 334; imshow(imA);
        subplot 335; imshow(imB);
        subplot 336; imshow(imF);
        bins = (0:255) / 255;
        subplot 337; hist(imA(:),  bins); axis tight;
        subplot 338; hist(imB(:),  bins); axis tight;
        subplot 339; hist(imF(:),  bins); axis tight;

    end;

end;


%%%%%%%%%%  PRIVATE FUNCTIONS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

methods (Static)

    function [Xf, E] = p_transferPoints1d(X, Xref, interpolation)
    % Private -- TRANSFER POINTS 1D     - perform 1D transfer, row by row.
    %
    %   Usage:
    %       Xf = Histogram.p_transferPoints1d(X, Xref);
    %       Xf = Histogram.p_transferPoints1d(X, Xref, interpolation);
    %       [Xf, E] = ...
    %
    %   Arguments:
    %       X    (dxN )     - a set of N vectors in R^d
    %       Xref (dxN')     - a set of N' vectors
    %       interpolation   - see 'METHOD' in 'help INTERP1'
    %                         default interpolation is 'NEAREST'
    %
    %   Return:
    %       Xf  (dxN)       - a set of N vectors in R^d
    %       E               - squared distance between the histograms
    %
    %   Note:
    %       histogram transfer is performed independently on the rows:
    %       the row Xf(i,:) is obtained as the permutation of Xref(i,:)
    %       minimizing the distance to X(i,:).

        % Size and indices
        [d,  N]  = size(X);
        [dd, NN] = size(Xref);
        indices = @(M) (0:M-1) / (M-1);

        % Optional arguments
        if (~exist('interpolation', 'var') || isempty(interpolation))
            interpolation = 'nearest';
        end;

        % Check arguments
        if (dd ~= d)
            error('Invalid dimension "d" for Xref');
        end;

        % Target histogram
        Xsort = sort(Xref, 2);
        Xvalues = interp1(indices(NN)', Xsort', indices(N)', interpolation)';

        % Ordered input values
        Xinput = X + 1e-12*rand(size(X));   % dequantization
        [osef,ranks] = sort(Xinput, 2);
        rowsInd = repmat((1:d)', 1, N);

        % Transfer histogram
        Xf = zeros(d, N);
        Xf(sub2ind([d, N], rowsInd, ranks)) = Xvalues;

        % Error
        E = sum((Xf(:) - X(:)) .^ 2);
    end;

    function [Xeq, E] = p_transferPointsPca(X, Xref, interpolation)
    % Private -- TRANSFER POINTS PCA    - approximate n-D transfer.
    %
    %   Usage:
    %       Xf = Histogram.p_transferPointsPca(X, Xref);
    %       Xf = Histogram.p_transferPointsPca(X, Xref, interpolation);
    %       [Xf, E] = ...
    %
    %   Arguments:
    %       X    (dxN )     - a set of N vectors in R^d
    %       Xref (dxN')     - a set of N' vectors
    %       interpolation   - see 'METHOD' in 'help INTERP1'
    %                         default interpolation is 'NEAREST'
    %
    %   Return:
    %       Xf  (dxN)       - a set of N vectors in R^d
    %       E               - squared distance between the histograms
    %
    %   Note:
    %       this is an approximation, better when the distribution (in R^d)
    %       of X and Xref are quite similar.
    %
    %   Reference:
    %       a better (but slower) approximation can be found in:
    %       "Wasserstein barycenter and its application to texture mixing",
    %       by J. Rabin, G. Peyre, J. Delon, and M. Bernot.

        % Default arguments
        if (~exist('interpolation', 'var'))
            interpolation = '';
        end;

        % Extract the PCA basis
        [base, diagmat] = eig(Xref * Xref');
        if (diagmat(3,3) > 1000 * diagmat(2,2))     % if monochrome
            [base, o] = eig(X * X');
        end;

        % Perform transfer in PCA basis
        Y    = base' * X;
        Yref = base' * Xref;
        [Yeq, E] = Histogram.p_transferPoints1d(Y, Yref, interpolation);
        Xeq = base * Yeq;
    end;

end;

end
