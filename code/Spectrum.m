classdef Spectrum
% SPECTRUM      - provide functions to manipulate Fourier spectrum.
%
%   See also:
%       Spectrum.demo
%       Spectrum.transfer
%       Spectrum.randomize
%       Spectrum.periodic
%
%   Contact:
%       Guillaume.Tartavel @ telecom-paristech.fr
%       Downloadable on perso.enst.fr/~tartavel/


properties (Constant)
    periodicPreProcessing = true;       % use periodic component for transfer
end;

methods (Static)

    function P = periodic(im)
    % PERIODIC      - perdiodic component of an image.
    %
    %   Usage:
    %       P = Spectrum.periodic(im);
    %
    %   Arguments:
    %       im  - an image
    %
    %   Return:
    %       P   - periodic component of the input image
    %
    %   Reference:
    %       the periodic+smooth decomposition is described in:
    %       "Periodic plus smooth image decomposition",
    %       by Lionel Moisan.

        % Laplacian over the border
        du = zeros(size(im));
        du(1,1:end,:)   = du(1,1:end,:)   + im(end,1:end,:) - im(1,1:end,:);
        du(end,1:end,:) = du(end,1:end,:) + im(1,1:end,:)   - im(end,1:end,:);
        du(1:end,1,:)   = du(1:end,1,:)   + im(1:end,end,:) - im(1:end,1,:);
        du(1:end,end,:) = du(1:end,end,:) + im(1:end,1,:)   - im(1:end,end,:);

        % Fourier processing
        L = Spectrum.p_fourierLaplacian(size(im));
        S = fft2(du) ./ L; % = FFT(smooth(im))
        S(1,1,:) = 0;
        P = im - ifft2(S, 'symmetric');
    end;

    function [imF, E] = transfer(im, imRef)
    % TRANSFER      - transfer the spectrum of an image.
    %
    %   Usage:
    %       imF = Spectrum.transfer(im, imRef);
    %       [imF, E] = ...
    %
    %   Argument:
    %       im      - an image
    %       imRef   - another image
    %
    %   Return:
    %       imF     - image with the same spectrum as 'imRef'
    %       E       - squared distance between the spectra
    %
    %   Note:
    %       the periodic component of 'imRef' is extracted first
    %
    %   Note:
    %       - for gray images with the same size:
    %           imF = arg min || im - imBis ||^2
    %         among all imBis with the same spectrum amplitude as imRef.
    %       - for color images with the same size:
    %         the difference of phases between channels are preserved.
    %       - for images with different sizes:
    %         the image is first cropped or windowed+padded.

        % Compute FFT
        ftIm = fft2(im);
        ftRef = fft2(Spectrum.p_extendImage(imRef, size(im)));

        % Compute new image FT
        innerProd = sum(ftIm .* conj(ftRef), 3);
        dephase = innerProd ./ (abs(innerProd) + eps);
        ftNew = bsxfun(@times, ftRef, dephase);

        % Inverse FFT
        imF = ifft2(ftNew, 'symmetric');
        E = sum((im(:) - imF(:)) .^2);
    end;

    function synth = randomize(im, outSize)
    % RANDOMIZE     - randomize the phases of the image Fourier tansform.
    %
    %   Usage:
    %       synth = Spectrum.randomize(im);
    %       synth = Spectrum.randomize(im, outSize);
    %
    %   Arguments:
    %       im                  - an image
    %       outSize = [M N]     - output size, default is the size of 'im'

        % Default argument
        if (~exist('outSize', 'var'))
            outSize = size(im);
        end;

        % Check arguments
        K = size(im, 3);
        switch (numel(outSize))
            case 1, outSize = [outSize outSize K];
            case 2, outSize = [outSize K];
            case 3, % nothing to be done
            otherwise, error('Invalid argument ''newSize''');
        end;

        % Compute the spectrum
        imPadded = Spectrum.p_extendImage(im, outSize);
        spectrum = fft2(imPadded);

        % Generate random phases (antisymmetric)
        phasesRaw = 2 * rand(outSize(1:2)) - 1;
        phases = phasesRaw - circshift(rot90(phasesRaw, 2), [1 1]);

        % Dephase the spectrum
        ft = bsxfun(@times, spectrum, exp(1i * pi * phases));
        synth = ifft2(ft, 'symmetric');
    end;

    function d2 = distance2(im, imRef, periodize)
    % DISTANCE 2    - compute the distance between the spectrum.
    %
    %   Usage:
    %       d2 = Spectrum.distance2(im, imRef);
    %       d2 = Spectrum.distance2(im, imRef, periodize);
    %
    %   Argument:
    %       im, imRef   - images of the same size
    %       periodize   - if true (default), images are made periodic first
    %
    %   Return:
    %       the (squared) distance between the spectrum
    %
    %   Note:
    %       this function works only for images of the same size
    %
    %   See also:
    %       Spectrum.transfer

        % Check size
        if (any(size(im) ~= size(imRef)))
            error('Spectrum distance only implemented for images of same size');
        end;

        % Periodize
        if (~exist('periodize', 'var') || isempty(periodize) || periodize)
            im = Spectrum.periodic(im);
            imRef = Spectrum.periodic(imRef);
        end;

        % Compute the distance
        ftIm = fft2(im);
        ftRef = fft2(imRef);
        innerProd = sum(ftIm .* conj(ftRef), 3);
        dephase = innerProd ./ (abs(innerProd) + eps);
        ftDiff = ftIm - bsxfun(@times, ftRef, dephase);
        d2 = sum(abs(ftDiff(:)).^2) / length(ftDiff(:));

    end;

    function demo()
    % DEMO      - launch a demo of the spectrum tools.
    %
    %   Usage:
    %       Spectrum.demo

        % Name of the images to be loaded
        imNameA = 'stone.png';
        imNameB = 'peppers.png';

        % Load images an gray versions
        imA = double(imread(imNameA)) / 255;
        imB = double(imread(imNameB)) / 255;


        % Compute periodic component and random phases
        imP = Spectrum.periodic(imA);
        imR = Spectrum.randomize(imA);
        imS = Spectrum.transfer(imB, imP);

        % Some functions
        getSpectrum = @(X) 20 * log10(abs(fftshift(fft2(rgb2gray(X)))));
        shiftImage = @(X) circshift(X, floor([size(X, 1) size(X, 2)] / 2));

        % Display the image
        clf;
        subplot(3,4,1);  imshow(imA); title('Original');
        subplot(3,4,9);  imshow(shiftImage(imA)); title('Circ shift');
        subplot(3,4,5);  imagesc(getSpectrum(imA));
            set(gca, 'XTick', []); set(gca, 'YTick', []);
        subplot(3,4,2);  imshow(imP); title('Periodic component');
        subplot(3,4,10); imshow(shiftImage(imP));
        subplot(3,4,6); imagesc(getSpectrum(imP));
            set(gca, 'XTick', []); set(gca, 'YTick', []);
        subplot(3,4,3);  imshow(imR); title('Random phases');
        subplot(3,4,7); imagesc(getSpectrum(imR));
            set(gca, 'XTick', []); set(gca, 'YTick', []);
        subplot(3,4,4);  imshow(imS); title('Spectrum transfer');
        subplot(3,4,12); imshow(imB); title('Target image');
        subplot(3,4,8); imagesc(getSpectrum(imS));
            set(gca, 'XTick', []); set(gca, 'YTick', []);
    end;

end;


%%%%%%%%%%  PRIVATE FUNCTIONS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

methods (Static)

    function L = p_fourierLaplacian(imSize)
    % Private -- FOURIER LAPLACIAN      - Fourier transf. of the Laplacian.
    %
    %   Usage:
    %       L = p_fourierLaplacian(imSize);
    %
    %   Arguments:
    %       imSize = [M N]      - dimension of an image
    %
    %   Return:
    %       L                   - DFT of the 2D Laplacian operator
    %                             over the set of images of size MxN

        % Check dimensions
        switch (numel(imSize))
            case 1, imSize = [imSize imSize 1];
            case 2, imSize = [imSize 1];
            case 3, % nothing to be done
            otherwise, error('Invalid argument ''imSize''');
        end;

        % Build the matrix
        p = (1:imSize(1))' - 1;     % column vector
        q = (1:imSize(2))  - 1;     % row vector
        cosP = cos(2 * pi * p / imSize(1));
        cosQ = cos(2 * pi * q / imSize(2));
        L1 = bsxfun(@plus, cosP, cosQ) * 2 - 4;
        L = repmat(L1, [1 1 imSize(3)]);
    end;

    function win = p_getWindow1d(N, alpha)
    % Private -- GET WINDOW 1D      - return a 1D window.
    %
    %   Usage:
    %       win = Spectrum.p_getWindow1d(N, alpha);
    %       [win, wsum] = Spectrum.p_getWindow1d(...);
    %       Spectrum.p_getWindow1d(...);
    %
    %   Arguments:
    %       N       - size of the window
    %       alpha   - roll-off parameter
    %
    %   Return:
    %       win     - a window of size N
    %       wsum    - sum of the coefficients of the window
    %
    %   Note:
    %       - the window is a raised cosine.
    %       - if no argument is provided, the window is displayed.

        % Dimension and function
        dN = floor(alpha * N);
        border = .5 - .5 * cos(pi * (1:dN) / (dN + 1));

        % Check size
        if (2 * dN > N)
            error('alpha must be less than 1/2');
        end;

        % Window
        win = ones(N, 1);
        win(1:dN) = border;
        win(end+1 - (1:dN)) = border;

        % Display
        if (nargout == 0)
            Nfft = 4*N;
            ft = 20 * log10(abs(fft(win, 2 * Nfft)));
            ft = ft - max(ft);
            clf;
            subplot 121; plot(1:N, win, '.-');
                title('Window');
            subplot 122; plot(0:Nfft-1, ft(1:Nfft));
                xlim([0 Nfft-1]);
                title('Fourier Transform');
        end;
    end;

    function imF = p_extendImage(im, newSize)
    % Private -- EXTEND IMAGE   - extend the image to match the new size.
    %
    %   Usage:
    %       imF = Spectrum.p_extendImage(im, newSize);
    %
    %   Arguments:
    %       im                  - an image
    %       newSize = [M N]     - new size
    %
    %   Return:
    %       imF     - obtained by cropping / padding 'im' to the new size.
    %
    %   Note:
    %       to avoid spetrum artifact, the image is
    %       - windowed before being zero-padded, or
    %       - made periodic after being cropped if 'periodicPreProcessing' is true.

        % Check dimensions
        [M N K] = size(im);
        switch (numel(newSize))
            case 1, newSize = [newSize newSize K];
            case 2, newSize = [newSize K];
            case 3, % nothing to be done
            otherwise, error('Invalid argument ''newSize''');
        end;

        % Deal with color
        if (newSize(3) == 1 && K ~= 1)
            im = rgb2gray(im);
        elseif (newSize(3) ~= 1 && K == 1)
            im = repmat(im, [1 1 newSize(3)]);
        elseif (newSize(3) ~= K)
            error('Incompatible number of channels');
        end;

        % Function to get first and last indices
        diffSize = [M N] - newSize(1:2);
        first = 1 + max(0, floor(diffSize / 2));
        last = [M N] - max(0, ceil(diffSize / 2));
        imCropped = im(first(1):last(1), first(2):last(2), :);

        % Keep the periodic component
        if (newSize(1) <= M && newSize(2) <= N && Spectrum.periodicPreProcessing)
            imF = Spectrum.periodic(imCropped);
        else

            % Compute a window (with L2 normalization)
            alpha = .1;
            winA = Spectrum.p_getWindow1d(size(imCropped, 1), alpha);
            winB = Spectrum.p_getWindow1d(size(imCropped, 2), alpha);
            win = bsxfun(@times, winA, winB');

            % Apply it (preserve mean and variance)
            meanValue = mean(mean(imCropped, 1), 2);
            imCentered = bsxfun(@minus, imCropped, meanValue);
            ratio = newSize(1) * newSize(2) / sum(winA.^2) / sum(winB.^2);
            imWinCentered = sqrt(ratio) * bsxfun(@times, imCentered, win);
            imWin = bsxfun(@plus, imWinCentered, meanValue);

            % Zero-pad it
            imF = repmat(meanValue, newSize(1:2));
            imF(1:size(win,1), 1:size(win,2), :) = imWin;
        end;

    end;

end;

end
