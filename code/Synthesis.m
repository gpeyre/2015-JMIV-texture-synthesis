classdef Synthesis
% SYNTHESIS 	- variational sparse texture synthesis algorithm.
%
%   Reference:
%       "Variational Sparse Texture Synthesis",
%       by G. Tartavel, Y. Gousseau, and G. Peyré.
%
%   See also:
%       Synthesis.demo
%       Synthesis.dictionary
%       Synthesis.launch
%
%   Contact:
%       Guillaume.Tartavel @ telecom-paristech.fr
%       Downloadable on perso.enst.fr/~tartavel/


properties (Constant)
    viewIterations = true;          % display images during the synthesis
    randomPatchOffsets = true;      % use random offset for the patches
    postProcessing = true;          % post-process using a spectrum transfer
    inputPatchesSubsampling = 1;    % subsampling of the exemplar patches
end;

methods (Static)

    function [D, im] = dictionary(dico, im, nScales, pSize, dicoSize, sparsity, nIter)
    % DICTIONARY    - learn or load a dictionary
    %
    %   Usage (loading):
    %       D = Synthesis.dictionary(dico, im);
    %       D = Synthesis.dictionary(dico, im, nScales);
    %       D = Synthesis.dictionary(dico, im, nScales, ...);
    %       [D,im] = Synthesis.dictionary(...);
    %
    %   Usage (learning):
    %       D = Syn...ry(dico, im, nScales, pSize, dicoSize, sparsity, nIter);
    %       D = Synthesis.dictionary([], im, ...);
    %       [D,im] = Synthesis.dictionary(...);
    %
    %   Arguments:
    %       dico        - dictionary file (.mat)
    %       im          - image or image filename
    %       nScales     - number of scales
    %       pSize       - patch size
    %       dicoSize    - size of the dictionary
    %       sparsity    - sparsity of the dictionary decomposition
    %       nIter       - number of iteration (per scale) for learning
    %
    %   Return:
    %       D           - dictionary D(:,:,s) at each scale s
    %       im          - the loaded image
    %
    %   Loading / learning:
    %       if dico is a file: it is loaded
    %       up to 3 arguments: dico is loaded
    %       above 3 arguments: the dictionary is learned and saved as dico
    %
    %   See also:
    %       Synthesis.launch

        % Load the image
        if (ischar(im))
            im = double(imread(im)) / 255;
        end;

        % Load or learn the dico
        s = Synthesis.p_signature(im);
        if (nargin <= 3 || exist(dico, 'file'))
            D = Synthesis.p_dicoLoad(dico, s);
            if (exist('nScales', 'var') && ~isempty(nScales))
                D = D(:,:,1:nScales);
            end;
        else
            Sparsity.checkToolbox;
            D = Synthesis.p_dicoLearn(im, nScales, pSize, dicoSize, sparsity, nIter);
            if (exist('dico', 'var') && ~isempty(dico))
                Synthesis.p_dicoSave(D, dico, s);
            end;
        end;
    end;

    function [synth, E] = launch(D, im, outSize, weights, sparsity, pSpacing, nIter)
    % LAUNCH    - launch the synthesis.
    %
    %   Usage:
    %       synth = Synthesis.launch(D, im, outSize, weights, sparsity, pSpacing, nIter);
    %       [synth, E] = Synthesis.launch(...);
    %
    %   Arguments:
    %       D           - dictionary, learned from the input image
    %       im          - input image
    %       outSize     - size of the synthesized image, as [M N]
    %       weights     - weights of the synthesis
    %                     struct. with fields 'histo', 'spectrum', 'patch'
    %       sparsity    - sparsity of the pathc decomposition
    %       pSpacing    - grid step for patches (1 = all patches)
    %       nIter       - number of iteration per scale
    %
    %   Return:
    %       synth       - the synthesized image
    %       E           - value of the cost function at each iteration/scale
    %
    %   Note:
    %       The energy may not be decreasing due to random patch offset.
    %       Set 'Synthesis.randomPatchOffset' to 'false' to disable them.
    %
    %   See also:
    %       Synthesis.dictionary

        % Initialization
        lastScale = ~isfield(weights, 'subscale');
        weights.subscale = true;
        [pSize,Nc] = Synthesis.p_dicoDims(D);

        % Perform multi-scale synthesis
        if (size(D, 3) == 1)
            tmp = min(numel(outSize), 2);
            init = rand([outSize(1) outSize(tmp) Nc]);
            E = [];
        else
            nextIm = Synthesis.p_downsample2(im);
            nextD = D(:,:,2:end);
            nextSize = ceil(outSize / 2);
            [nextSynth, E] = Synthesis.launch(nextD, nextIm, nextSize, weights, sparsity, pSpacing, nIter);
            init = Synthesis.p_upsample2(nextSynth);
        end;

        % Create initial image
        set(gcf, 'Name', sprintf('Synthesis at scale %d', size(D, 3)));
        pInit = Patches(init, pSize, pSpacing, true);
        F0 = Synthesis.p_getAtomFreqs(im, D(:,:,1), pSize, sparsity);

        % Perform synthesis
        [synth, Ek] = Synthesis.p_synthesize(D(:,:,1), im, pInit, weights, sparsity, F0, nIter);
        Ew = [weights.histo; weights.spectrum; weights.patch];
        E = [Ek*Ew E];

        % Post processing
        if (lastScale)
            if (Synthesis.postProcessing)
                synth = Spectrum.transfer(synth, im);
                synth = Histogram.transfer(synth, im);
            end;
            synth = min(max(synth, 0), 1);
        end;

        % Display the result

    end;

    function demo(imPath)
    % DEMO      - launch a demo of the synthesis.
    %
    %   Usage:
    %       Syntehsis.demo
    %       Synthesis.demo(imPath);
    %       synth = Synthesis.demo(...);
    %
    %   Argument:
    %       imPath  - path of an image of texture
    %
    %   Return:
    %       synth   - a synthesis
    %
    %   Warning:
    %       Dictionary is learn first, this may take some time.
    %       It is faster to save it.
    %
    %   See also:
    %       Synthesis.dictionary
    %       Synthesis.launch


        % Parameters
        outSize = [256 256];   	% Size of the output
        pSize = 12;           	% Size of the patches
        pSpacing = 4;           % Spacing of the patches
        sparsity = 4;           % Sparsity of the patch decomposition
        dicoSize = 256;       	% Size of the dicitonary (if learned)
        weights = struct('histo', 1, 'spectrum', 1, 'patch', 1);

        % Learning and synthesis
        if (~exist('imPath', 'var'))
            [D,im] = Synthesis.dictionary('sand.mat', 'sand.png');
        else
            [D,im] = Synthesis.dictionary([], imPath, 3, pSize, dicoSize, sparsity, 30);
        end;
        [synth, E] = Synthesis.launch(D, im, outSize, weights, sparsity, pSpacing, 30);

        % Result
        clf;
        subplot 221; imshow(im); title('Input image');
        subplot 222; imshow(synth); title('Synthesis');
        subplot 212; plot(E, '.-'); xlim([1 30]); ylim([0 max(E(:))]);
        title('Energy -- for each scale'); xlabel('Iterations');
    end;

end;


%%%%%%%%%% PRIVATE FUNCTIONS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

methods (Static)


    % Dictionary functions

    function s = p_signature(X, N)
    % Private -- SIGNATURE      - compute a signature for any array.
    %
    %   Usage:
    %       s = Synthesis.p_signature(X);
    %       s = Synthesis.p_signature(X, N);
    %
    %   Arguments:
    %       X   - an array
    %       N   - size of the signature, default is 2
    %
    %   Return:
    %       s   - a signature of X: its moments of order 0:N

        if (~exist('N', 'var') || isempty(N))
            N = 2;
        end;

        s = sum(bsxfun(@power, X(:), 0:N), 1);

    end;

    function D = p_dicoLearn(im, nScales, pSize, dicoSize, sparsity, nIter)
    % Private -- DICO LEARN     - learn a dictionary from an image.
    %
    %   Usage:
    %       D = Synthesis.p_dicoLearn(im, nScales, pSize, ...
    %               dicoSize, sparsity, nIter);
    %
    %   See also:
    %       Synthesis.dictionary

        % Learn the dictionary
        sub = Synthesis.inputPatchesSubsampling;
        P = Patches(im, pSize, sub).asMatrix;
        D = Sparsity.factorize(P, dicoSize, sparsity, nIter);

        % Deal with the next scales
        if (nScales > 1)
            imLow = Synthesis.p_downsample2(im);
            Dn = Synthesis.p_dicoLearn(imLow, nScales-1, pSize, dicoSize, sparsity, nIter);
            D = cat(3, D, Dn);
        end;
    end;

    function D = p_dicoLoad(filename, s)
    % Private -- DICO LOAD      - load a dictionary.
    %
    %   Usage:
    %       D = Synthesis.p_dicoLoad(filename);
    %       D = Synthesis.p_dicoLoad(filename, s);
    %
    %   Argument:
    %       filename    - name of the dictionary file
    %       s           - signature, used to check the loaded dico
    %
    %   Return:
    %       D           - the loaded dictionary

        % Check if the file exists
        if (~exist(filename, 'file'))
            error('Dictionary "%s" does not exist.', filename);
        end;

        % Load it
        vars = load(filename);
        if (~isfield(vars, 'D'))
            error('Dictionary "%s" is not valid.', filename);
        else
            D = vars.D;
        end;

        % Check the signature
        if (exist('s', 'var') && isfield(vars, 'signature'))
            N = min(numel(s), numel(vars.signature));
            err = s(1:N) - vars.signature(1:N);
            if (any(abs(err) > 1e-6 * abs(vars.signature)))
                [o,name,ext] = fileparts(filename);
                warning('FILE:wrongSignature', ...
                        'Dictionary "%s" has a wrong signature.', strcat(name, ext));
            end;
        end;

        % Normalize the dictionary
        D = bsxfun(@rdivide, D, sqrt(sum(D.^2, 1)));
    end

    function p_dicoSave(D, filename, s)
    % Private -- DICO SAVE      - save a dictionary
    %
    %   Usage:
    %       p_dicoSave(D, filename);
    %       p_dicoSave(D, filename, s);
    %
    %   Arguments:
    %       D           - the dictionary (array)
    %       filename    - output filename
    %       s           - signature

        vars = struct('D', D);
        if (exist('s', 'var'))
            vars.signature = s;
        end;
        save(filename, '-struct', 'vars');
    end;

    function [pSize, Nc, Natoms, Nscales] = p_dicoDims(D)
    % Private -- DICO DIMS  - get the dimensions of the dictionary.
    %
    %   Usage:
    %       [pSize, Nc, Natoms, Nscales] = Synthesis.p_dicoDims(D);
    %
    %   Arguments:
    %       D       - a dictionary
    %
    %   Return:
    %       pSize   - patch size
    %       Nc      - number of channels (1 = gray, 3 = colors)
    %       Natoms  - number of atoms in the dictionary
    %       Nscales - number of scales, for a multi-scale dictionary

        [o,Natoms,Nscales] = size(D);
        pDico = Patches([], D(:,1:2,1));
        pSize = pDico.Nyi;
        Nc = pDico.Nc;

    end;


    % Image functions

    function small = p_downsample2(im)
    % Private -- downsample an image.
    %
    %   Bilinear downsampling by a factor 2.
    %
    %   Usage:
    %       small = Synthesis.p_downsample2(im);

        aux = im(1:2:end-1,:,:) + im(2:2:end,:,:);
        small = (aux(:,1:2:end-1,:) + aux(:,2:2:end,:)) / 4;

    end;

    function big = p_upsample2(im)
    % Private -- upsample an image.
    %
    %   Bicubic periodic upsampling by a factor 2.
    %
    %   Usage:
    %       big = Synthesis.p_upsample2(im);

        [M,N,K] = size(im);
        cconvn = @(x, h) convn(padarray(x, size(h)-1, 'circular', 'pre'), h, 'valid');

        hr = [-1; 9; 9; -1] / 16;
        hl = [0; 1; 0; 0];

        aux = zeros(2*M, N, K);
        aux(1:2:end,:,:) = cconvn(im, hr);
        aux(2:2:end,:,:) = cconvn(im, hl);

        raw = zeros(2*M, 2*N, K);
        raw(:,1:2:end,:) = cconvn(aux, hr');
        raw(:,2:2:end,:) = cconvn(aux, hl');

        big = circshift(raw, [-3 -3]);
    end;

    function [imf, E, patches] = p_sparsifyPatches(im, patches, D, S, F)
    % Private -- SPARSIFY PATCHES   - approximate the image using sparsity.
    %
    %   Usage:
    %       imf = Synthesis.p_sparsifyPatches(im, patches, D, S, F)
    %       [imf, E] = ...
    %
    %   Arguments:
    %       im      - current image u
    %       patches - Patch object with the right dimensions
    %       D       - dictionary
    %       S       - sparsity
    %       F       - frequencies of atoms
    %
    %   Return:
    %       imf     - image u - P*(P(u) - DW) / Z
    %       E       - energy || P - DW || ^2 / Z    -- with the old P.

        % Get the patches and their weights
        patches.image = im;
        P = patches.asMatrix;
        Z = ceil(patches.Nyi / patches.Dp) ^ 2;
        [W, EZ] = Sparsity.decompose(P, D, S, F);

        % Reconstruct the image
        patches.asMatrix = P - D * W;
        imf = im - patches.image .* patches.weights / Z;
        E = EZ / Z;
    end;


    % Synthesis function

    function F0 = p_getAtomFreqs(im, D, pSize, sparsity)
    % Private -- GET ATOM FREQS     - get the atom frequencies for the image.
    %
    %   Usage:
    %       F0 = Synthesis.p_getAtomFreqs(im, D, pSize, sparsity);

        sub = Synthesis.inputPatchesSubsampling;
        pIm = Patches(im, pSize, sub);
        P = pIm.asMatrix;
        W0 = Sparsity.decompose(P, D, sparsity);
        F0 = sum(W0 ~= 0, 2);

    end;

    function [synth, E] = p_synthesize(D, im, pSynth, weights, S, F0, nIter)
    % Private -- SYNTHESIZE     - perform the synthesis of a scale.
    %
    %   Usage:
    %       synth = Synthesis.p_synthesize(D, im, weights, pSynth, S, F0, nIter);
    %       [synth, E] = ...
    %
    %   Arguments:
    %       D           - dictionary
    %       im          - input image
    %       pSynth      - current synthesis as a Patches object
    %       S           - sparsity of the pathc decomposition
    %       F0          - frequencies of each atoms of D
    %       nIter       - number of iteration
    %       weights     - weights of each term
    %           .histo
    %           .spectrum
    %           .patch
    %
    %   Return:
    %       synth       - the synthesized image
    %       E           - squared error before each iteration n,
    %                       E(n,:) = [E_hist, E_spec, E_patch]
    %
    %   See also:
    %       Synthesis.launch

        % Parameters
        wHis = weights.histo;
        wSpe = weights.spectrum;
        wPat = weights.patch;
        wSum = wHis + wSpe + wPat;
        F = pSynth.Nyp * pSynth.Nzp * S  * F0 / sum(F0);

        % Launch the synthesis
        E = zeros(nIter, 3);
        synth = pSynth.image;
        for n = 1:nIter

            % Random offset
            offset = randi(pSynth.Dp, [2 1]) * Synthesis.randomPatchOffsets;
            synth = circshift(synth, offset);

            % Sparse approximation, histogram and spectrum transfer
            [imPatch, E(n,3)] = Synthesis.p_sparsifyPatches(synth, pSynth, D, S, F);
            [imSpec, E(n,2)] = Spectrum.transfer(synth, im);  % = F*(|Fu|/|Fim| Fu)
            [imHist, E(n,1)] = Histogram.transfer(synth, im); % = im o sigma

            % Go to next step
            synth = (wPat * imPatch + wSpe * imSpec + wHis * imHist) / wSum;
            synth = circshift(synth, -offset);

            % Display the images while synthesizing
            if (Synthesis.viewIterations)
                v = @(x) imshow(min(max(x, 0), 1));
                subplot 221; v(synth);   title(sprintf('Step %d/%d', n, nIter));
                subplot 222; v(imPatch); title('Sparse coded');
                subplot 223; v(imSpec);  title('Spectrum transfered');
                subplot 224; v(imHist);  title('Histogram transfered');
                drawnow;
            end;
       end;
    end;

end;

end
