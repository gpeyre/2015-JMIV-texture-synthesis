classdef Sparsity
% SPARSITY  - provide tools for sparse decompositions.
%
%   See also:
%       Sparsity.demo
%       Sparsity.factorize
%       Sparsity.decompose
%
%   Contact:
%       Guillaume.Tartavel @ telecom-paristech.fr
%       Downloadable on perso.enst.fr/~tartavel/


methods (Static)

    function b = checkToolbox()
    % CHECK TOOLBOX     - check wether the SPAMS toolbox is installed.
    %
    %   The SPAMS toolbox provides a better learning algorithm.
    %   See http://spams-devel.gforge.inria.fr/
    %
    %   Usage:
    %       Sparsity.checkToolbox;
    %       b = Sparsity.checkToolbox;
    %
    %   Return:
    %           - display a warning if the toolbox is not found.
    %       b   - true if the toolbox is found.
    %
    %   See also:
    %       Sparsity.factorize

        b = logical(exist('mexTrainDL', 'file'));
        if (~b)
            warning('TOOLBOX:missing', 'SPAMS toolbox not found, sparse factorization will be slower.');
        end;
    end;

    function [D, W, e] = factorize(X, N, L0, niter)
    % FACTORIZE     - factorize a set of vector into a dicitonary.
    %
    %   Usage:
    %       D = Sparsity.factorize(X, N, L0, niter);
    %       [D, W] = Sparsity.factorize(...);
    %
    %   Arguments:
    %       X  (dxK)    - matrix to be approximated
    %       N           - size of the dictionary to be learned
    %       L0          - max number of non-zeros per column in W
    %       niter       - number of iteration of the algorithm
    %
    %   Return:
    %       D, W    - sparse factorization, approximated:
    %                   min || P - DW ||^2
    %                 among D with N rows
    %                   and W s.t. #{i : W(i,k) ~= 0} <= L0(k)
    %
    %   Reference:
    %       alternate minimization using KSVD, described in:
    %   	"Image denoising via learned dict. and sparse representation",
    %       by M. Elad and M. Aharon.

        % Parameters
        [d, K] = size(X);

        % Use SPAMS toolbox if available
        if (exist('mexTrainDL', 'file'))
            args = struct('mode', 3, 'lambda', L0, 'K', N, 'iter', K*niter/512, 'verbose', 0);
            D = mexTrainDL(X, args);
            D = bsxfun(@rdivide, D, sqrt(sum(D.^2, 1)));
            W = mexOMP(X, D, struct('L', L0, 'eps', 0));
            e = [];
            return;
        end;

        % Select N random patches for D
        I0 = randperm(K);
        D = X(:,I0(1:N));
        D = bsxfun(@rdivide, D, sqrt(sum(D.^2, 1)));

        % Monitor the errors
        getError2 = @(Y) sum(Y(:).^2);
        e = zeros(2 * niter + 1, 1);
        e(1) = getError2(X);

        % Alternate MP and KSVD
        h = waitbar(0, 'Dictionary: 0%');
        for it = 1:niter
            W = Sparsity.decompose(X, D, L0);
            e(2 * it) = getError2(X - D * W);
            [D, W] = Sparsity.p_KSVD(X, D, W);
            e(2 * it + 1) = getError2(X - D * W);
            waitbar(it/niter, h, sprintf('Dictionary: %.0f%%', 100*it/niter));
        end;
        close(h);
        W = Sparsity.decompose(X, D, L0);
    end;

    function [W, E] = decompose(X, D, L0, L0t)
    % DECOMPOSE     - decompose a set of vectors into a dictionary.
    %
    %   Usage:
    %       W = Sparsity.decompose(X, D, L0);
    %       W = Sparsity.decompose(X, D, [], L0t);
    %       W = Sparsity.decompose(X, D, L0, L0t);
    %       [W, E] = ...
    %
    %   Arguments:
    %       X   (dxK)   - matrix to be approximated
    %       D   (dxN)   - dictionary, used to approximate X as DW
    %       L0  (1|K)   - max number of non-zeros per column in W
    %       L0t (1|N)   - max number of non-zeros per row in W
    %
    %   Return:
    %       W   - sparse decomposition, approximated:
    %               W = arg min || P - DW' ||^2
    %             among W s.t. #{i : W(i,k) ~= 0} <= L0(k)
    %                      and #{k : W(i,k) ~= 0} <= L0t(i)
    %       E   - squared error || P - DW ||^2
    %
    %   Reference:
    %       the heuristic (with no L0t) is the MP + back-projection:
    %       "Matching pursuits with time-frequency dictionaries",
    %        by S. Mallat and Z. Zhang.

        % Check normalization of D
        Dnorm2 = sum(D.^2, 1);
        if (max(abs(Dnorm2 - 1)) > 1e-9)
            error('Dictionary D must be normalized');
        end;

        % Decomposition
        iter_factor = 1.5;
        if (~exist('L0t', 'var') || isempty(L0t))
            W = Sparsity.p_MP(X, D, L0, -iter_factor);
        else
            W = Sparsity.p_GMP(X, D, L0, L0t, -iter_factor);
        end;
        W = Sparsity.p_backProjection(X, D, W);

        % Error
        R = X - D * W;
        E = sum(R(:) .^ 2);
    end;

    function d2 = distance2(X, D, W)
    % DISTANCE 2    - compute the approximation error.
    %
    %   Usage:
    %       d2 = Sparsity.distance2(X, D, W);
    %
    %   Return:
    %       || X - DW ||^2, the approximation error
        R = X - D * W;
        d2 = sum(R(:).^2);
    end;

    function demo()
    % DEMO      - launch a demo of the sparsity tools.
    %
    %   Usage:
    %       Sparsity.demo

        % Choose some dimensions and constraints
        d = 64;
        K = 1000;
        N = 100;
        L0 = 8;
        L0t = L0 * K/N;     % uniform

        % Generate a synthetical example
        X = rand(d, K);
        D = rand(d, N);
        D = bsxfun(@rdivide, D, sqrt(sum(D.^2, 1)));

        % Decompose it
        W = Sparsity.decompose(X, D, L0);           % row constrained
        Wt = Sparsity.decompose(X, D, L0, L0t);     % row + col constrained

        % Display
        clf;
        err2 = sqrt(Sparsity.distance2(X, D, W));
        subplot 223; plot(1:N, sum(W~=0, 2), '.b');
        xlabel('Dictionary atoms'); xlim([1 N]);
        ylabel('Number of occurences'); ylim([0 max(sum(W~=0, 2))+.5]);
        subplot 221; hold on;
        plot([1 K], [L0 L0]+.25, '-r', 'LineWidth', 2);
        plot(1:K, sum(W~=0, 1), '.b');
        title(sprintf('L2 error: %.1f', sqrt(err2)));
        xlabel('ID of the approximated vectors'); xlim([1 K]);
        ylabel('Sparsity of the approximation'); ylim([0 L0+.5]);

        err2t = sqrt(Sparsity.distance2(X, D, Wt));
        subplot 224; hold on;
        plot([1 N], [L0t L0t]+.25, '-r', 'LineWidth', 2);
        plot(1:N, sum(Wt~=0, 2), '.b');
        xlabel('Dictionary atoms'); xlim([1 N]);
        ylabel('Number of occurences'); ylim([0 L0t+.5]);
        subplot 222; hold on;
        plot([1 K], [L0 L0]+.25, '-r', 'LineWidth', 2);
        plot(1:K, sum(Wt~=0, 1), '.b');
        title(sprintf('L2 error: %.1f', sqrt(err2t)));
        xlabel('ID of the approximated vectors'); xlim([1 K]);
        ylabel('Sparsity of the approximation'); ylim([0 L0+.5]);
        drawnow;

        % Learn a dictionary on an image
        patchSize = 8;
        N = 192;
        L0 = 4;
        pat = Patches('peppers.png', patchSize);
        im = pat.image;
        [D,W,e] = Sparsity.factorize(pat.asMatrix, N, L0, 30);

        % Reconstruct the image and a mosaic
        pat.asMatrix = D * W;
        dico = Patches([], patchSize, patchSize + 2);
        dico.asMatrix = bsxfun(@rdivide, abs(D), max(abs(D), [], 1));

        % Display it
        psnr = @(X) 20 * log10(sum(abs(X(:))));
        figure;
        subplot 222; imshow(dico.image); title('Dictionary');
        subplot 223; imshow(im); title('Input');
        subplot 224; imshow(pat.image);
        title(sprintf('Approximation, %.2f dB', psnr(im - pat.image)));
        subplot 221; plot(1:length(e), e / K, '.-'); title('Energy');
        xlabel('Iterations'); axis tight;
        ylabel('Squared error per patch'); ylim([0 1.1*max(e)/K]);
    end;

end;


%%%%%%%%%% PRIVATE FUNCTIONS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

methods (Static)

    function [U, VS, s] = p_fastSVD(X, K)
    % Private -- FAST SVD       - fast SVD for large matrices.
    %
    %   Note:
    %       this function is significantly faster than Matlab SVD function
    %       for matrices of size MxN with M << N
    %
    %   Usage:
    %       [U, VS] = Sparsity.p_fastSVD(X);
    %       [U, VS] = Sparsity.p_fastSVD(X, k);
    %       [U, VS, s] = ...
    %
    %   Argument:
    %       X       - a matrix of size MxN, factorized as U * VS'
    %       K       - number of singular values/vectors to be computed
    %                 must be at most M, default is M (exact decomposition)
    %
    %   Return:
    %       U       - matrix of size MxK, satisfying U'  * U   = eye(K)
    %       VS      - matrix of size NxK, satisfying VS' * VS  = diag(s.^2)
    %       s       - vector of K singular values of X, in decreasing order

        % Number of values of the decomposition
        if (~exist('K', 'var') || isempty(K))
            K = size(X, 1);
        elseif (K <= 0 || K > size(X, 1))
            error('invalid argument ''K''');
        end;

        % Get left eigen-vectors
        [U, SS] = eig(X * X');
        s = sqrt(diag(SS));

        % Order it in decreasing order
        k0 = size(X, 1) - K + 1;
        U = U(:,end:-1:k0);
        s = s(end:-1:k0);

        % Deduce VS
        VS = X' * U;
    end

    function [D, W] = p_KSVD(X, D, W)
    % Private -- KSVD   - update a dictionary using K-SVD.
    %
    %   Usage:
    %       [D, W] = Sparsity.p_KSVD(X, D, W);
    %
    %   Arguments:
    %       X  (dxK)    - matrix to be approximated
    %       D  (dxN)    - dictionary, used to approximate X as DW
    %       W  (NxK)    - coefficients of X in D
    %
    %   Return:
    %       D, W        - updated matrices
    %
    %   Reference:
    %       KSVD step presented in:
    %   	"Image denoising via learned dict. and sparse representation",
    %       by M. Elad and M. Aharon.

        % Dimensions
        [d ,K ] = size(X);
        [dd,N]  = size(D);
        [NN,KK] = size(W);
        if (d ~= dd || N ~= NN || K ~= KK)
            error('Invalid arguments, dimension mismatch');
        end;

        % Updates...
        for n = randperm(N)
            ns = ( (1:N) ~= n);
            ks = (W(n,:) ~= 0);

            % Case of an unused atom
            if (~any(ks))
                ks = true(size(ks));
            end;

            % Decompose residual
            R = X(:,ks) - D(:,ns) * W(ns,ks);
            [U,VS] = Sparsity.p_fastSVD(R, 1);

            % Update
            D(:,n) = U(:,1);
            W(n,ks) = VS(:,1);
        end;

    end

    function W = p_MP(X, D, L0, niter)
    % Private -- MP     - sparse decomposition using matching pursuit.
    %
    %   Usage:
    %       W = Sparsity.p_MP(X, D, L0);
    %       W = Sparsity.p_MP(X, D, L0, niter);
    %
    %   Arguments:
    %       X  (dxK)    - matrix to be approximated
    %       D  (dxN)    - dictionary, used to approximate X as DW
    %                     its columns must be normalized with L2 norm = 1
    %       L0 (1|K)    - sparsity constraint on W
    %       niter       - number of iteration, times max(L0) if negative
    %                     default: -1 = max(L0)
    %
    %   Return:
    %       W  (NxK)    - approximated weights:
    %                       W ~ arg min || P - D W' ||
    %                     among W' s.t. #{i : W'(i,k) != 0} <= L0(k)
    %
    %   Reference:
    %       the heuristic used is the Matching Pursuit algorithm:
    %       "Matching pursuits with time-frequency dictionaries",
    %        by S. Mallat and Z. Zhang.
    %
    %   See also:
    %       Sparsity.p_backProjection

        % Dimensions
        [d, N] = size(D);
        [dd,K] = size(X);

        % Optional arguments
        if (~exist('niter', 'var'))
            niter = max(L0);
        elseif (niter < 0)
            niter = ceil(-niter * max(L0));
        end;
        if (isscalar(L0))
            L0 = repmat(L0, [1 K]);
        else
            L0 = L0(:)';
        end;

        % Check dimensions
        if (d ~= dd)
            error('Incompatible dimensions for D and X.');
        elseif (length(L0) ~= K)
            error('Invalid dimension for L0');
        end;

        % Initialization
        W = zeros(N, K);
        L0_count = zeros(size(L0));
        innerProd = D' * X;
        dicoCovar = D' * D;

        % Run MP
        for step = 1:niter

            % Find non saturated constraints
            Istep = find(L0_count < L0);

            % Find best correlation
            [osef, Imax] = max(abs(innerProd(:,Istep)), [], 1);
            Iw = sub2ind([N, K], Imax, Istep);

            % Update L0 count
            L0_count(Istep) = L0_count(Istep) + (W(Iw) == 0);

            % Update weights
            dW = innerProd(Iw);
            W(Iw) = W(Iw) + dW;

            % Update inner products
            innerProd(:,Istep) = innerProd(:,Istep) ...
                - bsxfun(@times, dW, dicoCovar(:,Imax));

            % Stop if convergence
            if (isempty(Istep))
                break;
            end;
        end;
    end;

    function W = p_GMP(X, D, L0, L0t, niter)
    % Private -- GMP     - sparse decomposition.
    %
    %   Usage:
    %       W = Sparsity.p_MP(X, D, L0, L0t);
    %       W = Sparsity.p_MP(X, D, L0, L0t, niter);
    %
    %   Arguments:
    %       X   (dxK)   - matrix to be approximated
    %       D   (dxN)   - dictionary, used to approximate X as DW
    %                     its columns must be normalized with L2 norm = 1
    %       L0  (1|K)   - sparsity constraint on the columns of W
    %       L0t (1|N)   - sparsity constraint on the rows of W
    %       niter       - number of iteration
    %                     if negative, relative to min(sum(L0), sum(L0t))
    %                     default is -1 = min(sum(L0), sum(L0t))
    %
    %   Return:
    %       W  (NxK)    - approximated weights:
    %                       W ~ arg min || P - D W' ||
    %                     among W' s.t. #{i : W'(i,k) != 0} <= L0(k)
    %                     among W' s.t. #{k : W'(i,k) != 0} <= L0t(i)
    %
    %   See also:
    %       Sparsity.p_backProjection

        % Dimensions
        [d, N] = size(D);
        [dd,K] = size(X);

        % Optional arguments
        if (isscalar(L0))
            L0 = repmat(L0, [1 K]);
        else
            L0 = L0(:)';
        end;
        if (isscalar(L0t))
            L0t = repmat(L0t, [N 1]);
        else
            L0t = L0t(:);
        end;
        if (~exist('niter', 'var'))
            niter = min(sum(L0), sum(L0t));
        elseif (niter < 0)
            niter = ceil(-niter * min(sum(L0), sum(L0t)));
        end;

        % Check dimensions
        if (d ~= dd)
            error('Incompatible dimensions for D and X.');
        elseif (length(L0) ~= K)
            error('Invalid dimension for L0');
        elseif (length(L0t) ~= N)
            error('Invalid dimension for L0t');
        end;

        % Test if mex function exist
        if (~exist('mexGMP', 'file'))
            error('Cannot find compiled mex function ''mexGMP''');
        end;
        W = mexGMP(X, D, L0, L0t, niter);

    end;

    function W = p_backProjection(X, D, W)
    % Private -- BACK PROJECTION    - optimize the weights on its support.
    %
    %   Usage:
    %       Wf = Sparsity.p_backProjection(X, D, W);
    %
    %   Argument:
    %       X       - matrix to be approximated
    %       DW      - factorization approximating X
    %
    %   Return:
    %       Wf      - optimized weights:
    %                   Wf = arg min || X - D W' ||_2
    %                 among W' s.t. support(W') = support(W)

        K = size(W, 2);
        for k = 1:K
            I = (W(:,k) ~= 0);
            W(I,k) = D(:,I) \ X(:,k);
        end;
    end;

end;

end
