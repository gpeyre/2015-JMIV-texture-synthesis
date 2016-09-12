classdef Similar
% Similar   -- namespace for patch similarity analysis.

methods (Static)

    function [I, dist2] = getNN(Y, Xref, batchSize)
        % GET NN        -- get the NN in Xref of each column of Y.
        %
        %   Usage:
        %       [I, dist2] = Similar.getNN(Y, Xref);
        %       [I, dist2] = Similar.getNN(Y, Xref, batchSize);
        %
        %   Arguments:
        %       Y     (dxN)     -- N elements of size d
        %       Xref  (dxK)     -- K elements of size d
        %       batchSize       -- split Y into batches (lower memory req.)
        %       I     (1xN)     -- col. indice of the NN of each Y in Xref
        %       dist2  (1xN)    -- dist^2 between cols. of Y and X(:,I)

        % Compute all distances at once.
        if (~exist('batchSize', 'var'))
                D2 = bsxfun(@plus, sum(Y.^2,1), sum(Xref.^2,1)') ...
                    - 2 * Xref'*Y;
                [dist2, I] = min(D2, [], 1);
            return;
        end

        % Initialize
        N = size(Y, 2);
        I = zeros(1, N);
        dist2 = zeros(1, N);

        % Compute indices
        jstart = 1:batchSize:N;
        jend = jstart + batchSize - 1;
        jend(end) = N;

        % Batch processing
        for k = 1:length(jstart)
            [ja, jb] = deal(jstart(k), jend(k));
            [I(ja:jb), dist2(ja:jb)] = Similar.getNN(Y(:,ja:jb), Xref);
        end;
    end;

    function [cmap, cref] = buildCoordMap(imNN, cref)
    % BUILD COORD MAP   -- build a coordinate map from indices
    %
    %   Usage:
    %       cmap = buildCoordMap(imNN, cref);
    %       [cmap, cref] = buildCoordMap(imNN, [M N]);
    %
    %   Arguments:
    %       imNN  (PxQ)     -- indices of the NN
    %       cmap  (PxQx3)   -- coordinate map
    %       cref  (MxNx3)   -- reference coordinate map

        % Generate the reference map
        if (numel(cref) <= 3)
            idx = @(m, n) repmat((0:m-1)' / (m-1), 1, n);
            [M,N] = deal(cref(1), cref(min(2,end)));
            [I,J] = deal(idx(M, N), idx(N, M)');
            cref = cat(3, I, I .* J, J);
        end;

        % Generate the resulting map
        mat = reshape(cref, [], 3);
        cmap = reshape(mat(imNN(:),:), size(imNN, 1), size(imNN, 2), 3);
    end;

    function [imAreas, L] = buildCopyMap(imNN, dimsRef)
    % BUILD COPY MAP    -- build a copy map from indices.
    %
    %   Usage:
    %       imAreas = buildCopyMap(imNN, [M N]);
    %       [imAreas, L] = buildCopyMap(...);
    %
    %   Arguments:
    %       imNN    (PxQ)   -- indices of the NN
    %                MxN    -- dimension of the reference image
    %       L       (PxQ)   -- labeling of copied regions (1..K)
    %       imAreas (PxQ)   -- area of each region

        % Get dimensions
        [M,N] = deal(dimsRef(1), dimsRef(min(2,end)));
        [P,Q] = size(imNN);

        % Build labels
        [I,J] = ind2sub([M N], imNN);
        [I,J] = deal(bsxfun(@minus, I, (1:P)'), bsxfun(@minus, J, 1:Q));
        [I,J] = deal(I - min(I(:)), J -  min(J(:)));
        L = 1 + I + max(I(:)) .* J;

        % Label from 1 to K
        list = zeros(max(L(:)), 1);
        list(L) = 1;
        list(list ~= 0) = 1:sum(list);
        L = list(L);

        % Get area of each region
        stats = regionprops(L, 'Area');
        areas = cat(1, stats.Area);
        imAreas = areas(L);
    end;

    function h = hist(Y, bins)
    % HIST      -- display a normalized histogram
    %
    %   Usage:
    %       h = Similar.hist(Y, bins);
    %
    %   Arguments:
    %       Y       -- 1D observations
    %       bins    -- centers of the histogram's bins

        % Width of each bin
        dhi = diff(bins) / 2;
        dh = dhi([1 1:end]) + dhi([1:end end]);

        % Normalized histogram
        hi = hist(Y, bins);
        hp = hi / sum(hi);      % sum to 1
        h = hp ./ dh;           % area is 1

        % Deal with dirac in 0
        ho = hp(1);
        if (ho > .01)
            h(1) = 0;
        end;

        % Display
        bar(bins, h, 1);
        xlim(bins([1 end]));
        title(sprintf('Copy: %.1f%%', 100 * ho));

    end;

end;

end