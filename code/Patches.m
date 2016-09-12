classdef Patches
% PATCHES       - interface to handle the patches of an image.
%
%   This class is an interface to access the patches of an image.
%
%   To create the object:
%   >> obj = Patches(im, 12);               % patches of size 12
%   >> obj = Patches('lena.png', 8);        % load the image
%   >> obj = Patches(im, 8, 8);             % only distinct patches
%
%   To view the image:
%   >> imshow(obj.image);                   % view the image
%   >> obj.view();                          % same, with some options
%
%   To access a patches:
%   >> p   = obj.asMatrix(:,1);             % as a vector
%   >> pim = obj.asImages(:,:,:,1,1);       % as an image
%
%   Example of patch processing:
%   >> P = obj.asMatrix;                    % get patches as vectors
%   >> Q = median(P, 1);                    % compute the median value
%   >> obj.asMatrix = Q;                    % put them back into the object
%   >> imshow(obj.image);                   % view the resulting image
%
%   See also:
%       Patches         - the constructor
%       view            - to view the image
%       asImages        - to extract and update the patches
%                         dimension: Nyi x Nzi x Nc x Nyp x Nzp
%       asMatrix        - to extract and update the patches
%                         dimension: Nyi*Nzi*Nc x Nyp*Nzp
%       Ny,  Nz,  Nc    - dimensions of the image and its weights
%       Nyi, Nzi, Nc    - dimension of each patch
%       Nyp, Nzp        - number of patches in the image
%
%   Contact:
%       Guillaume.Tartavel @ telecom-paristech.fr
%       Downloadable on perso.enst.fr/~tartavel/


properties (Dependent)
    isPeriodic;         % if true, use periodic boundary conditions
    image;              % reconstructed image
    weights;            % weights in the reconstructed image
    asImages;           % patches values (as small images)
    asMatrix;           % patches values (as a matrix)
    indicesAsMatrix;    % patches indices (as a matrix)
    indicesAsImages;    % patches indices (as small images)
    Nc;                 % number of channels
    Ny;                 % height of the image
    Nz;                 % width of the image
    Nyi;                % height of a patch
    Nzi;                % width of a patch
    Nyp;                % number of patches, vertically
    Nzp;                % number of patches, horizontally
    Dp;                 % step between adjacent patches
end;


methods (Access=public)

    function this = Patches(image, psize, pspace, periodic)
    % PATCHES   - constructor.
    %
    %   Usage:
    %       obj = Patches(image, psize);
    %       obj = Patches(image, psize, pspace);
    %       obj = Patches(image, psize, pspace, periodic);
    %       obj = Patches([], matrix);
    %       obj = Patches([], matrix, border);
    %       obj = Patches([], matrix, border, normalize);
    %
    %   Arguments:
    %       image       - image array or image path
    %       psize       - size of the patches (e.g. 12 for 12x12 patches)
    %       pspace      - step between adjacent patches
    %                     default: 1 (overlapping patches)
    %                     note: 0 means pspace = psize (tiling of patches)
    %       periodic    - if true, the image is assumed to be periodic
    %                     default: false
    %       matrix      - a matrix of patches, converted into a mosaic
    %       border      - border around each patch of the mosaic
    %                     default: 1
    %       normalize   - if true, patches are rescaled to fit in -1..1
    %                     default: false

        if (~isempty(image))

            % Load the image
            if (ischar(image))
                this.p_image = double(imread(image)) / 255;
            else
                this.p_image = image;
            end;

            % Set p_delta
            if (~exist('pspace', 'var'))
                this.p_delta = 1;
            elseif (pspace == 0)
                this.p_delta = psize;
            else
                this.p_delta = pspace;
            end;

            % Initialize all other properties
            this.p_periodic = exist('periodic', 'var') && logical(periodic);
            [this.p_indices, this.p_image] = this.p_buildIndices(psize);
            % ... p_weights could be computed later
            this.p_weights = this.p_flatten(ones(psize^2*this.Nc, this.Nyp*this.Nzp));

        else

            % Optional attributes
            normalize = exist('periodic', 'var') && logical(periodic);
            if (~exist('pspace', 'var'))
                pspace = 1;
            end;

            % Get the patch dimensions
            asMatrix = psize;
            [d,Np] = size(asMatrix);
            [psize,Nc] = Patches.p_matrixToDimensions(asMatrix);

            % Normalize
            if (normalize)
                norms = max(abs(asMatrix), [], 1);
                norms(norms == 0) = 1;
                asMatrix = bsxfun(@rdivide, asMatrix, norms);
            end;

            % Find a good ratio between 1:1 and 1:2
            Nzp = ceil(sqrt(Np)) : floor(sqrt(2*Np));
            Nyp = ceil(Np ./ Nzp);
            Ndiff = Nyp .* Nzp - Np;
            [o,I] = min(Ndiff);

            % Create the image
            P = [asMatrix zeros(d, Ndiff(I))];
            imsize = [Nyp(I) Nzp(I)] * (psize + pspace) - pspace;
            im = zeros([imsize Nc]);
            this = Patches(im, psize, psize + pspace);
            this.image = this.p_flatten(P);

        end;

    end;

    function view(this, ax, clickFcn)
    % VIEW      - interactive view of the patches.
    %
    %   Usage:
    %       obj.view();
    %       obj.view(ax);
    %       obj.view(ax, true);         % display the clicked patches
    %       obj.view(ax, clickFcn);     % custom function on patch click
    %
    %   Arguments:
    %       ax          - output axes handle
    %                     default: gca
    %       chickFcn    - function of YZ = [yp, zp] called on click
    %                     default: []
    %                     note: if true, display the clicked patches

        % Draw in the axes
        if (~exist('ax', 'var') || isempty(ax))
            ax = gca;
        end;
        imshow(this.image, 'Parent', ax);
        title(ax, 'Image');
        drawings = get(ax, 'Children');

        % Click processing
        if (exist('clickFcn', 'var') && ~isempty(clickFcn))
            offset = (this.Nyi + 1) / 2;
            yzMax = [this.Nyp, this.Nzp];

            % Functions
            clickXY = @(h) get(get(h, 'Parent'), 'CurrentPoint');
            convertXY = @(XY) (XY(1,[2,1]) - offset) / this.Dp + 1;
            formatYZ = @(YZ) max(1, min(yzMax, round(YZ)));
            clickYZ = @(obj) formatYZ(convertXY(clickXY(obj)));

            % Default function
            if (islogical(clickFcn))
                figure;
                axpatch = gca;
                figure(get(ax, 'Parent'));
                clickFcn = @(X) imshow(...
                    imresize(this.asImages(:,:,:,X(1),X(2)), 4), ...
                    'Parent', axpatch);
                clickFcn([1 1]);
            end;

            % Set callback
            set(drawings(end), 'ButtonDownFcn', @(obj, evt) clickFcn(clickYZ(obj)));
        end;
    end;

end;


methods % getters, setters

    function out = get.isPeriodic(this)
        out = this.p_periodic;
    end;
    function out = get.image(this)
        out = this.p_image;
    end;
    function out = get.weights(this)
        out = this.p_weights;
    end;
    function out = get.asImages(this)
        out = this.p_image(this.indicesAsImages);
    end;
    function out = get.asMatrix(this)
        out = Patches.p_imagesToMatrix(this.p_image(this.p_indices));
    end;
    function out = get.indicesAsImages(this)
        out = this.p_indices;
    end;
    function out = get.indicesAsMatrix(this)
        out = Patches.p_imagesToMatrix(this.p_indices);
    end;
    function out = get.Nc(this)
        out = size(this.p_indices, 3);
    end;
    function out = get.Ny(this)
        out = size(this.p_image, 1);
    end;
    function out = get.Nz(this)
        out = size(this.p_image, 2);
    end;
    function out = get.Nyi(this)
        out = size(this.p_indices, 1);
    end;
    function out = get.Nzi(this)
        out = size(this.p_indices, 2);
    end;
    function out = get.Nyp(this)
        out = size(this.p_indices, 4);
    end;
    function out = get.Nzp(this)
        out = size(this.p_indices, 5);
    end;
    function out = get.Dp(this)
        out = this.p_delta;
    end;

    function this = set.image(this, image)
        if (isequal(size(this.image), size(image)))
            this.p_image = image;
        else
            error('image dimension cannot be changed');
        end;
    end;

    function this = set.asImages(this, asImages)
        this = this.p_reconstruct(Patches.p_imagesToMatrix(asImages));
    end;

    function this = set.asMatrix(this, asMatrix)
        this = this.p_reconstruct(asMatrix);
    end;

end;


%%%%%%%%%% PRIVATE SECTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

methods (Static)

    function asMatrix = p_imagesToMatrix(asImages)
    % Private -- IMAGES TO MATRIX  - patches from images shape to matrix.
    %
    %   Usage:
    %       asMatrix = Patches.imagesToMatrix(asImages);
    %
    %   Input:
    %       asImages    - 5D array of patches
    %                     size: Nyi x Nzi x Nc x Nyp x Nzp
    %
    %   Return:
    %       asMatrix    - matrix of patches
    %                     size: Nyi*Nzi*Nc, Nyp*Nzp

        if (ndims(asImages) > 2)
            [Nyi,Nzi,Nc,Nyp,Nzp] = size(asImages);
            asMatrix = reshape(asImages, [Nyi*Nzi*Nc, Nyp*Nzp]);
        else
            error('expected a 5D array');
        end;
    end;

    function [psize, Nc] = p_matrixToDimensions(asMatrix)
    % Private -- MATRIX TO DIMENSION    - find patch dimensions.
    %
    %   Usage:
    %       [psize, Nc] = Patches.matrixToDimensions(asMatrix);

        % Find Nc
        [d, Np] = size(asMatrix);
        c = [1 3];                      % candidates for Nc
        n = round(sqrt(d ./ c));        % candidates for psize
        I = find(n.^2 .* c == d, 1);    % valid candidate

        % Validate
        if (~isempty(I))
            [psize, Nc] = deal(n(I), c(I));
        else
            error('invalid dimension of the patches');
        end;
    end;

end;


properties (Access=private)
    p_periodic;     % boolean
    p_image;
    p_weights;
    p_indices;      % as images
    p_delta;        % grid step
end;


methods (Access=private)

    function [I, image] = p_buildIndices(this, psize)
    % Private -- BUILD INDICES  - build the indice array.
    %
    %   Note:
    %       if the tiling is not perfect, the image is cropped.

        % Number of patches
        Dp = this.p_delta;
        [Ny,Nz,Nc] = size(this.p_image);
        Nb = double(~this.p_periodic) * (psize - Dp);   % offset
        Np = floor(([Ny Nz] - Nb) / Dp);                % number of patches
        Nim = Nb + Np * Dp;                             % size of the image

        % Clamp the image
        o = ([Ny Nz] - Nim) / 2;
        oA = 1 + floor(o);
        oB = ceil(o);
        image = this.p_image(oA(1):end-oB(1), oA(2):end-oB(2), :);

        % Indices of the patches
        [Ny,Nz,Nc] = size(image);
        mkgrid = @(n, step, dim) ...
            reshape((0:n-1) * step, circshift([n 1 1 1 1], [1 dim-1]));
        % ... indices inside a patch
        Ic = 1 + mkgrid(Nc, Ny*Nz, 3);
        Iyi = mkgrid(psize, 1, 1);
        Izi = mkgrid(psize, Ny, 2);
        % ... indices on the grid
        Iyp = mkgrid(Np(1), this.p_delta, 4);
        Izp = mkgrid(Np(2), this.p_delta * Ny, 5);
        % ... merge all
        Iy = mod(bsxfun(@plus, Iyi, Iyp), Ny);
        Iz = mod(bsxfun(@plus, Izi, Izp), Ny*Nz);
        Iyc = bsxfun(@plus, Iy, Ic);
        I = bsxfun(@plus, Iyc, Iz);
    end;

    function flat = p_flatten(this, asMatrix)
    % Private -- FLATTEN    - flatten a patch matrix.
    %
    %   Note:
    %       * the sum of overlapping pixels is computed.
    %       * the size of the patches mustn't change.

        % Check dimensions
        Np = this.Nyp * this.Nzp;
        Ni = this.Nyi * this.Nzi;
        d = Ni * this.Nc;
        if (~isequal(size(asMatrix), [d Np]))
            error('patch matrix has invalid dimensions')
        end;

        % Merge the patches
        flat = zeros(size(this.image));
        I = reshape(this.indicesAsMatrix, Ni, []);
        M = reshape(asMatrix, Ni, []);
        for k = 1:Ni
            Ik = I(k,:);
            flat(Ik) = flat(Ik) + M(k,:);
        end;
    end;

    function obj = p_reconstruct(this, asMatrix)
    % Private -- RECONSTRUCT    - reconstruct the image from its patches.
    %
    %   Note:
    %       * weights are updated if needed.
    %       * the size of the image may change (if patch size changed).

        % Get dimensions
        if (this.Nyp * this.Nzp ~= size(asMatrix, 2))
            error('the number of patches must not change');
        end;
        Np = [this.Nyp this.Nzp];
        [psize,Nc] = Patches.p_matrixToDimensions(asMatrix);

        % New spacing between patches
        prevsize = this.Nyi;
        if (psize > prevsize)
            Dp = this.p_delta + psize - prevsize;       % stretching
        elseif (psize < this.p_delta)
            Dp = this.p_delta + prevsize - psize;       % preserve tiling
        else
            Dp = this.p_delta;                          % preserve spacing
        end;

        % New object if different size
        if (size(asMatrix, 1) == this.Nyi * this.Nzi * this.Nc)
            obj = this;
        else
            imsize = Np * Dp + (psize - Dp) * double(~this.isPeriodic);
            im = zeros([imsize Nc]);
            obj = Patches(im, psize, Dp, this.isPeriodic);
        end;

        % Merge all
        if (isempty(obj.p_weights))
            obj.p_weights = obj.p_flatten(ones(size(asMatrix)));
        end;
        obj.p_image = obj.p_flatten(asMatrix) ./ obj.weights;
        obj.p_image(obj.weights == 0) = 0;
    end;

end;

end
