% LAUNCH COPYMAP    - display statistics about patch nearest neighbors.
%
%   For each patch of each synthesized image (from different algorithms),
%   we find the nearest patch in the reference image.
%   We then display:
%   * a histogram of the distances between each patch and its nearest neighbor.
%   * a map showing the location of the nearest neighbor of each patch.
%   * a graph showing the percentage Y of pixels belonging to a region of
%       at least X pixels which has been verbatim-copied from the input.
%
%   You can set the parameters in the first part of this file.
%
%   Contact:
%       Guillaume.Tartavel @ telecom-paristech.fr
%       Downloadable on perso.enst.fr/~tartavel/

clf;
clear;
close all;


%% Parameters -- edit them!
patchSize = 3;
imOur = Patches('copymap/our/Sand-06.png', patchSize);
imRpn = Patches('copymap/RPN/Sand-06.png', patchSize);
imEfr = Patches('copymap/efros/Sand-06.png', patchSize);
imSim = Patches('copymap/sim/Sand-06.png', patchSize);
imRef = Patches('copymap/input/Sand-06.png', patchSize);


% Compute distances
[iEfr,distEfr] = Similar.getNN(imEfr.asMatrix, imRef.asMatrix, 1024);
[iOur,distOur] = Similar.getNN(imOur.asMatrix, imRef.asMatrix, 1024);
[iRpn,distRpn] = Similar.getNN(imRpn.asMatrix, imRef.asMatrix, 1024);
[iSim,distSim] = Similar.getNN(imSim.asMatrix, imRef.asMatrix, 1024);


%% Compute reference map
[o,cmap] = Similar.buildCoordMap(1, [imRef.Nyp imRef.Nzp]);

% Reshape indices as images
iiEfr = reshape(iEfr, imEfr.Nyp, imEfr.Nzp);
iiOur = reshape(iOur, imOur.Nyp, imOur.Nzp);
iiRpn = reshape(iRpn, imRpn.Nyp, imRpn.Nzp);
iiSim = reshape(iSim, imSim.Nyp, imSim.Nzp);

% Get areas of copy
aEfr = Similar.buildCopyMap(iiEfr, size(cmap));
aOur = Similar.buildCopyMap(iiOur, size(cmap));
aRpn = Similar.buildCopyMap(iiRpn, size(cmap));
aSim = Similar.buildCopyMap(iiSim, size(cmap));


%% View images
subplot 451; imshow(imRef.image); title('Input (exemplar)');
subplot 452; imshow(imEfr.image); title('Efros synthesis');
subplot 453; imshow(imOur.image); title('Our synthesis');
subplot 454; imshow(imRpn.image); title('Random Phase Noise');
subplot 455; imshow(imSim.image); title('Portilla/Simoncelli');

% View distances histograms
% X: distance to the nearest neighbor
% Y: amount of patches
bins = (0:.01:1) * max(sqrt(distOur));
subplot(4,5,17); Similar.hist(sqrt(distEfr), bins);
subplot(4,5,18); Similar.hist(sqrt(distOur), bins);
subplot(4,5,19); Similar.hist(sqrt(distRpn), bins);
subplot(4,5,20); Similar.hist(sqrt(distSim), bins);

% View coordinate map
subplot 456; imshow(cmap);
subplot 457; imshow(Similar.buildCoordMap(iiEfr, cmap));
subplot 458; imshow(Similar.buildCoordMap(iiOur, cmap));
subplot 459; imshow(Similar.buildCoordMap(iiRpn, cmap));
subplot(4,5,10); imshow(Similar.buildCoordMap(iiSim, cmap));
subplot 457; title('Location of NN');

% View areas of copy
% X: size of the regions
% Y: prct of pixels belonging to regions of size X+ verbatim-copied from input
icumhist = @(x) flipud(cumsum(flipud(accumarray(x(:), 1))));
amax = max(max(aOur(:)), max(aRpn(:)));
subplot(4,5,12); plot(icumhist(aEfr) / numel(aEfr)); xlim([1 amax]); ylim([0 1]);
subplot(4,5,13); plot(icumhist(aOur) / numel(aOur)); xlim([1 amax]); ylim([0 1]);
subplot(4,5,14); plot(icumhist(aRpn) / numel(aRpn)); xlim([1 amax]); ylim([0 1]);
subplot(4,5,15); plot(icumhist(aSim) / numel(aSim)); xlim([1 amax]); ylim([0 1]);
subplot(4,5,12); title('Amount of pixels in areas bigger than...');
