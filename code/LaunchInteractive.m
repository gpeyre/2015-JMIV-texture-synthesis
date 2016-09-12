% LAUNCH INTERACTIVE    - interactive script to launch a synthesis.
%
%   Interactive choice of parameters for the synthesis process.
%   To use the default values in brackets, just press 'Enter'.
%
%   See also:
%       LaunchExample
%
%   Contact:
%       Guillaume.Tartavel @ telecom-paristech.fr
%       Downloadable on perso.enst.fr/~tartavel/


%% Initializations

clc;
clear;
close all;

% Display information
disp('----- Interactive settings -----');
disp(' ');

if (~exist('mexTrainDL', 'file'))
    disp('Note: the SPAMS toolbox is not installed.');
    disp('The learning step is more efficient with it.');
    disp(' ');
end;

% Some functions
default = @(x, def) [x, def(1:1-length(x))];
promptPrefix = '>> ';
prompt = @(def) default(input(promptPrefix), def);
promptInt = @(low, def, high) min(max(round(prompt(def)), low), high);


%% Images

% Input
disp('Loading the input image...');
imfilter = {
    '*.png',            'Images as PNG (*.png)'; ...
    '*.jpg; *.jpeg',    'Images as JPEG (*.jpg, *.jpeg)'; ...
    '*.*',              'All Files (*.*)'};
[imname,impath] = uigetfile(imfilter, 'Choice of the Input Image');

% Check it
if (~ischar(imname) || ~ischar(impath))
    error('Interrupted');
else
    impath = fullfile(impath, imname);
end;
disp(['>> ', imname]);
[o,imbasename,imext] = fileparts(imname);


% Output size
disp(' ');
disp('Size of the synthesized image (width/height):');
disp('   def. is [256 256]');
value = input(promptPrefix);
if (isempty(value))
    outSize = [256 256];
else
    outSize = round(value(end:-1:1));
end;


% Output
disp(' ');
disp('Save the synthesized image as...');
disp('   if none, perform only the dictionary learning');
imfilter = {
    '*.png',            'Images as PNG (*.png)'; ...
    '*.jpg; *.jpeg',    'Images as JPEG (*.jpg, *.jpeg)'; ...
    '*.*',              'All Files (*.*)'};
[synthname,synthpath] = uiputfile(imfilter, 'Save the Synthesized Image as...', strcat(imbasename, '-out.png'));

% Check it
noSynth = ~ischar(synthname) || ~ischar(synthpath);
if (~noSynth)
    synthpath = fullfile(synthpath, synthname);
else
    synthname = [];
end;
disp(['>> ', synthname]);


%% Dictionary

% Load or learn?
disp(' ');
disp('Dictionary:')
disp('   0. learn it');
disp('  [1] load it');
loadDico = logical(prompt(1));


% Dictionary file
disp(' ');
if (loadDico)
    disp('Loading the dictionary...');
    [diconame,dicopath] = uigetfile('*.mat', 'Choice of the Dictionary', strcat(imbasename, '.mat'));
else
    disp('Saving the dictionary as...');
    [diconame,dicopath] = uiputfile('*.mat', 'Save the Dictioanry as...', strcat(imbasename, '.mat'));
end;

% Check it
if (~ischar(diconame) || ~ischar(dicopath))
    error('Interrupted');
else
    dicopath = fullfile(dicopath, diconame);
end;
disp(['>> ', diconame]);

% Load it
if (loadDico)
    [D,im] = Synthesis.dictionary(dicopath, impath);
    patches = Patches([], D(:,1:2,1));
    disp('> Loaded!');
elseif (exist(dicopath, 'file'))
    movefile(dicopath, strcat(dicopath, '.bak'));
    disp(['> Note: existing dictionary renamed as ', diconame, '.bak']);
end;


%% Other parameters

% Scales
disp(' ');
disp('Number of scales:');
if (~loadDico)
    disp('   1 ... [3] ...');
	nScales = promptInt(1, 3, Inf);
elseif (~noSynth)
    def = size(D, 3);
    disp(['   1 ... [', num2str(def), ']']);
	nScales = promptInt(1, def, def);
    D = D(:,:, 1:nScales);
else
    def = size(D, 3);
    disp(['>> ', num2str(def)]);
end;

% Weights
if (~noSynth)
    disp(' ');
    disp('Weights for the 3 terms [histogram, spectrum, patch]:');
    disp('   def. is [1 1 1]');
    wmat = input(promptPrefix);
    if (isempty(wmat))
        wmat = [1 1 1];
    elseif (numel(wmat) ~= 3)
        error('Weights must be a 3-element vector');
    end;
    weights = struct('histo', wmat(1), 'spectrum', wmat(2), 'patch', wmat(3));
    disp(['> histo: ' num2str(wmat(1)) ', spectrum: ' num2str(wmat(2)) ', patch: ' num2str(wmat(3))]);
end

% Patch size
disp(' ');
disp('Size of the patches:');
if (loadDico)
    pSize = patches.Nyi;
    disp(['>> ', num2str(pSize)]);
else
    disp('   2 ... [12] ...');
    pSize = promptInt(2, 12, Inf);
end;

% Patch spacing
if (~noSynth)
    disp(' ');
    disp('Subsampling of the patches during synthesis:');
    def = 4;
    disp(['   1 ... [', num2str(def), '] ... ', num2str(pSize)]);
    pSpacing = promptInt(1, def, pSize);
end;

% Sparsity
if (~noSynth || ~loadDico)
    disp(' ');
    disp('Sparsity:');
    disp('   1 ... [4] ...');
    sparsity = promptInt(1, 4, Inf);
end;

% Dico size
disp(' ');
disp('Size of the dictionary:');
if (loadDico)
    disp(['>> ', num2str(size(D, 2))]);
else
    def = 2 * pSize^2;
    disp(['   2 ... [', num2str(def), '] ...']);
    dicoSize = promptInt(2, def, Inf);
end;

% Iterations (learning)
if (~loadDico)
    disp(' ');
    disp('Number of iterations per scale for learning:');
    disp('   1 ... [30] ...');
    learnIter = promptInt(1, 30, Inf);
end;

% Iteration (synthesis)
if (~noSynth)
    disp(' ');
    disp('Number of iterations per scale for synthesis:');
    disp('   1 ... [30] ...');
    synthIter = promptInt(1, 30, Inf);
end;


%% Learning and synthesis

% Learning
if (~loadDico)
    disp(' ');
    disp('Learning the dictionary... this may take some time...');
    [D,im] = Synthesis.dictionary(dicopath, impath, nScales, ...
        pSize, dicoSize, sparsity, learnIter);
    disp('> Learned!');
end;

% Synthesis
if (~noSynth)
    disp(' ');
    disp('Synthesizing... this may take some time...');
    synth = Synthesis.launch(D, im, outSize, weights, sparsity, pSpacing, synthIter);
    imwrite(synth, synthpath);
    disp('> Synthesized!');
end;

%% Display

clf;
disp(' ');
disp('Displaying the image(s)...');
nPlot = 2 + ~noSynth;

% Input image
subplot(1,nPlot,2);
imshow(im);
title(sprintf('Input, %dx%d', size(im, 2), size(im, 1)));

% Synthesized image
if (~noSynth)
    subplot(1,nPlot,3);
    imshow(synth);
    title(sprintf('Synthesis, %dx%d', outSize(end), outSize(1)));
end;

% The dictionary
disp('Displaying the dictionary...');
dicoMosaic = Patches([], D(:,:,1), 1, true);
subplot(1,nPlot,1);
imshow(abs(dicoMosaic.image));
title(sprintf('Dictionary, patches %dx%d', pSize, pSize));
disp('> Done!');


%% End the session
disp(' ');
disp('----- Interactive session is over! -----');
