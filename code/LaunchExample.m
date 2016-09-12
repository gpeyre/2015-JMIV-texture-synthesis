% LAUNCH EXAMPLE    - editable script to launch a synthesis.
%
%   You can set the parameters in the first part of this file.
%   See the comments for each variable in the 'Parameters' section below.
%
%   See also:
%       LaunchInteractive
%
%   Contact:
%       Guillaume.Tartavel @ telecom-paristech.fr
%       Downloadable on perso.enst.fr/~tartavel/

clc;
clear;
close all;


%% Parameters -- edit them!

% Files
name = 'sand';                  % Name of the image (without extensions)
imExt = '.png';                 % Extention of the image (with the dot)
imDir = './';                   % Directory of the image
dicoDir = './';                 % Directory to save the dictionaries
synthDir = './';                % Directory to save the synthesis
synthName = 'sand-out';         % Name of the synthesized image (w/o ext)

% Synthesis
outSize = [256 256];            % Size of the output
nScales = 3;                    % Number of scales
pSize = 12;                     % Size of the patches
pSpacing = 4;                   % Spacing of the patches
sparsity = 4;                   % Sparsity of the patch decomposition
synthIter = 30;                 % # of iteration for synthesis (per scale)
weights = struct('histo', 1, 'spectrum', 1, 'patch', 1);

% Dictionary
loadDico = true;                % If true, load the dictionary if existing
dicoSize = 256;                 % Size of the dicitonary (if learned)
learnIter = 30;                 % # of iterations for learning (if learned)


%% Some initializations

dicoPath = fullfile(dicoDir, strcat(name, '.mat'));
imPath = fullfile(imDir, strcat(name, imExt));
synthPath = fullfile(synthDir, strcat(synthName, '.png'));


%% Load or learn the dicitonary

if (~loadDico)
    if (exist(dicoPath, 'file'))
        movefile(dicoPath, strcat(dicoPath, '.bak'));
        warning('FILE:AlreadyExists', 'Dictionary already exist, back-up copy created');
    end;
end;

[D,im] = Synthesis.dictionary(dicoPath, imPath, nScales, ...
        pSize, dicoSize, sparsity, learnIter);


%% Perform the synthesis

if (exist(synthPath, 'file'))
    fprintf('Image %s has already been synthesized\n', name);
    fprintf('    -> %s\n', synthPath);
else
    synth = Synthesis.launch(D, im, outSize, weights, sparsity, pSpacing, synthIter);
    imwrite(synth, synthPath);
    clf; imshow(synth); title('Synthesized image');
end;
