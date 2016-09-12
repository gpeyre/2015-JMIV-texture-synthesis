# 2015-JMIV-texture-synthesis

Source code of the article

G. Tartavel, Y. Gousseau, G. Peyré. [Variational Texture Synthesis with Sparsity and Spectrum Constraints](https://hal.archives-ouvertes.fr/hal-00881847/). Journal of Mathematical Imaging and Vision, 52(1), pp. 124–144, 2015.

![Examples of syntheses](videos/synthesis.png)


Description
-----------

This is an implementation of our texture synthesis algorithm.
It is an exemplar-base approach: the goal is to generate a new texture from a sample image.

More details about the algorithm, the variational approach and the sparsity constraints are available in our publications.

This implementation is in Matlab and is composed of a set of *.m files.
Documentation is available in Matlab for each function "FCN" by typing:
>> help FCN


Usage
-----

Compile the MEX functions.
You can use the Makefile provided in 'mex/':
$ cd mex/ && make

Launch Matlab and execute the editable script or the interactive script:
>> LaunchExample;
>> LaunchInteractive;

You can use the matlab 'help' command to explore the files:
>> help Synthesis

Note about the Image Processing toolbox in Matlab.
If you don't have it, replacement functions can be found in 'improc_toolbox/'.
You can for example add this directory to Matlab's search path:
>> addpath('improc_toolbox');

Note about the SPAMS toolbox.
The SPAMS toolbox provide a more efficient dictionary learning step.
It is available here:
* http://spams-devel.gforge.inria.fr/


File descriptions
-----------------

README.txt          -- this information file

*.png               -- some images to try the algorithm and the toolboxes
*.mat               -- some dictionary already learned (to save time)

LaunchCopymap.m     -- editable script to analyze the innovation capacity
LaunchExample.m     -- editable script to launch a synthesis
LaunchInteractive.m -- interactive script to launch a synthesis

Synthesis.m         -- implements our synthesis algorithm
Patches.m           -- interface for patches in an image
Histogram.m         -- provides histogram tools (histogram transfer)
Sparsity.m          -- provides sparse decomposition functions
Spectrum.m          -- provides Fourier spectrum tools (spectrum transfer)

mex/                -- contains the C++ implementation of some MEX functions
mex/Makefile        -- useful to compile/generate the MEX files

improc_toolbox/     -- contains a basic implementation of some function
                       provided in Matlab's Image Processing toolbox


Copyright (c) 2015 Guillaume Tartavel
