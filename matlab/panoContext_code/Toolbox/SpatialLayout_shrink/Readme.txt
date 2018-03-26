%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Spatial Layout
Varsha Hedau (vhedau2@uiuc.edu)
10/27/2010
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

This distribution contains data and code for estimating spatial layout of a scene,
as described in [1]. The software outputs a) A box layout of the room,
approximating it as a 3Dbox, and, b) The likelihood of each pixel in the image
of being one of the planar surfaces of the room, or object clutter. Please
refer to [1] for the detailed meaning of these outputs

IMPORTANT:
To use this software, YOU MUST CITE the following in any resulting publication:

[1] Varsha Hedau, Derek Hoiem, David Forsyth, “Recovering the Spatial 
    Layout of Cluttered Rooms,” in the Twelfth IEEE International Conference 
    on Computer Vision, 2009.

[2] Varsha Hedau, Derek Hoiem, David Forsyth, “Thinking Inside the Box: 
    Using Appearance Models and Context Based on Room Geometry,” ECCV 2010.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

LICENSE:

Copyright (C) 2011 Varsha Hedau, University of Illinois at Urbana Champaign.

This software is available for non-commercial use only. The software is 
available for general use by academic or non-profit, or government-sponsored 
researchers. This license does not grant the right to use this software or 
any derivation of it in a for-profit enterprise. Any redistribution of the 
software shall require a written authorization by the authors.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

CONTENTS:

- The directory 'Images_resized' contains the original images used for experiments
in [1]. 

- The directory 'Imsegs' contains segmentation images

- The directory 'LabelClassifiers' contains surface label classifiers and ground truth data

- The directory 'LearntClassifiers' contains box layout classifiers 

- The directory 'spatiallayoutcode' contains source code. Please see the 'demo_script.m'
for a sample run.

- 'traintestinds.mat' contain train-test split used in [1].


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

WHAT TO RUN?

Sample code for generating output on all test images used in [1] is
provided in demo_script.m.


 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

REFERENCES:

[1] Varsha Hedau, Derek Hoiem, David Forsyth, “Recovering the Spatial 
    Layout of Cluttered Rooms,” in the Twelfth IEEE International Conference 
    on Computer Vision, 2009.

[2] Varsha Hedau, Derek Hoiem, David Forsyth, “Thinking Inside the Box: 
    Using Appearance Models and Context Based on Room Geometry,” ECCV 2010.
