Sketch Token Toolbox     V0.95
==============================

This software package provides tools to extract contour-based mid-level features, and to extract contour segmentations from images.
This tool is highly efficient in speed while maintains high accuracy in contour detection.
Also, [1] shows that extracted mid-level features provide additional information for object and pedestrian detections.

Installation
------------
<ol>
<li>
Download Piotr's Image & Video Matlab Toolbox (http://vision.ucsd.edu/~pdollar/toolbox/doc/)<br>
    SketchTokens/toolbox/ should have channels, classify, filters, images, matlab, etc,.
</li>

<li>
Download Berkeley Segmentation Data Set and Benchmarks 500 (http://www.eecs.berkeley.edu/Research/Projects/CS/vision/grouping/resources.html)<br>
    SketchTokens/data/BSR/ should have BSDS500, bench, and documentation.
</li>

<li>
Pre-trained models can be downloaded from: http://people.csail.mit.edu/lim/lzd_cvpr2013/st_data.tgz
</li>

<li>
Look up stDemo.m for how to train and test our code
</li>
</ol>

References
----------
Please cite the following paper if you end up using the code:<br>
[1] Joseph J. Lim, C. Lawrence Zitnick, and Piotr Dollar. "Sketch Tokens: A Learned Mid-level Representation for Contour and and Object Detection," CVPR2013.

License
-------
Copyright 2013 Joseph Lim [lim@csail.mit.edu]

Please email me if you find bugs, or have suggestions or questions!

Licensed under the Simplified BSD License [see bsd.txt] <br>

<b>Note: There is a patent pending on the ideas presented in this work so this code should only be used for academic purposes.</b><br>
