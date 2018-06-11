# LayoutNet
Torch implementation for CVPR 18 [paper](https://arxiv.org/pdf/1803.08999.pdf): "LayoutNet: Reconstructing the 3D Room Layout from a Single RGB Image"

See sample [video](https://youtu.be/WDzYXRP6XDs) of 3D reconstruced layouts by our method.

<img src='figs/teasor.jpg' width=400>

## Prerequisites
- Linux
- NVIDIA GPU + CUDA CuDNN
- Torch 7

matio: https://github.com/tbeu/matio

- Matlab

## Data
- Download preprocessed (aligned to horizontal floor plane) training/validation/testing [data](https://drive.google.com/file/d/1vsIvZ5L-VT0sH-GgbUL1sRYiEHn2Jn3B/view?usp=sharing) to current folder

This includes the panoramas from both the panoContext dataset and our labeled stanford 2d-3d dataset.

- Download groundtruth [data](https://drive.google.com/file/d/1j91sz8Jt6Jsg198riA0ggz8Mjj4lSntx/view?usp=sharing) to current folder

This includes the groundtruth 2D posiiton of room corners in .mat format from the two dataset

## Pretrained model
- Download our pretrained [model](https://drive.google.com/file/d/1qqrKkT_nTN1RzjiLN92VvoB023ZoD28v/view?usp=sharing) to current folder

This includes the pretrained full approach on the panoContext dataset, the joint boudary and corner prediction branch, the single boundary prediction branch and the 3D layout box regressor.

## Image preprocess
We provide sample script in ./matlab/getManhattanAndAlign.m

## Train network
- To train our full approach:
```
th driver_pano_full.lua
```
Note that this loads the pretrained joint prediction branch and the 3D layout box regressor.

- To train the joint prediction branch of boudary and corner:
```
th driver_pano_joint.lua
```
Note that this loads the pretrained boundary prediction branch.

- To train the boudary prediction branch
```
th driver_pano_edg.lua
```
- To train the layout box regressor:
```
th driver_pano_box.lua
```

## Test network
- To test on our full approach:
```
th testNet_pano_full.lua
```
This saves predicted boundary, corner and 3D layout parameter in "result/" folder

## Optimization
- To Add Manhattan constraints and optimize for a better layout, open Matlab, then:
```
cd matlab
panoOptimization.m
```
This loads saved predictions from the network output and performs sampling.

## Evaluation

We provide the Matlab evaluation code for 3D IoU (compute3dOcc\_eval.m) and the generation of 2D layout label (getSegMask\_eval.m) for evaluating layout pixel accuracy.

## Extension to perspective images
Main training code in driver\_pano\_joint\_lsun.lua

## Citation
Please cite our paper for any purpose of usage.
```
@article{zou2018layoutnet,
  title={LayoutNet: Reconstructing the 3D Room Layout from a Single RGB Image},
  author={Zou, Chuhang and Colburn, Alex and Shan, Qi and Hoiem, Derek},
  journal={arXiv preprint arXiv:1803.08999},
  year={2018}
}
```
