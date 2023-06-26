# PWSPD for Clustering Data
This repository contains code accompanying the paper:
```
@article{mckenzie2019power,
  title={Power weighted shortest paths for clustering Euclidean data},
  author={Mckenzie, Daniel and Damelin, Steven},
  journal={arXiv preprint arXiv:1905.13345},
  year={2019}
}
```
Code to reproduce the experiments is contained in the folder ```BenchMarking```. 
Matlab code for using Power Weighted Shortest Path Distances (AKA Fermat distances) in spectral clustering. The modified version of Dijktra's algorithm proposed in the paper is implemented as ```Dijkstra_with_early_stopping```. Note you may have to compile some ```.mex``` files to get this to work. If you find this code useful please consider citing the paper mentioned above! 
