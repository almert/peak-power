# Peak-Power Capacity
### MATLAB codes for computing the maximum radius constraints on vector Gaussian channels

## Overview
This project contains the MATLAB codes for finding the maximum radius constraints such that the optimal distributions are supported on a single sphere. This is performed for vector Gaussian channels under two settings: (1) finding the input distribution that maximizes the mutual information, (2) finding the 'least-favorable prior distribution' that maximizes the minimum mean square error (MMSE).

The codes were used to generate the reported radius values in [[1](#citation)] (https://arxiv.org/abs/1804.08524).

## Usage
The codes for the optimal mutual information setting are under "RMax Mutual Information" and the codes for the maximum MMSE setting are under "RMax MMSE". We also included .mat files, which contain the results reported in the paper. 

The maximum radius for each setting has been shown to be the solution of an integral equation. The codes perform binary searches to find the solutions of these equations, where the integrals are evaluated via Monte Carlo methods. The equation for the mutual information setting can be found in Theorem 3 in the paper, while the equation for the MMSE setting can be found under the Discussion section.

Please cite [[1](#citation)] in your work when using these codes in your experiments.

## License
Princeton University

## Citation
```
[1] Dytso, Alex, Mert Al, H. Vincent Poor, and Shlomo Shamai. On the Capacity of the Peak Power Constrained Vector Gaussian Channel: An Estimation Theoretic Perspective. 2018. arXiv preprint arXiv:1804.08524.
```

BibTeX format:
```
@article{dysto2018capacity,
  title     = {On the Capacity of the Peak Power Constrained Vector Gaussian Channel: An Estimation Theoretic Perspective},
  author    = {Alex Dytso and
               Mert Al and
               H. Vincent Poor and
               Shlomo Shamai},
  journal   = {CoRR},
  volume    = {abs/1804.08524},
  year      = {2018},
  url       = {http://arxiv.org/abs/1804.08524},
}
```
