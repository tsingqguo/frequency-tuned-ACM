# frequency-tuned-ACM
This repository contains the demos of frequency-tuned ACM published in ICASSP 2017.

In the frequency domain, we propose a generalized Region based ACM that presents a new way to understand the essence of classical RACMs whose segmentation results are determined by a frequency filter to extract the proposed frequency boundary energy. Then, we introduce the difference of Gaussians as the optimal filter to exclude strong noise and intensity inhomogeneity effectively. Please refer to our project webside http://tsingqguo.gq/fbeproject.html for details.

'Demo.m' contains 3 demos. Please set 'demo' in the code as 1,2 or 3 to run each demo. fbe_acm is the proposed algorithm. The files suffixed by '.mat' contain parameters for each demo.

if you uses the provided codes, please cite the paper:
@inproceedings{Guo2017,
  author={Q. Guo and S. Sun and F. Dong and W. Feng and S. Ma and B. Z. Gao},
  title={Frequency-tuned Active Contour Model for Biomedical Image Segmentation},
  booktitle = {International Conference on Acoustics, Speech and Signal Processing},
  year = {2017}
}
