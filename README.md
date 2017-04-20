# SphericalSFM

This repository will contain Matlab implementations of the key steps in a spherical structure-from-motion pipeline. The code will be added to the repository over time as we clean and comment the code. The following list will describe the key functions and will be updated as we add more functionality:

1. SnP.m - This implements three solutions to the spherical-n-point problem. It is called with:
```matlab
    [ R_est,t_est ] = SnP( x,X,method,niter,gamma )
```
where x is a 3 by npoints matrix containing unit vectors corresponding to spherical image points, X is a 3 by npoints matrix containing 3D world points in correspondence with x, method is either 'hard', 'soft' or 'unconstrained' for the different methods, niter is the number of iterations of the hard or soft constrained methods, gamma is the weight of the soft constraint.

Reference
---------

If you use this code in your research, please cite the following paper:

H. Guan and W.A.P. Smith. Structure-from-motion in Spherical Video using the von Mises-Fisher Distribution. IEEE Transactions on Image Processing, Volume 26, Number 2, pp. 711-723, 2017.

Bibtex:

    @article{guan2017structure, 
        author={H. Guan and W. A. P. Smith}, 
        journal={IEEE Transactions on Image Processing}, 
        title={Structure-From-Motion in Spherical Video Using the von Mises-Fisher Distribution}, 
        year={2017}, 
        volume={26}, 
        number={2}, 
        pages={711-723}
    }
    
Dependencies
------------

If you use the soft constrained method in SnP.m then you need to install the CVX toolbox, available here: http://cvxr.com/cvx/
