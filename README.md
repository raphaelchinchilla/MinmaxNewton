# MinmaxNewton

The code in this repository implements the algorithm and numerical examples from our paper [Newton and interior-point methods for (constrained) nonconvex-nonconcave minmax optimization with stability guarantees](https://arxiv.org/abs/2205.08038).

The code is implemented in MATLAB and requires installing [TensCalc](https://github.com/hespanha/tenscalc) (and its dependencies). TensCalc is used as a backend to compute the symbolic differentiation and to generate optimize code, either MATLAB or C++.

The bulk of the algorithm is implemented in two files:
- `generate_tens_functions.m` parses the symbolic optimization problem and generates the code to compute the appropriate gradients and Hessians.
- `ip_newton_minmax.m` implements the interior-point (and consequently, Newton method) described in the paper, with the appropriate Hessian modifications.

The two files in the main folder of the repository, `benchmark.m` and `pursuit_evasion.m`, implement the numerical examples from our paper.
