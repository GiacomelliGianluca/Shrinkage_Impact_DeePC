# Title

Matlab implementation of data-driven control strategies as described in 
"Paper Title."

## Prerequisites

You will need:

- `Matlab`
- `Matlab statistics_toolbox` (see https://nl.mathworks.com/products/statistics.html)
- `CVX` (see https://cvxr.com/)
- `Hybrid Toolbox` (see https://cse.lab.imtlucca.it/~bemporad/hybrid/toolbox/)

## Installation

To clone this repository, see https://nl.mathworks.com/help/simulink/ug/clone-git-repository.html

## Example: 

The considered PieceWise Affine (PWA) system is described as follows

```math
\begin{align}
     &x_{t+1} = \begin{cases}
         -0.3 x_t + 1.4 u_t &\text{if}~x_t<0,\\
         0.9 x_t + 0.15 u_t &\text{if}~x_t\geq 0,
\end{cases}\\
 &y_t = x_t,
\end{align}
```
