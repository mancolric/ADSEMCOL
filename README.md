## ADSEMCOL

ADSEMCOL (Adaptive Spectral Element Method for Systems of Conservation Laws) is an in-house finite element code for the solution of systems of conservation laws. Its main features are:
- Continuous finite elements up to degree 10
- Anisotropic h-adaptation
- Adaptive, 3rd order, linearly implicit Runge--Kutta time stepping

Currently solves:
- Nonlinear scalar convection--diffusion problems
- Compressible Euler and Navier--Stokes equations

Under development (to be included in future updates):
- Hydrogen combustion
- Shallow water equations

The code is currently under development, i.e., it is not optimized and it can sometimes fail due to incorrect settings and not sufficiently robust implementation.

![](/test/figures/TriplePoint.png)

### Authors and contact

Manuel Colera-Rico (Universidad Politécnica de Madrid)\
Vít Dolejší (Charles University Prague)

m.colera@upm.es

### Cite as 

Colera, M., Dolejší, V., 2025. An adaptive spectral element method for systems of conservation laws. Under review in Journal of Computational Physics.

### Installation

Requires BAMG mesh generator, available via [FreeFem](https://freefem.org/), and [Julia](https://julialang.org/) language. 

In the main instalation directory, run `julia install` and follow the instructions to install all the necessary files and Julia packages. REVISAR

### Usage

The directories `ADSEMCOL/test/test_ConvectionDiffusion` and `ADSEMCOL/test/test_CompressibleFlow` contain several tests concerning the nonlinear convection-diffusion equation and the Euler and Navier--Stokes compressible flow equations.

Each test is defined as a function that can be executed from the Julia terminal. Normally, the function has two compolsory arguments: the resolution $h/p$ of the initial mesh (necessary to start the mesh adaptation process), and the discretization order $p$. After that, there are many optional parameters with default values which can be seen in the function header.

As an example, we consider the the triple-point test with an initial mesh size $h=0.10$, space discretization order $p=5$, space tolerance $\mathrm{Tol}_S=10^{-5}$, time tolerance $\mathrm{Tol}_T=10^{-5}$, artificial viscosity $\epsilon=3 \cdot 10^{-5}$ and final time $t_f=0.1$. We need to save the results every 0.01 time units. We execute

```julia

#change directory to test_CompressibleFlow:
cd(path-to-test_CompressibleFlow)

#includet is recommended to make changes in the file without needing to recompile:
using Revise 	
includet("TriplePoint.jl") 

#Run triple-point problem
TriplePoint(0.02, 5, TolS=1e-5, TolT=1e-3, epsilon=3e-5, tf=0.1, SaveRes=true, Deltat_SaveRes=0.01, SC=10)

```

Above, `SC` is the _study case_, i.e., a number associated with the problem. In this case, `SC=10`. The results are saved in the directory specified in the `ResUbi.jl` file (generated during the installation REVISAR) with the names `SC10_0.jld2` ($t=0$), `SC10_1.jld2` ($t=0.01$), `SC10_2.jld2` ($t=0.02$), ...

Detailed documentation will be uploaded in the future.
