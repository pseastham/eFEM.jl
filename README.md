# eFEM.jl

<!--[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://pseastham.github.io/eFEM.jl/stable)-->
<!--[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://pseastham.github.io/eFEM.jl/dev)-->
[![Build Status](https://travis-ci.com/pseastham/eFEM.jl.svg?branch=master)](https://travis-ci.com/pseastham/eFEM.jl)
[![Codecov](https://codecov.io/gh/pseastham/eFEM.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/pseastham/eFEM.jl)

Finite Element code in the [Julia language](https://julialang.org/) focused on fluid-dynamics applications.

This repository allows the use of Finite Elements discretizations to solve
common problems in fluid dynamics. Functionally, this package is designed to be used
with the larger [eFEMpart.jl](https://github.com/pseastham/eFEMpart.jl), which also includes the ability to simulate particles.

## Installation

eFEM.jl is an unregistered package. To install eFEM.jl using the REPL, type

`pkg> add https://github.com/pseastham/eFEM.jl`

To enter the Pkg environment from the REPL, type `]`.

## Meshes

For simple geometries (rectangles...), you can use the built-in 
geometry code but for more complicated geometries we suggest building your mesh with 
an external library (such as [GMSH](http://gmsh.info/)).

## Equations

As of right now, the following equations are solvable:

* Poisson's Equation (`:Poisson2D`)

<p align="center"><img src="/tex/a1e55dd0d6f8247d8b884e241419c34e.svg?invert_in_darkmode&sanitize=true" align=middle width=75.003885pt height=17.399144399999997pt/></p>

* Darcy's Equation (`:Darcy2D`)

<p align="center"><img src="/tex/3ba9ca5ab07d4c987d667c9f4956512c.svg?invert_in_darkmode&sanitize=true" align=middle width=118.8451539pt height=19.726228499999998pt/></p>

* Advection-Diffusion Equation (`:AdvDiff2D`)

<p align="center"><img src="/tex/50aaf8695606a64a2aba3412a4cd7ca3.svg?invert_in_darkmode&sanitize=true" align=middle width=178.72117724999998pt height=19.726228499999998pt/></p>

* Stokes' Equation (`:Stokes2D`)

<p align="center"><img src="/tex/f7e35892f79b733caf605eb9762d82c0.svg?invert_in_darkmode&sanitize=true" align=middle width=170.03593694999998pt height=19.726228499999998pt/></p>
<p align="center"><img src="/tex/efbfbcd0f130f2b91fea06b34868e681.svg?invert_in_darkmode&sanitize=true" align=middle width=66.2097216pt height=11.232861749999998pt/></p>

* Brinkman's Equation (`:Brinkman2D`)

<p align="center"><img src="/tex/07e57a540d72768f0e3d8ca41934ad8a.svg?invert_in_darkmode&sanitize=true" align=middle width=200.24691225pt height=19.726228499999998pt/></p>
<p align="center"><img src="/tex/efbfbcd0f130f2b91fea06b34868e681.svg?invert_in_darkmode&sanitize=true" align=middle width=66.2097216pt height=11.232861749999998pt/></p>

* Brinkman's Multiphase Equation (`:BrinkmanMP2D`)

<p align="center"><img src="/tex/8b86a228922df2e57a458e4cbd5379e8.svg?invert_in_darkmode&sanitize=true" align=middle width=186.47236739999997pt height=17.399144399999997pt/></p>
<p align="center"><img src="/tex/efbfbcd0f130f2b91fea06b34868e681.svg?invert_in_darkmode&sanitize=true" align=middle width=66.2097216pt height=11.232861749999998pt/></p>

All parameterized equations can be solvable with either constant or variable-in-space parameters. Additionally, Axisymmetric version of the Advection-Diffusion and Stokes equations are available with the Operator Types of `:AdvDiffAS` and `:StokesAS`, respectively. 

Equations in this `README` were generated by the GitHub app [TeXify](https://github.com/apps/texify)

## Auxiliary Information

Boundary conditions are treated intuitively, based on the mesh given. The functions `Dirichlet`, `Neumann`, and `Robin` allow assignment of boundaries to have certain boundary conditions, and the functions `Dirichlet`, `Neumann`, `Forcing` allow for the definition of the actual boundary conditions at these boundaries. See [examples](examples/) for how this is used in practice.

## Visualization

We export all solutions in a [legacy VTK format](https://www.vtk.org/VTK/img/file-formats.pdf). For visualizing these files, we suggest using [VisIt](https://wci.llnl.gov/simulation/computer-codes/visit/) or [ParaView](https://www.paraview.org/).

## Examples 

Check out the [examples folder](examples/) to see how to use our syntax.
