# Optimal Transport Metrics for Source Localization
This package provides the code to compute the Flat Metric, which is the metric used in the paper <a href="https://arxiv.org/abs/2010.13423">Optimal-transport-based metric for SMLM</a> for performance assessment in Single Molecule Localization Microscopy (SMLM).

The computation of the metric is done by solving an optimization problem formulated as a linear program thanks to a standard linear solver.

# Installation
The code is in Julia and to work properly requires:
<ul>
	<li><a href="https://julialang.org">Julia 1.0.5</a> (LTS release),</li>
	<li>the <a href="https://github.com/jump-dev/JuMP.jl">JuMP</a> package for optimization,</li>
	<li>the solver Mosek (and the package MosekTools). Beware that using Mosek requires a license. Alternatively, one can use the solver <a href="https://github.com/jump-dev/GLPK.jl">GLPK</a>,</li>
	<li>optional (for notebooks and plots): package <a href="https://github.com/timholy/Revise.jl">Revise</a> and <a href="https://github.com/JuliaPy/PyPlot.jl">PyPlot</a>.</li>
</ul>

# Examples

Chec

