# Cutting Planes
Selected topic in the course "Seminar Applied Optimization" in the spring semester 2024 at the university of Bern.

Cutting Planes are a method to solve a Mixed Integer Linear Program (MILP) of the form

$$
\begin{align}
\min_x \quad &c \cdot x \\
\text{s.t.}\quad &Ax \leq b \\
&Bx = d \\
&x \geq 0 \\
&x \in \mathbb{Z}^{n_1} \times \mathbb{R}^{n_2}
\end{align}
$$

by iteratively solving the problem relaxation (where the integral constraints are ignored) and adding additional inequalities (so called cutting planes) which
are valid for any feasible point but cut off the relaxed solution if it violates any integral constraints.


This project contains
- a c++ implementation of the Simplex Algorithm to solve LPs
```cpp
Hello World
```

- a c++ implementation of the Cutting Planes Method using Mixed Integer Gomory Cuts to solve MILPs
- a python program to visualize the Cutting Planes Method in 2D
