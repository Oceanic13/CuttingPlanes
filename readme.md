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
are valid for any feasible point but cut off the relaxed solution if it violates any integral constraints. This process is repeated until a feasible solution has been found.

Some textbooks for the interested:
- [John Wiley and Sons, Ltd - 2020 - Integer Programming](https://doi.org/10.1002/9781119606475.oth1)
- [Papadimitriou, Christos and Steiglitz, Kenneth - 1982 - Combinatorial Optimization: Algorithms and Complexity](https://doi.org/10.1109/TASSP.1984.1164450)


This project contains
- a c++ implementation of the Simplex Algorithm to solve LPs ([cpp](./CuttingPlanes/))
```cpp
#include "Utils.h"
#include "MixedIntegerLinearProgram.h"
#include "Solvers.h"
using MILP = MixedIntegerLinearProgram

// Make a MILP
// minimize -y s.t. 3x+2y <= 6 and -3x+2y <= 0 and x,y >= 0 and x,y integral
Vecd c(2); c << 0, -1; // objective
Matd A(2,2); A << 3, 2, -3, 2; // inequalities
Vecd b(2); b << 6, 0;
Matd B(0,2); // equalities
Vecd d(0);
int n1 = 2; // number of integral constraints
auto milp = MILP(c, A, b, B, d, n1);

auto solver = ToblexSolver(milp);
solver.solve();

Vec2d optimal_solution = solver.getOptimalSolution(); // (1, 1.5)
double optimal_value = solver.getOptimalValue(); // - 1.5

```

- a c++ implementation of the Cutting Planes Method using Mixed Integer Gomory Cuts to solve MILPs ([cpp](./CuttingPlanes/))
```cpp
#include "Utils.h"
#include "MixedIntegerLinearProgram.h"
#include "Solvers.h"
using MILP = MixedIntegerLinearProgram

// Make a MILP
// minimize -y s.t. 3x+2y <= 6 and -3x+2y <= 0 and x,y >= 0 and x,y integral
Vecd c(2); c << 0, -1; // objective
Matd A(2,2); A << 3, 2, -3, 2; // inequalities
Vecd b(2); b << 6, 0;
Matd B(0,2); // equalities
Vecd d(0);
int n1 = 2; // number of integral constraints
auto milp = MILP(c, A, b, B, d, n1);

auto solver = CuttingPlanesSolver(milp);
solver.solve(); // generates two cutting planes until the solution is integral

Vec2d optimal_solution = solver.getOptimalSolution(); // (1, 1)
double optimal_value = solver.getOptimalValue(); // - 1

```
- a python program to visualize the Cutting Planes Method in 2D ([python](./Visualizer.ipynb))
  1. Build the cpp project using cmake
  2. Run ([Visualizer.ipynb](./Visualizer.ipynb)) (specify the correct build path)
  3. Have fun! Mixed Integer Linear Programs can be stored in json format imported/exported in the cpp project as well as in the python project.
