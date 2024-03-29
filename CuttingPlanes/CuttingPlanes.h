#pragma once

#include <vector>
#include "Utils.h"
#include "soplex.h"

class CuttingPlanes
{
public:
    explicit CuttingPlanes();

    void solve(const Vec& c, const Mat& A, const Vec& b, const std::vector<bool>& I);

    void export_json(const std::string& filename);

private:
    soplex::SoPlex simplex;

    Vec c;
    Mat A;
    Vec b;
    uint n; // dimensionality
    uint m; // number of (original) inequalities
    uint n_cuts;
    std::vector<bool> I; // Integer constraints

    std::vector<Vec> simplex_solutions;

    /**
     * Returns the first variable index where the current solution violates an integer constraint
     * or -1 if the solution does not violate any integer constraints
     */
    int integer_constraint_violation(const Vec& x);

    /**
     * Initializes the simplex solver, setting up the system
     * min cx s.t. Ax <= b
     */
    void init_simplex(const Vec& c, const Mat& A, const Vec& b);

    /**
     * Solves the system
     * min cx s.t. Ax <= b
     */
    bool solve_simplex(Vec& primal);

    /**
     * Adds the linear constraint ax <= d to the system
     */
    void add_constraint(const Vec& a, const double& d);
};
