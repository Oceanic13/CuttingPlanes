#pragma once

#include <vector>
#include "Typedefs.h"

class CuttingPlanes
{
public:
    explicit CuttingPlanes(const Vec& c, const Mat& A, const Vec& b, const std::vector<bool>& I);

    void solve();

    void export_json(const std::string& filename);

private:
    uint n; // dimensionality
    uint m; // number of (original) inequalities
    Vec c; // Objective Function (cx)
    Mat A; // Inequality Matrix (Ax <= b), grows with each iteration
    Vec b; // Inequality Vector (Ax <= b), grows with each iteration
    std::vector<bool> I; // Integer constraints

    std::vector<Vec> simplex_solutions;

    /**
     * Adds the linear constraint ax <= d to the system
     */
    inline void add_constraint(const Vec& a, const double& d)
    {
        assert(A.rows() == b.size());
        uint k = A.rows();
        A.conservativeResize(k+1, Eigen::NoChange);
        A.row(k) = a;
        b.conservativeResize(k+1);
        b[k] = d;
    }
};
