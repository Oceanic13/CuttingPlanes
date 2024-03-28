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
    Vec c; // Objective Function (cx)
    Mat A; // Inequality Matrix (Ax <= b)
    Vec b; // Inequality Vector (Ax <= b)
    std::vector<bool> I; // Integer constraints

    std::vector<std::pair<Vec, double>> cutting_planes;

    std::vector<Vec> simplex_solutions;
};
