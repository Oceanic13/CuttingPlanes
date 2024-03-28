#include "CuttingPlanes.h"
#include <iomanip>
#include <iostream>
#include <fstream>

CuttingPlanes::CuttingPlanes(const Vec& c, const Mat& A, const Vec& b, const std::vector<bool>& I)
    : c(c), A(A), b(b), I(I)
{}

void CuttingPlanes::solve()
{
    cutting_planes.clear();
    simplex_solutions.clear();

    Vec t(2); t << 0, 1;
    cutting_planes.emplace_back(t, 1);
    t << 1, 1.5;
    simplex_solutions.push_back(t);

    std::cout << "Hello World" << std::endl;
}

void CuttingPlanes::export_json(const std::string& filename)
{
    json j;

    j["system"] = {
        {"n", c.size()},
        {"m", b.size()},
        {"c", c},
        {"b", b},
        {"A", A.rowwise()},
        {"I", I}
    };

    j["algo"] = {
        {"simplex", simplex_solutions},
        {"cuts", cutting_planes}
    };

    std::cout << std::setw(4) << j << std::endl;

    std::ofstream o(filename);
    o << std::setw(4) << j << std::endl;
    o.close();
}
