#include "CuttingPlanes.h"
#include <iomanip>
#include <iostream>
#include <fstream>
#include "SoPlexSolver.cpp";

CuttingPlanes::CuttingPlanes(const Vec& c, const Mat& A, const Vec& b, const std::vector<bool>& I)
    : n(c.size()), m(b.size()), c(c), A(A), b(b), I(I)
{}

void CuttingPlanes::solve()
{
    Vec x(n);

    //Vec a(2); a << 0, 1;
    //add_constraint(a, 1);

    auto s = soplex::SoPlexSolver();
    s.solve(c, A, b);


    //Vec t(2); t << 1, 1.5;
    //simplex_solutions.push_back(t);

    std::cout << "Hello World" << std::endl;
}

void CuttingPlanes::export_json(const std::string& filename)
{
    json j = {
        {"n", n},
        {"m", m},
        {"c", c},
        {"I", I},
        {"simplex", simplex_solutions},
        {"A", A.rowwise()},
        {"b", b}
    };

    std::cout << std::setw(4) << j << std::endl;

    std::ofstream o(filename);
    o << std::setw(4) << j << std::endl;
    o.close();
}
