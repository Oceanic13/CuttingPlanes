#include "CuttingPlanes.h"
#include "soplex/spxdefines.h"
#include <iomanip>
#include <iostream>
#include <fstream>

CuttingPlanes::CuttingPlanes()
{}

void CuttingPlanes::solve(const Vec& c, const Mat& A, const Vec& b, const std::vector<bool>& I)
{
    this->I = I;
    n = c.size();
    m = b.size();
    this->c = c;
    this->A = A;
    this->b = b;
    n_cuts = 0;

    Vec a(n);

    init_simplex(c, A, b);

    for (uint iter = 0; iter < 3; ++iter) {

        Vec primal(n);
        bool success = solve_simplex(primal);
        simplex_solutions.push_back(primal);

        int k = integer_constraint_violation(primal);

        if (k == -1) {
            std::cout << "Found Integer Solution: " << primal.transpose() << std::endl;
            break;
        }

        // TODO: Generate Gomory Cut

        if (iter == 0) {
            a << 0, 1;
            add_constraint(a, 1);
        } else if (iter == 1) {
            a << -1, 1;
            add_constraint(a, 0);
        }
    }

    /* write LP in .lp format */
    simplex.writeFileReal("/Users/tobiaskohler/Uni/CuttingPlanes/soplex.lp", NULL, NULL, NULL);
    simplex.writeBasisFile("/Users/tobiaskohler/Uni/CuttingPlanes/basis.lp", NULL, NULL, false);

    std::cout << "Hello World" << std::endl;
}

void CuttingPlanes::init_simplex(const Vec& c, const Mat& A, const Vec& b)
{
    simplex = soplex::SoPlex();

    /* set the objective sense */
    simplex.setIntParam(soplex::SoPlex::OBJSENSE, soplex::SoPlex::OBJSENSE_MINIMIZE);

    /* we first add variables */
    soplex::DSVector dummycol(0);
    for (uint i = 0; i < n; ++i) {
        simplex.addColReal( soplex::LPCol(c[i], dummycol, soplex::infinity, - soplex::infinity));
    }

    /* then constraints one by one */
    for (uint i = 0; i < m; ++i) {
        soplex::DSVector row(n);
        for (uint j = 0; j < n; ++j) {
            row.add(j, -A(i,j));
        }
        simplex.addRowReal( soplex::LPRow(-b[i], row,  soplex::infinity));
    }
}

bool CuttingPlanes::solve_simplex(Vec& primal)
{

    /* solve LP */
    soplex::SPxSolver::Status stat;
    soplex::DVector p(n);
    stat = simplex.optimize();

    /* get solution */
    if(stat == soplex::SPxSolver::OPTIMAL) {
        simplex.getPrimal(p);
        for (uint i = 0; i < n; ++i) {
            primal[i] = p[i];
        }
        return true;
    } else {
        std::cerr << "Error: SoPlex returned with status " << stat << ".\n";
        return false;
    }
}

void CuttingPlanes::add_constraint(const Vec& a, const double& d)
{
    n_cuts++;
    b.conservativeResize(m + n_cuts);
    A.conservativeResize(m + n_cuts, Eigen::NoChange);
    b[m + n_cuts - 1] = d;
    A.row(m + n_cuts - 1) = a;

    soplex::DSVector row(n);
    for (uint j = 0; j < n; ++j) {
        row.add(j, -a[j]);
    }
    simplex.addRowReal(soplex::LPRow(-d, row, soplex::infinity));
}

int CuttingPlanes::integer_constraint_violation(const Vec& x)
{
    for (uint i = 0; i < n; ++i) {
        if (I[i] && trunc(x[i]) != x[i]) {
            return i;
        }
    }
    return -1;
}

void CuttingPlanes::export_json(const std::string& filename)
{

    json j = {
        {"A", A.rowwise()},
        {"b", b},
        {"c", c},
        {"n", n},
        {"m", m},
        {"I", I},
        {"sols", simplex_solutions}
    };

    std::cout << std::setw(4) << j << std::endl;

    std::ofstream o(filename);
    o << std::setw(4) << j << std::endl;
    o.close();
}
