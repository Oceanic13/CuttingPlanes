#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include "LinearProgram.h"
#include "Utils.h"

namespace CP
{
template <class RMILPSolver>
class CuttingPlanes
{
public:
    explicit CuttingPlanes(LinearProgram& problem) : problem(problem), solver(RMILPSolver(problem)) {

    }

    /**
     * @brief Solves the Mixed Integer Linear Program min c^T*x s.t. Ax >= b and x >= 0 and x_i integer for all i s.t. I[i]
     * @param c Linear Objective Function of size n
     * @param A Inequality lhs (Ax >= b) of size m x n
     * @param b Inequality rhs (Ax >= b) of size m
     * @param I Integer Constraints of size n
     */
    void solve(const std::vector<bool>& I)
    {
        uint n = problem.dimension();
        this->I = I;
        n_cuts = 0;
        assert(I.size() == n);
        
        Vecd a(n);
        Vecd primal(n);

        for (uint iter = 0; iter < 3; ++iter) {
            
            std::cout << "CUTTING PLANES ITERATION " << iter << std::endl;

            bool success = solver.solve(primal);
            simplex_solutions.push_back(primal);


            int k = integer_constraint_violation(primal);
            if (k == -1) {break;}

            // TODO: Generate Gomory Cut

            if (iter == 0) {
                a << 0, 1;
                addConstraint(a, 1);
            } else if (iter == 1) {
                a << -1, 1;
                addConstraint(a, 0);
            }
        }

        solver.write_to_file();
    }

    void export_json(const std::string& filename)
    {
        Json j = {
            {"A", problem.constraintMatrix().rowwise()},
            {"b", problem.constraintVector()},
            {"c", problem.costCoefficients()},
            {"n", problem.dimension()},
            {"m", problem.nInequalities() - n_cuts},
            {"I", I},
            {"sols", simplex_solutions}
        };

        //std::cout << std::setw(4) << j << std::endl;

        std::ofstream o(filename);
        o << std::setw(4) << j << std::endl;
        o.close();
    }

private:
    RMILPSolver solver;
    
    LinearProgram& problem;
    uint n_cuts;
    std::vector<bool> I; // Integer constraints
    
    std::vector<Vecd> simplex_solutions;

    /**
     * Returns the first variable index where the current solution violates an integer constraint
     * or -1 if the solution does not violate any integer constraints
     */
    int integer_constraint_violation(const Vecd& x)
    {
        for (uint i = 0; i < x.size(); ++i) {
            if (I[i] && !isInt(x[i])) {
                return i;
            }
        }
        return -1;
    }

    /**
     * Adds the linear constraint ax >= d to the system
     */
    void addConstraint(const Vecd& a, const double& d)
    {
        problem.addConstraint(a, d);
        solver.addConstraint(a, d);
    }
};
}
