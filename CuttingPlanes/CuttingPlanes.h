#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include "MixedIntegerLinearProgram.h"
#include "ToblexSolver.h"
#include "Utils.h"

namespace CP
{
class CuttingPlanes
{
public:
    explicit CuttingPlanes(MixedIntegerLinearProgram& problem) : problem(problem), solver(ToblexSolver(problem)) {

    }

    /**
     * @brief Solves the Mixed Integer Linear Program min c^T*x s.t. Ax >= b and x >= 0 and x_i integer for all i s.t. I[i]
     * @param c Linear Objective Function of size n
     * @param A Inequality lhs (Ax >= b) of size m x n
     * @param b Inequality rhs (Ax >= b) of size m
     * @param I Integer Constraints of size n
     */
    void solve()
    {
        uint n = problem.dimension();
        n_cuts = 0;
        
        Vecd Ai(n);
        double bi;
        Vecd x(n);

        for (uint iter = 0; iter < 2; ++iter) {
            
            std::cout << "CUTTING PLANES ITERATION " << iter << std::endl;

            bool success = solver.solve(x);
            simplex_solutions.push_back(x);

            std::cout << solver << std::endl;

            if (solver.generateGomoryMixedIntegerCut(x, Ai, bi)) {
                addInequalityConstraint(Ai, bi);
            } else {
                break;
            }
            continue;


            // TODO: Generate Gomory Cut
        /*
            if (iter == 0) {
                a << 0, 1;
                addInequalityConstraint(a, 1);
            } else if (iter == 1) {
                a << -1, 1;
                addInequalityConstraint(a, 0);
            }*/
        }

        //solver.write_to_file();
    }

    void export_json(const std::string& filename)
    {
        Json j = {
            {"A", problem.inequaltyMatrix().rowwise()},
            {"b", problem.inequalityVector()},
            {"B", problem.equaltyMatrix().rowwise()},
            {"d", problem.equalityVector()},
            {"c", problem.costCoefficients()},
            {"nCuts", problem.dimension()},
            {"n1", problem.nIntegerConstraints()},
            {"sols", simplex_solutions}
        };

        //std::cout << std::setw(4) << j << std::endl;

        std::ofstream o(filename);
        o << std::setw(4) << j << std::endl;
        o.close();
    }

private:
    ToblexSolver solver;
    
    MixedIntegerLinearProgram& problem;
    uint n_cuts;
    
    std::vector<Vecd> simplex_solutions;

    /**
     * Adds the linear constraint ax >= d to the system
     */
    void addInequalityConstraint(const Vecd& Ai, const double& bi)
    {
        problem.addInequalityConstraint(Ai, bi);
        solver.addInequalityConstraint(Ai, bi);
    }
};
}
