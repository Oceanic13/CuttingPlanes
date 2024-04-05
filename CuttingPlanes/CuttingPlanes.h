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
    using MILP = MixedIntegerLinearProgram;
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
        solStat = ToblexSolver::NONOPTIMAL;
        cuts_coeffs.clear();
        cuts_rhs.clear();
        uint n = problem.dimension();
        
        Vecd Ai(n);
        double bi;
        Vecd x(n);

        for (uint iter = 0; iter < 100; ++iter) {

            solver.solve();

            if (solver.isInfeasible()) {
                solStat = ToblexSolver::INFEASIBLE;
                break;
            }

            double v = solver.getOptimalValue();
            x = solver.getOptimalSolution();

            simplex_solutions.push_back(x);
            simplex_values.push_back(v);

            //std::cout << solver << std::endl;

            if (solver.generateGomoryMixedIntegerCut(x, Ai, bi)) {
                cuts_coeffs.push_back(Ai);
                cuts_rhs.push_back(bi);
                addInequalityConstraint(Ai, bi);
            } else {
                solStat = ToblexSolver::OPTIMAL;
                optimal_solution = x;
                optimal_value = v;
                break;
            }
        }
    }

    void exportJson(const std::string& filename)
    {
        Json j;
        j["sols"] = simplex_solutions;
        j["cuts_coeffs"] = cuts_coeffs;
        j["cuts_rhs"] = cuts_rhs;

        std::ofstream o(filename);
        o << std::setw(4) << j << std::endl;
        o.close();
    }

    inline const uint numberOfCuts() {return cuts_rhs.size();}

    inline const Vecd& optimalSolution() {return optimal_solution;}

    inline const double optimalValue() {return optimal_value;}

    inline bool isInfeasible() {return solStat == ToblexSolver::INFEASIBLE;}

private:
    ToblexSolver solver;
    ToblexSolver::SolStatus solStat;
    
    MixedIntegerLinearProgram& problem;

    Vecd optimal_solution;
    double optimal_value;
    
    std::vector<Vecd> simplex_solutions;
    std::vector<double> simplex_values;
    std::vector<Vecd> cuts_coeffs;
    std::vector<double> cuts_rhs;

    /**
     * Adds the linear constraint ax >= d to the system
     */
    void addInequalityConstraint(const Vecd& Ai, const double& bi)
    {
        problem.addConstraint(Ai, bi, MILP::LEQ);
        solver.addConstraint(Ai, bi, MILP::LEQ);
    }
};
}
