#pragma once

#include "LinearProgram.h"
#include "Utils.h"

namespace CP
{
class ToblexSolver
{
public:
    explicit ToblexSolver(LinearProgram& problem)
        : problem(problem)
    {
    }

    bool solve(Vecd& primal, Vecd& dual)
    {
        solStat = NONOPTIMAL;
        // TODO: Check when solution is infeasible (ie Feasible set is empty)

        initTableau();

        if (problem.constraintVector().minCoeff() < 0) {
            // TODO   : Solve Phase 1 Problem to get feasible starting point
            std::cerr << "Phase 1 has not been implemented yet!" << std::endl;
            return false;

            // Setup Phase 1 Tableau, ie add columns for artificial vars (based on number of negative entries in b)
            // rows in which rhs < 0 -> multiply by -1 and add artificial variable (ie also add column)

            // Basis per row is artificial variable or slack variable if no artifical var is present

            // If a basis var appears in the obj. function, remove it using row op.
            removeBasicsInObjective();

            // Solve it using Simplex
            solveSimplex();

            // Once phase 1 is solved, drop columns of artificial variables and reset objective function row
            Matd newTableau;

            // If a basis var appears in the obj. function, remove it using row op.
            removeBasicsInObjective();

        }

        // Solve Simplex
        solveSimplex();

        std::cout << *this << std::endl;

        getPrimalSolution(primal);
        return false;
    }

    void addConstraint(const Vecd& a, const double& d)
    {
        // TODO: When just adding a new constraint, we should probably not reset the whole tableau
        // Make this smarter, ie as few changes as possible to the tableau to solve subsequent iterations of CP faster (hopefully)!

        return; // already handled since this solver refers to the LinearProgram problem.
        /*
        uint n = problem.dimension();
        assert(a.size() == n);

        // Add new row to tableau for ineq. and col for slack var
        uint nr = tableau.rows()+1;
        uint nc = tableau.cols()+1;
        tableau.conservativeResize(nr,nc);

        // move rhs col to the new last col and add rhs value to last entry
        tableau.block(0,nc-1,nr,1) = tableau.block(0,nc-2,nr,1);
        tableau(nr-1,nc-1) = d;

        // set column of new var to unit vector
        tableau.block(0,nc-2,nr,1).fill(0);
        tableau(nr-1,nc-2) = 1;

        // Set last row to new constraint
        tableau(nr-1,0) = 0;
        tableau.block(nr-1,1,1,n) = a.transpose();

        // add new slack var to basis
        B.push_back(n+(nr-1));

        solStat = NONOPTIMAL;
        */
    }

    void write_to_file()
    {
    }

    friend std::ostream& operator<< (std::ostream& out, const ToblexSolver& solver) {
        Vecd p = Vecd(solver.problem.dimension());
        solver.getPrimalSolution(p);
        out << "=================================================================================================" << std::endl;
        out << "Tableau:\n" << solver.tableau << std::endl;
        out << "Basis: " << solver.B << std::endl;
        out << "Status: " << solver.solStat << std::endl;
        out << "Primal Solution: " << p.transpose() << std::endl;
        out << "Optimal Value: " << solver.getOptimalValue() << std::endl;
        out << "=================================================================================================" << std::endl;
        return out;
    }

private:
    enum SolStatus {NONOPTIMAL=0, OPTIMAL=1, UNBOUNDED=2, INFEASIBLE=3};

    LinearProgram& problem;
    Matd tableau; // Simplex Tableau of size (nIneq+1) x (1+nVars+1)
    std::vector<int> B; // Basis Variables (1st entry will always be 0 for the cost variable z)
    SolStatus solStat;

    void initTableau()
    {
        // Setup Phase 2 Tableau
        uint m = problem.nInequalities();
        uint n = problem.dimension();
        tableau = Matd(m+1, 1+n+m+1);
        tableau(0,0) = 1; // 1
        tableau.block(1,0,m,1).fill(0); // 0
        tableau.block(0,1,1,n) = problem.costCoefficients().transpose(); // c_x
        tableau.block(0,n+1,1,m+1).fill(0); // c_slack
        tableau.block(1,1,m,n) = problem.constraintMatrix(); // A
        tableau.block(1,n+1,m,m) = Matd::Identity(m,m); // Id
        tableau.block(1,n+m+1,m,1) = problem.constraintVector(); // rhs b

        // Initial Basis is slack variables
        B.resize(m+1);
        B[0] = 0;
        for (uint i = 1; i < m+1; ++i) {B[i] = n+i;}
    }

    std::pair<int,int> pickPivot()
    {
        // Pick Pivot Column
        // If all coeffs in obj. are nonneg, sol is optimal, return
        int col = -1;
        double min = INFINITY;
        for (uint i = 1; i < tableau.cols()-1; ++i) {
            double t = tableau(0,i);
            if (t < 0 && t < min) {
                min = t;
                col = i;
            }
        }
        if (col == -1) {
            std::cout << "SOLUTION IS OPTIMAL" << std::endl;
            solStat = OPTIMAL;
            return std::make_pair(-1,-1);
        }

        // Pick Pivot Row
        // If all entries in col are nonpos, sol is unbounded, return
        int row = -1;
        min = INFINITY;
        for (uint i = 1; i < tableau.rows(); ++i) {
            double t = tableau(i,col);
            if (t > 0) {
                double b = tableau(i,tableau.cols()-1);
                t = b / t;
                if (t < min) {
                    min = t;
                    row = i;
                }
            }
        }
        if (row == -1) {
            std::cout << "SOLUTION IS UNBOUNDED" << std::endl;
            solStat = UNBOUNDED;
            return std::make_pair(col,-1);
        }

        solStat = NONOPTIMAL;
        return std::make_pair(col,row);
    }

    void doPivot(const int& col, const int& row)
    {
        double p = tableau(row, col);
        assert(p > 0);

        // Divide row by p to make pivot 1
        tableau.row(row) /= p;
        assert(std::abs(tableau(row,col)-1.) < 1e-9);
        tableau(row, col) = 1;

        // Make all other entries in column 0 by using elementary row operations
        for (uint i = 0; i < tableau.rows(); ++i) {
            if (i == row) {continue;}
            tableau.row(i) -= tableau(i,col) * tableau.row(row);
            assert(std::abs(tableau(i,col)) < 1e-9);
            tableau(i,col) = 0;
        }

        // Change Basis
        B[row] = col;
    }

    void solveSimplex()
    {
        // Assert rhs >= 0
        for (auto i = 1; i < tableau.rows(); ++i) {assert(tableau(i,tableau.cols()-1) >= 0);}

        while (solStat != OPTIMAL) {

            auto cr = pickPivot();
            //std::cout << tableau << std::endl;
            //std::cout << "Pivot " << cr.first << ", " << cr.second << std::endl;

            if (solStat != NONOPTIMAL) {
                return;
            }

            doPivot(cr.first, cr.second);
        }
    }

    void getPrimalSolution(Vecd& p) const {
        // TODO: Is it correct that if a variable is not basic (because n > m), it is 0?

        p.fill(0);
        for (auto row = 1u; row < tableau.rows(); ++row) {
            const auto& i = B[row];
            if (i >= 1 && i <= problem.dimension()) {
                p[i-1] = tableau(row, tableau.cols()-1);
            }
        }
    }

    double getOptimalValue() const {return tableau(0, tableau.cols()-1);}

    void removeBasicsInObjective()
    {
        // If a basis var appears in the obj. function, remove it using row op.
        for (auto row = 1u; row < tableau.rows(); ++row) {
            auto col = B[row];
            double c = tableau(0,col);
            if (c != 0) {
                assert(std::abs(tableau(row,col)-1) < 1e-9);
                tableau.row(0) -= c * tableau.row(row);
                assert(tableau(0,col) < 1e-9);
            }
        }
    }

    void resetObjectiveRow()
    {
        // Set the objective row to the original objective function (based on original vars with rhs=0)
        uint n = problem.dimension();
        tableau.block(0,1,1,n) = problem.costCoefficients().transpose(); // c_x
        tableau.block(0,n+1,1,tableau.cols()-n-1).fill(0); // c_slack

        std::cout << "Objective Row: " << tableau.row(0) << std::endl;
    }
};
}
