#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include "Utils.h"

template <class RMILPSolver>
class CuttingPlanes
{
public:
    explicit CuttingPlanes() {}

    /**
     * @brief Solves the Mixed Integer Linear Program min c^T*x s.t. Ax >= b and x >= 0 and x_i integer for all i s.t. I[i]
     * @param c Linear Objective Function of size n
     * @param A Inequality lhs (Ax >= b) of size m x n
     * @param b Inequality rhs (Ax >= b) of size m
     * @param I Integer Constraints of size n
     */
    void solve(const Vecd& c, const Matd& A, const Vecd& b, const std::vector<bool>& I)
    {
        this->I = I;
        this->n = c.size();
        this->m = b.size();
        this->c = c;
        this->A = A;
        this->b = b;
        n_cuts = 0;
        assert(A.rows() == m && A.cols() == n);
        assert(I.size() == n);
        
        Vecd a(n);
        Vecd primal(n);
        Vecd dual(1);

        solver.init(c, A, b);

        for (uint iter = 0; iter < 3; ++iter) {
            


            bool success = solver.solve(primal, dual);
            simplex_solutions.push_back(primal);


            int k = integer_constraint_violation(primal);
            if (k == -1) {break;}

            // TODO: Generate Gomory Cut

            if (iter == 0) {
                a << 0, -1;
                add_constraint(a, -1);
            } else if (iter == 1) {
                a << 1, -1;
                add_constraint(a, 0);
            }
        }

        solver.write_to_file();

        std::cout << "Hello World" << std::endl;
    }

    void export_json(const std::string& filename)
    {
        Json j = {
            {"A", A.rowwise()},
            {"b", b},
            {"c", c},
            {"n", n},
            {"m", m},
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
    
    Vecd c;
    Matd A;
    Vecd b;
    uint n; // dimensionality
    uint m; // number of (original) inequalities
    uint n_cuts;
    std::vector<bool> I; // Integer constraints
    
    std::vector<Vecd> simplex_solutions;

    /**
     * Returns the first variable index where the current solution violates an integer constraint
     * or -1 if the solution does not violate any integer constraints
     */
    int integer_constraint_violation(const Vecd& x)
    {
        for (uint i = 0; i < n; ++i) {
            if (I[i] && !isInt(x[i])) {
                return i;
            }
        }
        return -1;
    }

    /**
     * Adds the linear constraint ax >= d to the system
     */
    void add_constraint(const Vecd& a, const double& d)
    {
        n_cuts++;
        b.conservativeResize(m + n_cuts);
        A.conservativeResize(m + n_cuts, Eigen::NoChange);
        b[m + n_cuts - 1] = d;
        A.row(m + n_cuts - 1) = a;

        solver.add_constraint(a, d);
    }
};
