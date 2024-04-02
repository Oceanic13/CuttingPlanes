#pragma once

#include "MixedIntegerLinearProgram.h"
#include "Utils.h"

namespace CP
{
using MILP = MixedIntegerLinearProgram;

// TODO: Much better to have a sparse Matrix for the Tableau...
// TODO: We probably do not need to reset the entire tableau everytime
// Set the Tableau to the initial tableau if it exists and then manually add the new inequality to it

/**
 * My own Simplex Solver to have something that works and that I understand!
 */
class ToblexSolver
{
public:
    enum SolStatus {NONOPTIMAL=0, OPTIMAL=1, UNBOUNDED=2, INFEASIBLE=3};

    explicit ToblexSolver(MILP& problem) : problem(problem)
    {
    }

    bool solve(Vecd& primal)
    {
        solStat = NONOPTIMAL;
        initTableau();

        //std::cout << "Initial:\n" << *this << std::endl;

        if (nArts > 0) {// Solve Phase 1 Problem, ie it means x=0 is not a feasible starting point

            // Minimize sum of artificial variables
            T.row(0).fill(0);
            for (uint j = problem.dimension() + nSlacks; j < T.cols()-1; ++j) {
                T(0,j) = -1;
            }

            // If a basis var (a) appears in the obj. function, remove it using row op.
            removeBasicsInObjective();

            // Solve it using Simplex
            solveSimplex();

            // If values of artificial vars are not all 0, problem is infeasible
            if (!artificialVariablesAreAllZero()) {
                std::cerr << "PROBLEM IS INFEASIBLE" << std::endl;
                solStat = INFEASIBLE;
                return false;
            }

            // Once phase 1 is solved, drop columns of artificial variables and reset objective function row
            removeArtificialColumns();
            T.row(0).fill(0);
            T.block(0,0,1,problem.dimension()) = - problem.costCoefficients().transpose();

            // If a basis var (x) appears in the obj. function, remove it using row op.
            removeBasicsInObjective();

            solStat = NONOPTIMAL;
            //std::cout << "After Phase 1:\n" << *this << std::endl;
        }

        initialT = T;
        initialBasis = basis;

        solveSimplex();

        //std::cout << "After Phase 2:\n" << *this << std::endl;

        getPrimalSolution(primal);

        return true;
    }

    bool generateGomoryMixedIntegerCut(const Vecd& x, Vecd& Ai, double& bi)
    {
        // TODO: MAke right
        uint n = problem.dimension();
        uint n1 = problem.nIntegerConstraints();
        uint n2 = n-n1;

        // Find Integer violation (column)
        int col = -1;
        for (auto i = 0u; i < n1; ++i) {
            assert(problem.isIntegerConstrained(i));
            if (!isInt(x[i])) {
                col = i;
                break;
            }
        }
        if (col == -1) {return false;}

        //std::cout << "Variable " << col << " violates Integer Constraint!" << std::endl;

        // Find corresponding Basis row
        int row = -1;
        for (uint r = 1; r <= basis.size(); ++r) {
            if (basis[r-1] == col) {
                row = r;
                break;
            }
        }
        assert(row != -1);
        assert(T(row,col)==1);
        //std::cout << "In row " << row << std::endl;

        // Get row and rhs
        bi = br(row);
        double fi = getFrac(bi);
        assert(fi > 0 && fi < 1);
        bi = floor(bi);
        //Vecd varCoeffs(n);

        // Mixed Integer Rounding
        for (uint j = 0; j < n1; ++j) {
            double aij = T(row,j);
            double fij = getFrac(aij);
            if (fij <= fi) { // N1<=
                Ai[j] = floor(aij);
            } else { // N1>
                Ai[j] = floor(aij) + ((fij - fi) / (1. - fi));
            }
        }
        for (uint j = n1; j < n1+n2+nSlacks; ++j) {
            double aij = T(row,j);
            if (aij >= 0) {continue;}

            double coeff = aij / (1. - fi);

            if (j < n) { // N2- problem var
                Ai[j] = coeff;
            } else { // N2- slack var
                Vecd varCoeffs(n); double constTerm;
                getSlackVarCoeffs(j, varCoeffs, constTerm);
                //std::cout << "VarCoeffs for slack " << j << " are " << varCoeffs.transpose() << " + " << constTerm << std::endl;
                Ai += coeff * varCoeffs;
                bi -= coeff * constTerm;
            }
        }

        return true;
    }

    void getPrimalSolution(Vecd& p) const
    {
        // TODO: Is it correct that if a variable is not basic (because n > m), it is 0?

        p.fill(0);
        for (auto row = 1u; row < T.rows(); ++row) {
            const auto& i = basis[row-1];
            if (i >= 0 && i < problem.dimension()) {
                p[i] = T(row, T.cols()-1);
            }
        }
    }

    double getOptimalValue() const {
        return T(0, T.cols()-1);
    }

    const Matd& tableau() {return T;}

    bool isUnbounded() {return solStat==UNBOUNDED;}
    bool isOptimal() {return solStat==OPTIMAL;}
    bool isInfeasible() {return solStat==INFEASIBLE;}

    void addInequalityConstraint(const Vecd& Ai, const double& bi) {
        // TODO: Insteaf of rebuilding the entire tableau, store the initial tableau and add a new row there
    }

    void addEqualityConstraint(const Vecd& Bi, const double& di) {
        // TODO: Insteaf of rebuilding the entire tableau, store the initial tableau and add a new row there
    }

private:
    Matd T;
    Matd initialT;
    std::vector<int> initialBasis;

    MILP& problem;
    std::vector<int> basis;
    //std::vector<int> basisRow;
    uint nSlacks;
    uint nArts;
    SolStatus solStat;


    void initTableau()
    {
        const Vecd& c = problem.costCoefficients();
        const Matd& A = problem.inequaltyMatrix();
        const Vecd& b = problem.inequalityVector();
        const Matd& B = problem.equaltyMatrix();
        const Vecd& d = problem.equalityVector();

        uint n = problem.dimension();
        uint m = problem.nInequalities();
        uint p = problem.nEqualities();

        Matd A2(m, n);
        Vecd b2(m);
        Matd B2(p, n);
        Vecd d2(p);
        std::vector<Triplet> slacks; // col. rel. to 1st slack col.
        std::vector<Triplet> arts; // col. rel. to 1st artificial col.
        uint nRows = 1 + m + p;
        basis.resize(nRows-1);

        //std::cout << "Initializing Inequality rows..." << std::endl;

        for (uint i = 0; i < A.rows(); ++i) {
            // For each inequality constraint <= b with b >= 0
            // add a slack variable
            // e.g. 2x + y <= 2 becomes 2x + y + s = 2
            // For each inequality constraint <= b with b < 0
            // multiply by -1, subtract a slack variable and add an artificial variable
            // e.g. 2x - y <= -3 becomes -2x + y - s + a = 3
            if (b[i] >= 0) {
                A2.row(i) = A.row(i);
                slacks.emplace_back(i,slacks.size(), 1);
                b2[i] = b[i];
                basis[i]= slacks.size();
            } else {
                A2.row(i) = -A.row(i);
                slacks.emplace_back(i,slacks.size(), -1);
                arts.emplace_back(i,arts.size(), 1);
                b2[i] = -b[i];
                basis[i]= -arts.size();
            }
        }

        //std::cout << "Initializing Equality rows..." << std::endl;

        for (uint i = 0; i < B.rows(); ++i) {
            // For each equality constraint = d
            // multiply the equation by -1 if d < 0 and add an artificial variable if d >= 0
            // e.g. 3x - 2y = 4 becomes 3x - 2y + a = 4
            // 2x + y -7z = -15 becomes -2x - y + 7z + a = 15
            if (d[i] >= 0) {
                B2.row(i) = B.row(i);
                arts.emplace_back(i+m,arts.size(), 1);
                d2[i] = d[i];
                basis[i+m] = -arts.size();
            } else {
                B2.row(i) = -B.row(i);
                arts.emplace_back(i+m,arts.size(), 1);
                d2[i] = -d[i];
                basis[i+m] = -arts.size();
            }
        }

        //std::cout << "Initializing Columns..." << std::endl;

        // Setup Tableau
        nSlacks = slacks.size();
        nArts = arts.size();
        uint nCols = n + nSlacks + nArts + 1;
        T = Matd(nRows, nCols);
        //std::cout << "Tableau Size " << nRows << " x " << nCols << std::endl;
        T.fill(0);
        T.block(0,0,1,n) = - c.transpose();
        T.block(1,0,m,n) = A2;
        T.block(m+1,0,p,n) = B2;
        for (const Triplet& t : slacks) {T(t.row()+1,t.col()+n) = t.value();}
        for (const Triplet& t : arts) {T(t.row()+1,t.col()+n+nSlacks) = t.value();}
        T.block(1,nCols-1,m,1) = b2;
        T.block(1+m,nCols-1,p,1) = d2;

        // Initial Basis is artificial variable if it is present or else slack variable
        for (uint i = 0; i < basis.size(); ++i) {
            if (basis[i] > 0) { // basis is slack
                basis[i] += n;
            } else {assert(basis[i] < 0); // basis is artificial
                basis[i] = -basis[i] + n + nSlacks;
            }
            basis[i] -= 1;
        }
        //std::cout << "Basis " << basis << std::endl;
    }

    void solveSimplex()
    {
        // Assert rhs >= 0
        for (auto i = 1; i < T.rows(); ++i) {assert(T(i,T.cols()-1) >= 0);}

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

    std::pair<int,int> pickPivot()
    {
        // Pick Pivot Column
        // If all coeffs in obj. are nonpos, sol is optimal, return
        int col = -1;
        double max = -INFINITY;
        for (uint i = 0; i < T.cols()-1; ++i) {
            double t = T(0,i);
            if (t > 0 && t > max) {
                max = t;
                col = i;
            }
        }
        if (col == -1) {
            //std::cout << *this << std::endl;
            //std::cout << "SOLUTION IS OPTIMAL" << std::endl;
            solStat = OPTIMAL;
            return std::make_pair(-1,-1);
        }

        // Pick Pivot Row
        // If all entries in col are nonpos, sol is unbounded, return
        int row = -1;
        double min = INFINITY;
        for (uint i = 1; i < T.rows(); ++i) {
            double t = T(i,col);
            if (t > 0) {
                double b = br(i);
                t = b / t;
                if (t < min) {
                    min = t;
                    row = i;
                }
            }
        }
        if (row == -1) {
            //std::cout << "SOLUTION IS UNBOUNDED" << std::endl;
            solStat = UNBOUNDED;
            return std::make_pair(col,-1);
        }

        solStat = NONOPTIMAL;
        return std::make_pair(col,row);
    }

    void doPivot(const int& col, const int& row)
    {
        double p = T(row, col);
        assert(p > 0);

        // Divide row by p to make pivot 1
        T.row(row) /= p;
        assert(std::abs(T(row,col)-1.) < 1e-9);
        T(row, col) = 1;

        // Make all other entries in column 0 by using elementary row operations
        for (uint i = 0; i < T.rows(); ++i) {
            if (i == row) {continue;}
            T.row(i) -= T(i,col) * T.row(row);
            assert(std::abs(T(i,col)) < 1e-9);
            T(i,col) = 0;
        }

        // Change Basis
        basis[row-1] = col;
    }

    void removeBasicsInObjective()
    {
        // If a basis var appears in the obj. function, remove it using row op.
        for (auto row = 1u; row < T.rows(); ++row) {
            auto col = basis[row-1];
            //std::cout << "Basis of row " << row <<" is " << col << std::endl;
            double c = T(0,col);
            if (c != 0) {
                assert(std::abs(T(row,col)-1) < 1e-9);
                T.row(0) -= c * T.row(row);
                assert(T(0,col) < 1e-9);
            }
        }
    }

    bool artificialVariablesAreAllZero()
    {
        uint minArtIdx = problem.dimension() + nSlacks;
        uint maxArtIdx = minArtIdx + nArts;
        for (auto r = 1u; r < T.rows(); ++r) {
            auto c = basis[r-1];
            if (c >= minArtIdx && c < maxArtIdx) {
                if (br(r) != 0) {
                    return false;
                }
            }
        }
        return true;
    }

    void removeArtificialColumns() {
        removeCols(T, problem.dimension() + nSlacks, nArts);
        nArts = 0;
    }

    /**
     * Returns a slack variable s_j expressed as sum c_ij x_i + constTerm
     */
    void getSlackVarCoeffs(const int& j, Vecd& varCoeffs, double& constTerm) {
        uint n = problem.dimension();
        assert(j >= n && j < n+nSlacks);
        for (uint r = 1; r <= initialBasis.size(); ++r) {
            if (initialBasis[r-1] == j) {
                varCoeffs = - initialT.row(r).segment(0, n);
                constTerm = initialT(r, initialT.cols()-1);
                return;
            }
        }
        varCoeffs.setZero();
        constTerm = 0;
        return;
    }

    double& cj(const int& j) {return T(0,j);}
    double& Arj(const int& r, const int& j) {return T(r,j);}
    double& br(const int& r) {return T(r,T.cols()-1);}

    friend std::ostream& operator<< (std::ostream& out, const ToblexSolver& s) {
        Vecd p = Vecd(s.problem.dimension());
        s.getPrimalSolution(p);
        out << "=================================================================================================" << std::endl;
        out << "Initial Tableau:\n" << s.initialT << std::endl;
        out << "Initial Basis: " << s.initialBasis << std::endl;
        out << "Tableau:\n" << s.T << std::endl;
        out << "nVars, nSlack, nArtificial: " << s.problem.dimension() << ", " << s.nSlacks << ", " << s.nArts << std::endl;
        out << "Basis: " << s.basis << std::endl;
        out << "Status: " << s.solStat << std::endl;
        out << "Primal Solution: " << p.transpose() << std::endl;
        out << "Optimal Value: " << s.getOptimalValue() << std::endl;
        out << "=================================================================================================" << std::endl;
        return out;
    }

    friend std::ostream& operator<< (std::ostream& out, const SolStatus& s) {
        if (s == NONOPTIMAL) return out << "Nonoptimal";
        if (s == OPTIMAL) return out << "Optimal";
        if (s == UNBOUNDED) return out << "Unbounded";
        if (s == INFEASIBLE) return out << "Infeasible";
        return out;
    }
};
}
