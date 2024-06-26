#pragma once

#include "MixedIntegerLinearProgram.h"
#include "Utils.h"
#include <fstream>

namespace CP
{
using MILP = MixedIntegerLinearProgram;

// TODO: Move generateGomoryCut to Cutting Planes Solver

// TODO: Sparse Matrix instead of dense for tableau?

// TODO: Don't start from 0 again in each Cutting Planes iteration

// TODO: Smarter Cutting Planes generation: Don't just pick a random non-integer variable
// generally more control: Select Type of cut, e.g. specific for binary program?

// TODO: Handle unbounded problems correctly (for cutting planes)

/**
 * Two Phase Primal Simplex Algorithm to solve Linear Programs of the general form
 * minimize c*x
 * s.t. Ax <= b and Bx = d and x >= 0
 */
class ToblexSolver
{
    friend class CuttingPlanesSolver;
    friend class BranchBoundSolver;

public:
    enum SolStatus {NONOPTIMAL=0, OPTIMAL=1, UNBOUNDED=2, INFEASIBLE=3};

    explicit ToblexSolver(MILP& problem) :
        x(Vecd(problem.dimension())),
        nProblemVars(problem.dimension()),
        problem(problem)
    {

    }

    void solve()
    {
        initTableau();
        initialT = T;
        initialBasis = basis;

        // PHASE 1 SIMPLEX
        if (nArtificialVars > 0) { // Check if we need to solve the Phase 1 Problem first to find a feasible starting point other than x=0
            T.row(0).fill(0);
            T.row(0).segment(problem.dimension()+nSlackVars, nArtificialVars).fill(-1); // Set Sum of Artificial Variables as Function to Minimize
            removeBasicsInObjective(); // If a basis var (artificial) appears in the obj. function, remove it using row op.
            solveSimplex();
            if (!artificialVariablesAreAllZero()) {
                std::cerr << "PROBLEM IS INFEASIBLE" << std::endl;
                solStat = INFEASIBLE;
                return;
            }
            removeArtificialColumns();
            setProblemObjective(problem.costCoefficients());
            removeBasicsInObjective(); // If a basis var (decision) appears in the obj. function, remove it using row op.
        }

        // PHASE 2 SIMPLEX
        solveSimplex();

        // Get optimal solution
        x.setZero();
        for (auto r = 1u; r < T.rows(); ++r) {
            if (basis[r-1] < problem.dimension()) {
                x[basis[r-1]] = T(r,T.cols()-1);
            }
        }

        return;
    }


    ///
    /// \brief generateGomoryMixedIntegerCut
    /// \param x Current optimal solution
    /// \param Ai lhs coefficient of cut
    /// \param bi rhs constant of cut
    /// \return true iff cut was successfully generated
    ///
    bool generateGomoryMixedIntegerCut(const Vecd& x, Vecd& Ai, double& bi)
    {
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
        for (uint j = n1; j < n1+n2+nSlackVars; ++j) {
            double aij = T(row,j);
            if (aij >= 0) {continue;}

            double coeff = aij / (1. - fi);

            if (j < n) { // N2- problem var
                Ai[j] = coeff;
            } else { // N2- slack var
                Vecd varCoeffs(n); double constTerm;
                getSlackVarCoeffs(j, varCoeffs, constTerm); // resubstitute
                //std::cout << "VarCoeffs for slack " << j << " are " << varCoeffs.transpose() << " + " << constTerm << std::endl;
                Ai += coeff * varCoeffs;
                bi -= coeff * constTerm;
            }
        }

        for (uint j = 0; j < Ai.size(); ++j) {if (isInt(Ai[j])) {Ai[j] = std::round(Ai[j]);}}
        if (isInt(bi)) {bi = std::round(bi);}

        return true;
    }

    const Vecd& getOptimalSolution() const {return x;}
    double getOptimalValue() const {assert(isOptimal()); return T(0, T.cols()-1);}

    bool isUnbounded() const {return solStat==UNBOUNDED;}
    bool isOptimal() const {return solStat==OPTIMAL;}
    bool isInfeasible() const {return solStat==INFEASIBLE;}

    const Matd& tableau() {return T;}

    void addConstraint(const Vecd& Ai, const double& bi, const MILP::ConstraintType type = MILP::LEQ) {
        if (type == MILP::GEQ) {addConstraint(-Ai, -bi, MILP::LEQ);}

        //     // For each inequality constraint <= b with b >= 0
        //     // add a slack variable
        //     // e.g. 2x + y <= 2 becomes 2x + y + s = 2

        //     // For each inequality constraint <= b with b < 0
        //     // multiply by -1, subtract a slack variable and add an artificial variable
        //     // e.g. 2x - y <= -3 becomes -2x + y - s + a = 3

        //     // For each equality constraint = d
        //     // multiply the equation by -1 if d < 0 and add an artificial variable if d >= 0
        //     // e.g. 3x - 2y = 4 becomes 3x - 2y + a = 4
        //     // 2x + y -7z = -15 becomes -2x - y + 7z + a = 15

        uint n = Ai.size();
        uint nRows = T.rows()+1;
        T.conservativeResize(nRows, Eigen::NoChange);

        T.row(nRows-1).segment(0, n) = (bi>=0)? Ai : -Ai;
        T(nRows-1, T.cols()-1) = (bi>=0)? bi : -bi;

        if (type == MILP::EQ) {
            addArtificialCol(1);
        } else {
            assert(type == MILP::LEQ);
            if (bi >= 0) {
                addSlackCol(1);
            } else {
                addSlackCol(-1);
                addArtificialCol(1);
            }
        }
    }

private:
    Matd T;
    std::vector<int> basis;
    Matd initialT;
    std::vector<int> initialBasis;
    MILP& problem;
    uint nProblemVars;
    uint nSlackVars;
    uint nArtificialVars;
    Vecd x;
    SolStatus solStat;


    void initTableau()
    {
        nSlackVars = nArtificialVars = 0;
        const Vecd& c = problem.costCoefficients();
        const Matd& A = problem.inequaltyMatrix();
        const Vecd& b = problem.inequalityVector();
        const Matd& B = problem.equalityMatrix();
        const Vecd& d = problem.equalityVector();

        uint n = problem.dimension();
        uint m = problem.nInequalities();
        uint p = problem.nEqualities();

        basis.resize(m+p);

        T = Matd(1, n+1);
        setProblemObjective(c);
        for (uint i = 0; i < m; ++i) {addConstraint(A.row(i), b[i], MILP::LEQ);}
        for (uint i = 0; i < p; ++i) {addConstraint(B.row(i), d[i], MILP::EQ);}
    }

    void solveSimplex()
    {
        assert((T.col(T.cols()-1).segment(1,T.rows()-1).array() >= 0).all()); // Assert rhs >= 0
        solStat = NONOPTIMAL;

        int row;
        int col;
        while (solStat != OPTIMAL) {
            pickPivot(col, row);
            if (solStat != NONOPTIMAL) {return;}
            assert(T(row,col) > 0);
            gaussianPivot(T, row, col);

            // Change Basis
            basis[row-1] = col;
        }
    }

    void pickPivot(int& col, int& row)
    {
        // Pick Pivot Column
        // If all coeffs in obj. are nonpos, sol is optimal, return
        col = -1;
        double max = -INFINITY;
        for (uint i = 0; i < T.cols()-1; ++i) {
            double t = T(0,i);
            if (t > 0 && t > max) {
                max = t;
                col = i;
            }
        }
        if (col == -1) {
            solStat = OPTIMAL;
            return;
        }

        // Pick Pivot Row
        // If all entries in col are nonpos, sol is unbounded, return
        row = -1;
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
            solStat = UNBOUNDED;
            return;
        }

        solStat = NONOPTIMAL;
        return;
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
                //T.row(0) -= c * T.row(row);
                addRowToRow(T, 0, -c, row);
                assert(T(0,col) < 1e-9);
            }
        }
    }

    bool artificialVariablesAreAllZero()
    {
        uint minArtIdx = problem.dimension() + nSlackVars;
        uint maxArtIdx = minArtIdx + nArtificialVars;
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
        removeCols(T, problem.dimension() + nSlackVars, nArtificialVars);
        nArtificialVars = 0;
    }

    /**
     * Returns a slack variable s_j expressed as sum c_ij x_i + constTerm
     */
    void getSlackVarCoeffs(const int& j, Vecd& varCoeffs, double& constTerm) {
        uint n = problem.dimension();
        assert(j >= n && j < n+nSlackVars);
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

    // Inserts a slack column and sets the last enrty to s and sets the basis of the last row to the new column
    void addSlackCol(const int& s)
    {
        assert(s == -1 || s == 1);
        int r = T.rows()-1;
        assert(r >= 1);
        int c = problem.dimension() + nSlackVars++;
        insertEmptyCol(T, c);
        T(r, c) = s;

        // Correct column refs for artificial basics
        for (uint i = 0; i < basis.size(); i++) {
            if (basis[i] >= c) {
                basis[i]++;
            }
        }

        basis[r-1]= c;
    }

    // Inserts a artificial column and sets the last enrty to a and sets the basis of the last row to the new column
    void addArtificialCol(const int& a)
    {
        assert(a == 1);
        int r = T.rows()-1;
        assert(r >= 1);
        int c = problem.dimension() + nSlackVars + nArtificialVars++;
        assert(c == T.cols()-1);
        insertEmptyCol(T, c);
        T(r, c) = a;
        basis[r-1]= c;
    }

    void setProblemObjective(const Vecd& c)
    {
        assert(c.size() == problem.dimension());
        T.row(0).fill(0);
        T.row(0).segment(0,c.size()) = -c;
    }

    double& cj(const int& j) {return T(0,j);}
    double& Arj(const int& r, const int& j) {return T(r,j);}
    double& br(const int& r) {return T(r,T.cols()-1);}

    friend std::ostream& operator<< (std::ostream& out, const ToblexSolver& s) {
        out << "=================================================================================================" << std::endl;
        out << "Initial Tableau:\n" << s.initialT << std::endl;
        out << "Initial Basis: " << s.initialBasis << std::endl;
        out << "Tableau:\n" << s.T << std::endl;
        out << "nVars, nSlack, nArtificial: " << s.problem.dimension() << ", " << s.nSlackVars << ", " << s.nArtificialVars << std::endl;
        out << "Basis: " << s.basis << std::endl;
        out << "Status: " << s.solStat << std::endl;
        out << "Optimal Solution: " << s.getOptimalSolution().transpose() << std::endl;
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

class CuttingPlanesHistory {
    friend class CuttingPlanesSolver;
public:
    explicit CuttingPlanesHistory() {}

    void push_back_solution(Vecd x, double v)
    {
        simplex_solutions.push_back(x);
        simplex_values.push_back(v);
    }

    void push_back_cut(Vecd a, double b)
    {
        cuts_coeffs.push_back(a);
        cuts_rhs.push_back(b);
    }

    void clear()
    {
        simplex_solutions.clear();
        simplex_values.clear();
        cuts_coeffs.clear();
        cuts_rhs.clear();
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

    int n_cuts() {return cuts_coeffs.size();}

private:
    std::vector<Vecd> simplex_solutions;
    std::vector<double> simplex_values;
    std::vector<Vecd> cuts_coeffs;
    std::vector<double> cuts_rhs;
};

class CuttingPlanesSolver
{
public:
    explicit CuttingPlanesSolver(MixedIntegerLinearProgram& problem) : problem(problem), toblex(ToblexSolver(problem)) {

    }

    /// Solves the Mixed Integer Linear Program using Cutting Planes
    void solve(int max_iters=100)
    {
        solStat = ToblexSolver::NONOPTIMAL;
        history.clear();
        uint n = problem.dimension();

        Vecd Ai(n);
        double bi;
        Vecd x(n);

        for (uint iter = 0; iter < max_iters; ++iter) {

            toblex.solve();

            if (toblex.isInfeasible()) {
                solStat = ToblexSolver::INFEASIBLE;
                break;
            }

            double v = toblex.getOptimalValue();
            x = toblex.getOptimalSolution();

            std::cout << "Simplex: " << x.transpose() << " with value " << v << std::endl;

            history.push_back_solution(x, v);

            if (toblex.generateGomoryMixedIntegerCut(x, Ai, bi)) {
                history.push_back_cut(Ai, bi);
                addInequalityConstraint(Ai, bi);
                std::cout << "Add Cut " << Ai.transpose() << " <= " << bi << std::endl;
            } else {
                solStat = ToblexSolver::OPTIMAL;
                optimal_solution = x;
                optimal_value = v;
                std::cout << "Finished: " << x.transpose() << " with value " << v << std::endl;
                break;
            }
        }
    }


    CuttingPlanesHistory& getHistory()
    {
        return history;
    }

    inline const Vecd& optimalSolution() {return optimal_solution;}

    inline const double optimalValue() {return optimal_value;}

    inline bool isInfeasible() {return solStat == ToblexSolver::INFEASIBLE;}

private:
    MILP& problem;
    ToblexSolver toblex;
    ToblexSolver::SolStatus solStat;

    Vecd optimal_solution;
    double optimal_value;
    CuttingPlanesHistory history;

    /**
     * Adds the linear constraint ax <= b to the system
     */
    void addInequalityConstraint(const Vecd& Ai, const double& bi)
    {
        problem.addConstraint(Ai, bi, MILP::LEQ);
        toblex.addConstraint(Ai, bi, MILP::LEQ);
    }
};


// TODO: Sometime
class BranchBoundSolver
{
public:
    explicit BranchBoundSolver(MixedIntegerLinearProgram& problem)
    {
    }
};

}
