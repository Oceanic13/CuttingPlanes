#pragma once

#include "Utils.h"
#include <fstream>

namespace CP
{
class MixedIntegerLinearProgram
{
public:

    enum ConstraintType {LEQ=0, GEQ=1, EQ=2};

    /**
     * @brief MixedIntegerLinearProgram min c^Tx s.t. Ax <= b and Bx = d and x >= 0 and x_i integer for i = 1,...,n1 <= n
     * @param c Cost Coefficients
     * @param A Inequality Constraint Matrix
     * @param b Inequality Constraint Vector
     * @param n1 Number of integer constraints
     */
    explicit MixedIntegerLinearProgram(const Vecd& c, const Matd& A, const Vecd& b, const Matd& B, const Vecd& d, const int& n1)
        : n(c.size()), m(b.size()), p(d.size()), c(c), A(A), b(b), B(B), d(d), n1(n1)
    {
        assert(A.rows() == m && A.cols() == n);
        assert(B.rows() == p && B.cols() == n);
        assert(n1 <= n);
    }

    MixedIntegerLinearProgram(const Vecd& c) : MixedIntegerLinearProgram(c,Matd(0,c.size()),Vecd(0),Matd(0,c.size()),Vecd(0),c.size()) {}

    MixedIntegerLinearProgram(int n) : MixedIntegerLinearProgram(Vecd(n),Matd(0,n),Vecd(0),Matd(0,n),Vecd(0),n) {}

    MixedIntegerLinearProgram() : MixedIntegerLinearProgram(0) {}

    MixedIntegerLinearProgram(const std::string& filename)
    {
        importJson(filename);
    }

    /// Set the objective function to c*x
    void setObjective(const Vecd& c) {
        assert(c.size() == dimension());
        this->c = c;
    }

    /// Returns the Integer Linear Program minimize c*x s.t. Ax <= b
    static MixedIntegerLinearProgram inequalityILP(const Vecd& c, const Matd& A, const Vecd& b) {
        Vecd d(0);
        Matd B(0,c.size());
        return MixedIntegerLinearProgram(c, A, b, B, d, c.size());
    }

    /// Returns the Integer Linear Program minimize c*x s.t. Bx <= d
    static MixedIntegerLinearProgram equalityILP(const Vecd& c, const Matd& B, const Vecd& d) {
        Vecd b(0);
        Matd A(0,c.size());
        return MixedIntegerLinearProgram(c, A, b, B, d, c.size());
    }

    /// Adds an inequality or equality constraint to the problem
    void addConstraint(const Vecd& coeffs, const double& rhs, const ConstraintType type = LEQ)
    {
        if (type == EQ)
        {
            p++;
            d.conservativeResize(p);
            B.conservativeResize(p, Eigen::NoChange);
            d[p-1] = rhs;
            B.row(p-1) = coeffs;
        }
        else
        {
            m++;
            b.conservativeResize(m);
            A.conservativeResize(m, Eigen::NoChange);
            b[m-1] = (type==LEQ)? rhs : -rhs;
            A.row(m-1) = (type==LEQ)? coeffs : -coeffs;
        }
    }

    /// Checks if x is a feasible point for the problem relaxation
    /*
    inline bool isFeasible(const Vecd& x)
    {
        assert(x.size()==n);
        return (x.array()>=0).all() &&
               ((A*x).array() <= b.array()).all() &&
               (B*x).isApprox(d);
    }
*/

    inline const uint dimension() {return n;}

    inline const uint nInequalities() {return m;}

    inline const uint nEqualities() {return p;}

    inline const uint nIntegerConstraints() {return n1;}

    inline const bool isIntegerConstrained(const uint& i) {assert(i>=0&&i<n); return i < n1;}

    inline const bool isLP() {return n1==0;}

    inline const bool isILP() {return n1==n;}

    inline const Vecd& costCoefficients() {return c;}

    inline const Matd& inequaltyMatrix() {return A;}

    inline const Vecd& inequalityVector() {return b;}

    inline const Matd& equalityMatrix() {return B;}

    inline const Vecd& equalityVector() {return d;}

    void exportJson(const std::string& filename)
    {
        std::ofstream o(filename);
        o << std::setw(4) << asJson() << std::endl;
        o.close();
    }

    Json asJson()
    {
        Json j = {
            {"c", c},
            {"A", A.rowwise()},
            {"b", b},
            {"B", B.rowwise()},
            {"d", d},
            {"n", n},
            {"m", m},
            {"p", p},
            {"n1", n1},
        };
        return j;
    }

    void importJson(const std::string& filename)
    {
        Json j = Json::parse(std::ifstream(filename));
        n = j["n"];
        m = j["m"];
        p = j["p"];
        n1 = j["n1"];
        A = toMatd(j["A"]);
        b = toVecd(j["b"]);
        B = toMatd(j["B"]);
        d = toVecd(j["d"]);
        c = toVecd(j["c"]);
        A.conservativeResize(m, n);
        B.conservativeResize(p, n);
        b.conservativeResize(m);
        d.conservativeResize(p);
        c.conservativeResize(n);
    }

    friend std::ostream& operator<< (std::ostream& out, const MixedIntegerLinearProgram& lp) {
        out << "================ Mixed Integer Linear Program ================" << std::endl;
        //out << &lp << std::endl;
        out << "n, m, p, n1: " << lp.n << ", " << lp.m << ", " << lp.p << ", " << lp.n1 << std::endl;
        out << "c:\n" << lp.c.transpose() << std::endl;
        out << "A:\n" << lp.A << std::endl;
        out << "b:\n" << lp.b.transpose() << std::endl;
        out << "B:\n" << lp.B << std::endl;
        out << "d:\n" << lp.d.transpose() << std::endl;
        out << "==============================================================" << std::endl;
        return out;
    }

private:
    uint n; // dimension
    uint m; // number of inequalities
    uint p; // number of equalities
    Vecd c; // objective function (gradient)
    Matd A; // inequality matrix (Ax <= b)
    Vecd b; // inequalities rhs (Ax <= b)
    Matd B; // equality matrix (Bx = d)
    Vecd d; // equalities rhs (Bx = d)
    uint n1; // number of integer constraints (x1,...,xn1 are integers)
};

class SparseMixedIntegerLinearProgram
{

public:
    enum ConstraintType {LEQ=0, GEQ=1, EQ=2};

    SparseMixedIntegerLinearProgram(const SVecd& c, const SMatd& A, const SVecd& b, const SMatd& B, const SVecd& d, const uint& n1)
        : n(c.size()), m(b.size()), p(d.size()), c(c), A(A), b(b), B(B), d(d), n1(n1)
    {
        assert(A.rows() == m && A.cols() == n);
        assert(B.rows() == p && B.cols() == n);
        assert(n1 <= n);
    }

    /// MILP of dimension n with n1 integral constraints
    SparseMixedIntegerLinearProgram(uint n, uint n1) : SparseMixedIntegerLinearProgram(SVecd(n),SMatd(0,n),SVecd(0),SMatd(0,n),SVecd(0),n1) {}

    /// Sets the i-th coefficient of the objective gradient to ci
    void setObjective(int i, double ci) {
        this->c.coeffRef(i) = ci;
    }

    /// Adds an inequality or equality constraint to the problem
    void addConstraint(const SVecd& coeffs, const double& rhs, const ConstraintType type = LEQ)
    {
        if (type == EQ)
        {
            p++;
            d.conservativeResize(p);
            B.conservativeResize(p, n);

            for (SVecdIter it(coeffs); it; ++it) {
                B.coeffRef(p-1, it.index()) = it.value();
            }

            d.coeffRef(p-1) = rhs;
        }
        else
        {
            m++;
            b.conservativeResize(m);
            A.conservativeResize(m, n);

            for (SVecdIter it(coeffs); it; ++it) {
                A.coeffRef(m-1, it.index()) = it.value();
            }

            b.coeffRef(m-1) = (type==LEQ)? rhs : -rhs;
        }
    }

    inline const uint dimension() {return n;}

    inline const uint nInequalities() {return m;}

    inline const uint nEqualities() {return p;}

    inline const uint nIntegerConstraints() {return n1;}

    inline const bool isLP() {return n1==0;}

    inline const bool isILP() {return n1==n;}

    inline const SVecd& objectiveGradient() {return c;}

    inline const SMatd& inequaltyMatrix() {return A;}

    inline const SVecd& inequalityVector() {return b;}

    inline const SMatd& equalityMatrix() {return B;}

    inline const SVecd& equalityVector() {return d;}

    inline void prune()
    {
        A.prune(EPSILON);
        b.prune(EPSILON);
        B.prune(EPSILON);
        d.prune(EPSILON);
        c.prune(EPSILON);
    }

    friend std::ostream& operator<< (std::ostream& out, const SparseMixedIntegerLinearProgram& lp) {
        out << "============= Sparse Mixed Integer Linear Program ============" << std::endl;
        //out << &lp << std::endl;
        out << "n, m, p, n1: " << lp.n << ", " << lp.m << ", " << lp.p << ", " << lp.n1 << std::endl;
        out << "c:\n" << lp.c.transpose() << std::endl;
        out << "A (" << lp.A.rows() << " x " << lp.A.cols() << ")\n"  << lp.A << std::endl;
        out << "b:\n" << lp.b.transpose() << std::endl;
        out << "B (" << lp.B.rows() << " x " << lp.B.cols() << ")\n" << lp.B << std::endl;
        out << "d:\n" << lp.d.transpose() << std::endl;
        out << "==============================================================" << std::endl;
        return out;
    }


private:
    uint n; // dimension
    uint m; // number of inequalities
    uint p; // number of equalities
    SVecd c; // objective function (gradient)
    SMatd A; // inequality matrix (Ax <= b)
    SVecd b; // inequalities rhs (Ax <= b)
    SMatd B; // equality matrix (Bx = d)
    SVecd d; // equalities rhs (Bx = d)
    uint n1; // number of integer constraints (x1,...,xn1 are integers)
};
}
