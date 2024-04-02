#pragma once

#include "Utils.h"

namespace CP
{
class MixedIntegerLinearProgram
{
public:

    /**
     * @brief MixedIntegerLinearProgram min c^Tx s.t. Ax <= b and Bx = d and x >= 0 and x_i integer for i = 1,...,n1 <= n
     * @param c Cost Coefficients
     * @param A Inequality Constraint Matrix
     * @param b Inequality Constraint Vector
     * @param n1 Number of integer constraints
     */
    explicit MixedIntegerLinearProgram(const Vecd& c, const Matd& A, const Vecd& b, const Matd& B, const Vecd& d, const int& n1)
        : dim(c.size()), nIneqs(b.size()), nEqs(d.size()), c(c), A(A), b(b), B(B), d(d), n1(n1)
    {
        assert(A.rows() == nIneqs && A.cols() == dim);
        assert(B.rows() == nEqs && B.cols() == dim);
        assert(n1 <= dim);
    }

    MixedIntegerLinearProgram(int n) : MixedIntegerLinearProgram(Vecd(n),Matd(0,n),Vecd(0),Matd(0,n),Vecd(0),n) {}

    MixedIntegerLinearProgram() : MixedIntegerLinearProgram(0) {}

    void setObjective(const Vecd& c) {
        assert(c.size() == dimension());
        this->c = c;
    }

    static MixedIntegerLinearProgram inequalityILP(const Vecd& c, const Matd& A, const Vecd& b) {
        Vecd d(0);
        Matd B(0,c.size());
        return MixedIntegerLinearProgram(c, A, b, B, d, c.size());
    }

    static MixedIntegerLinearProgram equalityILP(const Vecd& c, const Matd& B, const Vecd& d) {
        Vecd b(0);
        Matd A(0,c.size());
        return MixedIntegerLinearProgram(c, A, b, B, d, c.size());
    }

    void addInequalityConstraint(const Vecd& Ai, const double& bi)
    {
        nIneqs++;
        b.conservativeResize(nIneqs);
        A.conservativeResize(nIneqs, Eigen::NoChange);
        b[nIneqs-1] = bi;
        A.row(nIneqs-1) = Ai;
    }

    void addEqualityConstraint(const Vecd& Bi, const double& di)
    {
        nEqs++;
        d.conservativeResize(nEqs);
        B.conservativeResize(nEqs, Eigen::NoChange);
        d[nEqs-1] = di;
        B.row(nEqs-1) = Bi;
    }

    inline const uint dimension() {return dim;}

    inline const uint nInequalities() {return nIneqs;}

    inline const uint nEqualities() {return nEqs;}

    inline const uint nIntegerConstraints() {return n1;}

    inline const bool isIntegerConstrained(const uint& i) {assert(i>=0&&i<dim); return i < n1;}

    inline const Vecd& costCoefficients() {return c;}

    inline const Matd& inequaltyMatrix() {return A;}

    inline const Vecd& inequalityVector() {return b;}

    inline const Matd& equaltyMatrix() {return B;}

    inline const Vecd& equalityVector() {return d;}

    friend std::ostream& operator<< (std::ostream& out, const MixedIntegerLinearProgram& lp) {
        out << "========================Linear Program========================" << std::endl;
        out << &lp << std::endl;
        out << "c:\n" << lp.c.transpose() << std::endl;
        out << "A:\n" << lp.A << std::endl;
        out << "b:\n" << lp.b.transpose() << std::endl;
        out << "B:\n" << lp.B << std::endl;
        out << "d:\n" << lp.d.transpose() << std::endl;
        out << "==============================================================" << std::endl;
        return out;
    }

private:
    uint dim;
    uint nIneqs;
    uint nEqs;
    Vecd c;
    Matd A;
    Vecd b;
    Matd B;
    Vecd d;
    uint n1;
};
}
