#pragma once

#include "Utils.h"

namespace CP
{
class LinearProgram
{
public:

    explicit LinearProgram(const Vecd& c, const Matd& A, const Vecd& b)
        : n(c.size()), m(b.size()), c(c), A(A), b(b)
    {
        assert(A.rows() == m && A.cols() == n);
    }

    LinearProgram() {}

    void addConstraint(const Vecd& a, const double& d)
    {
        m++;
        b.conservativeResize(m);
        A.conservativeResize(m, Eigen::NoChange);
        b[m-1] = d;
        A.row(m-1) = a;
    }

    inline const uint dimension() {return n;}

    inline const uint nInequalities() {return m;}

    inline const Vecd& costCoefficients() {return c;}

    inline const Matd& constraintMatrix() {return A;}

    inline const Vecd& constraintVector() {return b;}

    friend std::ostream& operator<< (std::ostream& out, const LinearProgram& lp) {
        out << "========================Linear Program========================" << std::endl;
        out << &lp << std::endl;
        out << "c:\n" << lp.c.transpose() << std::endl;
        out << "A:\n" << lp.A << std::endl;
        out << "b:\n" << lp.b.transpose() << std::endl;
        out << "==============================================================" << std::endl;
        return out;
    }

private:
    uint n;
    uint m;
    Vecd c;
    Matd A;
    Vecd b;

};
}
