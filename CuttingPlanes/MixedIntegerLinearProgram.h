#pragma once

#include "Utils.h"

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

    inline bool isFeasible(const Vecd& x)
    {
        assert(x.size()==n);
        return (x.array()>=0).all() &&
               ((A*x).array() <= b.array()).all() &&
               (B*x).isApprox(d);
    }

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

    inline const Matd& equaltyMatrix() {return B;}

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
    uint n;
    uint m;
    uint p;
    Vecd c;
    Matd A;
    Vecd b;
    Matd B;
    Vecd d;
    uint n1;
};
}
