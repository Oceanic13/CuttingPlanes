#pragma once

#include "Utils.h"
#include "alglib/stdafx.h"
#include "alglib/optimization.h"

class AlgLibSimplexSolver
{
public:
    explicit AlgLibSimplexSolver() {}

    void init(const Vec& c, const Mat& A, const Vec& b)
    {
        n = c.size();
        uint m = b.size();

        // Setup Min MP with n variables
        state = alglib::minlpstate();
        alglib::minlpcreate(n, state);

        // Set Objective Function
        alglib::minlpsetcost(state, toVec(c));

        // Set Box Constraints 0 <= x <= inf
        alglib::minlpsetbcall(state, 0, INFINITY);

        // Set Constraints Ax >= b
        alglib::real_2d_array a;
        a.setlength(m, n+1);
        for (auto i = 0; i < m; ++i) {
            for (auto j = 0; j < n; ++j) {
                a[i][j] = A(i, j);
            }
            a[i][n] = b[i];
        }
        alglib::integer_1d_array ct;
        ct.setlength(m);
        for (auto i = 0; i < m; ++i) {
            ct[i] = 1;
        }
        alglib::minlpsetlc(state, a, ct);
    }

    bool solve(Vec& primal, Vec& dual)
    {
        alglib::minlpoptimize(state);
        alglib::real_1d_array x;
        alglib::minlpreport r;
        alglib::minlpresults(state, x, r);

        for (auto i = 0; i < n; ++i) {
            primal[i] = x[i];
        }

        return true;
    }

    void add_constraint(const Vec& a, const double& d)
    {
        alglib::minlpaddlc2dense(state, toVec(a), d, INFINITY);
    }

    void write_to_file()
    {
    }

private:
    uint n;
    alglib::minlpstate state;

    alglib::real_1d_array toVec(const Vec& v)
    {
        alglib::real_1d_array r;
        r.setcontent(v.size(), v.data());
        return r;
    }
};
