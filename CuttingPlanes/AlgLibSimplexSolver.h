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

        alglib::minlpstate state;
        alglib::minlpcreate(n, state);
    }

    bool solve(Vec& primal, Vec& dual)
    {
        return false;
    }

    void add_constraint(const Vec& a, const double& d)
    {
    }

    void write_to_file()
    {
    }

private:
    uint n;
};
