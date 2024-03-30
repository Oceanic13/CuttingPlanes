#pragma once

#include "Utils.h"
#include "soplex.h"

class SoPlexSolver
{
public:
    explicit SoPlexSolver() {}
    
    void init(const Vecd& c, const Matd& A, const Vecd& b)
    {
        n = c.size(); // dimensionality
        uint m = b.size(); // number of inequalities
        soplex = soplex::SoPlex();

        // set the objective sense
        soplex.setIntParam(soplex::SoPlex::OBJSENSE, soplex::SoPlex::OBJSENSE_MINIMIZE);

        //we first add variables
        soplex::DSVector dummycol(0);
        for (uint i = 0; i < n; ++i) {
            soplex.addColReal( soplex::LPCol(c[i], dummycol, soplex::infinity, 0.0));
        }

        // then constraints one by one
        for (uint i = 0; i < m; ++i) {
            soplex::DSVector row(n);
            for (uint j = 0; j < n; ++j) {
                row.add(j, A(i,j));
            }
            soplex.addRowReal( soplex::LPRow(b[i], row,  soplex::infinity));
        }
    }
    
    bool solve(Vecd& primal, Vecd& dual)
    {
        // solve LP
        soplex::SPxSolver::Status stat;
        soplex::DVector p(n);
        stat = soplex.optimize();

        // get solution
        if(stat == soplex::SPxSolver::OPTIMAL) {
            soplex.getPrimal(p);
            for (uint i = 0; i < n; ++i) {
                primal[i] = p[i];
            }
            return true;
        } else {
            std::cerr << "Error: SoPlex returned with status " << stat << ".\n";
            return false;
        }
    }
    
    void add_constraint(const Vecd& a, const double& d)
    {
        soplex::DSVector row(n);
        for (uint j = 0; j < n; ++j) {
            row.add(j, a[j]);
        }
        soplex.addRowReal(soplex::LPRow(d, row,  soplex::infinity));
    }

    void write_to_file()
    {
        soplex.writeFileReal("/Users/tobiaskohler/Uni/CuttingPlanes/soplex.lp", NULL, NULL, NULL);
        soplex.writeBasisFile("/Users/tobiaskohler/Uni/CuttingPlanes/basis.lp", NULL, NULL, false);
    }

private:
    soplex::SoPlex soplex;
    uint n;

};
