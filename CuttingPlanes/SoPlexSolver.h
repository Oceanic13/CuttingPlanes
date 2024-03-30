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
        soplex::DVector d(1);
        stat = soplex.optimize();

        // get solution
        if(stat == soplex::SPxSolver::OPTIMAL) {
            soplex.getPrimal(p);
            soplex.getDual(d);
            for (uint i = 0; i < n; ++i) {primal[i] = p[i];}
            for (uint i = 0; i < 1; ++i) {dual[i] = d[i];}

            std::cout << "================================================================" << std::endl;
            uint m = soplex.numRows();
            assert(n == soplex.numCols());

            for (auto i = 0u; i < n; ++i) {
                soplex::DSVector col(m);
                soplex.getColVectorReal(i,col);
                std::cout << "Col " << i << ": " << col << std::endl;
            }

            for (auto i = 0u; i < m; ++i) {
                soplex::DSVector row(n);
                soplex.getRowVectorReal(i,row);
                std::cout << "Row " << i << ": " << row << std::endl;
            }

            soplex::DVector lhs(m);
            soplex.getLhsReal(lhs);
            std::cout << "LHS: " << lhs << std::endl;

            soplex::DVector rhs(m);
            soplex.getRhsReal(rhs);
            std::cout << "RHS: " << rhs << std::endl;

            Matd C(m, n);
            for (auto i = 0u; i < m; ++i) {
                for (auto j = 0u; j < n; ++j) {
                    C(i,j) = soplex.coefReal(i,j);
                }
            }
            std::cout << "Coeff:\n" << C << std::endl;


            VecT<int> b(m);
            soplex.getBasisInd(b.data());
            std::cout << "Basis Indices: " << b.transpose() << std::endl;

            std::cout << "Primal = " << primal.transpose() << std::endl;
            std::cout << "Dual = " << dual.transpose() << std::endl;
            std::cout << "================================================================" << std::endl;

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
        soplex.writeBasisFile("/Users/tobiaskohler/Uni/CuttingPlanes/basis.bas", NULL, NULL, false);
    }

private:
    soplex::SoPlex soplex;
    uint n;

};
