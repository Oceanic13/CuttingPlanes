#include "soplex.h"
#include "Typedefs.h"

namespace soplex {

class SoPlexSolver
{

public:
    SoPlexSolver() {}

    void solve(const Vec& c, const Mat& A, const Vec& b)
    {
        uint n = c.size();
        uint m = b.size();

        SoPlex mysoplex;

        /* set the objective sense */
        mysoplex.setIntParam(soplex::SoPlex::OBJSENSE, soplex::SoPlex::OBJSENSE_MINIMIZE);

        /* we first add variables */
        soplex::DSVector dummycol(0);
        for (uint i = 0; i < n; ++i) {
            mysoplex.addColReal(LPCol(c[i], dummycol, infinity, -infinity));
        }

        /* then constraints one by one */
        for (uint i = 0; i < m; ++i) {
            soplex::DSVector row(n);
            for (uint j = 0; j < n; ++j) {
                row.add(j, -A(i,j));
            }
            mysoplex.addRowReal(LPRow(-b[i], row, infinity));
        }

        /* write LP in .lp format */
        mysoplex.writeFileReal("/Users/tobiaskohler/Uni/CuttingPlanes/soplex.lp", NULL, NULL, NULL);

        /* solve LP */
        soplex::SPxSolver::Status stat;
        soplex::DVector prim(n);
        stat = mysoplex.optimize();

        /* get solution */
        if(stat == soplex::SPxSolver::OPTIMAL)
        {
            mysoplex.getPrimal(prim);
            std::cout << "LP solved to optimality.\n";
            std::cout << "Objective value is " << mysoplex.objValueReal() << ".\n";
            std::cout << "Primal solution is " << prim << ".\n";
        }
        else
        {
            std::cout << "Error: SoPlex returned with status " << stat << ".\n";
        }
    }
};
}
