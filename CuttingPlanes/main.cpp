#include "AlgLibSimplexSolver.h"
#include "CuttingPlanes.h"
#include "SoPlexSolver.h"

int main()
{
    Matd A(2,2); A << 3, 2, -3, 2;
    Vecd b(2); b << 6, 0;
    Vecd c(2); c << 0, -1;

    std::cout << A.cast<double>() << std::endl;
    
    auto cp = CuttingPlanes<SoPlexSolver>();
    cp.solve(c, -A, -b, {true,true});
    cp.export_json("/Users/tobiaskohler/Uni/CuttingPlanes/cp.json");

    return 0;
}
