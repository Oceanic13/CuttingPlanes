#include "AlgLibSimplexSolver.h"
#include "CuttingPlanes.h"
#include "MixedIntegerLinearProgram.h"
#include "ToblexSolver.h"
#include "SoPlexSolver.h"

using namespace CP;
using MILP = MixedIntegerLinearProgram;

int main()
{

    // Example from Book
    Matd A(2,2); A << 3, 2, -3, 2;
    Vecd b(2); b << 6, 0;
    Vecd c(2); c << 0, -1;
    auto milp = MILP::inequalityILP(c, A, b);

    A << 1, 0, 0, 1;
    b << 1.5, 1.5;
    c << -1, -1;
    auto milp2 = MILP::inequalityILP(c, A, b);
    
    auto cp = CuttingPlanes(milp);

    cp.solve();
    cp.export_json("/Users/tobiaskohler/Uni/CuttingPlanes/cp.json");

    return 0;
}
