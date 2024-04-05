#include "AlgLibSimplexSolver.h"
#include "CuttingPlanes.h"
#include "MixedIntegerLinearProgram.h"
#include "ToblexSolver.h"
#include "SoPlexSolver.h"

using namespace CP;
using MILP = MixedIntegerLinearProgram;

int main(int argc, char* argv[])
{
    if (argc > 1) {
        std::cout << argv[1] << std::endl;

        return 0;
    }

    // Example from Book
    Matd A(2,2); A << 3, 2, -3, 2;
    Vecd b(2); b << 6, 0;
    Vecd c(2); c << 0, -1;
    auto milp = MILP::inequalityILP(c, A, b);
    milp.exportJson("/Users/tobiaskohler/Uni/CuttingPlanes/milp.json");
    
    auto cp = CuttingPlanes(milp);
    cp.solve();
    cp.exportJson("/Users/tobiaskohler/Uni/CuttingPlanes/cp.json");

    return 0;
}
