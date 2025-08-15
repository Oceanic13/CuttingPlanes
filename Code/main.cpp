
//#include "CuttingPlanes.h"
#include "MixedIntegerLinearProgram.h"
#include "Solvers.h"

using namespace CP;
using MILP = MixedIntegerLinearProgram;

int main(int argc, char* argv[])
{
    if (argc == 3) {
        std::cout << "Cutting Planes" << std::endl;

        auto inFile = argv[1];
        auto outFile = argv[2];

        auto milp = MILP(inFile);

        std::cout << milp << std::endl;

        std::cout << "Cut me some slack!" << std::endl;

        auto cp = CuttingPlanesSolver(milp);
        cp.solve();
        cp.getHistory().exportJson(outFile);

        return 0;
    }

    std::cerr << "Running Example without args" << std::endl;

    // Example from Book p331
    Matd A(2,2); A << 4, 2, -3, 2;
    Vecd b(2); b << 6.5, -1;
    Vecd c(2); c << -0.5, -1;
    auto milp = MILP::inequalityILP(c, A, b);
    std::cout << milp << std::endl;
    milp.exportJson("/Users/tobiaskohler/Uni/CuttingPlanes/data/milp.json");
    
    auto cp = CuttingPlanesSolver(milp);
    cp.solve();
    auto& h = cp.getHistory();
    std::cout << h.n_cuts() << std::endl;
    std::cout << cp.optimalSolution().transpose() << std::endl;
    h.exportJson("/Users/tobiaskohler/Uni/CuttingPlanes/data/cp.json");

    return 0;
}
