
#include "CuttingPlanes.h"
#include "MixedIntegerLinearProgram.h"
#include "ToblexSolver.h"

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

        auto cp = CuttingPlanes(milp);
        cp.solve();
        cp.exportJson(outFile);

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

    // Example from Book p35
    /*
    milp = MILP(toVecd({-1, -2, -3}));
    milp.addConstraint(toVecd({1,2.5,1.75}), 8, MILP::LEQ);
    milp.addConstraint(toVecd({1,0,0}), 3.5, MILP::LEQ);
    milp.addConstraint(toVecd({0,0,1}), 5.25, MILP::LEQ);
    milp.addConstraint(toVecd({0,3,2}), 6, MILP::LEQ);
    milp.exportJson("/Users/tobiaskohler/Uni/CuttingPlanes/milp.json");
*/
    
    auto cp = CuttingPlanes(milp);
    cp.solve();
    std::cout << cp.numberOfCuts() << std::endl;
    std::cout << cp.optimalSolution().transpose() << std::endl;
    cp.exportJson("/Users/tobiaskohler/Uni/CuttingPlanes/data/cp.json");

    return 0;
}
