#include "AlgLibSimplexSolver.h"
#include "CuttingPlanes.h"
#include "LinearProgram.h"
#include "ToblexSolver.h"
#include "SoPlexSolver.h"

using namespace CP;

int main()
{
    Matd A(2,2); A << 3, 2, -3, 2;
    Vecd b(2); b << 6, 0;
    Vecd c(2); c << 0, -1;
    auto lpOptBookP331 = LinearProgram(c, A, b);

    c << -1, -1;
    A.resize(1,2); A << 1, -1;
    b.resize(1); b << 0;
    auto lpUnbounded = LinearProgram(c, A, b);
    
    auto cp = CuttingPlanes<ToblexSolver>(lpOptBookP331);

    cp.solve({true,true});
    cp.export_json("/Users/tobiaskohler/Uni/CuttingPlanes/cp.json");

    return 0;
}
