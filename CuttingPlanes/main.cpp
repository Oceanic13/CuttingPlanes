#include "CuttingPlanes.h"
#include "SoPlexSolver.h"

int main()
{
    Mat A(2,2); A << 3, 2, -3, 2;
    Vec b(2); b << 6, 0;
    Vec c(2); c << 0, -1;

    auto cp = CuttingPlanes<SoPlexSolver>();
    cp.solve(c, -A, -b, {true,true});
    cp.export_json("/Users/tobiaskohler/Uni/CuttingPlanes/test.json");

    return 0;
}
