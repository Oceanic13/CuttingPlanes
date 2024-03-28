
#include "CuttingPlanes.h"

int main()
{
    Mat A(3,2); A << 3, 2, -3, 2, 0, -1;
    Vec b(3); b << 6, 0, 0;
    Vec c(2); c << 0, -1;

    auto cp = CuttingPlanes(c, A, b, {true,true});
    cp.solve();
    cp.export_json("/Users/tobiaskohler/Uni/CuttingPlanes/test.json");

    return 0;
}
