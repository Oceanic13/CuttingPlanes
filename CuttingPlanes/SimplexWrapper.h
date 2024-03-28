#pragma once

#include "Typedefs.h"

class SimplexWrapper
{

public:
    explicit SimplexWrapper()
    {}

    void solve(Vec c, Mat A, Vec b);
};
