
#include "LinearProgram.h"
#include "ToblexSolver.h"
#include "Utils.h"
#include "gtest/gtest.h"

namespace CP
{

TEST(ToblexTest, OptBookP331Test)
{
    Matd A(2,2); A << 3, 2, -3, 2;
    Vecd b(2); b << 6, 0;
    Vecd c(2); c << 0, -1;
    Vecd x(2);
    auto lp = LinearProgram(c, A, b);
    auto solver = ToblexSolver(lp);
    solver.solve(x);

    Vecd expected(2); expected << 1, 1.5;

    ASSERT_TRUE(solver.isOptimal());
    ASSERT_EQ(expected, x);
}

TEST(ToblexTest, UnboundedTest)
{
    Matd A(1,2); A << 1, -1;
    Vecd b(1); b << 0;
    Vecd c(2); c << -1, -1;
    Vecd x(2);
    auto lp = LinearProgram(c, A, b);
    auto solver = ToblexSolver(lp);
    solver.solve(x);

    ASSERT_TRUE(solver.isUnbounded());
}

TEST(ToblexTest, SingletonTest)
{
    Matd A(1,2); A << 1, 1;
    Vecd b(1); b << 0;
    Vecd c(2);  c << 1, -2;
    Vecd x(2);
    auto lp = LinearProgram(c, A, b);
    auto solver = ToblexSolver(lp);
    solver.solve(x);

    Vecd expected(2); expected << 0, 0;

    ASSERT_TRUE(solver.isOptimal());
    ASSERT_EQ(expected, x);
}

int main(int argc, char** argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();

}
}
