
#include "MixedIntegerLinearProgram.h"
#include "ToblexSolver.h"
#include "Utils.h"
#include "gtest/gtest.h"

namespace CP
{

using MILP = MixedIntegerLinearProgram;

TEST(UtilsTest, GetIntTest)
{
    ASSERT_EQ(10, getInt(10));
    ASSERT_EQ(1, getInt(1.9));
    ASSERT_EQ(0, getInt(0.2));
    ASSERT_EQ(-5, getInt(-4.75));
    ASSERT_EQ(-1, getInt(-0.0001));
}

TEST(UtilsTest, GetFracTest)
{
    ASSERT_DOUBLE_EQ(0, getFrac(42));
    ASSERT_DOUBLE_EQ(0.9, getFrac(1.9));
    ASSERT_DOUBLE_EQ(0.2, getFrac(0.2));
    ASSERT_DOUBLE_EQ(0.25, getFrac(-4.75));
    ASSERT_DOUBLE_EQ(0.9999, getFrac(-0.0001));
}

TEST(UtilsTest, IsIntTest)
{
    ASSERT_TRUE(isInt(42));
    ASSERT_FALSE(isInt(-0.0034));
    ASSERT_TRUE(isInt(-19));
    ASSERT_FALSE(isInt(20.5));
    ASSERT_FALSE(isInt(10.00001));
}

TEST(UtilsTest, removeColsTest)
{
    Matd A(3, 5);
    A << 1, 2, 3, 4, 5,
        6, 7, 8, 9, 10,
        11, 12, 13, 14, 15;

    removeCols(A, 1, 2);

    Matd B(3,3);
    B << 1, 4, 5,
        6, 9, 10,
        11, 14, 15;
    ASSERT_EQ(B, A);
}

TEST(ToblexTest, OptBookP331Test)
{
    Matd A(2,2); A << 3, 2, -3, 2;
    Vecd b(2); b << 6, 0;
    Vecd c(2); c << 0, -1;
    Vecd x(2);
    auto lp = MILP::inequalityILP(c, A, b);
    auto solver = ToblexSolver(lp);
    solver.solve(x);

    Vecd expected(2); expected << 1, 1.5;

    ASSERT_TRUE(solver.isOptimal());
    ASSERT_EQ(expected, x);
    ASSERT_DOUBLE_EQ(-1.5, solver.getOptimalValue());
}

TEST(ToblexTest, UnboundedTest)
{
    Matd A(1,2); A << 1, -1;
    Vecd b(1); b << 0;
    Vecd c(2); c << -1, -1;
    Vecd x(2);
    auto lp = MILP::inequalityILP(c, A, b);
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
    auto lp = MILP::inequalityILP(c, A, b);
    auto solver = ToblexSolver(lp);
    solver.solve(x);

    Vecd expected(2); expected << 0, 0;

    ASSERT_TRUE(solver.isOptimal());
    ASSERT_EQ(expected, x);
    ASSERT_DOUBLE_EQ(0.0, solver.getOptimalValue());
}

TEST(ToblexTest, UnconstrainedProblemUnboundedTest)
{
    auto milp = MILP(5);
    Vecd c(5); c << 1, 3, -2, 0.5, -2.1; // vars 3, 5 to inf lead to arbitrarily small values
    milp.setObjective(c);

    Vecd x(5);

    auto s = ToblexSolver(milp);
    s.solve(x);

    ASSERT_TRUE(s.isUnbounded());
}

TEST(ToblexTest, UnconstrainedProblemWithOptimalSolutionTest)
{
    auto milp = MILP(5);
    Vecd c(5); c << 1, 3, 2, 0.5, 2.1; // since all cost coeffs are nonneg, the min. is eached by setting x=0
    milp.setObjective(c);

    Vecd x(5);
    Vecd expected(5); expected << 0, 0, 0, 0, 0;

    auto s = ToblexSolver(milp);
    s.solve(x);

    ASSERT_TRUE(s.isOptimal());
    ASSERT_EQ(expected, x);
    ASSERT_DOUBLE_EQ(0, s.getOptimalValue());
}

TEST(ToblexTest, ConflictingEqualityConstraintsTest)
{
    auto milp = MILP(2);
    Vecd c(2); c << 1, -1;
    milp.setObjective(c);

    Vecd Bi(2);

    Bi << 1, 0;
    milp.addEqualityConstraint(Bi, 3); // x = 3

    Bi << 1, 0;
    milp.addEqualityConstraint(Bi, 5); // x = 5

    Vecd x(2);

    auto s = ToblexSolver(milp);
    s.solve(x);

    ASSERT_TRUE(s.isInfeasible());
}

TEST(ToblexTest, Phase1OptimalTest1)
{
    auto milp = MILP(2);
    Vecd c(2); c << 2, 1;
    milp.setObjective(c);

    Vecd Ai(2); Ai << -1, -1;
    milp.addInequalityConstraint(Ai, -3);

    Vecd x(2);
    Vecd expected(2); expected << 0, 3;

    auto s = ToblexSolver(milp);
    s.solve(x);

    ASSERT_TRUE(s.isOptimal());
    ASSERT_EQ(expected, x);
    ASSERT_DOUBLE_EQ(3, s.getOptimalValue());
}

TEST(ToblexTest, Phase1UnboundedTest1)
{
    auto milp = MILP(2);
    Vecd c(2); c << -2, -1;
    milp.setObjective(c);

    Vecd Ai(2); Ai << -1, -1;
    milp.addInequalityConstraint(Ai, -3);

    Vecd x(2);

    auto s = ToblexSolver(milp);
    s.solve(x);

    ASSERT_TRUE(s.isUnbounded());
}

TEST(ToblexTest, SomeRandomProblemTest)
{
    Matd A(2,2); A <<
        0, 2,
        -1,2;
    Vecd b(2); b << 2, -3;
    Matd B(1,2); B << 1, 1;
    Vecd d(1); d << 4;
    Vecd c(2);  c << 2, -1;
    auto milp = MILP(c, A, b, B, d, 2);

    Vecd x(2);
    Vecd expected(2); expected << 11./3., 1./3.;

    ASSERT_EQ(2, milp.dimension());
    ASSERT_EQ(2, milp.nInequalities());
    ASSERT_EQ(1, milp.nEqualities());

    auto s = ToblexSolver(milp);
    s.solve(x);

    std::cout << s << std::endl;

    ASSERT_TRUE(s.isOptimal());
    ASSERT_TRUE((expected - x).squaredNorm() < 1e-9);
    ASSERT_DOUBLE_EQ(7.0, s.getOptimalValue());
}

int main(int argc, char** argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();

}
}
