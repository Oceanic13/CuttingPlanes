
#include "MixedIntegerLinearProgram.h"
#include "Solvers.h"
#include "Utils.h"
#include "gtest/gtest.h"

namespace CP
{

using MILP = MixedIntegerLinearProgram;
using SMILP = SparseMixedIntegerLinearProgram;

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

    auto lp = MILP::inequalityILP(c, A, b);
    auto solver = ToblexSolver(lp);
    solver.solve();

    Vecd expected(2); expected << 1, 1.5;

    ASSERT_TRUE(solver.isOptimal());
    ASSERT_EQ(expected, solver.getOptimalSolution());
    ASSERT_DOUBLE_EQ(-1.5, solver.getOptimalValue());
}

TEST(ToblexTest, UnboundedTest)
{
    Matd A(1,2); A << 1, -1;
    Vecd b(1); b << 0;
    Vecd c(2); c << -1, -1;

    auto lp = MILP::inequalityILP(c, A, b);
    auto solver = ToblexSolver(lp);
    solver.solve();

    ASSERT_TRUE(solver.isUnbounded());
}

TEST(ToblexTest, SingletonTest)
{
    Matd A(1,2); A << 1, 1;
    Vecd b(1); b << 0;
    Vecd c(2);  c << 1, -2;

    auto lp = MILP::inequalityILP(c, A, b);
    auto solver = ToblexSolver(lp);
    solver.solve();

    Vecd expected(2); expected << 0, 0;

    ASSERT_TRUE(solver.isOptimal());
    ASSERT_EQ(expected, solver.getOptimalSolution());
    ASSERT_DOUBLE_EQ(0.0, solver.getOptimalValue());
}

TEST(ToblexTest, UnconstrainedProblemUnboundedTest)
{
    auto milp = MILP(5);
    Vecd c(5); c << 1, 3, -2, 0.5, -2.1; // vars 3, 5 to inf lead to arbitrarily small values
    milp.setObjective(c);

    auto s = ToblexSolver(milp);
    s.solve();

    ASSERT_TRUE(s.isUnbounded());
}

TEST(ToblexTest, UnconstrainedProblemWithOptimalSolutionTest)
{
    auto milp = MILP(5);
    Vecd c(5); c << 1, 3, 2, 0.5, 2.1; // since all cost coeffs are nonneg, the min. is eached by setting x=0
    milp.setObjective(c);

    Vecd expected(5); expected << 0, 0, 0, 0, 0;

    auto s = ToblexSolver(milp);
    s.solve();

    ASSERT_TRUE(s.isOptimal());
    ASSERT_EQ(expected, s.getOptimalSolution());
    ASSERT_DOUBLE_EQ(0, s.getOptimalValue());
}

TEST(ToblexTest, ConflictingEqualityConstraintsTest)
{
    auto milp = MILP(2);
    Vecd c(2); c << 1, -1;
    milp.setObjective(c);

    Vecd Bi(2);

    Bi << 1, 0;
    milp.addConstraint(Bi, 3, MILP::EQ); // x = 3

    Bi << 1, 0;
    milp.addConstraint(Bi, 5, MILP::EQ); // x = 5

    auto s = ToblexSolver(milp);
    s.solve();

    ASSERT_TRUE(s.isInfeasible());
}

TEST(ToblexTest, Phase1OptimalTest1)
{
    auto milp = MILP(2);
    Vecd c(2); c << 2, 1;
    milp.setObjective(c);

    Vecd Ai(2); Ai << -1, -1;
    milp.addConstraint(Ai, -3, MILP::LEQ);

    Vecd expected(2); expected << 0, 3;

    auto s = ToblexSolver(milp);
    s.solve();

    ASSERT_TRUE(s.isOptimal());
    ASSERT_EQ(expected, s.getOptimalSolution());
    ASSERT_DOUBLE_EQ(3, s.getOptimalValue());
}

TEST(ToblexTest, Phase1UnboundedTest1)
{
    auto milp = MILP(2);
    Vecd c(2); c << -2, -1;
    milp.setObjective(c);

    Vecd Ai(2); Ai << -1, -1;
    milp.addConstraint(Ai, -3, MILP::LEQ);

    auto s = ToblexSolver(milp);
    s.solve();

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

    Vecd expected(2); expected << 11./3., 1./3.;

    ASSERT_EQ(2, milp.dimension());
    ASSERT_EQ(2, milp.nInequalities());
    ASSERT_EQ(1, milp.nEqualities());

    auto s = ToblexSolver(milp);
    s.solve();

    //std::cout << s << std::endl;

    ASSERT_TRUE(s.isOptimal());
    ASSERT_TRUE((expected - s.getOptimalSolution()).squaredNorm() < 1e-9);
    ASSERT_DOUBLE_EQ(7.0, s.getOptimalValue());
}

TEST(CuttingPlaneTest, OptBookP331Test)
{
    Matd A(2,2); A << 3, 2, -3, 2;
    Vecd b(2); b << 6, 0;
    Vecd c(2); c << 0, -1;

    auto milp = MILP::inequalityILP(c, A, b);
    auto solver = CuttingPlanesSolver(milp);
    solver.solve();

    Vecd expected(2); expected << 1, 1;

    ASSERT_EQ(2, solver.getHistory().n_cuts());
    ASSERT_EQ(expected, solver.optimalSolution());
    ASSERT_DOUBLE_EQ(-1, solver.optimalValue());
}

TEST(CuttingPlaneTest, IntegerProblemIsInfeasibleTest)
{
    auto milp = MILP(2);
    Vecd v(2);
    v << -1, -1; milp.setObjective(v);
    v << 1, 0; milp.addConstraint(v, 0.25, MILP::GEQ); // x >= 0.25
    v << 1, 0; milp.addConstraint(v, 0.75, MILP::LEQ); // x <= 0.75
    v << 0, 1; milp.addConstraint(v, 0.25, MILP::GEQ); // y >= 0.25
    v << 0, 1; milp.addConstraint(v, 0.75, MILP::LEQ); // y <= 0.75

    auto solver = CuttingPlanesSolver(milp);
    solver.solve();

    ASSERT_EQ(1, solver.getHistory().n_cuts());
    ASSERT_TRUE(solver.isInfeasible());
}

TEST(ToblexTest, SparseMILPTest)
{
    auto milp = SMILP(2, 2);
    milp.setObjective(1, -1);

    auto a = SVecd(2);
    a.coeffRef(0) = 3;
    a.coeffRef(1) = 2;
    milp.addConstraint(a, 6, SMILP::LEQ);

    a.coeffRef(0) = -3;
    a.coeffRef(1) = 2;
    milp.addConstraint(a, 0, SMILP::LEQ);
}

int main(int argc, char** argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();

}
}
