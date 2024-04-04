#pragma once

#include "soplex.h"
#include <nlohmann/json.hpp>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <vector>
#include <iterator>
#include <algorithm>

namespace CP
{
    constexpr double EPSILON = 1e-9;

    using Vecd = Eigen::VectorXd;
    using Matd = Eigen::MatrixXd;
    using SVecd = Eigen::SparseVector<double>;
    using SMatd = Eigen::SparseMatrix<double>;
    using Triplet = Eigen::Triplet<double>;
    using MatdBlock = Eigen::Block<double>;

    void removeCols(Eigen::MatrixXd& A, uint col, uint nColsToRemove)
    {
        uint nCols = A.cols();
        Eigen::MatrixXd B(A.rows(), nCols - nColsToRemove);
        B << A.leftCols(col), A.rightCols(nCols - (col + nColsToRemove));
        A = B;
    }

    void getNonZeroRows(const SMatd& A, const uint& col, std::vector<int>& nonZeroRows)
    {
        nonZeroRows.clear();
        for (int k = A.outerIndexPtr()[col]; k < A.outerIndexPtr()[col + 1]; ++k) {
            nonZeroRows.push_back(A.innerIndexPtr()[k]);
        }
    }

    void getNonZeroCols(const SMatd& A, const uint& row, std::vector<int>& nonZeroCols)
    {
        nonZeroCols.clear();
        for (SMatd::InnerIterator it(A,0); it; ++it) {
            nonZeroCols.push_back(it.col());
        }
    }

    void setRow(SMatd& A, const uint& row, const Vecd& v)
    {
        uint n = v.size();
        assert(A.cols() == n);
        for (auto j = 0u; j < n; ++j) {
            double a = v[j];
            if (a != 0.0) {
                A.insert(row, j) = a;
            }
        }
    }

    // R1 := R1 + c*R2
    void addRowToRow(SMatd& A, const uint& r1, const double& c, const uint& r2)
    {
        for (SMatd::InnerIterator it(A, r2); it; ++it) {
            A.coeffRef(r1, it.col()) += c * it.value();
        }
    }

    typedef nlohmann::json Json;

    inline constexpr int getInt(const double& d) {return std::floor(d);}

    inline constexpr double getFrac(const double& d) {return d - getInt(d);}

    inline constexpr bool isInt(const double& d) {return getFrac(d) < EPSILON;}

    template <typename T>
    std::ostream& operator<< (std::ostream& out, const std::vector<T>& v) {
        if (!v.empty()) {
            out << '[';
            std::copy(v.begin(), v.end(), std::ostream_iterator<T>(out, ", "));
            out << "\b\b]";
        } else {
            out << "[]";
        }
        return out;
    }
}

/*
# enable_testing()
# include(GoogleTest)
# find_package(Boost REQUIRED COMPONENTS system filesystem)
# add_executable(${PROJECT_NAME}-test
#         tests.cpp
#         )
*/
