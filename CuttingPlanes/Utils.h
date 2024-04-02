#pragma once

#include <nlohmann/json.hpp>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <vector>
#include <iterator>
#include <algorithm>

namespace CP
{
    using Vecd = Eigen::VectorXd;
    using Matd = Eigen::MatrixXd;
    using SVecd = Eigen::SparseVector<double>;
    using SMatd = Eigen::SparseMatrix<double>;
    using Triplet = Eigen::Triplet<double>;
    using MatdBlock = Eigen::Block<double>;

    void removeRow(Eigen::MatrixXd& matrix, unsigned int rowToRemove)
    {
        unsigned int numRows = matrix.rows()-1;
        unsigned int numCols = matrix.cols();

        if( rowToRemove < numRows )
            matrix.block(rowToRemove,0,numRows-rowToRemove,numCols) = matrix.block(rowToRemove+1,0,numRows-rowToRemove,numCols);

        matrix.conservativeResize(numRows, Eigen::NoChange);
    }

    void removeCols(Eigen::MatrixXd& A, uint col, uint nColsToRemove)
    {
        uint nCols = A.cols();
        Eigen::MatrixXd B(A.rows(), nCols - nColsToRemove);
        B << A.leftCols(col), A.rightCols(nCols - (col + nColsToRemove));
        A = B;
    }

    typedef nlohmann::json Json;

    inline constexpr int getInt(const double& d) {return std::floor(d);}

    inline constexpr double getFrac(const double& d) {return d - getInt(d);}

    inline constexpr bool isInt(const double& d) {return getFrac(d) < 1e-9;}

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
