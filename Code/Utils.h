#pragma once

//#include "soplex.h"
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
    using SMatd = Eigen::SparseMatrix<double,Eigen::RowMajor>;
    using Triplet = Eigen::Triplet<double>;
    using Triplets = std::vector<Triplet>;
    using MatdBlock = Eigen::Block<double>;
    using SVecdIter = SVecd::InnerIterator;
    using SMatdIter = SMatd::InnerIterator;
    //using SMatdRowIter = Eigen::SparseMatrix<double,Eigen::ColMajor>::InnerIterator;
    //using SMatdColIter = Eigen::SparseMatrix<double,Eigen::RowMajor>::InnerIterator;

    Vecd toVecd(const std::vector<double> v) {
        Vecd r(v.size());
        for (auto i = 0; i < v.size(); ++i) r[i] = v[i];
        return r;
    }

    Matd toMatd(const std::vector<std::vector<double>> v) {
        if (v.size() == 0) return Matd(0,0);
        auto nRows = v.size();
        auto nCols = v[0].size();
        Matd r(nRows, nCols);
        for (auto i = 0; i < nRows; ++i)
            for (auto j = 0; j < nCols; ++j)
                r(i,j) = v[i][j];
        return r;
    }

    void insertEmptyCol(Matd& A, const uint& c)
    {
        assert(c >= 0 && c < A.cols());
        Eigen::MatrixXd B(A.rows(), A.cols() + 1);
        B.block(0, 0, A.rows(), c) = A.leftCols(c);
        B.block(0, c + 1, A.rows(), A.cols() - c) = A.rightCols(A.cols() - c);
        A = std::move(B);
    }

    void removeCols(Matd& A, uint col, uint nColsToRemove)
    {
        uint nCols = A.cols();
        Eigen::MatrixXd B(A.rows(), nCols - nColsToRemove);
        B << A.leftCols(col), A.rightCols(nCols - (col + nColsToRemove));
        A = B;
    }

    // R1 := R1 + c*R2
    void addRowToRow(Matd& A, const uint& r1, const double& c, const uint& r2)
    {
        A.row(r1) += c * A.row(r2);
    }

    // Make A[r, c] = 1 (should not be 0) and A[i, c] = 0 by row ops
    // s.t. the c-th column of A will be the r-th unit vector
    void gaussianPivot(Matd& A, const uint& r, const uint& c)
    {
        const double p = A(r, c);
        assert(A(r,c) != 0);

        // Divide row by p to make pivot 1
        A.row(r) /= p;
        A(r, c) = 1;

        // Make all other entries in column 0 by using elementary row operations
        for (uint i = 0; i < A.rows(); ++i) {
            if (i == r) {continue;}
            addRowToRow(A, i, -A(i,c), r);
            A(i,c) = 0;
        }
    }

    typedef nlohmann::json Json;

    /// integral part of a number (rounded down)
    inline constexpr int getInt(const double& d) {return std::floor(d);}

    /// fractional part of a number
    inline constexpr double getFrac(const double& d) {return d - getInt(d);}

    inline constexpr bool isInt(const double& d) {double f = getFrac(d); return f<EPSILON || f>1-EPSILON;}

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

    // Sparse Matrix Operations
    //inline void pruneRow(SMatd& A, int r) {A.prune([&r](int i, int j, float) {return i!=r;});}
    //inline void pruneCol(SMatd& A, int c) {A.prune([&c](int i, int j, float) {return j!=c;});}
    //inline void fillRow(SMatd& A, int i, double d) {for (int j=0;j<A.cols();++j) {A.coeffRef(i,j)=d;}}
    //inline void fillCol(SMatd& A, int j, double d) {for (int i=0;i<A.rows();++i) {A.coeffRef(i,j)=d;}}
    //inline void setRow(SMatd& A, int i, const SVecd& a) {pruneRow(A, i); for (SVecdIter it(a); it; ++it) {A.coeffRef(i, it.index()) = it.value();}}
    //inline void setCol(SMatd& A, int j, const SVecd& a) {pruneCol(A, j); for (SVecdIter it(a); it; ++it) {A.coeffRef(it.index(), j) = it.value();}}
    inline void fillSVecd(SVecd& vec, double d) {for (int i = 0; i < vec.size(); ++i) {vec.coeffRef(i) = d;}}
}
