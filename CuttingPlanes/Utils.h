#pragma once

#include <nlohmann/json.hpp>
#include <Eigen/Dense>
#include <iostream>
#include <vector>
#include <iterator>
#include <algorithm>

typedef double Num;
typedef Eigen::VectorX<Num> Vec;
typedef Eigen::MatrixX<Num> Mat;
typedef Eigen::Block<Mat, Eigen::Dynamic, Eigen::Dynamic, false> MatBlock;

using json = nlohmann::json;

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
