#pragma once

#include <nlohmann/json.hpp>
#include <Eigen/Dense>
#include <iostream>
#include <vector>
#include <iterator>
#include <algorithm>

namespace CP
{
    template<typename T>
    using VecT = Eigen::VectorX<T>;
    template<typename T>
    using MatT = Eigen::MatrixX<T>;

    using Vecd = VecT<double>;
    using Matd = MatT<double>;

    using MatBlockd = Eigen::Block<Matd>;

    typedef nlohmann::json Json;

    inline constexpr bool isInt(const double& d) {return (d-trunc(d))<1e-9;}

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
