#pragma once

#include <nlohmann/json.hpp>
#include <Eigen/Dense>

typedef double Num;
typedef Eigen::VectorX<Num> Vec;
typedef Eigen::MatrixX<Num> Mat;

using json = nlohmann::json;
