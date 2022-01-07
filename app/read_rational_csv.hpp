#pragma once

#include <cpu_ti/type.hpp>

#include <vector>
#include <array>
#include <string>

std::vector<std::array<cpu_ti::Scalar, 3>>
read_rational_csv(const std::string &inputFileName, std::vector<bool> &results);
