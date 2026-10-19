#pragma once

#include <stdexcept>
#include <string>

namespace molsymm {

class Error : public std::runtime_error {
public:
  explicit Error(const std::string& message) : std::runtime_error(message) {}
};

} // namespace molsymm
