#pragma once

#include "molsymm/molsymm.hpp"

#include <cmath>
#include <set>
#include <sstream>
#include <string>
#include <vector>

namespace molsymm::test {

// Closure, conjugacy bookkeeping, orthogonality (allowing real-E norm 2h), and
// integer L=1 / L=2 reductions. Invokes `fail(msg)` on each problem.
template <typename Fail>
void check_character_table(const molsymm::Analysis& a, const std::string& label, Fail&& fail) {
  const auto& ops = a.operations;
  const auto& ct = a.character_table;
  const int h = static_cast<int>(ops.size());
  if (h != a.order) {
    fail(label + " order mismatch");
    return;
  }
  for (const auto& x : ops) {
    for (const auto& y : ops) {
      const molsymm::Mat3 p = x.matrix * y.matrix;
      bool found = false;
      for (const auto& w : ops) {
        if ((w.matrix - p).cwiseAbs().maxCoeff() < 5e-3) {
          found = true;
          break;
        }
      }
      if (!found) {
        fail(label + " not closed");
        return;
      }
    }
  }
  const int nir = static_cast<int>(ct.irrep_names.size());
  std::vector<int> csize(ct.class_labels.size(), 0);
  for (const auto& op : ops) {
    if (op.class_index < 0 || op.class_index >= static_cast<int>(ct.class_labels.size())) {
      fail(label + " class index out of range");
      return;
    }
    csize[static_cast<std::size_t>(op.class_index)]++;
  }
  {
    std::vector<int> cls(static_cast<std::size_t>(h), -1);
    int ncls = 0;
    for (int i = 0; i < h; ++i) {
      if (cls[static_cast<std::size_t>(i)] >= 0) {
        continue;
      }
      for (int gi = 0; gi < h; ++gi) {
        const molsymm::Mat3 c = ops[static_cast<std::size_t>(gi)].matrix *
                                ops[static_cast<std::size_t>(i)].matrix *
                                ops[static_cast<std::size_t>(gi)].matrix.transpose();
        for (int j = 0; j < h; ++j) {
          if ((ops[static_cast<std::size_t>(j)].matrix - c).cwiseAbs().maxCoeff() < 5e-3) {
            cls[static_cast<std::size_t>(j)] = ncls;
          }
        }
      }
      ++ncls;
    }
    if (ncls != static_cast<int>(ct.class_labels.size())) {
      fail(label + " has " + std::to_string(ct.class_labels.size()) + " classes, expected " +
           std::to_string(ncls));
    }
  }
  for (std::size_t c = 0; c < csize.size(); ++c) {
    if (c < ct.class_sizes.size() && ct.class_sizes[c] != csize[c]) {
      fail(label + " class size mismatch for " + ct.class_labels[c]);
    }
  }
  static const std::set<std::string> cubic{"T", "Th", "Td", "O", "Oh", "I", "Ih"};
  if (cubic.count(a.point_group) > 0 &&
      std::set<std::string>(ct.class_labels.begin(), ct.class_labels.end()).size() !=
          ct.class_labels.size()) {
    std::string all;
    for (const auto& l : ct.class_labels) {
      all += " " + l;
    }
    fail(label + " duplicate class labels:" + all);
  }
  auto chi = [&](int ir, const molsymm::Operation& op) {
    return ct.characters[static_cast<std::size_t>(ir)][static_cast<std::size_t>(op.class_index)];
  };
  std::vector<double> norm(static_cast<std::size_t>(nir), 0.0);
  for (int i = 0; i < nir; ++i) {
    for (int j = i; j < nir; ++j) {
      double s = 0.0;
      for (const auto& op : ops) {
        s += chi(i, op) * chi(j, op);
      }
      if (i == j) {
        norm[static_cast<std::size_t>(i)] = s;
        const double k = s / h;
        if (std::abs(k - 1.0) > 1e-3 && std::abs(k - 2.0) > 1e-3) {
          std::ostringstream os;
          os << label << " <" << ct.irrep_names[static_cast<std::size_t>(i)] << "|"
             << ct.irrep_names[static_cast<std::size_t>(i)] << "> = " << k << " h";
          fail(os.str());
          return;
        }
      } else if (std::abs(s) > 1e-3 * h) {
        std::ostringstream os;
        os << label << " <" << ct.irrep_names[static_cast<std::size_t>(i)] << "|"
           << ct.irrep_names[static_cast<std::size_t>(j)] << "> = " << s;
        fail(os.str());
        return;
      }
    }
  }
  for (int rep = 1; rep <= 2; ++rep) {
    int dim_sum = 0;
    for (int i = 0; i < nir; ++i) {
      double s = 0.0;
      for (const auto& op : ops) {
        const double tr = op.matrix.trace();
        const double tr2 = (op.matrix * op.matrix).trace();
        s += chi(i, op) * ((rep == 1) ? tr : 0.5 * (tr * tr + tr2) - 1.0);
      }
      const double mult = s / norm[static_cast<std::size_t>(i)];
      const double r = std::round(mult);
      if (std::abs(mult - r) > 1e-3 || r < 0) {
        std::ostringstream os;
        os << label << " L=" << rep << " multiplicity of "
           << ct.irrep_names[static_cast<std::size_t>(i)] << " = " << mult;
        fail(os.str());
        return;
      }
      dim_sum += static_cast<int>(r) * ct.irrep_dims[static_cast<std::size_t>(i)];
    }
    if (dim_sum != 2 * rep + 1) {
      fail(label + " L=" + std::to_string(rep) + " dims sum to " + std::to_string(dim_sum));
      return;
    }
  }
}

} // namespace molsymm::test
