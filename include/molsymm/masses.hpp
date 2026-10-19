#pragma once

namespace molsymm {

/// Standard atomic weight for atomic number Z (1..118). Returns 0 for unknown Z.
double atomic_mass(int atomic_number);

} // namespace molsymm
