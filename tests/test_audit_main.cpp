// Heavy combinatorial audits (synthetic + subgroup probe). Own ctest target so
// Debug/Release/coverage always run them with a long timeout; ASan skips via -E.
#include <exception>
#include <iostream>

void test_synthetic_audit();
void test_subgroup_probe();

int main() {
  try {
    test_synthetic_audit();
    test_subgroup_probe();
  } catch (const std::exception& ex) {
    std::cerr << "exception: " << ex.what() << "\n";
    return 1;
  }
  std::cout << "audit tests passed\n";
  return 0;
}
