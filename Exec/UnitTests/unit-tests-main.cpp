#include <iosfwd> // for string

#include "gtest/gtest_pred_impl.h" // for AddGlobalTestEnvironment, InitGoo...

#include "AmrexTestEnv.H" // for AmrexTestEnv

// Necessary as it's used in other source files
std::string inputs_name;

//! Global instance of the environment (for access in tests)
pelec_tests::AmrexTestEnv* utest_env = nullptr;

int
main(int argc, char** argv)
{
  ::testing::InitGoogleTest(&argc, argv);

  utest_env = new pelec_tests::AmrexTestEnv(argc, argv);
  ::testing::AddGlobalTestEnvironment(utest_env);

  return RUN_ALL_TESTS();
}
